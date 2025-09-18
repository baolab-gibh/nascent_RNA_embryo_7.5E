#!/usr/bin/env python3
# File: embryo_development_stage.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Nov 27, 2024
# Updated:

import json
import logging
import pickle
import shutil
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from copy import deepcopy
from pathlib import Path
from typing import List

import click
import anndata as adt
import numpy as npy
import polars as pls
import scanpy as scp
import matplotlib.pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as fun
from torch.optim import Adam
from torch.utils.data import DataLoader, Dataset

from sklearn.base import BaseEstimator
from sklearn.cluster import AffinityPropagation
from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline

from scipy.sparse import csr_matrix
from scipy.stats import pearsonr, spearmanr, kendalltau

from tqdm import tqdm


class LogManager(logging.Logger):
    def __init__(self, name, level=logging.INFO, logstream: bool = True, logfile: str = ""):
        super(LogManager, self).__init__(name)
        fmt = logging.Formatter(
            "{levelname: >8} | {asctime} | {name: ^20} | {message}", style="{", datefmt="%y-%m-%d,%H:%M:%S"
        )
        if logstream:
            self._add_handler(logging.StreamHandler(), level, fmt)

        if logfile:
            self._add_handler(logging.FileHandler(logfile), level, fmt)

    def _add_handler(self, hdl, lvl, fmt):
        hdl.setLevel(lvl)
        hdl.setFormatter(fmt)
        self.addHandler(hdl)


class ExpDataSet(Dataset):
    def __init__(self, X: pls.DataFrame, y=None):
        super(ExpDataSet, self).__init__()
        self.x_mat, self.y_vec = X, y

    def __len__(self):
        n_row, _ = self.x_mat.shape
        if self.y_vec is not None:
            assert n_row == self.y_vec.shape[0]
        return n_row

    def __getitem__(self, idx):
        x = self.x_mat[idx, :]
        if self.y_vec is not None:
            y = self.y_vec[idx]
            return x, y
        return x

    def __getitems__(self, idx: list):
        return self.__getitem__(idx)


class DatasetAligner(BaseEstimator):
    def __init__(
        self, max_missing_rate: float = 0.5, fill_value: float | int | str = 0, logman: LogManager = LogManager("DatasetAligner")
    ) -> None:
        self.feature_order, self.feature_status = None, None
        self.fill_value = fill_value
        self.max_missing_rate = max_missing_rate
        self.logman = logman

    def fit(self, X: pls.DataFrame, y=None):
        """Learn the feature order of the dataset."""
        self.feature_order = tuple(X.columns)
        self.feature_status = pls.DataFrame(
            dict(feature=self.feature_order, mean=X.mean().to_numpy().squeeze(), median=X.median().to_numpy().squeeze())
        )

        return self

    def transform(self, X: pls.DataFrame):
        """Align the dataset to the feature order of the fitted model."""
        if self.feature_order is None:
            raise ValueError("DatasetAligner must be fitted first.")

        # handle perfect match
        if tuple(X.columns) == self.feature_order:
            return X.to_numpy().astype("float32")

        # handle extra features
        extra_features = [x for x in X.columns if x not in self.feature_order]
        self.logman.warning(f"Found {len(extra_features)} extra features, e.g., {extra_features[:5]}")

        # handle missing features
        misses = [x for x in self.feature_order if x not in X.columns]
        n_misses, n_required = len(misses), len(self.feature_order)
        if n_misses / float(n_required) > self.max_missing_rate:
            self.logman.warning(f"Too many missing features: {n_misses} out of {n_required}. E.g. {misses[:5]}.")

        for per_mis in misses:
            if self.feature_status is None:
                X[per_mis] = 0
            elif isinstance(self.fill_value, (int, float)):
                X[per_mis] = self.fill_value
            elif isinstance(self.fill_value, str) and self.fill_value in ["mean", "median"]:
                X[per_mis] = self.feature_status.filter(feature=per_mis)[self.fill_value]
            else:
                self.logman.warning(f"Unknown fill value type: {type(self.fill_value)}. Using 0 instead.")
                X[per_mis] = 0

        # reorder and return.
        return X.select(self.feature_order).to_numpy().astype("float32")

    def fit_transform(self, X: pls.DataFrame, y=None):
        """Fit and transform the dataset."""
        self.fit(X)
        return self.transform(X)


class FeatureSelector(BaseEstimator):
    def __init__(
        self, method="pearson", threshold: float = 0.3, max_p_adj: float = 0.05, max_missing_rate: float = 0.75,
        missing_values: float = 0, logman: LogManager = LogManager("FeatureSelector")
    ):
        self.selected_features = None

        self.method = method
        self.threshold = threshold
        self.max_p_adj = max_p_adj

        self.max_missing_rate = max_missing_rate
        self.missing_values = missing_values

        self.logman = logman

    def fit(self, X: pls.DataFrame, y=None):
        n_row, n_col = X.shape
        lmr_cols = ( # low missing rate columns
            (X.with_columns(pls.all() == self.missing_values).sum() / float(n_row))
            .transpose(include_header=True)
            .filter(pls.col("column_0") <= self.max_missing_rate)
            .get_column("column")
            .to_list()
        )
        self.logman.info(f"Low missing rate columns: {len(lmr_cols)} out of {n_col} columns, e.g. {lmr_cols[:5]}")

        if y is None:
            self.selected_features = lmr_cols
        else:
            # if the y is available, to calculate the correlation and remove columns based on self.threshold, self.max_p_adj
            if len(y) != n_row:
                raise ValueError(f"The number of elments in y ({len(y)}) doesn't match the number of rows in X ({n_row}).")

            corr_method = {"spearman": spearmanr, "kendall": kendalltau, "pearson": pearsonr}.get(self.method, pearsonr)
            if self.method != "pearson":
                self.logman.warning("Unknown correlation method: {self.method}. Using pearson instead.")

            x_mat = X.select(lmr_cols)
            n_row, n_col = x_mat.shape

            corr_dict = {"feature_id": [], "p_value": [], "correlation": []}
            for per_col in lmr_cols:
                corr = corr_method(x_mat[:, per_col].to_numpy(), y)
                corr_dict["feature_id"].append(per_col)
                corr_dict["p_value"].append(corr.pvalue)
                corr_dict["correlation"].append(corr.statistic)

            self.selected_features = (
                pls.DataFrame(corr_dict)
                .sort("p_value")
                .with_row_index("rank")
                .with_columns(((pls.col("p_value") * pls.col("rank")) / pls.len()).alias("p_value_adj"))
                .filter(pls.col("correlation").abs() >= self.threshold, pls.col("p_value_adj") <= self.max_p_adj)
                .get_column("feature_id")
                .to_list()
            )

        return self

    def transform(self, X: pls.DataFrame):
        if self.selected_features is None:
            raise ValueError("FeatureSelector must be fitted first.")

        return X.select(self.selected_features)


class VAE(nn.Module):
    '''Variational Autoencoder.'''
    def __init__(
        self, in_dims: int, n_classes: int | None = None, hidden_dims: int = 32, latent_dims: int = 8
    ) -> None:
        super(VAE, self).__init__()
        self.in_dims = in_dims
        self.n_classes = n_classes
        self.hidden_dims = hidden_dims
        self.latent_dims = latent_dims

        self.enc_fc1 = nn.Linear(in_dims, hidden_dims)
        self.enc_fc2 = nn.Linear(hidden_dims, hidden_dims)
        self.enc_mu = nn.Linear(hidden_dims, latent_dims)
        self.enc_logvar = nn.Linear(hidden_dims, latent_dims)

        self.dec_fc1 = nn.Linear(latent_dims, hidden_dims)
        self.dec_fc2 = nn.Linear(hidden_dims, hidden_dims)
        self.dec_fc3 = nn.Linear(hidden_dims, in_dims)

        if n_classes is None:
            self.pre_fc1 = None
        else:
            self.pre_fc1 = nn.Linear(latent_dims, n_classes)

        self.relu = nn.ReLU()
        self.drop = nn.Dropout(0.5)
        self.sigmoid = nn.Sigmoid()
        self.leaky_relu = nn.LeakyReLU(0.05)

    @property
    def has_classifier(self):
        return self.pre_fc1 is not None

    def encode(self, x):
        '''Encode the input vector.'''
        x = self.leaky_relu(self.enc_fc1(x))
        x = self.leaky_relu(self.enc_fc2(x))
        x_enc = self.drop(x)
        mu = self.leaky_relu(self.enc_mu(x_enc))
        logvar = self.sigmoid(self.enc_logvar(x_enc))

        return x_enc, mu, logvar

    def latentz(self, mu, logvar):
        '''The latent layer.'''
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        z = eps * std + mu

        return z

    def decode(self, x):
        '''Decode the latent vector.'''
        x = self.leaky_relu(self.dec_fc1(x))
        x = self.leaky_relu(self.dec_fc2(x))
        x = self.drop(self.dec_fc3(x))

        return x

    def classify(self, x):
        '''Classify the latent vector.'''
        if self.n_classes is not None and self.pre_fc1 is not None:
            x = self.relu(self.pre_fc1(x))
            x = self.drop(x)
            y_hat = torch.round(x)
            y_prob = fun.softmax(x, dim=1)
            return y_hat, y_prob
        
        return None, None

    def forward(self, x):
        '''Propagate the input vector through the network.'''
        _, mu, logvar = self.encode(x)
        mu_hat = self.latentz(mu, logvar)
        x_dec = self.decode(mu_hat)
        y_pred = self.classify(mu_hat)

        return mu, logvar, mu_hat, x_dec, y_pred


class VAEMachine(BaseEstimator):
    def __init__(
        self, model=None, transform_into: str = "x_enc", logman: LogManager = LogManager("VAEMachine"), **kwargs
    ) -> None:
        super(VAEMachine, self).__init__()

        if model is None:
            self.in_dims = kwargs.get("in_dims", None)
            self.out_dims = kwargs.get("in_dims", None)
            self.n_classes = kwargs.get("n_classes", None)
            self.hidden_dims = kwargs.get("hidden_dims", 32)
            self.latent_dims = kwargs.get("latent_dims", 16)

            self.epochs = kwargs.get("epochs", 100)
            self.batch_size = kwargs.get("batch_size", 32)
            self.learning_rate = kwargs.get("learning_rate", 1e-3)

            self.model = VAE(self.in_dims, self.n_classes, self.hidden_dims, self.latent_dims)
        else:
            self.model = torch.load(model)

            self.in_dims = self.model.in_dims
            self.out_dims = self.model.in_dims
            self.n_classes = self.model.n_classes
            self.hidden_dims = self.model.hidden_dims
            self.latent_dims = self.model.latent_dims

        self.device = kwargs.get("device", "cpu")
        self.logman = logman
        self.transform_into = transform_into

    def get_params(self, deep=True):
        return {k: v for k, v in self.__dict__.items()}

    def set_params(self, **params):
        for p, v in params.items():
            setattr(self, p, v)
        return self

    def to(self, device):
        self.model.to(device)
        return self

    def fit(self, X, y=None):
        # Optimizers
        optimizer = Adam(self.model.parameters(), lr=self.learning_rate)

        # Loss functions
        mse_loss_fun = nn.MSELoss()
        kld_loss_fun = nn.KLDivLoss(reduction="batchmean", log_target=True)
        if y is not None and self.model.has_classifier:
            cel_loss_fun = nn.CrossEntropyLoss()
        else:
            cel_loss_fun = None

        # Turn on the training mode
        if not self.model.training: self.model.train()

        if y is not None:
            exp_dataset = ExpDataSet(X, y)
        else:
            exp_dataset = ExpDataSet(X)

        data_loader = DataLoader(exp_dataset, batch_size=self.batch_size, shuffle=True)
        for _ in tqdm(range(self.epochs)):
            for _, xy in enumerate(data_loader):
                if isinstance(xy, tuple):
                    per_x, per_y = xy
                else:
                    per_x, per_y = xy, None
                per_x = per_x.to(self.device)
                mu, _, mu_hat, x_dec, _ = self.model.to(self.device)(per_x.to(self.device)) # Forward pass

                # Losses
                loss_rec = mse_loss_fun(x_dec, per_x)
                loss_kld = kld_loss_fun(mu, mu_hat)
                loss = loss_rec + loss_kld
                if cel_loss_fun is not None:
                    loss_cls = cel_loss_fun(mu_hat, per_y)
                    loss += loss_cls

                loss.backward() # Backward pass
                optimizer.step() # Update
                optimizer.zero_grad()

        return self

    def transform(self, X):
        '''Transform the input.

        This method is designed for the sklearn.Pipeline. For the usage of single VAE, please use the `compress`,
        `generate` and `reconstruct` methods.
        '''
        with torch.no_grad():
            if self.model.training: self.model.eval()
            if self.transform_into == "x_enc":
                return self.compress(X)
            else:
                if self.transform_into != "x_dec":
                    self.logman.error(f"Invalid transform_into: {self.transform_into}, will transform into decoded X.")
                return self.reconstruct(X)

    def predict(self, X):
        '''Predict the output.'''
        y_pred, _ = self.classify(X)

        return y_pred

    def predict_proba(self, X):
        '''Predict the output probabilities.'''
        _, y_prob =  self.classify(X)

        return y_prob

    def compress(self, X):
        '''Compress the input vector to generate the latent vector.'''
        if self.model.training: self.model.eval()

        with torch.no_grad():
            x_enc, *_ = self.model.encode(X)

        return x_enc

    def generate(self, mu = None, logvar = None):
        '''Generate the output vector from the latent vector.'''
        if mu is None:
            mu = torch.randn((1, self.latent_dims))

        if logvar is None:
            logvar = torch.randn((1, self.latent_dims))

        if self.model.training: self.model.eval()

        with torch.no_grad():
            mu_hat = self.model.latentz(mu, logvar)
            x_dec = self.model.decode(mu_hat)

        return x_dec

    def reconstruct(self, X):
        '''Reconstruct the input vector.'''
        if self.model.training: self.model.eval()

        with torch.no_grad():
            *_, x_dec, _ = self.model.forward(X)

        return x_dec

    def classify(self, X):
        '''Classify the input vector.'''
        y_hat, y_prob = None, None
        if self.model.has_classifier:
            if self.model.training: self.model.eval()

            with torch.no_grad():
                y_hat, y_prob = self.model.classify(X)
                return y_hat, y_prob

        self.logman.warning("No classification layer.")
        return y_hat, y_prob


def preprocess(
    anndata: adt.AnnData, batch_key: str | None = None, min_cells: int = 3, min_genes: int = 200,
    max_mito_pct: float = 10, max_dbl_score: float = 0.2, n_comps: int = 50, resolution: float = 0.1,
    logman: LogManager = LogManager("Preprocess")
):
    '''Preprocess input data.'''
    anndata = anndata.copy()

    logman.debug("Making cell names unique ...")
    anndata.obs_names_make_unique()

    logman.debug("Creating QC metrics variables ...")
    anndata.var["mt"] = anndata.var_names.str.startswith(("mt-", "MT-", "Mt-"))
    anndata.var["ribo"] = anndata.var_names.str.startswith(("RPL", "RPS", "Rpl", "Rps"))
    anndata.var["hb"] = anndata.var_names.str.contains("^HB[^(P)]")
    scp.pp.calculate_qc_metrics(anndata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    logman.debug("Estimating doublets ...")
    scp.pp.scrublet(anndata, batch_key=batch_key)

    if anndata.X is not None and anndata.layers.get("raw_counts", None) is None:
        logman.debug("Save raw counts into layers")
        anndata.layers["raw_counts"] = deepcopy(anndata.X)

    logman.debug(f"Filtering cells and genes ...")
    scp.pp.filter_cells(anndata, min_genes=min_genes)
    scp.pp.filter_genes(anndata, min_cells=min_cells)

    logman.debug("Filtering mitochondrial genes, and doublets ...")
    kept = (anndata.obs["pct_counts_mt"] <= max_mito_pct) & (anndata.obs["doublet_score"] <= max_dbl_score)
    anndata = anndata[kept]

    logman.debug(f"Normalizing and log-transforming data ...")
    scp.pp.normalize_total(anndata, target_sum=1e6)
    scp.pp.log1p(anndata)

    logman.debug("Identifying highly variable genes ...")
    scp.pp.highly_variable_genes(anndata, batch_key=batch_key)

    logman.debug("Computing PCA ...")
    scp.pp.pca(anndata, n_comps=n_comps)

    logman.debug("Computing UMAP ...")
    scp.pp.neighbors(anndata)
    scp.tl.umap(anndata)

    logman.debug("Computing diffusion map ...")
    scp.tl.diffmap(anndata)

    logman.debug("Identifying clusters ...")
    scp.tl.louvain(anndata, resolution=resolution)

    return anndata


def affinity_propagation(tbl: pls.DataFrame, **kwargs):
    '''Affinity propagation clustering.'''
    group_by = kwargs.get("group_by", [])
    group_id = kwargs.get("group_id", "")
    min_cells_pmc = kwargs.get("min_cells_pmc", 1)

    if isinstance(group_by, str): group_by = [group_by]
    non_data_cols = ["cell_barcode"] + group_by
    data_cols = [x for x in tbl.columns if x not in non_data_cols]

    ap_keys = ["damping", "max_iter", "convergence_iter", "copy", "preference", "affinity", "verbose", "random_state"]
    ap_kwargs = {k: v for k, v in kwargs.items() if k in ap_keys}
    af = AffinityPropagation(**ap_kwargs, copy=True).fit(tbl[:, data_cols])

    n_rows, *_ = tbl.shape
    is_cls_center = [x in af.cluster_centers_indices_ for x in range(n_rows)]
    mc_raw_barcode = [f"{group_id}-MC{x:06d}" for x in af.labels_]
    cls_tbl = (
        tbl[:, non_data_cols]
        .with_columns(pls.Series(mc_raw_barcode).alias("mc_raw_barcode"), pls.Series(is_cls_center).alias("mc_raw_cc"))
        .with_columns(
            pls.when(pls.col("mc_raw_cc"))
            .then(pls.col("cell_barcode"))
            .otherwise(pls.col("cell_barcode").sort_by("mc_raw_cc").last().over("mc_raw_barcode")).alias("mc_raw_cc_cb")
        )
        .with_columns(pls.col("mc_raw_barcode").is_unique().alias("mc_raw_is_singleton"))
        .with_columns((pls.count().over("mc_raw_barcode") <= min_cells_pmc).alias("mc_raw_belongs_to_tiny_cluster"))
    )

    # Determine the new cluster of samples classified into tiny clusters
    tiny_cluster_sample = cls_tbl.filter(pls.col("mc_raw_belongs_to_tiny_cluster"))["cell_barcode"]
    normal_center = cls_tbl.filter((~pls.col("mc_raw_belongs_to_tiny_cluster") & pls.col("mc_raw_cc")))["cell_barcode"]
    selected_sample = pls.concat([tiny_cluster_sample, normal_center])
    new_clusters = (
        tbl.filter(pls.col("cell_barcode").is_in(selected_sample))
        .select(["cell_barcode"] + data_cols) .transpose(column_names="cell_barcode")
        .corr()
        .with_columns(selected_sample.alias("cell_barcode"))
        .unpivot(list(selected_sample), index="cell_barcode", variable_name="refer_barcode")
        .filter(pls.col("refer_barcode").is_in(normal_center), pls.col("cell_barcode").is_in(tiny_cluster_sample))
        .with_columns(
            pls.col("value").max().over("cell_barcode").alias("mc_regroup_corr"),
            pls.col("refer_barcode").sort_by("value").last().over("cell_barcode").alias("mc_regroup_cc_cb"))
        .filter(pls.col("mc_regroup_cc_cb") == pls.col("refer_barcode"))
        .select("cell_barcode", "mc_regroup_corr", "mc_regroup_cc_cb")
    )

    # Create clusters
    cls_tbl = (
        cls_tbl.join(new_clusters, on="cell_barcode", how="left")
        .with_columns(
            pls.when(pls.col("mc_regroup_corr").is_null())
            .then(pls.col("mc_raw_cc_cb"))
            .otherwise(pls.col("mc_regroup_cc_cb"))
            .alias("mc_regroup_cc_cb")
        ))
    new_b2m_dict = dict(cls_tbl.select("mc_raw_cc_cb", "mc_raw_barcode").rows())
    regrouped_barcodes = [new_b2m_dict[x] for x in cls_tbl["mc_regroup_cc_cb"]]
    cls_tbl = cls_tbl.with_columns(pls.Series(regrouped_barcodes).alias("mc_regroup_barcode"))

    return cls_tbl


def resample_mc_expression(
    adata: adt.AnnData, mc_adata: adt.AnnData, logman: LogManager = LogManager("Resampling"), **kwargs
):
    '''Estimate the expression of metacell by resampling.'''
    perms = kwargs.get("perms", 1000)
    n_samples = kwargs.get("n_samples", 100)

    if adata.X is None:
        raise ValueError("anndata.X is None.")

    sampling_params = {}
    if isinstance(perms, int):
        logman.debug(f"Will using {perms} samples")
        sampling_params = {"n": perms}
    else:
        if not isinstance(perms, float):
            logman.warning("n_perms should be a real number. Using perms=2.0 by default")
            perms = 2.0

        if perms >= 1:
            logman.debug("Up-sampling.")
        else:
            logman.debug("Down-sampling.")

        sampling_params = {"fraction": perms}

    obs_names = adata.obs_names.tolist()
    var_names = adata.var_names.tolist()
    perm_tab =  pls.DataFrame(mc_adata.obs).with_columns(pls.col("cell_barcode").str.split("|")).explode("cell_barcode")
    obs_info_tab = (
        pls.DataFrame(adata.obs)
        .with_columns(pls.Series(adata.obs_names).alias("cell_barcode"))
        .join(perm_tab, on="cell_barcode")
        .sort("cell_barcode", pls.col("cell_barcode").cast(pls.Enum(obs_names)))
        .with_row_index()
        .select(["index", "cell_barcode", "mc_regroup_barcode"])
        .sort("index")
    )

    resampling_exp_tab_list, resampling_info_tab_list = [], []
    for per_name, per_group in obs_info_tab.group_by("mc_regroup_barcode"):
        per_name = "-".join(per_name)
        rsmp_vec_list, resampling_idx_list, rsmp_mmb_list = [], [], []
        for resample_index in range(1, n_samples + 1):
            rs_tab = per_group.sample(**sampling_params)

            selected_indices = rs_tab["index"].to_list()
            rsmp_vec_list.append(adata.X[selected_indices, :].toarray().mean(axis=0).T)
            resampling_idx_list.append(f"RI{resample_index:06d}")
            rsmp_mmb_list.append(rs_tab["cell_barcode"].str.join("|").to_list())

        base_tab = pls.DataFrame({"rs_barcode": resampling_idx_list, "mc_regroup_barcode": per_name})
        per_rsmp_exp_tab = pls.concat([pls.DataFrame(rsmp_vec_list, schema=var_names), base_tab], how="horizontal")
        per_rsmp_info_tab = pls.concat([pls.DataFrame({"cell_barcodes": rsmp_mmb_list}), base_tab], how="horizontal")

        resampling_exp_tab_list.append(per_rsmp_exp_tab)
        resampling_info_tab_list.append(per_rsmp_info_tab)

    exp_tab = pls.concat(resampling_exp_tab_list).with_columns(
        (pls.col("mc_regroup_barcode") + "-" + pls.col("rs_barcode")).alias("cell_barcode")
    )
    info_tab = pls.concat(resampling_info_tab_list)
    return exp_tab, info_tab


def create_metacells(
    adata: adt.AnnData, group_by: str | None | List[str] = None, min_cells_pmc: int = 10, min_cells_pg: int = 100,
    assay: str = "X_pca", first_n_dims: int = 30, do_resampling: bool = True, n_artificats: int = 50,
    project_back: bool = False, logman: LogManager = LogManager("CreateMetaCells"),
):
    '''Create metacells from anndata object.'''
    adata = adata.copy()

    if assay == "X":
        if first_n_dims is not None:
            logman.warning("The parameter 'first_n_dims' has no effect when the assay is 'X'.")

        adata.uns["mc_min_cells_per_gene"] = min_cells_pg # The parameter is only available when the assay is X.
        kept, _ = scp.pp.filter_genes(adata, min_cells=min_cells_pg, inplace=False)

        if adata.X is not None and isinstance(adata.X, csr_matrix):
            x_array = deepcopy(adata[:, kept].X.toarray())
        elif adata.X is not None and isinstance(adata.X, np.ndarray):
            x_array = deepcopy(adata[:, kept].X)
        else:
            raise ValueError("The input data is not a sparse matrix or a dense matrix.")

        assay_tbl = pls.DataFrame(x_array).rename(lambda x: adata.var_names[kept][x], axis=1)
    else:
        if assay not in ["X_pca", "X_umap"]:
            logman.warning("Unknown assay. Use 'X_pca' instead. Available choices (case sensitive): X, X_pca, X_umap.")
            assay = "X_pca"

        assay_key = "UMAP" if assay == "X_umap" else "PC"
        _, available_dims = adata.obsm[assay].shape

        assert isinstance(first_n_dims, int) and first_n_dims > 0, "`first_n_dims` must be a positive integer."
        if first_n_dims > available_dims:
            logman.warning(f"Only {available_dims} available PCs. Using {available_dims} PCs instead.")
            first_n_dims = available_dims

        selected_dims = range(first_n_dims)
        assay_tbl = pls.DataFrame(adata.obsm[assay][:, selected_dims]).rename(lambda x: x.replace("column", assay_key))

    cell_barcode = pls.DataFrame().with_columns(cell_barcode = pls.Series(adata.obs_names))
    mtdata_tbl = pls.from_pandas(adata.obs.loc[:, group_by])
    comb_tbl = pls.concat([cell_barcode, mtdata_tbl, assay_tbl], how="horizontal")

    tmp_list = []
    if group_by is None or (isinstance(group_by, list) and len(group_by) == 0):
        tmp_list.append(affinity_propagation(comb_tbl, min_cells_pmc=min_cells_pmc, group_by=group_by))
    else:
        for group_id, group_data in comb_tbl.group_by(group_by):
            logman.debug(f"Working on group {group_id} ...")
            group_id = "-".join(group_id)
            tmp_list.append(
                affinity_propagation(group_data, min_cells_pmc=min_cells_pmc, group_by=group_by, group_id=group_id)
            )
    cls_tbl = pls.concat(tmp_list)

    # MC features, gene expression
    obs_names, var_names = adata.obs_names.tolist(), adata.var_names.tolist()
    x_mat = pls.DataFrame(adata.X.toarray(), schema=var_names).with_columns(pls.Series(obs_names).alias("cell_barcode"))
    mc_v_tbl = pls.DataFrame({"gene_id": var_names}).to_pandas().set_index("gene_id", drop=False)

    mc_o_tbl = (
        cls_tbl[["mc_regroup_barcode", "cell_barcode"]]
        .group_by("mc_regroup_barcode")
        .agg(pls.col("cell_barcode").unique().str.join("|"))
        .to_pandas()
        .set_index("mc_regroup_barcode", drop=False)
        .sort_index()
    )

    # Average expression matrix and variance
    mc_obs_names = mc_o_tbl.index.to_list()
    exp_tbl = (
        cls_tbl[:, ["cell_barcode", "mc_regroup_barcode"]]
        .join(x_mat, on="cell_barcode")
        .drop("cell_barcode")
        .group_by("mc_regroup_barcode")
    )

    mc_x_tbl = csr_matrix(
        exp_tbl.agg(pls.all().mean())
        .sort("mc_regroup_barcode", pls.col("mc_regroup_barcode").cast(pls.Enum(mc_obs_names)))
        .to_pandas()
        .set_index("mc_regroup_barcode")
    )

    mc_x_var_tbl = csr_matrix(
        exp_tbl.agg(pls.all().var())
        .sort("mc_regroup_barcode", pls.col("mc_regroup_barcode").cast(pls.Enum(mc_obs_names)))
        .to_pandas()
        .set_index("mc_regroup_barcode")
    )

    mc_anndata = adt.AnnData(X=mc_x_tbl, obs=mc_o_tbl, var=mc_v_tbl)
    mc_anndata.layers["mc_exp_var"] = mc_x_var_tbl
    mc_anndata.uns["mc_source_assay"] = assay
    mc_anndata.uns["mc_min_cells_per_metacell"] = min_cells_pmc

    if project_back:
        sub_cls_tbl = cls_tbl[["mc_regroup_barcode", "cell_barcode"]]
        if "X_pca" in adata.obsm:
            pca_tbl = (
                pls.concat([sub_cls_tbl, pls.DataFrame(adata.obsm["X_pca"])], how="horizontal")
                .group_by("mc_regroup_barcode")
                .agg(pls.all().exclude("cell_barcode").mean())
                .sort("mc_regroup_barcode", pls.col("mc_regroup_barcode").cast(pls.Enum(mc_obs_names)))
            )
            mc_anndata.obsm["X_pca"] = pca_tbl.drop("mc_regroup_barcode").to_numpy()

        if "X_umap" in adata.obsm:
            umap_tbl = (
                pls.concat([sub_cls_tbl, pls.DataFrame(adata.obsm["X_umap"])], how="horizontal")
                .group_by("mc_regroup_barcode")
                .agg(pls.all().exclude("cell_barcode").mean())
                .sort("mc_regroup_barcode", pls.col("mc_regroup_barcode").cast(pls.Enum(mc_obs_names)))
            )
            mc_anndata.obsm["X_umap"] = umap_tbl.drop("mc_regroup_barcode").to_numpy()

    if do_resampling: # Expression matrix by resampling
        non_data_cols = ["rs_barcode", "mc_regroup_barcode", "cell_barcode"]
        rsmp_exp_tbl, _ = resample_mc_expression(adata, mc_anndata, perms=0.5, n_samples=n_artificats)
        var_dict = {"gene_id":rsmp_exp_tbl.drop(non_data_cols).columns}
        rsmp_v_tbl = pls.DataFrame(var_dict).to_pandas().set_index("gene_id", drop=False)
        rsmp_o_tbl = rsmp_exp_tbl.select(non_data_cols).to_pandas().set_index("cell_barcode", drop=False)
        rsmp_x_tbl = rsmp_exp_tbl.drop(["rs_barcode", "mc_regroup_barcode"]).to_pandas().set_index("cell_barcode")
        rsmp_anndata = adt.AnnData(X=rsmp_x_tbl, obs=rsmp_o_tbl, var=rsmp_v_tbl)
        rsmp_anndata.uns["resampling_perm"] = 0.5
    else:
        rsmp_anndata = None

    return mc_anndata, rsmp_anndata, cls_tbl


def load_dataset(
    in_file: str | Path, y_col: str | None = None, covar_cols: list[str] | None = None, barcode_col: str = "barcode",
    test_frac: float = 0.25
):
    """Load dataset from file."""
    logman = LogManager("Load Dataset")
    if covar_cols is None: covar_cols = []
    if isinstance(in_file, str): in_file = Path(in_file)

    if in_file.suffix == ".h5ad":
        adata = scp.read_h5ad(in_file)

        # Labels
        obs_tbl = pls.from_pandas(adata.obs.reset_index(drop=False, names=barcode_col).copy())
        y_vec = obs_tbl[y_col].to_list() if y_col in obs_tbl.columns else None

        # Expression matrix
        exp_mat = pls.from_pandas(adata.to_df().reset_index(drop=False, names=barcode_col).copy())
        x_mat = pls.DataFrame(exp_mat.select(pls.all()))

        if covar_cols: # Add covariables to expression matrix
            cov_tbl = pls.DataFrame(obs_tbl.select(covar_cols))
            x_mat = x_mat.join(cov_tbl, on=barcode_col)
    else:
        if in_file.suffix != ".csv":
            logman.warning(f"Unknown file format {in_file.suffix}. Using .csv as default.")

        x_mat = pls.read_csv(in_file, has_header=True)
        if y_col in x_mat.columns:
            y_vec = x_mat[y_col].to_list()
            x_mat = x_mat.drop(y_col)
        else:
            y_vec = None

    x_mat = x_mat.drop(barcode_col)
    if y_vec is None:
        x_mat_train, x_mat_test = train_test_split(x_mat, test_size=test_frac)
        y_vec_train, y_vec_test = None, None
    else:
        x_mat_train, y_vec_train, x_mat_test, y_vec_test = train_test_split(x_mat, y_vec, test_size=test_frac)
        x_mat_test = x_mat_test

    return x_mat_train, y_vec_train, x_mat_test, y_vec_test


def train_model(
    X, y = None, n_epochs: int = 100, batch_size: int = 32, learning_rate: float = 1e-3, disable_gpu: bool = False
):
    """Train a VAE model using given data."""
    device = "cuda" if torch.cuda.is_available() and not disable_gpu else "cpu"
    n_classes = len(y) if y is not None else None

    preproc_steps = [("feature_selection", FeatureSelector()), ("aligner", DatasetAligner())]
    preproc_ppl = Pipeline(steps=preproc_steps)
    preproc_ppl.fit(X, y)
    x_mat = preproc_ppl.fit_transform(X)
    _, in_dims = x_mat.shape

    vae_machine = VAEMachine(
        in_dims=in_dims, n_classes=n_classes, epochs=n_epochs, batch_size=batch_size, learning_rate=learning_rate,
        device=device
    )
    vae_machine.fit(x_mat, y)

    return preproc_ppl, vae_machine


def evaluate_model(vae_model: VAEMachine, preproc: Pipeline, X, y_true=None):
    """Evaluate the model using confusion matrix."""
    x_mat = torch.Tensor(preproc.transform(X)) # Preprocess the given dataset.

    x_recon = vae_model.to("cpu").reconstruct(x_mat).detach() # Reconstruct the dataset by the model.
    recon_loss = fun.mse_loss(x_recon, x_mat).detach().item() # Estimate the reconstruction loss.

    precision = accuracy = None
    if y_true:
        y_pred = vae_model.to("cpu").predict(X)
        precision = precision_score(y_true, y_pred, average="macro")
        accuracy = accuracy_score(y_true, y_pred)
        recall = recall_score(y_true, y_pred, average="macro")

    return {"Reconstruction loss": recon_loss, "Precision": precision, "Accuracy": accuracy, "Recall": recall}


def predict_samples(X, model):
    pass


def prepare_working_folder(path: str | Path, force: bool = False):
    if isinstance(path, str): path = Path(path)
    if path.exists():
        if force:
            shutil.rmtree(path)
        else:
            raise FileExistsError(f"Found {path}. Using force=True to remove the old and create a new one.")
    path.mkdir(parents=True)


def plot_qc_matrics(adata, out_dir, **kwargs):
    '''Plot quality control metrics.'''

    fig = plt.figure(layout=None, figsize=(10, 9), tight_layout=True)
    gs = fig.add_gridspec(nrows=2, ncols=4)

    axe1 = fig.add_subplot(gs[0, 0])
    scp.pl.violin(adata, keys=["n_genes"], ax=axe1)

    axe2 = fig.add_subplot(gs[0, 1])
    scp.pl.violin(adata, keys=["total_counts"], ax=axe2)

    axe3 = fig.add_subplot(gs[0, 2])
    scp.pl.violin(adata, keys=["pct_counts_mt"], ax=axe3)

    axe4 = fig.add_subplot(gs[0, 3])
    scp.pl.violin(adata, keys=["doublet_score"], ax=axe4)

    axe5 = fig.add_subplot(gs[1, :2])
    scp.pl.pca(adata, color="louvain", ax=axe5)

    axe6 = fig.add_subplot(gs[1, 2:])
    scp.pl.umap(adata, color="louvain", ax=axe6)

    fig.savefig(out_dir / "qc_metrics.pdf")
    fig.clear()
    plt.close()


@click.group()
@click.version_option(version="0.1.0")
@click.option("-v", "--verbose", default=0, count=True, help="Verbose mode.")
def app(verbose):
    logging.basicConfig(level=logging.INFO if verbose >= 2 else logging.WARNING)
    logman = LogManager("Main")
    logman.debug("Test")


@app.command()
@click.argument("input", metavar="DIR|AnnData")
@click.option("-b", "--batch-key", metavar="STR", show_default=True, help="Batch key.")
@click.option("-a", "--resolution", default=0.1, type=float, metavar="FLOAT", show_default=True, help="Resolution for clustering.")
@click.option("-p", "--n-comps", default=50, type=int, metavar="INT", show_default=True, help="Number of principal components to use.")
@click.option("-c", "--min-cells", default=3, type=int, metavar="INT", show_default=True, help="A gene to be kept if it is detected at minimum number of cells.")
@click.option("-g", "--min-genes", default=200, type=int, metavar="INT", show_default=True, help="A cell to be kept if minimum number of genes are detected.")
@click.option("-m", "--max-mito-pct", default=10, type=float, metavar="FLOAT", show_default=True, help="Maximum percentage of mitochondrial reads.")
@click.option("-d", "--max-dbl-score", default=0.2, type=float, metavar="FLOAT", show_default=True, help="Doublet score threshold.")
@click.option("-F", "--force", default=False, is_flag=True, show_default=True, help="!!!BE CAREFUL!!!. Removing existing output directory and create new one.")
@click.option("-D", "--debug", default=False, is_flag=True, show_default=True, help="Debug mode.")
# @click.option("-r", "--required-var", multiple=True, metavar="STR", show_default=True, help="Required metadata.")
# @click.option("--use-highly-variable-genes", default=False, is_flag=True, show_default=True, help="Use highly variable genes.")
# @click.option("--n-highly-variable-genes", default=2000, type=int, metavar="INT", show_default=True, help="Number of highly variable genes to use.")
# @click.option("--template", metavar="DIR", show_default=True, help="Template directory.")
@click.option("-o", "--out-dir", default="Preprocess", metavar="DIR", show_default=True, help="Output directory.")
def preproc(
    input, batch_key, resolution, n_comps, min_cells, min_genes, max_mito_pct, max_dbl_score, force, debug, out_dir
):
    '''Preprocess input data.'''
    logman = LogManager("Preprocessing", level=logging.DEBUG if debug else logging.INFO)

    out_dir = Path(out_dir)
    prepare_working_folder(out_dir, force)

    input = Path(input)
    if input.is_dir():
        logman.debug(f"Reading 10x data from {input} ...")
        anndata = scp.read_10x_mtx(input)
    else:
        logman.debug(f"Reading AnnData from {input} ...")
        anndata = scp.read_h5ad(input)

    logman.debug("Preprocessing data ...")
    anndata = preprocess(anndata, batch_key, min_cells, min_genes, max_mito_pct, max_dbl_score, n_comps, resolution)

    logman.debug("Creating overview plots ...")
    plot_qc_matrics(anndata, out_dir)

    logman.debug("Writing processed AnnData to disk ...")
    scp.write(out_dir / "preprocessed.h5ad", anndata)

    logman.debug("Writing parameters into disk ...")
    with open(out_dir / "parameters.json", "w") as fhandle:
        json.dump({
            "batch_key": batch_key, "resolution": resolution, "n_comps": n_comps, "min_cells": min_cells,
            "min_genes": min_genes, "max_mito_pct": max_mito_pct, "max_dbl_score": max_dbl_score
        }, fhandle, indent=4)


@app.command()
@click.argument("in_file", metavar="ANNDATA")
@click.option("-a", "--assay", default="X_pca", metavar="STR", type=click.Choice(["X", "X_pca", "X_umap"]), show_default=True, help="PCA assay to use.")
@click.option("-n", "--n-components", default=30, metavar="INT", type=int, show_default=True, help="Number of components to use if the assay is 'X_pca' or 'X_umap'.")
@click.option("-g", "--group-by", multiple=True, metavar="STR", show_default=True, help="Group cells by metadata.")
@click.option("-f", "--selected-features", default=None, metavar="FILE", show_default=True, help="Selected features file.")
@click.option("-m", "--min-cells-pmc", default=5, metavar="INT", type=int, show_default=True, help="Minimum number of cells per metacell.")
@click.option("-G", "--min-cells-pg", default=100, metavar="INT", type=int, show_default=True, help="Minimum number of cells per gene.")
@click.option("-P", "--project-back", default=False, is_flag=True, show_default=True, help="Project meta-cells or resampled-cells back to original space.")
@click.option("-R", "--do-resampling", default=False, is_flag=True, show_default=True, help="Create fake cells by resampling.")
@click.option("-N", "--n-artifacts", default=50, metavar="INT", type=int, show_default=True, help="Number of fake cells to create.")
@click.option("-D", "--debug", default=False, is_flag=True, show_default=True, help="Debug mode.")
@click.option("-F", "--force", default=False, is_flag=True, show_default=True, help="!!!BE CAREFUL!!!. Removing existing output directory and create new one.")
@click.option("-o", "--out-dir", default="Metacell", metavar="DIR", type=str, show_default=True, help="Output directory.")
def metacell(
    in_file: str, assay: str, n_components: int, group_by: str | list[str], selected_features: str, min_cells_pmc: int,
    min_cells_pg: int, project_back: bool, do_resampling: bool, n_artifacts: int, debug: bool, force: bool,
    out_dir: str | Path
):
    '''Create meta-cells from single-cell RNA-seq data.'''
    logman = LogManager("Metacell", level=logging.DEBUG if debug else logging.INFO)

    out_dir = Path(out_dir)
    prepare_working_folder(out_dir, force=force)
    
    # Loading data and obtain its subset if selected_features is available.
    logman.debug("Loading data ...")
    raw_adata = scp.read_h5ad(in_file)
    if selected_features:
        logman.debug("To use selected features ...")
        with open(selected_features, "r") as fhandle:
            target_features = fhandle.read().splitlines()
        example_features = ",".join(target_features[:5])
        logman.debug(f"The first 5 features are: {example_features}")
        raw_adata = raw_adata[:, target_features]

    # TODO, project resampled cells back to original PCA and UMAP
    logman.debug("Creating metacells ...")
    if isinstance(group_by, str):
        group_by = [group_by]
    elif isinstance(group_by, tuple):
        group_by = list(group_by)
    else:
        raise ValueError("`-g/--group-by` must be str, tuple, or list.")

    mc_adata, rsmp_adata, cls_tbl = create_metacells(
        raw_adata, group_by, min_cells_pmc, min_cells_pg, assay, n_components, do_resampling, n_artifacts, project_back
    )

    # Save results to disk
    logman.debug("Writing results to disk ...")
    if isinstance(mc_adata, adt.AnnData):
        logman.debug("Writing metacells to disk ...")
        mc_adata.write_h5ad(out_dir / "meta_cells_by_clustering.h5ad")
    else:
        logman.debug("The metacells are not in AnnData format. Skipping ...")

    if isinstance(rsmp_adata, adt.AnnData):
        logman.debug("Writing resampled cells to disk ...")
        rsmp_adata.write_h5ad(out_dir / "artificial_cell_by_resampling.h5ad")
    else:
        logman.debug("The resampled cells are not in AnnData format. Skipping ...")

    if isinstance(cls_tbl, pls.DataFrame):
        logman.debug("Writing clustering results to disk ...")
        cls_tbl.write_csv(out_dir / "cluster_results.csv")
    else:
        logman.debug("The cluster results are not in polars.DataFrame format. Skipping ...")


@app.command()
@click.argument("in_file", metavar="ANNDATA|CSV")
# @click.argument("model", metavar="MODEL")
@click.option("-y", "--y-col", default="label", metavar="STR", show_default=True, help="Column name of the label.")
@click.option("-p", "--test-set-ratio", default=0.25, metavar="FLOAT", type=float, show_default=True, help="Train-test split ratio.")
@click.option("-e", "--n-epochs", default=100, metavar="INT", type=int, show_default=True, help="Number of epochs to train the model.")
@click.option("-b", "--batch-size", default=32, metavar="INT", type=int, show_default=True, help="Batch size.")
@click.option("-l", "--learning-rate", default=1e-3, metavar="FLOAT", type=float, show_default=True, help="Learning rate.")
@click.option("-G", "--disable-gpu", default=False, is_flag=True, show_default=True, help="Disable GPU.")
@click.option("-D", "--debug", default=False, is_flag=True, show_default=True, help="Debug mode.")
@click.option("-F", "--force", default=False, is_flag=True, show_default=True, help="!!!BE CAREFUL!!!. Removing existing output directory and create new one.")
@click.option("-o", "--out-dir", metavar="DIR", default="Train", show_default=True, help="Output directory.")
def train(
    in_file: str, test_set_ratio: float, y_col: str, n_epochs: int, batch_size: int, learning_rate: float,
    disable_gpu: bool, debug: bool, force: bool, out_dir: str | Path
):
    '''Train a model to predict embryo development stage from single-cell RNA-seq data.'''
    logman = LogManager("Train", level=logging.DEBUG if debug else logging.INFO)

    out_dir = Path(out_dir)
    prepare_working_folder(out_dir, force=force)
    
    logman.debug("Loading data ...")
    x_mat_train, y_vec_train, x_mat_test, y_vec_test = load_dataset(in_file, y_col, test_frac=test_set_ratio)

    logman.debug("Training model ...")
    preproc_ppl, vae_model = train_model(x_mat_train, y_vec_train, n_epochs, batch_size, learning_rate, disable_gpu)

    logman.debug("Saving model ...")
    torch.save(vae_model.model.state_dict(), out_dir / "vae_machine.pickle") # Save model to file: vae_machine.pickle

    logman.debug("Saving preprocessing pipeline ...")
    with open(out_dir / "preproc_pipeline.pickle", "wb") as f:
        pickle.dump(preproc_ppl, f)

    logman.debug("Evaluating model ...")
    test_matics = evaluate_model(vae_model, preproc_ppl, x_mat_test, y_vec_test)

    logman.debug("Saving model evaluation result to disk ...")
    with open(out_dir / "test_metrics.json", "w") as f:
        json.dump(test_matics, f, indent=4)


@app.command()
@click.argument("expr-mat", metavar="EXPRESSION_MATRIX")
@click.argument("model", metavar="MODEL")
@click.option("-o", "--out-dir", metavar="DIR", default="Predict", help="Final prediction to be saved.")
def predict(expr_mat, model, out_dir):
    '''Predict embryo development stage from gene expression matrix using a model trained by this tool.'''
    vae_machine = VAEMachine()
    vae_machine.model.load_state_dict(torch.load(model))
    vae_machine.predict(expr_mat)


if __name__ == '__main__':
    app(max_content_width=120)

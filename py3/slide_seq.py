#!/usr/bin/env python3
# File: kidney_spatial.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jul 24, 2025
# Updated:

'''A script to process Slide-seq results.'''

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)
# warnings.filterwarnings("ignore", category=NumbaWarning)

import copy
import logging
import math
import itertools
from pathlib import Path

import click
import anndata as adt
import numpy as npy
import pandas as pds
import polars as pls
import scanpy as scp

import torch
import torch.nn.functional as func

from scipy.sparse import csr_matrix, csc_matrix
from plotnine import ggplot, aes, labs, theme, theme_classic, geom_tile, scale_fill_continuous


# Logging manager.
class LogManager(logging.Logger):
    def __init__(self, name, level=logging.INFO, logstream: bool = True,
                 logfile: str = ""):
        super(LogManager, self).__init__(name)

        fmt = logging.Formatter("{levelname: >8} | {asctime} | {name: >20} | {message} ...", "%Y%m%d,%H:%M:%S", "{")
        if logstream:
            self._add_handler(logging.StreamHandler(), level, fmt)

        if logfile:
            self._add_handler(logging.FileHandler(logfile), level, fmt)

    def _add_handler(self, hdl, lvl, fmt):
        hdl.setLevel(lvl)
        hdl.setFormatter(fmt)
        self.addHandler(hdl)


def load_adata(
    in_dir: str | Path,
    logman: LogManager = LogManager("LoadAnnData")
):
    """Create an AnnData object from 10x read count matrix."""
    if isinstance(in_dir, str): in_dir = Path(in_dir)
    try:
        # load data
        logman.info("Load data from disk")
        logman.info(f"Resources: {in_dir}")
        tenx_dir, spatial_info = in_dir / 'matrix', in_dir / 'spatial'
        pos_tab = pls.read_parquet(spatial_info / 'tissue_positions.parquet').filter(pls.col('in_tissue') == 1)
        per_adt = scp.read_10x_mtx(tenx_dir)

        # add information of observations
        per_adt.obs['barcode'] = per_adt.obs.index
        per_adt.obsm["spatial"] = (pls.from_pandas(per_adt.obs.reset_index(names="Barcode"))
                                   .join(pos_tab, on="Barcode", how="left")
                                   .select(["row_in_bin", "col_in_bin"])
                                   .to_numpy())
        per_adt.obs['x'] = per_adt.obsm.get('spatial')[:, 0]
        per_adt.obs['y'] = per_adt.obsm.get('spatial')[:, 1]

        per_adt.raw = per_adt # keep the raw results.
    except FileNotFoundError as e:
        logman.critical(f"Failed to load 10x results for {in_dir} due to: {e}")
        per_adt = None

    return per_adt


def load_labeled_adata(
    total_dir: str | Path,
    nascent_dir: str | Path | None = None,
    min_genes: int = 100,
    min_beads: int = 3,
    extra_meta: dict | None = None,
    qc_patterns: dict | None = {"Mt": "^mt-", "Rpl": "^Rpl"},
    logman: LogManager = LogManager("LoadLabeledAnnData")
):
    """Load Labeled AnnData."""
    if isinstance(total_dir, str): total_dir = Path(total_dir)
    try:
        # load data
        adata = load_adata(total_dir)
        if isinstance(extra_meta, dict) and adata is not None:
            for key, val in extra_meta.items():
                adata.obs[key] = val
        if adata is None: return None

        if adata.layers is None: adata.layers = {}
        if "total" not in adata.layers: adata.layers["total"] = adata.X

        # add nascent read counts
        if nascent_dir is not None:
            nascent_adata = load_adata(nascent_dir)
            if nascent_adata is not None:
                logman.info("Add nascent reads to the object")
                total_vars, total_obs = adata.var_names.to_list(), adata.obs_names.to_list()
                total_spatial_pos = pls.from_pandas(adata.obs).select("barcode", "x", "y")

                nascent_vars = nascent_adata.var_names.to_list()
                missing_cols = [x for x in total_vars if x not in nascent_vars]
                n_nascent_cells, n_missing_cols = nascent_adata.shape[0], len(missing_cols)
                missing_cols = pls.DataFrame(npy.zeros((n_nascent_cells, n_missing_cols)), missing_cols)

                count_matrix = pls.from_numpy(
                    nascent_adata.X.toarray(), schema=nascent_vars
                ).with_columns(
                    pls.Series("barcode", nascent_adata.obs_names.to_list()),
                    pls.Series("x", nascent_adata.obs.x),
                    pls.Series("y", nascent_adata.obs.y)
                ).with_columns(
                    missing_cols
                ).join(
                    total_spatial_pos, on=["barcode", "x", "y"], how="right"
                ).select(["barcode"] + total_vars).to_pandas().set_index("barcode").loc[total_obs, :]

                adata.layers["nascent"] = csc_matrix(count_matrix)

        # filter
        logman.info(f"Filter cells and genes. min_genes = {min_genes}, min_beads = {min_beads}")
        scp.pp.filter_cells(adata, min_genes=min_genes)
        scp.pp.filter_genes(adata, min_cells=min_beads)
        assert adata.shape[0] >= 10, ValueError(
            f"Not enough beads left for analysis: observed {adata.shape[0]}, but expect at least 10"
        )

        # add QC metrics
        if isinstance(qc_patterns, dict) and len(qc_patterns) > 0:
            logman.info("Calculate QC metrics. " + ", ".join(qc_patterns.keys()))
            for key, val in qc_patterns.items():
                adata.var[key] = adata.var.index.str.contains(val)
            scp.pp.calculate_qc_metrics(adata, qc_vars=list(qc_patterns.keys()), inplace=True)

        return adata
    except FileNotFoundError as e:
        logman.critical(f"Failed to load 10x results for {total_dir}, {nascent_dir}, due to: {e}")

    return None


def remove_outliers(col: pds.Series, n_times: float = 2, method: str = "mean"):
    """Removing outliers using IQR for pandas columns"""
    if method == "quantile":
        q1, q3 = col.quantile(0.25), col.quantile(0.75)
        lower_bound, upper_bound = q1 - n_times * (q3 - q1), q3 + n_times * (q3 - q1)
    else:
        if method != "mean": print("Using default method: 'mean'")
        m, v = col.mean(), col.std()
        lower_bound, upper_bound = m - n_times * v, m + n_times * v
    return col.where((lower_bound <= col) & (col <= upper_bound))


def convolute_count_matrix(
    adata,
    kernel_size: int = 4,
    target_dim: int = 0,
    min_counts: int = 25,
    drop_singletons: bool = False,
    singleton_max_n_neigh: int = 2,
    keep_meta_info: str | list[str] = ["batch_id", "sample_id"],
    logman: LogManager = LogManager("ConvoluteCountMatrix")
) -> adt.AnnData:
    """Convolute the read counts matrix to improve the expression profile."""
    if target_dim % kernel_size != 0:
        raise ValueError(f"target_dim () should be divided by kernel_size ({kernel_size}) with remainder 0.")

    # Ensemble X array
    pos_list = adata.obs.loc[:, ['x', 'y']].apply(remove_outliers) # Removing outliers from the original data.
    (x_min, y_min), (x_max, y_max) = pos_list.loc[:, ['x', 'y']].min().to_list(), pos_list.loc[:, ['x', 'y']].max().to_list()
    width_raw = math.ceil((max(x_max - x_min, y_max - y_min) / kernel_size)) + 1

    logman.info("Remove outliers")
    layer_data, singleton_pos = {}, []
    for key in ["X", *list(adata.layers.keys())]:
        matrix = adata.X if key == "X" else adata.layers[key]
        if isinstance(matrix, (csc_matrix, csr_matrix)): matrix = matrix.toarray()
        assert isinstance(matrix, npy.ndarray), ValueError("matrix should be a numpy.ndarray.")

        logman.info(f"Convolute '{key}'")
        conv_fill_list, conv_list = [], []
        for per_values in matrix.T:
            # Transfer the vector of gene expression into a matrix.
            mat_cmp = torch.zeros((width_raw, width_raw)).unsqueeze(0)
            for val, (_, row) in zip(per_values, pos_list.iterrows()):
                x_c, y_c = row.to_list()
                x_c, y_c = (x_c - x_min) / kernel_size, (y_c - y_min) / kernel_size
                if npy.isnan(x_c) or npy.isnan(y_c): continue
                mat_cmp[0, int(x_c), int(y_c)] = torch.tensor(val)

            # Convolve
            per_mat_conv = func.avg_pool2d(mat_cmp, kernel_size, stride=kernel_size) * (kernel_size ** 2)
            if target_dim >= 5:
                per_mat_conv = func.adaptive_avg_pool2d(per_mat_conv, target_dim)
                per_mat_conv = (per_mat_conv / per_mat_conv.sum() * mat_cmp.sum()).ceil()
                width_conv = target_dim
            else:
                _, width_conv, _ = per_mat_conv.shape
            conv_list.append(copy.copy(per_mat_conv.reshape(1, width_conv**2)))

            # Fill empty holes
            mask = (per_mat_conv == 0).float()
            smoothed = func.avg_pool2d(per_mat_conv, 3, stride=1, padding=1)
            per_mat_conv = per_mat_conv * (1 - mask) + smoothed * mask
            conv_fill_list.append(per_mat_conv.reshape((1, width_conv**2))) # Reshape it into matrix

        conv_mat, conv_fill_mat = torch.stack(conv_list, dim=2).squeeze(), torch.stack(conv_fill_list, dim=2).squeeze()

        # Removing singletons
        if key == "X": # For now only using X to identify singletons
            logman.info(f"Identify singletons using '{key}' matrix")
            width_final = int(math.sqrt(conv_fill_mat.shape[0]))
            x_conv_sum = func.pad(conv_fill_mat.sum(1).reshape(width_final, width_final), (1, 1, 1, 1))
            for i, j in itertools.product(range(width_final), range(width_final)):
                center = x_conv_sum[i, j]
                if center == 0: continue
                window = x_conv_sum[i-1:i+2, j-1:j+2] # Extract 3x3 window
                is_singleton = 0 < (window > 0).sum().item() <= singleton_max_n_neigh # Check its neighbors
                if is_singleton: singleton_pos.append((i - 1, j - 1))

        logman.info(f"Update convolved matrix for '{key}'")
        layer_data.update({f"{key}_convolved": csc_matrix(copy.copy(conv_mat.numpy().astype("float32"))),
                           f"{key}_convolved_filled": csc_matrix(copy.copy(conv_fill_mat.numpy().astype("float32")))})

    # Meta information
    logman.info("Add meta information")
    width_final, _ = layer_data["X_convolved"].shape
    new_xy = list(itertools.product(range(int(math.sqrt(width_final))), repeat=2))
    barcode_len = math.floor(math.log2(width_final) / 2) + 1
    obs_info_dict = {
        "x": [i[0] for i in new_xy], "y": [i[1] for i in new_xy], "conv_kernel_size": kernel_size,
        "is_singleton": [(x, y) in singleton_pos for x, y in new_xy],
        "barcodes": list(set([''.join(x) for x in itertools.product("ATCG", repeat=barcode_len)]))[:width_final]}
    for per_var in keep_meta_info:
        required_meta_info = adata.obs.loc[:, per_var].drop_duplicates().to_list()
        if len(required_meta_info) == 1:
            obs_info_dict[per_var] = required_meta_info[0]

    logman.info("Create new AnnData from the processed results")
    obs_df = pds.DataFrame(obs_info_dict).set_index("barcodes", drop=False)
    var_df = pds.DataFrame({"gene_ids": adata.var_names}).set_index("gene_ids", drop=False)
    new_adata = adt.AnnData(
        layer_data["X_convolved_filled"], obs=obs_df, var=var_df, layers=layer_data,
        obsm={"X_spatial": obs_df.loc[:, ["x", "y"]].to_numpy()}
    )
    new_adata.uns["raw"] = adata.copy()

    if drop_singletons:
        logman.info(f"Remove singletons")
        new_adata = new_adata[new_adata.obs.is_singleton == False, :]

    scp.pp.filter_cells(adata, min_counts=min_counts)

    return new_adata


def get_spatial_field(
    adata,
    k_size: int = 0,
    as_type: str = "mat"
):
    """Create a matrix to contain beads in the space."""
    max_pos = adata.obs.loc[:, ['x', 'y']].max().max()
    min_pos = adata.obs.loc[:, ['x', 'y']].min().min()
    if k_size > 0:
        width_cmp = math.ceil((max_pos - min_pos + 1) / k_size)
    else:
        width_cmp = max_pos - min_pos + 1

    if as_type == "vec":
        spatial_matrix = npy.zeros(width_cmp ** 2)
    else:
        if as_type != "mat": print("[W]: Return mat by default.")
        spatial_matrix = npy.zeros([width_cmp, width_cmp])

    return spatial_matrix, width_cmp, min_pos, max_pos


def plot_features(
    adata, feature: str = "Total counts", save_to: str | Path = "total_counts.png", figsize=(7, 7), cmap_name="plasma"
):
    """Create an expression matrix in the x-y space of the slide-seq chip."""
    img, _, min_pos, _ = get_spatial_field(adata)

    if feature == "Total counts":
        count_table = adata.X.toarray().sum(1)
    elif feature == "Genes by counts":
        count_table = (adata.X.toarray() > 0).sum(1)
    elif feature in adata.var_names:
        count_table = adata[:, feature].X.toarray()
    elif feature in adata.layers.keys():
        count_table = adata.layers.get(feature).sum(1)
    else:
        raise ValueError(f"Not found {feature}")

    plot_tab = adata.obs.loc[:, ['x', 'y']].assign(expr=count_table)
    p = ggplot(plot_tab, aes("x", "y", fill="expr")) \
        + geom_tile() \
        + scale_fill_continuous(cmap_name=cmap_name, name=feature) \
        + labs(x = "spatial_1", y = "spatial_2") \
        + theme_classic() \
        + theme(legend_position="top")
    fig_width, fig_height = figsize
    p.save(save_to, width=fig_width, height=fig_height)


def create_meta_info(meta_info):
    """Create meta information from given string."""
    new_meta_info = {}
    for x in meta_info:
        key, *val = x.split(":")
        if len(val) < 1:
            raise ValueError(f"Failed to parse {x}")
        else:
            val = ";".join(val)
        new_meta_info[key] = val

    return new_meta_info


@click.command()
@click.argument("in_dir", type=click.Path(exists=True, dir_okay=True))
@click.argument("batch_id", type=str)
@click.argument("sample_id", type=str)
@click.option("-k", "--kernel-size", type=int, default=4, show_default=True, help="The kernel size to convolve the pixel matrix.")
@click.option("-m", "--min-counts", type=int, default=25, show_default=True, help="Minmal number of total read counts per bead.")
@click.option("-s", "--target-dim", type=int, default=0, show_default=True, help="The size of final output pixel matrix.")
@click.option("-D", "--drop-singletons", is_flag=True, help="Drop singletons.")
@click.option("-M", "--meta-info", type=str, multiple=True, show_default=True, help="Meta information for current sample. E.g., Layer:1")
@click.option("-o", "--out-dir", type=str, default=None, show_default=True, help="The prefix of output files.")
def main(in_dir: Path | str, batch_id: str, sample_id: str, **kwargs):
    control_samples = {
        "20250226_decoded_embryo": [ "OST110029", "OST110090", ],
        "20250403_decoded_embryo": [ "STO110111" ],
        "20250314_decoded_embryo": [ "OST110090" ],
    }

    all_samples = {
        # TODO update the T->C results using quality corrected pipeline.
        "20250314_decoded_embryo": {
            "successed_samples": [ "OST110091", "OST110093", ],
            "failed_samples": [ "OST110092", "OST110094", ],
        },
        "20250403_decoded_embryo": {
            "successed_samples": [ "STO110112", ],
            "failed_samples": [ "STO110113", ]
        },
        "20250423_decoded_embryo": {
            "successed_samples": [ "STO110107", "STO110108", "STO110109", "STO110114", "STO110115", ],
            "failed_samples": [ "STO110110", "STO110116", "STO110117", ]
        },
        "20250522_decoded_embryo" : {
            "successed_samples": [ "OST110079", "OST110080", "STO110135", "STO110136", "STO110137", ],
            "failed_samples": [ "STO110118", "OST110082", "OST110081", ],
        },
        "20250604_decoded_embryo": {
            "successed_samples": [ "STO110150", "STO110149", "STO110148", "STO110147", "STO110146", "STO110143", ],
            "failed_samples": [ "STO110145", "STO110144", ]
        },
        "20250728_decoded_embryo": {
            "successed_samples": [ "STO110151", "STO110152", "STO110153", "STO110154" ],
            "failed_samples": []
        }
    }

    # PROJECT_DIR = Path('~/Documents/projects/wp_vasaseq').expanduser()
    # base_in_dir = PROJECT_DIR / 'outputs/analysis/preprocessing/slide_seq_decoded'

    if batch_id in all_samples:
        if sample_id not in all_samples[batch_id]["successed_samples"]:
            raise ValueError(f"{sample_id} is not valid.")
    else:
        raise ValueError(f"{batch_id} is not valid.")

    in_dir = Path(in_dir).expanduser()

    kernel_size = kwargs["kernel_size"]
    target_dim = kwargs["target_dim"]
    min_counts = kwargs["min_counts"]
    drop_singletons = kwargs["drop_singletons"]

    meta_info = kwargs["meta_info"]
    meta_info = create_meta_info(meta_info)
    meta_info.update({"batch_id": batch_id, "sample_id": sample_id})

    out_dir = kwargs["out_dir"]
    if kwargs.get("out_dir") is None:
        out_dir = f"{batch_id}.{sample_id}"
    out_dir = Path(out_dir)

    if not out_dir.exists(): out_dir.mkdir(parents=True)

    try:
        total_dir = in_dir / batch_id / sample_id / f"06.binSegment/square_bin/{sample_id}_Raw"
        nascent_dir = in_dir / batch_id / sample_id / f"06.binSegment_nascent/square_bin/{sample_id}_Raw"
        adata = load_labeled_adata(total_dir, nascent_dir, extra_meta=meta_info)
        adata_conv = convolute_count_matrix(adata, kernel_size, target_dim, min_counts, drop_singletons)

        out_file = out_dir / f"k_size_{kernel_size}.min_counts_{min_counts}.target_dim_{target_dim}.total_counts.png"
        plot_features(adata_conv, save_to=out_file)

        out_file = out_dir / f"k_size_{kernel_size}.min_counts_{min_counts}.target_dim_{target_dim}.genes_by_counts.png"
        plot_features(adata_conv, feature="Genes by counts", save_to=out_file)

        out_file = out_dir / f"k_size_{kernel_size}.min_counts_{min_counts}.target_dim_{target_dim}.h5ad"
        adata_conv.write_h5ad(out_file)
    except Exception as e:
        print(batch_id, sample_id, e)
        return 1


if __name__ == "__main__":
    main()

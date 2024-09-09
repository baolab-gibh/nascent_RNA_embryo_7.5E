#!/usr/bin/env python3
# File: infer_embryo_stage.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jul 27, 2024
# Updated: Jul 27, 2024

# Built-in
import os
import logging
from argparse import ArgumentParser

# Third-party
import numpy as np
import joblib as jb
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.base import BaseEstimator
from sklearn.metrics import recall_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.pipeline import Pipeline
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectPercentile

import scanpy

import torch
import torch.nn as nn
import torch.nn.functional as func
import torch.distributions as dist
import torch.optim as optim
from torch.utils.data.dataloader import Dataset
from torch.utils.data.dataloader import DataLoader
from torchviz import make_dot
from torchsummary import summary as tr_summary
from torchsampler import ImbalancedDatasetSampler


# Log manager
class LogManager(logging.Logger):
    def __init__(self, name, level=logging.INFO, logstream: bool = True, logfile: str = ""):
        super(LogManager, self).__init__(name)
        fmt = logging.Formatter("{levelname: >8}|{asctime}|{name: >8}| {message}", style="{", datefmt="%Y%m%d,%H%M%S")
        if logstream:
            self._add_handler(logging.StreamHandler(), level, fmt)

        if logfile:
            self._add_handler(logging.FileHandler(logfile), level, fmt)

    def _add_handler(self, hdl, lvl, fmt):
        hdl.setLevel(lvl)
        hdl.setFormatter(fmt)
        self.addHandler(hdl)


class DatasetAligner(BaseEstimator):
    def __init__(self):
        self.required_cols_ = []

    def fit(self, X: pd.DataFrame, y=None):
        self.required_cols_ = X.columns.to_list()
        return self

    def transform(self, X: pd.DataFrame, y=None):
        if not isinstance(X, pd.DataFrame):
            raise ValueError("The X should be a pandas.DataFrame!")

        available_cols = X.columns.to_list()
        extra_cols = [x for x in available_cols if x not in self.required_cols_]
        miss_cols = [x for x in self.required_cols_ if x not in available_cols]
        if len(extra_cols) != 0 or len(miss_cols) != 0:
            raise ValueError("Input DataFrame columns don't match ones on which model was trained."
                             "The --meta-cols columns are removed before this step."
                             f"Missing columns: {miss_cols}"
                             f"Extra columns: {extra_cols}")
        
        return X.loc[:, self.required_cols_]

    def fit_transform(self, X: pd.DataFrame, y=None):
        self.fit(X)
        return self.transform(X)


class CategoryEncoder(BaseEstimator):
    """A category variable transformer applied for whole DataFrame."""
    def __init__(self):
        self.categories_ = {}

    def _transformer(self, series: pd.Series, inverse=False):
        trans_dict = self.categories_.get(series.name, {})
        if not inverse:
            trans_dict = {v: k for k, v in trans_dict.items()}
        if trans_dict:
            return series.apply(lambda x: trans_dict[x])

        return series

    def fit(self, X: pd.DataFrame, y=None):
        _is_cat = X.dtypes == "object"
        self.categories_ = (X.loc[:, _is_cat].reset_index(drop=True).apply(lambda x: x.drop_duplicates().reset_index(drop=True)).to_dict())

        return self

    def transform(self, X: pd.DataFrame, y=None, copy=None):
        return X.copy().apply(self._transformer) if copy else X.apply(self._transformer)

    def inverse_transform(self, X, y=None, copy=None):
        return X.copy().apply(self._transformer, inverse=True) if copy else X.apply(self._transformer, inverse=True)

    def fit_transform(self, X, y=None):
        self.fit(X)
        return self.transform(X)


class ExpDataSet(Dataset):
    """Load expression data."""
    def __init__(self, X, y=None):
        super(ExpDataSet, self).__init__()
        self.x_mat = X
        self.y_vec = y

    def __len__(self):
        return self.x_mat.shape[0]

    def __getitem__(self, idx: int):
        x_vec = torch.Tensor(self.x_mat[idx, :])

        if self.y_vec is None:
            y_vec = torch.tensor(torch.nan)
        else:
            y_vec = torch.tensor(self.y_vec[idx], dtype=torch.long)

        return x_vec, y_vec

    def __next__(self):
        for x in range(len(self)):
            yield self[x]

    def get_labels(self):
        return self.y_vec


# Metacell


# VAE
class Encoder(nn.Module):
    """Encoder."""
    def __init__(self, latent_dims, in_size=256):
        super(Encoder, self).__init__()
        self.ll_1 = nn.Linear(in_size, 128)
        self.ll_2 = nn.Linear(128, 64)
        self.ll_3 = nn.Linear(64, 32)
        self.ll_4 = nn.Linear(32, 16)
        self.ll_5 = nn.Linear(16, latent_dims)
        self.dp = nn.Dropout(0.5)

    def forward(self, x):
        x = func.relu(self.ll_1(x))
        x = func.relu(self.ll_2(x))
        x = func.relu(self.ll_3(x))
        x = func.relu(self.ll_4(x))
        x = self.dp(x)
        return self.ll_5(x)


class LatentZ(nn.Module):
    """Latent layer."""
    def __init__(self, latent_dims, n_labs):
        super(LatentZ, self).__init__()
        self.ll_mu = nn.Linear(latent_dims, latent_dims) # LD mean
        self.ll_lnvar = nn.Linear(latent_dims, latent_dims) # LD log(var)
        self.ll_pred = nn.Linear(latent_dims, n_labs) # Prediction

    @property
    def MINF32(self):
        return 1.175494e-34

    def forward(self, x):
        mu_exp = self.ll_mu(x) # Learn the expectation (E)
        lnvar = func.softplus(self.ll_lnvar(x))

        std = torch.exp(0.5 * lnvar) # Learn the stdev (SD)
        latent_dist = dist.Normal(mu_exp, std)
        mu_hat = latent_dist.rsample()

        y_pred = torch.softmax(self.ll_pred(mu_hat), 1) # Prediction using the resampled mu
        return mu_hat, mu_exp, y_pred


class Decoder(nn.Module):
    """Decoder."""
    def __init__(self, latent_dims, out_size=256):
        super(Decoder, self).__init__()
        self.ll_1 = nn.Linear(latent_dims, 16)
        self.ll_2 = nn.Linear(16, 32)
        self.ll_3 = nn.Linear(32, 64)
        self.ll_4 = nn.Linear(64, 128)
        self.ll_5 = nn.Linear(128, out_size)

    def forward(self, x):
        x = func.relu(self.ll_1(x))
        x = func.relu(self.ll_2(x))
        x = func.relu(self.ll_3(x))
        x = func.relu(self.ll_4(x))
        return self.ll_5(x)


class VAE(nn.Module):
    """Variational AutoEncoder."""
    def __init__(self, latent_dims, in_size, n_labs: int = 2):
        super(VAE, self).__init__()
        self.encoder = Encoder(latent_dims, in_size)
        self.latentz = LatentZ(latent_dims, n_labs)
        self.decoder = Decoder(latent_dims, in_size)

    def forward(self, x):
        x_enc = self.encoder(x)
        mu_hat, mu_exp, y_pred = self.latentz(x_enc)
        x_dec = self.decoder(mu_hat)
        return x_enc, mu_hat, mu_exp, y_pred, x_dec


class VAEMachine(BaseEstimator):
    """A VAE designed for sklearn.Pipeline."""
    def __init__(self, learning_rate=1e-6, latent_dims=4, n_epoch=25, bsize=64,
                 disable_cuda=False, log_per_n_epoch=20, label_order=[]):
        self.learning_rate = learning_rate
        self.latent_dims = latent_dims
        self.n_epoch = n_epoch
        self.bsize = bsize
        self.model = None

        self._loss_ce = nn.CrossEntropyLoss(reduction="mean")
        self._loss_kld = nn.KLDivLoss(reduction="batchmean")
        self._loss_mse = nn.MSELoss(reduction="mean")
        self._eval_vec = {"Loss": {"Reconstruction": [], "Prediction": []}, "Train": {"Accuracy": [], "Recall": [], "Precision": []}}

        self._disable_cuda = disable_cuda
        self._log_per_n_epoch = log_per_n_epoch
        self._label_encoder = LabelEncoder()
        self._label_order = label_order

    def _calc_rec_loss(self, x_dec, x_mat, mu_hat, mu_exp):
        return (self._loss_mse(x_dec, x_mat) + self._loss_kld(torch.log_softmax(mu_hat, dim=0), torch.softmax(mu_exp, dim=0)))

    def _calc_pre_loss(self, pred, true_lab):
        return self._loss_ce(pred, true_lab)

    def _calc_eval_par(self, y_pred, y_true, which="Train"):
        y_pred_lab = np.int8(y_pred.to("cpu").max(1)[1].detach().numpy())
        y_true_lab = np.int8(y_true.to("cpu").detach().numpy())
        reca = recall_score(y_true_lab, y_pred_lab, average="macro") # tp/p
        self._eval_vec[which]["Recall"].append(reca)
        accu = accuracy_score(y_true_lab, y_pred_lab) # (tp+tn)/(p+n)
        self._eval_vec[which]["Accuracy"].append(accu)
        prec = precision_score(y_true_lab, y_pred_lab, average="macro", zero_division=0) # tp/pp
        self._eval_vec[which]["Precision"].append(prec)

        return reca, accu, prec

    def _running_log(self, idx, **items):
        item_key = ["Epoch", "Reconstruction", "Prediction", "Accuracy", "Precision", "Recall", "Partition"]
        head_line, record_line = [], []
        for per_key in item_key:
            per_val = items.get(per_key)
            if per_val is None: continue
            if idx == 0: head_line.append(f"{per_key:>12.12}")

            if per_key in ["Epoch"]: record_line.append(f"{per_val+1:12d}")
            elif per_key in ["Partition"]: record_line.append(f"{per_val:>12.12}")
            elif per_key in ["Reconstruction", "Prediction"]: record_line.append(f"{per_val:12.3f}")
            else: record_line.append(f"{per_val:12.10f}")
        if idx == 0: print(",".join(head_line) + ", Device")
        print(",".join(record_line) + ", " + self.usedev)

    @property
    def get_labels(self):
        return self._label_encoder.classes_

    @property
    def eval_per_epoch(self):
        return self._eval_vec

    @property
    def optimizer(self):
        if self.model is None:
            raise ValueError("Model isn't trained yet, please use fit() first.")
        return optim.Adam(self.model.parameters(), lr=self.learning_rate)

    @property
    def usedev(self): # Using GPU
        if torch.cuda.is_available() and not self._disable_cuda:
            return "cuda:" + str(torch.cuda.current_device())
        return "cpu"

    def fit(self, X, y):
        """Fit the model to given dataset."""
        #How to backward two loss terms between which weights are crossovered.
        #https://discuss.pytorch.org/t/how-to-combine-multiple-criterions-to-a-loss-function/348/50
        sort_by = None # Sort the labels and encode it into integers.
        if self._label_order and isinstance(self._label_order, list):
            sort_by = lambda x: self._label_order.index(x)
        avail_labs = y.drop_duplicates().sort_values(ignore_index=True, key=sort_by)
        self._label_encoder.fit(avail_labs)
        y = self._label_encoder.transform(y)

        trainset = ExpDataSet(X, y)
        trainset_loader = DataLoader(trainset, sampler=ImbalancedDatasetSampler(trainset), batch_size=self.bsize)
        self.model = VAE(self.latent_dims, X.shape[1], len(avail_labs)).to(self.usedev)
        for pe in range(self.n_epoch):
            y_pred_pe, y_true_pe = torch.Tensor().to(self.usedev), torch.Tensor().to(self.usedev)
            preloss_pe = torch.tensor(0, dtype=torch.float, device=self.usedev)
            recloss_pe = torch.tensor(0, dtype=torch.float, device=self.usedev)

            for x_mat, y_vec in trainset_loader:
                x_mat, y_vec = x_mat.to(self.usedev), y_vec.to(self.usedev)
                *_, mu_hat, mu_exp, y_pred, x_dec = self.model(x_mat)

                self.optimizer.zero_grad()
                loss_pre = self._calc_pre_loss(y_pred, y_vec)
                preloss_pe += loss_pre.item()
                loss_pre.backward(retain_graph=True)
                loss_rec = self._calc_rec_loss(x_dec, x_mat, mu_hat, mu_exp)
                recloss_pe += loss_rec.item()
                loss_rec.backward()
                self.optimizer.step()

                y_pred_pe = torch.concat((y_pred_pe, y_pred))
                y_true_pe = torch.concat((y_true_pe, y_vec))

            recloss = recloss_pe.to("cpu").detach().item()
            self._eval_vec["Loss"]["Reconstruction"].append(recloss)
            preloss = preloss_pe.to("cpu").detach().item()
            self._eval_vec["Loss"]["Prediction"].append(preloss)

            tn_reca, tn_accu, tn_prec = self._calc_eval_par(y_pred_pe, y_true_pe)
            if pe == 0 or (pe + 1) % self._log_per_n_epoch == 0:
                item_dict = dict(Epoch=pe, Reconstruction=recloss, Prediction=preloss, Accuracy=tn_accu,
                                 Precision=tn_prec, Recall=tn_reca, Partition="Train")
                self._running_log(pe, **item_dict)

        return self
    
    def predict(self, X, y=None):
        if self.model is None:
            raise ValueError("Model isn't ready yet, please use fit() first.")

        if y is not None:
            y = self._label_encoder.transform(y)

        data_loader = DataLoader(ExpDataSet(X, y), batch_size=X.shape[0], shuffle=False)
        with torch.no_grad():
            for x_mat, y_true in data_loader:
                x_mat, y_true = x_mat.to(self.usedev), y_true.to(self.usedev)
                x_enc, mu_hat, mu_exp, all_prob, x_dec = self.model.to(self.usedev)(x_mat)

                pred_prob, pred_idx = torch.max(all_prob, dim=1)
                loss_rec, loss_pre = None, None
                if y is not None:
                    recloss = self._calc_rec_loss(x_dec, x_mat, mu_hat, mu_exp).to("cpu").numpy()
                    preloss = self._calc_pre_loss(all_prob, y_true).to("cpu").numpy()
                    tt_reca, tt_accu, tt_prec = self._calc_eval_par(all_prob, y_true)
                    item_dict = dict(Epoch=0, Reconstruction=recloss, Prediction=preloss,
                                     Accuracy=tt_accu, Precision=tt_prec, Recall=tt_reca,
                                     Partition="Predict")
                    self._running_log(0, **item_dict)

                x_enc = x_enc.to("cpu").numpy()
                mu_hat = mu_hat.to("cpu").numpy()
                mu_exp = mu_exp.to("cpu").numpy()
                x_dec = x_dec.to("cpu").numpy()
                pred_idx = pred_idx.to("cpu").numpy()
                pred_prob = pred_prob.to("cpu").numpy()
                all_prob = all_prob.to("cpu").numpy()
                pred_label = self._label_encoder.inverse_transform(pred_idx)

                yield (x_enc, mu_exp, mu_hat, x_dec, all_prob, pred_idx, pred_label, pred_prob, loss_rec, loss_pre)

    def predict_proba(self, X):
        return self.predict(X)

    def encode(self, X):
        return self.predict(X)


# Train, test, predict
def train():
    pass


def predict():
    pass


def get_cliopts():
    pass


def main():
    pass


if __name__ == "__main__":
    main()

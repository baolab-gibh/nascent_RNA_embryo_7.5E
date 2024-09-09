#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: gene_regulatory_network.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Nov 01, 2023
# Updated: Nov 01, 2023

import os
import pickle
from pathlib import Path

import scanpy as sc
import pandas as pd

from dask.distributed import Client, LocalCluster
from ctxcore.rnkdb import FeatherRankingDatabase
from arboreto.algo import grnboost2
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import df2regulons
from pyscenic.prune import prune2df
from pyscenic.aucell import aucell


overwrite = False
n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", 2))

condition_list = ["T0_RPMI", "T0_LPS", "T3m_RPMI", "T3m_LPS"]
target_genes = ["CD55", "SLFN5", "ADCY3"]
target_tfs = ["SPI1", "STAT1", "ESR1", "FOXK2", "EP300", "FOXP1", "CEBPA", "KLF5", "TBX3", "HLF", "RARA", "CTCF"]

proj_dir = Path("~/Documents/projects/wp_bcg_eqtl").expanduser()
feather_dir = proj_dir / "inputs/RcisTarget_db/database"
target_tf_dir = proj_dir / "inputs/ADASTRA/TF"

# TF used in the ASTFB DB
tf_names = [x.stem.split("_")[0] for x in target_tf_dir.glob("*.tsv")]

# loom files from scRNA-seq
loom_files = [proj_dir / f"outputs/pseudo_bulk/gene_regulatory_network/data/300bcg_{ts}.loom" for ts in condition_list]

# motif annotation
motif_ann_file = proj_dir / "inputs/RcisTarget_db/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

# motif database
db_fnames = [[p.name.split("_")[1], str(p)] for p in list(feather_dir.glob("*.rankings.feather"))]
dbs = [FeatherRankingDatabase(fname=p, name=n) for n, p in db_fnames]

# Gene regulatory network from scRNA-seq per condition, check condition_list for more details
for per_name, per_loom in zip(condition_list, loom_files):
    pbmc_ad = sc.read_loom(per_loom)

    sc.pp.filter_cells(pbmc_ad, min_genes=200)
    sc.pp.filter_genes(pbmc_ad, min_cells=3)

    tf_names = [x for x in tf_names if x in pbmc_ad.var_names]
    # Expression matrix of given genes and TFs, of given cell type
    ex_matrix = pd.DataFrame(pbmc_ad.X.toarray(), index=pbmc_ad.obs_names, columns=pbmc_ad.var_names).loc[:, target_genes + tf_names]

    # Identify adjacencies
    adjacency_save_to = proj_dir / f"outputs/pseudo_bulk/gene_regulatory_network/regulons/{per_name}.adjacencies.pkl"
    if not adjacency_save_to.exists() or overwrite:
        local_cluster = LocalCluster(n_workers=n_workers)
        my_client = Client(local_cluster)
        adjacencies = grnboost2(ex_matrix, tf_names=tf_names, client_or_address=my_client, seed=31415, verbose=True)
        my_client.close()
        local_cluster.close()
        with open(adjacency_save_to, "wb") as f:
            pickle.dump(adjacencies, f)
    else:
        adjacencies = pickle.load(open(adjacency_save_to, "rb"))

    # Obtain modules
    module_save_to = proj_dir / f"outputs/pseudo_bulk/gene_regulatory_network/regulons/{per_name}.modules.pkl"
    if not module_save_to.exists() or overwrite:
        modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
        with open(module_save_to, "wb") as f:
            pickle.dump(modules, f)
    else:
        modules = pickle.load(open(module_save_to, "rb"))

    # Create regulons
    regulon_save_to = proj_dir / f"outputs/pseudo_bulk/gene_regulatory_network/regulons/{per_name}.regulons.pkl"
    if not regulon_save_to.exists() or overwrite:
        module_df = prune2df(dbs, modules, str(motif_ann_file), num_workers=n_workers, client_or_address="custom_multiprocessing")
        regulons = df2regulons(module_df)
        with open(regulon_save_to, "wb") as f:
            pickle.dump(regulons, f)
    else:
        regulons = pickle.load(open(regulon_save_to, "rb"))

    # Create AUCell matrix
    auc_mtx = aucell(ex_matrix, regulons, num_workers=n_workers)
    auc_mtx_save_to = proj_dir / f"outputs/pseudo_bulk/gene_regulatory_network/regulons/{per_name}.auc_mtx.csv"
    if not auc_mtx_save_to.exists() or overwrite:
        auc_mtx.to_csv(auc_mtx_save_to)
    else:
        auc_mtx = pd.read_csv(auc_mtx_save_to)

print("-- Done ---")

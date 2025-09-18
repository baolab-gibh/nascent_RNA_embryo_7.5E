#!/usr/bin/env python3
# File: rna_velocity.py
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Jul 01, 2024
# Updated: Aug 06, 2024

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

import pickle
from pathlib import Path

# matplotlib: version 3.5.3, matplotlib may raise errors due to depreciation of register_cmap()
import pandas as pds
import scanpy as scp

from dask.distributed import Client, LocalCluster
from ctxcore.rnkdb import FeatherRankingDatabase
from arboreto.algo import grnboost2
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import df2regulons
from pyscenic.prune import prune2df
from pyscenic.aucell import aucell
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss # Recommendations: save the data to disk and plot the result by ggplot2 in R.


PROJECT_DIR = Path("~/Documents/projects/wp_vasaseq").expanduser()


#
## Analysis, regulon
#
overwrite = False
feather_dir = PROJECT_DIR / "inputs/reference/RcisTarget_db/database"
assert feather_dir.exists()

adata = scp.read_h5ad(PROJECT_DIR / "outputs/analysis/slide_seq/adata.h5ad")


db_fnames = [[p.name.split("_")[1], str(p)] for p in list(feather_dir.glob("*.rankings.feather"))]
dbs = [FeatherRankingDatabase(fname=p, name=n) for n, p in db_fnames]

motif_ann_file = PROJECT_DIR / "inputs/reference/RcisTarget_db/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
motif_ann_tab = pds.read_table(motif_ann_file)
tf_names = [x for x in motif_ann_tab.gene_name.drop_duplicates() if x in adata.var_names]

ex_matrix = pds.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)

# Estimate adjacencies
adjacency_save_to = PROJECT_DIR / "outputs/analysis/regulon/all_batches.adjacencies.pkl"
if not adjacency_save_to.exists() or overwrite:
    local_cluster = LocalCluster(n_workers=10)
    my_client = Client(local_cluster)
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, client_or_address=my_client, seed=31415, verbose=True)
    my_client.close()
    local_cluster.close()
    with open(adjacency_save_to, "wb") as f:
        pickle.dump(adjacencies, f)
else:
    adjacencies = pickle.load(open(adjacency_save_to, "rb"))
print("-- Estimating adjacencies DONE ---")

# Crate modules
module_save_to = PROJECT_DIR / "outputs/analysis/regulon/all_batches.modules.pkl"
if not module_save_to.exists() or overwrite:
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    with open(module_save_to, "wb") as f:
        pickle.dump(modules, f)
else:
    modules = pickle.load(open(module_save_to, "rb"))
print("-- Creating modules DONE ---")

# Create regulon
regulon_save_to = PROJECT_DIR / "outputs/analysis/regulon/all_batches.regulons.pkl"
if not regulon_save_to.exists() or overwrite:
    module_df = prune2df(dbs, modules, str(motif_ann_file), num_workers=10, client_or_address="custom_multiprocessing")
    regulons = df2regulons(module_df)
    with open(regulon_save_to, "wb") as f:
        pickle.dump(regulons, f)
else:
    regulons = pickle.load(open(regulon_save_to, "rb"))
print("-- Creating regulons DONE ---")

# Create AUCell matrix
auc_mtx_save_to = PROJECT_DIR / "outputs/analysis/regulon/all_batches.auc_mtx.csv"
if not auc_mtx_save_to.exists() or overwrite:
    auc_mtx = aucell(ex_matrix, regulons, num_workers=10)
    auc_mtx.to_csv(auc_mtx_save_to)
else:
    auc_mtx = pds.read_csv(auc_mtx_save_to, index_col = 0)
print("-- Creating AUCell matrix DONE ---")

# Visualization
selected_tf = auc_mtx.gt(0).sum().ge(30).pipe(lambda x: x[x].index)

## Regulon specificity scores
meta_info = adata.obs.loc[auc_mtx.index, ["Batches", "Regions", "Layers", "Cell_types"]].copy()
rss_regions_l0 = regulon_specificity_scores(auc_mtx.loc[:, selected_tf], meta_info.Cell_types)
rss_regions_l0.to_csv(PROJECT_DIR / "outputs/analysis/regulon/all_batches.rss_regions_l0.csv")
rss_regions_l1 = regulon_specificity_scores(auc_mtx.loc[:, selected_tf], meta_info.Regions)
rss_regions_l1.to_csv(PROJECT_DIR / "outputs/analysis/regulon/all_batches.rss_regions_l1.csv")
rss_layers = regulon_specificity_scores(auc_mtx.loc[:, selected_tf], meta_info.Layers)
rss_layers.to_csv(PROJECT_DIR / "outputs/analysis/regulon/all_batches.rss_layers.csv")

# Check 
top_tf = (rss_regions_l0.T.apply(lambda x: x.sort_values(ascending=False).head(50).index).stack().reset_index(drop=True).drop_duplicates().to_list())
auc_mtx_zscore = auc_mtx.loc[:, top_tf].apply(lambda x: (x - x.mean()) / x.std())
auc_mtx_zscore.to_csv(PROJECT_DIR / "outputs/analysis/regulon/all_batches.auc_mtx_zscore.rss_regions_l0.top_50.csv")

print("-- Done ---")

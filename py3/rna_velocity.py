#!/usr/bin/env python3
# File: rna_velocity.py
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Jul 01, 2024
# Updated: Aug 06, 2024

import warnings
from pathlib import Path

# matplotlib: version 3.5.3, matplotlib may raise errors due to depreciation of register_cmap()
import numpy as np
import pandas as pd
import polars as pl
import anndata as adt
import scanpy as scp
import dynamo as dyn
import scvelo as scv
import seaborn as sbn
from scipy import sparse
from scipy import cluster
from scipy import spatial
from dynamo.preprocessing import Preprocessor
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from sklearn.decomposition import PCA

mpl.rcParams["legend.loc"] = "lower right"

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=DeprecationWarning)

version = "version_3" # Including five batches, 240409_Lib_embryo, 240612_Lib_28region, 240620_Lib_38region, 240703_Lib_32region, 240710_Lib_37region


def convert_pos(adt, into="cell_type"):
    _pos = ['A', 'P', 'EA', 'EP', 'L', 'R', 'MA', 'MP', 'O']
    if into == "regions":
        return adt.obs.index.str.extract("[0-9]+([A-Z]{1,2})").loc[:, 0].tolist()
    else:
        if into == "pseudotime":
            _target = [3, 3, 1, 1, 3, 3, 2, 2, 0]
        else:
            _target = ["Ectoderm", "Ectoderm", "Endoderm", "Endoderm", "Ectoderm", "Ectoderm", "Mesoderm", "Mesoderm", "Other"]

        pseudo_time = dict(zip(_pos, _target))
        regions = [x[0] if x else "O" for x in adt.obs.index.str.findall("[0-9]+([A-Z]{1,2})")]
        return [pseudo_time[x] for x in regions]


def init_adata(ttl_adt, new_adt=None, new_as_base=False):
    adata = new_adata.copy() if new_as_base else ttl_adt.copy()
    if new_adt is not None:
        adata.layers["new"] = new_adt.X.copy()
        adata.layers["total"] = ttl_adt.X.copy()
    adata.layers["counts"] = adata.X.copy()
    adata.obs["Groups"] = convert_pos(adata, into="pseudotime")
    adata.obs["Regions"] = convert_pos(adata, into="regions")
    adata.obs["Cell_types"] = convert_pos(adata)
    adata.obs["Layer_region"] = adata.obs.index.str.extract("([0-9]+[A-Z]{1,2})").loc[:, 0].tolist()
    adata.obs["Sampling_dates"] = adata.obs.Batches.str.extract("^([0-9]+)_").loc[:, 0].tolist()

    return adata


def calc_ntr(adata, axis=None, min_tc=10, min_nc=3, max_tc=5000, max_nc=5000, use_bf=True, fillnan=0):
    has_total, has_new = "total" in adata.layers, "new" in adata.layers
    if not (has_new and has_total): raise KeyError("Missing 'new' or 'total' key in adata.layers data.")

    tobe_incl = True
    if min_tc is not None and max_tc is not None: # Check min total counts (min_tc) and max total counts (max_tc)
        tobe_incl += (min_tc <= adata.layers["total"]).A * (adata.layers["total"] <= max_tc).A
    if min_nc is not None and max_nc is not None: # Check min new counts (min_nc) and max new counts (max_nc)
        tobe_incl *= (min_nc <= adata.layers["new"]).A * (adata.layers["new"] <= max_nc).A
    if "pass_basic_filter" in adata.var and use_bf: # Using pass_basic_filter for variables (i.e., genes)
        tobe_incl *= adata.var["pass_basic_filter"].values
    if "pass_basic_filter" in adata.obs and use_bf: # Using pass_basic_filter for observations (i.e., samples)
        tobe_incl *= adata.obs["pass_basic_filter"].values

    if axis is None:
        ntr_mat = (adata.layers["new"] / adata.layers["total"]).A * tobe_incl
    elif axis in [0, 1]:
        ntr_mat = ((adata.layers["new"].A * tobe_incl).sum(axis) / (adata.layers["total"].A * tobe_incl).sum(axis))
    else:
        raise ValueError("axis should be None for per gene per sample, 1 for per gene, 0 for per sample")

    if fillnan is not None and isinstance(fillnan, (int, float)):
        ntr_mat = np.nan_to_num(ntr_mat, nan=fillnan)

    if axis is None:
        return sparse.csr_matrix(ntr_mat)

    return ntr_mat.squeeze()



dyn.configuration.set_figure_params('dynamo', background='black')
dyn.configuration.set_pub_style()


project_dir = Path("~/Documents/projects/wp_vasaseq").expanduser()
all_batches = ["240409_Lib_embryo", "240612_Lib_28region", "240620_Lib_38region", "240703_Lib_32region", "240710_Lib_37region", "240717_Lib_28region"]
gene_bkl = ["Gm42418", "Gm26917"]

ttl_adt_list, new_adt_list, spl_adt_list = {}, {}, {}
for per_batch in all_batches:
    ten_x_dir = project_dir / "outputs/analysis/preprocessing/quantification" / per_batch / "slamseq/10X"
    new_ad, ttl_ad = scp.read_10x_mtx(ten_x_dir / "new"), scp.read_10x_mtx(ten_x_dir / "total")
    new_ad = new_ad[:, new_ad.var.gene_ids.drop_duplicates()]
    ttl_ad = ttl_ad[:, ttl_ad.var.gene_ids.drop_duplicates()]
    new_ad.obs["Batches"], ttl_ad.obs["Batches"] = per_batch, per_batch
    ttl_ad.obs["Sampling_dates"] = ttl_ad.obs.Batches.str.extract("^([0-9]+)_").loc[:, 0].tolist()
    new_ad.obs["Sampling_dates"] = new_ad.obs.Batches.str.extract("^([0-9]+)_").loc[:, 0].tolist()

    ttl_adt_list[per_batch] = ttl_ad
    new_adt_list[per_batch] = new_ad

    spl_dir = project_dir / "outputs/analysis/preprocessing/velocity_counts" / per_batch / "velocyto"
    spl_ad = dyn.read_loom(spl_dir / "one_file_per_cell.velocyto_run.loom")
    spl_ad.obs["Batches"] = per_batch
    spl_ad.obs.index = [x[1] for x in spl_ad.obs.index.str.split("[:.]")]
    spl_ad.var["gene_ids"] = spl_ad.var.index
    spl_ad.var.index = spl_ad.var.Accession
    spl_ad = spl_ad[:, spl_ad.var.gene_ids.drop_duplicates().index]
    spl_ad.var.index = spl_ad.var.gene_ids
    spl_adt_list[per_batch] = spl_ad


# Nascent RNA information
ttl_adata = adt.concat(ttl_adt_list, merge="same")
ttl_adata.obs_names_make_unique()
ttl_adata = ttl_adata[~ttl_adata.obs.index.str.contains("NC_"), [x not in gene_bkl for x in ttl_adata.var.index]]

new_adata = adt.concat(new_adt_list, merge="same")
new_adata.obs_names_make_unique()
new_adata = new_adata[~new_adata.obs.index.str.contains("NC_"), [x not in gene_bkl for x in new_adata.var.index]]

labeled_adata = init_adata(ttl_adata, new_adata)


# Add corn plot coordinates
corn_coords = pd.read_csv(project_dir / "inputs/reference/corn_coordination/corn_axis.e_7_5.csv")
corn_coords["Layer_region"] = corn_coords.apply(lambda x: str(x["y_pos"]) + x["Regions_l1"], axis=1)
corn_coords["x_pos"] = corn_coords["x_pos"] * 1000
labeled_adata.obsm["X_corn"] = pd.merge(labeled_adata.obs.copy().reset_index(), corn_coords, how="left", on="Layer_region").set_index("index").loc[:, ["y_pos", "x_pos"]].to_numpy()


# Using scanpy to check transcriptional profile using total RNAs
adata_dict = {}
for per_cat in ["new", "total"]:
    if per_cat == "new":
        tmp_adata = init_adata(ttl_adata, new_adata, new_as_base=True)
    else:
        tmp_adata = init_adata(ttl_adata, new_adata)

    tmp_adata.var["NTR"] = calc_ntr(tmp_adata, 0) # total per gene
    tmp_adata.obs["NTR"] = calc_ntr(tmp_adata, 1) # total per sample
    tmp_adata.layers["NTR"] = calc_ntr(tmp_adata, None) # total per gene per sample
    scp.pp.filter_cells(tmp_adata, min_genes=100)
    scp.pp.filter_genes(tmp_adata, min_cells=3)
    scp.pp.normalize_total(tmp_adata)
    scp.pp.log1p(tmp_adata)
    scp.pp.combat(tmp_adata, key="Batches")
    scp.pp.highly_variable_genes(tmp_adata, n_top_genes=2000, batch_key="Batches")
    scp.tl.pca(tmp_adata, n_comps=50, use_highly_variable=True, svd_solver="arpack")
    scp.pp.neighbors(tmp_adata, n_pcs=10)
    scp.tl.umap(tmp_adata)
    scp.tl.dendrogram(tmp_adata, groupby="Cell_types")
    scp.tl.rank_genes_groups(tmp_adata, groupby="Cell_types", method="wilcoxon")
    
    adata_dict[per_cat] = tmp_adata

total_adata, nascent_adata = adata_dict["total"], adata_dict["new"]
fig = plt.figure(figsize=(9, 8))
axe1, axe2, axe3, axe4 = plt.subplot(321), plt.subplot(322), plt.subplot(312), plt.subplot(313)
with mpl.rc_context({"font.size": 8}):
    _ = scp.pl.umap(total_adata, color="Cell_types", title="UMAP by total RNA", size=200, edgecolor="0.5", legend_fontoutline=2, legend_loc="on data", palette="Set2", frameon=False, ax=axe1, show=False)
    _ = scp.pl.umap(nascent_adata, color="Cell_types", title="UMAP by new RNA", size=200, edgecolor="0.5", legend_fontoutline=2, legend_loc="on data", palette="Set2", frameon=False, ax=axe2, show=False)
    _ = scp.pl.rank_genes_groups_dotplot(total_adata, groupby="Cell_types", standard_scale="var", n_genes=10, ax=axe3, show=False)
    _ = scp.pl.rank_genes_groups_dotplot(nascent_adata, groupby="Cell_types", standard_scale="var", n_genes=10, var_group_rotation=0.5, ax=axe4, show=False)
plt.subplots_adjust(hspace=-0.15, wspace=0.1, left=0.075)
fig.savefig(project_dir / "outputs/analysis/velocity" / version / "scanpy.total_and_new.umap.pdf")
fig.clear()
plt.close(fig)


# Estimate the ntr ratio
ntr_adata = init_adata(ttl_adata, new_adata)
dyn.pp.filter_genes_by_outliers(ntr_adata, min_cell_s=3, min_cell_u=3, min_count_s=10, min_count_u=3)
ntr_adata.var["pass_basic_filter"] = ntr_adata.var["pass_basic_filter"] * (ntr_adata.X <= 10000).toarray().all(0)
# dyn.pp.filter_cells_by_outliers(ntr_adata, shared_count=None) # We don't filter out any samples for the current moment.
dyn.pp.calc_sz_factor(ntr_adata)
dyn.pp.cell_cycle_scores(ntr_adata) # Bulk-data, meanlingless.
ntr_adata.var["NTR"] = calc_ntr(ntr_adata, 0) # NTR per gene
ntr_adata.obs["NTR"] = calc_ntr(ntr_adata, 1) # NTR per sample
ntr_adata.layers["NTR"] = calc_ntr(ntr_adata, None) # NTR per gene per sample
selected_genes = ntr_adata.var["pass_basic_filter"]
selected_samples = ntr_adata.obs["pass_basic_filter"] if "pass_basic_filter" in ntr_adata.obs else [True] * len(ntr_adata)
ntr_adata_sub = ntr_adata[selected_samples, selected_genes].copy()
ntr_tab = pd.DataFrame(ntr_adata_sub.layers["NTR"].A, index=ntr_adata_sub.obs_names, columns=ntr_adata_sub.var_names) #.merge(ntr_adata_sub.obs, how="inner", left_index=True, right_index=True)


# Check PCA using NTR
pca = PCA(n_components=50)
pca.fit(ntr_tab.T)
pca.get_covariance()
pc1_vp, pc2_vp, *_ = pca.explained_variance_ratio_ * 100
pca_tab = pd.DataFrame(pca.components_.T, index=ntr_tab.index, columns=[f"PC{i}" for i in range(1, 51)])

fig, (axe1, axe2) = plt.subplots(1, 2, figsize=(10, 5))
with mpl.rc_context({"legend.loc": "upper left", "legend.borderpad": 1, "legend.handletextpad": 1}):
    _ = sbn.scatterplot(data=pca_tab, x="PC1", y="PC2", hue=ntr_adata.obs["Cell_types"], palette="deep", edgecolor="0", alpha=0.75, ax=axe1)
    _ = sbn.scatterplot(data=pca_tab, x="PC1", y="PC2", hue=ntr_adata.obs["Batches"], palette="Set2", edgecolor="0", alpha=0.75, ax=axe2)
_ = axe1.set_xlabel(f"PC1 ({pc1_vp:.1f}%)")
_ = axe1.set_ylabel(f"PC2 ({pc2_vp:.1f}%)")
_ = axe1.set_title("PCA plot by NTR, colored by cell types")
_ = axe1.spines["top"].set_visible(False)
_ = axe1.spines["right"].set_visible(False)
_ = axe2.set_xlabel(f"PC1 ({pc1_vp:.1f}%)")
_ = axe2.set_ylabel(f"PC2 ({pc2_vp:.1f}%)")
_ = axe2.set_title("PCA plot by NTR, colored by sequencing batches")
_ = axe2.spines["top"].set_visible(False)
_ = axe2.spines["right"].set_visible(False)
plt.subplots_adjust(left=0.075)
fig.savefig(project_dir / "outputs/analysis/velocity" / version / "ntr.pca.pdf")
fig.clear()
plt.close(fig)


# Cluster of samples based on NTR
batches = ['240409_Lib_embryo', '240612_Lib_28region', '240620_Lib_38region', '240703_Lib_32region', '240710_Lib_37region', '240717_Lib_28region']
cell_types = ["Ectoderm", "Endoderm", "Mesoderm"]
colors = mpl.cm.get_cmap("Dark2").colors[:3]
color_map = dict(zip(cell_types, colors))
network_colors = {k: color_map[v] for k, v in ntr_adata_sub.obs.loc[:, "Cell_types"].to_dict().items()}
network_colors = [network_colors[k] for k in pca_tab.index]

pca_corr = pca_tab.T.corr()
grid = sbn.clustermap(
    pca_corr, method="ward", center=0, cmap="RdBu_r", row_colors=network_colors, figsize=(5, 5), cbar_pos=(0, 0.075, 0.125, 0.025),
    cbar_kws={"label": "Correlation", "orientation": "horizontal"},
)
_ = grid.ax_heatmap.legend(
    handles=[mpatches.Patch(color=v, label=k) for k, v in color_map.items()], title="Cell type", loc="upper left",
    bbox_to_anchor=(-.35, 0.1), fontsize="large", title_fontsize="x-large", borderpad=1.25
)
grid.ax_col_dendrogram.remove()
grid.savefig(project_dir / "outputs/analysis/velocity" / version / "in_house.new_to_total_ratio.cluster_by_pca.pdf")





gene_dist = spatial.distance.pdist(ntr_tab.values.T)
gene_z = cluster.hierarchy.linkage(gene_dist, "ward")
sample_dist = spatial.distance.pdist(ntr_tab.values)
sample_z = cluster.hierarchy.linkage(sample_dist, "ward")

# Heatmap to show clusters
fig_size = (8, 6)
gridspec_kws = dict(width_ratios=[1, 8, 1], height_ratios=[1, 8, 1])
axe_keys = [["place holder 0", "gene level NTR", "place holder 1"], ['sample dendrogram', 'heatmap', "sample level NTR"], ['place holder 3', 'gene dendrogram', 'place holder 4']]
fig, axd = plt.subplot_mosaic(axe_keys, gridspec_kw=gridspec_kws, figsize=fig_size, layout="constrained")
with mpl.rc_context({"lines.linewidth": 0.5, "lines.color": "black"}):
    for k, axe in axd.items():
        axe.set_frame_on(False)
        axe.xaxis.set_visible(False)
        axe.yaxis.set_visible(False)
        if k == "sample dendrogram":
            _ = cluster.hierarchy.dendrogram(sample_z, ax=axe, orientation="left", link_color_func=lambda _: "black")
        elif k in ["heatmap", "sample level NTR"]:
            sample_order = ntr_tab.index[cluster.hierarchy.leaves_list(sample_z)].to_list()
            gene_order = ntr_tab.columns[cluster.hierarchy.leaves_list(gene_z)].to_list()
            ntr_tab_sub = ntr_tab.loc[sample_order, gene_order]
            if k == "sample level NTR":
                x_vals = ntr_tab_sub.values.T.mean(0)
                y_vals = ntr_tab_sub.index.tolist()
                _ = axe.plot(x_vals, y_vals)
                axe.set_ylim([1, 172])
                axe.xaxis.set_visible(True)
                axe.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
            else:
                _ = axe.pcolormesh(ntr_tab_sub.values, cmap="Blues", rasterized=True)
        elif k == "gene dendrogram":
            _ = cluster.hierarchy.dendrogram(gene_z, ax=axe, orientation="bottom", link_color_func=lambda _: "black")
    
fig.savefig(project_dir / "outputs/analysis/velocity" / version / "in_house.new_to_total_ratio.pdf")
fig.clear()
plt.close(fig)


# Splicing information
spl_adata = adt.concat(spl_adt_list, merge="same")
spl_adata.obs_names_make_unique()
spl_adata = spl_adata[~spl_adata.obs.index.str.contains("NC_"), :]
spliced_adata = spl_adata[:, labeled_adata.var.index.tolist()]
spliced_adata.obs["Groups"] = convert_pos(spliced_adata, into="pseudotime")
spliced_adata.obs["Regions"] = convert_pos(spliced_adata, into="regions")
spliced_adata.obs["Cell_types"] = convert_pos(spliced_adata)
spliced_adata.obs["Layer_region"] = spliced_adata.obs.index.str.extract("([0-9]+[A-Z]{1,2})").loc[:, 0].tolist()
spliced_adata.obs["Sampling_dates"] = spliced_adata.obs.Batches.str.extract("^([0-9]+)_").loc[:, 0].tolist()


#
## Velocity using splicing information
#
scv.pp.filter_and_normalize(spliced_adata)
scp.pp.pca(spliced_adata, n_comps=20)
scp.pp.neighbors(spliced_adata, n_pcs=20, n_neighbors=10)
scv.pp.moments(spliced_adata, n_pcs=0, n_neighbors=0)
scp.tl.umap(spliced_adata)
scv.tl.velocity(spliced_adata, mode='deterministic')
scv.tl.velocity_graph(spliced_adata)

fig, (axe1, axe2, axe3) = plt.subplots(1, 3, figsize=(12, 4.5))
scv.pl.velocity_embedding(spliced_adata, color="Cell_types", size=600, alpha=0.7, basis='umap', ax=axe3)
scv.pl.velocity_embedding_grid(spliced_adata, color="Cell_types", size=600, alpha=0.7, basis='umap', ax=axe2)
scv.pl.velocity_embedding_stream(spliced_adata, color="Cell_types", size=600, alpha=0.7, basis='umap', ax=axe1)
x_lim = axe3.get_xlim()
y_lim = axe3.get_ylim()
_ = axe2.set_xlim(x_lim)
_ = axe2.set_ylim(y_lim)
_ = axe1.set_xlim(x_lim)
_ = axe1.set_ylim(y_lim)
fig.savefig(project_dir / "outputs/analysis/velocity" / version / "scvelo.splicing_based.umap.pdf")
fig.clear()

# Save the anndata into h5ad format
spliced_adata.write_h5ad(project_dir / "outputs/analysis/velocity" / version / "scvelo.velocity_analysis.v2.h5ad")


#
## Velocity using labeling/unlabeling data
#
# # Using scanpy to check the batch effects
# labeled_adata.var["rbio"] = labeled_adata.var_names.str.startswith(("Rps", "Rpl"))
# scp.pp.filter_cells(labeled_adata, min_genes=100)
# scp.pp.filter_genes(labeled_adata, min_cells=3)
# scp.pp.normalize_total(labeled_adata)
# scp.pp.log1p(labeled_adata)
# scp.pp.highly_variable_genes(labeled_adata, n_top_genes=2000, batch_key="Batches")
# scp.tl.pca(labeled_adata)
# scp.pl.pca_variance_ratio(labeled_adata, n_pcs=50, log=True, save="")
# scp.pp.neighbors(labeled_adata)
# scp.tl.umap(labeled_adata)
# fig, axe = plt.subplots(1, 1, constrained_layout=True, figsize=(5, 5))
# scp.pl.umap(labeled_adata, color="Batches", size=40, ax=axe)
# fig.savefig(project_dir / "outputs/analysis/velocity" / version / "scanpy.check_batch_effects.umap.pdf")
# fig.clear()

# Preprocessing.
pp = Preprocessor(cell_cycle_score_enable=True)
pp.config_monocle_recipe(labeled_adata)
pp.preprocess_adata(labeled_adata, recipe="monocle", tkey="Groups", experiment_type="one-shot")
scp.pp.combat(labeled_adata, key="Batches")
pp.pca(labeled_adata)
dyn.tl.reduceDimension(labeled_adata, basis="pca")

# Velocities
dyn.tl.dynamics(labeled_adata)
dyn.tl.cell_velocities(labeled_adata, basis="umap", method='pearson', other_kernels_dict={'transform': 'sqrt'})
dyn.tl.cell_wise_confidence(labeled_adata, ekey = "X_total", vkey = "velocity_N")
dyn.tl.confident_cell_velocities(labeled_adata, group="Cell_types", ekey="total", vkey="velocity_N", lineage_dict={'Endoderm': 'Ectoderm'})
dyn.vf.VectorField(labeled_adata, basis='umap', vkey="velocity_N", M=1000, pot_curl_div=True)
dyn.vf.rank_velocity_genes(labeled_adata, vkey="velocity_N")

# Vector filed reconstruction
dyn.tl.cell_velocities(labeled_adata, basis="pca", method='pearson', other_kernels_dict={'transform': 'sqrt'})
dyn.vf.VectorField(labeled_adata, basis='pca', vkey="velocity_N", M=1000, pot_curl_div=True)
dyn.vf.speed(labeled_adata, basis='pca')
dyn.vf.curl(labeled_adata, basis='umap')
dyn.vf.divergence(labeled_adata, basis='pca')
dyn.vf.acceleration(labeled_adata, basis='pca')
dyn.vf.curvature(labeled_adata, basis='pca')

# Plot vector field energy alterations
fig, axe = plt.subplots(1, 1, constrained_layout=True, figsize=(10, 4))
p = dyn.pl.plot_energy(labeled_adata, basis='umap', fig=fig)
fig.savefig(project_dir / "outputs/analysis/velocity" / version / "dynamo.energy_and_energy_change_rate.pdf")
fig.clear()

# Velocity in the corn space
labeled_adata.obsm["X_corn"] = labeled_adata.obsm["X_umap"]
labeled_adata.obsm["X_corn_SparseVFC"] = labeled_adata.obsm["X_umap_SparseVFC"]
labeled_adata.obsm["velocity_corn"] = labeled_adata.obsm["velocity_umap"]
labeled_adata.obsm["velocity_corn_SparseVFC"] = labeled_adata.obsm["velocity_umap_SparseVFC"]

# Visualize the results in corn plot
fig, (axe1, axe2) = plt.subplots(1, 2, constrained_layout=True, figsize=(12, 5))
dyn.pl.topography(labeled_adata, basis='corn', x=1, y=0, color='Cell_types', pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, frontier=True, ekey="M_t", vkey="velocity_N", ax=axe1, save_show_or_return="return") # Plot topography
x_lim = axe1.get_xlim()
y_lim = axe1.get_ylim()
dyn.pl.streamline_plot(labeled_adata, basis="corn", x=1, y=0, color="Cell_types", pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, ekey="M_t", vkey="velocity_N", ax=axe2, save_show_or_return="return") # Plot velocity
axe.set_ylim(y_lim)
axe.set_xlim(x_lim)
fig.savefig(project_dir / "outputs/analysis/velocity" / version / "dynamo.corn_plot.pdf")
fig.clear()

# Vector field topography and velocity stream line
fig, (axe1, axe2) = plt.subplots(1, 2, constrained_layout=True, figsize=(12, 5))
dyn.pl.topography(labeled_adata, basis='umap', color='Cell_types', pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, frontier=True, ekey="M_n", vkey="velocity_N", ax=axe1, save_show_or_return="return", show_arrowed_spines=True) # Plot topography
x_lim = axe1.get_xlim()
y_lim = axe1.get_ylim()
_ = axe1.set_title("Vector field topography")
dyn.pl.streamline_plot(labeled_adata, basis="umap", color="Cell_types", pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, frontier=True, ekey="M_n", vkey="velocity_N", ax=axe2, save_show_or_return="return", show_arrowed_spines=True) # Plot velocity
_ = axe2.set_xlim(x_lim)
_ = axe2.set_ylim(y_lim)
_ = axe2.set_title("Velocity stream lines")
fig.savefig(project_dir / "outputs/analysis/velocity" / version / "dynamo.topography_and_streamline.pdf")
fig.clear()

# Other dynamics
fig, ((axe1, axe2), (axe3, axe4)) = plt.subplots(2, 2, constrained_layout=True, figsize=(12, 10))
dyn.pl.cell_wise_vectors(labeled_adata, basis='umap', color=['Cell_types'], ax=axe1, show_legend='on data', quiver_length=4, quiver_size=4, pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, show_arrowed_spines=False)
axe1.set_xlim(x_lim)
axe1.set_ylim(y_lim)
dyn.pl.grid_vectors(labeled_adata, color='divergence_pca', ax=axe2, quiver_length=4, quiver_size=4, pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, save_show_or_return='return')
axe2.set_xlim(x_lim)
axe2.set_ylim(y_lim)
dyn.pl.streamline_plot(labeled_adata, color='acceleration_pca', ax=axe3, pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, save_show_or_return='return')
axe3.set_xlim(x_lim)
axe3.set_ylim(y_lim)
dyn.pl.streamline_plot(labeled_adata, color='curvature_pca', ax=axe4, pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, save_show_or_return='return')
axe4.set_xlim(x_lim)
axe4.set_ylim(y_lim)
fig.savefig(project_dir / "outputs/analysis/velocity" / version / "dynamo.integrative_analysis.pdf")
fig.clear()


# Check marker genes
marker_genes = ["Sox2", "Sox17", "Mesp1", "Hnf4a"]
fig = plt.figure(figsize=(12, 6), frameon=False, clear=True, layout="tight")
axe1, axe2, axe3, axe4, axe5 = plt.subplot(243), plt.subplot(244), plt.subplot(247), plt.subplot(248), plt.subplot(121)

dyn.pl.streamline_plot(labeled_adata, color="Cell_types", basis="umap", pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, ax=axe5, save_show_or_return="return")
axe5.set_title("Velocity by dynamo")
axe5.set_xlim(x_lim)
axe5.set_ylim(y_lim)

for per_gene, per_axe in zip(marker_genes, [axe1, axe2, axe3, axe4]):
    dyn.pl.umap(labeled_adata, color=per_gene, pointsize=0.125, alpha=0.7, ax=per_axe, save_show_or_return="return")
    per_axe.set_title(f"Gene expression {per_gene}")
    per_axe.set_xlim(x_lim)
    per_axe.set_ylim(y_lim)
fig.savefig(project_dir / "outputs/analysis/velocity" / version / "dynamo.marker_gene_expression.pdf")


# In silico perturbation
gene = "Sox17"
dyn.pd.perturbation(labeled_adata, gene, [-100], emb_basis="umap")

fig, ((axe1, axe2), (axe3, axe4)) = plt.subplots(2, 2, figsize=(10, 10), layout="tight", clear=True, frameon=False)
dyn.pl.streamline_plot(labeled_adata, color="Cell_types", basis="umap", pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, ax=axe1, save_show_or_return="return")
axe1.set_title("Velocity before perturbation")
axe1.set_xlim(x_lim)
axe1.set_ylim(y_lim)
dyn.pl.umap(labeled_adata, color='umap_ddhodge_potential', pointsize=0.5, alpha=0.7, frontier=True, ax=axe2, save_show_or_return="return")
axe2.set_title("Dyanamics")
axe2.set_xlim(x_lim)
axe2.set_ylim(y_lim)
dyn.pl.streamline_plot(labeled_adata, color=gene, basis="umap", pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, ax=axe3, save_show_or_return="return")
axe3.set_title(f"Gene expression {gene}")
axe3.set_xlim(x_lim)
axe3.set_ylim(y_lim)
dyn.pl.streamline_plot(labeled_adata, color="Cell_types", basis="umap_perturbation", pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, ax=axe4, save_show_or_return="return")
axe4.set_title("Velocity after perturbation")
axe4.set_xlim(x_lim)
axe4.set_ylim(y_lim)
fig.savefig(project_dir / "outputs/analysis/velocity" / version / ("dynamo.perturbation_" + gene + ".pdf"))
fig.clear()


#
## Misc
#
fig, axe = plt.subplots(1, 1, figsize=(18, 18))
color_map = {"Ectoderm": "red", "Endoderm": "blue", "Mesoderm": "green"}
colors = labeled_adata.obs.Cell_types.apply(lambda x: color_map[x]).tolist()
axe.scatter(x = labeled_adata.obsm["X_umap"][:, 0], y = labeled_adata.obsm["X_umap"][:, 1], c=colors, s=500, alpha=0.5)
labels = [".".join([x, y]) for x, y in zip(labeled_adata.obs.Layer_region, labeled_adata.obs.Sampling_dates)]
for px, py, pt in zip(labeled_adata.obsm["X_umap"][:, 0], labeled_adata.obsm["X_umap"][:, 1], labels):
    axe.text(px, py, pt, alpha=0.75, fontsize="x-small", ha="center", va="center")
fig.savefig(project_dir / "outputs/analysis/velocity" / version / "dynamo.outliers.pdf")
fig.clear()

fig, axe = plt.subplots(1, 1, figsize=(5, 7))
color_map = {"Ectoderm": "red", "Endoderm": "blue", "Mesoderm": "green"}
colors = labeled_adata.obs.Cell_types.apply(lambda x: color_map[x]).tolist()
axe.scatter(x = labeled_adata.obsm["X_corn"][:, 1], y = labeled_adata.obsm["X_corn"][:, 0], c=colors, s=500, alpha=0.5)
labels = [".".join([x, y]) for x, y in zip(labeled_adata.obs.Layer_region, labeled_adata.obs.Sampling_dates)]
for px, py, pt in zip(labeled_adata.obsm["X_corn"][:, 1], labeled_adata.obsm["X_corn"][:, 0], labels):
    axe.text(px, py, pt, alpha=0.5, fontsize="x-small", ha="center", va="center")
fig.savefig(project_dir / "outputs/analysis/velocity" / version / "dynamo.corn_plot.sample_positions.pdf")
fig.clear()


# Save the anndata into h5ad format
labeled_adata.write_h5ad(project_dir / "outputs/analysis/velocity" / version / "dynamo.velocity_analysis.h5ad")
dyn.session_info()

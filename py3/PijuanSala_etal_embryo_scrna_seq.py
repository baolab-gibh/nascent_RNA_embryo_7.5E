#!/usr/bin/env python3
# File: PijuanSala_etal_embryo_scrna_seq.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Feb 12, 2025
# Updated:

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from pathlib import Path

import scvi
import torch
import anndata as adt
import polars as pls
import scanpy as scp
import seaborn as sbn

import matplotlib as mpl
import matplotlib.pyplot as plt


project_dir = Path("~/Documents/projects/wp_vasaseq/").expanduser()
torch.set_float32_matmul_precision("high")
stage_list = ['E6.5', 'E6.75', 'E7.0', 'E7.25', 'E7.5', 'E7.75', 'E8.0', 'E8.25', 'E8.5']


#
## Load raw read counts matrix and meta data from disk.
#
h5ad_file = project_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/raw_feature_bc_matrix.h5ad"
if h5ad_file.exists():
    adata = scp.read_h5ad(h5ad_file)
else:
    # FIXME: the features.csv.gz should have three columns: gene_symbols, gene_id, feature_type. The third column is
    # missing in the original data.
    # Load 10X data from disk.
    read_count_dir = project_dir / "inputs/PijuanSala_Nature_2019/10X/outs/raw_feature_bc_matrix"
    adata = scp.read_10x_mtx(read_count_dir)

    # Adding meta data.
    metadata_file = project_dir / "inputs/PijuanSala_Nature_2019/10X/meta.csv"
    metadata = pls.read_csv(metadata_file).to_pandas().assign(cell_barcode = adata.obs_names)
    metadata["sequencing.batch"] = metadata["sequencing.batch"].astype("category")
    adata.obs = metadata.set_index("cell_barcode")

    adata.write_h5ad(h5ad_file)
# adata.raw = adata.copy()
print("Loaded adata: {} cells x {} genes".format(*adata.shape))


#
## Preprocessing
#

h5ad_saveto = project_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/PijuanSala_etal_Nature_2019.processed.h5ad"
if h5ad_saveto.exists():
    adata = scp.read_h5ad(h5ad_saveto)
else:
    # Calculate QC matrics
    adata.var["mt"] = adata.var_names.str.contains("^mt-")
    adata.var["ribo"] = adata.var_names.str.contains("^Rp[ls]")
    scp.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=False)
    n_mt_genes, n_ribo_genes = adata.var["mt"].sum(), adata.var["ribo"].sum()
    print(f"Observed {n_mt_genes} MT genes and {n_ribo_genes} ribosomal genes")

    # Filtering genes and cells
    scp.pp.filter_cells(adata, min_genes=200)
    scp.pp.filter_genes(adata, min_cells=10)
    n_cells, n_genes = adata.shape
    print(f"Available gene and cells after filtering: {n_cells} cells x {n_genes} genes")

    # Normalization and scaling.
    scp.pp.normalize_total(adata, target_sum=1e6) # Using CPM
    scp.pp.log1p(adata)
    # scp.pp.scale(adata)

    # Identify highly variable genes
    scp.pp.highly_variable_genes(adata, batch_key="sequencing.batch")
    n_hvg = adata.var.highly_variable.sum()
    example_hvg = adata.var_names[adata.var.highly_variable][:5].str.cat(sep=", ")
    print(f"Found {n_hvg} HVGs, e.g. {example_hvg} ...")

    # Identify doublets, check the "doublet" in adata.obs
    # scp.pp.scrublet(adata). Not necessary

    # Keep cells with clear stage label and non-doublet.
    adata.obs.loc[:, "celltype_stage"] = adata.obs.celltype.astype(str) + "_" + adata.obs.stage.astype(str)
    selected_groups = (
        pls.from_pandas(adata.obs.loc[:, ["celltype_stage"]])
        .group_by("celltype_stage")
        .agg(pls.count())
        .filter(pls.col("count") >= 10)
        .get_column("celltype_stage")
        .to_list()
    )
    selected_cells = (
        (adata.obs.stage != "mixed_gastrulation") &
            (adata.obs.doublet == False) &
            (adata.obs.celltype != "NA") &
            adata.obs.celltype_stage.isin(selected_groups)
    )
    adata = adata[selected_cells, :]
    n_cells, _ = adata.shape
    print(f"After removing mixed gastrulation and doublets, {n_cells} cells remains.")

    # Order the stages
    adata.obs["stage"] = adata.obs.stage.cat.reorder_categories(stage_list, ordered=True)

    # Decomposition and visualization
    scp.pp.pca(adata, n_comps=50) # Create PCA
    scp.pp.neighbors(adata, n_neighbors=30, n_pcs=50) # Find neighbors
    scp.tl.umap(adata) # Create UMAP
    scp.tl.diffmap(adata, n_comps=30) # Create disfusion map

    fig, (axe1, axe2) = plt.subplots(2, 1, figsize=(14, 14), layout="constrained")
    scp.pl.umap(adata, color="celltype", size = 40, ax=axe1, title=None, add_outline=True)
    scp.pl.umap(adata, color="stage", size=40, ax=axe2, title=None, add_outline=True)
    fig.savefig(project_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/overview.umap_pca_diffmap.pdf")
    fig.clear()
    plt.close()

    fig, axe_list = plt.subplots(3, 3, figsize=(12, 12), layout="constrained", tight_layout=True, sharex=True, sharey=True)
    for idx, per_stage in enumerate(stage_list):
        n_row, n_col = int(idx / 3), idx % 3
        per_axe = axe_list[n_row, n_col]
        scp.pl.umap(adata[adata.obs.stage == per_stage, :], color="celltype", size=40, ax=per_axe, title=per_stage, add_outline=True)
        per_axe.legend([], [])
    fig.savefig(project_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/overview.umap_by_stage.pdf")
    fig.clear()
    plt.close()

    # Identify cell type by stage markers
    scp.tl.rank_genes_groups(adata, groupby="stage", method="wilcoxon", use_raw=False, pts=True, key_added="by_stage")
    scp.tl.rank_genes_groups(adata, groupby="celltype", method="wilcoxon", use_raw=False, pts=True, key_added="by_celltype")
    scp.tl.rank_genes_groups(adata, groupby="celltype_stage", method="wilcoxon", use_raw=False, pts=True, key_added="by_celltype_stage")

    # Save the adata in h5ad format
    adata.write_h5ad(h5ad_saveto)


#
## Create a model
#
deg_tab = pls.DataFrame(scp.get.rank_genes_groups_df(adata, "E6.5", key="by_stage"))
deg_tab.filter(pls.col("pct_nz_group") > 0.7, pls.col("pct_nz_reference") < 0.3)

with mpl.rc_context({"grid.linewidth": 0}):
    fig, axe = plt.subplots(1, 1, figsize=(6, 10), tight_layout=True)
    scp.pl.rank_genes_groups_matrixplot(
        adata, key="by_stage", n_genes=100, categories_order=stage_list, ax=axe, values_to_plot="logfoldchanges", cmap='RdYlBu_r',
        vmin=-4, vmax=4, min_logfoldchange=1, dendrogram=False, colorbar_title="Log(FC)", var_group_rotation=0, swap_axes=True
    )
    axe.set_xticks([])
    axe.set_xticklabels([])
    fig.savefig(project_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/overview.rank_genes_groups.pdf")
plt.clf()
plt.close()

deg_tab_list = []
for per_stage in stage_list:
    # per_stage = "E6.5"
    deg_this_stage_tab = (
        pls.from_pandas(scp.get.rank_genes_groups_df(adata, per_stage, key="by_stage"))
        .filter((pls.col("pvals_adj") < 0.05), (pls.col("logfoldchanges").abs() >= 1), ~pls.col("names").str.contains("Rik$"))
        .filter((pls.col("logfoldchanges").rank(descending=True) <= 100))
        .with_columns(pls.lit(per_stage).alias("group"))
        .select("names", "group", "logfoldchanges")
        .pivot("group", index="names", values="logfoldchanges")
    )
    selected_features = deg_this_stage_tab.get_column("names").unique().to_list()

    rest_stage_list = [x for x in stage_list if x != per_stage]
    deg_rest_stage_tab = (
        pls.from_pandas(scp.get.rank_genes_groups_df(adata, rest_stage_list, key="by_stage"))
        .filter(pls.col("names").is_in(selected_features))
        .select("names", "group", "logfoldchanges")
        .pivot("group", index="names", values="logfoldchanges")
    )

    deg_tab = deg_this_stage_tab.join(deg_rest_stage_tab, on="names", how="left").select("names", *stage_list).sort(per_stage, descending=False)
    deg_tab_list.append(deg_tab)
plt_tab = pls.concat(deg_tab_list, how="vertical").to_pandas().set_index("names")

fig, axe = plt.subplots(1, 1, figsize=(6, 10), tight_layout=True)
_ = sbn.heatmap(plt_tab, vmin=-7, vmax=7, center=0, cmap="Spectral_r", ax=axe, yticklabels=[])
fig.savefig(project_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/overview.rank_genes_groups_scatter.pdf")
plt.clf()
plt.close()

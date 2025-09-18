#!/usr/bin/env python3
# File: integration.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Apr 21, 2025
# Updated:

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from pathlib import Path

import scvi
import torch
import polars as pls
import pandas as pds
import anndata as adt
import scanpy as scp
import scvelo as scv
import seaborn as sbn

import matplotlib as mpl
import matplotlib.pyplot as plt

scvi.settings.seed = 0
plt.rcParams['figure.constrained_layout.use'] = True

def plots_velocity(adata, figsize=(8, 8), color_by="Cell_type", batch_key="sample_id", save_to=None):
    with mpl.rc_context({"figure.dpi": 300}):
        fig, ((axe1, axe2), (axe3, axe4)) = plt.subplots(2, 2, figsize=figsize, constrained_layout=True, tight_layout=True)

        if color_by not in adata.obs.columns:
            "color_by is not available in the meta data of adata.obs"
            color_by = None
        scv.pl.umap(adata, color=batch_key, ax=axe1)
        scv.pl.umap(adata, color=color_by, ax=axe2)
        scv.pl.velocity_embedding_stream(adata, color=color_by, ax=axe3)
        scv.pl.velocity_graph(adata, color=color_by, ax=axe4)

        title_maps = zip([axe1, axe2, axe3, axe4], ["Batches", f"Clusters by {color_by}", "Velocity streams", "Velocity graph"])
        for per_axe, per_title in title_maps:
            per_axe.set_title(per_title)

    if save_to is None: save_to = "velocity.pdf"
    fig.savefig(save_to)
    fig.clear()
    plt.close()


REGION_ORDER = ["EA", "MA", "A", "L", "R", "P", "MP", "EP"]
PROJECT_DIR = Path("~/Documents/projects/wp_vasaseq").expanduser()
ALL_BATCHES = [
    "240409_Lib_embryo", "240612_Lib_28region", "240620_Lib_38region", "240703_Lib_32region",
    "240710_Lib_37region", "240717_Lib_28region"
]

CELLTYPE_DICT = {
    'ExE endoderm': "Endoderm", 'Def. endoderm': "Endoderm", 'Parietal endoderm': "Endoderm", 'Visceral endoderm': "Endoderm",
    'ExE mesoderm': "Mesoderm", 'Mixed mesoderm': "Mesoderm", 'Intermediate mesoderm': "Mesoderm", 'Pharyngeal mesoderm': "Mesoderm", 'Nascent mesoderm': "Mesoderm", 'Paraxial mesoderm': "Mesoderm", 'Somitic mesoderm': "Mesoderm", 'Caudal Mesoderm': "Mesoderm",
    'ExE ectoderm': "Ectoderm", 'Surface ectoderm': "Ectoderm", 'Rostral neurectoderm': "Ectoderm", 'Caudal neurectoderm': "Ectoderm",
    'Caudal epiblast': "Epiblast", 'Epiblast': "Epiblast",
    'Primitive Streak': "Primitive Streak", 'Anterior Primitive Streak': "Primitive Streak",
    # 'Haematoendothelial progenitors', 'Blood progenitors 2', 'Blood progenitors 1',
    # 'PGC', 'Gut', 'Notochord', 'Allantois', 'Mesenchyme',
}

# Gene length
gtf_path = Path("~/Documents/projects/resources/references/gencode.vM29lift38.annotation.gtf").expanduser()
gene_len_list = (
    pls.read_csv(gtf_path, separator="\t", comment_prefix="#", has_header=False)
    .filter(pls.col("column_3") == "gene")
    .with_columns(pls.col("column_9").str.extract("gene_name \"(.*?)\";", 1).alias("gene_name"))
    .with_columns((pls.col("column_5") - pls.col("column_4")).alias("length"))
    .select("gene_name", "length")
)
gene_length_dict = dict(zip(gene_len_list["gene_name"], gene_len_list["length"]))
del gene_len_list


# GEO-VASA-SLAM by Pan et al., X 2025
# gsv_path = PROJECT_DIR / "outputs/analysis/preprocessing/geo_vasa_slam/240409_Lib_embryo.raw_data.h5ad"
# gsv_embryo = scp.read_h5ad(gsv_path)

# VASA-seq by Salmen et al., NBT, 2022
vasa_path = PROJECT_DIR / "outputs/analysis/public/Salmen_etal_NBT_2022/velocity_v2/raw.h5ad"
vasa_embryo = scp.read_h5ad(vasa_path)
vasa_embryo.obs_names_make_unique()
vasa_embryo.obs["barcode"] = vasa_embryo.obs_names
vasa_embryo.obs["celltype"] = "Unknown"
vasa_embryo.obs["batches"] = vasa_embryo.obs["sample_id"].astype(str) + "_vasaseq"
vasa_embryo.obs["sequencing_tech"] = "vasaseq"
vasa_embryo.var_names = vasa_embryo.var_names.str.extract(r"^.+_(.+)_.+$")[0].to_list()
vasa_embryo.var["length"] = [gene_length_dict.get(gene, 10) for gene in vasa_embryo.var_names]
nodup = vasa_embryo.var_names.duplicated()
vasa_embryo = vasa_embryo[:, ~nodup]
vasa_embryo.raw = vasa_embryo # Store the raw data


# Single-cell RNA-seq by Pijuan-Sala et al., Nature 2019
sc_path = PROJECT_DIR / "outputs/analysis/public/PijuanSala_etal_Nature_2019/anndata/PijuanSala_etal_Nature_2019.raw.h5ad"
sc_embryo = scp.read_h5ad(sc_path)
sc_embryo.obs["batches"] = sc_embryo.obs["sample"].astype(str) + "_scseq"
sc_embryo.obs["sequencing_tech"] = "scseq"
sc_embryo.var["length"] = [gene_length_dict.get(gene, 10) for gene in sc_embryo.var_names] # Adding gene length information
sc_embryo = sc_embryo[sc_embryo.obs.stage.isin(["E7.25", "E7.5", "E7.75"]), :].copy() # Keep only cells from E7.x
sc_embryo.raw = sc_embryo # Store the raw data


# Remove bias due to gene length
vasa_embryo.X = vasa_embryo.X / vasa_embryo.raw.var["length"].to_numpy() * vasa_embryo.var["length"].median()
sc_embryo.X = sc_embryo.X / sc_embryo.raw.var["length"].to_numpy() * sc_embryo.var["length"].median()


# Merging
embryo = scp.concat([vasa_embryo, sc_embryo], join="inner")
# del vasa_embryo
del sc_embryo

embryo.layers["counts"] = embryo.X.copy()
scp.pp.normalize_total(embryo, target_sum=1e4)
scp.pp.log1p(embryo)
embryo.raw = embryo  # keep full dimension safe
scp.pp.highly_variable_genes(embryo, flavor="seurat_v3", n_top_genes=2000, layer="counts", batch_key="sequencing_tech", subset=True,)


#
## Integration using scvi-tools
#
# model training
torch.set_float32_matmul_precision("high")
embryo.obs["celltype"] = embryo.obs["celltype"].apply(lambda x: "Unknown" if x == "Uknown" else x).to_list()

scvi.model.SCVI.setup_anndata(embryo, layer="counts", batch_key="sequencing_tech")
scvi_model = scvi.model.SCVI(embryo, n_layers=2, n_latent=30)
scvi_model.train()
embryo.obsm["X_scVI"] = scvi_model.get_latent_representation()

scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, adata=embryo, unlabeled_category="Unknown", labels_key="celltype")
scanvi_model.train(max_epochs=20, n_samples_per_label=100)
embryo.obsm["X_scANVI"] = scanvi_model.get_latent_representation(embryo)
embryo.obs["C_scANVI"] = scanvi_model.predict(embryo)

# Get latent representation and create UMAP
scp.pp.neighbors(embryo, use_rep="X_scVI")
scp.tl.umap(embryo, min_dist=0.3)
scp.tl.embedding_density(embryo, groupby="sequencing_tech")

# Visualization
fig, axe = plt.subplots(1, 1, figsize=(7, 5), tight_layout=True)
scp.pl.umap(embryo, color="sequencing_tech", size=50, ax=axe, alpha=0.75)
fig.savefig(PROJECT_DIR / "outputs/analysis/integration/embryo.umap.scvi_integration.by_sequencing_tech.pdf")
plt.close()

fig = scp.pl.embedding_density(embryo, groupby="sequencing_tech", alpha=0.75, return_fig=True)
fig.set_size_inches(13, 5)
fig.savefig(PROJECT_DIR / "outputs/analysis/integration/embryo.density.scvi_integration.by_sequencing_tech.pdf")
plt.close()

fig, axe = plt.subplots(1, 1, figsize=(13, 5))
scp.pl.umap(embryo, color="C_scANVI", size=50, ax=axe, alpha=0.5)
fig.savefig(PROJECT_DIR / "outputs/analysis/integration/embryo.umap.scanvi_annotations.pdf")
plt.close()


#
## Velocity after celltype annotations
#
# Adding cell type annotations
vasa_path = PROJECT_DIR / "outputs/analysis/public/Salmen_etal_NBT_2022/velocity_v2/raw.h5ad"
vasa_embryo = scp.read_h5ad(vasa_path)
vasa_embryo.obs_names_make_unique()
vasa_embryo.raw = vasa_embryo # Store the raw data

annotations = embryo.obs.query("sequencing_tech == 'vasaseq'").loc[:, ["C_scANVI"]].to_dict().get("C_scANVI")
vasa_embryo.obs["celltype"] = vasa_embryo.obs.index.map(annotations, na_action="ignore").fillna("Others") #[annotations.get(x, "Uknown") for x in vasa_embryo.obs.barcode]
vasa_embryo.obs["celltype_l0"] = vasa_embryo.obs.celltype.map(CELLTYPE_DICT, na_action="ignore").fillna("Others")
selected_cells = vasa_embryo.obs_names.drop_duplicates().to_list()
selected_genes = vasa_embryo.var_names.drop_duplicates().to_list()
vasa_embryo.obsm["X_pca_scVI"] = embryo[selected_cells, :].obsm["X_scVI"]
vasa_embryo.obsm["X_umap_scVI"] = embryo[selected_cells, :].obsm["X_umap"]

# velocity
scv.pp.filter_and_normalize(vasa_embryo, min_counts=3, min_counts_u=1, min_cells=3, min_cells_u=3)
scp.pp.highly_variable_genes(vasa_embryo, n_top_genes=2000)
scp.pp.pca(vasa_embryo, n_comps=50)
scp.pp.neighbors(vasa_embryo, n_pcs=50, n_neighbors=30)
scp.tl.umap(vasa_embryo)
scp.tl.leiden(vasa_embryo, flavor="igraph", n_iterations=10)
# scp.tl.rank_genes_groups(vasa_embryo, groupby="leiden", method="wilcoxon")

scv.pp.moments(vasa_embryo, n_pcs=None, n_neighbors=None)
scv.tl.velocity(vasa_embryo, mode="stochastic")
scv.tl.velocity_graph(vasa_embryo, n_jobs=10)
scv.tl.velocity_pseudotime(vasa_embryo)
scv.tl.velocity_embedding(vasa_embryo, basis="umap")

save_to = PROJECT_DIR / f"outputs/analysis/integration/velocity.colorby_celltype.png"
plots_velocity(vasa_embryo, figsize=(9, 8), color_by="celltype_l0", save_to=save_to)

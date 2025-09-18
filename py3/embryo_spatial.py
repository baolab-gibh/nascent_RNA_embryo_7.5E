#!/usr/bin/env python3
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)

from pathlib import Path
import yaml

import anndata as adt
import scanpy as scp
import squidpy as sqp

from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import polars as pls


plt.rcParams['figure.constrained_layout.use'] = True

CELLTYPE_DICT = {
    'ExE endoderm': "Endoderm", 'Def. endoderm': "Endoderm", 'Parietal endoderm': "Endoderm", 'Visceral endoderm': "Endoderm",
    'ExE mesoderm': "Mesoderm", 'Mixed mesoderm': "Mesoderm", 'Intermediate mesoderm': "Mesoderm", 'Pharyngeal mesoderm': "Mesoderm", 'Nascent mesoderm': "Mesoderm", 'Paraxial mesoderm': "Mesoderm", 'Somitic mesoderm': "Mesoderm", 'Caudal Mesoderm': "Mesoderm",
    'ExE ectoderm': "Ectoderm", 'Surface ectoderm': "Ectoderm", 'Rostral neurectoderm': "Ectoderm", 'Caudal neurectoderm': "Ectoderm",
    'Caudal epiblast': "Epiblast", 'Epiblast': "Epiblast",
    'Primitive Streak': "Primitive Streak", 'Anterior Primitive Streak': "Primitive Streak",
    # 'Haematoendothelial progenitors', 'Blood progenitors 2', 'Blood progenitors 1',
    # 'PGC', 'Gut', 'Notochord', 'Allantois', 'Mesenchyme',
}

proj_dir = Path('~/Documents/projects/wp_vasaseq').expanduser()
result_dir = proj_dir / 'outputs/analysis/preprocessing/slide_seq_decoded'

gtf_path = Path("~/Documents/projects/resources/references/gencode.vM29lift38.annotation.gtf").expanduser()
gene_len_list = (
    pls.read_csv(gtf_path, separator="\t", comment_prefix="#", has_header=False)
    .filter(pls.col("column_3") == "gene")
    .with_columns(pls.col("column_9").str.extract("gene_name \"(.*?)\";", 1).alias("gene_name"))
    .with_columns((pls.col("column_5") - pls.col("column_4")).alias("length"))
    .select("gene_name", "length")
)
gene_length_dict = dict(zip(gene_len_list["gene_name"], gene_len_list["length"]))


# Define celltype
in_file = proj_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/PijuanSala_etal_Nature_2019.processed.h5ad"
adata_pijuansala = scp.read_h5ad(in_file)
selected_cells = adata_pijuansala.obs.stage.isin(["E7.25", "E7.5"])
adata_pijuansala = adata_pijuansala[selected_cells, :].copy()
adata_pijuansala.obs["celltype_l0"] = adata_pijuansala.obs.celltype.map(CELLTYPE_DICT, na_action="ignore").fillna("Others")
adata_pijuansala.obs["Sequencing"] = "Single-cell"
scp.tl.dendrogram(adata_pijuansala, groupby="celltype", n_pcs=50)
scp.tl.embedding_density(adata_pijuansala, groupby="celltype_l0")


groups_to_show = ["Ectoderm", "Endoderm", "Mesoderm", "Epiblast", "Primitive Streak"]

# UMAP of original cell type by Pijuan-Sala et al. Nature 2019
fig = scp.pl.umap(adata_pijuansala, color="celltype", return_fig=True)
fig_saveto = proj_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/overview.E7x.umap.celltype_l1.png"
fig.set_size_inches(9, 4)
fig.savefig(fig_saveto)
plt.close()

# UMAP of cell type merged based on definitions by Pijuan-Sala et al., Nature 2019
fig = scp.pl.umap(adata_pijuansala, color="celltype_l0", return_fig=True)
fig_saveto = proj_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/overview.E7x.umap.celltype_l0.png"
fig.set_size_inches(5.5, 4)
fig.savefig(fig_saveto)
plt.close()

# UMAPs to show density of identifile cell types (L0)
fig = scp.pl.embedding_density(adata_pijuansala, key='umap_density_celltype_l0', group=groups_to_show, return_fig=True)
# fig_saveto = proj_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/overview.umap.embedding_density.E7x.pdf"
fig_saveto = proj_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/overview.E7x.umap.embedding_density.celltype_l0.png"
fig.set_size_inches(15, 7)
fig.savefig(fig_saveto)
plt.close()

# Identify markers genes per cell type
scp.tl.rank_genes_groups(adata_pijuansala, groupby="celltype_l0", key_added="by_celltype_l0")
deg_tab = pls.from_pandas(scp.get.rank_genes_groups_df(adata_pijuansala, group=groups_to_show, key="by_celltype_l0"))

# Matrixplot to show marker genes
marker_features = (
    deg_tab
    .filter(pls.col("pvals_adj") < 0.05, pls.col("logfoldchanges").abs() >= 2) 
    .with_columns(pls.col("logfoldchanges").abs().rank(descending=True).over("group").alias("rank")) 
    .filter(pls.col("rank") <= 5)
    .sort("group", pls.col("logfoldchanges").abs(), descending=True)
    .get_column("names")
    .to_list())

fig_size = (7, 12)
fig = scp.pl.rank_genes_groups_matrixplot(
    adata_pijuansala, groups=groups_to_show + ["Others"], n_genes=10, groupby="celltype_l0", key="by_celltype_l0",
    return_fig=True, dendrogram=True, swap_axes=True, figsize=fig_size, cmap="Spectral_r", edgecolors="none"
)
fig_saveto = proj_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/overview.E7x.matrixplot_by_celltype_l0.marker_genes.pdf"
fig.savefig(fig_saveto)
plt.close()


# Define marker genes from GEO-seq by Luyao Pan et al., X, 2025
adata_gvsseq = scp.read_h5ad(proj_dir / "outputs/analysis/velocity/version_3/dynamo.velocity_analysis.h5ad")

regions = adata_gvsseq.obs.Regions.unique().tolist()
cell_types = adata_gvsseq.obs.Cell_types.unique().tolist()
scp.tl.rank_genes_groups(adata_gvsseq, groupby="Regions", use_raw=False, pts=True, key_added="by_Regions")
scp.tl.rank_genes_groups(adata_gvsseq, groupby="Cell_types", use_raw=False, pts=True, key_added="by_Cell_types")
scp.tl.dendrogram(adata_gvsseq, groupby="Regions", n_pcs=50)


# Markter gene by logfoldchanges
fig_size = (7, 12)
fig = scp.pl.rank_genes_groups_matrixplot(
    adata_gvsseq, groups=regions, n_genes=10, min_logfoldchange=1, groupby="Regions", key="by_Regions",
    values_to_plot="logfoldchanges", dendrogram=True, vmin=-1.5, vmax=1.5,
    swap_axes=True, figsize=fig_size, cmap="Spectral_r", return_fig=True
)
fig_saveto = proj_dir / f"outputs/analysis/velocity/version_3/overview.all_regions.matrixplot_by_region.by_logfoldchanges.marker_genes.pdf"
fig.savefig(fig_saveto)
plt.close()


fig_size = (5, 12)
fig = scp.pl.rank_genes_groups_matrixplot(
    adata_gvsseq, groups=cell_types, n_genes=30, min_logfoldchange=1, groupby="Cell_types", key="by_Cell_types",
    values_to_plot="logfoldchanges", dendrogram=True, vmin=-3, vmax=3,
    swap_axes=True, figsize=fig_size, cmap="Spectral_r", return_fig=True
)
fig_saveto = proj_dir / f"outputs/analysis/velocity/version_3/overview.all_regions.matrixplot_by_cell_types.by_logfoldchanges.marker_genes.pdf"
fig.savefig(fig_saveto)
plt.close()

# Maker gene by mean expression
fig_size = (7, 12)
fig = scp.pl.rank_genes_groups_matrixplot(
    adata_gvsseq, groups=regions, n_genes=10, min_logfoldchange=1, groupby="Regions", key="by_Regions", dendrogram=True,
    swap_axes=True, figsize=fig_size, cmap="Spectral_r", return_fig=True
)
fig_saveto = proj_dir / f"outputs/analysis/velocity/version_3/overview.all_regions.matrixplot_by_region.by_mean_expression.marker_genes.pdf"
fig.savefig(fig_saveto)
plt.close()


fig_size = (5, 12)
fig = scp.pl.rank_genes_groups_matrixplot(
    adata_gvsseq, groups=cell_types, n_genes=30, min_logfoldchange=1, groupby="Cell_types", key="by_Cell_types", dendrogram=True,
    swap_axes=True, figsize=fig_size, cmap="Spectral_r", return_fig=True
)
fig_saveto = proj_dir / f"outputs/analysis/velocity/version_3/overview.all_regions.matrixplot_by_cell_types.by_mean_expression.marker_genes.pdf"
fig.savefig(fig_saveto)
plt.close()


# Integration between LuyaoPan_etal_x_2025 and PijuanSala_etal_Nature_2019
in_file = proj_dir / "outputs/analysis/velocity/version_3/all_batches.raw_data.h5ad"
adata_gvsseq = scp.read_h5ad(in_file)
adata_gvsseq.obs["Sequencing"] = "Geo-Vasa-Slam-seq"
adata_gvsseq.var["Length"] = adata_gvsseq.var.index.map(gene_length_dict).fillna(10000)
adata_gvsseq.raw = adata_gvsseq
scp.pp.normalize_total(adata_gvsseq, target_sum=1e6)
scp.pp.log1p(adata_gvsseq)

in_file = proj_dir / "outputs/analysis/development_stage/PijuanSala_etal_Nature_2019/PijuanSala_etal_Nature_2019.processed.h5ad"
adata_pijuansala = scp.read_h5ad(in_file)
selected_cells = adata_pijuansala.obs.stage.isin(["E7.5"])
adata_pijuansala = adata_pijuansala[selected_cells, :].copy()
adata_pijuansala.raw = adata_pijuansala
adata_pijuansala.obs["Cell_types"] = adata_pijuansala.obs.celltype.map(CELLTYPE_DICT, na_action="ignore").fillna("Others")
adata_pijuansala.obs["Sequencing"] = "Single-cell"
adata_pijuansala.obs["Batches"] = "PijuanSala_etal_Nature_2019"
adata_pijuansala.obs["Regions"] = "Single-cell"

adata_cmb = adt.concat([adata_pijuansala, adata_gvsseq], join="inner")
scp.pp.filter_cells(adata_cmb, min_genes=200)
scp.pp.filter_genes(adata_cmb, min_cells=3)
scp.pp.combat(adata_cmb, key="Sequencing")
scp.pp.highly_variable_genes(adata_cmb, batch_key="Sequencing")
scp.tl.pca(adata_cmb, n_comps=50, use_highly_variable=True, svd_solver="arpack")
scp.pp.neighbors(adata_cmb, n_pcs=50)
scp.tl.umap(adata_cmb)

regions = ["A", "L", "R", "P", "MA", "MP", "EA", "EP"]
regions_colors = dict(zip(regions, plt.cm.get_cmap("Set1").colors))
regions_colors.update({"Single-cell": "0.9"})

cell_types = ["Ectoderm", "Mesoderm", "Endoderm", "Epiblast", "Primitive Streak", "Others"]
cell_types_colors = dict(zip(cell_types, plt.cm.get_cmap("tab10").colors))

fig_saveto = proj_dir / "outputs/analysis/overview/overview.umap.PijuanSalaEtalNature2019_with_GeoVasaSlamSeq.pdf"
fig, (axe1, axe2) = plt.subplots(1, 2, figsize=(15, 6), tight_layout=True)
scp.pl.umap(adata_cmb, color="Cell_types", groups=cell_types, add_outline=True, size=150, palette=cell_types_colors, ax=axe1)
scp.pl.umap(adata_cmb, color="Regions", groups=regions + ["Single-cell"], add_outline=True, size=250, ax=axe2, palette=regions_colors)
fig.savefig(fig_saveto)
plt.close()


# Check the marker genes for Pijuan-Sala et al., Nature 2019
library_info = {
    "20250319": {
        "20250319Lib1": "OST110083", "20250319Lib2": "OST110084", "20250319Lib3": "OST110085",
        "20250319Lib4": "OST110086", "20250319Lib5": "OST110089"
    }
}

selected_features = {
    "PijuanSala_etal_Nature_2019": {
        "Ectoderm": ["Id2", "Gjb3", "Fabp3", "Elf5", "Gm1673", "Ddah1", "Cldn3", "S100a6", "Tfap2c", "T"],
        "Endoderm": ["Fth1", "Spink1", "Ctsh", "Cotl1", "Amn", "Apom", "Dab2", "Ttr", "Apoa1", "Pla2g12b"],
        "Epiblast": ["Pim2", "Pou5f1", "Slc7a3", "Dnmt3b", "Wfdc2", "Utf1", "Igfbp2", "Snrpn", "Tuba1a", "Gng3"],
        "PrimitiveStreak": ["Pou5f1", "Tdgf1", "Fgf8", "Igfbp2", "Snrpn", "Fst", "Tuba1a", "Npm3", "Eomes", "Pim2"],
        "Mesoderm": ["Ifitm1", "Mesp1", "Ccnd2", "Rbms1", "Fn1", "Hmgb1", "Fgf3", "Phlda2", "Dll3", "Pclaf"],
    },
    "LuyaoPan_etal_x_2025": {
        "Ectoderm": ["Sox3", "Igsf9b", "Akap12", "Pim2", "Map1b", "Nup210", "Jph4", "Pou5f1", "Slitrk5", "Nfasc", "Pdzd4", "Nrxn2"],
        "Endoderm": ["Slc2a3", "Emb", "Car4", "Krt8", "Krt18", "Slc16a1", "Trh", "Col4a1", "Tmprss2", "Col4a2", "Sox17", "Cpm", "Gsn"],
        "Mesoderm": ["Pdgfra", "Fn1", "Sh3pxd2a", "Foxc2", "Phlda2", "Tenm4", "Aplnr", "Lama1", "Pcdh19", "Dusp9", "Fbn2", "Adam19"],
    },
    "LuyaoPan_etal_x_2025.by_regions": {
        "A": ["Snrpn", "Pim2", "Rapgef5", "Plcxd1", "Sbk1", "Nrxn2", "Camkv", "Dnmt3b", "Slc7a3"], # "Hist3h2ba",
        "L": ["Cxcl12", "Cenpo", "Tbkbp1", "Syt11", "Alpl", "Mrps31", "Pdzd4", "Irx3", "Kif1a", "Cdkn1a"],
        "R": ["Cxcl12", "Sox3", "Prr5l", "Prr36", "Celsr3", "Irx3", "Klhl7", "Sema4d", "Slc35f2", "Usp44"],
        "P": ["Pou5f1", "Nkx1-2", "Fzd10", "Nefm", "Sema6a", "Apln", "Hoxb1", "Zfp428"], #"Cxx1a", "Cxx1b"
        "MA": ["Pdgfra", "Foxc2", "Fn1", "Lama1", "Dusp9", "Pkdcc", "Phlda2", "Epb41l3", "Fbn1", "Robo3"],
        "MP": ["Fn1", "Aplnr", "Pdgfra", "Pcdh19", "Igfbp4", "Sh3pxd2a", "Phlda2", "Ifitm1", "Adam19", "Cdc42ep4"],
        "EA": ["Slc2a3", "Emb", "Slc16a1", "Trh", "Col4a2", "Krt18", "Krt8", "Car4", "Cldn6", "Prss8"],
        "EP": ["Krt8", "Car4", "Slc39a8", "Slc2a3", "Tmprss2", "Krt18", "Emb", "Col4a1", "Trh", "Apoe"],
    },
}

overwrite = True
for bin_size in [10, 20, 50]:
    for batch_id, batch_info in library_info.items():
        for sample_id, chip_id in batch_info.items():
            # if sample_id != target_sample: continue

            try:
                in_file =  result_dir / f'{batch_id}_decoded_embryo/{sample_id}/07.analysis_bin{bin_size}/{chip_id}_bin{bin_size}.h5ad'
                adata_spatial = scp.read_h5ad(in_file)
            except Exception as e:
                adata_spatial = None
                print(e)
                continue

            for source, marker_features in selected_features.items():
                for celltype, feature_list in marker_features.items():
                    for per_feature in feature_list:
                        fig_saveto = proj_dir / 'outputs/analysis/spatial' / f'{batch_id}/{sample_id}/marker_genes/squidpy.spatial_map.bin{bin_size}.{source}.{celltype}.{per_feature}_expression.png'
                        if fig_saveto.exists() and not overwrite: continue
                        try:
                            fig_size = (10, 8)
                            fig, ((axe1, axe2), (axe3, axe4)) = plt.subplots(2, 2, figsize=fig_size, constrained_layout=True)
                            scp.pl.embedding(adata_spatial, basis="umap", color="cluster", ax=axe1)
                            scp.pl.embedding(adata_spatial, basis="umap", color=per_feature, ax=axe2)
                            scp.pl.embedding(adata_spatial, basis="spatial", color="cluster", ax=axe3)
                            scp.pl.embedding(adata_spatial, basis="spatial", color=per_feature, ax=axe4)
                            fig.savefig(fig_saveto)
                        except Exception as e:
                            print(e)
                        plt.close()
print("--- Done ---")


# Merging spatial and Geo-Vasa-Slam-seq
in_file = proj_dir / "outputs/analysis/velocity/version_3/all_batches.raw_data.h5ad"
adata_gvsseq = scp.read_h5ad(in_file)
adata_gvsseq.obs["Sequencing"] = "GEO-seq"
for bin_size in [10, 20, 50, 100]:
    for batch_id, batch_info in library_info.items():
        for sample_id, chip_id in batch_info.items():
            fig_saveto = proj_dir / f"outputs/analysis/spatial/{batch_id}/{sample_id}/overview.umap.normalized_fpkm.integrated_with_geo_seq.bin{bin_size}.pdf"

            #if fig_saveto.exists(): continue
            try:
                in_file =  result_dir / f'{batch_id}_decoded_embryo/{sample_id}/07.analysis_bin20/{chip_id}_bin{bin_size}.h5ad'
                if not in_file.exists():
                    print(f"{in_file} does not exist")
                    continue
                adata_spatial = scp.read_h5ad(in_file)
                adata_spatial.obs["Batches"] = sample_id
                adata_spatial.obs["Regions"] = "Unknown"
                adata_spatial.obs["Layers"] = "Unknown"
                adata_spatial.obs["Cell_types"] = "Unknown"
                adata_spatial.obs["Sequencing"] = "Slide-seq"
                adata_spatial.X = adata_spatial.layers["raw"]

                adata_cmb = adt.concat([adata_gvsseq, adata_spatial], join="inner")
                adata_cmb.raw = adata_cmb

                adata_cmb.var["Length"] = adata_cmb.var.index.map(gene_length_dict).fillna(10000)

                scp.pp.filter_cells(adata_cmb, min_genes=500)
                scp.pp.filter_genes(adata_cmb, min_cells=3)

                adata_cmb.X = adata_cmb.X / adata_cmb.var["Length"].to_numpy() * 1000000
                # (Number of Fragments / Gene Length (kb)) * 1000000
                # scp.pp.normalize_total(adata_cmb, target_sum=1e6)

                scp.pp.log1p(adata_cmb)
                scp.pp.combat(adata_cmb, key="Batches")
                scp.pp.highly_variable_genes(adata_cmb, n_top_genes=2000, batch_key="Batches")
                scp.tl.pca(adata_cmb, n_comps=50, use_highly_variable=True, svd_solver="arpack")
                scp.pp.neighbors(adata_cmb, n_pcs=10)
                scp.tl.umap(adata_cmb)

                fig, axe = plt.subplots(1, 1, figsize=(7, 5), tight_layout=True)
                scp.pl.umap(adata_cmb, color="Regions", ax=axe)
                fig.savefig(fig_saveto)
            except Exception as e:
                print(e)
            else:
                print(f"--- {sample_id}, {bin_size}. done ---")

            plt.close()
print("--- Done ---")



# Merge PijuanSala_etal_Nature_2019 and Slide-seq
for bin_size in [10, 20, 50, 100]:
    for batch_id, batch_info in library_info.items():
        for sample_id, chip_id in batch_info.items():
            fig_saveto = proj_dir / f"outputs/analysis/spatial/{batch_id}/{sample_id}/overview.umap.integrated_with_PijuanSala_etal.bin{bin_size}.pdf"

            # if fig_saveto.exists(): continue
            try:
                in_file =  result_dir / f'{batch_id}_decoded_embryo/{sample_id}/07.analysis_bin20/{chip_id}_bin{bin_size}.h5ad'
                if not in_file.exists():
                    print(f"{in_file} does not exist")
                    continue
                adata_spatial = scp.read_h5ad(in_file)
                adata_spatial.obs["sample"] = sample_id
                adata_spatial.obs["celltype_l0"] = "Unknown"
                adata_spatial.obs["Sequencing"] = "Slide-seq"
                adata_spatial.X = adata_spatial.layers["raw"]
                scp.pp.normalize_total(adata_spatial, target_sum=1e6)
                scp.pp.log1p(adata_spatial)

                adata_cmb = adt.concat([adata_pijuansala, adata_spatial], join="inner")
                adata_cmb.raw = adata_cmb

                scp.pp.filter_cells(adata_cmb, min_genes=500)
                scp.pp.filter_genes(adata_cmb, min_cells=3)
                scp.pp.combat(adata_cmb, key="Sequencing")
                scp.pp.highly_variable_genes(adata_cmb, n_top_genes=2000, batch_key="Sequencing")
                scp.tl.pca(adata_cmb, n_comps=50, use_highly_variable=True, svd_solver="arpack")
                scp.pp.neighbors(adata_cmb, n_pcs=50)
                scp.tl.umap(adata_cmb)

                fig, axe = plt.subplots(1, 1, figsize=(7, 5), tight_layout=True)
                scp.pl.umap(adata_cmb, color="celltype_l0", ax=axe)
                fig.savefig(fig_saveto)
                print(f"--- {sample_id}, {bin_size}. done ---")

                del adata_cmb, adata_spatial
            except Exception as e:
                print(f"--- {sample_id}, {bin_size}. failed ---", e)


            plt.close()
print("--- Done ---")


# library information

feature_list = {
    "A": [  "Snrpn",   "Pim2", "Rapgef5", "Plcxd1",   "Sbk1", "Nrxn2", "Camkv", "Dnmt3b", "Slc7a3", "Sox2ot"], #, "Hist3h2ba" 
    "P": [ "Pou5f1", "Nkx1-2",   "Fzd10",   "Nefm", "Sema6a",  "Apln", "Hoxb1", "Zfp428", "T", "Mixl1", "Mesp1", ], #, "Cxx1b"
}

library_info_path = proj_dir / "scripts/snakemake/slide_seq.configuration.yaml"
with open(library_info_path) as inhand:
    library_info = yaml.load(inhand, Loader=yaml.FullLoader)

obj_list = []
bin_size = 10
for batch_id, batch_info in library_info.get("sequencing_batches").items():
    for sample_id in batch_info.get('samples').keys():
        try:
            in_file =  result_dir / f'{batch_id}/{sample_id}/07.analysis_bin{bin_size}_with_image/{sample_id}_bin{bin_size}.h5ad'
            if not in_file.exists():
                print(f"{in_file} does not exist")
                continue
            adata_spatial = scp.read_h5ad(in_file)
            adata_spatial.obs["sample"] = sample_id
            adata_spatial.obs["celltype_l0"] = "Unknown"
            adata_spatial.obs["Sequencing"] = "Slide-seq"
            obj_list.append(adata_spatial)
            # del adata_cmb, adata_spatial
        except Exception as e:
            print(f"--- {batch_id}, {sample_id}, {bin_size}. failed due to {e} ---")

print("---done---")
embryo = adt.concat(obj_list, join="inner")
embryo.obs_names_make_unique()
embryo = embryo[embryo.obs.total_counts >= 2500, :]

embryo.obsm["PerX_pca"] = embryo.obsm["X_pca"]
embryo.obsm["PerX_umap"] = embryo.obsm["X_umap"]
embryo.obsm["PerX_tsne"] = embryo.obsm["X_tsne"]

scp.pp.filter_cells(embryo, min_genes=200)
scp.pp.filter_genes(embryo, min_cells=3)
scp.pp.highly_variable_genes(embryo, n_top_genes=2000)
scp.pp.pca(embryo, n_comps=50, use_highly_variable=True, svd_solver="arpack")
scp.pp.neighbors(embryo)
scp.tl.umap(embryo)

fig, axe = plt.subplots(1, 1, figsize=(7, 5), tight_layout=True)
scp.pl.embedding(embryo, basis="umap", color="sample", ax=axe)
fig.savefig(proj_dir / "outputs/analysis/spatial/overview.umap.20250423_embryo.pdf")

fig, axe = plt.subplots(1, 1, figsize=(7, 5), tight_layout=True)
scp.pl.embedding(embryo, basis="PerX_umap", color="sample", ax=axe)
fig.savefig(proj_dir / "outputs/analysis/spatial/overview.raw_umap.20250423_embryo.pdf")

fig, axe = plt.subplots(1, 1, figsize=(7, 5), tight_layout=True)
scp.pl.embedding(embryo, basis="spatial", color="sample", ax=axe)
fig.savefig(proj_dir / "outputs/analysis/spatial/overview.spatial.20250423_embryo.pdf")

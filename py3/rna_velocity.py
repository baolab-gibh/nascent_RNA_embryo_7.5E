#!/usr/bin/env python3
# File: rna_velocity.py
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Jul 01, 2024
# Updated: Aug 06, 2024

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

import os
import pickle
import itertools
from pathlib import Path

# matplotlib: version 3.5.3, matplotlib may raise errors due to depreciation of register_cmap()
import numpy as npy
import pandas as pds
import polars as pls
import anndata as adt
# import squidpy as sqp
import scanpy as scp
import dynamo as dyn
import scvelo as scv
import seaborn as sbn
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl

from dynamo.preprocessing import Preprocessor
from scipy import sparse
from scipy import cluster
from scipy import spatial
from sklearn.decomposition import PCA


# Global settings
# mpl.rcParams["legend.loc"] = "lower right"
# mpl.rcParams["legend.borderpad"] = 0.1
# mpl.rcParams["legend.frameon"] = False
# mpl.rcParams["legend.scatterpoints"] = 5
# 
# dyn.configuration.set_figure_params(True, background='black')
# dyn.configuration.set_pub_style()


VERSION = "version_3" # Including five batches, 240409_Lib_embryo, 240612_Lib_28region, 240620_Lib_38region, 240703_Lib_32region, 240710_Lib_37region
REGION_ORDER = ["EA", "MA", "A", "L", "R", "P", "MP", "EP"]
PROJECT_DIR = Path("~/Documents/projects/wp_vasaseq").expanduser()
ALL_BATCHES = [
    "240409_Lib_embryo", "240612_Lib_28region", "240620_Lib_38region", "240703_Lib_32region",
    "240710_Lib_37region", "240717_Lib_28region"
]


def convert_pos(adt, into="cell_type"):
    '''Convert sample name into cell type, regions, or layers.'''
    _pos = ['A', 'P', 'EA', 'EP', 'L', 'R', 'MA', 'MP', 'O']
    if into == "regions":
        return adt.obs.index.str.extract("[0-9]+([A-Z]{1,2})").loc[:, 0].tolist()
    elif into == "layers":
        return adt.obs.index.str.extract("(NC_[0-9]+|[0-9]+)[A-Z]{1,2}").loc[:, 0].tolist()
    else:
        if into == "pseudotime":
            #_target = [3, 3, 1, 1, 3, 3, 2, 2, 0]
            _target = [1, 1, 1, 1, 1, 1, 1, 1, 0]
        else:
            _target = ["Ectoderm", "Ectoderm", "Endoderm", "Endoderm", "Ectoderm", "Ectoderm", "Mesoderm", "Mesoderm", "Other"]

        pseudo_time = dict(zip(_pos, _target))
        regions = [x[0] if x else "O" for x in adt.obs.index.str.findall("[0-9]+([A-Z]{1,2})")]
        return [pseudo_time[x] for x in regions]


def add_obs_ntr(adata, base_dir: Path):
    '''Add new-to-total ratio (NTR) to `adata.obs`.'''
    all_ntrfile_list = {x: (base_dir / x) / (x + ".nascent_rna.ntrstat.tsv") for x in ALL_BATCHES}

    all_ntrtab_list = []
    for key, val in all_ntrfile_list.items():
        per_ntr_tab = pls.read_csv(val, separator="\t").with_columns(pls.lit(key).alias("Batches"))
        all_ntrtab_list.append(per_ntr_tab)

    required_cols = ["Sample", "T->C", "ntr_lower", "ntr", "ntr_upper"]
    ntrtab = (
        pls.concat(all_ntrtab_list, how="vertical")
        .with_columns(
            pls.col("Condition").map_elements(
                lambda x: x.replace("L_", "R_") if "L_" in x else (x.replace("R_", "L_") if "R_" in x else x)))
        .with_columns((pls.col("Condition") + "-" + pls.col("Batches")).alias("Sample"))
        .select("Sample", "Type", "Condition", "T->C", "ntr_lower", "ntr", "ntr_upper")
    )

    ntrtab_list = [adata.obs.copy()]
    for per_type in ["Exonic", "Intronic"]:
        per_tab = (
            ntrtab.filter(pls.col("Type") == per_type, ~pls.col("Condition").str.contains("NC_"))
            .select(required_cols)
            .to_pandas()
            .set_index("Sample")
            .rename(columns=lambda x: per_type + "_" + x)
        )
        ntrtab_list.append(per_tab)

    return pds.concat(ntrtab_list, axis=1).copy()


def add_x_ntr(adata: scp.AnnData, base_dir: Path, min_rc: int = 10):
    """Add NTR of each gene of each sample for `AnnData.X`."""
    base_dir = PROJECT_DIR / "outputs/analysis/preprocessing/geo_vasa_slam/nascent_rna"
    all_ntrfile_list = {x: (base_dir / x) / (x + ".nascent_rna.tsv.gz") for x in ALL_BATCHES}

    all_ntrtab_list = []
    required_cols = ["Symbol", "Sample", "NTR"]
    selected_features = adata.var_names.unique().tolist()
    for key, val in all_ntrfile_list.items():
        per_ntr_tab = (
            pls.read_csv(val, separator="\t")
            .filter(pls.col("Symbol").is_in(selected_features))
            .select("Symbol", pls.col("^.* MAP|alpha|beta|Readcount$"))
            .unpivot(index="Symbol", value_name="Values", variable_name="Sample")
            .filter(~pls.col("Sample").str.starts_with("NC_"))
            .with_columns(pls.col("Sample").str.extract("^(.*) (MAP|Readcount|alpha|beta)$", 2).alias("Type"))
            .with_columns(pls.col("Sample").str.extract("^(.*) (MAP|Readcount|alpha|beta)$", 1) + "-" + key)
            .fill_nan(0)
            .pivot("Type", index=["Symbol", "Sample"], values="Values", aggregate_function="mean")
            .with_columns(
                Kept = (pls.col("Readcount") > min_rc) &
                ((0.5 <= pls.col("alpha")) & (pls.col("alpha") <= 1500)) &
                ((0.5 <= pls.col("beta")) & (pls.col("beta") <= 1500)))
            .with_columns(NTR = pls.when(pls.col("Kept")).then(pls.col("MAP")).otherwise(0.0))
            .select(required_cols)
        )

        all_ntrtab_list.append(per_ntr_tab)
    
    return sparse.csr_matrix(
        pls.concat(all_ntrtab_list, how="vertical")
        .with_columns(
            pls.col("Sample").map_elements(
                lambda x: x.replace("L_", "R_") if "L_" in x else (x.replace("R_", "L_") if "R_" in x else x)
            ))
        .pivot("Symbol", index="Sample", values="NTR", aggregate_function="mean")
        .select("Sample", *selected_features).fill_null(0).fill_nan(0)
        .to_pandas().set_index("Sample").loc[adata.obs_names, adata.var_names].to_numpy()
    )


def add_var_ntr(adata, min_tc=10, min_nc=3, max_tc=5000, max_nc=5000, use_bf=True, fillnan=0):
    '''Add NTR to `adata.var`.'''
    has_total, has_new = "total" in adata.layers, "new" in adata.layers
    if not (has_new and has_total): raise KeyError("Missing 'new' or 'total' key in adata.layers data.")

    tobe_incl = True
    if min_tc is not None and max_tc is not None: # Check min total counts (min_tc) and max total counts (max_tc)
        tobe_incl *= (min_tc <= adata.layers["total"]).A * (adata.layers["total"] <= max_tc).A
    if min_nc is not None and max_nc is not None: # Check min new counts (min_nc) and max new counts (max_nc)
        tobe_incl *= (min_nc <= adata.layers["new"]).A * (adata.layers["new"] <= max_nc).A
    if "pass_basic_filter" in adata.var and use_bf: # Using pass_basic_filter for variables (i.e., genes)
        tobe_incl *= adata.var["pass_basic_filter"].values
    if "pass_basic_filter" in adata.obs and use_bf: # Using pass_basic_filter for observations (i.e., samples)
        tobe_incl *= adata.obs["pass_basic_filter"].values

    ntr_mat = ((adata.layers["new"].A * tobe_incl).sum(0) / (adata.layers["total"].A * tobe_incl).sum(0))
    if fillnan is not None and isinstance(fillnan, (int, float)):
        ntr_mat = npy.nan_to_num(ntr_mat, nan=fillnan)

    return adata.var.copy().assign(NTR=ntr_mat)


def init_adata(ttl_adt, new_adt=None, new_as_base=False, ntr_base_dir=None):
    '''Create an AnnData from new and total reads AnnData'''
    adata = new_adt.copy() if new_as_base and new_adt is not None else ttl_adt.copy()
    if new_adt is not None:
        adata.layers["new"] = new_adt.X.copy()
        adata.layers["total"] = ttl_adt.X.copy()
    adata.layers["counts"] = adata.X.copy()
    adata.obs["Groups"] = convert_pos(adata, into="pseudotime")
    adata.obs["Regions"] = convert_pos(adata, into="regions")
    adata.obs["Layers"] = convert_pos(adata, into="layers")
    adata.obs["Cell_types"] = convert_pos(adata)
    adata.obs["Layer_region"] = adata.obs.index.str.extract("([0-9]+[A-Z]{1,2})").loc[:, 0].tolist()
    adata.obs["Sampling_dates"] = adata.obs.Batches.str.extract("^([0-9]+)_").loc[:, 0].tolist()

    if ntr_base_dir is not None:
        adata.layers["NTR"] = add_x_ntr(adata, ntr_base_dir)
        adata.obs = add_obs_ntr(adata, ntr_base_dir)
        adata.var = add_var_ntr(adata)

    return adata


def fetch_group_degs(adata, group_by, top_by="logfoldchanges", min_logfc=1, max_pval_adj=0.05, **kwargs) -> pls.DataFrame:
    '''Fetch marker genes for each comparison'''
    group_name = adata.obs[group_by].unique().tolist()
    deg_tab_list = []
    for qry_group, ref_group in itertools.combinations(group_name, 2):
        tmp_adata = scp.tl.rank_genes_groups(adata, groupby=group_by, groups=[qry_group], reference=ref_group, copy=True, **kwargs)
        deg_tab = (pls.DataFrame(scp.get.rank_genes_groups_df(tmp_adata, qry_group, key="rank_genes_groups"))
                   .with_columns(pls.lit(qry_group).alias("query_group"), pls.lit(ref_group).alias("reference_group")))

        deg_tab_list.append(deg_tab)

    deg_tab_all = (pls.concat(deg_tab_list, how="vertical")
                   .filter(pls.col("pvals_adj") < max_pval_adj, pls.col("logfoldchanges").abs() > min_logfc)
                   .with_columns(pls.col(top_by).abs().neg().rank().over("query_group", "reference_group").alias("rank"))
                   .sort("reference_group", "rank", maintain_order=True))

    return deg_tab_all


def plot_feature_in_corn(ttl_adata, new_adata, feature, fig_size=(6, 3), save_to=None, **kwargs):
    '''Plot feature expression in corn plot, new vs total'''
    if "s" not in kwargs: kwargs["s"] = 100
    if "palette" not in kwargs: kwargs["palette"] = "Blues"

    # Basic data
    meta_data = pls.DataFrame(ttl_adata.obs).with_columns(pls.Series(ttl_adata.obs_names))
    pos_tab = pls.DataFrame(ttl_adata.obsm["X_corn"]).rename({"column_0": "layer_axis", "column_1": "region_axis"})
    
    # Total expression
    ttl_exp_tab = pls.DataFrame(ttl_adata[:, feature].X.toarray()).rename({"column_0": feature})
    ttl_plot_tab = (pls.concat([meta_data, pos_tab, ttl_exp_tab], how="horizontal").group_by("layer_axis", "region_axis")
                    .agg(pls.count().alias("N_samples"), pls.col(feature).min().alias(feature)))

    # Nascent expression
    new_exp_tab = pls.DataFrame(new_adata[:, feature].X.toarray()).rename({"column_0": feature})
    new_plot_tab = (pls.concat([meta_data, pos_tab, new_exp_tab], how="horizontal").group_by("layer_axis", "region_axis")
                    .agg(pls.count().alias("N_samples"), pls.col(feature).min().alias(feature)))

    # New to total ratio
    ntr_tab = pls.DataFrame(ttl_adata[:, feature].layers["NTR"].toarray()).rename({"column_0": feature})
    ntr_plot_tab = (pls.concat([meta_data, pos_tab, ntr_tab], how="horizontal").group_by("layer_axis", "region_axis")
                    .agg(pls.count().alias("N_samples"), pls.col(feature).min().alias(feature)))

    max_exp = max(ttl_plot_tab[feature].max(), new_plot_tab[feature].max())
    min_exp = min(ttl_plot_tab[feature].min(), new_plot_tab[feature].min())

    # New to total expression
    meta_tab = pls.DataFrame(ttl_adata.obs).select("Sampling_dates", "Regions", "Layers", "Cell_types")
    nt_expr_tab = pls.DataFrame({ f"ttl_{feature}": ttl_exp_tab[feature], f"new_{feature}": new_exp_tab[feature] })
    expr_plot_tab = pls.concat([meta_tab, nt_expr_tab], how="horizontal")

    # Plots
    with mpl.rc_context({"legend.markerscale": 0.5, "legend.frameon": True, "legend.handletextpad": 1.0, "legend.borderpad": 1.0}):
        fig = plt.figure(figsize=fig_size, tight_layout=True, layout="constrained")
        grid_spec = fig.add_gridspec(ncols=3, nrows=5)
        axe1 = fig.add_subplot(grid_spec[:3, :])
        axe2 = fig.add_subplot(grid_spec[3:, 0])
        axe3 = fig.add_subplot(grid_spec[3:, 1], sharex=axe2)
        axe4 = fig.add_subplot(grid_spec[3:, 2], sharex=axe3)

        _ = sbn.scatterplot(expr_plot_tab, x=f"ttl_{feature}", y=f"new_{feature}", style="Cell_types", hue="Regions", ax=axe1, edgecolor=None)
        _ = axe1.axline((0, 0), slope=1, color="k", linestyle="dashed", linewidth=0.25)
        _ = sbn.scatterplot(ttl_plot_tab, x="region_axis", y="layer_axis", hue=feature, hue_norm=(min_exp, max_exp), ax=axe2, **kwargs)
        _ = sbn.scatterplot(new_plot_tab, x="region_axis", y="layer_axis", hue=feature, hue_norm=(min_exp, max_exp), ax=axe3, **kwargs)
        _ = sbn.scatterplot(ntr_plot_tab, x="region_axis", y="layer_axis", hue=feature, hue_norm=(0, 1), ax=axe4, **kwargs)

        # x tick labels
        x_tick_pos = [-10.490, -7.824, -4.863, -1.605, 1.605, 4.863, 7.824, 10.490]
        x_tick_lab = ["EA", "MA", "A", "L", "R", "P", "MP", "EP"]
        ta_iters = zip(["Total vs nascent RNA", "Total RNA", "Nascent RNA", "New-to-total ratio"], [axe1, axe2, axe3, axe4])
        for title, per_axe in ta_iters:
            if title == "Total vs nascent RNA":
                per_axe.set_xlabel("Total RNA expression")
                per_axe.set_ylabel("Nascent RNA expression")
                per_axe.axis("equal")
                per_axe.set_aspect("equal")
                per_axe.spines["right"].set_visible(False)
                per_axe.spines["top"].set_visible(False)
            else:
                if title == "Total RNA":
                    per_axe.set_ylabel("Layer")
                    per_axe.set_yticks(list(range(1, 17)))
                else:
                    per_axe.set_ylabel(None)
                    per_axe.set_yticks([])

                per_axe.xaxis.tick_top()
                per_axe.set_xticks(x_tick_pos, x_tick_lab)
                per_axe.set_xlabel(title)
                per_axe.margins(0.075, 0.05)
            per_axe.set_title(None)
            per_axe.legend(title=feature, loc="upper left", bbox_to_anchor=(1.0, 1.0)) # Adjust legend

    if save_to is None: return fig
    fig.savefig(save_to)
    fig.clear()
    plt.close(fig)


def project_exp_into_corn(adata, feature, assay="X", layer_axis="Layers", region_axis="Regions", null_to=0.0):
    """Project the expression of given feature into corn coordinations."""
    global REGION_ORDER
    if isinstance(feature, str):
        feature = [feature]
    elif not isinstance(feature, list):
        raise TypeError("`feature` must be str or list")

    meta_cols = ["cell_barcode", layer_axis, region_axis]
    meta_tab = pls.DataFrame(adata.obs.copy()).with_columns(pls.Series(adata.obs_names).alias("cell_barcode")).select(meta_cols)

    if assay == "X":
        val_tab = adata[:, feature].X
    elif assay in ["new", "total", "NTR"]:
        val_tab = adata[:, feature].layers[assay].toarray()
    else:
        raise ValueError(f"Unknown assay: {assay}")
    col_names = {f"column_{i}": g for i, g in enumerate(feature)}
    val_tab = pls.DataFrame(val_tab).rename(col_names)
    corn_coords = pls.DataFrame(adata.obsm["X_corn"]).rename({"column_0": "layer_axis", "column_1": "region_axis"})
    obs_tab = pls.concat([meta_tab, val_tab, corn_coords], how="horizontal").select(meta_cols + ["layer_axis", "region_axis"] + feature)

    return (
        obs_tab.pivot(region_axis, index=layer_axis, values=feature, aggregate_function="mean")
        .sort(layer_axis, descending=False).select(REGION_ORDER).fill_null(null_to)
    )


def load_data():
    global PROJECT_DIR, ALL_BATCHES, REGION_ORDER
    gene_bkl = ["Gm42418", "Gm26917"]

    ttl_adt_list, new_adt_list, spl_adt_list = {}, {}, {}
    for per_batch in ALL_BATCHES:
        ten_x_dir = PROJECT_DIR / "outputs/analysis/preprocessing/geo_vasa_slam/quantification" / per_batch / "slamseq/10X"
        new_ad, ttl_ad = scp.read_10x_mtx(ten_x_dir / "new"), scp.read_10x_mtx(ten_x_dir / "total")
        new_ad = new_ad[:, new_ad.var.gene_ids.drop_duplicates()]
        new_ad.obs_names = new_ad.obs_names.to_series().str.replace("([0-9]+)L", "\\1Z", regex=True).str.replace("([0-9]+)R", "\\1L", regex=True).str.replace("([0-9]+)Z", "\\1R", regex=True).add("-" + per_batch)
        ttl_ad = ttl_ad[:, ttl_ad.var.gene_ids.drop_duplicates()]
        ttl_ad.obs_names = ttl_ad.obs_names.to_series().str.replace("([0-9]+)L", "\\1Z", regex=True).str.replace("([0-9]+)R", "\\1L", regex=True).str.replace("([0-9]+)Z", "\\1R", regex=True).add("-" + per_batch)
        new_ad.obs["Batches"], ttl_ad.obs["Batches"] = per_batch, per_batch
        ttl_ad.obs["Sampling_dates"] = ttl_ad.obs.Batches.str.extract("^([0-9]+)_").loc[:, 0].tolist()
        new_ad.obs["Sampling_dates"] = new_ad.obs.Batches.str.extract("^([0-9]+)_").loc[:, 0].tolist()

        ttl_adt_list[per_batch] = ttl_ad
        new_adt_list[per_batch] = new_ad

        spl_dir = PROJECT_DIR / "outputs/analysis/preprocessing/geo_vasa_slam/velocity_counts" / per_batch / "velocyto"
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

    ntr_base_dir = PROJECT_DIR / "outputs/analysis/preprocessing/geo_vasa_slam/nascent_rna"
    labeled_adata = init_adata(ttl_adata, new_adata, ntr_base_dir=ntr_base_dir)

    # Sliced data
    spl_adata = adt.concat(spl_adt_list, merge="same")

    # Add corn plot coordinates
    corn_coords = pds.read_csv(PROJECT_DIR / "inputs/reference/corn_coordination/corn_axis.e_7_5.csv")
    corn_coords["Layer_region"] = corn_coords.apply(lambda x: str(x["y_pos"]) + x["Regions_l1"], axis=1)
    corn_coords["x_pos"] = corn_coords["x_pos"]
    labeled_adata.obsm["X_corn"] = pds.merge(labeled_adata.obs.copy().reset_index(), corn_coords, how="left", on="Layer_region").set_index("index").loc[:, ["y_pos", "x_pos"]].to_numpy()

    return ttl_adata, new_adata, labeled_adata, spl_adata


def umap_xy_lims(adata: scp.AnnData, ext_ratio: float = 0.1):
    if "X_umap" in adata.obsm:
        x_lim_min, x_lim_max = adata.obsm["X_umap"][:, 0].min(), labeled_adata.obsm["X_umap"][:, 0].max()
        y_lim_min, y_lim_max = adata.obsm["X_umap"][:, 1].min(), labeled_adata.obsm["X_umap"][:, 1].max()
    else:
        return (None, None), (None, None)

    x_lim_span = (x_lim_max - x_lim_min) * ext_ratio
    x_lim_min, x_lim_max = x_lim_min - x_lim_span, x_lim_max + x_lim_span
    y_lim_span = (y_lim_max - y_lim_min) * ext_ratio
    y_lim_min, y_lim_max = y_lim_min - y_lim_span, y_lim_max + y_lim_span
    return (x_lim_min, x_lim_max), (y_lim_min, y_lim_max)


print("--- Define functions done ---")


#
## Analysis, velocity
#
# Load all data. ttl_adata for total read counts, new_adata for nascent read counts, labeled_adata for both.
labeled_adata_path = PROJECT_DIR / "outputs/analysis/geo_vasa_slam/velocity" / VERSION / "all_batches.raw_data.h5ad"
if labeled_adata_path.exists():
    labeled_adata = scp.read_h5ad(labeled_adata_path)
else:
    ttl_adata, new_adata, labeled_adata, spl_adata = load_data()
    labeled_adata.write_h5ad(labeled_adata_path)
print("--- loading data done ---")

# Corn map to check duplicated samples. Supplementary Figure xx
corn_coord_save_to = PROJECT_DIR / "outputs/analysis/overview/corn_map.samples_per_region.csv"
if not corn_coord_save_to.exists():
    corn_plot_tab = (
        pls.DataFrame(labeled_adata.obsm["X_corn"])
        .rename({"column_0": "layer_axis", "column_1": "region_axis"})
        .with_columns(pls.Series(labeled_adata.obs.index).alias("Samples"))
        .with_columns(pls.Series(labeled_adata.obs.Batches).alias("Batches"))
        .with_columns(pls.Series(labeled_adata.obs.Regions).alias("Regions"))
        .with_columns(pls.Series(labeled_adata.obs.Layer_region).alias("Layer_region"))
        .group_by("Layer_region")
        .agg(pls.count().alias("N_samples"), pls.col("layer_axis").mean().alias("layer_axis"), pls.col("region_axis").mean().alias("region_axis"))
        .with_columns(pls.col("Layer_region").str.extract("(^[0-9]+)").alias("Layer"), pls.col("Layer_region").str.extract("([A-Z]+$)").alias("Region")))
    corn_plot_tab.write_csv()
# Check r script for plotting
print("--- Samples per corn postion: done ---")


# Using scanpy to check transcriptional profile using total RNAs
adata_dict = {}
ntr_base_dir = PROJECT_DIR / "outputs/analysis/preprocessing/geo_vasa_slam/nascent_rna"
for per_cat in ["new", "total"]:
    new_as_base = per_cat == "new"
    tmp_adata = init_adata(ttl_adata, new_adata, new_as_base, ntr_base_dir)
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
    scp.tl.rank_genes_groups(tmp_adata, groupby="Cell_types", method="wilcoxon", key_added="One_vs_rest")

    adata_dict[per_cat] = tmp_adata

total_adata, nascent_adata = adata_dict["total"], adata_dict["new"]
total_adata.obsm["X_corn"], nascent_adata.obsm["X_corn"] = labeled_adata.obsm["X_corn"].copy(), labeled_adata.obsm["X_corn"].copy()


# UMAP by total and nascent RNA
umap_axe_kwargs = dict(frameon=False, edgecolor="0.5", legend_loc="on data", palette="Set2", show=False, size=300, color="Cell_types")
matrixplot_axe_kwargs = dict(var_group_rotation=0.5, show=False, n_genes=15, cmap="Spectral_r")

fig = plt.figure(figsize=(8, 10), layout="constrained", tight_layout=True)
axe1, axe2, axe3, axe4 = plt.subplot(421), plt.subplot(422), plt.subplot(412), plt.subplot(413)
with mpl.rc_context({"font.size": 8}):
    _ = scp.pl.umap(total_adata, title="UMAP by total RNA", ax=axe1, **umap_axe_kwargs)
    _ = scp.pl.umap(nascent_adata, title="UMAP by new RNA", ax=axe2, **umap_axe_kwargs)
    _ = scp.pl.rank_genes_groups_matrixplot(total_adata, groupby="Cell_types", key="One_vs_rest", ax=axe3, **matrixplot_axe_kwargs)
    _ = scp.pl.rank_genes_groups_matrixplot(nascent_adata, groupby="Cell_types", key="One_vs_rest", ax=axe4, **matrixplot_axe_kwargs)
plt.subplots_adjust(hspace=0.2, wspace=0.1, left=0.05, right=0.95, bottom=0.1, top=0.9)
fig.savefig(PROJECT_DIR / "outputs/analysis/overview/scanpy.total_and_new.umap_and_mark_gene.pdf")
fig.clear()
plt.close(fig)

umap_axe_kwargs = dict(frameon=False, edgecolor="0.5", legend_loc="on data", palette="Set2", show=False, size=300, color="Regions")
matrixplot_axe_kwargs = dict(var_group_rotation=0.5, show=False, n_genes=15, cmap="Spectral_r")
fig = plt.figure(figsize=(8, 10), layout="constrained", tight_layout=True)
axe1, axe2, axe3, axe4 = plt.subplot(421), plt.subplot(422), plt.subplot(412), plt.subplot(413)
with mpl.rc_context({"font.size": 8}):
    _ = scp.pl.umap(total_adata, title="UMAP by total RNA", ax=axe1, **umap_axe_kwargs)
    _ = scp.pl.umap(nascent_adata, title="UMAP by new RNA", ax=axe2, **umap_axe_kwargs)
    _ = scp.pl.rank_genes_groups_matrixplot(total_adata, groupby="Cell_types", key="One_vs_rest", ax=axe3, **matrixplot_axe_kwargs)
    _ = scp.pl.rank_genes_groups_matrixplot(nascent_adata, groupby="Cell_types", key="One_vs_rest", ax=axe4, **matrixplot_axe_kwargs)
plt.subplots_adjust(hspace=0.2, wspace=0.1, left=0.05, right=0.95, bottom=0.1, top=0.9)
fig.savefig(PROJECT_DIR / "outputs/analysis/overview/scanpy.total_and_new.umap_and_mark_gene.by_regions.pdf")
fig.clear()
plt.close(fig)


# Violin plots to show DEGs
ttl_deg_tab = fetch_group_degs(total_adata, "Cell_types", top_by="scores").with_columns((pls.col("query_group") + ".vs." + pls.col("reference_group")).alias("comparison"))
total_markers = ttl_deg_tab.filter(pls.col("names").str.contains("Rik$|^Gm").not_()).top_k(20, by="rank", reverse=True)["names"].unique().to_list()
gg = scp.pl.matrixplot(total_adata, total_markers, groupby='Regions', figsize=(4, 5), dendrogram=True, return_fig=True, cmap="Spectral_r", swap_axes=True)
_ = gg.add_totals().style(edge_color='black')
save_to = PROJECT_DIR / "outputs/analysis/overview/scanpy.total.degs_by_total_reads.heatmap.pdf"
gg.savefig(save_to, bbox_inches="tight")


# ------
fig, ((axe1, axe2, axe3), (axe4, axe5, axe6)) = plt.subplots(2, 3, figsize=(16 * 0.75, 8 * 0.75), sharey=True, layout="constrained")
with mpl.rc_context({"font.size": 8}):
    _ = scp.pl.rank_genes_groups_violin(total_adata, groups='Ectoderm', gene_names=total_markers, key="One_vs_rest", ax=axe1)
    _ = scp.pl.rank_genes_groups_violin(total_adata, groups='Endoderm', gene_names=total_markers, key="One_vs_rest", ax=axe2)
    _ = scp.pl.rank_genes_groups_violin(total_adata, groups='Mesoderm', gene_names=total_markers, key="One_vs_rest", ax=axe3)
    _ = scp.pl.rank_genes_groups_violin(nascent_adata, groups='Ectoderm', gene_names=total_markers, key="One_vs_rest", ax=axe4)
    _ = scp.pl.rank_genes_groups_violin(nascent_adata, groups='Endoderm', gene_names=total_markers, key="One_vs_rest", ax=axe5)
    _ = scp.pl.rank_genes_groups_violin(nascent_adata, groups='Mesoderm', gene_names=total_markers, key="One_vs_rest", ax=axe6)
    for idx, per_axe in enumerate([axe1, axe2, axe3, axe4, axe5, axe6]):
        per_axe.set_xlabel("")
        if idx > 2:
            per_axe.set_ylabel("Expression (Nascent)")
        else:
            per_axe.set_ylabel("Expression (Total)")
plt.subplots_adjust(hspace=0.4, wspace=0.05, left=0.05, right=0.95, bottom=0.1, top=0.9)
fig.savefig(PROJECT_DIR / "outputs/analysis/overview/scanpy.total_marker_genes.violin_plot.pdf")
fig.clear()
plt.close(fig)


# Violin plots to show DEGs by nascent data
new_deg_tab = fetch_group_degs(nascent_adata, "Cell_types", top_by="scores").with_columns((pls.col("query_group") + ".vs." + pls.col("reference_group")).alias("comparison"))
new_markers = new_deg_tab.filter(pls.col("names").str.contains("Rik$|^Gm").not_()).top_k(20, by="rank", reverse=True)["names"].unique().to_list()
gg = scp.pl.matrixplot(nascent_adata, new_markers, groupby='Regions', figsize=(4, 5), dendrogram=True, return_fig=True, cmap="Spectral_r", swap_axes=True)
_ = gg.add_totals().style(edge_color='black')
save_to = PROJECT_DIR / "outputs/analysis/overview/scanpy.total.degs_by_new_reads.heatmap.pdf"
gg.savefig(save_to, bbox_inches="tight")


fig, ((axe1, axe2, axe3), (axe4, axe5, axe6)) = plt.subplots(2, 3, figsize=(16 * 0.75, 8 * 0.75), sharey=True, layout="constrained")
with mpl.rc_context({"font.size": 8}):
    _ = scp.pl.rank_genes_groups_violin(total_adata, groups='Ectoderm', gene_names=new_markers, key="One_vs_rest", ax=axe1)
    _ = scp.pl.rank_genes_groups_violin(total_adata, groups='Endoderm', gene_names=new_markers, key="One_vs_rest", ax=axe2)
    _ = scp.pl.rank_genes_groups_violin(total_adata, groups='Mesoderm', gene_names=new_markers, key="One_vs_rest", ax=axe3)
    _ = scp.pl.rank_genes_groups_violin(nascent_adata, groups='Ectoderm', gene_names=new_markers, key="One_vs_rest", ax=axe4)
    _ = scp.pl.rank_genes_groups_violin(nascent_adata, groups='Endoderm', gene_names=new_markers, key="One_vs_rest", ax=axe5)
    _ = scp.pl.rank_genes_groups_violin(nascent_adata, groups='Mesoderm', gene_names=new_markers, key="One_vs_rest", ax=axe6)
    for idx, per_axe in enumerate([axe1, axe2, axe3, axe4, axe5, axe6]):
        per_axe.set_xlabel("")
        if idx > 2:
            per_axe.set_ylabel("Expression (Nascent)")
        else:
            per_axe.set_ylabel("Expression (Total)")
plt.subplots_adjust(hspace=0.4, wspace=0.05, left=0.05, right=0.95, bottom=0.1, top=0.9)
fig.savefig(PROJECT_DIR / "outputs/analysis/overview/scanpy.new_marker_genes.violin_plot.pdf")
fig.clear()
plt.close(fig)


# Show the expression of given gene in corn plot
overwrite = True
marker_features = ["Foxa1", "Dnmt3b", "Neat1", "Sox17", "Twist1", "T", "Gata6", "Gata3", "Tenm4", "Dlc1", "Smarcd3", "Foxc1", "Meis1", "Sox17", "Sox2", "Cpt1b"]
marker_features = ["Wnt3", "Nodal", "Lefty2", "Bmp4", "Axin2", "Fgfr1", "Fgf4", "Fgf8", "Tbx6", "Dll1", "Sox2", "Mesp1", "Mesp2", "Net1", "Rhoa", "T"]
for per_feature in marker_features:
    save_to = PROJECT_DIR / "outputs/analysis/overview" / f"scanpy.total_and_new.corn_plot.{per_feature}.pdf"
    if not save_to.exists() or overwrite:
        plot_feature_in_corn(total_adata, nascent_adata, per_feature, fig_size=(7, 7), save_to=save_to, linewidth=0, palette="Reds")
print("--- Done ---")


# Estimate the ntr ratio
ntr_adata = init_adata(ttl_adata, new_adata)
dyn.pp.filter_genes_by_outliers(ntr_adata, min_cell_s=3, min_cell_u=3, min_count_s=10, min_count_u=3)
ntr_adata.var["pass_basic_filter"] = ntr_adata.var["pass_basic_filter"] * (ntr_adata.X <= 10000).toarray().all(0)
dyn.pp.calc_sz_factor(ntr_adata)
dyn.pp.cell_cycle_scores(ntr_adata) # Bulk-data, meanlingless.
ntr_adata.var["NTR"] = add_var_ntr(ntr_adata) # NTR per gene
ntr_adata.obs["NTR"] = add_obs_ntr(ntr_adata, ntr_base_dir) # NTR per sample
ntr_adata.layers["NTR"] = add_x_ntr(ntr_adata, ntr_base_dir) # NTR per gene per sample
selected_genes = ntr_adata.var["pass_basic_filter"]
selected_samples = ntr_adata.obs["pass_basic_filter"] if "pass_basic_filter" in ntr_adata.obs else [True] * len(ntr_adata)
ntr_adata_sub = ntr_adata[selected_samples, selected_genes].copy()
ntr_tab = pds.DataFrame(ntr_adata_sub.layers["NTR"].A, index=ntr_adata_sub.obs_names, columns=ntr_adata_sub.var_names) #.merge(ntr_adata_sub.obs, how="inner", left_index=True, right_index=True)


# Check PCA using NTR
pca = PCA(n_components=50)
pca.fit(ntr_tab.T)
pca.get_covariance()
pc1_vp, pc2_vp, *_ = pca.explained_variance_ratio_ * 100
pca_tab = pds.DataFrame(pca.components_.T, index=ntr_tab.index, columns=[f"PC{i}" for i in range(1, 51)])

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
fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "ntr.pca.pdf")
fig.clear()
plt.close(fig)


# Cluster of samples based on NTR
batches = ['240409_Lib_embryo', '240612_Lib_28region', '240620_Lib_38region', '240703_Lib_32region', '240710_Lib_37region', '240717_Lib_28region']
cell_types = ["Ectoderm", "Endoderm", "Mesoderm"]
colors = mpl.colormaps.get_cmap("Dark2").colors[:3]
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
grid.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "in_house.new_to_total_ratio.cluster_by_pca.pdf")


# NTR clustering
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
    
fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "in_house.new_to_total_ratio.pdf")
fig.clear()
plt.close(fig)


#
## Splicing information
#
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
fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "scvelo.splicing_based.umap.pdf")
fig.clear()

# Save the anndata into h5ad format
spliced_adata.write_h5ad(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "scvelo.velocity_analysis.v2.h5ad")


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
_, _, labeled_adata, _ = load_data()
del labeled_adata.layers["counts"]
del labeled_adata.layers["NTR"]
pp = Preprocessor(cell_cycle_score_enable=True)
pp.config_seurat_recipe(labeled_adata)
pp.preprocess_adata_seurat(labeled_adata, tkey="Groups", experiment_type="mix_std_stm")
scp.pp.combat(labeled_adata, key="Batches")
pp.pca(labeled_adata)

# Calculate dynamics and reduce dimension
dyn.tl.dynamics(labeled_adata)
dyn.tl.reduceDimension(labeled_adata)

# Velocity in UMAP space
dyn.tl.cell_velocities(labeled_adata, basis="umap", method='pearson', other_kernels_dict={'transform': 'sqrt'})
dyn.vf.VectorField(labeled_adata, basis='umap', velocity_key="velocity_N", M=1000, pot_curl_div=True)

dyn.tl.cell_wise_confidence(labeled_adata, ekey = "X_total", vkey = "velocity_N")
dyn.tl.confident_cell_velocities(labeled_adata, group="Cell_types", ekey="X_total", vkey="velocity_N", lineage_dict={'Endoderm': 'Ectoderm'})
dyn.vf.rank_velocity_genes(labeled_adata, vkey="velocity_N")

# Velocity in PCA space
dyn.tl.cell_velocities(labeled_adata, basis="pca", method='pearson', other_kernels_dict={'transform': 'sqrt'})
dyn.vf.VectorField(labeled_adata, basis='pca', vkey="velocity_N", M=1000, pot_curl_div=True)

# Velocity in corn space
dyn.tl.cell_velocities(labeled_adata, basis="corn", method='pearson', other_kernels_dict={'transform': 'sqrt'})
dyn.vf.VectorField(labeled_adata, basis='corn', velocity_key="velocity_N", M=1000, pot_curl_div=True)

# Calculate parameters in new constructed vector fields
dyn.vf.curl(labeled_adata, basis='umap')
dyn.vf.speed(labeled_adata, basis='pca')
dyn.vf.divergence(labeled_adata, basis='pca')
dyn.vf.acceleration(labeled_adata, basis='pca')
dyn.vf.curvature(labeled_adata, basis='pca')

# Plot vector field energy alterations
fig, axe = plt.subplots(1, 1, constrained_layout=True, figsize=(10, 4))
p = dyn.pl.plot_energy(labeled_adata, basis='umap', fig=fig)
fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "dynamo.energy_and_energy_change_rate.pdf")
fig.clear()

# Define the lims of UMAP related figures
umap_x_lim, umap_y_lim = umap_xy_lims(labeled_adata)

# Vector field topography and velocity stream line
general_kwargs = dict(pointsize=0.75, s_kwargs_dict={"alpha": 0.7}, frontier=True, ekey="M_n", vkey="velocity_N", save_show_or_return="return", show_arrowed_spines=True)
fig, ((axe1, axe2), (axe3, axe4)) = plt.subplots(2, 2, constrained_layout=True, tight_layout=True, figsize=(12, 10))
_ = dyn.pl.topography(labeled_adata, color='Cell_types', ax=axe1, **general_kwargs) # Plot topography
_ = dyn.pl.streamline_plot(labeled_adata, color="Cell_types", ax=axe2, **general_kwargs) # Plot velocity
_ = dyn.pl.streamline_plot(labeled_adata, color="Regions", ax=axe3, **general_kwargs) # Plot velocity
_ = dyn.pl.streamline_plot(labeled_adata, basis="corn", color="Regions", ax=axe4, x=1, y=0, **general_kwargs) # Plot velocity
for idx, per_axe in enumerate([axe1, axe2, axe3, axe4]):
    if idx == 3:
        axe4.set_xlim((-20, 20))
        axe4.set_ylim((-1, 18))
    else:
        per_axe.set_xlim(umap_x_lim)
        per_axe.set_ylim(umap_y_lim)
# fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "dynamo.topography_and_streamline.pdf")
fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "dynamo.topography_and_streamline.png")
fig.clear()

# Other dynamics
general_kwargs = dict(quiver_length=4, quiver_size=4, pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, save_show_or_return="return")
fig, ((axe1, axe2), (axe3, axe4)) = plt.subplots(2, 2, constrained_layout=True, tight_layout=True, figsize=(12, 10))
_ = dyn.pl.streamline_plot(labeled_adata, color='divergence_pca', ax=axe1, **general_kwargs)
_ = dyn.pl.streamline_plot(labeled_adata, color='speed_pca', ax=axe2, **general_kwargs)
_ = dyn.pl.streamline_plot(labeled_adata, color='acceleration_pca', ax=axe3, **general_kwargs)
_ = dyn.pl.streamline_plot(labeled_adata, color='curvature_pca', ax=axe4, **general_kwargs)
for per_axe in [axe1, axe2, axe3, axe4]:
    per_axe.set_ylim(umap_y_lim)
    per_axe.set_xlim(umap_x_lim)
#fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "dynamo.integrative_analysis.pdf")
fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "dynamo.integrative_analysis.png")
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
fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "dynamo.marker_gene_expression.pdf")


# In silico perturbation
gene = "Mid1"
dyn.pd.perturbation(labeled_adata, gene, [-10000], emb_basis="umap")

fig, ((axe1, axe2), (axe3, axe4)) = plt.subplots(2, 2, figsize=(10, 10), layout="tight")
dyn.pl.streamline_plot(labeled_adata, color="Regions", basis="umap", pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, ax=axe1, save_show_or_return="return")
axe1.set_title("Velocity before perturbation")
dyn.pl.umap(labeled_adata, color='umap_ddhodge_potential', pointsize=0.5, alpha=0.7, frontier=True, ax=axe2, save_show_or_return="return")
axe2.set_title("Dyanamics")
dyn.pl.streamline_plot(labeled_adata, color=gene, pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, ax=axe3, save_show_or_return="return")
axe3.set_title(f"Gene expression {gene}")
dyn.pl.streamline_plot(labeled_adata, color="Regions", basis="umap_perturbation", pointsize=0.5, s_kwargs_dict={"alpha": 0.7}, ax=axe4, save_show_or_return="return")
axe4.set_title("Velocity after perturbation")
for per_axe in [axe1, axe2, axe3, axe4]:
    per_axe.set_ylim(y_lim)
    per_axe.set_xlim(x_lim)
# fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / ("dynamo.perturbation_" + gene + ".pdf"))
fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / ("dynamo.perturbation_" + gene + ".png"))
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
fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "dynamo.outliers.pdf")
fig.clear()

fig, axe = plt.subplots(1, 1, figsize=(5, 7))
color_map = {"Ectoderm": "red", "Endoderm": "blue", "Mesoderm": "green"}
colors = labeled_adata.obs.Cell_types.apply(lambda x: color_map[x]).tolist()
axe.scatter(x = labeled_adata.obsm["X_corn"][:, 1], y = labeled_adata.obsm["X_corn"][:, 0], c=colors, s=500, alpha=0.5)
labels = [".".join([x, y]) for x, y in zip(labeled_adata.obs.Layer_region, labeled_adata.obs.Sampling_dates)]
for px, py, pt in zip(labeled_adata.obsm["X_corn"][:, 1], labeled_adata.obsm["X_corn"][:, 0], labels):
    axe.text(px, py, pt, alpha=0.5, fontsize="x-small", ha="center", va="center")
fig.savefig(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "dynamo.corn_plot.sample_positions.pdf")
fig.clear()


# Save the anndata into h5ad format
labeled_adata.write_h5ad(PROJECT_DIR / "outputs/analysis/velocity" / VERSION / "dynamo.velocity_analysis.h5ad")
dyn.session_info()

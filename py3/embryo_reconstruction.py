#!/usr/bin/env python3
# File: embryo_reconstruction.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jul 01, 2025
# Updated:

# https://www.nature.com/articles/s41467-023-41482-5
# https://doi.org/10.1186/s13059-024-03347-y

import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

from typing import Literal

import copy
import math
import logging
import itertools
from pathlib import Path
from functools import partial

import anndata as adt
import polars as pls
import scanpy as scp
import pandas as pds

import sklearn as skl
from sklearn import cluster
# from sklearn.cross_decomposition import CCA
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise

import numpy as npy
from numpy import random
from numpy.random import randint

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

from scipy.sparse import csc_matrix, csr_matrix
from scipy.stats import spearmanr, pearsonr


CELL_CYCLE_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1",
    "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2",
    "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1",
    "CHAF1B", "BRIP1", "E2F8", "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
    "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11",
    "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE",
    "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
]

BSNH = "Biospatial_neighborhood"


class SciColorMaps:
    def __init__(self, name: str):
        self._name = name
        self._color_map = self._define_color_map(name)

    @property
    def name(self):
        return self._name

    def __call__(self, n=None):
        return self._color_map

    def _define_color_map(self, name):
        if name.lower() in ["npg"]:
            return ListedColormap([
                "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148",
                "#B09C85"
            ], name=name)
        elif name.lower() in ["nejm"]:
            return ListedColormap([
                "#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97"
            ], name=name)
        elif name.lower() in ["lancet"]:
            return ListedColormap([
                "#00468B", "#EDOOOO", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#ADOO2A", "#ADB6B6", "#1B1919"
            ], name=name)
        elif name.lower() in ["jama"]:
            return ListedColormap([
                "#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "#6A6599", "#80796B"
            ], name=name)
        elif name.lower() in ["bmj"]:
            return ListedColormap([
                "#2A6EBB", "#FOABOO", "#C50084", "#7D5CC6", "#E37222", "#69BE28", "#00B2A9", "#CD202C", "#747678"
            ], name=name)
        else:
            return ListedColormap([
                "red", "orange", "yellow", "pink", "green", "blue", "purple", "brown", "grey", "black", "white"
            ], name=name)

    def register_color(self):
        if self._name not in mpl.colormaps:
            mpl.colormaps.register(self._color_map, name=self._name)


class LogManager(logging.Logger):
    def __init__(self, name, level=logging.INFO, logstream: bool = True,
                 logfile: str = ""):
        super(LogManager, self).__init__(name)

        fmt = logging.Formatter("{levelname: >8} | {asctime} | {name: >30} | {message} ...", "%Y%m%d,%H:%M:%S", "{")
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
    file_pattern: str = "*.h5ad",
    logman: LogManager = LogManager("LoadAnnData")
):
    """Create an AnnData object from 10x read count matrix."""
    if isinstance(in_dir, str): in_dir = Path(in_dir)

    in_files = list(in_dir.glob(file_pattern))
    if len(in_files) == 0:
        logman.critical(f"No files found for {in_dir}")
        return None

    adata_list = []
    for per_adata in in_files:
        if per_adata.is_file():
            adata_list.append(scp.read_h5ad(per_adata))
        elif per_adata.is_dir():
            adata_list.append(scp.read_10x_mtx(per_adata))
        else:
            logman.warning(f"{per_adata} is not a file or a dir.")

    if len(adata_list) == 0:
        logman.critical(f"No AnnData objects found for {in_dir}")
        return None

    return adt.concat(adata_list, join="outer")


def cluster_cells(
    adata: adt.AnnData,
    by: str = "umap",
    suffix: str | None = None,
    n_comps: int | None = 2,
    marker_features: list[str] | None = None,
    logman: LogManager = LogManager("ClusterCells"),
    **kwargs
):
    if by.lower() in ["spatial", "x_spatial"]:
        # Data
        if "X_spatial" in adata.obsm:
            positions = adata.obsm["X_spatial"]
        else:
            raise ValueError(f"Spatial data not found in {adata}.")
    elif by.lower() in ["umap", "x_umap", "pca", "x_pca", "feature", "features"]:
        # Data
        if by.lower() in ["features", "feature"]:
            if marker_features is None:
                raise ValueError("Marker features must be provided when by='feature'.")
            positions = adata[:, marker_features].X
        elif by.lower() in ["pca", "x_pca"]:
            if "X_pca" in adata.obsm:
                positions = adata.obsm["X_pca"][:, :n_comps]
            else:
                raise ValueError(f"PCA data not found in {adata}.")
        else:
            if by.lower() not in ["umap", "x_umap"]:
                logman.warning(f"Unknown method: {by}. Using default method: umap.")
            if "X_umap" in adata.obsm:
                positions = adata.obsm["X_umap"][:, :n_comps]
            else:
                raise ValueError(f"UMAP data not found in {adata}.")
    else:
        raise ValueError(f"Unknown method: {by}")

    # Model
    damping = kwargs.get("damping", 0.5)
    max_iter = kwargs.get("max_iter", 200)
    random_state = kwargs.get("random_state", 0)
    affinity = kwargs.get("affinity", "euclidean")
    convergence_iter = kwargs.get("convergence_iter", 15)
    model = cluster.AffinityPropagation(
        damping=damping, max_iter=max_iter, convergence_iter=convergence_iter, affinity=affinity,
        random_state=random_state
    )

    # Clustering results
    labels = model.fit_predict(positions)

    cluster_name = f"Cluster" if suffix is None else f"Cluster_{suffix}"
    return pls.DataFrame({"Barcode": adata.obs_names, cluster_name: labels}).with_columns(pls.col(cluster_name).cast(pls.String))


def find_biospatial_neighbors(
    adata: adt.AnnData,
    slice_key: str | None = None,
    in_place: bool = True,
    logman: LogManager = LogManager("FindNeighbors")
):
    """Group cells into neighborhoods accounting for biologically and spatial distances.

    The biologcal distances are determined by the expression profiles of each cell, while the spatial distances are
    distantiated by the spatial positions of each cell.
    """
    cluster_tbl_list = []
    all_groups = adata.obs.loc[:, slice_key].unique().tolist()
    for per_group in all_groups:
        # Select observations for each layer
        selected_obs = adata.obs.loc[adata.obs.loc[:, slice_key] == per_group, :].index
        sub_adata = adata[selected_obs, :]

        # Cluster observations by biological profiles.
        logman.info(f"Finding biological neighbors for {per_group}")
        nh_biological = cluster_cells(sub_adata, by="umap", n_comps=2, suffix="bio")

        # Spatial neighbors are defined as cells with similar spatial positions
        logman.info(f"Finding spatial neighbors for {per_group}")
        nh_spatial = cluster_cells(sub_adata, by="spatial", suffix="spatial")

        # Define neighborhood based on spatial and biological neighbors.
        neighborhoods = nh_biological.join(
            nh_spatial, how="outer", on="Barcode"
        ).with_columns(
            (pls.col('Cluster_bio') + "-" + pls.col('Cluster_spatial')).alias(BSNH)
        ).select('Barcode', BSNH)

        cluster_tbl_list.append(neighborhoods)

    clusters = pls.concat(cluster_tbl_list)
    if in_place:
        adata.uns[BSNH] = clusters
        adata.obs[BSNH] = clusters[BSNH]
        return

    return clusters


def define_marker_features(
    adata: adt.AnnData,
    marker_features: list[str] | str | None,
    marker_feature_key: str = "marker_features",
    group_by: str | None = None,
    rank_method: Literal["wilcoxon", "t-test"] | None = "wilcoxon",
    deg_pval_adj: float | None = 0.05,
    deg_top_n: int | None = 10,
    deg_log2fc: float | None = 1,
    logman: LogManager = LogManager("DefineMarkerFeatures")
):
    logman.info("Defining marker features ...")
    if marker_features is not None:
        if isinstance(marker_features, str):
            marker_features = [marker_features]
    elif marker_feature_key in adata.var:
        marker_features = adata.var.pipe(lambda x: x[marker_feature_key][x[marker_feature_key]].index.tolist())
    elif deg_pval_adj is not None and deg_top_n is not None and deg_log2fc is not None and group_by is not None:
        per_adata = scp.tl.rank_genes_groups(adata, groupby=group_by, method=rank_method, copy=True)
        if per_adata is None: raise ValueError("Failed to find marker features.")

        deg_df = scp.get.rank_genes_groups_df(per_adata, group=None)
        if deg_df is None: raise ValueError("Failed to find marker features.")

        marker_features = (pls.from_pandas(deg_df)
                           .filter( pls.col("pvals_adj") < deg_pval_adj, pls.col("logfoldchanges").abs() >= deg_log2fc)
                           .group_by("group").head(deg_top_n).get_column("names").unique().to_list())
    else:
        raise ValueError("Marker features must be provided.")

    # If no min_marker_expr or min_marker_features are provided, return None
    adata.var[marker_feature_key] = adata.var_names.isin(marker_features)


def compute_cluster_distances(
    adata: adt.AnnData,
    slice_key: str,
    reference_slice_name: str,
    query_slice_name: str,
    query_slice_angle: float,
    neighborhood_key: str | None = BSNH,
    neighborhood: pls.DataFrame | None = None,
    minimal_neighborhood_size: int | None = 10,
    marker_feature_key: str | None = "marker_features",
    layer_key: str | None = None,
    logman: LogManager = LogManager("ComputeClusterDistances"),
    **kwargs
):
    logman.info(f"Computing distance matrix for {reference_slice_name} and {query_slice_name}")

    # Ensure the slice information is available in the adata.obs
    if slice_key not in adata.obs:
        raise ValueError(f"Slice {slice_key} not found in {adata}.")

    # Ensure every provided slice is in the adata.obs[slice_key]
    for per_slice in [reference_slice_name, query_slice_name]:
        if per_slice not in adata.obs[slice_key].unique():
            raise ValueError(f"Slice {per_slice} not found in {adata}.")

    # Neighborhood information, either from the adata.obs or given by the `neighborhood`
    if neighborhood_key is None and neighborhood is None:
        raise ValueError("Either neighborhood_key or neighborhoods must be provided.")
    elif neighborhood is None:
        neighborhood = pls.from_pandas(adata.obs.reset_index(names="Barcode"))
    elif neighborhood_key is None:
        assert all([x in neighborhood.columns for x in ["Barcode", slice_key, neighborhood_key]])
    else:
        raise ValueError("Both neighborhood_key and neighborhoods cannot be provided.")

    # Define neighborhood to work on
    if neighborhood is None: raise ValueError("neighborhood is None!")
    candidate_nhs = (neighborhood
                     .select("Barcode", slice_key, neighborhood_key)
                     .group_by(slice_key, neighborhood_key)
                     .agg(pls.count().alias("Count"))
                     .filter(pls.col("Count") >= minimal_neighborhood_size))
    reference_slice_nhs = candidate_nhs.filter(pls.col(slice_key) == reference_slice_name).get_column(neighborhood_key)
    query_slice_nhs = candidate_nhs.filter(pls.col(slice_key) == query_slice_name).get_column(neighborhood_key)

    # The data matrix to work on, either the X or any matrix in adata.layers
    if layer_key is None:
        data_mat = adata.X.toarray() if isinstance(adata.X, (csr_matrix, csc_matrix)) else adata.X
    elif layer_key in adata.layers:
        layer_data = adata.layers[layer_key]
        data_mat = layer_data.toarray() if isinstance(layer_data, (csr_matrix, csc_matrix)) else layer_data
    else:
        raise ValueError(f"Layer {layer_key} not found in {adata}.")

    neighborhood_distance_dict = {reference_slice_name: [], query_slice_name: [], "Distance": []}
    marker_feature_index, *_ = npy.where(adata.var[marker_feature_key])
    data_mat = data_mat[:, marker_feature_index]
    for ref_nhs, qry_nhs in itertools.product(reference_slice_nhs, query_slice_nhs):
        reference_slice_idx, *_ = npy.where(adata.obs[neighborhood_key] == ref_nhs)
        query_slice_idx, *_ = npy.where(adata.obs[neighborhood_key] == qry_nhs)
        reference_slice_mat, query_slice_mat = data_mat[reference_slice_idx, :], data_mat[query_slice_idx, :]

        rotate_slice(query_slice_mat, query_slice_angle)

        distance = pairwise.pairwise_distances(reference_slice_mat, query_slice_mat).mean()

        # collection results
        neighborhood_distance_dict[reference_slice_name].append(ref_nhs)
        neighborhood_distance_dict[query_slice_name].append(qry_nhs)
        neighborhood_distance_dict["Distance"].append(distance)

    return pls.DataFrame(neighborhood_distance_dict)


def rotate_slice(positions, angle, rotate_center=(0, 0), scale=1, shift=(0, 0)):
    theta = npy.deg2rad(angle)
    rotation_mat = npy.array([[npy.cos(theta), -npy.sin(theta)], [npy.sin(theta),  npy.cos(theta)]])
    shifted = (positions - npy.array(rotate_center)) * scale + shift
    rotated = shifted @ rotation_mat.T + npy.array(rotate_center)

    return rotated


def define_adjustments(
    adata: adt.AnnData,
    slice_key: str,
    slice_order: list[str],
    neighborhood_key: str | None = "Biospatial_neighborhood",
    minimal_anchor_neighborhoods: int = 10,
    corner_stone_slice: str | None = None,
    logman: LogManager = LogManager("DefineAdjustments"),
    **kwargs,
):
    logman.info("Defining adjustments ...")
    for slice_idx, reference_slice in enumerate(slice_order):
        if slice_idx == len(slice_order) - 1: break
        query_slice = slice_order[slice_idx + 1]

        for degree in range(1, 360, 10): # rotate the slice and calculate the distance between anchor neighborhoods.
            nh_dist_tbl = compute_cluster_distances(
                adata,  slice_key, reference_slice, query_slice, degree, neighborhood_key=neighborhood_key, **kwargs
            ).sort("Distance").head(minimal_anchor_neighborhoods)

            nh_dist_tbl.select([reference_slice, query_slice]).get_columns()


def adjust_slice(
    adata: adt.AnnData,
    logman: LogManager = LogManager("RotateSlice"),
    **kwargs,
):
    logman.info(f"Adjusiting slice")
    if "X_spatial" in adata.obsm.keys():
        positions = copy.deepcopy(adata.obsm["X_spatial"])
    elif 'x' in adata.obs.columns and 'y' in adata.obs.columns:
        positions = copy.deepcopy(adata.obs.loc[:, ['x', 'y']])
    else:
        raise ValueError(f"Spatial data not found in {adata}.")


def align_slices(
    adata: adt.AnnData,
    slice_key: str,
    marker_features: list[str] | None = None,
    deg_pval_adj: float | None = 1,
    deg_top_n: int | None = 10,
    logman: LogManager = LogManager("AlignSlices"),
):
    logman.info("Aligning slices...")
    define_marker_features(adata, marker_features, group_by=slice_key, deg_pval_adj=deg_pval_adj, deg_top_n=deg_top_n)
    find_biospatial_neighbors(adata, slice_key=slice_key)
    # compute_cluster_distances(adata, BSNH)
    define_adjustments(adata, slice_key=slice_key, slice_order=adata.obs[slice_key].unique().tolist())
    adjust_slice(adata)


def obtain_data_vec(
    adata: adt.AnnData,
    feature: str = "total counts",
    layer: str = "X",
    logman: LogManager = LogManager("ObtainDataVec")
):
    """Obtain data vector for a feature from X or a layer in an AnnData object."""
    count_vec = None
    if feature == "total counts": # total counts from X
        if isinstance(adata.X, (csr_matrix, csc_matrix)):
            count_vec = adata.X.toarray().sum(1)
        elif isinstance(adata.X, npy.ndarray):
            count_vec = adata.X.sum(1)
    elif feature in adata.layers: # Sum from a layer
        layer_data = adata.layers[feature]
        if isinstance(layer_data, (csr_matrix, csc_matrix)):
            count_vec = layer_data.toarray().sum(1)
        elif isinstance(layer_data, npy.ndarray):
            count_vec = layer_data.sum(1)
    elif feature in adata.var_names: # Any feature in var_names, data of a layer or X.
        feature_idx = npy.where(adata.var_names == feature)[0]
        if layer == "X":
            if isinstance(adata.X, (csr_matrix, csc_matrix)):
                count_vec = adata.X.toarray()[:, feature_idx]
            elif isinstance(adata.X, npy.ndarray):
                count_vec = adata.X[:, feature_idx]
        else:
            layer_data = adata.layers[layer]
            if isinstance(layer_data, (csr_matrix, csc_matrix)):
                count_vec = layer_data.toarray()[:, feature_idx]
            elif isinstance(layer_data, npy.ndarray):
                count_vec = layer_data[:, feature_idx]
    else:
        logman.critical(f"Unknown feature: {feature}. Valid ones: `total counts`, features in `var_names` or `layers`")

    if count_vec is not None:
        count_vec = count_vec.squeeze()

    return count_vec


def visualize(
    adata: adt.AnnData,
    feature: str = "total counts",
    layer: str = "X",
    color_map: str = "tab20",
    color_by: str = "feature",
    save_to: str | Path = "plots.pdf",
    gif_elev: int = 20,
    gif_frames: int = 72,
    gif_fps: int = 10,
    use_relative_z: bool = True,
    v_jitter: bool | float = True,
    layout: str = "stacked"
):
    '''Plot gene expression in spatial.'''
    count_vec = obtain_data_vec(adata, feature, layer)
    if count_vec is None:
        print(f"Unknown type of for {feature}")
        return None

    meta_data = pls.from_pandas(adata.obs)
    if use_relative_z:
        meta_data = meta_data.with_columns(pls.col('z').rank('dense'))

    data_mat = meta_data.with_columns(
        pls.Series(count_vec).alias(feature),
        (pls.col('x') - (pls.col('x').median().over('sample_id') - pls.col('x').mean())).alias('x_adj'),
        (pls.col('y') - (pls.col('y').median().over('sample_id') - pls.col('y').mean())).alias('y_adj'),)

    x_adj, y_adj = data_mat['x_adj'], data_mat['y_adj']
    if v_jitter:
        jitter_width = 0.25 if isinstance(v_jitter, bool) else v_jitter
        z_adj = data_mat['z'] + random.normal(0, jitter_width, data_mat.height)
    else:
        z_adj = data_mat['z']

    # Plotting
    fig = plt.figure()
    axe = fig.add_subplot(projection="3d")

    # Decide the color of feature
    if color_by == "feature":
        feature_color = data_mat[feature]
        _color_map = color_map
    else:
        if color_by not in ["z", "layer", "layers"]:
            print("Unknown way to color samples, using default: layers")

        layers = meta_data['z'].unique().sort().to_list()
        colors = plt.cm.get_cmap(color_map, len(layers)).colors
        color_x_layer = {str(layer): plt.cm.colors.to_hex(color) for (layer, color) in zip(layers, colors)}

        feature_color = data_mat.with_columns(pls.col("z").cast(pls.String).replace(color_x_layer).alias("color"))['color']
        _color_map = None

    # Initial scatters
    plot = axe.scatter(x_adj, y_adj, z_adj, c=feature_color, cmap=_color_map, alpha=0.5)
    fig.colorbar(plot, ax=axe, label=f"{feature}")

    # Create the plot, either animation or static plot.
    save_to = Path(save_to)
    if save_to.suffix == ".gif":
        if layout == "stacked":
            def update_stacked(frame, **kwargs):
                _gif_elev = kwargs["elev"]
                axe.view_init(elev=_gif_elev, azim=frame) # Change azimuth angle
                return axe,

            gif_frames = 360
            update_func = partial(update_stacked, elev=gif_elev)
        elif layout == "one_by_one" and color_by == "feature":
            def update_one_by_one(frame, **kwargs):
                n_ttl_frames = kwargs["n_ttl_frames"]
                color_map = kwargs["color_map"]
                feature_data = kwargs["feature_data"]
                layer_idx, x_adj, y_adj, z_adj = kwargs['layer_idx'], kwargs['x_adj'], kwargs['y_adj'], kwargs['z_adj']
                clear_existing = kwargs["clear_existing"]

                # Determine layer
                n_layers = len(layer_idx.unique())
                if frame <= n_ttl_frames / 2:
                    current_layer = min(math.ceil(2 * frame / n_ttl_frames * n_layers), n_layers)
                else:
                    current_layer = min(math.ceil(2 * (1 - frame / n_ttl_frames) * n_layers), n_layers)

                # Obtain data for current layer
                current_data = layer_idx == current_layer
                x_adj_ = x_adj.filter(current_data)
                y_adj_ = y_adj.filter(current_data)
                z_adj_ = z_adj.filter(current_data)
                feature_data_ = feature_data.filter(current_data)
                vmin, vmax = feature_data.min(), feature_data.max()

                # Plotting and define the canvas
                if clear_existing: axe.clear()
                axe.scatter(x_adj_, y_adj_, z_adj_, c=feature_data_, cmap=color_map, vmin=vmin, vmax=vmax, alpha=0.5)
                #axe.set_title(layer_idx)
                axe.set_xlim((x_adj.min(), x_adj.max()))
                axe.set_ylim((y_adj.min(), y_adj.max()))
                axe.set_zlim((z_adj.min(), z_adj.max()))

                return axe,

            update_func = partial(
                update_one_by_one, n_ttl_frames=gif_frames, feature_data=feature_color, color_map=color_map,
                x_adj=x_adj, y_adj=y_adj, z_adj=z_adj, layer_idx=data_mat['z'], clear_existing=True,)
        else:
            fig.clear()
            plt.close()
            print(f"Unknown way to update the frames: {layout}")
            return None

        ani = FuncAnimation(fig, update_func, repeat=True, frames=gif_frames, interval=1000)
        writer = PillowWriter(fps=gif_fps, metadata=dict(artist='Me'), bitrate=1800)
        ani.save(save_to, writer=writer)
    else:
        fig.savefig(save_to)

    fig.clear()
    plt.close()


PROJECT_DIR = Path("~/Documents/projects/wp_vasaseq").expanduser()

npg = SciColorMaps("npg").register_color()
jama = SciColorMaps("jama").register_color()

# 
sample_info = {
    "OST110091": { "slice_layer":  0, "skip":  True, "batch_id": "20250314_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "OST110093": { "slice_layer":  0, "skip":  True, "batch_id": "20250314_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110112": { "slice_layer":  0, "skip":  True, "batch_id": "20250403_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },

    "STO110107": { "slice_layer": 16, "skip": False, "batch_id": "20250423_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110108": { "slice_layer": 14, "skip": False, "batch_id": "20250423_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110109": { "slice_layer": 11, "skip": False, "batch_id": "20250423_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110114": { "slice_layer":  6, "skip": False, "batch_id": "20250423_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110115": { "slice_layer":  4, "skip": False, "batch_id": "20250423_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },

    "OST110079": { "slice_layer": 17, "skip": False, "batch_id": "20250522_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "OST110080": { "slice_layer": 15, "skip": False, "batch_id": "20250522_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110135": { "slice_layer": 10, "skip": False, "batch_id": "20250522_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110136": { "slice_layer":  8, "skip": False, "batch_id": "20250522_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110137": { "slice_layer":  6, "skip": False, "batch_id": "20250522_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },

    "STO110143": { "slice_layer": 16, "skip": False, "batch_id": "20250604_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110146": { "slice_layer": 10, "skip": False, "batch_id": "20250604_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110147": { "slice_layer":  8, "skip": False, "batch_id": "20250604_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110148": { "slice_layer":  6, "skip": False, "batch_id": "20250604_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110149": { "slice_layer":  4, "skip": False, "batch_id": "20250604_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110150": { "slice_layer":  2, "skip": False, "batch_id": "20250604_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },

    "STO110151": { "slice_layer": 10, "skip": False, "batch_id": "20250728_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110152": { "slice_layer":  8, "skip": False, "batch_id": "20250728_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110153": { "slice_layer":  6, "skip": False, "batch_id": "20250728_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
    "STO110154": { "slice_layer":  4, "skip": False, "batch_id": "20250728_decoded_embryo", "input_file": "k_size_4.min_counts_25.target_dim_0.h5ad" },
}

selected_samples = [

    #"OST110079", # original 17, use as 10
    #"STO110143", # original 16, use as  9
    "STO110107", # original 16, use as  9
    "OST110080", # original 15, use as  8
    #"STO110108", # original 14, use as  7
    "STO110109", # original 11, use as  6
    #"STO110135", # original 10, use as  5
    "STO110146", # original 10, use as  5
    "STO110147", # original  8, use as  4
    "STO110148", # original  6, use as  3
    "STO110149", # original  4, use as  2
    #"STO110115", # original  4, use as  2
    "STO110150", # original  2, use as  1
]

selected_samples = [
    "STO110107", "STO110108", "STO110109", "STO110114", "STO110115",
    "OST110079", "OST110080", "STO110135", "STO110136", "STO110137",
    "STO110143", "STO110146", "STO110147", "STO110148", "STO110149", "STO110150",
    "STO110151", "STO110152", "STO110153", "STO110154",
]


# Load and merge dataset
resolution = 0.35
reload = False
work_dir = PROJECT_DIR / "outputs/analysis/spatial/reconstruction"
processed_out_adata_file = work_dir / f"embryo_reconstruction.processed.resolution_{resolution}.h5ad"
if processed_out_adata_file.exists():
    adata = scp.read_h5ad(processed_out_adata_file)
else:
    raw_out_adata_file = work_dir / "embryo_reconstruction.raw.h5ad"
    if raw_out_adata_file.exists():
        adata = scp.read_h5ad(raw_out_adata_file)
    else:
        adata_list = {}
        adata_folder = work_dir / "preprocess"
        for ii, sample_id in enumerate(selected_samples):
            batch_id = sample_info[sample_id]['batch_id']
            per_adata_file = adata_folder / f'{batch_id}.{sample_id}' / sample_info[sample_id]['input_file']
            if sample_info[sample_id]["skip"]: continue
            if sample_id not in selected_samples: continue

            per_adata = scp.read_h5ad(per_adata_file)
            per_adata.X = per_adata.X.astype("float64")
            for per_layer in per_adata.layers:
                per_adata.layers[per_layer] = per_adata.layers[per_layer].astype("float64")
            per_adata.obs["z"] = sample_info[sample_id]['slice_layer']
            adata_list[sample_id] = per_adata

        # Concat data
        adata = scp.concat(adata_list, join="outer", label="batch")
        adata.obs.index.name = None
        adata.obs_names_make_unique()
        adata.write_h5ad(raw_out_adata_file)

    # Basic filtering
    scp.pp.filter_cells(adata, min_genes=200)
    scp.pp.filter_cells(adata, min_counts=500)
    scp.pp.filter_genes(adata, min_cells=100)
    scp.pp.filter_genes(adata, min_counts=20)

    # Normalization
    scp.pp.normalize_total(adata, target_sum=1e6, exclude_highly_expressed=True, layers="all")
    scp.pp.log1p(adata)
    # scp.pp.scale(adata)

    # Define highly variable genes and PCA
    scp.pp.highly_variable_genes(adata, n_top_genes=3000)
    scp.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver="arpack")

    # Removing batch effects
    #scp.external.pp.harmony_integrate(adata, key="batch_id", basis="X_pca")
    scp.external.pp.harmony_integrate(adata, key="sample_id", basis="X_pca")

    # Clusters by raw
    scp.pp.neighbors(adata, use_rep="X_pca", key_added="neighbors_raw")
    scp.tl.umap(adata, neighbors_key="neighbors_raw")
    adata.obsm["X_umap_raw"] = adata.obsm["X_umap"].copy()
    scp.tl.embedding_density(adata, basis="umap_raw", groupby='sample_id')
    scp.tl.embedding_density(adata, basis="umap_raw", groupby='batch_id')
    scp.tl.leiden(adata, resolution, neighbors_key="neighbors_raw", key_added=f"leiden_{resolution}_raw")
    scp.tl.louvain(adata, resolution, neighbors_key="neighbors_raw", key_added=f"louvain_{resolution}_raw")

    # Clusters by harmony
    scp.pp.neighbors(adata, use_rep="X_pca_harmony")
    scp.tl.umap(adata)
    scp.tl.embedding_density(adata, groupby='sample_id')
    scp.tl.embedding_density(adata, groupby='batch_id')
    scp.tl.leiden(adata, resolution, key_added=f"leiden_{resolution}")
    scp.tl.louvain(adata, resolution, key_added=f"louvain_{resolution}")

    # Visualization
    umap_plot_params = { "size": 80, "alpha": 0.5, "show": False }
    fig, ((axe1, axe2), (axe3, axe4), (axe5, axe6)) = plt.subplots(3, 2, figsize=(16, 12), tight_layout=True)
    _ = scp.pl.umap(adata, color="batch_id", palette="npg", ax=axe1, **umap_plot_params)
    _ = scp.pl.umap(adata, color="sample_id", palette="tab20", ax=axe2, **umap_plot_params)
    _ = scp.pl.umap(adata, color="z", ax=axe3, **umap_plot_params)
    _ = scp.pl.umap(adata, color=f"louvain_{resolution}", palette="jama", ax=axe4, **umap_plot_params)
    _ = scp.pl.umap(adata, color="n_genes", color_map="Blues", ax=axe5, **umap_plot_params)
    _ = scp.pl.umap(adata, color="n_counts", color_map="Greens", ax=axe6, **umap_plot_params)
    for per_axe, title in zip([axe1, axe2, axe3, axe4, axe5, axe6], ["Batch", "Sample", "Slice", f"Louvain {resolution}", "No. of genes", "Read counts"]):
        per_axe.set_title(title, fontsize=12)
    fig.savefig(work_dir / f"plots/umap.properties.pdf")
    fig.clear()

    for key in ["sample", "batch"]: # Density to show samples
        fig = scp.pl.embedding_density(adata, key=f"umap_density_{key}_id", return_fig=True)
        fig.savefig(work_dir / f"plots/per_{key}_density.pdf")
        fig.clear()
    plt.close()

    # Save to disk
    adata.write_h5ad(processed_out_adata_file)




# Define marker genes per clusters
align_slices(adata, "sample_id")

if False:
# Integartion with public dataset, PijuanSala et. al., Nature, 2019
    single_cell_adata_file = PROJECT_DIR / "outputs/analysis/public/PijuanSala_etal_Nature_2019/anndata/PijuanSala_etal_Nature_2019.processed.h5ad"
    single_cell_adata = scp.read_h5ad(single_cell_adata_file)
    single_cell_adata.obs["barcodes"] = single_cell_adata.obs_names
    single_cell_adata.obs["sample_id"] = single_cell_adata.obs["sample"]
    single_cell_adata.obs["batch_id"] = "PijuanSala_etal_Nature_2019"

    spatial_adata_file = work_dir / f"embryo_reconstruction.processed.resolution_{resolution}.h5ad"
    spatial_adata = scp.read_h5ad(spatial_adata_file)
    spatial_adata.obs["celltype"] = "Unknown"

    comb_adata = adt.concat({"spatial": spatial_adata, "single_cell": single_cell_adata}, join="inner", label="source")
    comb_adata.obs.batch_id = comb_adata.obs.batch_id.str.replace("PijuanSala_etal_Nature_2019", "PNat2019", regex=False)
    comb_adata.obs.batch_id = comb_adata.obs.batch_id.str.replace("_decoded_embryo", "", regex=False)
    comb_adata.obs.sample_id = comb_adata.obs.sample_id.astype("category")

    scp.pp.highly_variable_genes(comb_adata, n_top_genes=3000)
    scp.pp.pca(comb_adata, n_comps=50, use_highly_variable=True, svd_solver="arpack")

    # Removing batch effects
    #scp.external.pp.harmony_integrate(adata, key="batch_id", basis="X_pca")
    scp.external.pp.harmony_integrate(comb_adata, key="sample_id", basis="X_pca")

    # Clusters by raw
    scp.pp.neighbors(comb_adata, use_rep="X_pca", key_added="neighbors_raw")
    scp.tl.umap(comb_adata, neighbors_key="neighbors_raw")
    comb_adata.obsm["X_umap_raw"] = comb_adata.obsm["X_umap"].copy()
    scp.tl.embedding_density(comb_adata, basis="umap_raw", groupby='sample_id')
    scp.tl.embedding_density(comb_adata, basis="umap_raw", groupby='batch_id')
    scp.tl.leiden(comb_adata, resolution, neighbors_key="neighbors_raw", key_added=f"leiden_{resolution}_raw")
    scp.tl.louvain(comb_adata, resolution, neighbors_key="neighbors_raw", key_added=f"louvain_{resolution}_raw")

    # Clusters by harmony
    scp.pp.neighbors(comb_adata, use_rep="X_pca_harmony")
    scp.tl.umap(comb_adata)
    scp.tl.embedding_density(comb_adata, groupby='sample_id')
    scp.tl.embedding_density(comb_adata, groupby='batch_id')
    scp.tl.leiden(comb_adata, resolution, key_added=f"leiden_{resolution}")
    scp.tl.louvain(comb_adata, resolution, key_added=f"louvain_{resolution}")

    fig, (axe1, axe2, axe3) = plt.subplots(1, 3, figsize=(21, 4), layout="constrained")
    _ = scp.pl.embedding(comb_adata, basis="umap_raw", color="batch_id", ax=axe1, title="Raw")
    _ = scp.pl.embedding(comb_adata, basis="umap", color="batch_id", ax=axe2, title="Harmonized")
    _ = scp.pl.embedding(comb_adata, basis="umap", color="celltype", ax=axe3, title="Celltype(Harmonized)")
    fig.savefig(work_dir / f"plots/umap.integration_with_public_single_cell_data.pdf")
    fig.clear()

# adata.layers["NTR"] = csr_matrix(adata.layers["nascent_convolved_filled"].toarray() / (adata.layers["total_convolved_filled"].toarray() + 1))
# visualize(adata, "total counts", color_map="Spectral_r", save_to=f"version_2/plots.total_counts.stacked.gif")
# visualize(adata, "total counts", color_map="Spectral_r", save_to=f"version_2/plots.total_counts.one_by_one.gif", layout="one_by_one")
# visualize(adata, "total counts", color_by="layers", save_to=f"version_2/plots.total_counts.pdf", layout="stacked")
# visualize(adata, "nascent_convolved_filled", color_map="Spectral_r", save_to=f"version_2/plots.nascent_counts.stacked.gif")
# visualize(adata, "NTR", color_map="Spectral_r", save_to=f"version_2/plots.new_to_toal_ratio.stacked.gif")

# scp.pp.normalize_per_cell(adata, layers='all', )
# feature_list = [
#     # "Dll1", "Nkx1-2", "Pyy", "Ripply2", "Sox17", # Missing in the data.
#     "Cdx1", "Hhex", "Irx5", "Mesp1", "Pou5f1", "Rapgef5", "Sbk1", "Sox2ot", "Supt3", "T", "Tbx6", "Wnt3a",
# ]
# for feature in feature_list:
#     #visualize(adata, feature, color_map="Spectral_r", save_to=f"version_2/plots.{feature}.stacked.pdf")
#     #visualize(adata, feature, color_map="Spectral_r", save_to=f"version_2/plots.{feature}.stacked.gif")
#     #visualize(adata, feature, color_map="Spectral_r", save_to=f"version_2/plots.{feature}.one_by_one.gif", layout="one_by_one", gif_fps=5)
#     #visualize(adata, feature, "nascent_convolved_filled", color_map="Spectral_r", save_to=f"version_2/plots.{feature}.nascent_expr.one_by_one.gif", layout="one_by_one", gif_fps=5)
#     print(f"----- {feature} done -----")

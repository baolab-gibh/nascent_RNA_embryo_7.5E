#!/usr/bin/env python3
# File: Salmen_etal_rna_velocity.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jan 09, 2025
# Updated:

from pathlib import Path
import warnings

import click
import pandas as pds
import scvelo as scv
import scanpy as scp
import anndata as adt

import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy import sparse

warnings.filterwarnings("ignore", category=DeprecationWarning)

ext_list = ("total.UFICounts.tsv.gz", "uniaggGenes_spliced.UFICounts.tsv.gz", "uniaggGenes_unspliced.UFICounts.tsv.gz")
samples_list = (
    "GSM5369504_E7.5-1_i4", "GSM5369505_E7.5-2_i5", "GSM5369506_E7.5-3_i14", "GSM5369507_E7.5-4_i6",
    "GSM5369508_E7.5-5_i15", "GSM5369509_E7.5-6_i31", "GSM5369511_E7.5-7_i32"
)


def read_ufi_data(ttl_path, spl_path, uspl_path, **kwargs):
    '''Load spliced and non-spliced read counts from the disk and create a AnnData from scvelo.'''
    ttl_tab = pds.read_table(ttl_path, index_col=0)
    spl_tab = pds.read_table(spl_path, index_col=0)
    uspl_tab = pds.read_table(uspl_path, index_col=0)

    features = list(set(ttl_tab.index) & set(spl_tab.index) & set(uspl_tab.index))
    ttl_tab, spl_tab, uspl_tab = ttl_tab.loc[features, :], spl_tab.loc[features, :], uspl_tab.loc[features, :]
    row_matches = all([a == b == c for a, b, c in zip(ttl_tab.index, spl_tab.index, uspl_tab.index)])
    col_matches = all([a == b == c for a, b, c in zip(ttl_tab.columns, spl_tab.columns, uspl_tab.columns)])
    if not col_matches or not row_matches:
        raise ValueError("The columns or rows of the three files don't match")

    var_list = ttl_tab.index.to_frame().reset_index(drop=True).rename(columns={0: "gene_ids"}).set_index("gene_ids", drop=False)
    obs_list = ttl_tab.columns.to_frame().reset_index(drop=True).rename(columns={0: "cellbarcodes"}).set_index("cellbarcodes")
    for meta_key, meta_val in kwargs.items():
        obs_list[meta_key] = meta_val

    adata = adt.AnnData(sparse.csr_matrix(ttl_tab.to_numpy(dtype="int32").T), var=var_list, obs=obs_list)
    adata.layers["spliced"] = sparse.csr_matrix(spl_tab.to_numpy(dtype="int32").T)
    adata.layers["unspliced"] = sparse.csr_matrix(uspl_tab.to_numpy(dtype="int32").T)

    return adata


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


def preprocess(sample_sheet, out_dir: Path, over_write: bool, n_pcs: int, cluster_resolution: float, batch_key: str):
    """Preprocess the data and save to h5ad file."""
    h5ad_save_to = out_dir / "velocity.h5ad"
    if h5ad_save_to.exists() and not over_write:
        adata = adt.read_h5ad(h5ad_save_to)
    else:
        adata_list = list()
        with open(sample_sheet) as ss:
            for per_run in ss:
                if per_run.startswith("#"): continue
                sample_id, ttl_path, spl_path, uspl_path, *_ = per_run.strip().split(",")
                per_adata = read_ufi_data(ttl_path, spl_path, uspl_path, sample_id=sample_id, stage="E7.5")
                adata_list.append(per_adata)
                
        adata = adt.concat(adata_list)
        adata.write_h5ad(out_dir / "raw.h5ad")
        scv.pp.filter_and_normalize(adata, min_counts=3, min_counts_u=1, min_cells=3, min_cells_u=3)

        adata.var["mt"] = adata.var_names.str.contains("_mt-")
        adata.var["ribo"] = adata.var_names.str.contains("_Rp[sl]")
        adata.var["hb"] = adata.var_names.str.contains("_Hb[^(p)]")
        scp.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

        scp.pp.highly_variable_genes(adata, n_top_genes=5000, batch_key=batch_key)
        scp.pp.pca(adata, n_comps=n_pcs)
        scp.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=30)
        scp.tl.umap(adata)
        scp.tl.leiden(adata, flavor="igraph", n_iterations=10, resolution=cluster_resolution)
        scp.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

        # Save the adata in h5ad format
        adata.write_h5ad(out_dir / "velocity.h5ad")

    return adata


@click.command()
@click.argument("sample_sheet", type=click.Path(exists=True))
@click.option("-n", "--n-cpus", default=1, show_default=True, type=int, help="Number of CPUs used to create velocity graph.")
@click.option("-c", "--color-by", default="Cell_type", show_default=True, type=str, help="Which column to be used to color the UMAP.")
@click.option("-b", "--batch-key", default="sample_id", show_default=True, type=str, help="Which column to be used to indicate batches.")
@click.option("-r", "--cluster-resolution", default=1, show_default=True, type=float, help="Resolution for clustering.")
@click.option("-p", "--n-pcs", default=50, show_default=True, type=int, help="Number of PCs used to create velocity graph.")
@click.option("--fig-format", default="png", show_default=True, type=str, help="Format of the figure.")
@click.option("-f", "--over-write", default=False, is_flag=True, help="Overwrite existing h5ad file.")
@click.option("-o", "--out-dir", default="Velocity", show_default=True, type=Path, help="Output directory.")
def main(
    sample_sheet, out_dir: Path, over_write: bool, n_cpus: int, color_by: str, n_pcs: int, cluster_resolution: float,
    fig_format: str, batch_key: str
):
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = preprocess(sample_sheet, out_dir, over_write, n_pcs, cluster_resolution, batch_key)

    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata, n_jobs=n_cpus)
    scv.tl.velocity_pseudotime(adata)
    scv.tl.velocity_embedding(adata, basis="umap")

    save_to = out_dir / f"velocity.colorby_{color_by}.{fig_format}"
    plots_velocity(adata, color_by=color_by, batch_key=batch_key, save_to=save_to)


if __name__ == "__main__":
    main()

#!/usr/bin/env Rscript
# File: slide_seq.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 12, 2024
# Updated:
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)

  library(Seurat)
})

project_dir <- "~/Documents/projects/wp_vasaseq"


beads_obj <- file.path(project_dir, "outputs/analysis/undecoded_with_polya/spNEB.Solo.out/Gene/filtered") %>%
  Read10X() %>%
  CreateSeuratObject(min.cell = 10, min.features = 5)
p <- VlnPlot(beads_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
figure_save_to <- file.path(project_dir, "outputs/analysis/undecoded_with_polya/plots/spNEB.beads_as_cell.vlnplot.pdf")
ggsave(p, filename = figure_save_to, width = 8.5, height = 6)


beads_obj <- NormalizeData(beads_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
beads_obj <- FindVariableFeatures(beads_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
beads_obj <- ScaleData(beads_obj, verbose = FALSE)
beads_obj <- RunPCA(beads_obj, features = VariableFeatures(object = beads_obj), verbose = FALSE)
beads_obj <- FindNeighbors(beads_obj, verbose = FALSE)
beads_obj <- FindClusters(beads_obj, resolution = 0.5, verbose = FALSE)
beads_obj <- RunUMAP(beads_obj, reduction = "pca", verbose = FALSE)

p <- DimPlot(beads_obj, reduction = "umap", label = TRUE, group.by = "Regions", pt.size = 10) + NoLegend()
umap_save_to <- file.path(project_dir, "outputs/analysis/undecoded_with_polya/plots/spNEB.beads_as_cell.umap.pdf")
ggsave(p, filename = umap_save_to, width = 8.5, height = 6)

# 
mtx <- file.path(project_dir, "outputs/analysis/undecoded_with_polya/spNEB.Solo.out/Gene/raw/matrix.mtx.gz") %>%
  fread(col.names = c("feature_id", "sample_id", "umi_count")) %>%
  dplyr::filter(umi_count < 10000)

plot_tab <- mtx %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarize(total_umi_count = sum(umi_count)) %>%
  dplyr::arrange(desc(total_umi_count)) %>%
  dplyr::mutate(barcode_index = 1:dplyr::n())

p <- ggplot() +
  geom_point(aes(x = barcode_index, y = total_umi_count), plot_tab) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(x = NULL, y = "Total UMI count")

umi_count_save_to <- file.path(project_dir, "outputs/analysis/undecoded_with_polya/plots/spNEB.beads_as_cell.raw.umi_count.pdf")
ggsave(p, filename = umi_count_save_to, width = 8.5, height = 6)

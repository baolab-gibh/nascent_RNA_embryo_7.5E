#!/usr/bin/env Rscript
# File: enhancer_rna.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 06, 2024
# Updated:


suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(RColorBrewer)

  library(DESeq2)
  library(Seurat)
  library(ggsci)
  library(patchwork)
  library(harmony)
})

rbgc <- brewer.pal(n = 11, name = "RdBu")
npg_3 <- ggsci::pal_npg()(3)

all_batches <- c("240409_Lib_embryo", "240612_Lib_28region", "240620_Lib_38region", "240703_Lib_32region", "240710_Lib_37region", "240717_Lib_28region")
regions_l0 <- c("Ectoderm", "Ectoderm", "Endoderm", "Endoderm", "Ectoderm", "Ectoderm", "Mesoderm", "Mesoderm", "Control") %>%
  purrr::set_names(c('A', 'P', 'EA', 'EP', 'L', 'R', 'MA', 'MP', "NC"))

project_dir <- "~/Documents/projects/wp_vasaseq"


#
## eRNA
#
# Enhancer table
if (FALSE) {
  enhancer_tab <- file.path(project_dir, "outputs/analysis/preprocessing/geo_seq/quantification") %>%
    list.files(pattern = ".enhancer_rna.read_counts.txt$", full.names = TRUE, recursive = TRUE) %>%
    purrr::discard(~stringr::str_detect(.x, "\\.bk")) %>%
    lapply(function(x) {
      batch <- stringr::str_split(x, "/", simplify = TRUE) %>% purrr::keep(~stringr::str_detect(.x, "^[0-9]+_Lib_")) %>% `[`(1)
      fread(x) %>%
        dplyr::rename_with(.col = dplyr::starts_with("/"), .fn = ~(stringr::str_split(.x, "/", simplify = TRUE) %>% purrr::keep(~stringr::str_detect(.x, ".merged.rmdup.bam$")) %>% stringr::str_replace("merged.rmdup.bam$", batch)))
    }) %>%
    purrr::reduce(dplyr::inner_join, by = c("Geneid", "Chr", "Start", "End", "Strand", "Length")) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate()

  features_list <- enhancer_tab %>% dplyr::pull("Geneid") %>% data.frame(., .)
  barcode_id <- enhancer_tab %>% colnames() %>% purrr::keep(~stringr::str_detect(.x, "_Lib_"))
  matrix <- enhancer_tab %>%
    dplyr::select(-c(Chr, Start, End, Strand, Length)) %>%
    tidyr::pivot_longer(-Geneid, names_to = "Barcodes", values_to = "Counts") %>%
    dplyr::mutate(Geneid_idx = match(Geneid, features_list[,1]), Barcodes_idx = match(Barcodes, barcode_id)) %>%
    dplyr::select(Geneid_idx, Barcodes_idx, Counts) %>%
    dplyr::bind_rows(data.frame(Geneid_idx = nrow(features_list), Barcodes_idx = length(barcode_id), Counts = dim(.)[1]), .) %>%
    apply(1, paste, collapse = "\t")

  out_dir <- file.path(project_dir, "outputs/analysis/enhancer_rna/enhancer_10X")

  write.table(features_list, file.path(out_dir, "features.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  # writeLines(features_list, file.path(out_dir, "features.tsv"))
  writeLines(barcode_id, file.path(out_dir, "barcodes.tsv"))
  writeLines(c("%%MatrixMarket matrix coordinate integer general", "%", matrix), file.path(out_dir, "matrix.mtx"))

  # Enhancer RNA, enc or ENC
  enc_embryo <- Read10X(out_dir) %>% CreateSeuratObject(project = "EnhancerRNA", min.cell = 20, min.features = 5)
  tar_samples <- colnames(enc_embryo) %>% purrr::discard(~stringr::str_detect(.x, "NC_"))
  enc_embryo <- enc_embryo[, tar_samples]
  enc_embryo@meta.data <- enc_embryo@meta.data %>% as.data.frame() %>%
    dplyr::mutate(Sample = rownames(.)) %>%
    dplyr::mutate(
      Layers = stringr::str_extract(Sample, "^([0-9]+)[A-Z]+_", group = 1),
      Region_l1 = stringr::str_extract(Sample, "^[0-9]+([A-Z]+)_", group = 1),
      Region_l0 = dplyr::case_when(
        Region_l1 == "A" ~ "Ectoderm",
        Region_l1 == "P" ~ "Ectoderm",
        Region_l1 == "L" ~ "Ectoderm",
        Region_l1 == "R" ~ "Ectoderm",
        Region_l1 == "EA" ~ "Endoderm",
        Region_l1 == "EP" ~ "Endoderm",
        Region_l1 == "MA" ~ "Mesoderm",
        Region_l1 == "MP" ~ "Mesoderm",
        TRUE ~ "Control"
      )) %>%
    dplyr::mutate(Region_l0 = factor(Region_l0, levels = c("Ectoderm", "Mesoderm", "Endoderm", "Control"))) %>%
    dplyr::mutate(Region_l1 = factor(Region_l1, levels = c("A", "P", "L", "R", "MA", "MP", "EA", "EP")))

  enc_embryo <- NormalizeData(enc_embryo, normalization.method = "LogNormalize", scale.factor = 10000)
  enc_embryo <- FindVariableFeatures(enc_embryo, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
  enc_embryo <- ScaleData(enc_embryo, features = rownames(enc_embryo), verbose = FALSE)
  enc_embryo <- RunPCA(enc_embryo, features = VariableFeatures(enc_embryo), verbose = FALSE)
  enc_embryo <- FindNeighbors(enc_embryo, dims = 1:10, verbose = FALSE)
  enc_embryo <- FindClusters(enc_embryo, resolution = 0.8, verbose = FALSE)
  enc_embryo <- RunUMAP(enc_embryo, dims = 1:10, verbose = FALSE)

  p <- DimPlot(enc_embryo, reduction = "pca", label = TRUE, repel = TRUE, pt.size = 2) + NoLegend()
  p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/pca.pdf")
  ggsave(p, filename = p_save_to, width = 4, height = 4)

  p <- DimPlot(enc_embryo, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 2) + NoLegend()
  p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/umap.pdf")
  ggsave(p, filename = p_save_to, width = 4, height = 4)

  p <- DimPlot(enc_embryo, reduction = "umap", group.by = "Region_l0", label = TRUE, repel = TRUE, pt.size = 2)
  p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/umap.Region_l0.pdf")
  ggsave(p, filename = p_save_to, width = 5, height = 4)

  p <- DimPlot(enc_embryo, reduction = "umap", group.by = "Region_l1", label = TRUE, repel = TRUE, pt.size = 2)
  p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/umap.Region_l1.pdf")
  ggsave(p, filename = p_save_to, width = 5, height = 4)

  p <- DimPlot(enc_embryo, reduction = "pca", group.by = "Region_l1", label = TRUE, repel = TRUE, pt.size = 2)
  p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/pca.Region_l1.pdf")
  ggsave(p, filename = p_save_to, width = 5, height = 4)

  p <- DimPlot(enc_embryo, reduction = "pca", group.by = "Region_l0", label = TRUE, repel = TRUE, pt.size = 2)
  p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/pca.Region_l0.pdf")
  ggsave(p, filename = p_save_to, width = 5, height = 4)


  # Find markers
  Idents(enc_embryo) <- "Region_l0"
  DefaultAssay(enc_embryo) <- "RNA"
  de_marker_tab <- FindAllMarkers(enc_embryo, logfc.threshold = 0.001)

  marker_features <- de_marker_tab %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(desc(avg_log2FC)) %>% dplyr::pull(gene) %>% unique()
  marker_features <- de_marker_tab %>% dplyr::filter(gene %in% marker_features, cluster == "Ectoderm") %>% dplyr::arrange(desc(avg_log2FC)) %>% dplyr::pull(gene)
  p <- DotPlot(enc_embryo, features = marker_features, cluster.idents = TRUE) +
    scale_color_gradient2(low = rbgc[11], mid = rbgc[6], high = rbgc[1], name = "Avg. Exp.") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()
  p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/dotplot.marker.pdf")
  ggsave(p, filename = p_save_to, width = 6, height = 12)


  for (example_erna_id in marker_features) {
    p <- FeaturePlot(enc_embryo, features = example_erna_id, pt.size = 2, order = TRUE) + scale_color_gradient(low = rbgc[5], high = rbgc[1])
    p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots", paste0("example_erna.umap.", stringr::str_replace_all(example_erna_id, "[:-]", "_"), ".pdf"))
    ggsave(p, filename = p_save_to, width = 4, height = 4)

    p <- FeaturePlot(enc_embryo, reduction = "pca", features = example_erna_id, pt.size = 2, order = TRUE) + scale_color_gradient(low = rbgc[5], high = rbgc[1])
    p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots", paste0("example_erna.pca.", stringr::str_replace_all(example_erna_id, "[:-]", "_"), ".pdf"))
    ggsave(p, filename = p_save_to, width = 4, height = 4)

    p_regionl0 <- VlnPlot(enc_embryo, assay = "RNA", features = example_erna_id, group.by = "Region_l0") +
      scale_fill_jama() +
      labs(x = NULL) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))
    p_regionl1 <- VlnPlot(enc_embryo, assay = "RNA", features = example_erna_id, group.by = "Region_l1") +
      scale_x_discrete(labels = c("A", "P", "L", "R", "MA", "MP", "EA", "EP")) +
      labs(title = NULL) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

    p <- p_regionl0 + p_regionl1 + plot_layout(guides = "collect")
    p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots", paste0("example_erna.vlnplot.", stringr::str_replace_all(example_erna_id, "[:-]", "_"), ".pdf"))
    ggsave(p, filename = p_save_to, width = 8, height = 4)
  }
}
# Failed due to too few read counts.


#
## Non-coding transcriptional units
#
noncoding_tab <- file.path(project_dir, "outputs/analysis/preprocessing/geo_seq/quantification") %>%
  list.files(pattern = ".non_coding_transcripts.counts.txt$", full.names = TRUE, recursive = TRUE) %>%
  purrr::discard(~stringr::str_detect(.x, "\\.bk")) %>%
  lapply(function(x) {
    batch <- stringr::str_split(x, "/", simplify = TRUE) %>% purrr::keep(~stringr::str_detect(.x, "^[0-9]+_Lib_")) %>% `[`(1)
    fread(x) %>%
      dplyr::rename_with(.col = dplyr::starts_with("/"), .fn = ~(stringr::str_split(.x, "/", simplify = TRUE) %>% purrr::keep(~stringr::str_detect(.x, ".merged.rmdup.bam$")) %>% stringr::str_replace("merged.rmdup.bam$", batch)))
  }) %>%
  purrr::reduce(dplyr::inner_join, by = c("Geneid", "gene_id", "gene_name", "Chr", "Start", "End", "Strand", "Length")) %>%
  dplyr::as_tibble()



# DESeq2 method
noncoding_count_data <- noncoding_tab %>%
  dplyr::select(-c(Chr, Start, End, Strand, Length, gene_id, gene_name, dplyr::starts_with("NC"))) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("Geneid") %>% as.matrix()

noncoding_col_data <- data.frame(Sample = colnames(noncoding_count_data)) %>%
  dplyr::mutate(Batch = stringr::str_extract(Sample, "^.*\\.(.*$)", group = 1)) %>%
  dplyr::mutate(Layer = stringr::str_extract(Sample, "^([0-9]+)[A-Z]+_", group = 1)) %>%
  dplyr::mutate(Region_l1 = stringr::str_extract(Sample, "^[0-9]+([A-Z]+)_", group = 1)) %>%
  dplyr::mutate(Region_l0 = dplyr::case_when(
    Region_l1 == "A" ~ "Ectoderm",
    Region_l1 == "P" ~ "Ectoderm",
    Region_l1 == "L" ~ "Ectoderm",
    Region_l1 == "R" ~ "Ectoderm",
    Region_l1 == "EA" ~ "Endoderm",
    Region_l1 == "EP" ~ "Endoderm",
    Region_l1 == "MA" ~ "Mesoderm",
    Region_l1 == "MP" ~ "Mesoderm",
    TRUE ~ "Control"
  ))

noncoding_dds <- DESeqDataSetFromMatrix(countData = noncoding_count_data, colData = noncoding_col_data, design = ~Batch)
noncoding_keep <- (rowSums(counts(noncoding_dds) >= 10) >= 5)
noncoding_dds <- noncoding_dds[noncoding_keep, ]
noncoding_dds <- DESeq(noncoding_dds)


# Pseudo single cell, by Seurat
features_list <- noncoding_tab %>% dplyr::pull("Geneid") %>% data.frame(., .)
barcode_id <- noncoding_tab %>% colnames() %>% purrr::keep(~stringr::str_detect(.x, "_Lib_"))
matrix <- noncoding_tab %>%
  dplyr::select(-c(Chr, Start, End, Strand, Length, gene_id, gene_name)) %>%
  tidyr::pivot_longer(-Geneid, names_to = "Barcodes", values_to = "Counts") %>%
  dplyr::mutate(Geneid_idx = match(Geneid, features_list[,1]), Barcodes_idx = match(Barcodes, barcode_id)) %>%
  dplyr::select(Geneid_idx, Barcodes_idx, Counts) %>%
  dplyr::bind_rows(data.frame(Geneid_idx = nrow(features_list), Barcodes_idx = length(barcode_id), Counts = dim(.)[1]), .) %>%
  apply(1, paste, collapse = "\t")

out_dir <- file.path(project_dir, "outputs/analysis/enhancer_rna/transcriptional_unit_10X")

write.table(features_list, file.path(out_dir, "features.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
# writeLines(features_list, file.path(out_dir, "features.tsv"))
writeLines(barcode_id, file.path(out_dir, "barcodes.tsv"))
writeLines(c("%%MatrixMarket matrix coordinate integer general", "%", matrix), file.path(out_dir, "matrix.mtx"))

# Transcriptional units, TU or tu
tu_embryo <- Read10X(out_dir) %>% CreateSeuratObject(project = "TranscriptionalUnits", min.cell = 20, min.features = 5)
tar_samples <- colnames(tu_embryo) %>% purrr::discard(~stringr::str_detect(.x, "NC_"))
tu_embryo <- tu_embryo[, tar_samples]
tu_embryo@meta.data <- tu_embryo@meta.data %>% as.data.frame() %>%
  dplyr::mutate(Sample = rownames(.)) %>%
  dplyr::mutate(
    Batch = stringr::str_extract(Sample, "^.*\\.(.*$)", group = 1),
    Layers = stringr::str_extract(Sample, "^([0-9]+)[A-Z]+_", group = 1),
    Region_l1 = stringr::str_extract(Sample, "^[0-9]+([A-Z]+)_", group = 1),
    Region_l0 = dplyr::case_when(
      Region_l1 == "A" ~ "Ectoderm",
      Region_l1 == "P" ~ "Ectoderm",
      Region_l1 == "L" ~ "Ectoderm",
      Region_l1 == "R" ~ "Ectoderm",
      Region_l1 == "EA" ~ "Endoderm",
      Region_l1 == "EP" ~ "Endoderm",
      Region_l1 == "MA" ~ "Mesoderm",
      Region_l1 == "MP" ~ "Mesoderm",
      TRUE ~ "Control"
    )) %>%
  dplyr::mutate(Region_l0 = factor(Region_l0, levels = c("Ectoderm", "Mesoderm", "Endoderm", "Control"))) %>%
  dplyr::mutate(Region_l1 = factor(Region_l1, levels = c("A", "P", "L", "R", "MA", "MP", "EA", "EP")))

tu_embryo <- NormalizeData(tu_embryo, normalization.method = "LogNormalize")
tu_embryo <- FindVariableFeatures(tu_embryo, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
tu_embryo <- ScaleData(tu_embryo, features = rownames(tu_embryo), verbose = FALSE)
tu_embryo <- RunPCA(tu_embryo, features = VariableFeatures(tu_embryo), verbose = FALSE)
tu_embryo <- RunHarmony(tu_embryo, "Batch")
tu_embryo <- FindNeighbors(tu_embryo, dims = 1:10, verbose = FALSE)
tu_embryo <- FindClusters(tu_embryo, resolution = 0.8, verbose = FALSE)
tu_embryo <- RunUMAP(tu_embryo, dims = 1:10, verbose = FALSE)


p <- DimPlot(tu_embryo, reduction = "harmony", label = TRUE, repel = TRUE, pt.size = 2) + NoLegend()
p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/transcriptional_units.pca.pdf")
ggsave(p, filename = p_save_to, width = 4, height = 4)

p <- DimPlot(tu_embryo, reduction = "harmony", label = TRUE, repel = TRUE, pt.size = 2) + NoLegend()
p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/transcriptional_units.umap.pdf")
ggsave(p, filename = p_save_to, width = 4, height = 4)

p <- DimPlot(tu_embryo, reduction = "harmony", group.by = "Region_l0", label = TRUE, repel = TRUE, pt.size = 2)
p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/transcriptional_units.umap.Region_l0.pdf")
ggsave(p, filename = p_save_to, width = 5, height = 4)

p <- DimPlot(tu_embryo, reduction = "harmony", group.by = "Region_l1", label = TRUE, repel = TRUE, pt.size = 2)
p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/transcriptional_units.umap.Region_l1.pdf")
ggsave(p, filename = p_save_to, width = 5, height = 4)

p <- DimPlot(tu_embryo, reduction = "pca", group.by = "Region_l1", label = TRUE, repel = TRUE, pt.size = 2)
p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/transcriptional_units.pca.Region_l1.pdf")
ggsave(p, filename = p_save_to, width = 5, height = 4)

p <- DimPlot(tu_embryo, reduction = "pca", group.by = "Region_l0", label = TRUE, repel = TRUE, pt.size = 2)
p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/transcriptional_units.pca.Region_l0.pdf")
ggsave(p, filename = p_save_to, width = 5, height = 4)


# TU markers
Idents(tu_embryo) <- "Region_l0"
DefaultAssay(tu_embryo) <- "RNA"
de_marker_tab <- FindAllMarkers(tu_embryo, logfc.threshold = 0.001)

marker_features <- de_marker_tab %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(desc(avg_log2FC)) %>% dplyr::pull(gene) %>% unique()
avg_exp_tab <- AverageExpression(tu_embryo, assays = "RNA", group.by = "Region_l0", features = marker_features)[["RNA"]] %>%
  as.data.frame() %>%
  (function(tab) { hclust_res <- dist(tab) %>% hclust(method = "ward.D2"); tab[hclust_res$order,] }) %>%
  tibble::rownames_to_column("Geneid") %>%
  dplyr::mutate(Geneid = stringr::str_extract(Geneid, "(ENSMUST\\d+)\\.", group = 1)) %>%
  dplyr::mutate(Geneid = forcats::fct_inorder(Geneid)) %>%
  tidyr::pivot_longer(cols = -Geneid, names_to = "Region_l0", values_to = "avg_exp")

p <- ggplot(avg_exp_tab) +
  geom_tile(aes(x = Region_l0, y = Geneid, fill = avg_exp)) +
  scale_fill_gradient2(low = rbgc[11], mid = rbgc[6], high = rbgc[1], name = "Avg. Exp.", limits = c(0, 10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank())
p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots/transcriptional_units.heatmap.marker.pdf")
ggsave(p, filename = p_save_to, width = 6, height = 12)


for (example_tu_id in marker_features[1:5]) {
  p <- FeaturePlot(tu_embryo, features = example_tu_id, pt.size = 2, order = TRUE) + scale_color_gradient(low = rbgc[5], high = rbgc[1])
  p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots", paste0("example_erna.umap.", stringr::str_extract(example_tu_id, "(ENSMUST[0-9]+)\\.", group = 1), ".pdf"))
  ggsave(p, filename = p_save_to, width = 4, height = 4)

  p <- FeaturePlot(tu_embryo, reduction = "pca", features = example_tu_id, pt.size = 2, order = TRUE) + scale_color_gradient(low = rbgc[5], high = rbgc[1])
  p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots", paste0("example_erna.pca.", stringr::str_extract(example_tu_id, "(ENSMUST[0-9]+)\\.", group = 1), ".pdf"))
  ggsave(p, filename = p_save_to, width = 4, height = 4)

  p_regionl0 <- VlnPlot(tu_embryo, assay = "RNA", features = example_tu_id, group.by = "Region_l0") +
    scale_fill_jama() +
    labs(x = NULL) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))
  p_regionl1 <- VlnPlot(tu_embryo, assay = "RNA", features = example_tu_id, group.by = "Region_l1") +
    scale_x_discrete(labels = c("A", "P", "L", "R", "MA", "MP", "EA", "EP")) +
    labs(title = NULL) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  p <- p_regionl0 + p_regionl1 + plot_layout(guides = "collect")
  p_save_to <- file.path(project_dir, "outputs/analysis/enhancer_rna/plots", paste0("example_erna.vlnplot.", stringr::str_extract(example_tu_id, "(ENSMUST[0-9]+)\\.", group = 1), ".pdf"))
  ggsave(p, filename = p_save_to, width = 8, height = 4)
}

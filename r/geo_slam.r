#!/usr/bin/env Rscript
# File: geo_slam.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Apr 13, 2025
# Updated: 

suppressPackageStartupMessages({
  library(ggsci)
  library(ggrepel)
  library(patchwork)
  library(tidyverse)
  library(data.table)
  library(RColorBrewer)

  library(DESeq2)
  library(RcppML)
  library(Seurat)
  library(GenomicRanges)

})

project_dir <- "~/Documents/projects/wp_vasaseq"
rbgc <- brewer.pal(n = 11, name = "RdBu")




valid_barcodes <- file.path(project_dir, "outputs/analysis/preprocessing/geo_slam/GEO_SLAM_20250327/barcodes/GEO_SLAM_20250327.sample_barcode_mapping.txt") %>%
  fread(header = FALSE, col.names = c("Sample", "Barcode"))

matrix_dir <- file.path(project_dir, "outputs/analysis/preprocessing/geo_slam/GEO_SLAM_20250327/alignment/GEO_SLAM_20250327.Solo.out/Gene/raw")
embryo <- Read10X(matrix_dir) %>% CreateSeuratObject(min.cell = 10, min.features = 5)

selected_cells <- colnames(embryo) %>% purrr::keep(~.x %in% valid_barcodes$Barcode)
embryo <- embryo[, selected_cells]

saveRDS(embryo, file = file.path(project_dir, "outputs/analysis/geo_slam/overview/20250327_embryo.raw.rds"))

embryo <- NormalizeData(embryo, normalization.method = "LogNormalize", scale.factor = 10000)
embryo <- FindVariableFeatures(embryo, selection.method = "vst", nfeatures = 500)
embryo <- ScaleData(embryo)
embryo <- RunPCA(embryo, features = VariableFeatures(object = embryo), npcs = 20)
embryo <- FindNeighbors(embryo)
embryo <- FindClusters(embryo, resolution = 0.5)
embryo <- RunUMAP(embryo, dims = 1:20)

p_umap <- DimPlot(embryo, reduction = "umap", label = TRUE)
figure_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview/20250327_embryo.umap.pdf")
ggsave(p_umap, filename = figure_save_to, width = 6, height = 6)


# UMI per cell
umi_counts <- file.path(project_dir, "outputs/analysis/preprocessing/geo_slam/GEO_SLAM_20250327/alignment/GEO_SLAM_20250327.Solo.out/Gene/UMIperCellSorted.txt") %>%
  fread(header = FALSE, col.names = c("UMI_count")) %>%
  dplyr::mutate(UMI_rank = 1:dplyr::n())

umi_counts %>%
  dplyr::mutate(potential_cb = ifelse(UMI_rank <= 44, "Cell-associated", "Non-cell-associated")) %>%
  dplyr::group_by(potential_cb) %>%
  dplyr::summarize(mean_UMI = mean(UMI_count), median_UMI = median(UMI_count), sd_UMI = sd(UMI_count), Total = sum(UMI_count))
umi_counts$UMI_count %>% sum
n_cell_assocaited_umi <- 101599 / 371298


p <- ggplot() +
  geom_point(aes(x = UMI_rank, y = UMI_count), umi_counts) +
  geom_text_repel(aes(x=UMI_rank, y=UMI_count, label=UMI_rank), max.overlaps = Inf, umi_counts %>% head(44), size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(x = paste0("Log10(Rank by UMI-counts)", " (n = ", nrow(umi_counts), ")"), y = paste0("Log10(UMI counts)"))
p_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview/20250327_embryo.UMI_count.pdf")
ggsave(p, filename = p_save_to, width = 7, height = 7)



# HGV
top10 <- head(VariableFeatures(embryo), 10)
plot1 <- VariableFeaturePlot(embryo)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
p_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview/20250327_embryo.top10_HVG.pdf")
ggsave(plot = plot1 + plot2, filename = p_save_to, width = 14, height = 7)



# 0321, Luyao
valid_barcodes <- file.path(project_dir, "outputs/analysis/preprocessing/geo_slam/batch_20250321_byPLY/barcodes/batch_20250321_byPLY.sample_barcode_mapping.txt") %>%
  fread(header = FALSE, col.names = c("Sample", "Barcode"))

matrix_dir <- file.path(project_dir, "outputs/analysis/preprocessing/geo_slam/batch_20250321_byPLY/alignment/batch_20250321_byPLY.Solo.out/Gene/raw")
selected_cells <- valid_barcodes$Barcode

embryo <- Read10X(matrix_dir) %>% CreateSeuratObject(min.cell = 10, min.features = 5)
embryo <- embryo[, selected_cells]

embryo <- NormalizeData(embryo, normalization.method = "LogNormalize", scale.factor = 10000)
embryo <- FindVariableFeatures(embryo, selection.method = "vst", nfeatures = 500)
embryo <- ScaleData(embryo)
embryo <- RunPCA(embryo, features = VariableFeatures(object = embryo), npcs = 10)
embryo <- FindNeighbors(embryo)
embryo <- FindClusters(embryo, resolution = 0.5)
embryo <- RunUMAP(embryo, dims = 1:10, n.neighbors = 5)

p_umap <- DimPlot(embryo, reduction = "umap", label = TRUE)
figure_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview/20250321_PLY_kidey.umap.pdf")
ggsave(p_umap, filename = figure_save_to, width = 6, height = 6)

top10 <- head(VariableFeatures(embryo), 10)
plot1 <- VariableFeaturePlot(embryo)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
p_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview/20250321_PLY_kidey.top10_HVG.pdf")
ggsave(plot = plot1 + plot2, filename = p_save_to, width = 14, height = 7)

# PLY, 20250321
umi_counts <- file.path(project_dir, "outputs/analysis/preprocessing/geo_slam/batch_20250321_byPLY/alignment/batch_20250321_byPLY.Solo.out/Gene/UMIperCellSorted.txt") %>%
  fread(header = FALSE, col.names = c("UMI_count")) %>%
  dplyr::mutate(UMI_rank = 1:dplyr::n())

umi_counts %>%
  dplyr::mutate(potential_cb = ifelse(UMI_rank <= 44, "Cell-associated", "Non-cell-associated")) %>%
  dplyr::group_by(potential_cb) %>%
  dplyr::summarize(mean_UMI = mean(UMI_count), median_UMI = median(UMI_count), sd_UMI = sd(UMI_count), Total = sum(UMI_count))

umi_counts$UMI_count %>% sum
n_cell_assocaited_umi <- 9549223 / 10301924


p <- ggplot() +
  geom_point(aes(x = UMI_rank, y = UMI_count), umi_counts) +
  geom_text_repel(aes(x=UMI_rank, y=UMI_count, label=UMI_rank), max.overlaps = Inf, min.segment.length = 0, umi_counts %>% head(16), size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(x = paste0("Log10(Rank by UMI-counts)", " (n = ", nrow(umi_counts), ")"), y = paste0("Log10(UMI counts)"))
p_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview/20250321_PLY_kidey.UMI_count.pdf")
ggsave(p, filename = p_save_to, width = 7, height = 7)


# LDJ 20250321
umi_counts <- file.path(project_dir, "outputs/analysis/preprocessing/geo_slam/batch_20250321_byLDJ/alignment/batch_20250321_byLDJ.Solo.out/Gene/UMIperCellSorted.txt") %>%
  fread(header = FALSE, col.names = c("UMI_count")) %>%
  dplyr::mutate(UMI_rank = 1:dplyr::n())

umi_counts %>%
  dplyr::mutate(potential_cb = ifelse(UMI_rank <= 44, "Cell-associated", "Non-cell-associated")) %>%
  dplyr::group_by(potential_cb) %>%
  dplyr::summarize(mean_UMI = mean(UMI_count), median_UMI = median(UMI_count), sd_UMI = sd(UMI_count), Total = sum(UMI_count))
umi_counts$UMI_count %>% sum
n_cell_assocaited_umi <- 2992923 / 3177512

p <- ggplot() +
  geom_point(aes(x = UMI_rank, y = UMI_count), umi_counts) +
  geom_text_repel(aes(x=UMI_rank, y=UMI_count, label=UMI_rank), max.overlaps = Inf, min.segment.length = 0, umi_counts %>% head(8), size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(x = paste0("Log10(Rank by UMI-counts)", " (n = ", nrow(umi_counts), ")"), y = paste0("Log10(UMI counts)"))
p_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview/20250321_LDJ_kidey.UMI_count.pdf")
ggsave(p, filename = p_save_to, width = 7, height = 7)


gene_nr <- tibble::tribble(
  ~Sample, ~Counts,
  "10A",313,
  "10EA",214,
  "10EP",1729,
  "10L",43,
  "10MA",326,
  "10MP",486,
  "10P",176,
  "10R",354,
  "12A",420,
  "12EA",340,
  "12EP",2125,
  "12L",304,
  "12MA",795,
  "12MP",1009,
  "12P",201,
  "12R",420,
  "15A",655,
  "15EA",184,
  "15EP",1610,
  "15L",120,
  "15MA",966,
  "15MP",794,
  "15P",303,
  "15R",162,
  "2A",269,
  "2EA",156,
  "2EP",217,
  "2P",165,
  "4A",947,
  "4EA",81,
  "4EP",386,
  "4P",550,
  "6A",360,
  "6EA",77,
  "6EP",322,
  "6MA",446,
  "6MP",277,
  "6P",323,
  "8A",454,
  "8EA",301,
  "8EP",1065,
  "8MA",128,
  "8MP",265,
  "8P",100,
) %>%
  dplyr::mutate(Region = stringr::str_extract(Sample, "([A-Z]+)", group = 1))

p <- ggplot(data = gene_nr) +
  geom_violin(aes(x = Region, y = Counts, group = Region)) +
  geom_point(aes(x = Region, y = Counts, color = Region), position = position_jitterdodge(jitter.width = 0.5), size = 2) +
  theme_classic()
p_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview/20250327_embryo.gene_nr.pdf")
ggsave(p, filename = p_save_to, width = 7, height = 4)


# Count details
count_details <- file.path(project_dir, "outputs/analysis/preprocessing/geo_slam/SLAM_20250321_byLDJ/nascent_rna/all_samples/SLAM_20250321_byLDJ.mismatchdetails.tsv") %>%
  fread(header = TRUE) %>%
  dplyr::mutate(Mutation = paste0(Genomic, "->", Read)) %>%
  dplyr::mutate(Mutation = forcats::fct_reorder(Mutation, Mismatches))

# Coverage and mismatches per position (<=35), separated by category (Exonic vs Intronic) 
p1 <- ggplot(count_details %>% dplyr::filter(Category == "Exonic", Position <= 35)) +
  geom_bar(aes(x = Position, y = Coverage, fill = Mutation), stat = "identity", position = "dodge") +
  scale_color_d3(palette = "category20") +
  theme_classic() + theme(legend.position = "top") + guides(fill = guide_legend(ncol = 1)) + labs(title = "Exonic coverage")

p2 <- ggplot(count_details %>% dplyr::filter(Category == "Exonic", Position <= 35)) +
  geom_bar(aes(x = Position, y = Mismatches, fill = Mutation), stat = "identity", position = "dodge") +
  scale_color_d3(palette = "category20") +
  theme_classic() + theme(legend.position = "top") + guides(fill = guide_legend(ncol = 1)) + labs(title = "Exonic mismatches")

p3 <- ggplot(count_details %>% dplyr::filter(Category == "Intronic", Position <= 35)) +
  geom_bar(aes(x = Position, y = Coverage, fill = Mutation), stat = "identity", position = "dodge") +
  scale_color_d3(palette = "category20") +
  theme_classic() + theme(legend.position = "top") + guides(fill = guide_legend(ncol = 1)) + labs(title = "Intronic coverage")

p4 <- ggplot(count_details %>% dplyr::filter(Category == "Intronic", Position <= 35)) +
  geom_bar(aes(x = Position, y = Mismatches, fill = Mutation), stat = "identity", position = "dodge") +
  scale_color_d3(palette = "category20") +
  theme_classic() + theme(legend.position = "top") + guides(fill = guide_legend(ncol = 1)) + labs(title = "Intronic mismatches")

p <- (p1 / p2 / p3 / p4) + plot_layout(guides = "collect")
p_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview/20250321_LDJ_embryo.count_details.pdf")
ggsave(p, filename = p_save_to, width = 15, height = 7)


# Percentage of exonic and intronic reads
p_tab <- count_details %>%
  dplyr::filter(Category %in% c("Intronic", "Exonic")) %>%
  dplyr::group_by(Position, Category) %>%
  dplyr::mutate(Percentages = Mismatches / sum(Mismatches) * 100)

p1 <- ggplot(p_tab %>% dplyr::filter(Category == "Exonic")) +
  geom_bar(aes(x = Position, y = Percentages, fill = Mutation, group = Mutation), position = "stack", stat = "identity") +
  scale_color_d3(palette = "category20") +
  theme_classic() + theme(legend.position = "top") + guides(fill = guide_legend(ncol = 1)) + labs(title = "Exonic mismatche percentages")

p2 <- ggplot(p_tab %>% dplyr::filter(Category == "Intronic")) +
  geom_bar(aes(x = Position, y = Percentages, fill = Mutation, group = Mutation), position = "stack", stat = "identity") +
  scale_color_d3(palette = "category20") +
  theme_classic() + theme(legend.position = "top") + guides(fill = guide_legend(ncol = 1)) + labs(title = "Exonic mismatche percentages")
p <- (p1 / p2) + plot_layout(guides = "collect")
p_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview/20250321_LDJ_embryo.mismatch_percentages.pdf")
ggsave(p, filename = p_save_to, width = 15, height = 4)


# Count details
working_dir <- file.path(project_dir, "outputs/analysis/geo_slam/problems/smart_seq_3prime_utr")
cds_count_detail <- file.path(working_dir, "cds/cds_reads.nascent_rna.mismatchdetails.tsv") %>% fread() %>% dplyr::mutate(Read_source = "CDS")
utr_count_detail <- file.path(working_dir, "utr/utr_reads.nascent_rna.mismatchdetails.tsv") %>% fread() %>% dplyr::mutate(Read_source = "UTR")
exon_count_detail <- file.path(working_dir, "exon/exon_reads.nascent_rna.mismatchdetails.tsv") %>% fread() %>% dplyr::mutate(Read_source = "Exon")
intron_count_detail <- file.path(working_dir, "intron/intron_reads.nascent_rna.mismatchdetails.tsv") %>% fread() %>% dplyr::mutate(Read_source = "Intron")
total_count_detail <- file.path(project_dir, "outputs/analysis/preprocessing/geo_slam/SLAM_20250321_byLDJ/nascent_rna/per_sample/lps_6h_r3/lps_6h_r3.mismatchdetails.tsv") %>% fread() %>% dplyr::mutate(Read_source = "Total")

count_details <- dplyr::bind_rows(cds_count_detail, utr_count_detail, exon_count_detail, intron_count_detail, total_count_detail) %>%
  dplyr::filter(Category %in% c("Exonic", "Intronic")) %>%
  dplyr::mutate(Mutation = paste0(Genomic, "->", Read), Mutation_rates = Mismatches / Coverage) %>%
  dplyr::mutate(Mutation = forcats::fct_reorder(Mutation, Mismatches)) %>%
  dplyr::mutate(Read_source = factor(Read_source, levels = c("Total", "UTR", "CDS", "Exon", "Intron")))

p <- ggplot(count_details) +
  geom_line(aes(x = Position, y = Mutation_rates, color = Mutation)) +
  facet_wrap(~ Read_source + Category, ncol = 2, scales = "free_y") +
  theme_classic()
p_save_to <- file.path(working_dir, "20250321_LDJ_kidey.count_details.by_mapped_regions.pdf")
ggsave(p, filename = p_save_to, width = 10, height = 10)


p1 <- ggplot(count_details %>% dplyr::filter(Read_source == "CDS", 20 < Position, Position <= 35)) +
  geom_bar(aes(x = Position, y = Coverage, fill = Mutation), stat = "identity", position = "dodge") +
  facet_wrap(~ Category, ncol = 2, scales = "free_y") +
  theme_classic()
p2 <- ggplot(count_details %>% dplyr::filter(Read_source == "CDS", 20 < Position, Position <= 30)) +
  geom_bar(aes(x = Position, y = Mismatches, fill = Mutation), stat = "identity", position = "dodge") +
  facet_wrap(~ Category, ncol = 2, scales = "free_y") +
  theme_classic()
p <- p1 / p2
p_save_to <- file.path(working_dir, "20250321_LDJ_kidey.count_details.by_mapped_regions.cds.pdf")
ggsave(p, filename = p_save_to, width = 20, height = 10)

p1 <- ggplot(count_details %>% dplyr::filter(Read_source == "UTR", Position > 20, Position <= 30)) +
  geom_bar(aes(x = Position, y = Coverage, fill = Mutation), stat = "identity", position = "dodge") +
  facet_wrap(~ Category, ncol = 2, scales = "free_y") +
  theme_classic()
p2 <- ggplot(count_details %>% dplyr::filter(Read_source == "UTR", Position > 20, Position <= 30)) +
  geom_bar(aes(x = Position, y = Mismatches, fill = Mutation), stat = "identity", position = "dodge") +
  facet_wrap(~ Category, ncol = 2, scales = "free_y") +
  theme_classic()
p <- (p1 / p2) + plot_layout(guides = "collect")
p_save_to <- file.path(working_dir, "20250321_LDJ_kidey.count_details.by_mapped_regions.utr.pdf")
ggsave(p, filename = p_save_to, width = 20, height = 10)

p1 <- ggplot(count_details %>% dplyr::filter(Read_source == "Total", Position > 20, Position < 30)) +
  geom_bar(aes(x = Position, y = Coverage, fill = Mutation), stat = "identity", position = "dodge") +
  facet_wrap(~ Category, ncol = 2, scales = "free_y") +
  theme_classic()
p2 <- ggplot(count_details %>% dplyr::filter(Read_source == "Total", Position > 20, Position < 30)) +
  geom_bar(aes(x = Position, y = Mismatches, fill = Mutation), stat = "identity", position = "dodge") +
  facet_wrap(~ Category, ncol = 2, scales = "free_y") +
  theme_classic()
p <- (p1 / p2) + plot_layout(guides = "collect")
p_save_to <- file.path(working_dir, "20250321_LDJ_kidey.count_details.by_mapped_regions.total.pdf")
ggsave(p, filename = p_save_to, width = 20, height = 10)

p1 <- ggplot(count_details %>% dplyr::filter(Read_source == "Exon", Position > 20, Position < 30)) +
  geom_bar(aes(x = Position, y = Coverage, fill = Mutation), stat = "identity", position = "dodge") +
  facet_wrap(~ Category, ncol = 2, scales = "free_y") +
  theme_classic()
p2 <- ggplot(count_details %>% dplyr::filter(Read_source == "Exon", Position > 20, Position < 30)) +
  geom_bar(aes(x = Position, y = Mismatches, fill = Mutation), stat = "identity", position = "dodge") +
  facet_wrap(~ Category, ncol = 2, scales = "free_y") +
  theme_classic()
p <- (p1 / p2) + plot_layout(guides = "collect")
p_save_to <- file.path(working_dir, "20250321_LDJ_kidey.count_details.by_mapped_regions.exon.pdf")
ggsave(p, filename = p_save_to, width = 20, height = 10)


# 20250419 Smart-seqV2
# 1. Barcodes per sample
batch_id <- "GEO_SLAM_20250419_smartseq2"
batch_id <- "GEO_SLAM_20250420_vasaseq"

in_file <- file.path(project_dir, "outputs/analysis/preprocessing/geo_slam", batch_id, "alignment", paste0(batch_id, ".Solo.out"), "GeneFull/UMIperCellSorted.txt")
bc_tab <- fread(in_file) %>% dplyr::mutate(Rank = 1:n()) %>% dplyr::select(Counts = V1, Rank)
p1 <- ggplot(bc_tab) + geom_point(aes(x = Rank, y = Counts)) +
  geom_text_repel(aes(x = Rank, y = Counts, label=Rank), data = bc_tab %>% head(6), min.segment.length = 0) +
  scale_y_log10() +
  theme_classic() +
  labs(x = "Rank by UMI counts", y = "UMI counts")
p1_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview", paste0(batch_id, ".umi_counts_by_sample.rank_ordered.pdf"))
ggsave(p1_save_to, plot = p1, width = 6, height = 5)

bc_tab_top6 <- bc_tab %>% head(6)
p2 <- ggplot(bc_tab_top6) +
  geom_bar(aes(x = 1, y = Counts, fill = as.factor(Rank)), stat = "identity", position = "stack") +
  scale_fill_npg() +
  labs(fill = "Rank", x = NULL) +
  coord_polar('y') +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p2_save_to <- file.path(project_dir, "outputs/analysis/geo_slam/overview", paste0(batch_id, ".umi_counts_by_sample.top_6.pdf"))
ggsave(p2_save_to, plot = p2, width = 7, height = 7)

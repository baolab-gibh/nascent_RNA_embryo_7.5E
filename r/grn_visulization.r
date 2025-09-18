#!/usr/bin/env Rscript
# File: grn_visulization.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Mar 18, 2025
# Updated:


options(bitmapType = "cairo")
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  
  library(Seurat)
  library(patchwork)
  library(ggsci)
  library(ggrepel)

  library(RColorBrewer)
  library(ComplexHeatmap)
  library(dendextend)
  library(circlize)
  library(clusterProfiler)
  library(simplifyEnrichment)
  library(enrichplot)
})


sigmoid_rev <- function(x) -log(1 / x - 1) # Reverse of sigmoid
sigmoid <- function(x) 1 / (1 + exp(-x)) # Sigmoid function


rbgc <- brewer.pal(n = 11, name = "RdBu")
version <- "version_3"

project_dir <- "~/Documents/projects/wp_vasaseq"

regions_l0 <- c("Ectoderm", "Ectoderm", "Endoderm", "Endoderm", "Ectoderm", "Ectoderm", "Mesoderm", "Mesoderm") %>%
  purrr::set_names(c('A', 'P', 'EA', 'EP', 'L', 'R', 'MA', 'MP'))

all_batches <- c("240409_Lib_embryo", "240612_Lib_28region", "240620_Lib_38region", "240703_Lib_32region", "240710_Lib_37region", "240717_Lib_28region")



#
## Regulon activities
#
read_source <- "total_reads"
read_source <- "new_reads"


# Regulon specificity score rank plot
rss_tab <- file.path(project_dir, 'outputs/analysis/regulon', paste0('all_batches.', read_source, '.rss_regions_l0.csv')) %>%
  fread() %>%
  dplyr::select(Regions_l0 = V1, dplyr::everything()) %>%
  tidyr::pivot_longer(-Regions_l0, values_to = "RSS", names_to = "Regulon") %>%
  dplyr::group_by(Regions_l0) %>%
  dplyr::arrange(-RSS) %>%
  dplyr::mutate(Regulon_rank = dplyr::row_number())
lbl_tab <- rss_tab %>% dplyr::filter(Regulon_rank <= 5)

p <- ggplot() +
  geom_point(aes(x = Regulon_rank, y = RSS), data = rss_tab, color = "black", size = 0.5) +
  geom_point(aes(x = Regulon_rank, y = RSS, color = Regions_l0), data = lbl_tab, size = 3, alpha = 0.5) +
  geom_text_repel(aes(x = Regulon_rank, y = RSS, label = Regulon), data = lbl_tab, max.overlaps = Inf, min.segment.length = 0) +
  facet_wrap(~Regions_l0, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none")
psaveto <- file.path(project_dir, "outputs/analysis/regulon/plots", paste0('all_batches.', read_source, '.rss_regions_l0.regulon_spcificity_rank.pdf'))
ggsave(psaveto, plot = p, width = 9, height = 5)


# AUC z-score heatmap, top 50 per regulon
auc_mtx_zscore <- file.path(project_dir, 'outputs/analysis/regulon', paste0('all_batches.', read_source, '.auc_mtx_zscore.rss_regions_l0.top_50.csv')) %>%
  fread() %>%
  as.data.frame() %>%
  tibble::column_to_rownames("Cell") %>%
  as.matrix() %>%
  t()

cls_mtd <- "ward.D2"
row_cls_k <- 10

row_hcl <- stats::hclust(stats::dist(auc_mtx_zscore), method = cls_mtd)
hcl_tree <- stats::cutree(row_hcl, k = row_cls_k)

sample_info <- auc_mtx_zscore %>% colnames() %>% stringr::str_split(pattern = "-", simplify = TRUE) %>% as.data.frame() %>% 
  dplyr::select(Sample = V1, Batch = V2) %>%
  dplyr::mutate(Cell_id = paste0(Sample, "-", Batch)) %>%
  dplyr::mutate(Regions_l1 = stringr::str_extract(Sample, pattern = "^[0-9]+([ALRPEM]+)_[ATCG]+$", group = 1)) %>%
  dplyr::mutate(Regions_l0 = regions_l0[Regions_l1]) %>%
  dplyr::mutate(Layer = stringr::str_extract(Sample, pattern = "^([0-9]+)[ALRPEM]+_[ATCG]+$", group = 1) %>% as.integer())

selected_samples <- sample_info %>% dplyr::pull(Regions_l0, Cell_id)

# row_dend <- as.dendrogram(row_hcl, k = row_cls_k)
# group_labels <- data.frame(Labels = labels(row_dend), Cluster = hcl_tree[labels(row_dend)]) %>% dplyr::pull(Cluster) %>% unique() %>% as.character()
# row_dend <- row_dend %>% color_branches(row_dend, k = row_cls_k, groupLabels = group_labels)
# col_dend <- as.dendrogram(hclust(dist(t(auc_mtx_zscore)), method = cls_mtk), k = row_cls_k)
col_dend <- cluster_within_group(auc_mtx_zscore, selected_samples[colnames(auc_mtx_zscore)])

egg_list <- data.frame(Regulon = names(hcl_tree), Regulon_cluster = as.character(hcl_tree)) %>%
  dplyr::mutate(TF_name = stringr::str_extract(Regulon, "([A-Za-z0-9]+?)\\(\\+\\)", group = 1)) %>%
  dplyr::group_by(Regulon_cluster) %>%
  dplyr::summarize(Region_TF = list(TF_name)) %>%
  dplyr::pull(Region_TF, Regulon_cluster) %>%
  lapply(function(x) enrichGO(x, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "ALL"))

label_texts <- egg_list %>% lapply(function(x) {
  nr_genes <- x@gene %>% length %>% `/`(5) %>% ceiling()
  texts <- dplyr::filter(x@result, ONTOLOGY == "BP") %>% dplyr::slice_min(p.adjust, n = nr_genes, with_ties = FALSE) %>% dplyr::pull(ID)
  if (length(texts) == 0) {
    return("None")
  }
  return(texts)
})


boundary_val <- 2.
auc_mtx_zscore_trunc <- auc_mtx_zscore %>% (function(x) { x[x < -boundary_val] <- -boundary_val; x[x > +boundary_val] <- +boundary_val; x })
color_fun <- colorRamp2(c(min(auc_mtx_zscore_trunc), 0, max(auc_mtx_zscore_trunc)), c(rbgc[11], rbgc[6], rbgc[1]))

set.seed(31415926)
file.path(project_dir, "outputs/analysis/regulon/plots", paste0("all_batches.", read_source, ".regulon_activities.heatmap.pdf")) %>%
  pdf(width = 5, height = 7)
top_ann <- HeatmapAnnotation(Regions = selected_samples, annotation_name_side = "left", annotation_legend_param = list(Regions = list(direction = "horizontal"))) # Column annotations
right_ann <- rowAnnotation(textbox = anno_textbox(hcl_tree, label_texts, gp = gpar(fontsize = 7, color = "black"),  max_width = unit(35, "mm"), word_wrap = TRUE, add_new_line = TRUE))
# right_ann <- rowAnnotation(foo = anno_mark(at = c(199), labels = rownames(rgl_auc_tab)[c(199)], labels_gp = gpar(fontsize = 8))) # Row annotation to highligh regulons

Heatmap(auc_mtx_zscore_trunc, name = "Reg. Act.", col = color_fun,
  row_split = hcl_tree,
  # cluster_rows = row_dend,
  # clustering_method_rows = cls_mtd,
  cluster_columns = col_dend,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = top_ann,
  right_annotation = right_ann,
  heatmap_legend_param = list(direction = "horizontal")
) %>% draw(merge_legend = TRUE, heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off()


# Functional enrichment
for (which_cluster in names(egg_list)) {
  sub_egg <- pairwise_termsim(egg_list[[which_cluster]])
  p <- emapplot(sub_egg, showCategory = 20, layout = "kk", cex_category = 1.5)
  emap_save_to <- file.path(project_dir, "outputs/analysis/regulon/plots", paste0("all_batches.", read_source, ".regulon_cluster.go_erichemnt.cluster_", which_cluster, ".pdf"))
  ggsave(emap_save_to, plot = p, width = 8, heigh = 8)
}


# Plotting corn plots
corn_plot_tab <- file.path(project_dir, "outputs/analysis/overview/corn_map.samples_per_region.csv") %>% fread() %>%
  dplyr::mutate(N_samples = as.factor(N_samples))


x_tick_pos <- c(-10.490, -7.824, -4.863, -1.605, 1.605, 4.863, 7.824, 10.490)
x_tick_lab <- c("EA", "MA", "A", "L", "R", "P", "MP", "EP")
p <- ggplot(corn_plot_tab) +
  geom_point(aes(x = region_axis, y = layer_axis, color = N_samples), size = 10) +
  scale_x_continuous(breaks = x_tick_pos, labels = x_tick_lab, expand = c(.1, .1), position = "top") +
  scale_y_continuous(breaks = 1:17) +
  scale_color_npg() +
  labs(x = "Region", y = "Layer", color = "N Samples") +
  theme_classic() +
  theme(axis.ticks.x.bottom = element_blank(), axis.text.x.bottom = element_blank())
psaveto <- file.path(project_dir, "outputs/analysis/overview/corn_map.samples_per_region.pdf")
ggsave(psaveto, plot = p, width = 4.25, height = 5.5)

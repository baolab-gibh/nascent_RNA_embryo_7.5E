#!/usr/bin/env Rscript

#
## Install packages
#
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::version()
# # If your bioconductor version is previous to 4.0, see the section bellow
# 
# ## Required
# BiocManager::install(c("AUCell", "RcisTarget"))
# BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost
# 
# ## Optional (but highly recommended):
# # To score the network on cells (i.e. run AUCell):
# BiocManager::install(c("zoo", "mixtools"))
# devtools::install_version("rbokeh", version = '0.5.2')
# 
# # For various visualizations and perform t-SNEs:
# BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# 
# # To support paralell execution (not available in Windows):
# BiocManager::install(c("doMC", "doRNG"))
# 
# # To export/visualize in http://scope.aertslab.org
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
# devtools::install_github("aertslab/SCENIC@v1.1.2")
# BiocManager::install("Singlecellexperiment")
# devtools::install_version("Seurat", version = '4.2.0')


options(bitmapType = "cairo")
suppressPackageStartupMessages({
  library(GENIE3)
  library(RcisTarget)
  library(AUCell)
  library(SCENIC)
  
  library(Seurat)
  library(harmony)
  library(tidyverse)
  library(data.table)
  library(SingleCellExperiment)
  library(SCopeLoomR)
  library(patchwork)
  library(ggsci)

  library(RColorBrewer)
  library(ComplexHeatmap)
  library(dendextend)
  library(circlize)
  library(clusterProfiler)
  # library(DOSE)
  library(simplifyEnrichment)
  library(enrichplot)
})


sigmoid_rev <- function(x) -log(1 / x - 1) # Reverse of sigmoid
sigmoid <- function(x) 1 / (1 + exp(-x)) # Sigmoid function



# runSCENIC_4_aucell_binarize in-house version
rm(list=ls()); gc()

rbgc <- brewer.pal(n = 11, name = "RdBu")
version <- "version_3"

project_dir <- "~/Documents/projects/wp_vasaseq"

regions_l0 <- c("Ectoderm", "Ectoderm", "Endoderm", "Endoderm", "Ectoderm", "Ectoderm", "Mesoderm", "Mesoderm") %>%
  purrr::set_names(c('A', 'P', 'EA', 'EP', 'L', 'R', 'MA', 'MP'))


all_batches <- c("240409_Lib_embryo", "240612_Lib_28region", "240620_Lib_38region", "240703_Lib_32region", "240710_Lib_37region", "240717_Lib_28region")
embryo_saveto <- file.path(project_dir, "outputs/analysis/gene_expression", version, "mouse_embryo.rds")
if (file.exists(embryo_saveto)) {
  cat("[I]: Loading from the disk ...\n")
  embryo <- readRDS(embryo_saveto)
} else {
  cat("[I]: Dumping into the disk ...\n")
  obj_list <- lapply(all_batches, function(per_batch) {
    obj <- file.path(project_dir, "outputs/analysis/preprocessing/geo_seq/quantification", per_batch, "featureCounts/10X") %>%
      Read10X() %>%
      CreateSeuratObject(project = per_batch, min.cell = 10, min.features = 100)

    obj[["Batches"]] <- per_batch
    obj[["Regions"]] <- colnames(obj) %>% stringr::str_remove("_[ATCG]+$")
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 10, features = VariableFeatures(object = obj), verbose = FALSE)

    # Cluster by UMAP using gene expression
    tmp_obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
    tmp_obj <- FindClusters(tmp_obj, resolution = 0.5, verbose = FALSE)
    tmp_obj <- RunUMAP(tmp_obj, reduction = "pca", dims = 1:10, n.neighbors = 25, verbose = FALSE)
    # tmp_obj <- RunTSNE(tmp_obj, reduction = "pca", dims = 1:10, verbose = FALSE) # Perplexity is too large
    p <- DimPlot(tmp_obj, reduction = "umap", label = TRUE, group.by = "Regions", pt.size = 10) + NoLegend()
    umap_save_to <- file.path(project_dir, "outputs/analysis/gene_expression", version, paste0(per_batch, ".umap.pdf"))
    ggsave(p, filename = umap_save_to, width = 8.5, height = 6)

    obj
  }) %>%
    purrr::set_names(all_batches)

  rm("embryo"); gc()
  set.seed(3141592)
  embryo <- merge(obj_list[[1]], obj_list[2:length(obj_list)]) %>%
    (function(obj) {
      tar_cells <- colnames(obj) %>% purrr::discard(~stringr::str_detect(.x, "NC_"))
      obj[, tar_cells]
    }) %>%
    NormalizeData(normalization.method = "LogNormalize")

  VariableFeatures(embryo) <- obj_list %>%
    lapply(function(obj) {
      FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3500, verbose = FALSE) %>% VariableFeatures()
    }) %>%
    Reduce(intersect, .) %>% unique()

  embryo <- embryo %>% ScaleData(verbose = FALSE) %>%
    RunPCA(features = VariableFeatures(embryo), npcs = 15, verbose = FALSE) %>%
    RunHarmony("Batches")

  embryo <- FindNeighbors(embryo, reduction = "harmony", dims = 1:15, verbose = FALSE)
  embryo <- FindClusters(embryo, resolution = 0.5, verbose = FALSE)
  embryo <- RunUMAP(embryo, reduction = "harmony", dims = 1:15, verbose = FALSE)
  embryo <- RunTSNE(embryo, reduction = "harmony", dims = 1:15, verbose = FALSE) # Perplexity is too large

  embryo[["Treatments"]] <- ifelse(stringr::str_detect(unlist(embryo[["Regions"]]), "NC_"), "No4sU", "4sU")
  embryo[["Layers"]] <- embryo[["Regions"]] %>% unlist() %>% stringr::str_extract("[0-9]+") %>% as.integer() %>% purrr::imap(~ceiling(.x / 3) - (.x > 15)) %>% unlist()
  embryo[["Regions_l1"]] <- embryo[["Regions"]] %>% unlist() %>% stringr::str_extract("[ALMEPR]+$")
  embryo[["Regions_l0"]] <- embryo[["Regions"]] %>% unlist() %>% stringr::str_extract("[0-9]([A-Z]+)$", group = 1) %>% `[`(regions_l0, .) %>% purrr::set_names(NULL)

  p1 <- DimPlot(embryo, reduction = "umap", group.by = "Batches", pt.size = 2.5) + scale_color_jama() + labs(x = NULL, color = "Batches") # Check batch effects
  # p2 <- DimPlot(embryo, reduction = "umap", group.by = "Treatments", pt.size = 2.5) + scale_color_npg() + labs(x = NULL, y = NULL, color = "Treatments") # Check control vs treatment
  p2 <- DimPlot(embryo, reduction = "umap", group.by = "Layers", pt.size = 2.5) + labs(x = NULL, y = NULL, color = "Layers") # Check layers
  p3 <- DimPlot(embryo, reduction = "umap", group.by = "Regions_l1", pt.size = 2.5) + scale_color_igv() + labs(color = "Regions L1") # Check regions L1
  p4 <- DimPlot(embryo, reduction = "umap", group.by = "Regions_l0", pt.size = 2.5) + labs(y = NULL) + labs(color = "Regions L0") # Check regions L0

  p <- (p1 + p2 + p3 + p4) + plot_layout(guides = "collect")
  umap_save_to <- file.path(project_dir, "outputs/analysis/gene_expression", version, "all_batches.umap.check_groups.pdf")
  ggsave(p, filename = umap_save_to, width = 10.5, heigh = 8)

  p <- DimPlot(embryo, reduction = "tsne", group.by = "Regions_l1", pt.size = 2.5)
  umap_save_to <- file.path(project_dir, "outputs/analysis/gene_expression", version, "all_batches.tsne.color_by_region_l1.pdf")
  ggsave(p, filename = umap_save_to, width = 7.5, height = 6)

  # Save
  saveRDS(embryo, embryo_saveto)
}



#
## Regulons
#
sce <- embryo %>% as.SingleCellExperiment()
cell_info <- colData(sce)
loom_obj_save_to <- file.path(project_dir, "outputs/analysis/regulon", version, "mouse_embryo.loom")
if (file.exists(loom_obj_save_to)) {
  cat("[I]: Loading from the disk ...\n")
  loom <- open_loom(loom_obj_save_to)
  expr_mat <- get_dgem(loom)
  close_loom(loom)
} else {
  cat("[I]: Dumping into the disk ...\n")
  expr_mat <- logcounts(sce)
  loom <- build_loom(loom_obj_save_to, dgem=expr_mat)
  loom <- add_embedding(loom, embryo@reductions$umap@cell.embeddings, "umap")
  loom <- add_embedding(loom, embryo@reductions$pca@cell.embeddings, "pca")
  loom <- add_row_attr(loom, "cellInfo", rownames(cell_info))
  loom <- add_row_attr(loom, "colVars", cell_info$Regions)
  loom <- add_row_attr(loom, "Regions_l1", cell_info$Regions_l1)
  loom <- add_row_attr(loom, "Treatments", cell_info$Treatments)
  loom <- add_row_attr(loom, "Layers", cell_info$Layers)
  expr_mat <- get_dgem(loom)
  close_loom(loom)
}


# SCENIC seetings
setwd(file.path(project_dir, "outputs/analysis/regulon", version)); getwd()

re_init <- FALSE
scenic_options_save_to <- file.path(project_dir, "outputs/analysis/regulon", version, "scenic_options.rds")
if (re_init) {
  org <- "mgi"
  dataset_title <- "Mouse embryo"
  db_dir <- file.path(project_dir, "outputs/analysis/regulon/cis_target_db")
  dbs <- c("500bp" = "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather", "10kb" = "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
  scenicOptions <- initializeScenic(org = org, dbDir = db_dir, dbs = dbs, datasetTitle = dataset_title, nCores = 3)
  scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
  scenicOptions@settings$seed <- 31415926
  scenicOptions@settings$defaultTsne$perpl <- 20
  scenicOptions@settings$defaultTsne$dims <- 20
  saveRDS(cell_info, file=getDatasetInfo(scenicOptions, "cellInfo"))
  saveRDS(scenicOptions, file = scenic_options_save_to)

  # Co-expression network
  gene_kept <- geneFiltering(expr_mat, scenicOptions, 3 * .01 * ncol(expr_mat), ncol(expr_mat) *.01)
  expr_mat_flt <- expr_mat[gene_kept, ]
  runCorrelation(expr_mat_flt, scenicOptions)
  runGenie3(expr_mat_flt, scenicOptions)

  # Build and score the GRN
  expr_mat_log <- expr_mat
  runSCENIC_1_coexNetwork2modules(scenicOptions)
  runSCENIC_2_createRegulons(scenicOptions)
  runSCENIC_3_scoreCells(scenicOptions, expr_mat_log)
  runSCENIC_4_aucell_binarize(scenicOptions, exprMat = expr_mat_log)

  # Save the results into loom
  scenicOptions@fileNames$output["loomFile",] <- "output/mouse_embryo_SCENIC.loom"
  export2scope(scenicOptions, expr_mat)

  # Save the session again
  saveRDS(scenicOptions, scenic_options_save_to)
} else {
  scenicOptions <- readRDS(scenic_options_save_to)
  scenicLoom <- open_loom(scenicOptions@fileNames$output["loomFile",])
}


# t-SNE to show regulon activities
tSNE_fileName <- tsneAUC(scenicOptions, aucType = "Binary", filePrefix = getSettings(scenicOptions, "tSNE_filePrefix"), onlyHighConf = FALSE)
tSNE <- readRDS(tSNE_fileName)

sub <- paste0("t-SNE on ", tSNE$type)
tsne_tbl <- as.data.frame(tSNE$Y) %>% tibble::rownames_to_column("cell") %>%
  dplyr::left_join(cell_info %>% as.data.frame() %>% tibble::rownames_to_column("cell"), by = "cell") %>%
  dplyr::select(cell, tsne1, tsne2, Regions) %>%
  dplyr::mutate(by_layers = stringr::str_extract(Regions, "[0-9]+")) %>%
  dplyr::mutate(by_regions = stringr::str_extract(Regions, "[A-Za-z]+"))

p <- ggplot(tsne_tbl) +
  geom_point(aes(x = tsne1, y = tsne2, color = Regions), size = 10, alpha = 0.75) +
  geom_text(aes(x = tsne1, y = tsne2, label = Regions), size = 3) +
  labs(title = sub) +
  theme_classic() +
  theme(legend.position = "none")
tsne_save_to <- file.path(project_dir, "outputs/analysis/regulon", version, "tsne_embryo.all_samples.pdf")
ggsave(p, filename = tsne_save_to, width = 9, height = 8)

p <- ggplot(tsne_tbl) +
  geom_point(aes(x = tsne1, y = tsne2, color = by_layers), size = 10, alpha = 0.75) +
  geom_text(aes(x = tsne1, y = tsne2, label = Regions), size = 3) +
  labs(title = sub) +
  theme_classic() +
  theme(legend.position = "none")
tsne_save_to <- file.path(project_dir, "outputs/analysis/regulon", version, "tsne_embryo.all_samples.color_by_layers.pdf")
ggsave(p, filename = tsne_save_to, width = 9, height = 8)

p <- ggplot(tsne_tbl) +
  geom_point(aes(x = tsne1, y = tsne2, color = by_regions), size = 10, alpha = 0.75) +
  geom_text(aes(x = tsne1, y = tsne2, label = Regions), size = 3) +
  labs(title = sub) +
  theme_classic() +
  theme(legend.position = "none")
tsne_save_to <- file.path(project_dir, "outputs/analysis/regulon", version, "tsne_embryo.all_samples.color_by_regions.pdf")
ggsave(p, filename = tsne_save_to, width = 9, height = 8)


#
## Regulon activities
#
loom_file <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(loom_file)

regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)


cell_info <- readRDS(scenicOptions@inputDatasetInfo$cellInfo)
selected_samples <- cell_info %>% as.data.frame() %>%
  tibble::rownames_to_column("Cell_barcodes") %>%
    dplyr::mutate(Regions_l0 = dplyr::case_when(
      Regions_l1 %in% c("L", "A", "P", "R") ~ "Ectoderm",
      Regions_l1 %in% c("MA", "MP") ~ "Mesoderm",
      Regions_l1 %in% c("EA", "EP") ~ "Endoderm"
  )) %>%
  dplyr::mutate(Regions_l0 = factor(Regions_l0, levels = c("Ectoderm", "Mesoderm", "Endoderm"))) %>%
  dplyr::filter(!stringr::str_detect(Cell_barcodes, "NC_")) %>%
  dplyr::arrange(Regions_l0, Regions_l1) %>%
  dplyr::pull(Regions_l0, Cell_barcodes)

rgl_auc_subset <- regulonsAUC[onlyNonDuplicatedExtended(rownames(regulonsAUC)), names(selected_samples)] %>%
  getAUC() %>% t() %>% scale(scale = TRUE, center = TRUE) %>% t()
rgl_auc_subset[rgl_auc_subset >= 2] <- 2
rgl_auc_subset[rgl_auc_subset <= -2] <- -2

color_fun <- colorRamp2(c(min(rgl_auc_subset), 0, max(rgl_auc_subset)), c(rbgc[11], rbgc[6], rbgc[1]))

cls_mtd <- "ward.D2"
row_cls_k <- 10

row_hcl <- stats::hclust(stats::dist(rgl_auc_subset), method = cls_mtd)
hcl_tree <- stats::cutree(row_hcl, k = row_cls_k)

# row_dend <- as.dendrogram(row_hcl, k = row_cls_k)
# group_labels <- data.frame(Labels = labels(row_dend), Cluster = hcl_tree[labels(row_dend)]) %>% dplyr::pull(Cluster) %>% unique() %>% as.character()
# row_dend <- row_dend %>% color_branches(row_dend, k = row_cls_k, groupLabels = group_labels)
# col_dend <- as.dendrogram(hclust(dist(t(rgl_auc_subset)), method = cls_mtk), k = row_cls_k)
col_dend <- cluster_within_group(rgl_auc_subset, selected_samples[colnames(rgl_auc_subset)])

egg_list <- data.frame(Regulon = names(hcl_tree), Regulon_cluster = as.character(hcl_tree)) %>%
  dplyr::mutate(TF_name = stringr::str_extract(Regulon, "([A-Za-z0-9]+?)[_ ]", group = 1)) %>%
  dplyr::group_by(Regulon_cluster) %>%
  dplyr::summarize(Region_TF = list(TF_name)) %>%
  dplyr::pull(Region_TF, Regulon_cluster) %>%
  lapply(function(x) enrichGO(x, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "ALL"))

label_texts <- egg_list %>% lapply(function(x) {
  nr_genes <- x@gene %>% length %>% `/`(5) %>% ceiling()
  x@result %>%
    dplyr::filter(ONTOLOGY == "BP") %>%
    dplyr::slice_min(p.adjust, n = nr_genes, with_ties = FALSE) %>%
    dplyr::pull(ID)
})


set.seed(31415926)
file.path(project_dir, "outputs/analysis/regulon", version, "plots/regulon_activities.heatmap.pdf") %>% pdf(width = 5, height = 8)
top_ann <- HeatmapAnnotation(Regions = selected_samples, annotation_name_side = "left", annotation_legend_param = list(Regions = list(direction = "horizontal"))) # Column annotations
right_ann <- rowAnnotation(textbox = anno_textbox(hcl_tree, label_texts, gp = gpar(fontsize = 7, color = "black"),  max_width = unit(35, "mm"), word_wrap = TRUE, add_new_line = TRUE))
# right_ann <- rowAnnotation(foo = anno_mark(at = c(199), labels = rownames(rgl_auc_subset)[c(199)], labels_gp = gpar(fontsize = 8))) # Row annotation to highligh regulons

Heatmap(rgl_auc_subset, name = "Reg. Act.", col = color_fun,
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
  emap_save_to <- file.path(project_dir, "outputs/analysis/regulon", version, "plots", paste0("regulon_cluster.go_erichemnt.cluster_", which_cluster, ".pdf"))
  ggsave(emap_save_to, plot = p, width = 12, heigh = 12)
}

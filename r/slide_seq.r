#!/usr/bin/env Rscript
# File: slide_seq.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 12, 2024
# Updated:
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggbreak)

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


#
nascent_rna_tab <- file.path(project_dir, "outputs/analysis/preprocessing/slide_seq/nascent_rna/20241022_encoded_with_polya_embryo/20241022_encoded_with_polya_embryo.nascent_rna.tsv.gz") %>%
  fread()
plot_tab <- nascent_rna_tab %>%
  dplyr::select(Gene, Symbol, dplyr::ends_with("Readcount")) %>%
  tidyr::pivot_longer(-c(Gene, Symbol), names_to = "Sample", values_to = "Readcount") %>%
  dplyr::mutate(Sample = stringr::str_remove(Sample, " Readcount")) %>%
  dplyr::mutate(Sample = factor(Sample, levels = c("167NC", "166OIL", "167TFEA")))

p <- ggplot() +
  geom_violin(data = plot_tab, aes(x = Sample, y = Readcount + 1, fill = Sample)) +
  geom_boxplot(data = plot_tab, aes(x = Sample, y = Readcount + 1), outlier.size = 0.50, width = 0.25) +
  scale_y_log10() +
  facet_wrap(~Sample, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "Readcount")
plot_save_to <- file.path(project_dir, "outputs/analysis/preprocessing/slide_seq/nascent_rna/plots/20241022_encoded_with_polya_embryo.nascent_rna.boxplot.pdf")
ggsave(p, filename = plot_save_to, width = 6, height = 3)


# Nascent RNA control
embryo <- file.path(project_dir, "outputs/analysis/preprocessing/slide_seq/read_alignments/20241022_encoded_with_polya_embryo/167NC/167NC.single_end.Solo.out/Gene/raw") %>%
  Read10X() %>%
  CreateSeuratObject(min.cell = 0, min.features = 0)

count_tab <- embryo[["RNA"]]@counts %>% as.data.frame() %>% rowSums()
count_tab <- file.path(project_dir, 'outputs/analysis/preprocessing/slide_seq/read_count/20241022_encoded_with_polya_embryo/167NC/167NC.txt') %>%
  fread() %>%
  dplyr::select(gene_id = Geneid, read_count = 7)

gene_biotype_tab <- file.path(project_dir, "outputs/references/genomic_features/mus_musculus.90.gene_biotypes.gene_id.tsv") %>%
  fread(header = FALSE, col.names = c("gene_id.ver", "biotype")) %>%
  dplyr::mutate(gene_id = stringr::str_extract(gene_id.ver, '(^ENSMUSG[0-9]+)[_.]', group = 1))
other_biotypes <- c(
  "polymorphic_pseudogene", "processed_pseudogene", "processed_transcript", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene",
  "transcribed_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "sense_overlapping", "ribozyme")

plot_tab <- count_tab %>%
  dplyr::left_join(gene_biotype_tab, by = "gene_id") %>%
  dplyr::filter(!biotype %in% c(other_biotypes, "Mt_tRNA", "Mt_rRNA"), read_count >= 3) %>%
  dplyr::group_by(biotype) %>%
  dplyr::summarize(n_counts = sum(read_count)) %>%
  dplyr::mutate(perc = n_counts / sum(n_counts) * 100) %>%
  dplyr::mutate(biotype = dplyr::if_else(is.na(biotype), "Others", biotype)) %>%
  dplyr::mutate(biotype = forcats::fct_reorder(biotype, perc))

p <- ggplot() +
  geom_bar(aes(x = biotype, y = perc, fill=biotype), stat = "identity", data = plot_tab) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_save_to <- file.path(project_dir, "outputs/analysis/preprocessing/slide_seq/read_alignments/20241022_encoded_with_polya_embryo/167NC/biotype_percentages.with_rRNA.pdf")
ggsave(plot_save_to, plot = p, width = 5, height = 5)


# T->C mutations per position on the reads
samples_id <- c("Sample161A", "Sample161B", "Sample163A", "Sample163B", "Sample168A", "Sample168B")

t2c_pos_tab <- file.path(project_dir, "outputs/analysis/preprocessing/slide_seq/nascent_rna/20241104_encoded_with_polya_embryo") %>%
  list.files(pattern = "*.mismatchdetails.tsv$", recursive = TRUE, full.names = TRUE) %>%
  purrr::discard(~stringr::str_detect(.x, "all_samples")) %>%
  lapply(function(p) {
    sample_id <- stringr::str_split(p, "/", simplify = TRUE) %>% purrr::keep(~stringr::str_detect(.x, "^Sample[0-9]+[AB]$"))
    fread(p) %>% dplyr::mutate(Sample = sample_id)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(Category = dplyr::case_when(
    Category %in% c("Exonic", "ExonicSense", "ExonicAntisense") ~ "Exonic",
    Category %in% c("Intronic", "IntronicSense", "IntronicAntisense") ~ "Intronic",
  )) %>%
  dplyr::group_by(Sample, Category, Genomic, Read, Position) %>%
  dplyr::summarize(Coverage = sum(Coverage), Mismatches = sum(Mismatches)) %>%
  dplyr::mutate(MismatchRatio = Mismatches / Coverage, Genomic = paste0(Category, "-", Genomic))

p <- ggplot(t2c_pos_tab) +
  geom_line(aes(x = Position, y = MismatchRatio, color = Read)) +
  facet_wrap(Sample ~ Genomic, scales = "free_y", ncol = 8) +
  theme_classic()
psaveto <- file.path(project_dir, "outputs/analysis/preprocessing/slide_seq/nascent_rna/20241104_encoded_with_polya_embryo/plots/t2c_mismatch_ratio.per_sample.pdf")
ggsave(psaveto, plot = p, width = 20, height = 9)


t2c_pos_tab <- file.path(project_dir, "outputs/analysis/preprocessing/slide_seq/nascent_rna/20241104_encoded_with_polya_embryo") %>%
  list.files(pattern = "*.mismatchdetails.tsv$", recursive = TRUE, full.names = TRUE) %>%
  purrr::keep(~stringr::str_detect(.x, "all_samples")) %>%
  fread() %>%
  dplyr::mutate(Category = dplyr::case_when(
    Category %in% c("Exonic", "ExonicSense", "ExonicAntisense") ~ "Exonic",
    Category %in% c("Intronic", "IntronicSense", "IntronicAntisense") ~ "Intronic",
  )) %>%
  dplyr::group_by(Category, Genomic, Read, Position) %>%
  dplyr::summarize(Coverage = sum(Coverage), Mismatches = sum(Mismatches)) %>%
  dplyr::mutate(MismatchRatio = Mismatches / Coverage, Genomic = paste0(Category, "-", Genomic))

p <- ggplot(t2c_pos_tab) +
  geom_line(aes(x = Position, y = MismatchRatio, color = Read)) +
  facet_wrap(. ~ Genomic, scales = "free_y", ncol = 8) +
  theme_classic()
psaveto <- file.path(project_dir, "outputs/analysis/preprocessing/slide_seq/nascent_rna/20241104_encoded_with_polya_embryo/plots/t2c_mismatch_ratio.all_sample.pdf")
ggsave(psaveto, plot = p, width = 20, height = 4)


# # AWK script, should be removed
# 
# gtf_file=~/Documents/projects/resources/Gencode/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M29/GRCm38_mapping/gencode.vM29lift38.annotation.gtf.gz
# awk -F$'\t' -f- <<'EOF' <(zcat $gtf_file) >| ~/mus_musculus.90.gene_biotypes.tsv
# $3 ~ /^gene$/ {
#   split($9, arr, ";");
#   for(i in arr) {
#     if (arr[i]~/gene_id/) {split(arr[i], name, "\""); gene_id=name[2]}
#     else if (arr[i]~/gene_type/) {split(arr[i],type,"\""); gene_type=type[2]}
#   }
#   print gene_id"\t"gene_type
#   next
# }
# EOF

#!/usr/bin/env Rscript
# File: plots.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jul 22, 2024
# Updated: Aug 06, 2024

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggsci)
  library(GenomicRanges)
  library(RColorBrewer)
  library(DESeq2)
})

project_dir <- "~/Documents/projects/wp_vasaseq"
rbgc <- brewer.pal(n = 11, name = "RdBu")

all_batches <- c("240409_Lib_embryo", "240612_Lib_28region", "240620_Lib_38region", "240703_Lib_32region", "240710_Lib_37region", "240717_Lib_28region")
regions_l0 <- c("Ectoderm", "Ectoderm", "Endoderm", "Endoderm", "Ectoderm", "Ectoderm", "Mesoderm", "Mesoderm", "Control") %>%
  purrr::set_names(c('A', 'P', 'EA', 'EP', 'L', 'R', 'MA', 'MP', "NC"))


# Statistics of nascent RNAs
newrna_stats <- file.path(project_dir, 'outputs/analysis/preprocessing/geo_seq/nascent_rna', all_batches) %>%
  lapply(function(x) list.files(x, pattern = "nascent_rna.ntrstat.tsv", full.names = TRUE, recursive = TRUE))

newrna_ntr_tab <- newrna_stats %>%
  lapply(function(x) {
    path_vec <- stringr::str_split(x, "/", simplify = TRUE)
    batch_id <- path_vec[length(path_vec) - 1]
    fread(x) %>% dplyr::mutate(Batch = batch_id)
  }) %>%
  Reduce(rbind, .) %>%
  dplyr::mutate(Type = factor(Type, levels = c("Intronic", "Exonic"))) %>%
  dplyr::mutate(Sample = paste0(Condition, "_", Batch)) %>%
  dplyr::mutate(Sample = forcats::fct_reorder2(Sample, Type, `T->C`)) %>%
  dplyr::mutate(Treatment = dplyr::if_else(stringr::str_detect(Condition, "^NC_"), "Control", "With_4sU")) %>%
  dplyr::mutate(Regions_l1 = stringr::str_extract(Condition, "^[0-9]{1,2}([A-Z]{1,2})_", group = 1) %>% purrr::map_chr(~dplyr::if_else(is.na(.x), "NC", .x))) %>%
  dplyr::mutate(Regions_l0 = regions_l0[Regions_l1])


# Nascent RNA ratio plot
ntr_table <- newrna_ntr_tab %>%
  dplyr::select(Type, Sample, Batch, ntr, ntr_upper, ntr_lower, Treatment, Regions_l1, Regions_l0) %>%
  tidyr::pivot_wider(names_from = "Type", values_from = c("ntr", "ntr_lower", "ntr_upper")) %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = c("With_4sU", "Control")))

ntr_plot <- ggplot() +
  geom_abline(linetype="dotted") +
  geom_point(aes(x = ntr_Exonic, y = ntr_Intronic, shape = Treatment, color = Regions_l0), alpha = 0.75, size = 4, data = ntr_table) +
  scale_color_npg() +
  labs(x = "Exonic reads", y = "Intronic reads", color = "Regions", shape = "Treatment") +
  lims(x = c(-0.01, 1.01), y = c(-0.01, 1.01)) +
  theme_classic() +
  theme(legend.position = "right")
ntr_plot_saveto <- file.path(project_dir, 'outputs/analysis/overview/nascent_rna.new_to_total_ratio.color_by_regions.pdf')
ggsave(ntr_plot_saveto, plot = ntr_plot, width = 7.5, height = 7)

ntr_plot <- ggplot() +
  geom_abline(linetype="dotted") +
  geom_point(aes(x = ntr_Exonic, y = ntr_Intronic, shape = Treatment, color = Batch), alpha = 0.75, size = 4, data = ntr_table) +
  scale_color_npg() +
  labs(x = "Exonic reads", y = "Intronic reads", color = "Batches", shape = "Treatment") +
  lims(x = c(-0.01, 1.01), y = c(-0.01, 1.01)) +
  theme_classic() +
  theme(legend.position = "right")
ntr_plot_saveto <- file.path(project_dir, 'outputs/analysis/overview/nascent_rna.new_to_total_ratio.color_by_Batch.pdf')
ggsave(ntr_plot_saveto, plot = ntr_plot, width = 8, height = 7)


# T->C mutations ratio
gmean_tab <- newrna_ntr_tab %>%
  dplyr::filter(!stringr::str_detect(Sample, "^NC_")) %>%
  dplyr::group_by(Type) %>%
  dplyr::summarize(geomic_mean = exp(sum(log(`T->C`)) / n())) %>%
  dplyr::mutate(gmean_label = paste0(as.character(round(geomic_mean, 4) * 100), "%"))

t2c_plot <- ggplot() +
  geom_point(aes(x = Sample, y = `T->C`, color = Type, shape = Treatment), size = 2, data=newrna_ntr_tab) +
  geom_path(aes(x = Sample, y = `T->C`, group = Sample), data=newrna_ntr_tab) +
  geom_hline(aes(yintercept = geomic_mean, color = Type), linetype = "dotted", data = gmean_tab) +
  geom_text(aes(y = geomic_mean, x = 5, label = gmean_label, color = Type), hjust = 0, size = 6, data = gmean_tab) +
  labs(x = "Sample", y = "T->C ratio (Grand-SLAM)", shape = "Treatment", color = "Genomic region") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top")
t2c_plot_saveto <- file.path(project_dir, 'outputs/analysis/overview/nascent_rna.t_to_c_ratio.pdf')
ggsave(t2c_plot_saveto, plot = t2c_plot, width = 15, height = 5)

t2c_plot <- ggplot() +
  geom_point(aes(x = Sample, y = `T->C`, color = Type, shape = Treatment), size = 2, data=newrna_ntr_tab) +
  geom_path(aes(x = Sample, y = `T->C`, group = Sample), data=newrna_ntr_tab) +
  geom_hline(aes(yintercept = geomic_mean, color = Type), linetype = "dotted", data = gmean_tab) +
  scale_color_npg() +
  # geom_text(aes(y = geomic_mean, x = 25, label = gmean_label, color = Type), hjust = 0, size = 3, data = gmean_tab) +
  facet_grid(~Batch, scales = "free_x") +
  labs(x = "Sample", y = "T->C ratio (Grand-SLAM)", shape = "Treatment", color = "Genomic region") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top")
t2c_plot_saveto <- file.path(project_dir, 'outputs/analysis/overview/nascent_rna.t_to_c_ratio.by_batch.pdf')
ggsave(t2c_plot_saveto, plot = t2c_plot, width = 15, height = 5)



# Gene body coverage
quality_reports <- file.path(project_dir, 'outputs/analysis/preprocessing/geo_seq/alignment_statistics', all_batches) %>%
  lapply(function(x) list.files(x, pattern = "coverage_profile_along_genes_", full.names = TRUE, recursive = TRUE))

coverage_tab <- quality_reports %>%
  unlist() %>%
  lapply(function(x) {
    path_vec <- stringr::str_split(x, "/", simplify = TRUE)
    category <- path_vec[length(path_vec)] %>% stringr::str_extract("total|low|high")
    sample_id <- path_vec[length(path_vec) - 2] %>% stringr::str_remove("_[ATCG]+$")
    batch_id <- path_vec[length(path_vec) - 3]

    fread(x) %>% dplyr::mutate(Batch = batch_id, Sample = sample_id, Category = category)
  }) %>%
  Reduce(rbind, .) %>%
  dplyr::select(Position = `#Transcript position`, Coverage = `Transcript coverage profile`, Batch, Sample, Category)

plot_tab <- coverage_tab %>%
  dplyr::filter(Category == "total") %>%
  dplyr::group_by(Batch, Sample) %>%
  dplyr::mutate(Percentage = Coverage / sum(Coverage) * 100) %>%
  dplyr::ungroup()

cov_plot_pan <- ggplot(plot_tab) +
  geom_line(aes(x = Position, y = Coverage, group = Sample, color = Batch)) +
  facet_wrap(~Batch, ncol = 2, scales = "free_y") +
  scale_color_npg() +
  labs(x = "Position (scale into [0, 100))", y = "Coverage") +
  theme_classic() +
  theme(legend.position = "none")

coverage_plot_saveto <- file.path(project_dir, "outputs/analysis/overview", "all_batches.total_genes.coverage_along_gene.pdf")
ggsave(cov_plot_pan, filename = coverage_plot_saveto, width = 10.5, height = 7)


# Genebody coverage 2019 Peng et.al. Nature
quality_reports <- file.path(project_dir, 'outputs/analysis/preprocessing/2019_peng_etal_nature/mapping_reports') %>%
  lapply(function(x) list.files(x, pattern = "coverage_profile_along_genes_", full.names = TRUE, recursive = TRUE))

coverage_tab <- quality_reports %>%
  unlist() %>%
  lapply(function(x) {
    path_vec <- stringr::str_split(x, "/", simplify = TRUE)
    category <- path_vec[length(path_vec)] %>% stringr::str_extract("total|low|high")
    sample_id <- path_vec[length(path_vec) - 2] %>% stringr::str_remove("_[ATCG]+$")

    fread(x) %>% dplyr::mutate(Sample = sample_id, Category = category)
  }) %>%
  Reduce(rbind, .) %>%
  dplyr::select(Position = `#Transcript position`, Coverage = `Transcript coverage profile`, Sample, Category)


plot_tab <- coverage_tab %>%
  dplyr::filter(Category == "total") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(Percentage = Coverage / sum(Coverage) * 100)

cov_plot_peng <- ggplot(plot_tab) +
  geom_line(aes(x = Position, y = Coverage, color = Sample)) +
  # scale_y_log10() +
  labs(x = "Position (scale into [0, 100))", y = "Coverage") +
  theme_classic() +
  theme(legend.position = "None")

coverage_plot_saveto <- file.path(project_dir, "outputs/analysis/overview/2019_peng_etal_nature.total_genes.coverage_along_gene.pdf")
ggsave(cov_plot_peng, filename = coverage_plot_saveto, width = 8, height = 6)


# Exon/intron/intergenic reads proportions
reads_dist_tab <- file.path(project_dir, "outputs/analysis/preprocessing/alignment_statistics/rnaseq_qc_results.all.txt") %>%
  fread() %>%
  dplyr::select(Batch = V1, Sample = V2, Region = V3, Counts = V4) %>%
  dplyr::mutate(Regions_l1 = stringr::str_extract(Sample, "^[0-9]+([A-Z]+)_", group = 1)) %>%
  dplyr::mutate(Regions_l1 = purrr::map_chr(Regions_l1, ~dplyr::if_else(is.na(.x), "NC", .x))) %>%
  dplyr::mutate(Regions_l0 = regions_l0[Regions_l1]) %>%
  dplyr::mutate(Treatment = dplyr::if_else(stringr::str_detect(Sample, "^NC_"), "Control", "With_4sU")) %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = c("With_4sU", "Control"))) %>%
  tidyr::pivot_wider(names_from = Region, values_from = Counts) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(total = exonic + intronic + intergenic + overlapping_exon) %>%
  dplyr::mutate(dplyr::across(c(exonic, intronic, intergenic, overlapping_exon), ~.x / total * 100, .names = "{.col}_perc"))

reads_dist_plot <- ggplot() +
  geom_point(aes(x = exonic_perc, y = intronic_perc, size = intergenic_perc, color = overlapping_exon_perc, shape = Treatment), data = reads_dist_tab, alpha = 0.75) +
  scale_color_gradient(low = rbgc[5], high = rbgc[1]) +
  labs(x = "Reads on exons (%)", y = "Reads on introns (%)", color = "Reads on SJ (%)", size = "Other reads (%)") +
  theme_classic() +
  theme(legend.position = "right")

reads_dist_plot_saveto <- file.path(project_dir, "outputs/analysis/overview/alignment_statistics.reads_distributions.pdf")
ggsave(reads_dist_plot_saveto, plot = reads_dist_plot, width = 8.5, height = 7)

reads_dist_plot_by_region <- reads_dist_plot + facet_wrap(~Regions_l0)
reads_dist_plot_by_region_saveto <- file.path(project_dir, "outputs/analysis/overview/alignment_statistics.reads_distributions.by_regions.pdf")
ggsave(reads_dist_plot_by_region_saveto, plot = reads_dist_plot_by_region, width = 8.5, height = 7)

reads_dist_plot_by_batch <- reads_dist_plot + facet_wrap(~Batch)
reads_dist_plot_by_batch_saveto <- file.path(project_dir, "outputs/analysis/overview/alignment_statistics.reads_distributions.by_batchs.pdf")
ggsave(reads_dist_plot_by_batch_saveto, plot = reads_dist_plot_by_batch, width = 12, height = 7)


# Exon/intron/intergenic reads proportions, 2019 Peng et al. Nature
reads_dist_tab <- file.path(project_dir, 'outputs/analysis/preprocessing/2019_peng_etal_nature/mapping_reports/rnaseq_qc_results.all.txt') %>%
  fread() %>%
  dplyr::mutate(Percentage = stringr::str_remove(Percentage, "%$") %>% as.numeric()) %>%
  tidyr::pivot_wider(names_from = Region, values_from = c(Percentage, Counts))

reads_dist_plot <- ggplot() +
  geom_point(aes(x = Percentage_exonic, y = Percentage_intronic, size = Percentage_intergenic, color = Percentage_overlapping_exon), data = reads_dist_tab, alpha = 0.75) +
  scale_color_gradient(low = rbgc[5], high = rbgc[1]) +
  labs(x = "Reads on exons (%)", y = "Reads on introns (%)", color = "Reads on SJ (%)", size = "Other reads (%)") +
  theme_classic() +
  theme(legend.position = "right")
reads_dist_plot_saveto <- file.path(project_dir, "outputs/analysis/overview/2019_peng_etal_nature.alignment_statistics.reads_distributions.pdf")
ggsave(reads_dist_plot_saveto, plot = reads_dist_plot, width = 8.5, height = 7)


# Gene expression
selected_gene_type <- c(
  "protein_coding",
  "IG_V_gene", "IG_C_gene", "IG_J_gene", "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_D_gene", "IG_D_gene",
  "lncRNA", "miRNA", "snRNA", "TEC", "snoRNA", "scaRNA", "scRNA", "sRNA"
)
dds_obj_path <- file.path(project_dir, "outputs/analysis/gene_expression/deseq2.dds_object.rds")
norm_expr_tab_path <- file.path(project_dir, "outputs/analysis/gene_expression/deseq2.normalized_expression_table.csv")
if (all(file.exists(dds_obj_path, norm_expr_tab_path))) {
  norm_mat <- fread(norm_expr_tab_path) %>% as.data.frame() %>% tibble::rownames_to_column("Geneid") %>% as.matrix()
  dds <- readRDS(dds_obj)
} else {
  count_tab_paths <- file.path(project_dir, 'outputs/analysis/preprocessing/geo_seq/quantification', all_batches, "featureCounts/count_tables") %>%
    lapply(function(x) list.files(x, pattern = ".all_transcripts.read_counts.txt$", full.names = TRUE, recursive = TRUE))

  readcount_table <- count_tab_paths %>%
    lapply(function(x) {
      path_vec <- stringr::str_split(x, "/", simplify = TRUE)
      batch_id <- path_vec[length(path_vec) - 3]
      fread(x) %>% dplyr::rename_with(~paste0(basename(.x), ".", batch_id) %>% stringr::str_remove(".merged.rmdup.bam"), .cols = dplyr::starts_with("/home"))
    }) %>%
    Reduce(dplyr::left_join, .)

  readcount_mat <- readcount_table %>%
    dplyr::select(-c(Chr, Start, End, Strand, Length, gene_id, gene_name)) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Geneid") %>%
    as.matrix()

  non_samples <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", "gene_id", "gene_name")
  selected_sample_info <- data.frame(Sample = colnames(readcount_table)) %>% dplyr::filter(!Sample %in% non_samples) %>%
    dplyr::mutate(Treatment = dplyr::if_else(stringr::str_detect(Sample, "^NC_"), "Control", "With_4sU")) %>%
    dplyr::mutate(Regions_l1 = stringr::str_extract(Sample, "^[0-9]{1,2}([A-Z]{1,2})_", group = 1) %>% purrr::map_chr(~dplyr::if_else(is.na(.x), "NC", .x))) %>%
    dplyr::mutate(Regions_l0 = regions_l0[Regions_l1]) %>%
    tidyr::separate(Sample, into = c("Sample", "Batch"), sep = "\\.")

  # Create DESeq object
  dds <- DESeqDataSetFromMatrix(countData = readcount_mat, colData = selected_sample_info, design = ~1)

  # Filtering
  smallestGroupSize <- 3
  keep <- (rowSums(counts(dds) >= 10) >= smallestGroupSize) & (rowMax(counts(dds)) <= 10000)
  dds <- dds[keep,]

  # Basic processing
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  saveRDS(dds, dds_obj_path)


  # Check FPKM
  new_rowranges <- readcount_table %>%
    dplyr::select(Geneid, Chr, Start, End, Strand) %>%
    dplyr::filter(Geneid %in% names(rowRanges(dds))) %>%
    GenomicRanges::makeGRangesListFromDataFrame(split.field = "Geneid")
  rowRanges(dds) <- new_rowranges[names(rowRanges(dds))]

  fpkm_tab <- fpkm(dds) %>% as.data.frame() %>%
    tibble::rownames_to_column("Geneid") %>%
    tidyr::pivot_longer(-Geneid, values_to = "FPKM", names_to = "Sample_x_batch") %>%
    tidyr::pivot_wider(names_from = "Geneid", values_from = "FPKM") %>%
    tidyr::separate(Sample_x_batch, into = c("Sample", "Batch"), sep = "\\.", remove = FALSE) %>%
    dplyr::mutate(Treatment = dplyr::if_else(stringr::str_detect(Sample, "^NC_"), "Control", "With_4sU")) %>%
    dplyr::mutate(Regions_l1 = stringr::str_extract(Sample, "^[0-9]{1,2}([A-Z]{1,2})_", group = 1) %>% purrr::map_chr(~dplyr::if_else(is.na(.x), "NC", .x))) %>%
    dplyr::mutate(Regions_l0 = regions_l0[Regions_l1]) %>%
    dplyr::select(Sample_x_batch, Sample, Batch, Treatment, Regions_l0, Regions_l1, dplyr::starts_with("ENSMUST"))

  fpkm_gradients <- fpkm_tab %>%
    dplyr::rowwise() %>%
    dplyr::mutate(FPKM_0 = sum(dplyr::c_across(dplyr::starts_with("ENSMUST")) > 0)) %>%
    dplyr::mutate(FPKM_0.1 = sum(dplyr::c_across(dplyr::starts_with("ENSMUST")) > 0.1)) %>%
    dplyr::mutate(FPKM_1 = sum(dplyr::c_across(dplyr::starts_with("ENSMUST")) > 1)) %>%
    dplyr::mutate(FPKM_5 = sum(dplyr::c_across(dplyr::starts_with("ENSMUST")) > 5)) %>%
    dplyr::mutate(FPKM_10 = sum(dplyr::c_across(dplyr::starts_with("ENSMUST")) > 10)) %>%
    dplyr::mutate(FPKM_20 = sum(dplyr::c_across(dplyr::starts_with("ENSMUST")) > 20)) %>%
    dplyr::mutate(FPKM_40 = sum(dplyr::c_across(dplyr::starts_with("ENSMUST")) > 40)) %>%
    dplyr::mutate(FPKM_80 = sum(dplyr::c_across(dplyr::starts_with("ENSMUST")) > 80)) %>%
    dplyr::mutate(FPKM_160 = sum(dplyr::c_across(dplyr::starts_with("ENSMUST")) > 160)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::starts_with("ENSMUST"))

  fpkm_gradients_plot_tab <- fpkm_gradients %>%
    tidyr::pivot_longer(dplyr::starts_with("FPKM"), names_to = "FPKM_gradient", values_to = "Counts") %>%
    dplyr::mutate(FPKM_gradient = factor(FPKM_gradient, levels = c("FPKM_0", "FPKM_0.01", "FPKM_0.1", "FPKM_1", "FPKM_5", "FPKM_10", "FPKM_20", "FPKM_40", "FPKM_80", "FPKM_160")))

  # fpkm_gradients_plot_tab %>% fwrite(file.path(project_dir, "outputs/analysis/gene_expression/deseq2.fpkm_gradients.tsv"))

  fpkm_gradients_plot <- ggplot() +
    geom_line(aes(x = FPKM_gradient, y = Counts, group = Sample_x_batch, color = Batch), alpha = 0.75, data = fpkm_gradients_plot_tab) +
    geom_point(aes(x = FPKM_gradient, y = Counts, group = Sample_x_batch, color = Batch), alpha = 0.75, data = fpkm_gradients_plot_tab) +
    labs(x = "FPKM gradients", y = "Log10(# of genes)", color = "Batches") +
    scale_y_log10() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  fpkm_gradients_plot_saveto <- file.path(project_dir, "outputs/analysis/gene_expression/deseq2.fpkm_gradients.line_plot.pdf")
  ggsave(fpkm_gradients_plot_saveto, plot = fpkm_gradients_plot, width = 8, height = 5)

  fpkm_gradients_plot <- ggplot() +
    geom_line(aes(x = FPKM_gradient, y = Counts, group = Sample_x_batch, color = Batch), alpha = 0.75, data = fpkm_gradients_plot_tab) +
    geom_point(aes(x = FPKM_gradient, y = Counts, group = Sample_x_batch, color = Batch), alpha = 0.75, data = fpkm_gradients_plot_tab) +
    labs(x = "FPKM gradients", y = "Log10(# of genes)", color = "Batches") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  fpkm_gradients_plot_saveto <- file.path(project_dir, "outputs/analysis/gene_expression/deseq2.fpkm_gradients.no_log10.line_plot.pdf")
  ggsave(fpkm_gradients_plot_saveto, plot = fpkm_gradients_plot, width = 8, height = 5)

  # Number of genes per biotype
  transcript_biotype_tab <- file.path(project_dir, "outputs/references/genomic_features/mus_musculus.90.wto_rRNA.transcripts_biotypes.tsv") %>%
    fread(header = FALSE, col.names = c("Transcript_id", "Biotype"))
  other_biotypes <- c(
    "polymorphic_pseudogene", "processed_pseudogene", "processed_transcript", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene",
    "transcribed_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "sense_overlapping", "ribozyme")

  expressed_genes_per_biotype <- fpkm_tab %>%
    tidyr::pivot_longer(dplyr::starts_with("ENSMUST"), values_to = "FPKM", names_to = "Transcript_id") %>%
    dplyr::left_join(transcript_biotype_tab, by = "Transcript_id") %>%
    dplyr::mutate(Biotype = dplyr::case_when(Biotype %in% other_biotypes ~ "Others", TRUE ~ Biotype)) %>%
    dplyr::group_by(Batch, Sample, Sample_x_batch, Treatment, Regions_l0, Regions_l1, Biotype) %>%
    dplyr::summarize(
      FPKM_0 = sum(FPKM > 0), FPKM_0.01 = sum(FPKM > 0.01), FPKM_0.1 = sum(FPKM > 0.1), FPKM_1 = sum(FPKM > 1),
      FPKM_5 = sum(FPKM > 5), FPKM_10 = sum(FPKM > 10), FPKM_20 = sum(FPKM > 20), FPKM_40 = sum(FPKM > 40),
      FPKM_80 = sum(FPKM > 80), FPKM_160 = sum(FPKM > 160)
    ) %>%
    dplyr::ungroup()

  egpb_plot_tab <- expressed_genes_per_biotype %>% dplyr::select(Batch, Sample, FPKM_0, Biotype) %>%
    tidyr::pivot_wider(names_from = Biotype, values_from = FPKM_0) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(All = sum(dplyr::c_across(-c(Batch, Sample)))) %>%
    dplyr::mutate(dplyr::across(-c(Batch, Sample, All), ~.x / All * 100)) %>%
    dplyr::select(-All) %>%
    tidyr::pivot_longer(-c(Batch, Sample), names_to = "Biotype", values_to = "Percentage") %>%
    dplyr::mutate(Biotype = forcats::fct_reorder(Biotype, Percentage)) %>%
    dplyr::filter(Biotype != "protein_coding", Biotype != "Others")

  egpb_plot <- ggplot() +
    geom_col(aes(x = Sample, y = Percentage, fill = Biotype), data = egpb_plot_tab, position = "stack") +
    facet_wrap(~Batch, ncol = 2, scales = "free_x") +
    scale_fill_npg() +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  egpb_plot_saveto <- file.path(project_dir, "outputs/analysis/gene_expression/deseq2.expressed_genes_per_biotype.pdf")
  ggsave(egpb_plot_saveto, plot = egpb_plot, width = 10)



  # Normalization by varianceStablizingTransformation and removing batch effects by limma::removeBatchEffect
  # Note the vsd was not updated, i.e., with batch effects
  vsd <- vst(dds, blind=FALSE)
  model_matrix <- model.matrix(~1, colData(vsd))
  norm_mat <- assay(vsd) %>% limma::removeBatchEffect(batch=vsd$Batch, design=model_matrix)
  norm_mat %>% as.data.frame() %>% tibble::rownames_to_column("Geneid") %>% fwrite(norm_expr_tab_path)

  # Check the variance stablized 
  p_tab <- norm_mat %>% as.data.frame() %>% tibble::rownames_to_column("Geneid") %>%
    tidyr::pivot_longer(-Geneid, names_to = "Sample", values_to = "Expression_VST") %>%
    dplyr::mutate(Sample_x_batch = Sample) %>%
    tidyr::separate(Sample, into = c("Sample", "Batch"), sep = "\\.") %>%
    dplyr::group_by(Sample_x_batch) %>%
    dplyr::mutate(Mean_exp_vst = median(Expression_VST)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Sample_x_batch = forcats::fct_reorder(Sample_x_batch, Mean_exp_vst))

  p <- ggplot(p_tab) +
    geom_violin(aes(y = Expression_VST, x = Sample_x_batch, fill = Batch)) +
    geom_boxplot(aes(y = Expression_VST, x = Sample_x_batch, fill = Batch), outlier.shape = NA, width = 0.5) +
    scale_fill_npg() +
    labs(x = "Samples", y = "Normalized expression per gene") +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  file.path(project_dir, "outputs/analysis/gene_expression/deseq2.normalized_expression.pdf") %>%
    ggsave(plot = p, width = 20, height = 6.5)

  # Check the distribution at low dimentions
  pca <- norm_mat %>% t() %>% prcomp()
  pca_tab <- pca$x %>% as.data.frame() %>% tibble::rownames_to_column("Sample_x_batch") %>%
    tidyr::separate(Sample_x_batch, into = c("Sample", "Batch"), remove = FALSE, sep = "\\.") %>%
    dplyr::mutate(Treatment = dplyr::if_else(stringr::str_detect(Sample, "^NC_"), "Control", "With_4sU")) %>%
    dplyr::mutate(Regions_l1 = stringr::str_extract(Sample, "^[0-9]{1,2}([A-Z]{1,2})_", group = 1) %>% purrr::map_chr(~dplyr::if_else(is.na(.x), "NC", .x))) %>%
    dplyr::mutate(Regions_l0 = regions_l0[Regions_l1]) %>%
    dplyr::mutate(Regions_l0 = factor(Regions_l0, levels = c("Endoderm", "Mesoderm", "Ectoderm", "Control")))

  pca_plot <- ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = Regions_l0), data = pca_tab) +
    theme_classic()
  pca_plot_saveto <- file.path(project_dir, "outputs/analysis/gene_expression/deseq2.gene_expression.color_by_regions_l0.pca_plot.pdf")
  ggsave(pca_plot_saveto, plot = pca_plot, width = 5.5, height = 5)

  pca_plot <- ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = Regions_l1), data = pca_tab) +
    theme_classic()
  pca_plot_saveto <- file.path(project_dir, "outputs/analysis/gene_expression/deseq2.gene_expression.color_by_regions_l1.pca_plot.pdf")
  ggsave(pca_plot_saveto, plot = pca_plot, width = 5.5, height = 5)

  pca_plot <- ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = Batch), data = pca_tab) +
    theme_classic()
  pca_plot_saveto <- file.path(project_dir, "outputs/analysis/gene_expression/deseq2.gene_expression.color_by_batch.pca_plot.pdf")
  ggsave(pca_plot_saveto, plot = pca_plot, width = 5.5, height = 5)
}

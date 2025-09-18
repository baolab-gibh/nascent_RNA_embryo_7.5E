#!/usr/bin/env Rscript
# File: slam_seq.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Nov 09, 2024
# Updated:

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggsci)
})

project_dir <- "~/Documents/projects/wp_vasaseq"


gene_biotype <- list(
  Ribozyme = c("ribozyme"), 
  rRNA = c("rRNA", "rRNA_pseudogene"),
  Protein_coding = c("protein_coding"),
  # Mt_RNA = c("Mt_rRNA", "Mt_tRNA"),
  Non_coding_rna = c("lncRNA", "miRNA", "misc_RNA", "sRNA", "scRNA", "scaRNA", "snRNA", "snoRNA"),
  IG_or_TR= c("IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_J_gene", "IG_J_pseudogene", "IG_V_gene", "IG_V_pseudogene", "IG_pseudogene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_V_gene", "TR_V_pseudogene"),
  Uknnown_or_pseudogene = c("processed_pseudogene", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "vault_RNA", "TEC", "artifact")
)



count_tab <- file.path(project_dir, "outputs/analysis/preprocessing/slam_seq/read_count/20241031_slamseq") %>%
  list.files(recursive = TRUE, full.names = TRUE, pattern = ".txt$") %>%
  lapply(function(x) {
    path_vec <- stringr::str_split(x, "/", simplify = TRUE)
    batch_id <- path_vec[length(path_vec) - 1]
    fread(x, sep = "\t") %>% dplyr::mutate(Batch = batch_id) %>% dplyr::rename_with(function(x) return("Readcounts"), .cols = dplyr::starts_with("/home"))
  }) %>%
  Reduce(rbind, .) %>%
  dplyr::mutate(gene_type = dplyr::case_when(
    gene_type %in% gene_biotype[["Ribozyme"]] ~ "Ribozyme",
    gene_type %in% gene_biotype[["rRNA"]] ~ "rRNA",
    gene_type %in% gene_biotype[["Protein_coding"]] ~ "Protein_coding",
    gene_type %in% gene_biotype[["Non_coding_rna"]] ~ "Non_coding_rna",
    gene_type %in% gene_biotype[["IG_or_TR"]] ~ "IG_or_TR",
    gene_type %in% gene_biotype[["Uknnown_or_pseudogene"]] ~ "Uknnown_or_pseudogene",
    T ~ gene_type,
  )) %>%
  dplyr::group_by(Batch, gene_type) %>%
  dplyr::summarise(Readcounts = sum(Readcounts)) %>%
  dplyr::group_by(Batch) %>%
  dplyr::mutate(Proportion = Readcounts / sum(Readcounts), Percentage = Proportion * 100) %>%
  dplyr::mutate(gene_type = factor(gene_type, levels = c("Mt_rRNA", "Mt_tRNA", "Ribozyme", "rRNA", "Non_coding_rna", "IG_or_TR", "Protein_coding", "Uknnown_or_pseudogene")))


p <- ggplot(count_tab) +
  geom_bar(aes(x = Batch, y = Percentage, fill = gene_type), stat = "identity", position = "stack") +
  labs(title = "Gene species percentage by read counts") +
  scale_fill_npg() +
  theme_classic()
psaveto <- file.path(project_dir, "outputs/analysis/overview/slam_seq/plots/slam_seq.read_count.pdf")
ggsave(p, file = psaveto, width = 8, height = 6)

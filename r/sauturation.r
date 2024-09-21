#!/usr/bin/env Rscript
# File: sauturation.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Sep 12, 2024
# Updated:

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(DESeq2)
  library(GenomicRanges)
})

project_dir <- "~/Documents/projects/wp_vasaseq"

all_batches <- c("240409_Lib_embryo", "240612_Lib_28region", "240620_Lib_38region", "240703_Lib_32region", "240710_Lib_37region", "240717_Lib_28region")


sauturation_fpkm_gradient_tab_path <- file.path(project_dir, "outputs/analysis/preprocessing/geo_seq/sauturation/sauturation_fpkm_gradients.tsv")
if (!file.exists(sauturation_fpkm_gradient_tab)) {
  sauturation_fpkm_gradient_tab <- file.path(project_dir, "outputs/analysis/preprocessing/geo_seq/sauturation") %>%
    list.files(pattern = "read_counts.*.txt$", full.names = TRUE, recursive = TRUE) %>%
    split(sort(rep(all_batches, 10))) %>%
    lapply(function(path_vec) {
      path_vec %>% lapply(function(p) {
        batch_vec <- p %>% stringr::str_split("/") %>% unlist() %>% purrr::keep(~stringr::str_detect(.x, "_Lib_"))
        percent <- p %>% stringr::str_extract("read_counts\\.([0-9]+)\\.txt", group = 1) %>% as.numeric()

        readcount_table <- fread(p) %>%
          dplyr::select(-c(Geneid, Length)) %>%
          dplyr::rename_with(~stringr::str_extract(.x, "[0-9A-Z]+_[ATCG]+"), dplyr::starts_with("/")) %>%
          dplyr::mutate(Batch = batch_vec, Percent = percent) %>%
          dplyr::group_by(Batch, Percent, gene_name, gene_id) %>%
          dplyr::summarise(Chr = dplyr::first(Chr), Start = min(Start), End = max(End), Strand = dplyr::first(Strand), across(where(is.integer), ~sum(.x))) %>%
          dplyr::ungroup()

        count_data <- readcount_table %>%
          dplyr::select(-c(Batch, Percent, gene_name, Chr, Start, End, Strand)) %>%
          as.data.frame() %>%
          tibble::column_to_rownames("gene_id") %>%
          as.matrix()

        col_data <- data.frame(Sample = colnames(count_data)) %>%
          dplyr::mutate(Treatment = dplyr::if_else(stringr::str_detect(Sample, "^NC_"), "Control", "With_4sU"))

        dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~1)

        smallestGroupSize <- 3
        keep <- (rowSums(counts(dds) >= 10) >= smallestGroupSize) & (rowMax(counts(dds)) <= 10000)
        dds <- dds[keep,]

        # Basic processing
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds)

        new_rowranges <- readcount_table %>%
          dplyr::select(gene_id, Chr, Start, End, Strand) %>%
          dplyr::filter(gene_id %in% names(rowRanges(dds))) %>%
          GenomicRanges::makeGRangesListFromDataFrame(split.field = "gene_id")
        rowRanges(dds) <- new_rowranges[names(rowRanges(dds))]

        fpkm_tab <- fpkm(dds) %>%
          as.data.frame() %>%
          tibble::rownames_to_column("gene_id") %>%
          tidyr::pivot_longer(-gene_id, values_to = "FPKM", names_to = "Sample") %>%
          tidyr::pivot_wider(names_from = "gene_id", values_from = "FPKM") %>%
          dplyr::mutate(Treatment = dplyr::if_else(stringr::str_detect(Sample, "^NC_"), "Control", "With_4sU")) %>%
          dplyr::mutate(Regions_l1 = stringr::str_extract(Sample, "^[0-9]{1,2}([A-Z]{1,2})_", group = 1) %>% purrr::map_chr(~dplyr::if_else(is.na(.x), "NC", .x))) %>%
          dplyr::mutate(Batch = batch_vec, Percent = percent) %>%
          dplyr::select(Sample, Batch, Percent, Treatment, Regions_l1, dplyr::starts_with("ENSMUSG")) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(FPKM_0 = sum(dplyr::c_across(dplyr::starts_with("ENSMUSG")) > 0)) %>%
          dplyr::mutate(FPKM_0.1 = sum(dplyr::c_across(dplyr::starts_with("ENSMUSG")) > 0.1)) %>%
          dplyr::mutate(FPKM_1 = sum(dplyr::c_across(dplyr::starts_with("ENSMUSG")) > 1)) %>%
          dplyr::mutate(FPKM_5 = sum(dplyr::c_across(dplyr::starts_with("ENSMUSG")) > 5)) %>%
          dplyr::mutate(FPKM_10 = sum(dplyr::c_across(dplyr::starts_with("ENSMUSG")) > 10)) %>%
          dplyr::mutate(FPKM_20 = sum(dplyr::c_across(dplyr::starts_with("ENSMUSG")) > 20)) %>%
          dplyr::mutate(FPKM_40 = sum(dplyr::c_across(dplyr::starts_with("ENSMUSG")) > 40)) %>%
          dplyr::mutate(FPKM_80 = sum(dplyr::c_across(dplyr::starts_with("ENSMUSG")) > 80)) %>%
          dplyr::mutate(FPKM_160 = sum(dplyr::c_across(dplyr::starts_with("ENSMUSG")) > 160)) %>%
          dplyr::ungroup() %>%
          dplyr::select(-dplyr::starts_with("ENSMUSG")) %>%
          tidyr::pivot_longer(dplyr::starts_with("FPKM"), names_to = "FPKM_gradient", values_to = "Counts") %>%
          dplyr::mutate(FPKM_gradient = factor(FPKM_gradient, levels = c("FPKM_0", "FPKM_0.01", "FPKM_0.1", "FPKM_1", "FPKM_5", "FPKM_10", "FPKM_20", "FPKM_40", "FPKM_80", "FPKM_160")))
      }) %>%
      dplyr::bind_rows()
    }) %>%
    dplyr::bind_rows()
  sauturation_fpkm_gradient_tab %>% fwrite(sauturation_fpkm_gradient_tab_path)
} else {
  sauturation_fpkm_gradient_tab <- fread(sauturation_fpkm_gradient_tab_path)
}



p <- sauturation_fpkm_gradient_tab %>% dplyr::filter(FPKM_gradient %in% c("FPKM_1", "FPKM_5", "FPKM_10")) %>%
  ggplot() +
  geom_point(aes(x = Percent, y = Counts, group = Sample, color = Regions_l1), alpha = 0.5) +
  geom_line(aes(x = Percent, y = Counts, group = Sample), alpha = 0.5) +
  facet_grid(Batch ~ FPKM_gradient, scales = "free") +
  theme_bw()

p_saveto <- file.path(project_dir, "outputs/analysis/preprocessing/geo_seq/sauturation/sauturation.fpkm_gradients.pdf")
ggsave(p_saveto, width = 8, height = 7)

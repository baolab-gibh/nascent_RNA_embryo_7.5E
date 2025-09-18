#!/usr/bin/env Rscript
# File: explore_pijuan_etal.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com

library(data.table)
library(tidyverse)
library(Seurat)

# available_dev_stage <- c(
#   "embryonic day 6.5", "embryonic day 6.75", "embryonic day 7.0", "embryonic day 7.25", "embryonic day 7.5", "embryonic day 7.75",
#   "embryonic day 8.0", "embryonic day 8.25", "embryonic day 8.5", "mixed"
# )

selected_dev_stage <- c("embryonic day 7.0", "embryonic day 7.25", "embryonic day 7.5", "embryonic day 7.75", "embryonic day 8.0")

in_file <- "/home/zzhang/Documents/projects/wp_vasaseq/inputs/PijuanSala_Nature_2019/E-MTAB-6967.sdrf.txt"

# Check avialable columns
fread(in_file, nrows = 2) %>% colnames()

# Overvieiw sampls
tab <- fread(in_file) %>%
  dplyr::select(
    Source_Name = `Source Name`, developmental_stage = `Characteristics[developmental stage]`, Assay_Name = `Assay Name`,
    number_of_embryos_in_pool = `Parameter Value[number of embryos in pool]`, technical_replicate_group = `Comment[technical replicate group]`
  )
tab %>% dplyr::select(Source_Name, number_of_embryos_in_pool) %>% dplyr::distinct() %>% dplyr::summarise(total_samples = n(), total_embryos = sum(number_of_embryos_in_pool))

tab %>% dplyr::filter(!developmental_stage %in% c("embryonic day 8.25"), !Source_Name %in% c("Sample 11")) %>%
  dplyr::select(Source_Name, developmental_stage, number_of_embryos_in_pool) %>%
  dplyr::group_by(Source_Name) %>%
  dplyr::summarise(developmental_stage = unique(developmental_stage), n=n(), number_of_embryos_in_pool = unique(number_of_embryos_in_pool)) %>%
  dplyr::summarise(total_samples = n(), total_embryos = sum(number_of_embryos_in_pool))

tab %>% dplyr::filter(!developmental_stage %in% c("embryonic day 8.25"), !Source_Name %in% c("Sample 11")) %>% dplyr::pull(Source_Name) %>% unique()

tab %>% dplyr::filter(developmental_stage %in% selected_dev_stage) %>% tidyr::separate(Assay_Name, remove = FALSE, into = c("x", "y", "index", "z", "a"), sep = "_")
tab %>% dplyr::select(developmental_stage) %>% unique()

tab %>% dplyr::filter(developmental_stage %in% "embryonic day 7.0") %>% dplyr::select(Source_Name, developmental_stage) %>% table()



tab %>% dplyr::select(Source_Name, developmental_stage) %>% table
tab %>% dplyr::select(Source_Name) %>% table()
# tab %>% dplyr::filter(`Characteristics[developmental stage]` != "embryonic day 8.25")

tab %>% dplyr::select(`Source Name`) %>% table()

embryo <- Read10X("/home/zzhang/Documents/projects/wp_vasaseq/inputs/PijuanSala_Nature_2019/10X/outs/raw_feature_bc_matrix") %>%
  CreateSeuratObject(project = "embryo", min.cells = 3, min.features = 200)

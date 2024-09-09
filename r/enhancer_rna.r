#!/usr/bin/env Rscript
# File: enhancer_rna.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 06, 2024
# Updated:

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)

  library(DESeq2)
})

rbgc <- brewer.pal(n = 11, name = "RdBu")

all_batches <- c("240409_Lib_embryo", "240612_Lib_28region", "240620_Lib_38region", "240703_Lib_32region", "240710_Lib_37region", "240717_Lib_28region")
regions_l0 <- c("Ectoderm", "Ectoderm", "Endoderm", "Endoderm", "Ectoderm", "Ectoderm", "Mesoderm", "Mesoderm", "Control") %>%
  purrr::set_names(c('A', 'P', 'EA', 'EP', 'L', 'R', 'MA', 'MP', "NC"))

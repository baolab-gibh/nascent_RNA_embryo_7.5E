#!/usr/bin/env Rscript
# File: corn_coordination.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jul 22, 2024
# Updated: Jul 22, 2024

options(bitmapType = "cairo")
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggsci)
})


project_dir <- "~/Documents/projects/wp_vasaseq"

# Sigmod function
sigmod <- function(x) 1.0 / (1.0 + exp(3) ^ (-sqrt(x)))
shift_val <- sigmod(1:2) %>% (function(vc) c(min(vc), max(vc))) %>% quantile(probs = 0.075) %>% `*`(0.9975)
shrink <- function(y_pos, shift_val) {
  sigmod(1:y_pos) %>% (function(vc) c(min(vc), max(vc))) %>% quantile(probs = c(0.025, 0.065, 0.100, 0.135)) %>% as.vector() %>% `-`(shift_val)
}

# Dyanmics by corn plot
corn_axis <- data.frame(y_pos = sort(rep(2:16, 4))) %>%
  dplyr::group_by(y_pos) %>%
  dplyr::reframe(x_pos = shrink(y_pos, shift_val)) %>%
  dplyr::bind_rows(dplyr::filter(., y_pos == 2) %>% dplyr::slice_min(x_pos, n = 2) %>% dplyr::mutate(y_pos = 1)) %>%
  dplyr::bind_rows(dplyr::mutate(., x_pos = -x_pos)) %>%
  dplyr::mutate(Layers = paste0("L", y_pos), Layers = forcats::fct_reorder(Layers, y_pos)) %>%
  dplyr::group_by(Layers) %>%
  dplyr::arrange(x_pos) %>%
  dplyr::mutate(Regions_l1 = {
    group_name <- dplyr::cur_group()
    new_dtfm <- dplyr::cur_data()
    if (group_name == "L1") {
      new_dtfm %>% dplyr::mutate(Regions_l1 = c("A", "L", "R", "P"))
    } else {
      new_dtfm %>% dplyr::mutate(Regions_l1 = c("EA", "MA", "A", "L", "R", "P", "MP", "EP"))
    }
  }) %>%
  tidyr::unnest(Regions_l1, names_sep = ".") %>%
  dplyr::select(y_pos, x_pos, Layers, Regions_l1 = Regions_l1.Regions_l1) %>%
  dplyr::mutate(Regions_l1 = factor(Regions_l1, levels = c("EA", "MA", "A", "L", "R", "P", "MP", "EP"))) %>%
  dplyr::mutate(Regions_l0 = dplyr::case_when(
    Regions_l1 %in% c("A", "L", "R", "P") ~ "Ectoderm",
    Regions_l1 %in% c("EA", "EP") ~ "Endoderm",
    Regions_l1 %in% c("MA", "MP") ~ "Mesoderm",
  ))

corn_axis %>% fwrite(file.path(project_dir, "inputs/reference/corn_coordination/corn_axis.e_7_5.csv"))

p <- ggplot() +
  geom_point(aes(x = x_pos, y = Layers, fill = Regions_l1), corn_axis %>% dplyr::filter(!Regions_l1 %in% c("MA", "MP")), shape = 21, color = "black", size = 8) +
  geom_point(aes(x = x_pos, y = Layers, fill = Regions_l1), corn_axis %>% dplyr::filter(Regions_l1 %in% c("MA", "MP")), shape = 23, color = "black", size = 8) +
  scale_fill_npg() +
  labs(x = NULL, y = "Layers") +
  lims(x = c(-0.0075, 0.0075)) +
  theme_classic() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line())
corn_axis_plot_saveto <- file.path(project_dir, "inputs/reference/corn_coordination/corn_axis.e_7_5.pdf")
ggsave(corn_axis_plot_saveto, plot = p, width = 4, height = 4.5)

#!/usr/bin/env Rscript

# Plot nascent RNA analysis results from Grand-SLAM

suppressPackageStartupMessages({
  # Plotting
  library(ggplot2)
  library(patchwork)

  # Data manipulation
  library(data.table)
  library(dplyr)

  # Utilities
  library(magrittr)
  library(optparse)
})


#' Main
main <- function() {
  parser <- optparse::OptionParser(usage = "%prog [options] <work_dir>")
  parser <- optparse::add_option(
    parser, c("-o", "--save-to"), metavar="FILE", type = "character", default = "output.pdf", help = "Save to. [Default: %default]"
  )
  parser <- optparse::add_option(
    parser, c("-c", "--ctrl-sample"), metavar = "CHAR", type = "character", default = NULL,
    help = "Control sample name(s). Multiple samples should be separated by comma. [Default: %default]"
  )
  parser <- optparse::add_option(
    parser, c("-W", "--fig-width"), metavar="INT", type = "integer", default = 9, help = "Figure width. [Default: %default]"
  )
  parser <- optparse::add_option(
    parser, c("-H", "--fig-height"), metavar="INT", type = "integer",
    help = "Figure width. If it is not provided, it will be adjusted based on the number of samples. [Default: %default]"
  )

  # Positional arguments
  arguments <- optparse::parse_args(parser, positional_arguments = 1, convert_hyphens_to_underscores = TRUE)
  work_dir <- arguments$args

  # Optional arguments
  options <- arguments$options
  save_to <- options$save_to
  ctrl_sample <- options$ctrl_sample
  fig_width <- options$fig_width
  fig_height <- options$fig_height

  ctrl_sample <- ctrl_sample %>% strsplit(",") %>% unlist()

  # Plot NTR for each sample
  ntr_overview_tab <- file.path(work_dir, "all_samples") %>%
    list.files("*.ntrstat.tsv$", recursive = TRUE, full.names = TRUE) %>%
    data.table::fread()

  ntr_case_tab <- ntr_overview_tab %>% dplyr::filter(! Condition %in% ctrl_sample)
  ntr_ctrl_tab <- ntr_overview_tab %>% dplyr::filter(Condition %in% ctrl_sample)
  p_ntr_overview <- ggplot2::ggplot(ntr_case_tab) +
    ggplot2::geom_line(
      aes(x = Type, y = ntr, group = Condition), color = "grey", linetype = "dotted",
      position = ggplot2::position_dodge(0.1)
    ) +
    ggplot2::geom_pointrange(
      aes(x = Type, y = ntr, ymin = ntr_lower, ymax = ntr_upper, color = Condition), size = 0.25,
      position = ggplot2::position_dodge(0.1)
    ) +
    ggplot2::geom_pointrange(
      aes(x = Type, y = ntr, ymin = ntr_lower, ymax = ntr_upper), data = ntr_ctrl_tab, color = "black", size = 0.25,
      position = ggplot2::position_dodge(0.1)
    ) +
    ggplot2::geom_line(
      aes(x = Type, y = ntr, group = Condition), data = ntr_ctrl_tab, color = "grey", linetype = "dotted",
      position = ggplot2::position_dodge(0.1)
    ) +
    #ggsci::scale_color_npg() +
    ggplot2::labs(x = "Read source", y = "New to total ratio (NTR)", color = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "top")

  # Mismatch along position for each sample
  per_sample_files <- file.path(work_dir, "per_sample") %>%
    list.files("*.mismatchdetails.tsv$", recursive = TRUE, full.names = TRUE)
  mm_position_tab <- per_sample_files %>%
    lapply(function(p) {
      sample_id <- p %>% stringr::str_split("/", simplify = TRUE) %>% `[`(length(.) - 1)
      data.table::fread(p) %>% dplyr::mutate(Sample_id = sample_id)
    }) %>%
    Reduce(rbind, .) %>%
    dplyr::mutate(Mismatch_rate = round(Mismatches / Coverage, 6)) %>%
    dplyr::mutate(Mismatch_type = paste0(Genomic, "-", Read)) %>%
    dplyr::mutate(Genomic = factor(Genomic, levels = c("T", "A", "C", "G")))

  # Overview of all mismatches
  mm_overview_tab <- mm_position_tab %>%
    dplyr::filter(Category %in% c("Exonic", "Intronic")) %>%
    dplyr::group_by(Sample_id, Genomic, Read) %>%
    dplyr::summarize(Coverage = sum(Coverage), Mismatches = sum(Mismatches)) %>%
    dplyr::mutate(Mismatch_rate = Mismatches / Coverage) %>%
    dplyr::mutate(Mismatch_type = paste0(Genomic, "->", Read))

  p_mm_overview <- ggplot2::ggplot(mm_overview_tab) +
    ggplot2::geom_line(aes(x = Mismatch_type, y = Mismatch_rate, group = Sample_id), linetype = "dotted", color = "grey") +
    ggplot2::geom_point(aes(x = Mismatch_type, y = Mismatch_rate, color = Sample_id, group = Sample_id)) +
    # ggsci::scale_color_npg() +
    ggplot2::labs(x = "Conversion", y = "Conversion rate", color = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "top")

  # T->X rate per position from exonic reads
  exonic_mm_tab <- mm_position_tab %>% dplyr::filter(Genomic == "T", Category %in% c("Exonic"))
  p_mm_exonic_pos <- ggplot2::ggplot(exonic_mm_tab) +
    ggplot2::geom_line(aes(x = Position, y = Mismatch_rate, color = Read)) +
    ggh4x::facet_nested_wrap(Sample_id ~ ., ncol = 1, scales = "free_y") +
    ggplot2::labs(x = "Position", y = "T->C rate (Exonic)", color = "T->A/C/G") +
    ggplot2::theme_classic()

  # T->X rate per position from intronic reads
  intronic_tab <- mm_position_tab %>% dplyr::filter(Genomic == "T", Category %in% c("Intronic"))
  p_mm_intronic_pos <- ggplot2::ggplot(intronic_tab) +
    ggplot2::geom_line(aes(x = Position, y = Mismatch_rate, color = Read)) +
    ggh4x::facet_nested_wrap(. ~ Sample_id, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = "Position", y = "T->C rate (Intronic)", color = "T->A/C/G") +
    ggplot2::theme_classic()

  # Combine all plots
  n_treat_samples <- length(per_sample_files) - length(ctrl_sample)
  mismatch_panel_height <- ceiling(max(1, n_treat_samples / 3))
  p <- ((p_ntr_overview | p_mm_overview) + patchwork::plot_layout(ncol = 2, widths = c(2, 3))) /
    ((p_mm_exonic_pos | p_mm_intronic_pos) + patchwork::plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "top")) +
    patchwork::plot_layout(ncol = 1, heights = c(1, mismatch_panel_height))

  if (is.null(fig_height)) {
    n_samples <- ntr_overview_tab %>% dplyr::pull(Condition) %>% unique() %>% length()
    fig_height <- n_samples + 3.5
  }

  ggplot2::ggsave(p, file = save_to, width = fig_width, height = fig_height)
}


main()

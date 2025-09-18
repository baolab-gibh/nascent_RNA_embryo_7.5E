#!/usr/bin/env python3
# File: split_bam_by_bin.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Apr 14, 2025
# Updated: 

import click
import pysam
import polars as pls


def load_barcode_info(raw_bc_file, binned_bc_file, min_umi_count=500):
    '''Load barcode information.'''
    raw_bc_tab = pls.read_csv(raw_bc_file, separator="\t").rename({"Barcode": "Raw_barcode"})
    binned_bc_tab = pls.read_csv(binned_bc_file, separator="\t").select("Barcode").rename({"Barcode": "Binned_barcode"})

    valid_barcodes = dict(
        pls.concat([raw_bc_tab, binned_bc_tab], how="horizontal")
        .with_columns(pls.col("count").sum().over("Binned_barcode").alias("Binned_bc_count"))
        .filter(pls.col("Binned_bc_count") >= min_umi_count)
        .select(pls.col("Raw_barcode"), pls.col("Binned_barcode"))
        .unique()
        .group_by(pls.col("Binned_barcode"))
        .agg(raw_barcodes=pls.col("Raw_barcode"))
        .iter_rows()
    )

    return valid_barcodes


raw_bc_file = "./ST110272_A1_Raw/ST110272_A1_Raw_count_detail.txt.gz"
binned_bc_file = "./ST110272_A1_Bin20/ST110272_A1_Bin20_count_detail.txt.gz"
valid_bc = load_barcode_info(raw_bc_file, binned_bc_file)


def subset_bam_by_tag(in_bam, tag, value):
    with pysam.AlignmentFile(in_bam, "rb") as in_hand:
        for per_read in in_hand:
            per_read.get_tag(tag)


def split_bam_by_bin(in_bam, out_dir, n_bins):
    pass


def main():
    pass


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# File: dedup_by_umi.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Mar 16, 2025
# Updated:

import pysam
from multiprocessing import Pool

import click


def group_reads(read_list):
    return [], []


def dedup_by_umi(read_list, method: str = "random", random_seed: int = 42):
    return []


@click.command()
@click.argument('in_bam', type=click.Path(exists=True), help="Input BAM file")
@click.argument('out_bam', type=click.Path(), help="Output BAM file")
@click.option("-@", "--n-cpus", type=int, default=1, help="Number of CPUs to use for deduplication")
@click.option("-s", "--chunk-size", type=int, default=10000, show_default=True, help="Chunk size per split.")
@click.option("-i", "--insert-size", type=int, default=150, show_default=True, help="Insert size.")
@click.option("-g", "--group-tag", metavar="STR", multiple=True, show_default=True, default=None, help="Group reads by TAG. E.g., CB")
@click.option("-m", "--dedup-method", default="merge", show_default=True, type=click.Choice(["merge", "longest", "left", "right", "random"]), help="Deduplication method.")
def main(in_bam, out_bam, n_cpus, chunk_size: int = 10000):
    with (
        pysam.AlignmentFile(in_bam, "rb", threads=n_cpus) as inhandle,
        pysam.AlignmentFile(out_bam, "wb", threads=n_cpus) as outhandle
    ):
        read_list = []
        for idx, read in enumerate(inhandle):
            read_list.append(read)
            if idx % chunk_size == 0:
                dup_reads, nondup_reads = group_reads(read)
                dedup_reads = dedup_by_umi(dup_reads)

                for r in dedup_reads + nondup_reads:
                    outhandle.write(r)
                read_list = []

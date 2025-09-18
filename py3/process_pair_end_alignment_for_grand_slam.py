#!/usr/bin/env python3
# File: process_pair_end_alignment_for_grand_slam.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jul 30, 2025

import click
import pysam


CIGAR_OPS = "MIDNSHP=X"
TRANSLATOR = str.maketrans("ATCGNatcgn^", "TAGCNtagcn^")


@click.command()
@click.argument("in_bam_file", type=click.Path(exists=True, readable=True), metavar="<in.bam>")
@click.argument("out_bam_file", type=click.Path(exists=False, writable=True), metavar="<out.bam>")
@click.option("-@", "--threads", type=int, default=1, show_default=True, help="number of CPUs for compression or decompression.")
def main(in_bam_file, out_bam_file, **kwargs):
    threads = kwargs["threads"]
    d_threads, c_threads = max(1, int(threads / 2)), max(1, int(threads / 2))
    with pysam.AlignmentFile(in_bam_file, "r", threads=d_threads) as in_bam, \
         pysam.AlignmentFile(out_bam_file, "wb", template=in_bam, header=in_bam.header, threads=c_threads) as out_bam:
        for per_read in in_bam.fetch():
            per_read.flag |= (2047 - 64)
            per_read.flag &= (2047 - 128)
            if per_read.flag & 16 == 1:
                if per_read.query_name is not None:
                    per_read.query_name += "-1"

                if per_read.query_sequence is not None:
                    per_read.query_sequence = per_read.query_sequence[::-1].translate(TRANSLATOR)

                if per_read.query_qualities is not None:
                    per_read.query_qualities = per_read.query_qualities[::-1]

                if per_read.flag is not None:
                    per_read.flag ^= 16
            else:
                if per_read.query_name is not None:
                    per_read.query_name += "-2"

            out_bam.write(per_read)


if __name__ == "__main__":
    main(max_content_width=100)

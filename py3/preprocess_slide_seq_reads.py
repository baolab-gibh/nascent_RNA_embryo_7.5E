#!/usr/bin/env python3
# File: preprocess_slide_seq_reads.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jul 17, 2024
# Updated: Jul 17, 2024

# The library structure. Read 1: BC1(4bp)+CGACTCACTACAGGG+BC2(4bp)TCGGTGACACGATCG+BC3(4bp)+10umi TTT*

import sys
import gzip
import pathlib

from argparse import ArgumentParser

parser = ArgumentParser(description="Preprocess slide-seq reads")
parser.add_argument("-1", "--read1", type=str, required=True, help="input fastq file")
parser.add_argument("-2", "--read2", type=str, required=True, help="input fastq file")
parser.add_argument("-o", "--output", type=str, required=True, help="Output fastq file")
parser.add_argument("-w", "--white-list", type=str, required=True, help="White list file")

args = parser.parse_args()

with (
    gzip.open(args.read1, "rt") as f1,
    gzip.open(args.read2, "rt") as f2,
    gzip.open(args.output, "wb") as fout,
    open(args.white_list, "r") as fwl,
):
    read_id = None
    read_seq_1, read_seq_2 = None, None
    read_qual_1, read_qual_2 = None, None
    bc_wl = list(set([x.strip() for x in fwl]))

    for idx, (line1, line2) in enumerate(zip(f1, f2)):
        if idx % 4 == 0:
            read_id = str(line2).strip().split(" ").pop(0)
        elif idx % 4 == 1:
            read_seq_1 = str(line1).strip()
            read_seq_2 = str(line2).strip()
        elif idx % 4 == 2:
            continue
        elif idx % 4 == 3:
            read_qual_2 = str(line2).strip()
            if read_seq_1 and read_qual_2 and read_id and read_qual_2:
                if len(read_seq_1) <= 52: continue

                bc1, bc2, bc3 = str(read_seq_1[:4]), str(read_seq_1[19:23]), str(read_seq_1[38:42])
                if not all([bc1 in bc_wl, bc2 in bc_wl, bc3 in bc_wl]): continue

                umi = str(read_seq_1[42:52])
                barcode = bc1 + bc2 + bc3
                fout.write(f"{read_id}:CB_{barcode}:UMI_{umi}\n{read_seq_2}\n+\n{read_qual_2}\n".encode())

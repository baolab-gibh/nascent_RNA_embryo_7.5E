#!/usr/bin/env python3
# File: preprocess_slide_seq_reads.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jul 17, 2024
# Updated: Jul 17, 2024

# The library structure. Read 1: BC1(4bp)+CGACTCACTACAGGG+BC2(4bp)TCGGTGACACGATCG+BC3(4bp)+10umi TTT*

import re
import pgzip
from pathlib import Path
from argparse import ArgumentParser
from multiprocessing import Pool

import tqdm

BARCODE_PATTERN = re.compile("([A-Z])(\\d+)")
TRANS_DICT = str.maketrans("ACGTN", "TGCAN")

def decode_sequence(ttl_seq, pattern):
    """Fetch the barcode and umi information from the read 1."""
    barcode, umiseq, linker_seq_list, poly_t_list = "", "", [], []
    for ele, count in pattern:
        count = int(count)
        ele_seq = ttl_seq[:count]
        if ele == "C":
            barcode += ele_seq
        elif ele == "L":
            linker_seq_list.append(ele_seq)
        elif ele == "U":
            umiseq += ele_seq
        elif ele == "T":
            poly_t_list.append(ele_seq)
        ttl_seq = ttl_seq[count:]

    return barcode, umiseq, linker_seq_list, poly_t_list, ttl_seq


def check_barcode(qry_bc_lst: list[str] | str, ref_bc_lst: list[str] | None = None, n_miss: int = 0):
    """Check if the query barcodes are in the reference barcodes."""
    if ref_bc_lst is None or len(ref_bc_lst) == 0:
        return True

    if isinstance(qry_bc_lst, str):
        qry_bc_lst = [qry_bc_lst]

    if n_miss == 0:
        for qry_bc in qry_bc_lst:
            if qry_bc not in ref_bc_lst:
                return False
    else:
        for qry_bc in qry_bc_lst:
            for ref_bc in ref_bc_lst:
                if qry_bc in ref_bc:
                    return True
                obs_miss = sum(c1 != c2 for c1, c2 in zip(ref_bc, qry_bc))
                if obs_miss <= n_miss:
                    return True
        return False

    return True


def process_chunk(in_pairs):
    """Function to process a chunk of lines"""
    global n_workers, bc_whitelist, bc_pattern
    valid_bc = True
    pr_idx_pool, pr_r1_pool, pr_r2_pool = [], [], []
    read_id, idx_seq, r1_seq, r2_seq, index_qual, r1_qual, r2_qual = "", "", "", "", "", "", ""

    for idx, (line1, line2) in enumerate(in_pairs):
        line1, line2 = line1.strip(), line2.strip()

        if idx % 4 == 0: # Read ID
            read_id = line1.split(" ").pop(0)
            if "/1" in read_id or "/2" in read_id: read_id = read_id[:-2]
        elif idx % 4 == 1: # Read sequence
            barcode, umi_seq, *_, r1_seq = decode_sequence(line1, bc_pattern)
            valid_bc = check_barcode(barcode, bc_whitelist)
            idx_seq, r1_seq, r2_seq = barcode + umi_seq, r1_seq.translate(TRANS_DICT)[::-1], line2
        elif idx % 4 == 3: # Read quality
            barcode_qual, umi_qual, *_, r1_qual = decode_sequence(line1, bc_pattern)
            index_qual, r1_qual, r2_qual = barcode_qual + umi_qual, r1_qual[::-1], line2

            non_empty = all([x != "" for x in [idx_seq, r1_seq, r2_seq, index_qual, r1_qual, r2_qual]])
            if valid_bc and non_empty:
                pr_idx_pool.append(f"{read_id}/1\n{idx_seq}\n+\n{index_qual}\n")
                pr_r1_pool.append(f"{read_id}/2\n{r1_seq}\n+\n{r1_qual}\n")
                pr_r2_pool.append(f"{read_id}/2\n{r2_seq}\n+\n{r2_qual}\n")

            # Reset variables
            valid_bc, read_id, idx_seq, r1_seq, r2_seq, index_qual, r1_qual, r2_qual = True, "", "", "", "", "", "", ""

    return pr_idx_pool, pr_r1_pool, pr_r2_pool


def processor(
    in_r1: Path, in_r2: Path, out_r1: Path, out_r2: Path, out_idx: Path, n_workers: int = 4, chunk_size: int = 100000,
    batch_size: int = 5000, gz_threads: int = 2
):
    """Process the fastq file."""
    with (
        pgzip.open(in_r1, "rt", thread=gz_threads) as r1_hand, pgzip.open(in_r2, "rt", thread=gz_threads) as r2_hand,
        pgzip.open(out_idx, "wt", thread=gz_threads) as wi_hand,
        pgzip.open(out_r1, "wt", thread=gz_threads) as w1_hand, pgzip.open(out_r2, "wt", thread=gz_threads) as w2_hand,
    ):
        chunk_size, batch_size = chunk_size * 4, batch_size * 4
        pair_vec, pair_batch, pair_chunk = [], [], []
        for idx, per_pair in tqdm.tqdm(enumerate(zip(r1_hand, r2_hand))):
            if idx % batch_size == 0 and (idx + 1) >= batch_size:
                pair_chunk.append(pair_vec)
                pair_vec = []

            if idx % chunk_size == 0 and (idx + 1) >= chunk_size:
                with Pool(n_workers) as pool:
                    results = pool.map(process_chunk, pair_chunk)
                pool.join()

                for idx_pc, r1_pc, r2_pc in results:
                    wi_hand.writelines(idx_pc)
                    w1_hand.writelines(r1_pc)
                    w2_hand.writelines(r2_pc)
                pair_chunk = []

            pair_vec.append(per_pair)
        else:
            pair_chunk.append(pair_vec)
            with Pool(n_workers) as pool:
                results = pool.map(process_chunk, pair_chunk)
            pool.join()

            for idx_pc, r1_pc, r2_pc in results:
                wi_hand.writelines(idx_pc)
                w1_hand.writelines(r1_pc)
                w2_hand.writelines(r2_pc)

        del pair_chunk, pair_batch, pair_vec


def main():
    """Preprocess slide-seq reads."""
    global bc_whitelist, bc_pattern, n_workers, rev_cmp_r1, rev_cmp_r2

    parser = ArgumentParser(description="Preprocess slide-seq reads")
    parser.add_argument("read1", help="Input fastq file R1 (barcode and UMI). Required")
    parser.add_argument("read2", help="Input fastq file R2 (molecules). Required")
    parser.add_argument("-c", "--chunk-size", metavar="INT", type=int, default=2000000, help="Size per chunk. Default: %(default)s")
    parser.add_argument("-b", "--batch-size", metavar="INT", type=int, default=100000, help="Size per Batch. Default: %(default)s")
    parser.add_argument("-n", "--sample-name", metavar="STR", default="sample", help="Sample name. Default: %(default)s")
    parser.add_argument("-p", "--pattern", metavar="STR", default="C4L15C4L15C4U10T18", help="Pattern of R1 reads. Default: %(default)s")
    parser.add_argument("-w", "--white-list", metavar="PATH", default=None, help="White list file. Default: %(default)s")
    parser.add_argument("--rev-cmp-r1", action="store_true", default=False, help="Reverse complement R1 reads. Default: %(default)s")
    parser.add_argument("--rev-cmp-r2", action="store_true", default=False, help="Reverse complement R2 reads. Default: %(default)s")
    parser.add_argument("--white-list-sep", metavar="STR", default=",", help="Field separator in white list file. Default: %(default)s")
    parser.add_argument("--barcode-col", metavar="INT", type=int, default=3, help="Barcode's column index (1-based). Default: %(default)s")
    parser.add_argument("-N", "--n-workers", metavar="INT", type=int, default=1, help="Nr. of parallel workers. Default: %(default)s")
    parser.add_argument("-Z", "--gzip-threads", metavar="INT", type=int, default=1, help="Nr. of pgzip threads. Default: %(default)s")
    parser.add_argument("-o", "--out-dir", metavar="PATH", default="./", help="Path to store the processed FASTQ. Default: %(default)s")

    args = parser.parse_args()

    chunk_size, batch_size, sample_name, white_list = args.chunk_size, args.batch_size, args.sample_name, args.white_list
    gz_threads, n_workers = args.gzip_threads, args.n_workers
    rev_cmp_r1, rev_cmp_r2 = args.rev_cmp_r1, args.rev_cmp_r2
    in_r1, in_r2, out_dir  = Path(args.read1), Path(args.read2), Path(args.out_dir)

    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    bc_pattern = BARCODE_PATTERN.findall(args.pattern)
    bc_whitelist = []
    if white_list:
        white_list_sep = args.white_list_sep
        barcode_col = args.barcode_col - 1
        with open(white_list, "r") as fwl:
            bc_whitelist = list(set([x.strip().split(white_list_sep)[barcode_col] for x in fwl]))
    
    out_idx, out_r1, out_r2 = out_dir / f"{sample_name}.I1.fastq.gz", out_dir / f"{sample_name}.R1.fastq.gz", out_dir / f"{sample_name}.R2.fastq.gz"
    processor(in_r1, in_r2, out_r1, out_r2, out_idx, n_workers, chunk_size, batch_size, gz_threads)


if __name__ == "__main__":
    main()

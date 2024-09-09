#!/usr/bin/env python3
# File: t2c_counter.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jul 17, 2024
# Updated: Jul 17, 2024

import os, logging, pathlib
from multiprocessing import Pool
from multiprocessing.pool import ApplyResult
from argparse import ArgumentParser

import pysam
import pandas as pd


BASE_REVCOMP_DICT = {"a": "t", "c": "g", "t": "a", "g": "c", "n": "n"}


class LogManager(logging.Logger):
    def __init__(self, name, level=logging.INFO, logstream: bool = True, logfile: str = ""):
        super(LogManager, self).__init__(name)
        fmt = logging.Formatter("{levelname: >8}|{asctime}|{name: >8}| {message}", style="{", datefmt="%Y%m%d,%H%M%S")
        if logstream:
            self._add_handler(logging.StreamHandler(), level, fmt)

        if logfile:
            self._add_handler(logging.FileHandler(logfile), level, fmt)

    def _add_handler(self, hdl, lvl, fmt):
        hdl.setLevel(lvl)
        hdl.setFormatter(fmt)
        self.addHandler(hdl)


def get_cliopts():
    parser = ArgumentParser(description="A program to count Ti/Tv on each read.")
    parser.add_argument("bam_file_list", nargs="+", metavar="BAMFILE")

    parser.add_argument("-r", "--region", metavar="CHROM[:START[-END]]", default=None, help="Only count a given region, all by default (None). Default: %(default)s")
    parser.add_argument("-R", "--region-file", metavar="FILE", default=None, help="File contains regions to be processed in BED/GTF/GFF format. Default: %(default)s")
    parser.add_argument("-g", "--gff-file", metavar="FILE", default=None, help="Genomic feature files to group the Ti/Tv mutants. Default: %(default)s")
    parser.add_argument("-b", "--black-list", metavar="FILE", default=None, help="Positions should be removed. A file in VCF or BED format. Default: %(default)s")
    parser.add_argument("-@", "--n-cpus", metavar="INT", default=1, type=int, help="Number of parallel jobs. Default: %(default)s")

    parser.add_argument("-e", "--excl-flag", metavar="INT", default=1024, type=int, help="Flags to be removed. Default: %(default)s")
    parser.add_argument("-i", "--incl-flag", metavar="INT", default=0, type=int, help="Flags to be included. Default: %(default)s")
    parser.add_argument("-q", "--min-quality", metavar="INT", default=20, type=int, help="Minmum quality. Default: %(default)s")

    parser.add_argument("-O", "--one-sample-per-bam", action="store_true", help="Indicating one sample per BAM file.")
    parser.add_argument("-o", "--out-file", metavar="FILE", default="titv_counts.txt", help="Output file. Default: %(default)s")

    return parser.parse_args()


def crate_blacklist(blfile: str):
    blfile_path = pathlib.Path(blfile)

    if blfile.endswith(""):
        itvl = pysam.VariantFile(blfile)
    elif blfile.endswith(".bed"):
        itvl = pysam.TabixFile(blfile)


def count_titv(bamfile, region=None, black_list=None, excl_flag=1540, min_qual=20, max_mm=10, trim3end=10, trim5end=10, logman: LogManager = LogManager("TiTvCounter")):
    """Count transitions (Ti) and tranversions (Tv) on each read."""
    worker_pid = os.getpid()
    logman.info(f"Worker {worker_pid} is running ... ")
    black_list = [] if black_list is None else black_list
    reads_reports = {
        "Overlapped reads": 0, "Overlapped SE reads": 0, "Overlapped PE reads": 0, "Removed (excl-flag)": 0,
        "Removed (too many mismatches)": 0, "Skipped (no mismatches)": 0
    }
    titv_pool = []

    bam = pysam.AlignmentFile(bamfile, mode="r", duplicate_filehandle=True)
    reads_reports["Total reads"] = bam.mapped + bam.unmapped
    reads_reports["Mapped reads"], reads_reports["Unmapped reads"] = bam.mapped, bam.unmapped
    for per_rec in bam.fetch(region=region):
        # Collect basic statistis for the current bam
        reads_reports["Overlapped reads"] += 1
        reads_reports["Overlapped PE reads"] += per_rec.is_paired
        reads_reports["Overlapped SE reads"] += not per_rec.is_paired

        # Remove reads base on the flag parameter.
        if per_rec.flag & excl_flag > 0:
            reads_reports["Removed (excl-flag)"] += 1
            continue

        # Remove reads without any mismatch.
        if per_rec.get_tag("NM") < 1:
            reads_reports["Skipped (no mismatches)"] += 1
            continue

        refer_id, query_id = per_rec.reference_name, per_rec.query_name
        query_len, query_seq, query_qual = per_rec.query_length, per_rec.query_sequence, per_rec.query_qualities

        nm_trimmed = 0 # Count number of mismatches after trimming the reads.
        per_read_titv_pool = [query_id, refer_id]
        aligned_pairs = per_rec.get_aligned_pairs(with_seq=True) # Obtain matched pairs.
        aligned_pairs = aligned_pairs[trim5end:query_len - trim3end] # Trim the 3/5-prime ends.
        is_reverse = per_rec.is_reverse

        for query_pos, refer_pos, refer_base in aligned_pairs:
            if refer_base is None or refer_base not in "atcgn": continue # Skip non-mismatches
            if (refer_id, refer_pos) in black_list: continue # Skip mismatches in black list, e.g., SNPs
            if query_qual[query_pos] <= min_qual: continue # Skip mismatches due to low base quality

            query_base = query_seq[query_pos]
            if is_reverse:
                refer_base = BASE_REVCOMP_DICT.get(refer_base, "n")
                query_base = BASE_REVCOMP_DICT.get(query_base.lower(), "N").upper()
            base_trans = f"{refer_base}>{query_base}"
            per_read_titv_pool.append([refer_pos, base_trans, query_qual[query_pos]])
            nm_trimmed += 1

        if 0 < nm_trimmed <= max_mm: # Skipping reads due to too many mismatches after trimming. TODO: using threshold
            titv_pool.append(per_read_titv_pool)

    try:
        bam.close()
    except Exception as e:
        logman.warn(f"Failed to close bam file, worker pid ({worker_pid})")
        logman.warn(f"{e}")

    reads_reports["Overlapped PE reads"] *= 2
    logman.info(f"Worker {worker_pid} was done")
    return bamfile, reads_reports, titv_pool


def collect_results(reslist, one_sample_per_bam=True, logman: LogManager = LogManager("ResultCollector")):
    titv_res_list, stat_res_list, new_rna_list = [], [], []
    all_trans = [f"{x}>{y}" for x in "atcgn" for y in "ATCGN" if x.upper() != y]
    bam_file_list = []
    for perres in reslist:
        bam_file, stat_list, titv_list = [None] * 3

        try:
            bam_file, stat_list, titv_list = perres.get()
        except Exception as e:
            logman.info(f"{e}")

        if titv_list is None or stat_list is None: continue

        temp_titv_count_dict = {per_trans:{} for per_trans in all_trans}
        tmp_newrna_dict = {}
        for query_id, refer_id, *per_titv in titv_list:
            n_t2c, t2c_pos_list = 0, []
            for refer_pos, titv, _ in per_titv:
                cur_idx = (refer_id, refer_pos)
                if cur_idx in temp_titv_count_dict[titv]:
                    temp_titv_count_dict[titv][cur_idx] += 1
                else:
                    temp_titv_count_dict[titv][cur_idx] = 1

                if titv == "t>C":
                    n_t2c += 1
                    t2c_pos_list.append(":".join([refer_id, str(refer_pos)]))

            if n_t2c > 0:
                tmp_newrna_dict[query_id] = "|".join(t2c_pos_list)
        
        new_rna_list.append(pd.DataFrame({"read_id": tmp_newrna_dict.keys(), "positions": tmp_newrna_dict.values()}).assign(Datasource=bam_file))
        titv_res_list.append(pd.DataFrame(temp_titv_count_dict).assign(Datasource=bam_file).fillna(0))
        stat_res_list.append(pd.DataFrame({"Statistics": stat_list.keys(), "Value": stat_list.values()}).assign(Datasource=bam_file))
        bam_file_list.append(bam_file)

    titv_count_tab = pd.concat(titv_res_list).convert_dtypes()
    stat_res_tab = pd.concat(stat_res_list).convert_dtypes()
    new_rna_tab = pd.concat(new_rna_list).convert_dtypes()

    if not one_sample_per_bam:
        data_source = "|".join(bam_file_list)
        titv_count_tab.loc[:, "Datasource"] = data_source
        stat_res_tab.loc[:, "Datasource"] = data_source

    return stat_res_tab, titv_count_tab, new_rna_tab


def main(logman: LogManager = LogManager("Main")):
    opts = get_cliopts()
    n_cpus = opts.n_cpus
    region = opts.region
    one_sample_per_bam = opts.one_sample_per_bam

    out_file = opts.out_file
    bam_file_list = opts.bam_file_list

    black_list = opts.black_list
    logman.info("Creating black list from")

    logman.info("Collecting Ti/Tv for each read/segment ...")
    kwargs = {"region": region, "black_list": black_list}
    pool = Pool(n_cpus)
    results = [pool.apply_async(count_titv, args=(bamfile,), kwds=kwargs) for bamfile in bam_file_list]
    pool.close()
    pool.join()

    logman.info("Saving results into disk ...")
    stat_tab, titv_tab, new_rna_tab = collect_results(results, one_sample_per_bam)
    titv_tab.to_csv(out_file)


if __name__ == "__main__":
    main()

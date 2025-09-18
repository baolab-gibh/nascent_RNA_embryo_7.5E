#!/usr/bin/env python3
# File: extract_new_reads.py
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Dec 30, 2024
# Updated:
"""A script to extract new reads from 4sU labeled sequencing results."""


import re
import pathlib
import logging
import functools
import copy
# from typing import Self

import click
import pysam
import pandas as pds
import polars as pls


MD_TAG_REGEX = re.compile("[0-9]+|[A-Z]|\\^[A-Z]+")
CIGAR_OPS_REGEX = re.compile("[0-9]+[MIDNSHP=X]")

CIGAR_OPS = "MIDNSHP=X"
TRANSLATOR = str.maketrans("ATCGNatcgn^", "TAGCNtagcn^")

OUT_COLUMNS = [
    "ReadIdx", "GroupBy", "ReadSkipCode", "MappedChrom", "MappedStart", "MappedEnd", "Position", "RefBase", "AltBase",
    "AltQual", "BaseSkipCode", "IsT2C"
]


class LogManager(logging.Logger):
    def __init__(self, name, level=logging.INFO, logstream: bool = True, logfile: str = ""):
        super(LogManager, self).__init__(name)
        fmt = logging.Formatter(
            "{levelname: >8} | {asctime} | {name: ^20} | {message}", style="{", datefmt="%y-%m-%d,%H:%M:%S"
        )
        if logstream:
            self._add_handler(logging.StreamHandler(), level, fmt)

        if logfile:
            self._add_handler(logging.FileHandler(logfile), level, fmt)

    def _add_handler(self, hdl, lvl, fmt):
        hdl.setLevel(lvl)
        hdl.setFormatter(fmt)
        self.addHandler(hdl)


class SamAlignedSegment:
    def __init__(
        self, segment: pysam.AlignedSegment, header: pysam.AlignmentHeader | None = None, min_qual: int = -1,
        by_umi: bool = True, phred: int = 33, logman: LogManager = LogManager("Read")
    ) -> None:
        self._logman = logman
        if (segment.query_sequence is None or segment.query_qualities is None or segment.cigarstring is None 
            or segment.cigartuples is None or segment.get_tag("MD") is None
        ): return

        self._read = segment
        self._header = header

        self._chrom_id = segment.reference_id if segment.reference_id else -1
        self._chrom = segment.reference_name if segment.reference_name else ''
        self._start = segment.reference_start if segment.reference_start else -1
        self._end = segment.reference_end if segment.reference_end else -1

        self._min_qual = min_qual
        self._phred = phred

        self._by_umi = by_umi
        self._is_inverted = False

        if segment.has_tag("CR"):
            self._cb = str(segment.get_tag("CR"))
        else:
            self._cb = None
        if segment.has_tag("UR"):
            self._umi = str(segment.get_tag("UR"))
        else:
            self._cb = None

        self._md_string = self._read.get_tag("MD")
        self._md_list = MD_TAG_REGEX.findall(self._md_string)

        self._mm_info = []
        self._collect_mismatches()

    def __str__(self) -> str:
        example_bases = ''
        if self._read.query_sequence:
            example_bases = self._read.query_sequence[:10] if self._read.query_sequence else ''

        example_quals = ''
        if self._read.query_qualities:
            example_quals = ''.join([chr(x + self._phred) for x in self._read.query_qualities[:10]])

        return f"{self._chrom}:{self._start}-{self._end}: {example_bases}... {example_quals}..."

    def __expr__(self):
        return self.__str__()

    def __eq__(self, other):
        return self._chrom == other._chrom and self._start == other._start

    def __gt__(self, other):
        return self._chrom > other._chrom if self._chrom != other._chrom else self._start > other._start

    def __ge__(self, other):
        return self == other or self > other

    def __lt__(self, other):
        return self._chrom < other._chrom if self._chrom != other._chrom else self._start < other._start

    def __le__(self, other):
        return self == other or self < other

    def __getitem__(self, i):
        '''Subsetting a Read.'''
        if self._read.query_sequence is None or self._read.query_qualities is None: return None

        if isinstance(i, slice):
            return self._read.query_sequence[i], ''.join([chr(x + self._phred) for x in self._read.query_qualities[i]])

        return None, None

    def __add__(self, other): # TODO
        if True: return

        # Update flag
        # Merge qualities
        # Merge sequences. Mutations, deletions, insertions.
        # Update CIGAR string
        # Update MD string
        # return self, other

    def __invert__(self):
        """Obtain the inverse complement of the read.

        The operation will deep-copy the self._read, which is time and memory cost.
        """
        new_read = copy.deepcopy(self._read)
        new_quality = new_read.query_qualities[::-1]
        new_read.query_sequence = new_read.query_sequence[::-1].translate(TRANSLATOR)
        new_read.query_qualities = new_quality

        # Update the flag
        new_read.flag ^= 16

        # Invert CIGAR string
        cigar_tuples = new_read.cigartuples if new_read.cigartuples else []
        if len(cigar_tuples) > 1:
            new_cigar = ""
            for ops, length in cigar_tuples:
                new_cigar = f"{length}{CIGAR_OPS[ops]}" + new_cigar
            new_read.cigarstring = new_cigar

        # Invert MD tag
        md_tag = str(new_read.get_tag("MD")) if new_read.has_tag("MD") else ""
        if md_tag:
            new_md_tag = ""
            for ele in self._md_list:
                new_md_tag = ele + new_md_tag if ele.isdigit() else ele.translate(TRANSLATOR) + new_md_tag

            if new_md_tag: new_read.set_tag("MD", new_md_tag, "Z", True)

        my_read = SamAlignedSegment(new_read, self._header, self._min_qual, self._by_umi, self._phred)
        my_read._is_inverted = False if my_read._is_inverted else True

        return my_read

    def _collect_mismatches(self):
        insert_list = []
        if "I" in self._read.cigarstring: # We have to deal with insertions due to MD tag does not contains insertions.
            qry_bases, qry_quals, qry_pos = [], [], 0
            for op, length in self._read.cigartuples:
                if op in [0, 1, 4]:
                    if op == 1: insert_list.append((qry_pos, length))
                    qry_pos += length
                elif op in [2, 3, 5]:
                    continue
                else:
                    raise ValueError(f"Unsupported CIGAR op: {op}")
            if qry_quals: qry_quals = array.array("B", qry_quals)

        ref_md_pos = self._read.reference_start
        qry_bases = list(self._read.query_sequence)
        qry_quals = self._read.query_qualities

        new_md_list, qry_md_pos = [], 0
        for token in self._md_list:
            if token.isdigit():
                x_int = int(token)
                qry_md_pos += x_int
                ref_md_pos += x_int

                if len(new_md_list) == 0:
                    new_md_list.append(x_int)
                elif isinstance(new_md_list[-1], int):
                    new_md_list[-1] += x_int
                else:
                    new_md_list.append(x_int)
            elif token in 'ATCGNatcgn':
                offset = 0
                for insert_pos, insert_len in insert_list:
                    if qry_md_pos >= insert_pos:
                        offset += insert_len
                qry_base, qry_qual = qry_bases[qry_md_pos + offset], qry_quals[qry_md_pos + offset]
                adjusted = qry_qual <= self._min_qual

                if adjusted: # Adjuste the current base if it's sequencing quality is too low
                    qry_bases[qry_md_pos] = token
                    qry_quals[qry_md_pos] = self._min_qual

                    if len(new_md_list) == 0:
                        new_md_list.append(1)
                    elif isinstance(new_md_list[-1], int):
                        new_md_list[-1] += 1
                    else:
                        new_md_list.append(token)
                else:
                    new_md_list.append(token)

                qry_md_pos += 1
                ref_md_pos += offset + 1
                ref_pos, qry_pos, ref_base = ref_md_pos, qry_md_pos, token
                self._mm_info.append([ref_pos, ref_base, qry_pos, qry_base, qry_qual, adjusted])
            elif "^" in token:
                ref_md_pos += len(token) - 1
                new_md_list.append(token)
            else:
                raise ValueError(f"Unknown MD operation: {token}")

        if new_md_list:
            new_md_list = [str(x) for x in new_md_list]
            new_md_string = "".join(new_md_list)
            self._read.set_tag("MD", new_md_string, "Z", True)
            self._read.query_sequence = "".join(qry_bases).replace("-", "")
            self._read.query_qualities = qry_quals

    def collect_mismatches(self, trim3p: int = 10, trim5p: int = 10, min_qual: int = 20):
        mm_info = []
        if not isinstance(self._mm_info, list) or len(self._mm_info) <= 0: return mm_info
        for ref_pos, ref_base, qry_pos, qry_base, qry_qual, adjusted in self._mm_info:
            if adjusted: continue
            if trim5p >= qry_pos or qry_pos >= self._read.query_length - trim3p: continue
            if qry_qual <= min_qual: continue
            mm_info.append((ref_pos, ref_base, qry_pos, qry_base, qry_qual))

        return mm_info

    def is_nascent(
        self, excl_flag: int, incl_flag: int, excl_unpaired: bool = False, max_mm: int = 5, min_qual: int = 20,
        bl_list: list | dict | None = None, bl_binsize: int = 100000, trim3end: int = 10, trim5end: int = 10
    ):
        assert isinstance(self._read, pysam.AlignedSegment)
        read_skip_code = 0
        if self._read.flag & excl_flag != 0: read_skip_code |= 1 # Skip due to containing non-required flag
        if excl_unpaired and self._read.flag == 0: read_skip_code |= 2 # Skip due to unpaired reads
        if self._read.flag > 0 and self._read.flag & incl_flag == 0: read_skip_code |= 4 # Skip due to missing required flags

        try:
            nm_counts = self._read.get_tag("NM")
            nm_counts = int(nm_counts)
        except:
            nm_counts = 0
        if nm_counts < 1: read_skip_code |= 8 # Remove reads without mismatches.
        if nm_counts > max_mm: read_skip_code |= 16 # Remove reads with too many mismatches

        qry_seq = self._read.query_sequence
        if qry_seq is None: read_skip_code |= 32 # Remove reads without base sequence

        qry_quals = self._read.query_qualities
        if qry_quals is None: read_skip_code |= 64 # Remove reads without base quality

        # Walk through the matches if the read is valid.
        is_valid_new, t2c_count, per_read_mismatches = False, 0, []
        if read_skip_code == 0:
            ref_id, ref_start = self._read.reference_name, self._read.reference_start # Reference sequence information
            qry_len, qry_seq = self._read.query_length, self._read.query_sequence # Query sequence information
            min_pos, max_pos = trim5end, qry_len - trim3end # Define genomic interval to walk through

            mm_info = self.collect_mismatches(trim3end, trim5end, min_qual)
            for ref_pos, ref_base, qry_pos, qry_base, qry_qual in mm_info:
                base_skip_code = 0
                if qry_pos < min_pos: base_skip_code |= 1 # Skip the 3/5-prime ends, code by 1
                if qry_pos >= max_pos: base_skip_code |= 2 # Skip the 3/5-prime ends, code by 10
                if ref_base is None: base_skip_code |= 4 # Skip due to short deletion, coded by 100
                if ref_base not in "ATCGNatcgn": base_skip_code |= 8 # Skip due to non-canonical bases, coded by 1000
                if qry_qual <= min_qual: base_skip_code |= 16 # Skip mismatches due to low mapping quality, coded by 10000
                if in_blacklist(bl_list, ref_id, ref_pos, bl_binsize): base_skip_code |= 32 # Skip the base if it is in the black list, coded by 100000

                is_t2c = (ref_base in "Aa" and qry_base in "Gg") or (ref_base in "Tt" and qry_base in "Cc")
                if is_t2c and base_skip_code == 0: t2c_count += 1
                per_read_mismatches.append([ref_pos, ref_base, qry_pos, qry_base, qry_qual])
            is_valid_new = 0 < t2c_count

        return is_valid_new, int(read_skip_code), per_read_mismatches


class SamAlignedSegmentPool:
    def __init__(
        self, group_by: str | None = None, logman: LogManager = LogManager("AlignedSegmentPool")
    ) -> None:
        self.logman = logman
        self.aligned_segment_list = []
        self._group_by = group_by

    def __len__(self):
        return len(self.aligned_segment_list)

    def __iter__(self):
        return iter(self.aligned_segment_list)

    def __getitem__(self, k):
        return self.aligned_segment_list[k]

    def __contains__(self, o):
        pass

    def dedup(self, how: str = "random", group_by: str | None = None):
        if how not in ["random", "longest", "merge"]:
            self.logman.error(f"Invalid dedup method: {how}. Using default method: random.")
            how = "random"
        group_by = group_by if self._group_by is None else self._group_by

    def append(self, aligned_segment: SamAlignedSegment) -> None:
        self.aligned_segment_list.append(aligned_segment)

    def pop(self) -> SamAlignedSegment:
        return self.aligned_segment_list.pop()

    def sort(self, **kwargs):
        self.aligned_segment_list.sort(**kwargs)
        return self

    def process(self):
        for per_segment in self.aligned_segment_list:
            per_segment.collect_mismatches()


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


def create_blacklist(blfile_list, bin_size):
    """Create a blacklist."""
    black_list = {}

    for blfile in blfile_list:
        if blfile is None: continue
        if not isinstance(blfile, pathlib.Path): blfile = pathlib.Path(blfile)

        if blfile.suffix == ".vcf" or "".join(blfile.suffixes[:-2]) == ".vcf.gz":
            with pysam.VariantFile(str(blfile)) as itvl_list:
                for rec in itvl_list:
                    chrom_key, start_key, pos = rec.chrom, int(rec.pos / bin_size), rec.pos
                    if chrom_key in black_list:
                        if start_key in black_list[chrom_key]:
                            black_list[chrom_key][start_key].append(pos)
                        else:
                            black_list[chrom_key][start_key] = [pos]
                    else:
                        black_list[chrom_key] = {start_key: [pos]}
        elif blfile.suffix == ".bed" or "".join(blfile.suffixes[:-2]) == ".bed.gz":
            with pysam.TabixFile(str(blfile), parser=pysam.asBed()) as itvl_list:
                for rec in itvl_list.fetch():
                    chrom_key, start, end = rec.contig, rec.start, rec.end
                    pos_list = list(range(start, end))
                    start_key = int(start / bin_size)

                    if chrom_key in black_list:
                        if start_key in black_list[chrom_key]:
                            black_list[chrom_key][start_key].extend(pos_list)
                        else:
                            black_list[chrom_key][start_key] = pos_list
                    else:
                        black_list[chrom_key] = {start_key: pos_list}
        else:
            raise ValueError("Unrecognized blacklist format: %s" % blfile)

    return black_list


def in_blacklist(blacklist, chrom, pos, bin_size):
    if blacklist is None: return False
    if chrom in blacklist:
        pos_key = int(pos / bin_size)
        if pos_key in blacklist[chrom]:
            return pos in blacklist[chrom][pos_key]
    return False


def decode_skip_reason(code, mode="read"):
    """Decode skipping reasons from code."""
    reasons = []
    if mode == "read":
        if code & 1 == 1: reasons.append("Containing non-required flag")
        if code & 2 == 2: reasons.append("Is unpaired reads")
        if code & 4 == 4: reasons.append("Missing required flags")
        if code & 8 == 8: reasons.append("Missing mismatches")
        if code & 16 == 16: reasons.append("Too many mismatches")
        if code & 32 == 32: reasons.append("Missing base sequence")
        if code & 64 == 48: reasons.append("Missing base quality")
    else:
        if code & 1 == 1: reasons.append("3-prime ends")
        if code & 2 == 2: reasons.append("5-prime ends")
        if code & 4 == 4: reasons.append("Short deletion")
        if code & 8 == 8: reasons.append("Non-canonical bases")
        if code & 16 == 16: reasons.append("Low mapping quality")
        if code & 32 == 32: reasons.append("In black list")

    return reasons


@click.command()
@click.argument("in_bam_file", metavar="<in.bam>")
@click.option("-o", "--old-reads-out-file", metavar="FILE", show_default=True, default=None, type=str, help="Save old reads to.")
@click.option("-n", "--new-reads-out-file", metavar="FILE", show_default=True, default=None, type=str, help="Save new reads to.")
@click.option("-i", "--incl-flag", metavar="INT", show_default=True, default=4095, type=int, help="Include reads with this flag.")
@click.option("-e", "--excl-flag", metavar="INT", show_default=True, default=0, type=int, help="Exclude reads with this flag.")
@click.option("-g", "--group-tag", metavar="STR", multiple=True, show_default=True, default=[None], help="Group cells by metadata.")
@click.option("-b", "--black-list", metavar="FILE", multiple=True, show_default=True, default=[None], help="Black list file.")
@click.option("-s", "--black-list-bin-size", metavar="INT", show_default=True, default=100000, type=int, help="Black list bin size.")
@click.option("-r", "--region", metavar="STR", show_default=True, default=None, help="Region to process.")
# @click.option("-d", "--delete-duplicates", metavar="BOOL", show_default=True, default=False, is_flag=True, help="Delete duplicate reads.")
# @click.option("-D", "--dedup-by", metavar="STR", show_default=True, multiple=True, default=["CR"], type=click.Choice(["CR", "BR"]), help="Deduplication method.")
@click.option("-m", "--max-mismatches", metavar="INT", show_default=True, default=4, type=int, help="Maximum number of mismatches.")
@click.option("-q", "--min-qualities", metavar="INT", show_default=True, default=25, type=int, help="Minimum quality score.")
@click.option("-@", "--n-cpus", metavar="INT", show_default=True, default=1, help="Number of CPUs.")
@click.option("--trim-3-end", metavar="INT", show_default=True, default=3, type=int, help="Trim 3' ends.")
@click.option("--trim-5-end", metavar="INT", show_default=True, default=2, type=int, help="Trim 5' ends.")
@click.option("--excl-unpaired", metavar="BOOL", show_default=True, default=False, is_flag=True, help="Exclude unpaired reads.")
@click.option("-O", "--save-read-info", metavar="FILE", show_default=True, default=None, help="Output mismatches table.")
def main(in_bam_file, **kwargs):
    """Extract nascent reads."""
    old_reads_out_file = kwargs["old_reads_out_file"]
    new_reads_out_file = kwargs["new_reads_out_file"]

    excl_flag = kwargs["excl_flag"]
    incl_flag = kwargs["incl_flag"]

    group_tag = kwargs["group_tag"]
    black_list_file = kwargs["black_list"]
    black_list_bin_size = kwargs["black_list_bin_size"]
    region = kwargs["region"]

    n_cpus = kwargs["n_cpus"]
    n_in_cpus = max(1, int(n_cpus / 3))
    n_out_cpus = max(1, int((n_cpus - n_in_cpus) / 2))

    trim_5_end = kwargs["trim_5_end"]
    trim_3_end = kwargs["trim_3_end"]
    excl_unpaired = kwargs["excl_unpaired"]
    max_mismatches = kwargs["max_mismatches"]
    min_qualities = kwargs["min_qualities"]

    # dedup = kwargs["delete_duplicates"]
    # dedup_by = kwargs["dedup_by"]

    save_read_info = kwargs["save_read_info"]

    # Create black list
    black_list = create_blacklist(black_list_file, black_list_bin_size)

    read_info_list = []
    with pysam.AlignmentFile(in_bam_file, mode="r", threads=n_in_cpus, duplicate_filehandle=True) as in_bam:
        if new_reads_out_file is not None and old_reads_out_file is not None:
            new_out_bam = pysam.AlignmentFile(new_reads_out_file, mode="wb", threads=n_out_cpus, template=in_bam)
            old_out_bam = pysam.AlignmentFile(old_reads_out_file, mode="wb", threads=n_out_cpus, template=in_bam)
        else:
            new_out_bam, old_out_bam = None, None

        read_pool = in_bam.fetch() if region is None else in_bam.fetch(region=region)
        for per_read in read_pool:
            header = in_bam.header
            per_segment = SamAlignedSegment(per_read, header, min_qual=min_qualities)

            if new_out_bam is not None and old_out_bam is not None:
                is_new_read, read_skip_code, per_read_mismatches = per_segment.is_nascent(
                    incl_flag, excl_flag, excl_unpaired, max_mismatches, min_qualities, black_list,
                    black_list_bin_size, trim_3_end, trim_5_end
                )

                if is_new_read:
                    new_out_bam.write(per_segment._read)
                else:
                    old_out_bam.write(per_segment._read)

                if save_read_info:
                    qry_idx = per_read.query_name
                    ref_idx, ref_start, ref_end = per_read.reference_name, per_read.reference_start, per_read.reference_end
                    group_by = "NO_GROUP_TAG" if len(group_tag) == 0 else ";".join([str(per_read.get_tag(t)) for t in group_tag])

                    base_info = [qry_idx, group_by, read_skip_code, ref_idx, ref_start, ref_end]

                    if isinstance(per_read_mismatches, list) and len(per_read_mismatches) > 0:
                        read_info_list.extend([base_info + per_rec for per_rec in per_read_mismatches])
                    else:
                        read_info_list.append(base_info + [""] * 4)

        if new_out_bam: new_out_bam.close()
        if old_out_bam: old_out_bam.close()

    if save_read_info:
        col_keys = pds.Index(OUT_COLUMNS)
        all_mismatches_tab = pds.DataFrame(read_info_list, columns=col_keys)
        all_mismatches_tab.to_csv(save_read_info, index=False)


if __name__ == "__main__":
    main(max_content_width=100)

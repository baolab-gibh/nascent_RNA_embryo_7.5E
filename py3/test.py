# File: test.py

from os import wait
import copy
import time
import re
import logging
from typing import Self
import array
import tqdm

import pysam
from collections import OrderedDict

# MD_REGEX = re.compile("[0-9]+(([A-Z]|\\^[A-Z]+)[0-9]+)*")
MD_TAG_REGEX = re.compile("[0-9]+|[A-Z]|\\^[A-Z]+")
TRANSLATOR = str.maketrans("ATCGNatcgn^", "TAGCNtagcn^")
CIGAR_OPS = "MIDNSHP=X"
CIGAR_OPS_REGEX = re.compile("[0-9]+[MIDNSHP=X]")


def in_blacklist(blacklist, chrom, pos, bin_size):
    if blacklist is None: return False
    if chrom in blacklist:
        pos_key = int(pos / bin_size)
        if pos_key in blacklist[chrom]:
            return pos in blacklist[chrom][pos_key]
    return False


class LogManager(logging.Logger):
    def __init__(self, name, level=logging.INFO, logstream: bool = True, logfile: str = ""):
        super(LogManager, self).__init__(name)
        fmt = logging.Formatter("{levelname: >8}|{asctime}|{name: >12}| {message}", style="{", datefmt="%Y%m%d,%H%M%S")
        if logstream:
            self._add_handler(logging.StreamHandler(), level, fmt)

        if logfile:
            self._add_handler(logging.FileHandler(logfile), level, fmt)

        self._level = level

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
        self._cb = str(segment.get_tag("CR"))
        self._umi = str(segment.get_tag("UR"))
        self._is_inverted = False

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

    def __gt__(self, other: Self):
        return self._chrom > other._chrom if self._chrom != other._chrom else self._start > other._start

    def __ge__(self, other: Self):
        return self == other or self > other

    def __lt__(self, other: Self):
        return self._chrom < other._chrom if self._chrom != other._chrom else self._start < other._start

    def __le__(self, other: Self):
        return self == other or self < other

    def __getitem__(self, i):
        '''Subsetting a Read.'''
        if self._read.query_sequence is None or self._read.query_qualities is None: return None

        if isinstance(i, slice):
            return self._read.query_sequence[i], ''.join([chr(x + self._phred) for x in self._read.query_qualities[i]])

        return None, None

    def __add__(self, other: Self):
        if True: return
        if self._chrom != other._chrom: return None
        if self._read.get_overlap(other._start, other._end) == 0: return None
        if self._by_umi and other._by_umi and (self._cb != other._cb or self._umi != other._umi): return None

        l_bases, l_quals = self._read.query_sequence, self._read.query_qualities
        r_bases, r_quals = other._read.query_sequence, other._read.query_qualities

        l_offset, r_offset = self._read.reference_start, other._read.reference_start

        # CIGAR operations
        cigar_ops_list, l_qry_pos, r_qry_pos = [], 0, 0
        for ops, length in self._read.cigartuples:
            cigar_ops_list.append((l_qry_pos, l_qry_pos + l_offset, length, CIGAR_OPS[ops], "CIGAR", "L"))
            l_qry_pos += length

        for ops, length in other._read.cigartuples:
            cigar_ops_list.append((r_qry_pos, r_qry_pos + r_offset, length, CIGAR_OPS[ops], "CIGAR", "R"))
            r_qry_pos += length

        cigar_ops_list.sort(key=lambda x: x[1])

        # MD tags
        l_qry_pos, r_qry_pos = 0, 0
        md_ops_list = []
        for token in self._md_list:
            if token.isdigit():
                token_int = int(token)
                md_ops_list.append((l_qry_pos, l_qry_pos + l_offset, token_int, "M", "MD", "L"))
                l_qry_pos += token_int
            elif "^" in token:
                md_ops_list.append((l_qry_pos, l_qry_pos + l_offset, len(token) - 1, "D", "MD", "L"))
                continue
            elif token in 'ATCGNatcgn':
                md_ops_list.append((l_qry_pos, l_qry_pos + l_offset, 1, token, "MD", "L"))
                l_qry_pos += 1
            else:
                raise ValueError(f"Unknown MD token {token}")

        for token in other._md_list:
            if token.isdigit():
                token_int = int(token)
                md_ops_list.append((r_qry_pos, r_qry_pos + r_offset, token_int, "M", "MD", "R"))
                r_qry_pos += token_int
            elif "^" in token:
                md_ops_list.append((r_qry_pos, r_qry_pos + r_offset, len(token) - 1, "D", "MD", "R"))
                continue
            elif token in 'ATCGNatcgn':
                md_ops_list.append((r_qry_pos, r_qry_pos + r_offset, 1, token, "MD", "R"))
                r_qry_pos += 1
            else:
                raise ValueError(f"Unknown MD token {token}")

        md_ops_list.sort(key=lambda x: x[1])

        # Update flag
        # Merge qualities
        # Merge sequences. Mutations, deletions, insertions.
        # Update CIGAR string
        # Update MD string

        return self, other

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

                if adjusted:
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
        if not isinstance(self._mm_info, list) or len(self._mm_info) <= 0: return None

        mm_info = []
        for ref_pos, ref_base, qry_pos, qry_base, qry_qual, adjusted in self._mm_info:
            if adjusted: continue
            if trim5p >= qry_pos or qry_pos >= self._read.query_length - trim3p: continue
            if qry_qual <= min_qual: continue
            mm_info.append((ref_pos, ref_base, qry_pos, qry_base, qry_qual))

        return mm_info

    def is_nascent(
        self, excl_flag: int, incl_flag: int, excl_unpaired: bool = False, max_mm: int = 5, min_qual: int = 20,
        bl_list: list | None = None, bl_binsize: int = 100000, trim3end: int = 10, trim5end: int = 10
    ):
        assert isinstance(self._read, pysam.AlignedSegment)
        read_skip_code = 0
        if self._read.flag & excl_flag != 0: read_skip_code |= 1 # Skip due to containing non-required flag
        if excl_unpaired and self._read.flag == 0: read_skip_code |= 2 # Skip due to unpaired reads
        if self._read.flag > 0 and self._read.flag & incl_flag == 0: read_skip_code |= 4 # Skip due to missing required flags

        try:
            nm_counts = self._read.get_tag("nM")
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


with pysam.AlignmentFile('test.bam', 'rb') as bam:
    header = bam.header
    start = time.time()
    for i, x in enumerate(bam):
        if x.cigarstring is None: continue
        r2 = SamAlignedSegment(x, header, min_qual=30, by_umi=False)
        if r2.is_nascent(0, 4095)[0]: print("yes")
    end = time.time()
    print(f"Time used: {end-start}")

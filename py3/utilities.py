import logging

class LogManager(logging.Logger):
    def __init__(self, name, level=logging.INFO, logstream: bool = True, logfile: str = ""):
        super(LogManager, self).__init__(name)
        fmt = logging.Formatter("{levelname: >8}|{asctime}|{name: >12}| {message}", style="{", datefmt="%Y%m%d,%H%M%S")
        if logstream:
            self._add_handler(logging.StreamHandler(), level, fmt)

        if logfile:
            self._add_handler(logging.FileHandler(logfile), level, fmt)

    def _add_handler(self, hdl, lvl, fmt):
        hdl.setLevel(lvl)
        hdl.setFormatter(fmt)
        self.addHandler(hdl)


class GenomicInterval:
    """Genomic intervals"""
    def __init__(self, chrom: str, start: int, end: int, strand: str = "*", logman: LogManager = LogManager("GenomicInterval")):
        if start > end:
            raise ValueError(f"Start must be less than end! {start}, {end}")
        self.logman, self.strand, self.chrom, self. start, self.end = logman, strand, chrom, start, end

    def _parse_other(self, other):
        if isinstance(other, GenomicInterval):
            return other
        elif isinstance(other, (tuple, list)):
            if len(other) == 4:
                return GenomicInterval(other[0], other[1], other[2], other[3])
            return GenomicInterval(other[0], other[1], other[2])
        elif isinstance(other, dict):
            if "strand" in other:
                return GenomicInterval(other["chrom"], other["start"], other["end"], other["strand"])
            return GenomicInterval(other["chrom"], other["start"], other["end"])
        elif isinstance(other, str):
            chrom, start_end, *strand = other.split(":")
            start, end = start_end.split("-")
            if len(strand):
                return GenomicInterval(chrom, int(start), int(end), strand[0])
            return GenomicInterval(chrom, int(start), int(end))
        else:
            raise ValueError("Unrecognized interval type: %s" % type(other))

    def __setattr__(self, key, value):
        if key not in ["logman", "chrom", "start", "end", "strand"]:
            raise ValueError(f"Unrecognized attribute: {key}")

        if key in ["start", "end"]:
            if not isinstance(value, int):
                raise ValueError(f"{key} must be an integer! {value}")
            if value < 0:
                raise ValueError(f"{key} must be >= 0! {value}")

        if key in ["strand"]:
            if not isinstance(value, str):
                raise ValueError(f"strand must be a string! {value}")
            if value not in ["-", "+", "*"]:
                raise ValueError(f"Unrecognized strand: {value}. Available: [-, +, *]")

        super().__setattr__(key, value)

    def __str__(self):
        return f"Genomic Region: {self.chrom}:{self.start}-{self.end}:{self.strand}"

    def __repr__(self):
        return f"Genomic Region: {self.chrom}:{self.start}-{self.end}:{self.strand}"

    def __eq__(self, other):
        other = self._parse_other(other)
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end

    def __le__(self, other):
        other = self._parse_other(other)
        return self.chrom == other.chrom and self.start <= other.start

    def __lt__(self, other):
        other = self._parse_other(other)
        return self.chrom == other.chrom and self.start < other.start

    def __ge__(self, other):
        other = self._parse_other(other)
        return self.chrom == other.chrom and self.start >= other.start

    def __gt__(self, other):
        other = self._parse_other(other)
        return self.chrom == other.chrom and self.start > other.start

    def __contains__(self, other):
        other = self._parse_other(other)
        return self.chrom == other.chrom and self.start >= other.start and self.end <= other.start

    def __and__(self, other):
        other = self._parse_other(other)
        if self.chrom != other.chrom:
            self.logman.warning("Can't merge intervals on different chromosomes.")
            return None

        if self.end < other.start or self.start > other.end:
            self.logman.warning("Can't merge intervals that don't overlap.")
            return None

        return GenomicInterval(self.chrom, max(self.start, other.start), min(self.end, other.end))

    def __or__(self, other):
        other = self._parse_other(other)
        if self & other:
            return GenomicInterval(self.chrom, min(self.start, other.start), max(self.end, other.end))
        return None

    def __xor__(self, other):
        other = self._parse_other(other)

        overlap = self & other
        if overlap is None: return None
        
        if overlap.start == self.start and overlap.end == self.end:
            return None
        elif overlap.start == self.start:
            return GenomicInterval(self.chrom, overlap.end, self.end)
        elif overlap.end == self.end:
            return GenomicInterval(self.chrom, self.start, overlap.start)
        else:
            return GenomicInterval(self.chrom, self.start, overlap.start), GenomicInterval(self.chrom, overlap.end, self.end)

    def __rshift__(self, shift):
        '''Move the interval to 3' end.'''
        shift = shift if self.strand in ["+", "*"] else -shift
        return GenomicInterval(self.chrom, self.start + shift, self.end + shift, self.strand)

    def __irshift__(self, shift):
        '''Move the interval to 3' end. Modify the interval in place'''
        shift = shift if self.strand in ["+", "*"] else -shift
        self.start += shift
        self.end += shift

    def __lshift__(self, shift):
        '''Move the interval to 5' end.'''
        shift = shift if self.strand in ["+", "*"] else -shift
        return GenomicInterval(self.chrom, self.start - shift, self.end - shift, self.strand)

    def __ilshift__(self, shift):
        '''Move the interval to 5' end. Modify the interval in place'''
        shift = shift if self.strand in ["+", "*"] else -shift
        self.start -= shift
        self.end -= shift

    def __len__(self):
        '''Return the length of the interval'''
        return self.end - self.start

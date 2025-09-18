#!/usr/bin/env awk
# File: collect_read_alignment_overview.awk
# Author: Zhenhua Zhang

function get_digit(string, type) {
  if (type == "p") { return match(string, /([0-9.]+%$)/, DIGIT_STR) ? DIGIT_STR[1] : "NA" }
  else if (type == "f") { return match(string, /([0-9.]+$)/, DIGIT_STR) ? DIGIT_STR[1] : "NA" }
  else if (type == "n") { return match(string, /([0-9]+$)/, DIGIT_STR) ? DIGIT_STR[1] : "NA" }
  return "NA"
}

BEGIN {
  kstr = "FileName NrIR AvgLenIR NrUniMR PcUniMR NrUniSplMR RtMmBaseMR RtDelBaseMR AvgDelLenMR RtInsBaseMR AvgInsLenMR"
  kstr = kstr" NrMulLocMR PcMulLocMR NrTMLociMR PcTMLociMR NrMmUMR PcMmUMR NrShortUMR PcShortUMR NrOtherUMR PcOtherUMR"
  split(kstr, klist, " ")
  for (x in klist) {
    if (x == length(klist)) {
      print klist[x]
    } else {
      printf klist[x]","
    }
  }
}

FNR == 1 { RESULT["FileName"] = FILENAME } # Input file name
$0 ~ /Number of input reads/ { RESULT["NrIR"] = get_digit($0, "n") } # Input reads
$0 ~ /Average input read length/ { RESULT["AvgLenIR"] = get_digit($0, "f") } # Input Reads
$0 ~ /Uniquely mapped reads number/ { RESULT["NrUniMR"] = get_digit($0, "n") } # Uniquely mapped reads
$0 ~ /Uniquely mapped reads % / { RESULT["PcUniMR"] = get_digit($0, "p") } # Uniquely mapped reads
$0 ~ /Number of splices:/ { RESULT["NrUniSplMR"] = get_digit($0, "n") } # Uniquely mapped over splicing junction
$0 ~ /Mismatch rate per base, %/ { RESULT["RtMmBaseMR"] = get_digit($0, "p") } # Mismatches of uniquely mapped reads
$0 ~ /Deletion rate per base/ { RESULT["RtDelBaseMR"] = get_digit($0, "p") } # Mismatches of uniquely mapped reads
$0 ~ /Deletion average length/ { RESULT["AvgDelLenMR"] = get_digit($0, "f") } # Mismatches of uniquely mapped reads
$0 ~ /Insertion rate per base/ { RESULT["RtInsBaseMR"] = get_digit($0, "p") } # Mismatches of uniquely mapped reads
$0 ~ /Insertion average length/ { RESULT["AvgInsLenMR"] = get_digit($0, "f") } # Mismatches of uniquely mapped reads
$0 ~ /Number of reads mapped to multiple loci/ { RESULT["NrMulLocMR"] = get_digit($0, "n") } # Reads mapped to multiple loci
$0 ~ /% of reads mapped to multiple loci/ { RESULT["PcMulLocMR"] = get_digit($0, "p") } # Reads mapped to multiple loci
$0 ~ /Number of reads mapped to too many loci/ { RESULT["NrTMLociMR"] = get_digit($0, "n") } # Reads mapped to multiple loci
$0 ~ /% of reads mapped to too many loci/ { RESULT["PcTMLociMR"] = get_digit($0, "p") } # Reads mapped to multiple loci
$0 ~ /Number of reads unmapped: too many mismatches/ { RESULT["NrMmUMR"] = get_digit($0, "n") } # Unmapped reads
$0 ~ /% of reads unmapped: too many mismatches/ { RESULT["PcMmUMR"] = get_digit($0, "p") } # Unmapped reads
$0 ~ /Number of reads unmapped: too short/ { RESULT["NrShortUMR"] = get_digit($0, "n") } # Unmapped reads
$0 ~ /% of reads unmapped: too short/ { RESULT["PcShortUMR"] = get_digit($0, "p") } # Unmapped reads
$0 ~ /Number of reads unmapped: other/ { RESULT["NrOtherUMR"] = get_digit($0, "n") } # Unmapped reads
$0 ~ /% of reads unmapped: other/ { RESULT["PcOtherUMR"] = get_digit($0, "p") } # Unmapped reads

ENDFILE {
  for (x in klist) {
    if (x == length(klist)) {
      print RESULT[klist[x]] 
    } else {
      printf RESULT[klist[x]]"," 
    }
  }
}

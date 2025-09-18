#!/usr/bin/env bash


find /home/zzhang/Documents/projects/wp_vasaseq/inputs/2024* -name *2.fq.gz \
  | grep -e NC_L03_2 -e 161B -e 167NC -e BL0805_L03_R2 \
  | xargs -I% -P 2 bash -c 'zcat % | awk -v NAME=% "\$0 ~ /AAGCAGTGGTATCAA/ { COUNT++ } END {print NR, COUNT, NAME}"' % \
  >| reads_with_TSO.AAGCAGTGGTATCAA.txt

AAGCAGTGGTATCAA
FileName	Total	TSOCounts	SampleName
20241022	323550968	39220935	167NC
20241104	235988488	40794509	161B
20240813	967005300	150462623	test-A-BL0805_L03_R2.fq.gz
20241112	417223724	72627290	Lib-1-NC_L03_2.fq.gz


find /home/zzhang/Documents/projects/wp_vasaseq/inputs/2024* -name *2.fq.gz \
  | grep -e NC_L03_2 -e 161B -e 167NC -e BL0805_L03_R2 \
  | xargs -I% -P 2 bash -c 'zcat % | awk -v NAME=% "\$0 ~ /AGATCGGAAGAGCGTCGTGT/ { COUNT++ } END {print NAME, NR, COUNT}"' % \
  >| reads_with_TSO.AGATCGGAAGAGCGTCGTGT.txt

FileName                  , Total    , TSOCounts, SampleName
167NC                     , 323550968,  71347799, 20241022
161B                      , 235988488,  32307962, 20241104
test-A-BL0805_L03_R2.fq.gz, 967005300,  27106781, 20240813
Lib-1-NC_L03_2.fq.gz      , 417223724,  95662186, 20241112



for x in $(find /home/zzhang/Documents/projects/wp_vasaseq/inputs/2024* -name *2.fq.gz | grep -e NC_L03_2 -e 161B -e 167NC -e BL0805_L03_R2); do
  name=$(basename $x)
  awk -v NAME=${name%%_*} -f- <<'EOF' <(zcat $x)
  BEGIN {if (NAME == "test-A-BL0805") {print "SampleName,TotalReads,TSOReads,MultiTSOReads,SeqPrepReads,SeqPrepTSOReads"} }
NR % 4 != 2 { next }

{ 
  is_tso_reads = $1 ~ /AAGCAGTGGTATCAA/ ? 1 : 0
  is_seqprep_reads = $1 ~ /AGATCGGAAGAGCGTCGTGT/ ? 1 : 0
  is_mult_tso_reads = $1 ~ /AAGCAGTGGTATCAA.*AAGCAGTGGTATCAA/ ? 1 : 0

  if (is_tso_reads == 1) { TSOReads++ }
  if (is_seqprep_reads == 1) { SeqPrepReads++ }
  if (is_tso_reads == 1 && is_seqprep_reads == 1) { SeqPrepTSOReads++ }
  if (is_mult_tso_reads) { MultiTSOReads++ }

  Total++
}

END { OFS=","; print NAME,Total,TSOReads,MultiTSOReads,SeqPrepReads,SeqPrepTSOReads }
EOF
done

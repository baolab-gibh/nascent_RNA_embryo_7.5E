#!/bin/bash
# File: explore_slide_seq_results.sh
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Jan 02, 2025
# Updated:

project_dir=~/Documents/projects/wp_vasaseq

if [[ ${USER} == 'zhzhang_gibh' ]]; then
  source ~/tools/miniconda3/bin/activate
  mamba activate sm732 # Snakemake v7.32.4
elif [[ ${USER} =~ 'zzhang' && -d ${project_dir}/scripts/.env ]]; then
  source ${project_dir}/scripts/.env/bin/activate
fi

cd ${project_dir}/temps || return 1


# Testing the results using new trimming strategy
sample_id=2-PLY

out_dir=${project_dir}/temps/test_1
fastq_dir=${project_dir}/outputs/analysis/preprocessing/slide_seq/fastq/20241217_encoded_with_polya_embryo
fastq_r1=${fastq_dir}/${sample_id}/${sample_id}.R1.fq.gz
fastq_r2=${fastq_dir}/${sample_id}/${sample_id}.R2.fq.gz

# Reformat fastq and trimming both R1 and R2
fastq_fmt_r1=${out_dir}/${sample_id}.paired_end.R1.fq.gz
awk -v IDX_COL=8 -v SEP=":" -f- <<'EOF' <(zcat ${fastq_r1}) | gzip > ${fastq_fmt_r1}
NR % 4 == 1 {
  split($1, tmp_array, SEP)
  read_name = tmp_array[1]
  for (i in tmp_array) { if ( i == IDX_COL || i == 1 ) { continue } else { read_name = read_name":"tmp_array[i] } }
  print read_name"/"tmp_array[IDX_COL]
  next
}
{ print }
EOF

fastq_fmt_r2=${out_dir}/${sample_id}.paired_end.R2.fq.gz
awk -v IDX_COL=8 -v SEP=":" -f- <<'EOF' <(zcat ${fastq_r2}) | gzip > ${fastq_fmt_r2}
NR % 4 == 1 {
  split($1, tmp_array, SEP)
  read_name = tmp_array[1]
  for (i in tmp_array) { if ( i == IDX_COL || i == 1 ) { continue } else { read_name = read_name":"tmp_array[i] } }
  print read_name"/"tmp_array[IDX_COL]
  next
}
{ print }
EOF

fastq_clean_r1=${out_dir}/${sample_id}.paired_end.clean.R1.fq.gz
fastq_clean_r2=${out_dir}/${sample_id}.paired_end.clean.R2.fq.gz
cutadapt -j 10 -n 10 -e 0.4 -q 20 -m 22:40 --trim-n \
  -G AAGCAGTGGTATCAACGCAGAGATC -G GGAAGAGCGTCGTGTAGGGAAAGAG \
  -A AAAAAAAAAA -A TTTTTTTTTT \
  -o ${fastq_clean_r1} -p ${fastq_clean_r2} \
  ${fastq_fmt_r1} ${fastq_fmt_r2}

# Mapping reads to genomes by STAR (STAR-solo)
star_genome_dir=~/Documents/projects/resources/references/STAR_GRCm38
genomic_feature_file=${project_dir}/outputs/references/genome/mus_musculus.90.gtf
STAR --runMode alignReads \
  --twopassMode Basic \
  --soloType CB_UMI_Simple --soloCBwhitelist None \
  --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 10 \
  --genomeDir ${star_genome_dir} \
  --sjdbGTFfile ${genomic_feature_file} \
  --readFilesIn ${fastq_clean_r2} ${fastq_clean_r1} \
  --readFilesCommand zcat \
  --runThreadN 10 \
  --alignEndsType EndToEnd \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NM NH HI nM AS CR UR CB UB GX GN sS sQ sM MD \
  --outFileNamePrefix ${out_dir}/${sample_id}.paired_end.
samtools index ${out_dir}/${sample_id}.paired_end.Aligned.sortedByCoord.out.bam

# Mark duplications
raw_bam_file=${out_dir}/${sample_id}.paired_end.Aligned.sortedByCoord.out.bam
mkdup_bam_file=${out_dir}/${sample_id}.paired_end.star_solo.possorted.markdup.bam
samtools markdup -Sc --no-PG --include-fails --barcode-tag sS --threads 10 ${raw_bam_file} ${mkdup_bam_file}
samtools index ${mkdup_bam_file}

# Overview of read alignments. Works for *.Log.final.out and *.Solo.out/Gene/Summary.csv
log_final_out=${out_dir}/${sample_id}.paired_end.Log.final.out
star_solo_summary=${out_dir}/${sample_id}.paired_end.Solo.out/Gene/Summary.csv
total_reads=$(zgrep -c ^@ ${fastq_r1})
awk -v TOTAL_READS=${total_reads} -f- <<'EOF' ${log_final_out} ${star_solo_summary} | csvtk pretty
function get_digit(string, type) {
  if (type == "p") { return match(string, /([0-9.]+%$)/, DIGIT_STR) ? DIGIT_STR[1] : "NA" }
  else if (type == "f") { return match(string, /([0-9.]+$)/, DIGIT_STR) ? DIGIT_STR[1] : "NA" }
  else if (type == "n") { return match(string, /([0-9]+$)/, DIGIT_STR) ? DIGIT_STR[1] : "NA" }
  return "NA"
}

BEGIN {
  kstr = "TotalReads NrIR AvgLenIR NrUniMR PcUniMR RtMmBaseMR NrMulLocMR PcMulLocMR NrTMLociMR PcTMLociMR NrMmUMR"
  kstr = kstr" PcMmUMR NrShortUMR PcShortUMR NrOtherUMR PcOtherUMR NrChReads PcChReads NrBarcodes NrUniInCellMR"
  kstr = kstr" NrMeanUniPB NrMedianUniPB NrUMIsPB NrMeanUMIPB NrMedianUMIPB NrMeanGenePB NrMedianGenePB NrTtlGenes"
  split(kstr, klist, " ")
  RESULT["TotalReads"] = TOTAL_READS
  for (x in klist) { if (x == length(klist)) { print klist[x] } else { printf klist[x]"," } }
}

# Summary from Log.Final.out by STAR
$0 ~ /Number of input reads/ { RESULT["NrIR"] = get_digit($0, "n") } # Input reads
$0 ~ /Average input read length/ { RESULT["AvgLenIR"] = get_digit($0, "f") } # Input Reads
$0 ~ /Uniquely mapped reads number/ { RESULT["NrUniMR"] = get_digit($0, "n") } # Uniquely mapped reads
$0 ~ /Uniquely mapped reads % / { RESULT["PcUniMR"] = get_digit($0, "p") } # Uniquely mapped reads
$0 ~ /Mismatch rate per base, %/ { RESULT["RtMmBaseMR"] = get_digit($0, "p") } # Mismatches of uniquely mapped reads
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
$0 ~ /Number of chimeric reads/ { RESULT["NrChReads"] = get_digit($0, "n") }
$0 ~ /% of chimeric reads/ { RESULT["PcChReads"] = get_digit($0, "p") }

# Summary from Summary.csv by STAR-solo
$0 ~ /Estimated Number of Cells/ { RESULT["NrBarcodes"] = get_digit($0, "n") }
$0 ~ /Unique Reads in Cells Mapped to Gene/ { RESULT["NrUniInCellMR"] = get_digit($0, "n") }
$0 ~ /Mean Reads per Cell/ { RESULT["NrMeanUniPB"] = get_digit($0, "f") }
$0 ~ /Median Reads per Cell/ { RESULT["NrMedianUniPB"] = get_digit($0, "f") }
$0 ~ /UMIs in Cells/ { RESULT["NrUMIsPB"] = get_digit($0, "n") }
$0 ~ /Mean UMI per Cell/ { RESULT["NrMeanUMIPB"] = get_digit($0, "f") }
$0 ~ /Median UMI per Cell/ { RESULT["NrMedianUMIPB"] = get_digit($0, "f") }
$0 ~ /Mean Gene per Cell/ { RESULT["NrMeanGenePB"] = get_digit($0, "f") }
$0 ~ /Median Gene per Cell/ { RESULT["NrMedianGenePB"] = get_digit($0, "f") }
$0 ~ /Total Gene Detected/ { RESULT["NrTtlGenes"] = get_digit($0, "n") }

END { for (x in klist) { if (x == length(klist)) { print RESULT[klist[x]] } else { printf RESULT[klist[x]]"," } } }
EOF


# Quantification of reads
mkdup_bam_file=${out_dir}/${sample_id}.paired_end.star_solo.markdup.bam
feature_counts_tbl=${out_dir}/${sample_id}.paired_end.FeatureCounts.txt
featureCounts -T 12 --extraAttributes gene_name --ignoreDup -a ${genomic_feature_file} -o ${feature_counts_tbl} ${mkdup_bam_file}


# Overview of read alignments, works for multiple *.Log.final.out
log_final_out=${out_dir}/${sample_id}.paired_end.Log.final.out
awk -v TOTAL_READS=${total_reads} -f- <<'EOF' ${log_final_out} | csvtk pretty
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
  for (x in klist) { if (x == length(klist)) { print klist[x] } else { printf klist[x]"," } }
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

ENDFILE { for (x in klist) { if (x == length(klist)) { print RESULT[klist[x]] } else { printf RESULT[klist[x]]"," } } }
EOF

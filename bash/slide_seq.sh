#!/bin/bash
# File: slide_seq.sh
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Aug 11, 2024
# Updated:

project_dir=/home/zzhang/Documents/projects/wp_vasaseq
using_celatlas=true

# Load celatlas_spatial environment
if [[ -e ${project_dir}/scripts/.celatlas_env && $using_celatlas == 'true' ]]; then
  source "${project_dir}/scripts/.celatlas_env/bin/activate"
fi


# Preprocess slide-seq reads. I.e., remove read pairs without UMI and copy UMI to read 2 file which is also the output.
# python -m pdb ${project_dir}/scripts/py3/preprocess_slide_seq_reads.py \
# 1. Do we really need to discard read pairs that read 1 does not have proper UMI/BC information
# 2. The way to remove duplications?
if [[ $using_celatlas == 'true' ]]; then
  celatlas_spatial rna barcode \
    --chemistry BBV2.4 \
    --sample test-A-A \
    --fq1 ${project_dir}/inputs/20240813_encoded_E7.5toE7.75_MouseEmbryo/test-A-AL0805_L04_R1.fq.gz \
    --fq2 ${project_dir}/inputs/20240813_encoded_E7.5toE7.75_MouseEmbryo/test-A-AL0805_L04_R2.fq.gz \
    --outdir ${project_dir}/outputs/analysis/20240813_encoded_E7.5toE7.75_MouseEmbryo/test-A-A \
    --whitelist ${project_dir}/inputs/undecoded_with_polya/barcode.txt \
    --gzip
else
  python ${project_dir}/scripts/py3/preprocess_slide_seq_reads.py \
    -1 ${project_dir}/inputs/undecoded_with_polya/spNEB_1.fastq.gz \
    -2 ${project_dir}/inputs/undecoded_with_polya/spNEB_2.fastq.gz \
    -w ${project_dir}/inputs/undecoded_with_polya/barcode.txt \
    -o ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.with_umi_bc.fastq.gz
fi

# QC and preprocessing using fastp
fastp -w 4 -A \
  -i ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.with_umi_bc.fastq.gz \
  -o ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.with_umi_bc.clean.fastq.gz \
  --html ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.with_umi_bc.html \
  --json ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.with_umi_bc.json \
  --low_complexity_filter --complexity_threshold 20

# Remove poly A
cutadapt -j 4 -n 6 -m 30 -e 2 -a 'T{10}' -a 'A{10}' \
  -o ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.with_umi_bc.clean.rm_polyx.fastq.gz  \
  ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.with_umi_bc.clean.fastq.gz 


# Split the reads into barcode/UMI and molecule reads
awk -f- \
  -v OUT_R1=${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.R1.fastq \
  -v OUT_R2=${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.R2.fastq \
  <<'EOF' <(zcat ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.with_umi_bc.clean.rm_polyx.fastq.gz)
NR % 4 == 1 {
  read_name = $0; split(read_name, a, ":");
  cb = gensub("CB_", "", "g", a[8]); umi = gensub("UMI_", "", "g", a[9]); next
}
NR % 4 == 2 { read_bases = $0; next }
NR % 4 == 3 { next }
NR % 4 == 0 {
  read_qual = $0
  OFS = "\n"
  print read_name >> OUT_R1
  print cb""umi >> OUT_R1
  print "+" >> OUT_R1
  print "JJJJJJJJJJJJJJJJJJJJJJ" >> OUT_R1

  print read_name >> OUT_R2
  print read_bases >> OUT_R2
  print "+" >> OUT_R2
  print read_qual >> OUT_R2
}
EOF

# Compress the UMI/CB fastq and molecule fastq
gzip ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.R{1,2}.fastq

zgrep -Eo ':CB_[ATCG]+' ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.R1.fastq.gz | sort -u | wc -l

# Run STAR
# 2,280,342 mapped reads (read 1)
#   918,069 mapped barcodes
#       310 maximum number of reads per barcode
STAR --runMode alignReads \
  --twopassMode Basic \
  --soloType CB_UMI_Simple --soloCBwhitelist None \
  --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 10 \
  --genomeDir ${project_dir}/outputs/references/genome/STAR_build \
  --sjdbGTFfile ${project_dir}/outputs/references/genome/mus_musculus.90.gtf \
  --readFilesIn ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.R{2,1}.fastq.gz \
  --readFilesCommand zcat \
  --runThreadN 4 \
  --alignIntronMax 50000 \
  --alignEndsType EndToEnd \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
  --outSAMstrandField intronMotif \
  --outFileNamePrefix ${project_dir}/outputs/analysis/undecoded_with_polya/spNEB.

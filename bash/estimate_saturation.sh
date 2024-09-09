#!/bin/bash
# File: estimate_saturation.sh
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Aug 30, 2024
# Updated:

project_dir=~/Documents/projects/wp_vasaseq



fastq_dir=/mnt/WDRed/vasa_seq/inputs/
sauturation_dir=${project_dir}/outputs/analysis/preprocessing/geo_seq/sauturation
for per_batch in 240409_Lib_embryo 240612_Lib_28region 240620_Lib_38region 240703_Lib_32region 240710_Lib_37region 240717_Lib_28region; do
  work_dir=${sauturation_dir}/${per_batch}
  if [[ ! -d ${work_dir} ]]; then echo "Creating ${work_dir}" && mkdir -p ${work_dir}; fi

  in_fastq_r1=${fastq_dir}/${per_batch}/00.mergeRawFq/${per_batch}/${per_batch}_raw_1.fq.gz
  in_fastq_r2=${fastq_dir}/${per_batch}/00.mergeRawFq/${per_batch}/${per_batch}_raw_2.fq.gz
  if [[ ! -f ${in_fastq_r1} || ! -f ${in_fastq_r2} ]]; then echo "Missing ${in_fastq_r1} or ${in_fastq_r2}"; break; fi

  # All read ids
  zgrep '^@' ${in_fastq_r1} | sed 's#/1$##g; s#^@##g' | gzip > ${work_dir}/read_id.total.txt.gz
  total_reads=$(zcat "${work_dir}/read_id.total.txt.gz" | wc -l)

  for per_frac in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do
    n_reads=$(echo "${per_frac} * ${total_reads}" | bc | cut -f1 -d.)

    # Subset reads of a fraction by random selection.
    shuf -n "${n_reads}" "${work_dir}/read_id.total.txt.gz" > "${work_dir}/read_id.${per_frac}.txt"

    for per_bam in "${project_dir}"/outputs/analysis/preprocessing/geo_seq/read_alignments/"${per_batch}"/*/*.merged.bam; do
      if [[ ! -f ${per_bam} ]]; then echo "Missing ${per_bam}"; break; fi
      per_region=$(echo "${per_bam}" | xargs basename | cut -f1 -d.)
      samtools view "${per_bam}" | cut -f 1 -d: | gzip > "${work_dir}/read_id.all_mapped.${per_region}.txt.gz"
      zgrep -wf "${work_dir}/read_id.${per_frac}.txt" "${work_dir}/read_id.all_mapped.${per_region}.txt.gz" > "${work_dir}/read_id.${per_frac}.${per_region}.txt"
    done
  done

  break
done

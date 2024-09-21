#!/bin/bash
# File: estimate_saturation.sh
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Aug 30, 2024
# Updated:

project_dir=~/Documents/projects/wp_vasaseq


fastq_dir=/mnt/WDRed/vasa_seq/inputs/
sauturation_dir=${project_dir}/outputs/analysis/preprocessing/geo_seq/sauturation
# for per_batch in 240612_Lib_28region 240620_Lib_38region 240703_Lib_32region 240710_Lib_37region 240717_Lib_28region; do
for per_batch in 240409_Lib_embryo 240612_Lib_28region 240620_Lib_38region 240703_Lib_32region 240710_Lib_37region 240717_Lib_28region; do
  work_dir=${sauturation_dir}/${per_batch}
  if [[ ! -d ${work_dir} ]]; then echo "Creating ${work_dir}" && mkdir -p ${work_dir}; fi

  mkdir -p ${work_dir}/temp

  in_fastq_r1=${fastq_dir}/${per_batch}/00.mergeRawFq/${per_batch}/${per_batch}_raw_1.fq.gz
  in_fastq_r2=${fastq_dir}/${per_batch}/00.mergeRawFq/${per_batch}/${per_batch}_raw_2.fq.gz
  if [[ ! -f ${in_fastq_r1} || ! -f ${in_fastq_r2} ]]; then echo "Missing ${in_fastq_r1} or ${in_fastq_r2}"; break; fi

  # All read ids
  all_read_ids=${work_dir}/read_id.total.txt.gz 
  if [[ ! -f ${all_read_ids} ]]; then
    pigz -dc -p 2 ${in_fastq_r1} | grep '^@' | sed 's#/1$##g; s#^@##g' | pigz -p 2 > ${all_read_ids}
  fi
  total_reads=$(zcat ${work_dir}/read_id.total.txt.gz | wc -l)

  for per_perc in 10 20 30 40 50 60 70 80 90 100; do
    if [[ -f ${work_dir}/read_counts.${per_perc}.txt ]]; then echo "Skipping ${per_perc}"; continue; fi
    if [[ ${per_perc} -ne 100 ]]; then
      n_reads=$(echo "scale=2; ${per_perc} / 100 * ${total_reads}" | bc | cut -f1 -d.)

      # Subset reads of a fraction by random selection.
      zcat ${work_dir}/read_id.total.txt.gz | head -n ${n_reads} > ${work_dir}/temp/read_id.${per_perc}.txt
      # awk 'NR==FNR{seen[$0]++;next} seen[FNR]==1{print}' <(shuf -i 1-${total_reads} -n ${n_reads}) <() > ${work_dir}/temp/read_id.${per_perc}.txt
      # zcat ${work_dir}/read_id.total.txt.gz | shuf -n ${n_reads} > ${work_dir}/temp/read_id.${per_perc}.txt

      for per_bam in ${project_dir}/outputs/analysis/preprocessing/geo_seq/read_alignments/${per_batch}/*/*.merged.bam; do
        if [[ ! -f ${per_bam} ]]; then echo "Missing ${per_bam}"; break; fi
        per_region=$(echo "${per_bam}" | xargs basename | cut -f1 -d.)

        # mkdir -p ${work_dir}/${per_region}
        # samtools view ${per_bam} | cut -f 1 | split -d -a 4 -l 2500 - ${work_dir}/${per_region}/read_id.mapped.${per_region}.
        samtools view ${per_bam} | cut -f 1 > ${work_dir}/temp/read_id.mapped.${per_region}.txt
        awk 'NR==FNR { split($1,tmp_array,":"); ids[tmp_array[1]]=$1 } NR!=FNR { if ($1 in ids) {print(ids[$1])} else {next} }' \
          ${work_dir}/temp/read_id.mapped.${per_region}.txt ${work_dir}/temp/read_id.${per_perc}.txt \
          | sort -u \
          > ${work_dir}/temp/read_id.${per_perc}.${per_region}.txt

        # Subset BAM file by selected read ids
        samtools view -N ${work_dir}/temp/read_id.${per_perc}.${per_region}.txt -bo ${work_dir}/temp/read_id.${per_perc}.${per_region}.bam ${per_bam}

        # Report
        echo ${per_batch},${per_perc},${per_region},${n_reads},${total_reads} >> ${sauturation_dir}/sauturation_report.txt
      done
    else
      for per_bam in ${project_dir}/outputs/analysis/preprocessing/geo_seq/read_alignments/${per_batch}/*/*.merged.bam; do
        cp ${per_bam} ${work_dir}/temp/
      done
    fi

    # Count mapped reads
    featureCounts -p -t transcript -g transcript_id --extraAttributes gene_id,gene_name,gene_biotype -T 2 \
      -a ${project_dir}/outputs/references/genomic_features/mus_musculus.90.wto_rRNA.transcripts.gtf \
      -o ${work_dir}/read_counts.${per_perc}.txt ${work_dir}/temp/*.bam

    rm -fr ${work_dir}/temp/*
  done

  rm -fr ${all_read_ids} && touch ${all_read_ids}
done

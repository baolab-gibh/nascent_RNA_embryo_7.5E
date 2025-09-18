#!/bin/bash
# File: estimate_gene_body_coverave.sh

proj_dir=~/Documents/projects/wp_vasaseq
in_gtf=/home/zzhang/Documents/projects/wp_vasaseq/outputs/references/genome/mus_musculus.90.gtf

# GEO VASA SLAM seq
alignment_dir=${proj_dir}/outputs/analysis/preprocessing/geo_vasa_slam_seq/read_alignments
statistics_dir=${proj_dir}/outputs/analysis/preprocessing/geo_vasa_slam_seq/alignment_statistics
batch=240717_Lib_28region
for batch in 240409_Lib_embryo 240612_Lib_28region 240620_Lib_38region 240703_Lib_32region 240710_Lib_37region 240717_Lib_28region; do
  for sample in $(ls --color=never ${alignment_dir}/${batch}); do
    in_bam=${alignment_dir}/${batch}/${sample}/${sample}.merged.bam
    out_dir=${statistics_dir}/${batch}/${sample}

    rm -fr ${out_dir} && mkdir -p ${out_dir}
    samtools view -h -F 1 -O SAM -@ 4 ${in_bam} | samtools sort -O BAM -@ 4 -o ${out_dir}/${sample}.single_ended.bam
    java -jar ~/tools/QoRTs/QoRTs.jar QC --addFunctions calcDetailedGeneCounts \
      --maxReadLength 250 --singleEnded --generatePlots ${out_dir}/${sample}.single_ended.bam ${in_gtf} ${out_dir}
    rm -f ${out_dir}/${sample}.single_ended.bam 
  done
done


# 2019 Peng et al, Nature
batch=2019_peng_etal_nature
alignment_dir=${proj_dir}/outputs/analysis/preprocessing/2019_peng_etal_nature/read_alignments
statistics_dir=${proj_dir}/outputs/analysis/preprocessing/2019_peng_etal_nature/mapping_reports
for sample in $(ls --color=never ${alignment_dir}/); do
  in_bam=${alignment_dir}/${sample}/${sample}.rmdup.bam
  out_dir=${statistics_dir}/${sample}

  rm -fr ${out_dir} && mkdir -p ${out_dir}
  java -jar ~/tools/QoRTs/QoRTs.jar QC --addFunctions calcDetailedGeneCounts --maxReadLength 250 --generatePlots ${in_bam} ${in_gtf} ${out_dir}
done

#!/usr/bin/env bash
# File: estimate_mapping_distribution.sh
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Apr 17, 2025

app_image=~/tools/containers/images/vasaslam_seq.sif

project_dir=~/Documents/projects/wp_vasaseq/
working_dir=${project_dir}/outputs/analysis/geo_slam/problems/smart_seq_3prime_utr

bam_file=${project_dir}/outputs/analysis/preprocessing/geo_slam/SLAM_20250321_byLDJ/alignment_persample/lps_6h_r3.star_solo.rmdup.bam
gtf_file=~/Documents/projects/resources/references/gencode.vM29lift38.annotation.gtf
count_detail=${project_dir}/outputs/analysis/preprocessing/geo_slam/SLAM_20250321_byLDJ/nascent_rna/per_sample/lps_6h_r3/lps_6h_r3.mismatchdetails.tsv
nascent_rna_counts=${project_dir}/outputs/analysis/preprocessing/geo_slam/SLAM_20250321_byLDJ/nascent_rna/per_sample/lps_6h_r3/lps_6h_r3.tsv.gz
no4su_control_cit_file=${project_dir}/outputs/analysis/preprocessing/geo_slam/SLAM_20250321_byLDJ/nascent_rna/cit_files/per_sample/

# Genomic regions
utr_regions=${working_dir}/utr_regions.bed
cds_regions=${working_dir}/cds_regions.bed
gene_regions=${working_dir}/gene_regions.bed
exon_regions=${working_dir}/exon_regions.bed
intron_regions=${working_dir}/intron_regions.bed

awk -F$'\t' '$3=="UTR" {OFS="\t"; print $1,$4-1,$5}' ${gtf_file} > ${utr_regions}
awk -F$'\t' '$3=="CDS" {OFS="\t"; print $1,$4-1,$5}' ${gtf_file} > ${cds_regions}
awk -F$'\t' '$3=="gene" {OFS="\t"; print $1,$4-1,$5}' ${gtf_file} > ${gene_regions}
awk -F$'\t' '$3 != "gene" && $3 != "transcript" {OFS="\t"; print $1,$4-1,$5}' ${gtf_file} > ${exon_regions}
apptainer exec ${app_image} bedtools subtract -a ${gene_regions} -b ${exon_regions} > ${intron_regions}

# Read counts
n_ttl_reads=$(apptainer exec ${app_image} samtools view -c ${bam_file})
n_utr_reads=$(apptainer exec ${app_image} samtools view -c --region-file ${utr_regions} ${bam_file})
n_cds_reads=$(apptainer exec ${app_image} samtools view -c --region-file ${cds_regions} ${bam_file})
n_exon_reads=$(apptainer exec ${app_image} samtools view -c --region-file ${exon_regions} ${bam_file})
n_intron_reads=$(apptainer exec ${app_image} samtools view -c --region-file ${intron_regions} ${bam_file})

echo "Ttl: ${n_ttl_reads}" "UTR: ${n_utr_reads}" "CDS: ${n_cds_reads}" "Exon: ${n_exon_reads}" "Intron: ${n_intron_reads}"

for x in utr cds exon intron; do
  echo "[I]: Working on ${x}"
  mkdir -p ${working_dir}/${x}
  apptainer exec ${app_image} samtools view --region-file ${working_dir}/${x}_regions.bed -O BAM -o ${working_dir}/${x}/${x}_reads.bam ${bam_file}
  apptainer exec ${app_image} samtools index ${working_dir}/${x}/${x}_reads.bam
  apptainer exec ${app_image} gedi -e Bam2CIT -p ${working_dir}/${x}/${x}_reads.cit ${working_dir}/${x}/${x}_reads.bam
  apptainer exec ${app_image} gedi -e CorrectCIT -p ${working_dir}/${x}/${x}_reads.cit ${working_dir}/${x}/${x}_reads.corrected.cit
  apptainer exec ${app_image} gedi -e MergeCIT -p ${working_dir}/${x}/${x}_reads.corrected.with_ctrl.cit ${working_dir}/${x}/${x}_reads.corrected.cit ${working_dir}/nc_no4suctrl_r1.corrected.cit
  apptainer exec ${app_image} gedi -e ReadCount -p ${working_dir}/${x}/${x}_reads.corrected.with_ctrl.cit
  apptainer exec ${app_image} gedi -e Slam -genomic mus_musculus.90.wto_rRNA -reads ${working_dir}/${x}/${x}_reads.corrected.with_ctrl.cit \
    -prefix ${working_dir}/${x}/${x}_reads.nascent_rna -introns -minEstimateReads 5000 -trim3p 30 -trim5p 10 -progress -plot -D -full -allGenes \
    -no4sUpattern nc_no4suctrl
  echo "[I]: Finished ${x}"
done

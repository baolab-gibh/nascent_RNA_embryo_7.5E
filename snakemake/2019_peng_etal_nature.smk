#!/usr/bin/env snakemake
# File: vasalam_ppl.smk
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 17, 2024
# Updated: Jul 09, 2024

# TODO:
# 1. For gedi -e Slam, adding 

import sys, csv, pathlib
from snakemake.utils import min_version

# Ensure the snakemake knows how to handle moduliazations.
min_version('6.0')


project_dir = pathlib.Path("~/Documents/projects/wp_vasaseq").expanduser()

# Inputs
genomic_feature_file = project_dir / 'outputs/references/genome/mus_musculus.90.gtf' 
genomic_feature_without_rRNA = project_dir / 'outputs/references/genome/mus_musculus.90.wto_rRNA.gtf'
star_genome_dir = project_dir / 'outputs/references/genome/STAR_build'

sample_list_file = project_dir / 'inputs/2019_peng_etal_nature/selected_samples.txt'
batch = "2019_peng_etal_nature"

# Outputs
sra_dir = project_dir / 'outputs/analysis/preprocessing' / batch / 'sra'
fq_dir = project_dir / 'outputs/analysis/preprocessing' / batch / 'fastqs'
qc_dir = project_dir / 'outputs/analysis/preprocessing' / batch / 'quality_control'
ra_dir = project_dir / 'outputs/analysis/preprocessing' / batch / 'read_alignments'
qm_dir = project_dir / 'outputs/analysis/preprocessing' / batch / 'mapping_reports'

# Created wildcards
with open(sample_list_file, "r") as fh:
  sample_list = [x.strip() for x in fh]
final_outputs = expand(qm_dir / '{sample_id}/{sample_id}.done', sample_id=sample_list)


wildcard_constraints:
  sample_id = r'\w+', group_id = r'\w+'


localrules:
  all,
  s00_index_genome_star,
  s01_download_from_sra,
  s02_dump_fastq,
  s03_quality_control,
  s04_align_reads,
  s05_remove_duplicates,
  s06_estimate_mapping_quality,


rule all:
  input: final_outputs


rule s00_index_genome_star:
  input:
    fasta_file = project_dir / 'outputs/references/genome/mus_musculus.90.fasta',
    genomic_feature_file = project_dir / 'outputs/references/genome/mus_musculus.90.gtf',
  output:
    star_genome_dir = directory(star_genome_dir),
  resources:
    cpus_per_task = 6
  shell:
    '''
    STAR --runMode genomeGenerate --runThreadN {resources.cpus_per_task} --genomeFastaFiles mus_musculus.90.fasta --genomeDir {output.star_genome_dir}
    '''


rule s01_download_from_sra:
  input:
    samples = sample_list_file
  output:
    srr_file = sra_dir / '{sample_id}'
  params:
    out_dir = lambda _, output: Path(output(0)).parent,
  resources:
    cpus_per_task = 4
  shell:
    '''
    mkdir -p {params.out_dir}
    wget -cO {output.srr_file} 'https://sra-pub-run-odp.s3.amazonaws.com/sra/{wildcards.sample_id}/{wildcards.sample_id}'
    '''


rule s02_dump_fastq:
  input:
    srr_file = sra_dir / '{sample_id}',
  output:
    raw_r1 = fq_dir / '{sample_id}/{sample_id}_1.fastq.gz',
    raw_r2 = fq_dir / '{sample_id}/{sample_id}_2.fastq.gz',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  resources:
    cpus_per_task = 4
  shell:
    '''
    mkdir -p {params.out_dir}
    fastq-dump --split-files --gzip -O {params.out_dir} {input.srr_file}
    '''


rule s03_quality_control:
  input:
    raw_r1 = fq_dir / '{sample_id}/{sample_id}_1.fastq.gz',
    raw_r2 = fq_dir / '{sample_id}/{sample_id}_2.fastq.gz',
  output:
    clean_r1 = qc_dir / '{sample_id}/{sample_id}.clean.R1.fq.gz',
    clean_r2 = qc_dir / '{sample_id}/{sample_id}.clean.R2.fq.gz',
    json_report = qc_dir / '{sample_id}/{sample_id}.json',
    html_report = qc_dir / '{sample_id}/{sample_id}.html',
  resources:
    cpus_per_task = 4
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    fastp -w {resources.cpus_per_task} \
        -i {input.raw_r1} -I {input.raw_r2} \
        -o {output.clean_r1} -O {output.clean_r2} \
        --html {output.html_report} --json {output.json_report} \
        --correction --disable_length_filtering \
        --cut_front --cut_tail --cut_mean_quality 30 \
        --trim_poly_g --trim_poly_x --poly_x_min_len 10 \
        --umi --umi_loc read1 --umi_len 8 --umi_delim ' ' \
        --low_complexity_filter --complexity_threshold 50

    fastqc -q -t {resources.cpus_per_task} -o {params.out_dir} {input.raw_r1} {input.raw_r2}
    fastqc -q -t {resources.cpus_per_task} -o {params.out_dir} {output.clean_r1} {output.clean_r2}
    '''


rule s04_align_reads:
  input:
    star_genome_dir = star_genome_dir,
    genomic_feature_file = genomic_feature_file,
    clean_r1 = qc_dir / '{sample_id}/{sample_id}.clean.R1.fq.gz',
    clean_r2 = qc_dir / '{sample_id}/{sample_id}.clean.R2.fq.gz',
  output:
    bam_file = ra_dir / '{sample_id}/{sample_id}.bam',
    bai_file = ra_dir / '{sample_id}/{sample_id}.bam.bai',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    out_prefix = lambda wc: ra_dir / wc.sample_id / (wc.sample_id + '.'),
  resources:
    cpus_per_task = 4
  shell:
    '''
    mkdir -p {params.out_dir}
    echo "Aligning PE ends ..."
    STAR --runMode alignReads \
      --twopassMode Basic \
      --genomeDir {input.star_genome_dir} \
      --sjdbGTFfile {input.genomic_feature_file} \
      --readFilesIn {input.clean_r1} {input.clean_r2} \
      --readFilesCommand zcat \
      --runThreadN {resources.cpus_per_task} \
      --alignIntronMax 50000 \
      --alignEndsType EndToEnd \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterMismatchNmax 10 \
      --outFilterMatchNminOverLread 0.5 \
      --outFilterScoreMinOverLread 0.5 \
      --outSAMattributes All \
      --outSAMstrandField intronMotif \
      --outFileNamePrefix {params.out_prefix}

    mv {params.out_prefix}Aligned.sortedByCoord.out.bam {output.bam_file}
    samtools index -b -@ {resources.cpus_per_task} {output.bam_file}
    '''


rule s05_remove_duplicates:
  input:
    bam_file = ra_dir / '{sample_id}/{sample_id}.bam',
    bai_file = ra_dir / '{sample_id}/{sample_id}.bam.bai',
  output:
    bam_file = ra_dir / '{sample_id}/{sample_id}.rmdup.bam',
    bai_file = ra_dir / '{sample_id}/{sample_id}.rmdup.bam.bai',
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.out_dir}
    umi_tools dedup --stdin={input.bam_file} --stdout={output.bam_file} --paired --ignore-umi
    samtools index {output.bam_file}
    '''


rule s06_estimate_mapping_quality:
  input:
    bam_file = ra_dir / '{sample_id}/{sample_id}.rmdup.bam',
    bai_file = ra_dir / '{sample_id}/{sample_id}.rmdup.bam.bai',
  output:
    qualimap_report = qm_dir / '{sample_id}/{sample_id}.done',
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.out_dir}
    qualimap --java-mem-size=8G rnaseq --sorted -bam {input.bam_file} -gtf {genomic_feature_without_rRNA} -outdir {params.out_dir}
    touch {output.qualimap_report}
    '''

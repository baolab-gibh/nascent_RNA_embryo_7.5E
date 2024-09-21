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

house_keeping_genes = ['ACTIN', 'GAPDH', 'S16', 'U1', "beta-actin"]

# Ensure the snakemake knows how to handle moduliazations.
min_version('6.0')


def obtain_pefqs(wc):
  wc.sample_id


def include_control_samples(wc):
  global ctrl_samples, sample_info_list, batches
  if any([x in ctrl_samples for x in sample_info_list]): return []
  if batches == '240612_Lib_28region': return []
  return expand(ct_dir / '{per_ctrl}/{per_ctrl}.corrected.cit', per_ctrl=ctrl_samples)


project_dir = pathlib.Path("~/Documents/projects/wp_vasaseq").expanduser()

ctrl_samples = ["NC_I9_CAGATC", "NC_I10_AGCACT", "NC_O9_CACGTA", "NC_O10_ATGCTC"]

# Inputs
genomic_feature_file = project_dir / 'outputs/references/genome/mus_musculus.90.gtf' 
genomic_feature_without_rRNA = project_dir / 'outputs/references/genome/mus_musculus.90.wto_rRNA.gtf'

transcript_regions = project_dir / 'outputs/references/genomic_features/mus_musculus.90.wto_rRNA.transcripts.gtf'
nc_transcript_regions = project_dir / 'outputs/references/genomic_features/mus_musculus.90.wto_rRNA.noncoding_transcripts.gtf'
erna_regions_expanded = project_dir / 'outputs/references/genomic_features/mus_musculus.90.wto_rRNA.overlapped_with_fantom5_and_ensembl.gtf'
transcript_with_exon_regions = project_dir / 'outputs/references/genomic_features/mus_musculus.90.wto_rRNA.transcripts_with_exons.gtf'

star_genome_dir = project_dir / 'outputs/references/genome/STAR_build'
gedi_index = project_dir / 'outputs/references/genome/GEDI.wto_rRNA'

batches = "20240813_encoded_with_polya_embryo" # Done

# Outputs
fastq_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/fastq' / batches
decode_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/decode' / batches 
quality_control_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/quality_control' / batches 
read_alignment_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/read_alignments' / batches 
nascent_rna_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/nascent_rna' / batches 
mapping_quality_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/alignment_statistics' / batches 

# Created wildcards
sample_info_list = ["test-A-AL0805", "test-A-BL0805"]
final_outputs = [nascent_rna_dir / (batches + '.nascent_rna.done'), mapping_quality_dir / (batches + '.mapping_quality.done')]


# wildcard_constraints:
#   sample_id = r'\w+', group_id = r'\w+'


localrules:
  all,
  s00_index_genome_grand_slam,
  s00_index_genome_star,
  s01_decode_spatial_barcode,
  s02_check_quality,
  s03_remove_polyx,
  s04_split_into_paired_end,
  s05_align_reads_and_quantify_expression,
  s06_estimate_mapping_quality,
  s07_estimate_mapping_quality_done,
  s06_convert_bam_to_cit,
  s07_merge_cits,
  s08_estimate_nascent_rna,


rule all:
  input: final_outputs


rule s00_index_genome_grand_slam:
  input:
    star_genome_dir = star_genome_dir
  output:
    'gedi_index.done'
  shell:
    '''
    gedi -e IndexGenome -s mus_musculus.90.fasta -g mus_musculus.90.gtf -p -nobowtie -nostar -nokallisto
    touch {output}
    '''


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
    STAR --runMode genomeGenerate \
        --runThreadN {resources.cpus_per_task} \
        --sjdbGTFfile {input.genomic_feature_file} \
        --genomeFastaFiles {input.fasta_file} \
        --genomeDir {output.star_genome_dir}
    '''


# rule s01_decode_spatial_barcode:
#   input:
#     fastq_r1 = fastq_dir / '{sample_id}/{sample_id}.R1.fq.gz',
#     fastq_r2 = fastq_dir / '{sample_id}/{sample_id}.R1.fq.gz',
#     white_list = fastq_dir / 'white_list.txt'
#   output:
#     fastq_r1 = decode_dir / '{sample_id}/{sample_id}.R1.fq.gz',
#     fastq_r2 = decode_dir / '{sample_id}/{sample_id}.R2.fq.gz',
#   params:
#     output_dir = lambda _, output: Path(output[0]).parent,
#   resources:
#     cpus_per_task = 6
#   shell:
#     '''
#     celatlas_spatial rna barcode \
#       --chemistry BBV2.4 \
#       --mode strna \
#       --fq1 {input.fastq_r1} \
#       --fq2 {input.fastq_r2} \
#       --whitelist {input.white_list} \
#       --threads {resources.cpus_per_task} \
#       --outdir {params.output_dir} \
#       --gzip
#     '''


rule s01_decode_spatial_barcode:
  input:
    fastq_r1 = fastq_dir / '{sample_id}/{sample_id}.R1.fq.gz',
    fastq_r2 = fastq_dir / '{sample_id}/{sample_id}.R2.fq.gz',
    barcode_file = project_dir / 'inputs/reference/slide_seq.barcodes.txt'
  output:
    fastq_r2 = decode_dir / '{sample_id}/{sample_id}.with_barcode.fq.gz',
  params:
    output_dir = lambda _, output: Path(output[0]).parent,
  resources:
    cpus_per_task = 6
  shell:
    '''
    python {project_dir}/scripts/py3/preprocess_slide_seq_reads.py \
      -1 {input.fastq_r1} \
      -2 {input.fastq_r2} \
      -w {input.barcode_file} \
      -o {output.fastq_r2}
    '''


rule s02_check_quality:
  input:
    fastq_r2 = decode_dir / '{sample_id}/{sample_id}.with_barcode.fq.gz',
  output:
    json_report = quality_control_dir / '{sample_id}/{sample_id}.with_barcode.json',
    html_report = quality_control_dir / '{sample_id}/{sample_id}.with_barcode.html',
    fastq_r2 = quality_control_dir / '{sample_id}/{sample_id}.with_barcode.clean.fq.gz',
  resources:
    cpus_per_task = 4
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    fastp -w {resources.cpus_per_task} \
        -i {input.fastq_r2} \
        -o {output.fastq_r2} \
        --html {output.html_report} \
        --json {output.json_report} \
        --disable_length_filtering \
        --cut_front --cut_tail --cut_mean_quality 30 \
        --trim_poly_g --trim_poly_x --poly_x_min_len 10 \
        --low_complexity_filter --complexity_threshold 20

    fastqc -q -t {resources.cpus_per_task} -o {params.out_dir} {input.fastq_r2}
    fastqc -q -t {resources.cpus_per_task} -o {params.out_dir} {output.fastq_r2}
    '''


rule s03_remove_polyx:
  input:
    fastq_r2 = quality_control_dir / '{sample_id}/{sample_id}.with_barcode.clean.fq.gz',
  output:
    quality_control_dir / '{sample_id}/{sample_id}.with_barcode.clean.rmpolyx_fastqc.html',
    quality_control_dir / '{sample_id}/{sample_id}.with_barcode.clean.rmpolyx_fastqc.zip',
    fastq_r2 = quality_control_dir / '{sample_id}/{sample_id}.with_barcode.clean.rmpolyx.fq.gz',
  resources:
    cpus_per_task = 2
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    # For merged reads, removing poly-x (A/T) for both end
    cutadapt -j 4 -n 6 -m 30 -e 2 -a 'T{{10}}' -a 'A{{10}}' -o {output.fastq_r2} {input.fastq_r2}
    fastqc -q -t {resources.cpus_per_task} -o {params.out_dir} {output.fastq_r2}
    '''


rule s04_split_into_paired_end:
  input:
    fastq_r2 = quality_control_dir / '{sample_id}/{sample_id}.with_barcode.clean.rmpolyx.fq.gz',
  output:
    fastq_r1 = quality_control_dir / '{sample_id}/{sample_id}.with_barcode.clean.rmpolyx.paired_end.R1.fq.gz',
    fastq_r2 = quality_control_dir / '{sample_id}/{sample_id}.with_barcode.clean.rmpolyx.paired_end.R2.fq.gz',
  params:
    fastq_r1 = lambda _, output: Path(output["fastq_r1"]).with_suffix(""),
    fastq_r2 = lambda _, output: Path(output["fastq_r2"]).with_suffix(""),
  shell:
    '''
    awk -f- -v OUT_R1={params.fastq_r1} -v OUT_R2={params.fastq_r2} <<'EOF' <(zcat {input.fastq_r2})
NR % 4 == 1 {{ read_name=$0; split(read_name, a, ":"); cb=gensub("CB_", "", "g", a[8]); umi=gensub("UMI_", "", "g", a[9]); next }}
NR % 4 == 2 {{ read_bases = $0; next }}
NR % 4 == 3 {{ next }}
NR % 4 == 0 {{
  read_qual = $0
  OFS = "\\n"
  print read_name >> OUT_R1
  print cb""umi >> OUT_R1
  print "+" >> OUT_R1
  print "JJJJJJJJJJJJJJJJJJJJJJ" >> OUT_R1

  print read_name >> OUT_R2
  print read_bases >> OUT_R2
  print "+" >> OUT_R2
  print read_qual >> OUT_R2
}}
EOF
  gzip {params.fastq_r1} {params.fastq_r2}
    '''


rule s05_align_reads_and_quantify_expression:
  input:
    star_genome_dir = star_genome_dir,
    genomic_feature_file = genomic_feature_file,
    fastq_r1 = quality_control_dir / '{sample_id}/{sample_id}.with_barcode.clean.rmpolyx.paired_end.R1.fq.gz',
    fastq_r2 = quality_control_dir / '{sample_id}/{sample_id}.with_barcode.clean.rmpolyx.paired_end.R2.fq.gz',
  output:
    bam_file = read_alignment_dir / '{sample_id}/{sample_id}.star_solo.bam',
    bai_file = read_alignment_dir / '{sample_id}/{sample_id}.star_solo.bam.bai',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    pe_out_prefix = lambda wc: read_alignment_dir / wc.sample_id / (wc.sample_id + '.paired_end.'),
    se_out_prefix = lambda wc: read_alignment_dir / wc.sample_id / (wc.sample_id + '.single_end.'),
  resources:
    cpus_per_task = 4
  shell:
    '''
    mkdir -p {params.out_dir}
    echo "Aligning SE ends ..."
    STAR --runMode alignReads \
      --twopassMode Basic \
      --soloType CB_UMI_Simple --soloCBwhitelist None \
      --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 10 \
      --genomeDir {input.star_genome_dir} \
      --sjdbGTFfile {input.genomic_feature_file} \
      --readFilesIn {input.fastq_r2} {input.fastq_r1}  \
      --readFilesCommand zcat \
      --runThreadN 4 \
      --alignIntronMax 50000 \
      --alignEndsType EndToEnd \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM MD \
      --outSAMstrandField intronMotif \
      --outFileNamePrefix {params.se_out_prefix}

    mv {params.se_out_prefix}Aligned.sortedByCoord.out.bam {output.bam_file}
    rm -fr {params.se_out_prefix}_STAR{{genome,pass1,tmp}}

    echo "Creating index ..."
    samtools index -b -@ {resources.cpus_per_task} {output.bam_file}
    '''


rule s06_estimate_mapping_quality:
  input:
    bam_file = read_alignment_dir / '{sample_id}/{sample_id}.star_solo.bam',
    bai_file = read_alignment_dir / '{sample_id}/{sample_id}.star_solo.bam.bai',
  output:
    qualimap_report = mapping_quality_dir / '{sample_id}/qualimapReport.html'
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.out_dir}
    qualimap --java-mem-size=8G rnaseq --sorted -bam {input.bam_file} -gtf {genomic_feature_without_rRNA} -outdir {params.out_dir}
    '''


rule s07_estimate_mapping_quality_done:
  input:
    expand(mapping_quality_dir / '{sample_id}/qualimapReport.html', sample_id=sample_info_list),
  output:
    mapping_quality_dir / (batches + '.mapping_quality.done')
  shell:
    ''' touch {output} '''


rule s06_convert_bam_to_cit:
  input:
    bam_file = read_alignment_dir / '{sample_id}/{sample_id}.star_solo.bam',
    bai_file = read_alignment_dir / '{sample_id}/{sample_id}.star_solo.bam.bai',
  output:
    cit_file = nascent_rna_dir / '{sample_id}/{sample_id}.cit',
    cit_corrected_file = nascent_rna_dir / '{sample_id}/{sample_id}.corrected.cit',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    gedi -e Bam2CIT -p {output.cit_file} {input.bam_file}
    gedi -e CorrectCIT -p {output.cit_file} {output.cit_corrected_file}
    '''


rule s07_merge_cits:
  input:
    cit_corrected_files = expand(nascent_rna_dir / '{sample_id}/{sample_id}.corrected.cit', sample_id=sample_info_list),
  output:
    merged_cit = nascent_rna_dir / (batches + '.all_samples.cit'),
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    gedi -e MergeCIT -p {output.merged_cit} {input.cit_corrected_files}
    gedi -e ReadCount -p {output.merged_cit} 
    mv {output.merged_cit}.metadata.json {output.merged_cit}.metadata.json.bk
    cat {output.merged_cit}.metadata.json.bk | jq '.[] | {{conditions: [.[] | {{total: .total, name: .name | split("/") | .[-2]}}]}}' > {output.merged_cit}.metadata.json
    '''


rule s08_estimate_nascent_rna:
  input:
    merged_cit = nascent_rna_dir / (batches + '.all_samples.cit'),
  output:
    expand(
      nascent_rna_dir / (batches + '.nascent_rna.{per_ext}'),
      per_ext=[
        'binom.tsv', 'binomEstimated.tsv', 'binomFreq.png', 'binomNewProportion.png', 'binomOverlap.tsv',
        'double.pdf', 'doublehit.tsv', 'exonintron.pdf', 'ext.tsv', 'mismatchdetails.tsv', 'mismatches.pdf',
        'mismatches.tsv', 'mismatchpos.pdf', 'mismatchposzoomed.pdf', 'ntr.png', 'ntrstat.tsv', 'param',
        'rates.tsv', 'runtime', 'single_new.rates.png', 'single_old.rates.png', 'strandness', 'tsv.gz'
      ]),
    job_done = nascent_rna_dir / (batches + '.nascent_rna.done'),
  params:
    out_prefix = nascent_rna_dir / (batches + '.nascent_rna'),
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    gedi -e Slam \
        -genomic mus_musculus.90.wto_rRNA \
        -reads {input.merged_cit} \
        -prefix {params.out_prefix} \
        -minEstimateReads 5000 -trim3p 50 -trim5p 30 \
        -progress -plot -D -full \
        -no4sUpattern NC_

    touch {output.job_done}
    '''

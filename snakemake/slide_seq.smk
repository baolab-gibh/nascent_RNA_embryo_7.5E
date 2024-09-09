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

barch_list = ['E200022670-IAA_mix']
batch_list = [
  '240409_Lib_embryo', '240612_Lib_28region', '240620_Lib_38region', '240703_Lib_32region', '240710_Lib_37region', '240717_Lib_28region'
]
batch_list = [
  '20240813_encoded_AL085', '20240813_encoded_BL085'
]


batches = "240710_Lib_37region" # Done
if batches == '20240813_encoded_AL085':
  mp_r1 = project_dir / 'inputs' / batches / 'test-A-AL0805_L04_R1.fq.gz'
  mp_r1 = project_dir / 'inputs' / batches / 'test-A-AL0805_L04_R2.fq.gz'
  # sample_info_tab = project_dir / 'outputs/analysis/preprocessing/barcodes' / (batches + '.sample_barcode.txt')
elif batches == '20240813_encoded_BL085':
  mp_r1 = project_dir / 'inputs' / batches / 'test-A-BL0805_L04_R1.fq.gz'
  mp_r1 = project_dir / 'inputs' / batches / 'test-A-BL0805_L04_R2.fq.gz'
  # sample_info_tab = project_dir / 'outputs/analysis/preprocessing/barcodes' / (batches + '.sample_barcode.txt')
else:
  raise ValueError(f"Unknown batche {batches}")


# Outputs
dm_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/decode' / batches 
qc_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/quality_control' / batches 
ra_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/read_alignments' / batches 
qt_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/quantification' / batches 
nr_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/nascent_rna' / batches 
mq_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/alignment_statistics' / batches 
vc_dir = project_dir / 'outputs/analysis/preprocessing/slide_seq/velocity_counts' / batches 

# Created wildcards
# sample_info_list = []
# with open(sample_info_tab) as handle:
#   for line in handle:
#     if line.startswith('#'): continue
#     _sample_id, barcod = line.strip().split('\t')
#     sample_info_list.append(_sample_id)

final_outputs = [
  nr_dir / (batches + '.nascent_rna.tsv.gz'),
  qt_dir / 'featureCounts' / 'count_tables' / (batches + '.non_coding_transcripts.counts.txt'),
  qt_dir / 'featureCounts' / 'count_tables' / (batches + '.enhancer_rna.read_counts.txt'),
  qt_dir / 'featureCounts' / 'count_tables' / (batches + '.all_transcripts.read_counts.txt'),
  qt_dir / 'featureCounts' / '10X' / 'genes.tsv',
  qt_dir / 'featureCounts' / '10X' / 'matrix.mtx',
  qt_dir / 'featureCounts' / '10X' / 'barcodes.tsv',
  qt_dir / 'slamseq' / '10X' / 'new' / 'barcodes.tsv',
  qt_dir / 'slamseq' / '10X' / 'new' / 'matrix.mtx',
  qt_dir / 'slamseq' / '10X' / 'new' / 'genes.tsv',
  qt_dir / 'slamseq' / '10X' / 'total' / 'barcodes.tsv',
  qt_dir / 'slamseq' / '10X' / 'total' / 'matrix.mtx',
  qt_dir / 'slamseq' / '10X' / 'total' / 'genes.tsv',
  vc_dir / 'velocyto' / 'one_file_per_cell.velocyto_run.loom',
  mq_dir / 'mapping_quality_summarization.done',
]


wildcard_constraints:
  sample_id = r'\w+', group_id = r'\w+'


localrules:
  all,
  s00_index_genome_grand_slam,
  s00_index_genome_star,
  s01_decode,
  s02_check_quality,
  s03_remove_polyx,
  s04_align_reads,
  s05_remove_duplicates,
  s06_estimate_mapping_quality,
  s07_estimate_mapping_quality_done,
  s06_add_cell_barcode,
  s07_quantify_spliced_reads,
  s06_quantify_transcripts,
  s06_quantify_enhancer_rna,
  s06_quantify_noncoding_transcripts,
  s07_convert_transcripts_readcounts_into_10X,
  s06_convert_bam_to_cit,
  s07_merge_cits,
  s08_estimate_nascent_rna,
  s09_convert_nrna_readcounts_into_10X,


rule all:
  input: final_outputs


rule s00_index_genome_grand_slam:
  input:
    star_genome_dir = star_genome_dir
  output:
    'tmp'
  shell:
    '''
    gedi -e IndexGenome -s mus_musculus.90.fasta -g mus_musculus.90.gtf -p -nobowtie -nostar -nokallisto
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


rule s01_decode:
  input:
    read_1 = mp_r1,
    read_2 = mp_r2,
  output:
    dm_file_list = dm_dir / 'demultiplexed_file_list.txt'
  params:
    output_dir = lambda _, output: Path(output[0]).parent,
  resources:
    cpus_per_task = 6
  shell:
    '''
    celatlas_spatial rna barcode \
      --chemistry BBV2.4 \
      --mode strna \
      --fq1 {input.read_1} \
      --fq2 {input.read_2} \
      --whitelist {input.white_list} \
      --threads {resources.cpus_per_task} \
      --outdir {params.output_dir} \
      --gzip
    '''


rule s02_check_quality:
  input:
    dm_file_list = dm_dir / 'demultiplexed_file_list.txt'
  output:
    qc_dir / '{sample_id}/{sample_id}.clean_merged_fastqc.html',
    qc_dir / '{sample_id}/{sample_id}.clean_merged_fastqc.zip',
    qc_dir / '{sample_id}/{sample_id}.clean_unmerged.R1_fastqc.html',
    qc_dir / '{sample_id}/{sample_id}.clean_unmerged.R1_fastqc.zip',
    qc_dir / '{sample_id}/{sample_id}.clean_unmerged.R2_fastqc.html',
    qc_dir / '{sample_id}/{sample_id}.clean_unmerged.R2_fastqc.zip',
    json_report = qc_dir / '{sample_id}/{sample_id}.json',
    html_report = qc_dir / '{sample_id}/{sample_id}.html',
    clean_merged_fq = qc_dir / '{sample_id}/{sample_id}.clean_merged.fq.gz',
    clean_unmerged_r1 = qc_dir / '{sample_id}/{sample_id}.clean_unmerged.R1.fq.gz',
    clean_unmerged_r2 = qc_dir / '{sample_id}/{sample_id}.clean_unmerged.R2.fq.gz',
  resources:
    cpus_per_task = 4
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    read -r raw_r1 raw_r2 <<<$(grep -w {wildcards.sample_id} {input.dm_file_list})
    fastp -w {resources.cpus_per_task} \
        -i ${{raw_r1}} -I ${{raw_r2}} \
        -o {output.clean_unmerged_r1} -O {output.clean_unmerged_r2} \
        --html {output.html_report} --json {output.json_report} \
        --merge --overlap_diff_percent_limit 50 --merged_out {output.clean_merged_fq} \
        --correction --disable_length_filtering \
        --cut_front --cut_tail --cut_mean_quality 30 \
        --trim_poly_g --trim_poly_x --poly_x_min_len 10 \
        --umi --umi_loc read1 --umi_len 8 --umi_delim ' ' \
        --low_complexity_filter --complexity_threshold 50

    fastqc -q -t {resources.cpus_per_task} -o {params.out_dir} ${{raw_r1}} ${{raw_r2}}
    fastqc -q -t {resources.cpus_per_task} -o {params.out_dir} {output.clean_merged_fq} {output.clean_unmerged_r1} {output.clean_unmerged_r2}
    '''


rule s03_remove_polyx:
  input:
    clean_merged_fq = qc_dir / '{sample_id}/{sample_id}.clean_merged.fq.gz',
    clean_unmerged_r1 = qc_dir / '{sample_id}/{sample_id}.clean_unmerged.R1.fq.gz',
    clean_unmerged_r2 = qc_dir / '{sample_id}/{sample_id}.clean_unmerged.R2.fq.gz',
  output:
    qc_dir / '{sample_id}/{sample_id}.clean_merged.rmpolyx_fastqc.html',
    qc_dir / '{sample_id}/{sample_id}.clean_merged.rmpolyx_fastqc.zip',
    qc_dir / '{sample_id}/{sample_id}.clean_unmerged.rmpolyx.R1_fastqc.html',
    qc_dir / '{sample_id}/{sample_id}.clean_unmerged.rmpolyx.R1_fastqc.zip',
    qc_dir / '{sample_id}/{sample_id}.clean_unmerged.rmpolyx.R2_fastqc.html',
    qc_dir / '{sample_id}/{sample_id}.clean_unmerged.rmpolyx.R2_fastqc.zip',
    clean_merged_fq = qc_dir / '{sample_id}/{sample_id}.clean_merged.rmpolyx.fq.gz',
    clean_unmerged_r1 = qc_dir / '{sample_id}/{sample_id}.clean_unmerged.rmpolyx.R1.fq.gz',
    clean_unmerged_r2 = qc_dir / '{sample_id}/{sample_id}.clean_unmerged.rmpolyx.R2.fq.gz',
  resources:
    cpus_per_task = 2
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    # For merged reads, removing poly-x (A/T) for both end
    zcat {input.clean_merged_fq} \
      | cutadapt -j {resources.cpus_per_task} -n 6 -m 60 -e 2 -g '^T{{2}}' -g '^T{{4}}' -g '^T{{6}}' -g '^T{{8}}' -g '^T{{16}}' -g '^A{{2}}' -g '^A{{4}}' -g '^A{{6}}' -g '^A{{8}}' -g '^A{{16}}' - \
      | cutadapt -j {resources.cpus_per_task} -n 6 -m 60 -e 2 -a 'A{{1}}$' -a 'A{{2}}$' -a 'A{{4}}$' -a 'A{{6}}$' -a 'A{{8}}$' -a 'A{{16}}$' -a 'T{{1}}$' -a 'T{{2}}$' -a 'T{{4}}$' -a 'T{{6}}$' -a 'T{{8}}$' -a 'T{{16}}$' - \
      | sed -e 's#\(/[12]\) \([ATCGNatcgn]\+\) \(merged_[0-9]\+_[0-9]\+\)$#:\\3:\\2\\1#' \
      | seqtk seq -r \
      | gzip > {output.clean_merged_fq}

    # For unmerged reads, removing poly-x (A/T) for both end pair-wisely
    cutadapt -j {resources.cpus_per_task} -n 6 -m 60 -e 2 \
      -g '^T{{4}}' -g '^T{{6}}' -g '^T{{8}}' -g '^T{{16}}' -g '^A{{4}}' -g '^A{{6}}' -g '^A{{8}}' -g '^A{{16}}' \
      -A 'T{{10}}' -A 'A{{10}}' \
      -o {output.clean_unmerged_r1} \
      -p {output.clean_unmerged_r2} \
      {input.clean_unmerged_r1} {input.clean_unmerged_r2}

    zcat {output.clean_unmerged_r1} | sed -e 's#\(/[12]\) \([ATCGNatcgn]\+\)$#:\\2\\1#' | seqtk seq -r | gzip > {output.clean_unmerged_r1}.tmp.gz
    mv -f {output.clean_unmerged_r1}.tmp.gz {output.clean_unmerged_r1}

    zcat {output.clean_unmerged_r2} | sed -e 's#\(/[12]\) \([ATCGNatcgn]\+\)$#:\\2\\1#' | seqtk seq -r | gzip > {output.clean_unmerged_r2}.tmp.gz
    mv -f {output.clean_unmerged_r2}.tmp.gz {output.clean_unmerged_r2}

    fastqc -q -t {resources.cpus_per_task} -o {params.out_dir} {output.clean_merged_fq} {output.clean_unmerged_r1} {output.clean_unmerged_r2}

    for per_file in {input}; do
      rm -f ${{per_file}} && touch ${{per_file}}
    done
    touch {output}
    '''


rule s04_align_reads:
  input:
    star_genome_dir = star_genome_dir,
    genomic_feature_file = genomic_feature_file,
    clean_merged_fq = qc_dir / '{sample_id}/{sample_id}.clean_merged.rmpolyx.fq.gz',
    clean_unmerged_r1 = qc_dir / '{sample_id}/{sample_id}.clean_unmerged.rmpolyx.R1.fq.gz',
    clean_unmerged_r2 = qc_dir / '{sample_id}/{sample_id}.clean_unmerged.rmpolyx.R2.fq.gz',
  output:
    bam_file = ra_dir / '{sample_id}/{sample_id}.merged.bam',
    bai_file = ra_dir / '{sample_id}/{sample_id}.merged.bam.bai',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    pe_out_prefix = lambda wc: ra_dir / wc.sample_id / (wc.sample_id + '.paired_end.'),
    se_out_prefix = lambda wc: ra_dir / wc.sample_id / (wc.sample_id + '.single_end.'),
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
      --readFilesIn {input.clean_unmerged_r1} {input.clean_unmerged_r2} \
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
      --outFileNamePrefix {params.pe_out_prefix}

    echo "Aligning SE ends ..."
    STAR --runMode alignReads \
      --twopassMode Basic \
      --genomeDir {input.star_genome_dir} \
      --sjdbGTFfile {input.genomic_feature_file} \
      --readFilesIn {input.clean_merged_fq} \
      --readFilesCommand zcat \
      --runThreadN {resources.cpus_per_task} \
      --alignIntronMax 50000 \
      --alignEndsType EndToEnd \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterMismatchNmax 10 \
      --outFilterScoreMinOverLread 0.5 \
      --outFilterMatchNminOverLread 0.5 \
      --outSAMattributes All \
      --outSAMstrandField intronMotif \
      --outFileNamePrefix {params.se_out_prefix}

    echo "Merging PE and SE alignments ..."
    samtools merge --no-PG \
      -o {output.bam_file} \
      {params.pe_out_prefix}Aligned.sortedByCoord.out.bam {params.se_out_prefix}Aligned.sortedByCoord.out.bam 

    echo "Creating index ..."
    samtools index -b -@ {resources.cpus_per_task} {output.bam_file}

    for per_file in {input.clean_merged_fq} {input.clean_unmerged_r1} {input.clean_unmerged_r2}; do
      rm -f ${{per_file}} && touch ${{per_file}}
    done
    touch {output}
    '''


rule s05_remove_duplicates:
  input:
    bam_file = ra_dir / '{sample_id}/{sample_id}.merged.bam',
    bai_file = ra_dir / '{sample_id}/{sample_id}.merged.bam.bai',
  output:
    bam_file = ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam',
    bai_file = ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam.bai',
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.out_dir}
    umi_tools dedup --stdin={input.bam_file} --stdout={output.bam_file} --umi-separator=':'
    samtools index {output.bam_file}
    '''


rule s06_estimate_mapping_quality:
  input:
    bam_file = ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam',
    bai_file = ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam.bai',
  output:
    qualimap_report = mq_dir / '{sample_id}/qualimapReport.html'
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.out_dir}
    qualimap --java-mem-size=8G rnaseq \
      --sorted \
      -bam {input.bam_file} \
      -gtf {genomic_feature_without_rRNA} \
      -outdir {params.out_dir}
    '''


rule s07_estimate_mapping_quality_done:
  input:
    expand(mq_dir / '{sample_id}/qualimapReport.html', sample_id=sample_info_list),
  output:
    mq_dir / 'mapping_quality_summarization.done'
  shell:
    ''' touch {output} '''


rule s06_quantify_transcripts:
  input:
    transcript_regions = transcript_regions,
    bam_file = expand(ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam', sample_id=sample_info_list),
    bai_file = expand(ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam.bai', sample_id=sample_info_list),
  output:
    count_tab = qt_dir / 'featureCounts/count_tables' / (batches + '.all_transcripts.read_counts.txt'),
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  resources:
    cpus_per_task = 6
  shell:
    '''
    mkdir -p {params.out_dir}
    featureCounts -p -t transcript -g transcript_id --extraAttributes gene_id,gene_name -T {resources.cpus_per_task} \
        -a {input.transcript_regions} -o {output.count_tab} {input.bam_file}
    '''


rule s07_convert_transcripts_readcounts_into_10X:
  input:
    count_tab = qt_dir / 'featureCounts/count_tables' / (batches + '.all_transcripts.read_counts.txt'),
  output:
    genes_tsv = qt_dir / 'featureCounts/10X/genes.tsv',
    matrix_mtx = qt_dir / 'featureCounts/10X/matrix.mtx',
    barcodes_tsv = qt_dir / 'featureCounts/10X/barcodes.tsv',
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  resources:
    cpus_per_task = 6
  shell:
    '''
    mkdir -p {params.out_dir}

    grep -m 1 -e '^Geneid' {input.count_tab} | cut -f9- | tr '\\t' '\\n' | xargs -I% bash -c 'basename $(dirname %)' > {output.barcodes_tsv}
    grep -v -e '^# Program' -e '^Geneid' {input.count_tab} | awk -F$'\\t' '{{OFS="\\t"; print $8, $8}}' > {output.genes_tsv}
    n_ft=$(wc -l {output.genes_tsv} | cut -f1 -d' ')
    n_bc=$(wc -l {output.barcodes_tsv} | cut -f1 -d' ')
    awk -v N_FT=${{n_ft}} -v N_BC=${{n_bc}} -f- <<'EOF' {input.count_tab} > {output.matrix_mtx}
BEGIN  {{ OFS = "\\t"; print "%%MatrixMarket matrix coordinate integer general"; print "%"; print N_FT,N_BC,N_FT * N_BC }}
NR > 2 {{ OFS = "\\t"; i = 9; while (i <= NF) {{ print NR-2,i-8,$i; i++ }} }}
EOF
    '''


rule s06_quantify_noncoding_transcripts:
  input:
    nc_transcript_regions = nc_transcript_regions,
    bam_file = expand(ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam', sample_id=sample_info_list),
    bai_file = expand(ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam.bai', sample_id=sample_info_list),
  output:
    count_tab = qt_dir / 'featureCounts/count_tables' / (batches + '.non_coding_transcripts.counts.txt'),
  resources:
    cpus_per_task = 6
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.out_dir}
    featureCounts -p -t transcript -g transcript_id --extraAttributes gene_id,gene_name -T {resources.cpus_per_task} \
        -a {input.nc_transcript_regions} -o {output.count_tab} {input.bam_file}
    '''


rule s06_quantify_enhancer_rna:
  input:
    erna_regions_expanded = erna_regions_expanded,
    bam_file = expand(ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam', sample_id=sample_info_list),
    bai_file = expand(ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam.bai', sample_id=sample_info_list),
  output:
    count_tab = qt_dir / 'featureCounts/count_tables' / (batches + '.enhancer_rna.read_counts.txt'),
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  resources:
    cpus_per_task = 6
  shell:
    '''
    mkdir -p {params.out_dir}
    featureCounts -p -t enhancer -g enhancer_region -T {resources.cpus_per_task} \
        -a {input.erna_regions_expanded} -o {output.count_tab} {input.bam_file}
    '''


rule s06_add_cell_barcode:
  input:
    bam_file = ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam',
    bai_file = ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam.bai',
  output:
    bam_file = ra_dir / '{sample_id}/{sample_id}.merged.rmdup.add_cb_ub.bam',
    bai_file = ra_dir / '{sample_id}/{sample_id}.merged.rmdup.add_cb_ub.bam.bai',
  params:
    bc = lambda wc: wc.sample_id
  resources:
    cpus_per_task = 2
  shell:
    '''
    samtools view -h -@ {resources.cpus_per_task} {input.bam_file} \
        | awk '$0~/^@/{{print;next}} {{split($1,NM,":");UMI=NM[length(NM)];print $0"\\tCB:Z:"{params.bc}"\\tUB:Z:"UMI}}' \
        | samtools view -@ {resources.cpus_per_task} -O BAM -o {output.bam_file}
    samtools index -@ {resources.cpus_per_task} {output.bam_file}
    '''


rule s07_quantify_spliced_reads:
  input:
    bam_files = expand(ra_dir / '{sample_id}/{sample_id}.merged.rmdup.add_cb_ub.bam', sample_id=sample_info_list),
    bai_files = expand(ra_dir / '{sample_id}/{sample_id}.merged.rmdup.add_cb_ub.bam.bai', sample_id=sample_info_list),
    gtf_file = transcript_with_exon_regions,
  output:
    out_loom = vc_dir / 'velocyto' / 'one_file_per_cell.velocyto_run.loom',
  params:
    out_dir = directory(vc_dir / 'velocyto'),
  shell:
    '''
    mkdir -p {params.out_dir}
    velocyto run -c -o {params.out_dir} {input.bam_files} {input.gtf_file}
    ls {params.out_dir}/*.loom | head -1 | xargs -I% mv % {output.out_loom}
    '''


rule s06_convert_bam_to_cit:
  input:
    bam_file = ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam',
    bai_file = ra_dir / '{sample_id}/{sample_id}.merged.rmdup.bam.bai',
  output:
    cit_file = nr_dir / '{sample_id}/{sample_id}.cit',
    cit_corrected_file = nr_dir / '{sample_id}/{sample_id}.corrected.cit',
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
    cit_corrected_files = expand(nr_dir / '{sample_id}/{sample_id}.corrected.cit', sample_id=sample_info_list),
    control_cit = include_control_samples,
  output:
    merged_cit = nr_dir / (batches + '.all_samples.cit'),
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    gedi -e MergeCIT -p {output.merged_cit} {input.cit_corrected_files} {input.control_cit}
    gedi -e ReadCount -p {output.merged_cit} 
    mv {output.merged_cit}.metadata.json {output.merged_cit}.metadata.json.bk
    cat {output.merged_cit}.metadata.json.bk | jq '.[] | {{conditions: [.[] | {{total: .total, name: .name | split("/") | .[-2]}}]}}' > {output.merged_cit}.metadata.json
    '''


rule s08_estimate_nascent_rna:
  input:
    merged_cit = nr_dir / (batches + '.all_samples.cit'),
  output:
    nr_dir / (batches + '.nascent_rna.binom.tsv'),
    nr_dir / (batches + '.nascent_rna.binomEstimated.tsv'),
    nr_dir / (batches + '.nascent_rna.binomFreq.png'),
    nr_dir / (batches + '.nascent_rna.binomNewProportion.png'),
    nr_dir / (batches + '.nascent_rna.binomOverlap.tsv'),
    nr_dir / (batches + '.nascent_rna.double.pdf'),
    nr_dir / (batches + '.nascent_rna.doublehit.tsv'),
    nr_dir / (batches + '.nascent_rna.exonintron.pdf'),
    nr_dir / (batches + '.nascent_rna.ext.tsv'),
    nr_dir / (batches + '.nascent_rna.mismatchdetails.tsv'),
    nr_dir / (batches + '.nascent_rna.mismatches.pdf'),
    nr_dir / (batches + '.nascent_rna.mismatches.tsv'),
    nr_dir / (batches + '.nascent_rna.mismatchpos.pdf'),
    nr_dir / (batches + '.nascent_rna.mismatchposzoomed.pdf'),
    nr_dir / (batches + '.nascent_rna.ntr.png'),
    nr_dir / (batches + '.nascent_rna.ntrstat.tsv'),
    nr_dir / (batches + '.nascent_rna.param'),
    nr_dir / (batches + '.nascent_rna.rates.tsv'),
    nr_dir / (batches + '.nascent_rna.runtime'),
    nr_dir / (batches + '.nascent_rna.single_new.rates.png'),
    nr_dir / (batches + '.nascent_rna.single_old.rates.png'),
    nr_dir / (batches + '.nascent_rna.strandness'),
    nr_dir / (batches + '.nascent_rna.tsv.gz'),
  params:
    out_prefix = nr_dir / (batches + '.nascent_rna'),
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    gedi -e Slam \
        -genomic mus_musculus.90.wto_rRNA \
        -reads {input.merged_cit} \
        -prefix {params.out_prefix} \
        -minEstimateReads 5000 -trim3p 30 -trim5p 10 \
        -progress -plot -D -full \
        -no4sUpattern NC_
    '''


rule s09_convert_nrna_readcounts_into_10X:
  input:
    nr_dir / (batches + '.nascent_rna.tsv.gz'),
  output:
    new_barcodes = qt_dir / "slamseq/10X/new/barcodes.tsv",
    new_matrix = qt_dir / "slamseq/10X/new/matrix.mtx",
    new_genes = qt_dir / "slamseq/10X/new/genes.tsv",
    total_barcodes = qt_dir / "slamseq/10X/total/barcodes.tsv",
    total_matrix = qt_dir / "slamseq/10X/total/matrix.mtx",
    total_genes = qt_dir / "slamseq/10X/total/genes.tsv",
  params:
    total_count_dir = qt_dir / "slamseq/10X/total",
    new_count_dir = qt_dir / "slamseq/10X/new",
  shell:
    '''
    mkdir -p {params.total_count_dir} {params.new_count_dir}; rm -f {output}

    zgrep -m1 Readcount {input} | tr "\\t" "\\n" | grep -E "Readcount$" | sed 's/ Readcount$//g' | tee {output.total_barcodes} > {output.new_barcodes}
    zcat {input} | cut -f2 -d$'\\t' | sed -n '2,$p' | xargs -I% echo -e %"\\t"% | tee {output.total_genes} > {output.new_genes}

    n_ft=$(wc -l {output.total_genes} | cut -f1 -d' ')
    n_bc=$(wc -l {output.total_barcodes} | cut -f1 -d' ')

    awk -v N_FT=${{n_ft}} -v N_BC=${{n_bc}} -F$'\\t' -f- <<'EOF' <(zcat {input})
function ceil(x, y) {{ y = int(x); return (x > y ? y + 1 : y)}}
BEGIN  {{
  OFS = "\\t"
  print "%%MatrixMarket matrix coordinate integer general" >> "{output.total_matrix}"
  print "%" >> "{output.total_matrix}"
  print N_FT,N_BC,N_FT * N_BC >> "{output.total_matrix}"

  print "%%MatrixMarket matrix coordinate integer general" >> "{output.new_matrix}"
  print "%" >> "{output.new_matrix}"
  print N_FT,N_BC,N_FT * N_BC >> "{output.new_matrix}"

  close("{output.total_matrix}")
  close("{output.new_matrix}")
}}

NR == 1 {{
  ii = 1; iii = 1
  while (ii <= NF) {{
    if ($ii ~ / Readcount/) {{
      sample_id = gensub(" Readcount$", "", "g", $ii); readcount_idx[sample_id] = ii
      sample_list[iii] = sample_id; iii++
    }} else if ($ii ~ / Conversions/) {{
      sample_id = gensub(" Conversions$", "", "g", $ii); conversion_idx[sample_id] = ii
    }} else if ($ii ~ / MAP$/) {{
      sample_id = gensub(" MAP$", "", "g", $ii); new2total_idx[sample_id] = ii
    }}
    ii++
  }}
}}

NR >= 2 {{
  OFS = "\\t"
  for (iiii in sample_list) {{
    sample_id = sample_list[iiii]
    new_to_total_ratio = $new2total_idx[sample_id]
    total_count = ceil($readcount_idx[sample_id])
    new_count = ceil(total_count * new_to_total_ratio)
    print NR - 1, iiii, total_count >> "{output.total_matrix}"
    print NR - 1, iiii, new_count >> "{output.new_matrix}"
  }}

  close("{output.total_matrix}"); close("{output.new_matrix}")
}}
EOF
    '''

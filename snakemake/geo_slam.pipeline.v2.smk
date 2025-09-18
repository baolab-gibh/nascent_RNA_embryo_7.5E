#!/usr/bin/env snakemake
# File: geo_slam.pipeline.smk
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 17, 2024
# Updated: Apr 02, 2025


import os
import sys
import csv
import yaml
from pathlib import Path
from snakemake.utils import min_version
from snakemake.io import Wildcards

min_version('6.0') # Ensure the snakemake knows how to handle moduliazations.

TRANS_DICT = str.maketrans("ATGCN", "TACGN")


# Adapters
truseq_index_5_prime = "AATGATACGGCGACCACCGAGATCTACAC" # R1
truseq_index_3_prime = "ATCTCGTATGCCGTCTTCTGCTTG" # R1
transposase_adapter_i5 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
transposase_adapter_i7 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC" # R1

truseq_index_5_prime_revcmp = "GTGTAGATCTCGGTGGTCGCCGTATCATT"
truseq_index_3_prime_revcmp = "CAAGCAGAAGACGGCATACGAGAT"
transposase_adapter_i5_revcmp = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"
transposase_adapter_i7_revcmp = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"


# # read 1 adapters
# transposase_adapter_i7, truseq_index_5_prime
# transposase_adapter_i7, truseq_index_3_prime
# -a transposase_adapter_i7
# 
# # read 2 adapters
# transposase_adapter_i5_revcmp, truseq_index_5_prime_revcmp
# transposase_adapter_i5_revcmp, truseq_index_3_prime_revcmp
# -A truseq_index_5_prime_revcmp

LIBRARY_COMPONENTS = "".join([
    truseq_index_5_prime, "(?P<i5_index>[ACGT]{8})", transposase_adapter_i5, "(?P<cdna_sequence>[ACGT]+)",
  transposase_adapter_i7, "(?P<i7_index>[ACGT]{8})",   truseq_index_3_prime,
])


class RuleIO():
  """A hard coded class to determine all outputs and specific inputs for some rule.

  A initiation will take a working directory and a list of rule to determine.
  """
  def __init__(self, wk_dir: Path, rule: list[str] | str, info_db: dict | None = None):
    self._rule = rule
    self._wk_dir = wk_dir
    self._info_db = info_db
  def __call__(self, wc: Wildcards | None = None):
    if wc is None and isinstance(self._rule, list): # Determine required outputs.
      outs = []
      for per_batch, batch_info_dict in self._info_db.items():
        for per_sample, sample_info_dict in batch_info_dict["samples"].items():
          if per_sample in batch_info_dict["parameters"]["control_samples"]: continue
          if "quality_control" in self._rule: # s02_quality_control
            outs.append(self._wk_dir / f'{per_batch}/{per_sample}/quality_control/{per_sample}.clean.R1.fastq.gz')
            outs.append(self._wk_dir / f'{per_batch}/{per_sample}/quality_control/{per_sample}.clean.R2.fastq.gz')
          if "alignment" in self._rule: # s03_align_reads
            outs.append(self._wk_dir / f'{per_batch}/{per_sample}/alignment/{per_sample}.bam')
            outs.append(self._wk_dir / f'{per_batch}/{per_sample}/alignment/{per_sample}.bam.bai')
            outs.append(self._wk_dir / f'{per_batch}/{per_sample}/alignment/{per_sample}.final.Log.out')
          if "alignment_overview" in self._rule: # s04_overview_of_alignment
            outs.append(self._wk_dir / f'{per_batch}/{per_sample}/alignment/{per_sample}.overview_of_alignment.csv')
          if "alignment_reports" in self._rule: # s04_estimate_mapping_quality
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/mapping_reports/qualimap/qualimap.done')
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/mapping_reports/qorts/qorts.done')
          if "quantify_enhancer_rna" in self._rule: # s07_quantify_enhancer_rna
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/quantification/{per_sample}.enhancer_rna.read_counts.txt')
          if "quantify_genes" in self._rule: # s07_quantify_genes
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/quantification/{per_sample}.genes.read_counts.txt')
          if "nascent_rna_per_treatment" in self._rule: # s09_estimate_nascent_rna_per_treatment
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/nascent_rna/{per_sample}.tsv.gz')
      return sorted(outs)
    elif isinstance(wc, Wildcards) and isinstance(self._rule, str): # Determine inputs based on wildcards and rule name.
      batch_id, sample_id = getattr(wc, "batch_id", None), getattr(wc, "sample_id", None)
      if self._rule in ["s01_normalize_fastq_name"]:
        return dict(
          fastq_r1 = self._info_db.get(batch_id).get("samples").get(sample_id).get("fastqs").get("fastq_r1"),
          fastq_r2 = self._info_db.get(batch_id).get("samples").get(sample_id).get("fastqs").get("fastq_r2"),
        )
      if self._rule in ["s02_quality_control"]:
        return dict(
          fastq_r1 = self._wk_dir / f'{batch_id}/{sample_id}/fastq/{sample_id}.raw.R1.fastq.gz',
          fastq_r2 = self._wk_dir / f'{batch_id}/{sample_id}/fastq/{sample_id}.raw.R2.fastq.gz',
        )
      if self._rule == "s03_align_reads":
        return dict(
          fastq_r1 = self._wk_dir / f'{batch_id}/{sample_id}/quality_control/{sample_id}.clean.R1.fastq.gz',
          fastq_r2 = self._wk_dir / f'{batch_id}/{sample_id}/quality_control/{sample_id}.clean.R2.fastq.gz',
        )
      if self._rule in ["s04_estimate_mapping_quality"]:
        return dict(
          bam_file = self._wk_dir / f'{batch_id}/{sample_id}/alignment/{sample_id}.bam',
          bai_file = self._wk_dir / f'{batch_id}/{sample_id}/alignment/{sample_id}.bam.bai',
        )
      if self._rule == "s04_overview_of_alignment":
        return dict(
          fastq_raw_r1 = self._wk_dir / f'{batch_id}/{sample_id}/fastq/{sample_id}.raw.R1.fastq.gz',
          fastq_clean_r1 = self._wk_dir / f'{batch_id}/{sample_id}/quality_control/{sample_id}.clean.R1.fastq.gz',
          log_final_out = self._wk_dir / f'{batch_id}/{sample_id}/alignment/{sample_id}.final.Log.out',
        )
      if self._rule in ["s05_remove_duplicates"]:
        return dict(
          bam_file = self._wk_dir / f'{batch_id}/{sample_id}/alignment/{sample_id}.bam',
          bai_file = self._wk_dir / f'{batch_id}/{sample_id}/alignment/{sample_id}.bam.bai',
        )
      if self._rule in ["s06_correct_base_qualities"]:
        return dict(
          bam_file = self._wk_dir / f'{batch_id}/{sample_id}/alignment/{sample_id}.rmdup.bam',
          bai_file = self._wk_dir / f'{batch_id}/{sample_id}/alignment/{sample_id}.rmdup.bam.bai',
        )
      if self._rule in ["s07_quantify_enhancer_rna", "s07_quantify_genes", "s07_convert_bam_to_cit"]:
        return dict(
          bam_file = self._wk_dir / f'{batch_id}/{sample_id}/alignment/{sample_id}.rmdup.bam',
          bai_file = self._wk_dir / f'{batch_id}/{sample_id}/alignment/{sample_id}.rmdup.bam.bai'
        )
      if self._rule in ["s08_merge_cits_per_treatment"]:
        ctrl_cits = []
        for per_batch, sample_info_dict in self._info_db.items():
          if per_batch != batch_id: continue
          for ctrl_sample in sample_info_dict.get("parameters", {}).get("control_samples", []):
            ctrl_cits.append(self._wk_dir / f'{batch_id}/{ctrl_sample}/cit_files/{ctrl_sample}.corrected.cit')
        return dict(
          ctrl_cits = ctrl_cits,
          sample_cit = self._wk_dir / f'{batch_id}/{sample_id}/cit_files/{sample_id}.corrected.cit'
        )
      if self._rule in ["s09_estimate_nascent_rna_per_treatment"]:
        return dict(
          merged_cit = self._wk_dir / f'{batch_id}/{sample_id}/cit_files/{sample_id}.corrected.with_ctrl.cit'
        )
    else:
      raise ValueError(f"Not supported with wc: {type(wc)} and `self._rule`: {type(self._rule)}")


def load_configs(fpath: Path):
  """Load sample information from given input file."""
  with open(fpath, "r") as fhand:
    db = yaml.load(fhand, Loader=yaml.FullLoader)
    if len(db) == 0:
      raise ValueError(f"No sample information found in {fpath}")
  return db


def get_parameters(batch, which):
  """Get parameters from the config file."""
  global configurations
  meta_params = configurations.get("parameters", {})
  batch_info = configurations.get("sequencing_batches", {})
  batch_params = batch_info.get(batch, {}).get("parameters", {}).get(which, None)
  return batch_params if batch_params else meta_params.get(which, "")


#
## Inputs
#
config_file = Path("~/Documents/projects/wp_vasaseq/scripts/snakemake/geo_slam.configuration.v2.yaml").expanduser()
configurations = load_configs(config_file)
parameters = configurations.get("parameters")
sequencing_info = configurations.get("sequencing_batches")

working_dir = Path(configurations["working_dir"]).expanduser()
app_image = Path(configurations["app_image"]).expanduser()

# Inputs
genome_sequence_file = Path(parameters.get("genome_sequence_file")).expanduser()
genomic_feature_file = Path(parameters.get("genomic_feature_file")).expanduser()
erna_feature_file = Path(parameters.get("erna_feature_file")).expanduser()

star_index_dir = Path(parameters.get("star_index_dir")).expanduser()
gedi_index_dir = Path(parameters.get("gedi_index_dir")).expanduser()

required_outputs = parameters.get("required_outputs", [])
rule_io = RuleIO(working_dir, required_outputs, sequencing_info)
final_outputs = rule_io()


#
## Settings
#
wildcard_constraints:
  batch_id = r'\w+', sample_id = r'\w+'


#
## Pipeline rules
#
rule all:
  input: final_outputs


rule s00_index_genome_grand_slam:
  input:
    genome_sequence_file = genome_sequence_file, genomic_feature_file = genomic_feature_file,
  output:
    gedi_index_dir = directory(gedi_index_dir),
  shell:
    '''
    apptainer exec -B /mnt/backup {app_image} gedi \
      -e IndexGenome -s mus_musculus.90.fasta -g mus_musculus.90.gtf -p -nobowtie -nostar -nokallisto
    '''


rule s00_index_genome_star:
  input:
    genome_sequence_file = genome_sequence_file, genomic_feature_file = genomic_feature_file,
  output:
    star_index_dir = directory(star_index_dir),
  resources:
    n_cpus = 10
  shell:
    '''
    apptainer exec -B /mnt/backup {app_image} STAR \
      --runMode genomeGenerate \
      --runThreadN {resources.n_cpus} \
      --sjdbGTFfile {input.genomic_feature_file} \
      --genomeFastaFiles {input.genome_sequence_file} \
      --genomeDir {output.star_index_dir}
    '''


rule s01_normalize_fastq_name:
  input:
    unpack(RuleIO(working_dir, "s01_normalize_fastq_name", sequencing_info))
  output:
    fastq_r1 = working_dir / '{batch_id}/{sample_id}/fastq/{sample_id}.raw.R1.fastq.gz',
    fastq_r2 = working_dir / '{batch_id}/{sample_id}/fastq/{sample_id}.raw.R2.fastq.gz',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  resources:
    n_cpus = 10
  shell:
    '''
    mkdir -p {params.out_dir}

    pigz -dc -p {resources.n_cpus} {input.fastq_r1} | pigz -p {resources.n_cpus} > {output.fastq_r1}
    pigz -dc -p {resources.n_cpus} {input.fastq_r2} | pigz -p {resources.n_cpus} > {output.fastq_r2}
    apptainer exec -B /mnt/backup {app_image} fastqc \
      -q -t {resources.n_cpus} -o {params.out_dir} {output.fastq_r1} {output.fastq_r2}
    '''


rule s02_quality_control:
  input:
    unpack(RuleIO(working_dir, "s02_quality_control", sequencing_info))
  output:
    fastq_r1 = working_dir / '{batch_id}/{sample_id}/quality_control/{sample_id}.clean.R1.fastq.gz',
    fastq_r2 = working_dir / '{batch_id}/{sample_id}/quality_control/{sample_id}.clean.R2.fastq.gz',
  resources:
    n_cpus = 5
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    adapter_sequence = lambda wc: get_parameters(wc.batch_id, "adapter_sequence")
  shell:
    '''
    mkdir -p {params.out_dir}

    # Removing adapter, polyX. Trimming low quality bases. Removing 'N'. Ensure minimum read length (80).
    apptainer exec -B /mnt/backup {app_image} cutadapt \
      --quiet -j {resources.n_cpus} -n 10 -m 100 -q 0 -Q 25 -e 1 --max-n 0 --trim-n \
      -a {transposase_adapter_i7} -A {truseq_index_5_prime_revcmp} \
      -o {output.fastq_r1} -p {output.fastq_r2} \
      {input.fastq_r1} {input.fastq_r2}

    # Check quality
    apptainer exec -B /mnt/backup {app_image} fastqc -q -t {resources.n_cpus} -o {params.out_dir} {output.fastq_r1} {output.fastq_r2}

    # Reverse complement the R2 reads.
    apptainer exec -B /mnt/backup {app_image} seqtk seq -r {output.fastq_r2} \
      | apptainer exec -B /mnt/backup {app_image} gzip > {params.out_dir}/reverse_complemented.R2.fastq.gz
    mv -f {params.out_dir}/reverse_complemented.R2.fastq.gz {output.fastq_r2}
    '''


rule s03_align_reads:
  input:
    unpack(RuleIO(working_dir, "s03_align_reads", sequencing_info)), star_index_dir = star_index_dir,
  output:
    bam_file = working_dir / '{batch_id}/{sample_id}/alignment/{sample_id}.bam',
    bai_file = working_dir / '{batch_id}/{sample_id}/alignment/{sample_id}.bam.bai',
    log_final_out = working_dir / '{batch_id}/{sample_id}/alignment/{sample_id}.final.Log.out',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    out_prefix = lambda wc: working_dir / wc.batch_id / wc.sample_id / 'alignment' / f'{wc.sample_id}'
  resources:
    n_cpus = 10
  shell:
    '''
    ulimit -n 4096
    mkdir -p {params.out_dir}

    # Align reads to the reference genome.
    apptainer exec -B /mnt/backup {app_image} STAR \
      --runMode alignReads \
      --twopassMode Basic \
      --genomeDir {input.star_index_dir} \
      --readFilesIn {input.fastq_r1},{input.fastq_r2} \
      --readFilesCommand zcat \
      --runThreadN {resources.n_cpus} \
      --alignEndsType EndToEnd \
      --alignIntronMax 50000 \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NM NH HI nM AS MD \
      --outFilterMismatchNmax 5 \
      --outFileNamePrefix {params.out_prefix}.primary.

    # Obtain reads per strand
    apptainer exec -B /mnt/backup {app_image} samtools fastq -F 16 {params.out_prefix}.primary.Aligned.sortedByCoord.out.bam \
      | apptainer exec -B /mnt/backup {app_image} seqtk seq -r - \
      | apptainer exec -B /mnt/backup {app_image} gzip > {params.out_dir}/forward.fq.gz
    apptainer exec -B /mnt/backup {app_image} samtools fastq -f 16 {params.out_prefix}.primary.Aligned.sortedByCoord.out.bam \
      | apptainer exec -B /mnt/backup {app_image} gzip > {params.out_dir}/reverse.fq.gz

    # Realign all reads, including reads mapped to forward and revers strands.
    apptainer exec -B /mnt/backup {app_image} STAR \
      --runMode alignReads \
      --twopassMode Basic \
      --genomeDir {input.star_index_dir} \
      --readFilesIn {params.out_dir}/reverse.fq.gz,{params.out_dir}/forward.fq.gz \
      --readFilesCommand zcat \
      --runThreadN {resources.n_cpus} \
      --alignEndsType EndToEnd \
      --alignIntronMax 50000 \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NM NH HI nM AS MD \
      --outFilterMismatchNmax 5 \
      --outFileNamePrefix {params.out_prefix}.final.

    # Create soft links to the final alignments.
    ln -s {params.out_prefix}.final.Aligned.sortedByCoord.out.bam {output.bam_file}
    apptainer exec -B /mnt/backup {app_image} samtools index -@ {resources.n_cpus} {output.bam_file}
    '''


rule s04_overview_of_alignment:
  input:
    unpack(RuleIO(working_dir, "s04_overview_of_alignment", sequencing_info))
  output:
    overview_of_alignment = working_dir / '{batch_id}/{sample_id}/alignment/{sample_id}.overview_of_alignment.csv',
  resources:
    n_cpus = 10
  run:
    n_total_reads, n_clean_reads = -1, -1
    fastq_raw_r1, fastq_clean_r1 = Path(input.fastq_raw_r1), Path(input.fastq_clean_r1)
    if not fastq_raw_r1.exists() or fastq_raw_r1.stat().st_size == 0:
      n_total_reads = -1
    else:
      try:
        n_total_reads = os.popen(f"pigz -dcp {resources.n_cpus} {fastq_raw_r1} | grep . -c").read()
      except Exception as e:
        print("[E]: Failed to get number of RAW reads, due to", e)

    if not fastq_clean_r1.exists() or fastq_clean_r1.stat().st_size == 0:
      n_clean_reads = -1
    else:
      try:
        n_clean_reads = os.popen(f"pigz -dcp {resources.n_cpus} {fastq_clean_r1} | grep . -c").read()
      except Exception as e:
        print("[E]: Failed to get number of CLEAN reads, due to", e)

    required_statistics = [
      "Uniquely mapped reads number",
      "Uniquely mapped reads %",
      "Average mapped length",
      "Number of splices: Total",
      "Number of splices: Annotated (sjdb)",
      "Number of splices: Non-canonical",
      "Number of reads mapped to multiple loci",
      "% of reads mapped to multiple loci",
      "Number of reads unmapped: too many mismatches",
      "Number of reads unmapped: too short",
      "Number of reads unmapped: other",
      "Estimated Number of Cells",
      "Unique Reads in Cells Mapped to Gene",
      "Mean Reads per Cell",
      "Median Reads per Cell",
      "UMIs in Cells",
      "Mean UMI per Cell",
      "Median UMI per Cell",
      "Mean Gene per Cell",
      "Median Gene per Cell",
      "Total Gene Detected",
    ]
    report = ["Number of RAW reads", n_total_reads], ["Number of CLEAN reads", n_clean_reads]
    with (
      open(input.log_final_out, "r") as log_hand,
      open(output.overview_of_alignment, "w") as out_hand
    ):
      out_hand.write(f"Batch ID,{wildcards.batch_id}\n")
      out_hand.write(f"Number of RAW reads,{n_total_reads}\n")

      for line in log_hand:
        if "|" not in line: continue
        *key_list, _, value = [x for x in line.strip().split(" ") if x != ""]
        key = ' '.join(key_list)

        if key in required_statistics:
          out_hand.write(f"{key},{value}\n")

      for line in sum_hand:
        key, value = line.strip().split(",")
        if key in required_statistics:
          out_hand.write(f"{key},{value}\n")


rule s04_estimate_mapping_quality:
  input:
    unpack(RuleIO(working_dir, "s04_estimate_mapping_quality", sequencing_info)),
    gtf_file = genomic_feature_file,
    genome_sequence_file = genome_sequence_file
  output:
    qorts_done = working_dir / '{batch_id}/{sample_id}/mapping_reports/qorts/qorts.done',
    qorts_dir = directory(working_dir / '{batch_id}/{sample_id}/mapping_reports/qorts'),
    qualimap_done = working_dir / '{batch_id}/{sample_id}/mapping_reports/qualimap/qualimap.done',
    qualimap_dir = directory(working_dir / '{batch_id}/{sample_id}/mapping_reports/qualimap'),
  shell:
    '''
    mkdir -p {output.qualimap_dir} {output.qorts_dir}

    apptainer exec -B /mnt/backup {app_image} qualimap --java-mem-size=16G rnaseq \
      --sorted -bam {input.bam_file} -gtf {input.gtf_file} -outdir {output.qualimap_dir}

    apptainer exec -B /mnt/backup {app_image} qorts QC \
      --genomeFA {input.genome_sequence_file} \
      --addFunctions mismatchEngine,cigarMatch,referenceMatch,writeDocs,calcDetailedGeneCounts \
      --maxReadLength 250 --singleEnded --generatePlots \
      {input.bam_file} {input.gtf_file} {output.qorts_dir}

    touch {output.qorts_done} {output.qualimap_done}
    '''


rule s05_remove_duplicates:
  input:
    unpack(RuleIO(working_dir, "s05_remove_duplicates", sequencing_info)),
  output:
    bam_file = working_dir / '{batch_id}/{sample_id}/alignment/{sample_id}.rmdup.bam',
    bai_file = working_dir / '{batch_id}/{sample_id}/alignment/{sample_id}.rmdup.bam.bai',
  resources:
    n_cpus = 10
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} samtools sort -n -@ 2 {input.bam_file} \
      | apptainer exec -B /mnt/backup {app_image} samtools fixmate -mc --no-PG -@ 2 - - \
      | apptainer exec -B /mnt/backup {app_image} samtools sort --no-PG -@ 2 \
      | apptainer exec -B /mnt/backup {app_image} samtools markdup --no-PG -@ 2 - - \
      | apptainer exec -B /mnt/backup {app_image} samtools view -F 1024 -o {output.bam_file}
    apptainer exec -B /mnt/backup {app_image} samtools index {output.bam_file}

    # Clean up
    echo rm -f {input.bam_file} {input.bai_file}
    echo touch {input.bam_file} {input.bai_file} {output.bam_file} {output.bai_file}
    '''


rule s06_correct_base_qualities:
  input:
    unpack(RuleIO(working_dir, "s06_correct_base_qualities", sequencing_info))
  output:
    bam_file = working_dir / '{batch_id}/{sample_id}/alignment/{sample_id}.rmdup.base_corrected.bam',
    bai_file = working_dir / '{batch_id}/{sample_id}/alignment/{sample_id}.rmdup.base_corrected.bam.bai',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  resources:
    n_cpus = 10
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} extract_new_reads.py \
      -@ {resources.n_cpus} -o {params.out_dir}/old_reads.bam -n {params.out_dir}/new_reads.bam {input.bam_file}

    apptainer exec -B /mnt/backup {app_image} samtools cat {params.out_dir}/{{old_reads,new_reads}}.bam \
      | apptainer exec -B /mnt/backup {app_image} samtools sort -@ {resources.n_cpus} -o {output.bam_file}
    samtools index {output.bam_file}

    rm -f {params.out_dir}/{{old_reads,new_reads}}.bam
    '''


rule s07_quantify_enhancer_rna:
  input:
    unpack(RuleIO(working_dir, "s07_quantify_enhancer_rna", sequencing_info)), genomic_regions = erna_feature_file,
  output:
    count_tab = working_dir / '{batch_id}/{sample_id}/quantification/{sample_id}.enhancer_rna.read_counts.txt',
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  resources:
    n_cpus = 8
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} featureCounts \
      -p -t enhancer -g enhancer_region -T {resources.n_cpus} \
      -a {input.genomic_regions} -o {output.count_tab} \
      {input.bam_file}
    '''


rule s07_quantify_genes:
  input:
    unpack(RuleIO(working_dir, "s07_quantify_genes", sequencing_info)), genomic_regions = genomic_feature_file,
  output:
    count_tab = working_dir / '{batch_id}/{sample_id}/quantification/{sample_id}.genes.read_counts.txt',
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  resources:
    n_cpus = 8
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} featureCounts \
      -p -t gene -g gene_name --extraAttributes gene_id -T {resources.n_cpus} \
      -a {input.genomic_regions} -o {output.count_tab} {input.bam_file}
    '''


rule s07_convert_bam_to_cit:
  input:
    unpack(RuleIO(working_dir, "s07_convert_bam_to_cit", sequencing_info)),
  output:
    cit_file = working_dir / '{batch_id}/{sample_id}/cit_files/{sample_id}.cit',
    cit_corr_file = working_dir / '{batch_id}/{sample_id}/cit_files/{sample_id}.corrected.cit',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} gedi -e Bam2CIT {output.cit_file} {input.bam_file}
    apptainer exec -B /mnt/backup {app_image} gedi -e CorrectCIT {output.cit_file} {output.cit_corr_file}

    mv {output.cit_corr_file}.metadata.json {output.cit_corr_file}.metadata.json.bk
    jq '.[] | {{conditions: [.[] | {{total: .total, name: "{wildcards.sample_id}" }}]}}' \
      < {output.cit_corr_file}.metadata.json.bk \
      > {output.cit_corr_file}.metadata.json

    echo rm -f {output.cit_file}
    echo touch {output.cit_file} {output.cit_corr_file}
    '''


rule s08_merge_cits_per_treatment:
  input:
    unpack(RuleIO(working_dir, "s08_merge_cits_per_treatment", sequencing_info))
  output:
    merged_cit = working_dir / '{batch_id}/{sample_id}/cit_files/{sample_id}.corrected.with_ctrl.cit',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    external_controls = lambda wc: get_parameters(wc.batch_id, "external_controls")
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} \
      gedi -e MergeCIT {output.merged_cit} {input.ctrl_cits} {input.sample_cit} {params.external_controls}
    apptainer exec -B /mnt/backup {app_image} gedi -e ReadCount {output.merged_cit} 
    '''


rule s09_estimate_nascent_rna_per_treatment:
  input:
    unpack(RuleIO(working_dir, "s09_estimate_nascent_rna_per_treatment", sequencing_info))
  output:
    nascent_rna = working_dir / '{batch_id}/{sample_id}/nascent_rna/{sample_id}.tsv.gz',
  params:
    out_prefix = lambda _, output: output.nascent_rna.strip(".tsv.gz"),
    out_dir = lambda _, output: Path(output[0]).parent,
    ctrl_pattern = lambda wc: get_parameters(wc.batch_id, "control_pattern")
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} \
      gedi -e Slam -genomic GEDI_GRCm38 -reads {input.merged_cit} -prefix {params.out_prefix} \
      -minEstimateReads 1000 -trim3p 30 -trim5p 20 -plot -D -full -allGenes \
      -no4sUpattern {params.ctrl_pattern}
    '''

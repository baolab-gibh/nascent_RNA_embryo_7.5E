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
      for per_batch, sample_info_dict in self._info_db.items():
        if "quality_control" in self._rule: # s03_quality_control
          outs.append(self._wk_dir / f'{per_batch}/quality_control/{per_batch}.decoded.clean.I1.fastq.gz')
          outs.append(self._wk_dir / f'{per_batch}/quality_control/{per_batch}.decoded.clean.R1.fastq.gz')
          outs.append(self._wk_dir / f'{per_batch}/quality_control/{per_batch}.decoded.clean.I2.fastq.gz')
          outs.append(self._wk_dir / f'{per_batch}/quality_control/{per_batch}.decoded.clean.R2.fastq.gz')
        if "alignment" in self._rule: # s04_align_reads
          outs.append(self._wk_dir/f'{per_batch}/alignment/{per_batch}.star_solo.bam')
        if "alignment_overview" in self._rule: # s05_overview_of_alignment
          outs.append(self._wk_dir / f'{per_batch}/alignment/{per_batch}.overview_of_alignment.csv')
        if "alignment_reports" in self._rule: # s05_estimate_mapping_quality
          outs.append(self._wk_dir/f'{per_batch}/mapping_reports/qualimap/qualimap.done')
          outs.append(self._wk_dir/f'{per_batch}/mapping_reports/qorts/qorts.done')
        if "nascent_rna" in self._rule: # s10_estimate_nascent_rna
          outs.append(self._wk_dir/f'{per_batch}/nascent_rna/all_samples/{per_batch}.tsv.gz')
        # Per sample outputs
        for per_sample, sample_info in sample_info_dict.get("samples").items():
          if "quantify_enhancer_rna" in self._rule: # s08_quantify_enhancer_rna
            outs.append(self._wk_dir/f'{per_batch}/quantification/per_sample/{per_sample}.enhancer_rna.read_counts.txt')
          if "quantify_genes" in self._rule:
            outs.append(self._wk_dir/f'{per_batch}/quantification/per_sample/{per_sample}.genes.read_counts.txt')
          if "nascent_rna_per_treatment" in self._rule: # s10_estimate_nascent_rna_per_treatment
            if per_sample not in sample_info_dict["parameters"]["control_samples"]:
              outs.append(self._wk_dir/f'{per_batch}/nascent_rna/per_sample/{per_sample}/{per_sample}.tsv.gz')
      return sorted(outs)
    elif isinstance(wc, Wildcards) and isinstance(self._rule, str): # Determine inputs based on wildcards and rule name.
      batch_id, sample_id = getattr(wc, "batch_id", None), getattr(wc, "sample_id", None)
      if self._rule in ["s01_normalize_fastq_name"]:
        return dict(
          fastq_r1= self._info_db.get(batch_id).get("fastqs").get("fastq_r1"),
          fastq_r2= self._info_db.get(batch_id).get("fastqs").get("fastq_r2"),
        )
      if self._rule in ["s02_preproc_multiplexed_reads"]:
        return dict(
          fastq_r1=self._wk_dir / f'{batch_id}/fastq/{batch_id}.raw.R1.fastq.gz',
          fastq_r2=self._wk_dir / f'{batch_id}/fastq/{batch_id}.raw.R2.fastq.gz'
        )
      if self._rule in ["s03_quality_control"]:
        return dict(
          fastq_i1=self._wk_dir / f'{batch_id}/decode/{batch_id}.decoded.I1.fastq.gz',
          fastq_r1=self._wk_dir / f'{batch_id}/decode/{batch_id}.decoded.R1.fastq.gz',
          fastq_i2=self._wk_dir / f'{batch_id}/decode/{batch_id}.decoded.I2.fastq.gz',
          fastq_r2=self._wk_dir / f'{batch_id}/decode/{batch_id}.decoded.R2.fastq.gz',
        )
      if self._rule == "s04_align_reads":
        return dict(
          fastq_i1=self._wk_dir / f'{batch_id}/quality_control/{batch_id}.decoded.clean.I1.fastq.gz',
          fastq_r1=self._wk_dir / f'{batch_id}/quality_control/{batch_id}.decoded.clean.R1.fastq.gz',
          fastq_i2=self._wk_dir / f'{batch_id}/quality_control/{batch_id}.decoded.clean.I2.fastq.gz',
          fastq_r2=self._wk_dir / f'{batch_id}/quality_control/{batch_id}.decoded.clean.R2.fastq.gz',
        )
      if self._rule == "s05_overview_of_alignment":
        return dict(
          fastq_raw_r1=self._wk_dir / f'{batch_id}/fastq/{batch_id}.raw.R1.fastq.gz',
          fastq_clean_r1=self._wk_dir / f'{batch_id}/quality_control/{batch_id}.decoded.clean.R1.fastq.gz',
          log_final_out = self._wk_dir / f'{batch_id}/alignment/{batch_id}.Log.final.out',
          star_solo_summary = self._wk_dir / f'{batch_id}/alignment/{batch_id}.Solo.out/Gene/Summary.csv',
        )
      if self._rule in ["s05_estimate_mapping_quality"]:
        return dict(
          bam_file=self._wk_dir / f'{batch_id}/alignment/{batch_id}.star_solo.bam',
          bai_file=self._wk_dir / f'{batch_id}/alignment/{batch_id}.star_solo.bam.bai'
        )
      if self._rule in ["s05_split_bam_by_bc",]:
        return dict(
          bam_file=self._wk_dir / f'{batch_id}/alignment/{batch_id}.star_solo.bam',
          bai_file=self._wk_dir / f'{batch_id}/alignment/{batch_id}.star_solo.bam.bai',
          valid_barcode_file=self._wk_dir / f'{batch_id}/barcodes/{batch_id}.valid_barcodes.txt',
          sample_barcode_file=self._wk_dir / f'{batch_id}/barcodes/{batch_id}.sample_barcode_mapping.txt',
        )
      if self._rule in ["s06_remove_duplicates"]:
        return dict(
          bam_file=self._wk_dir / f'{batch_id}/alignment_persample/{sample_id}.star_solo.bam',
          bai_file=self._wk_dir / f'{batch_id}/alignment_persample/{sample_id}.star_solo.bam.bai',
        )
      if self._rule in ["s07_correct_base_qualities"]:
        return dict(
          bam_file=self._wk_dir / f'{batch_id}/alignment_persample/{sample_id}.star_solo.rmdup.bam',
          bai_file=self._wk_dir / f'{batch_id}/alignment_persample/{sample_id}.star_solo.rmdup.bam.bai',
        )
      if self._rule in ["s08_quantify_enhancer_rna", "s08_quantify_genes", "s08_convert_bam_to_cit"]:
        return dict(
          bam_file=self._wk_dir / f'{batch_id}/alignment_persample/{sample_id}.star_solo.rmdup.base_corrected.bam',
          bai_file=self._wk_dir / f'{batch_id}/alignment_persample/{sample_id}.star_solo.rmdup.base_corrected.bam.bai'
        )
      if self._rule in ["s09_merge_cits"]:
        input_files = []
        for per_batch, sample_info_dict in self._info_db.items():
          if batch_id != per_batch: continue
          for per_sample in sample_info_dict["samples"].keys():
            input_files.append(self._wk_dir / f'{batch_id}/cit_files/per_sample/{per_sample}.corrected.cit')
        return input_files
      if self._rule in ['s10_estimate_nascent_rna']:
        return dict(
          merged_cit=self._wk_dir / f'{batch_id}/cit_files/merged/{batch_id}.all_samples.cit'
        )
      if self._rule in ["s09_merge_cits_per_treatment"]:
        ctrl_cits = []
        for per_batch, sample_info_dict in self._info_db.items():
          if per_batch != batch_id: continue
          for ctrl_sample in sample_info_dict.get("parameters", {}).get("control_samples", []):
            ctrl_cits.append(self._wk_dir / f'{batch_id}/cit_files/per_sample/{ctrl_sample}.corrected.cit')
        return dict(
          ctrl_cits=ctrl_cits,
          sample_cit=self._wk_dir / f'{batch_id}/cit_files/per_sample/{sample_id}.corrected.cit'
        )
      if self._rule in ["s10_estimate_nascent_rna_per_treatment"]:
        return dict(merged_cit=self._wk_dir / f'{batch_id}/cit_files/merged/{sample_id}.corrected.with_ctrl.cit')
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
config_file = Path("~/Documents/projects/wp_vasaseq/scripts/snakemake/geo_slam.configuration.yaml").expanduser()
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

# barcode_pattern = parameters.get("barcode_pattern")

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
    apptainer exec -B /mnt/backup {app_image} gedi -e IndexGenome \
      -s mus_musculus.90.fasta -g mus_musculus.90.gtf -p -nobowtie -nostar -nokallisto
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
    apptainer exec -B /mnt/backup {app_image} STAR --runMode genomeGenerate \
      --runThreadN {resources.n_cpus} \
      --sjdbGTFfile {input.genomic_feature_file} \
      --genomeFastaFiles {input.genome_sequence_file} \
      --genomeDir {output.star_index_dir}
    '''


rule s01_create_barcodes:
  input:
    config_file,
  output:
    sample_barcode_file = working_dir / '{batch_id}/barcodes/{batch_id}.sample_barcode_mapping.txt',
    valid_barcode_file = working_dir / '{batch_id}/barcodes/{batch_id}.valid_barcodes.txt',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  run:
    params.out_dir.mkdir(parents=True, exist_ok=True)

    batch_id = wildcards.batch_id
    batch_info = sequencing_info.get(batch_id)
    with open(output.sample_barcode_file, "w") as out_1, open(output.valid_barcode_file, "w") as out_2:
      for sample_id, sample_info in batch_info.get("samples").items():
        barcode = sample_info.get("barcodes")
        out_1.write(f"{sample_id}\t{barcode}\n")
        out_2.write(f"{barcode}\n")


rule s01_normalize_fastq_name:
  input:
    unpack(RuleIO(working_dir, "s01_normalize_fastq_name", sequencing_info))
  output:
    fastq_r1 = working_dir / '{batch_id}/fastq/{batch_id}.raw.R1.fastq.gz',
    fastq_r2 = working_dir / '{batch_id}/fastq/{batch_id}.raw.R2.fastq.gz',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  resources:
    n_cpus = 10
  shell:
    '''
    mkdir -p {params.out_dir}

    ln -s {input.fastq_r1} {output.fastq_r1}; ln -s {input.fastq_r2} {output.fastq_r2}
    apptainer exec -B /mnt/backup {app_image} fastqc -q -t {resources.n_cpus} -o {params.out_dir} \
      {output.fastq_r1} {output.fastq_r2}
    '''


rule s02_preproc_multiplexed_reads:
  input:
    unpack(RuleIO(working_dir, "s02_preproc_multiplexed_reads", sequencing_info))
  output:
    fastq_i1 = working_dir / '{batch_id}/decode/{batch_id}.decoded.I1.fastq.gz',
    fastq_r1 = working_dir / '{batch_id}/decode/{batch_id}.decoded.R1.fastq.gz',
    fastq_i2 = working_dir / '{batch_id}/decode/{batch_id}.decoded.I2.fastq.gz',
    fastq_r2 = working_dir / '{batch_id}/decode/{batch_id}.decoded.R2.fastq.gz',
  resources:
    n_cpus = 10
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    barcode_pattern = lambda wc: get_parameters(wc.batch_id, "barcode_pattern")
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} preprocess_barcoded_reads.py \
      -o {params.out_dir} -n {wildcards.batch_id}.decoded -N 4 -Z 10 -p {params.barcode_pattern} \
      {input.fastq_r1} {input.fastq_r2}

    ln -s {output.fastq_i1} {output.fastq_i2}
    apptainer exec -B /mnt/backup {app_image} fastqc \
      -q -t {resources.n_cpus} -o {params.out_dir} {output.fastq_i1} {output.fastq_r1} {output.fastq_r2}
    '''


rule s03_quality_control_r1: # For read 1.
  input:
    unpack(RuleIO(working_dir, "s03_quality_control", sequencing_info))
  output:
    fastq_i1 = working_dir / '{batch_id}/quality_control/{batch_id}.decoded.clean.I1.fastq.gz',
    fastq_r1 = working_dir / '{batch_id}/quality_control/{batch_id}.decoded.clean.R1.fastq.gz',
  resources:
    n_cpus = 5
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    adapter_sequence = lambda wc: get_parameters(wc.batch_id, "adapter_sequence_rev_cmp")
  shell:
    '''
    mkdir -p {params.out_dir}

    # Removing adapter, polyX. Trimming low quality bases. Removing 'N'. Ensure minimum read length (80).
    apptainer exec -B /mnt/backup {app_image} cutadapt \
      --quiet -j {resources.n_cpus} -n 10 -m 14:50 -q 0 -Q 25 -e 1 --max-n 0 --trim-n \
      -A 'A{{6}}' -G {params.adapter_sequence} \
      -A 'A{{4}}$' -A 'A{{2}}$' -A 'A{{1}}$' \
      -A 'G{{4}}$' -A 'G{{2}}$' -A 'G{{1}}$' \
      -o {params.out_dir}/001.I1.fq.gz -p {params.out_dir}/001.R2.fq.gz \
      {input.fastq_i1} {input.fastq_r1}

    # Removing 5bp on per end
    apptainer exec -B /mnt/backup {app_image} cutadapt \
      --quiet -j {resources.n_cpus} -m 14:60 -U 5 -U -5 \
      -o {output.fastq_i1} -p {output.fastq_r1} \
      {params.out_dir}/001.I1.fq.gz {params.out_dir}/001.R2.fq.gz
    apptainer exec -B /mnt/backup {app_image} fastqc -q -t {resources.n_cpus} -o {params.out_dir} {output.fastq_r1}

    # Clean up
    rm -f {params.out_dir}/001.I1.fq.gz {params.out_dir}/001.R2.fq.gz
    echo rm -f {input.fastq_i1} {input.fastq_r1}
    echo touch {input.fastq_i1} {input.fastq_r1} {output.fastq_r1} {output.fastq_i1}
    '''


rule s03_quality_control_r2: # For read 2
  input:
    unpack(RuleIO(working_dir, "s03_quality_control", sequencing_info))
  output:
    fastq_i2 = working_dir / '{batch_id}/quality_control/{batch_id}.decoded.clean.I2.fastq.gz',
    fastq_r2 = working_dir / '{batch_id}/quality_control/{batch_id}.decoded.clean.R2.fastq.gz',
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
      --quiet -j {resources.n_cpus} -n 10 -m 14:50 -q 0 -Q 25 -e 1 --max-n 0 --trim-n \
      -A 'A{{10}}' -A {params.adapter_sequence} \
      -A 'A{{4}}$' -A 'A{{2}}$' -A 'A{{1}}$' \
      -A 'G{{4}}$' -A 'G{{2}}$' -A 'G{{1}}$' \
      -o {params.out_dir}/002.I1.fq.gz -p {params.out_dir}/002.R2.fq.gz \
      {input.fastq_i2} {input.fastq_r2}

    # Removing 5bp on per end
    apptainer exec -B /mnt/backup {app_image} cutadapt \
      --quiet -j {resources.n_cpus} -m 14:60 -U 5 -U -5 \
      -o {output.fastq_i2} -p {output.fastq_r2} \
      {params.out_dir}/002.I1.fq.gz {params.out_dir}/002.R2.fq.gz
    apptainer exec -B /mnt/backup {app_image} fastqc -q -t {resources.n_cpus} -o {params.out_dir} {output.fastq_r2}

    # Clean up
    rm -f {params.out_dir}/002.I1.fq.gz {params.out_dir}/002.R2.fq.gz
    echo rm -f {input.fastq_i2} {input.fastq_r2}
    echo touch {input.fastq_i2} {input.fastq_r2} {output.fastq_r2} {output.fastq_i2}
    '''


rule s04_align_reads:
  input:
    unpack(RuleIO(working_dir, "s04_align_reads", sequencing_info)), star_index_dir = star_index_dir,
  output:
    bam_file = working_dir / '{batch_id}/alignment/{batch_id}.star_solo.bam',
    bai_file = working_dir / '{batch_id}/alignment/{batch_id}.star_solo.bam.bai',
    log_final_out = working_dir / '{batch_id}/alignment/{batch_id}.Log.final.out',
    star_solo_summary = working_dir / '{batch_id}/alignment/{batch_id}.Solo.out/Gene/Summary.csv',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    out_prefix = lambda wc: working_dir / wc.batch_id / 'alignment' / (wc.batch_id + ".")
  resources:
    n_cpus = 10
  priority: 100
  shell:
    '''
    ulimit -n 4096
    mkdir -p {params.out_dir}

    echo -e Aligning using Solo mode ...
    apptainer exec -B /mnt/backup {app_image} STAR --runMode alignReads \
      --twopassMode Basic \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist None \
      --soloCBstart 1 --soloCBlen 6 \
      --soloUMIstart 7 --soloUMIlen 8 \
      --soloFeatures Gene GeneFull SJ Velocyto \
      --genomeDir {input.star_index_dir} \
      --readFilesIn {input.fastq_r2},{input.fastq_r1} {input.fastq_i2},{input.fastq_i1} \
      --readFilesCommand zcat \
      --runThreadN {resources.n_cpus} \
      --alignEndsType EndToEnd \
      --alignIntronMax 50000 \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NM NH HI nM AS CR UR CB UB GX GN MD \
      --outFilterMismatchNmax 10 \
      --outFileNamePrefix {params.out_prefix}

    echo -e Creating index ...
    ln -s {params.out_prefix}Aligned.sortedByCoord.out.bam {output.bam_file}
    samtools index -@ {resources.n_cpus} {output.bam_file}
    '''


rule s05_overview_of_alignment:
  input:
    unpack(RuleIO(working_dir, "s05_overview_of_alignment", sequencing_info))
  output:
    overview_of_alignment = working_dir / '{batch_id}/alignment/{batch_id}.overview_of_alignment.csv',
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
      "Uniquely mapped reads number", "Uniquely mapped reads %", "Average mapped length",
      "Number of splices: Total", "Number of splices: Annotated (sjdb)", "Number of splices: Non-canonical",
      "Number of reads mapped to multiple loci", "% of reads mapped to multiple loci",
      "Number of reads unmapped: too many mismatches", "Number of reads unmapped: too short",
      "Number of reads unmapped: other", "Estimated Number of Cells",
      "Unique Reads in Cells Mapped to Gene", "Mean Reads per Cell", "Median Reads per Cell",
      "UMIs in Cells", "Mean UMI per Cell", "Median UMI per Cell",
      "Mean Gene per Cell", "Median Gene per Cell", "Total Gene Detected"
    ]
    report = ["Number of RAW reads", n_total_reads], ["Number of CLEAN reads", n_clean_reads]
    with (
      open(input.log_final_out, "r") as log_hand,
      open(input.star_solo_summary, "r") as sum_hand,
      open(output.overview_of_alignment, "w") as out_hand
    ):
      out_hand.write(f"Batch ID,{wildcards.batch_id}\n")
      out_hand.write(f"Number of RAW reads,{n_total_reads}\n")

      for line in log_hand:
        if "|" not in line: continue

        line_list = [x for x in line.strip().split(" ") if x != ""]
        *key_list, _, value = line_list
        key = ' '.join(key_list)

        print(line)
        if key in required_statistics:
          out_hand.write(f"{key},{value}\n")

      for line in sum_hand:
        key, value = line.strip().split(",")
        if key in required_statistics:
          out_hand.write(f"{key},{value}\n")


rule s05_estimate_mapping_quality:
  input:
    unpack(RuleIO(working_dir, "s05_estimate_mapping_quality", sequencing_info)),
    gtf_file = genomic_feature_file,
    genome_sequence_file = genome_sequence_file
  output:
    qorts_done = working_dir / '{batch_id}/mapping_reports/qorts/qorts.done',
    qualimap_done = working_dir / '{batch_id}/mapping_reports/qualimap/qualimap.done',
    qorts_dir = directory(working_dir / '{batch_id}/mapping_reports/qorts'),
    qualimap_dir = directory(working_dir / '{batch_id}/mapping_reports/qualimap'),
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


rule s05_split_bam_by_bc:
  input:
    unpack(RuleIO(working_dir, "s05_split_bam_by_bc", sequencing_info)),
  output:
    bam_file = working_dir / '{batch_id}/alignment_persample/{sample_id}.star_solo.bam',
    bai_file = working_dir / '{batch_id}/alignment_persample/{sample_id}.star_solo.bam.bai',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    sample_barcode = lambda wc: sequencing_info.get(wc.batch_id, {}).get("samples", {}).get(wc.sample_id, {}).get("barcodes", ""),
  resources:
    n_cpus = 10
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} samtools view \
      -d CR:{params.sample_barcode} -O BAM -o {output.bam_file} -@ {resources.n_cpus} {input.bam_file}
    apptainer exec -B /mnt/backup {app_image} samtools index -@ {resources.n_cpus} {output.bam_file}
    '''


rule s06_remove_duplicates:
  input:
    unpack(RuleIO(working_dir, "s06_remove_duplicates", sequencing_info)),
  output:
    bam_file = working_dir / '{batch_id}/alignment_persample/{sample_id}.star_solo.rmdup.bam',
    bai_file = working_dir / '{batch_id}/alignment_persample/{sample_id}.star_solo.rmdup.bam.bai',
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


rule s07_correct_base_qualities:
  input:
    unpack(RuleIO(working_dir, "s07_correct_base_qualities", sequencing_info))
  output:
    bam_file = working_dir / '{batch_id}/alignment_persample/{sample_id}.star_solo.rmdup.base_corrected.bam',
    bai_file = working_dir / '{batch_id}/alignment_persample/{sample_id}.star_solo.rmdup.base_corrected.bam.bai',
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
    apptainer exec -B /mnt/backup {app_image} samtools index {output.bam_file}

    rm -f {params.out_dir}/{{old_reads,new_reads}}.bam
    '''


rule s08_quantify_enhancer_rna:
  input:
    unpack(RuleIO(working_dir, "s08_quantify_enhancer_rna", sequencing_info)), genomic_regions = erna_feature_file,
  output:
    count_tab = working_dir / '{batch_id}/quantification/per_sample/{sample_id}.enhancer_rna.read_counts.txt',
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  resources:
    n_cpus = 8
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} featureCounts \
      -t enhancer -g enhancer_region \
      -T {resources.n_cpus} -a {input.genomic_regions} -o {output.count_tab} {input.bam_file}
    '''


rule s08_quantify_genes:
  input:
    unpack(RuleIO(working_dir, "s08_quantify_genes", sequencing_info)), genomic_regions = genomic_feature_file,
  output:
    count_tab = working_dir / '{batch_id}/quantification/per_sample/{sample_id}.genes.read_counts.txt',
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  resources:
    n_cpus = 8
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} featureCounts \
      -t gene -g gene_name --extraAttributes gene_id -T {resources.n_cpus} \
      -a {input.genomic_regions} -o {output.count_tab} {input.bam_file}
    '''


rule s08_convert_bam_to_cit:
  input:
    unpack(RuleIO(working_dir, "s08_convert_bam_to_cit", sequencing_info)),
  output:
    cit_file = working_dir / '{batch_id}/cit_files/per_sample/{sample_id}.cit',
    cit_corr_file = working_dir / '{batch_id}/cit_files/per_sample/{sample_id}.corrected.cit',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} gedi -e Bam2CIT -p {output.cit_file} {input.bam_file}
    apptainer exec -B /mnt/backup {app_image} gedi -e CorrectCIT -p {output.cit_file} {output.cit_corr_file}

    mv {output.cit_corr_file}.metadata.json {output.cit_corr_file}.metadata.json.bk
    jq '.[] | {{conditions: [.[] | {{total: .total, name: "{wildcards.sample_id}" }}]}}' \
      < {output.cit_corr_file}.metadata.json.bk \
      > {output.cit_corr_file}.metadata.json

    echo rm -f {output.cit_file}
    echo touch {output.cit_file} {output.cit_corr_file}
    '''


rule s09_merge_cits:
  input:
    unpack(RuleIO(working_dir, "s09_merge_cits", sequencing_info))
  output:
    merged_cit = working_dir / '{batch_id}/cit_files/merged/{batch_id}.all_samples.cit',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    external_controls = lambda wc: get_parameters(wc.batch_id, "external_controls")
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} gedi -e MergeCIT -p {output.merged_cit} {input} {params.external_controls}
    apptainer exec -B /mnt/backup {app_image} gedi -e ReadCount -p {output.merged_cit}
    '''


rule s10_estimate_nascent_rna:
  input:
    unpack(RuleIO(working_dir, "s10_estimate_nascent_rna", sequencing_info))
  output:
    nascent_rna = working_dir / '{batch_id}/nascent_rna/all_samples/{batch_id}.tsv.gz',
  params:
    # out_prefix = lambda wc, output: Path(output[0]).parent / f'{wc.batch_id}.nascent_rna',
    out_prefix = lambda _, output: output.nascent_rna.strip(".tsv.gz"),
    out_dir = lambda _, output: Path(output[0]).parent,
    ctrl_pattern = lambda wc: get_parameters(wc.batch_id, "control_pattern"),
  priority: 100
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} gedi -e Slam \
      -genomic GEDI_GRCm38 -reads {input.merged_cit} -prefix {params.out_prefix} \
      -minEstimateReads 1000 -trim3p 30 -trim5p 20 -progress -plot -D -full -allGenes \
      -no4sUpattern {params.ctrl_pattern}
    '''


rule s09_merge_cits_per_treatment:
  input:
    unpack(RuleIO(working_dir, "s09_merge_cits_per_treatment", sequencing_info))
  output:
    merged_cit = working_dir / '{batch_id}/cit_files/merged/{sample_id}.corrected.with_ctrl.cit',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    external_controls = lambda wc: get_parameters(wc.batch_id, "external_controls")
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} gedi -e MergeCIT -p {output.merged_cit} \
      {input.ctrl_cits} {input.sample_cit} {params.external_controls}
    apptainer exec -B /mnt/backup {app_image} gedi -e ReadCount -p {output.merged_cit} 
    '''


rule s10_estimate_nascent_rna_per_treatment:
  input:
    unpack(RuleIO(working_dir, "s10_estimate_nascent_rna_per_treatment", sequencing_info))
  output:
    nascent_rna = working_dir / '{batch_id}/nascent_rna/per_sample/{sample_id}/{sample_id}.tsv.gz',
  params:
    out_prefix = lambda _, output: output.nascent_rna.strip(".tsv.gz"),
    out_dir = lambda _, output: Path(output[0]).parent,
    ctrl_pattern = lambda wc: get_parameters(wc.batch_id, "control_pattern")
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} gedi -e Slam \
      -genomic GEDI_GRCm38 -reads {input.merged_cit} -prefix {params.out_prefix} \
      -minEstimateReads 1000 -trim3p 30 -trim5p 20 -progress -plot -D -full -allGenes \
      -no4sUpattern {params.ctrl_pattern}
    '''

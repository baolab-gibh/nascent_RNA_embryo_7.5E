#!/usr/bin/env snakemake -f
# File: celatlas_spatial.pipeline.smk
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Mar 11, 2025
# Updated: Apr 02, 2025

import os
import sys
import csv
import yaml
from pathlib import Path
from snakemake.utils import min_version
from snakemake.io import Wildcards

min_version('6.0') # Ensure the snakemake knows how to handle moduliazations.


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
      global bin_size
      for per_batch, sample_info_dict in self._info_db.items():
        for per_sample, sample_info in sample_info_dict.get("samples").items():
          # Per batch outputs
          if "normalize_fastq_name" in self._rule:
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/00.fastqs/{per_sample}_1.fq.gz')
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/00.fastqs/{per_sample}_2.fq.gz')
          if "quality_control" in self._rule:
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/01.barcode/stat.txt')
          if "alignment" in self._rule:
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/03.star/{per_sample}_Aligned.sortedByCoord.out.bam')
          if "spatial" in self._rule:
            has_image = sample_info.get("image_file", None) is not None
            analysis_dir = f"07.analysis_bin{bin_size}_with_image" if has_image else f"07.analysis_bin{bin_size}"
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/{analysis_dir}/Bioinfodata')
            for ext in ["markers.tsv", "markers_raw.tsv", "tsne_coord.tsv", "bin10.h5ad", "bin20.h5ad", "bin50.h5ad", "bin100.h5ad"]:
              outs.append(self._wk_dir/f'{per_batch}/{per_sample}/{analysis_dir}/{per_sample}_{ext}')
          if "alignment_reports" in self._rule:
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/04.mapping_reports/{per_sample}_qorts/qorts.done'),
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/04.mapping_reports/{per_sample}_qorts'),
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/04.mapping_reports/{per_sample}_qualimap/qualimap.done'),
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/04.mapping_reports/{per_sample}_qualimap'),
          # Per sample outputs
          if "quantification_erna" in self._rule:
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/04.quantification/per_sample/enhancer_rna.read_counts.txt')
          if "nascent_rna" in self._rule:
            outs.append(self._wk_dir/f'{per_batch}/{per_sample}/06.binSegment_nascent/{per_sample}_Barcodes_tissue_positions.csv')
          if "nascent_rna_per_treatment" in self._rule:
            if per_sample not in sample_info_dict["parameters"]["control_samples"]:
              outs.append(self._wk_dir/f'{per_batch}/{per_sample}/07.nascent_rna/{per_sample}.tsv.gz')
      return sorted(outs)
    elif isinstance(wc, Wildcards) and isinstance(self._rule, str): # Determine inputs based on wildcards and rule name.
      per_batch, per_sample = getattr(wc, "batch_id", None), getattr(wc, "sample_id", None)
      # Standarded Celatlas Spatial pipeline.
      if self._rule in ["s00_normalize_fastq"]:
        return dict(
          fastq_r1=self._info_db[per_batch]["samples"][per_sample]["fastqs"]["fastq_r1"],
          fastq_r2=self._info_db[per_batch]["samples"][per_sample]["fastqs"]["fastq_r2"],
        )
      if self._rule in ["s01_decode_fastq"]:
        return dict(
          fastq_r1 = self._wk_dir / f'{per_batch}/{per_sample}/00.fastqs/{per_sample}_1.fq.gz',
          fastq_r2 = self._wk_dir / f'{per_batch}/{per_sample}/00.fastqs/{per_sample}_2.fq.gz',
          barcode_file=self._info_db[per_batch]["samples"][per_sample]["barcode_file"],
        )
      if self._rule in ["s02_trimming_fastq"]:
        return dict(
          fastq_r1=self._wk_dir/f'{per_batch}/{per_sample}/01.barcode/{per_sample}_1.fq.gz',
          fastq_r2=self._wk_dir/f'{per_batch}/{per_sample}/01.barcode/{per_sample}_2.fq.gz',
        )
      if self._rule in ["s03_align_reads"]:
        return dict(fastq_r2=self._wk_dir/f'{per_batch}/{per_sample}/02.cutadapt/{per_sample}_clean_2.fq.gz',)
      if self._rule in ["s04_count_reads"]:
        return dict(
          bam_file=self._wk_dir/f'{per_batch}/{per_sample}/03.star/{per_sample}_Aligned.sortedByCoord.out.bam',
          bai_file=self._wk_dir/f'{per_batch}/{per_sample}/03.star/{per_sample}_Aligned.sortedByCoord.out.bam.bai',
        )
      if self._rule in ["s04_count_reads_nascent"]:
        return dict(
          bam_file=self._wk_dir/f'{per_batch}/{per_sample}/03.star/{per_sample}_Aligned.sortedByCoord.out.nascent.bam',
          bai_file=self._wk_dir/f'{per_batch}/{per_sample}/03.star/{per_sample}_Aligned.sortedByCoord.out.nascent.bam.bai',
        )
      if self._rule in ["s05_generate_count_details"]:
        return dict(
          bam_file=self._wk_dir/f'{per_batch}/{per_sample}/04.featureCounts/{per_sample}_nameSorted.bam',
          count_table = self._wk_dir/f'{per_batch}/{per_sample}/04.featureCounts/{per_sample}',
          count_summary_table = self._wk_dir / f'{per_batch}/{per_sample}/04.featureCounts/{per_sample}.summary',
        )
      if self._rule in ["s05_generate_count_details_nascent"]:
        return dict(
          bam_file=self._wk_dir/f'{per_batch}/{per_sample}/04.featureCounts_nascent/{per_sample}_nameSorted.bam',
          count_table = self._wk_dir/f'{per_batch}/{per_sample}/04.featureCounts_nascent/{per_sample}',
          count_summary_table = self._wk_dir / f'{per_batch}/{per_sample}/04.featureCounts_nascent/{per_sample}.summary',
        )
      if self._rule in ["s06_bin_segment"]:
        return dict(count_detail_tab = self._wk_dir/f'{per_batch}/{per_sample}/05.count/{per_sample}_count_detail.txt',)
      if self._rule in ["s06_bin_segment_nascent"]:
        return dict(count_detail_tab = self._wk_dir/f'{per_batch}/{per_sample}/05.count_nascent/{per_sample}_count_detail.txt',)
      if self._rule in ["s07_analysis"]:
        return dict(square_bin=self._wk_dir/f'{per_batch}/{per_sample}/06.binSegment/square_bin',)
      if self._rule in ["s06_bin_segment_with_image"]:
        return dict(
          image_file = self._info_db[per_batch]["samples"][per_sample]["image_file"],
          count_detail_tab = self._wk_dir/f'{per_batch}/{per_sample}/05.count/{per_sample}_count_detail.txt',
        )
      if self._rule in ["s07_analysis_with_image"]:
        return dict(square_bin=self._wk_dir/f'{per_batch}/{per_sample}/06.binSegment_image/square_bin',)

      # Nascent RNA pipeline.
      if self._rule in ["s04_estimate_mapping_quality"]:
        return dict(
          bam_file=self._wk_dir/f'{per_batch}/{per_sample}/03.star/{per_sample}_Aligned.sortedByCoord.out.bam',
          bai_file=self._wk_dir/f'{per_batch}/{per_sample}/03.star/{per_sample}_Aligned.sortedByCoord.out.bam.bai',
        )
      if self._rule in ["s04_remove_duplicates"]:
        return dict(
          bam_file=self._wk_dir/f'{per_batch}/{per_sample}/03.star/{per_sample}_Aligned.sortedByCoord.out.bam',
          bai_file=self._wk_dir/f'{per_batch}/{per_sample}/03.star/{per_sample}_Aligned.sortedByCoord.out.bam.bai',
        )
      if self._rule in ["s05_correct_base_qualities"]:
        return dict(
          bam_file=self._wk_dir/f'{per_batch}/{per_sample}/04.rmdup/{per_sample}.rmdup.out.bam',
          bai_file=self._wk_dir/f'{per_batch}/{per_sample}/04.rmdup/{per_sample}.rmdup.out.bam.bai',
        )
      if self._rule in ["s06_convert_bam_to_cit"]:
        return dict(
          bam_file=self._wk_dir/f'{per_batch}/{per_sample}/04.rmdup/{per_sample}.rmdup.base_corrected.out.bam',
          bai_file=self._wk_dir/f'{per_batch}/{per_sample}/04.rmdup/{per_sample}.rmdup.base_corrected.out.bam.bai',
        )
      # if self._rule in ["s06_merge_cits"]:
      #   input_files = []
      #   for per_batch, sample_info_dict in self._info_db.items():
      #     if per_batch != per_batch: continue
      #     for per_sample in sample_info_dict.get("samples").keys():
      #       input_files.append(self._wk_dir/f'{per_batch}/cit_files/per_sample/{per_sample}.corrected.cit')
      #   return input_files
      # if self._rule in ['s07_estimate_nascent_rna']:
      #   return dict(merged_cit=self._wk_dir/f'{per_batch}/cit_files/merged/{per_batch}.all_samples.cit')
      if self._rule in ["s07_merge_cits_per_treatment"]:
        ctrl_cits = []
        for cur_sample in self._info_db[per_batch]["parameters"]["control_samples"]:
          ctrl_cits.append(self._wk_dir/f'{per_batch}/{cur_sample}/06.cit_files/{cur_sample}.corrected.cit')
        return dict(
          ctrl_cits=ctrl_cits,
          sample_cit=self._wk_dir / f'{per_batch}/{per_sample}/06.cit_files/{per_sample}.corrected.cit'
        )
      if self._rule in ["s08_estimate_nascent_rna_per_treatment"]:
        return dict(merged_cit=self._wk_dir/f'{per_batch}/{per_sample}/06.cit_files/{per_sample}.corrected.with_ctrl.cit',)
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
config_file = Path("~/Documents/projects/wp_vasaseq/scripts/snakemake/slide_seq.configuration.yaml").expanduser()
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
bin_size = parameters.get("bin_size")

required_outputs = parameters.get("required_outputs", [])
rule_io = RuleIO(working_dir, required_outputs, sequencing_info)
final_outputs = rule_io()


wildcard_constraints:
  batch_id = r'\w+', sample_id = r'\w+', per_sample = r'\w+'


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


rule s00_create_barcodes:
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
    batch_info = sequencing_info[batch_id]
    with open(output.sample_barcode_file, "w") as out_1, open(output.valid_barcode_file, "w") as out_2:
      for sample_id, sample_info in batch_info.get("samples").items():
        barcode = sample_info["barcode"]
        out_1.write(f"{sample_id}\t{barcode}\n")
        out_2.write(f"{barcode}\n")


rule s00_normalize_fastq:
  input:
    unpack(RuleIO(working_dir, "s00_normalize_fastq", sequencing_info))
  output:
    fastq_r1 = working_dir / '{batch_id}/{sample_id}/00.fastqs/{sample_id}_1.fq.gz',
    fastq_r2 = working_dir / '{batch_id}/{sample_id}/00.fastqs/{sample_id}_2.fq.gz',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    soft_list = True,
  resources:
    n_cpus = 10
  shell:
    '''
    mkdir -p {params.out_dir}
    pigz -dc -p 5 {input.fastq_r1} | pigz -p 5 > {output.fastq_r1}
    pigz -dc -p 5 {input.fastq_r2} | pigz -p 5 > {output.fastq_r2}
    apptainer exec -B /mnt/backup {app_image} \
      fastqc -q -t {resources.n_cpus} -o {params.out_dir} {output.fastq_r1} {output.fastq_r2}
    '''


rule s01_decode_fastq:
  input: unpack(RuleIO(working_dir, "s01_decode_fastq", sequencing_info))
  output:
    stat_report = working_dir / '{batch_id}/{sample_id}/01.barcode/stat.txt',
    fastq_r1 = working_dir / '{batch_id}/{sample_id}/01.barcode/{sample_id}_1.fq.gz',
    fastq_r2 = working_dir / '{batch_id}/{sample_id}/01.barcode/{sample_id}_2.fq.gz',
  resources:
    n_cpus = 10
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    chemistry = lambda wc: get_parameters(wc.batch_id, "chemistry"),
    chemistry_pattern = lambda wc: get_parameters(wc.batch_id, "chemistry_pattern"),
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna barcode \
      --sample {wildcards.sample_id} \
      --chemistry {params.chemistry} \
      --whitelist {input.barcode_file} \
      --pattern {params.chemistry_pattern} \
      --thread {resources.n_cpus} \
      --mode strna \
      --lowNum 2 \
      --output_R1 \
      --resume \
      --gzip \
      --fq1 {input.fastq_r1} \
      --fq2 {input.fastq_r2} \
      --outdir {params.out_dir}
    '''


rule s02_trimming_fastq:
  input:
    unpack(RuleIO(working_dir, "s02_trimming_fastq", sequencing_info))
  output:
    log_out = working_dir / '{batch_id}/{sample_id}/02.cutadapt/cutadapt.log',
    stat_report = working_dir / '{batch_id}/{sample_id}/02.cutadapt/stat.txt',
    fastq_r2 = working_dir / '{batch_id}/{sample_id}/02.cutadapt/{sample_id}_clean_2.fq.gz',
  resources:
    n_cpus = 10
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    insert_size = lambda wc: parameters["insert_size"],
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna cutadapt \
      --sample {wildcards.sample_id} \
      --thread {resources.n_cpus} \
      --overlap 10 \
      --minimum_length 60 \
      --nextseq_trim 0 \
      --gzip \
      --cutadapt_param '-a AGATCGGAAGAG -a 'T{{5}}' -g AAGCAGTGGTATCAACGC -e 1 -n 10 --trim-n -q 20' \
      --insert {params.insert_size}  \
      --fq {input.fastq_r2} \
      --outdir {params.out_dir}

    echo "FASTQ quality after trimming ..." # Output FASTQ quality
    apptainer exec -B /mnt/backup {app_image} fastqc -q -t {resources.n_cpus} -o {params.out_dir} {output.fastq_r2}
    '''


rule s03_align_reads:
  input:
    unpack(RuleIO(working_dir, "s03_align_reads", sequencing_info))
  output:
    stat_report = working_dir / '{batch_id}/{sample_id}/03.star/stat.txt',
    log_out = working_dir / '{batch_id}/{sample_id}/03.star/{sample_id}_Log.final.out',
    bam_file = working_dir / '{batch_id}/{sample_id}/03.star/{sample_id}_Aligned.sortedByCoord.out.bam',
    bai_file = working_dir / '{batch_id}/{sample_id}/03.star/{sample_id}_Aligned.sortedByCoord.out.bam.bai',
    nascent_bam_file = working_dir / '{batch_id}/{sample_id}/03.star/{sample_id}_Aligned.sortedByCoord.out.nascent.bam',
    nascent_bai_file = working_dir / '{batch_id}/{sample_id}/03.star/{sample_id}_Aligned.sortedByCoord.out.nascent.bam.bai',
    old_bam_file = working_dir / '{batch_id}/{sample_id}/03.star/{sample_id}_Aligned.sortedByCoord.out.old.bam',
    old_bai_file = working_dir / '{batch_id}/{sample_id}/03.star/{sample_id}_Aligned.sortedByCoord.out.old.bam.bai',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    out_prefix = lambda wc: f'{working_dir}/{wc.batch_id}/{wc.sample_id}/03.star/{wc.sample_id}_',
  resources:
    n_cpus = 10
  shell:
    '''
    mkdir -p {params.out_dir}

    echo "Aligning reads ..."
    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna star \
      --sample {wildcards.sample_id} \
      --thread {resources.n_cpus} \
      --genomeDir {star_index_dir} \
      --outFilterMultimapNmax 1 \
      --STAR_param "--outSAMattributes NH HI AS NM MD nM --alignEndsType EndToEnd" \
      --starMem 30 \
      --fq {input.fastq_r2} \
      --outdir {params.out_dir}

    echo "Creating index ..."
    apptainer exec -B /mnt/backup {app_image} samtools index -b -@ {resources.n_cpus} {output.bam_file}

    echo "Extracting nascent reads ..."
    apptainer exec -B /mnt/backup {app_image} python /opt/tools/scripts/extract_new_reads.py \
      -@ {resources.n_cpus} \
      -n {output.nascent_bam_file} \
      -o {output.old_bam_file} \
      {output.bam_file}

    samtools index -@ {resources.n_cpus} {output.nascent_bam_file}
    samtools index -@ {resources.n_cpus} {output.old_bam_file}
    '''


rule s04_count_reads_nascent:
  input:
    unpack(RuleIO(working_dir, "s04_count_reads_nascent", sequencing_info))
  output:
    stat_report = working_dir / '{batch_id}/{sample_id}/04.featureCounts_nascent/stat.txt',
    bam_file = working_dir / '{batch_id}/{sample_id}/04.featureCounts_nascent/{sample_id}_nameSorted.bam',
    count_table = working_dir / '{batch_id}/{sample_id}/04.featureCounts_nascent/{sample_id}',
    count_summary_table = working_dir / '{batch_id}/{sample_id}/04.featureCounts_nascent/{sample_id}.summary',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    feature_type = lambda wc: parameters["feature_type"],
  resources:
    n_cpus = 10
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna featureCounts \
      --sample {wildcards.sample_id} \
      --thread {resources.n_cpus} \
      --gtf_type {params.feature_type} \
      --genomeDir {star_index_dir} \
      --featureCounts_param '-s 1 ' \
      --input {input.bam_file} \
      --outdir {params.out_dir}
    '''


rule s05_generate_count_details_nascent:
  input:
    unpack(RuleIO(working_dir, "s05_generate_count_details_nascent", sequencing_info))
  output:
    stat_report = working_dir / '{batch_id}/{sample_id}/05.count_nascent/stat.txt',
    count_table = working_dir / '{batch_id}/{sample_id}/05.count_nascent/{sample_id}_counts.txt',
    count_detail_tab = working_dir / '{batch_id}/{sample_id}/05.count_nascent/{sample_id}_count_detail.txt',
    raw_bc_mat = directory(working_dir / '{batch_id}/{sample_id}/05.count_nascent/{sample_id}_raw_feature_bc_matrix'),
    filtered_bc_mat = directory(working_dir / '{batch_id}/{sample_id}/05.count_nascent/{sample_id}_filtered_feature_bc_matrix'),
  params:
    expected_cell_num = lambda wc: get_parameters(wc.batch_id, "expected_cell_num"),
    out_dir = lambda _, output: Path(output[0]).parent
  resources:
    n_cpus = 10
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna count \
      --sample {wildcards.sample_id} \
      --thread {resources.n_cpus} \
      --genomeDir {star_index_dir} \
      --expected_cell_num {params.expected_cell_num} \
      --cell_calling_method auto \
      --force_cell_num None \
      --bam {input.bam_file} \
      --outdir {params.out_dir}
    '''


rule s06_bin_segment_nascent:
  input:
    unpack(RuleIO(working_dir, "s06_bin_segment_nascent", sequencing_info))
  output:
    stat_report = working_dir / '{batch_id}/{sample_id}/06.binSegment_nascent/stat.txt',
    bin_statistics = working_dir / '{batch_id}/{sample_id}/06.binSegment_nascent/bin_statistics.html',
    bc_tissue_positions = working_dir / '{batch_id}/{sample_id}/06.binSegment_nascent/{sample_id}_Barcodes_tissue_positions.csv',
    images = directory(working_dir / '{batch_id}/{sample_id}/06.binSegment_nascent/images'),
    square_bin = directory(working_dir / '{batch_id}/{sample_id}/06.binSegment_nascent/square_bin'),
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    pixel_size = lambda wc: get_parameters(wc.batch_id, "pixel_size"),
    model_file = lambda wc: get_parameters(wc.batch_id, "model_file"),
    segment_method = lambda wc: get_parameters(wc.batch_id, "segment_method"),
    chip_data_dir = lambda wc: Path(sequencing_info[wc.batch_id]["samples"][wc.sample_id]["barcode_file"]).parent,
  resources:
    n_cpus = 10
  priority: 100
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna binSegment \
      --sample {wildcards.sample_id} \
      --thread {resources.n_cpus} \
      --genomeDir {star_index_dir} \
      --model {params.model_file} \
      --segment \
      --pixel-size {params.pixel_size} \
      --input {params.chip_data_dir} \
      --method {params.segment_method} \
      --count \
      --count_detail {input.count_detail_tab} \
      --outdir {params.out_dir}
    '''


rule s04_count_reads:
  input:
    unpack(RuleIO(working_dir, "s04_count_reads", sequencing_info))
  output:
    stat_report = working_dir / '{batch_id}/{sample_id}/04.featureCounts/stat.txt',
    bam_file = working_dir / '{batch_id}/{sample_id}/04.featureCounts/{sample_id}_nameSorted.bam',
    count_table = working_dir / '{batch_id}/{sample_id}/04.featureCounts/{sample_id}',
    count_summary_table = working_dir / '{batch_id}/{sample_id}/04.featureCounts/{sample_id}.summary',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    feature_type = lambda wc: parameters["feature_type"],
  resources:
    n_cpus = 10
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna featureCounts \
      --sample {wildcards.sample_id} \
      --thread {resources.n_cpus} \
      --gtf_type {params.feature_type} \
      --genomeDir {star_index_dir} \
      --featureCounts_param '-s 1 ' \
      --input {input.bam_file} \
      --outdir {params.out_dir}
    '''


rule s05_generate_count_details:
  input:
    unpack(RuleIO(working_dir, "s05_generate_count_details", sequencing_info))
  output:
    stat_report = working_dir / '{batch_id}/{sample_id}/05.count/stat.txt',
    count_table = working_dir / '{batch_id}/{sample_id}/05.count/{sample_id}_counts.txt',
    count_detail_tab = working_dir / '{batch_id}/{sample_id}/05.count/{sample_id}_count_detail.txt',
    raw_bc_mat = directory(working_dir / '{batch_id}/{sample_id}/05.count/{sample_id}_raw_feature_bc_matrix'),
    filtered_bc_mat = directory(working_dir / '{batch_id}/{sample_id}/05.count/{sample_id}_filtered_feature_bc_matrix'),
  params:
    expected_cell_num = lambda wc: get_parameters(wc.batch_id, "expected_cell_num"),
    out_dir = lambda _, output: Path(output[0]).parent
  resources:
    n_cpus = 10
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna count \
      --sample {wildcards.sample_id} \
      --thread {resources.n_cpus} \
      --genomeDir {star_index_dir} \
      --expected_cell_num {params.expected_cell_num} \
      --cell_calling_method auto \
      --force_cell_num None \
      --bam {input.bam_file} \
      --outdir {params.out_dir}
    '''


rule s06_bin_segment:
  input:
    unpack(RuleIO(working_dir, "s06_bin_segment", sequencing_info))
  output:
    stat_report = working_dir / '{batch_id}/{sample_id}/06.binSegment/stat.txt',
    bin_statistics = working_dir / '{batch_id}/{sample_id}/06.binSegment/bin_statistics.html',
    bc_tissue_positions = working_dir / '{batch_id}/{sample_id}/06.binSegment/{sample_id}_Barcodes_tissue_positions.csv',
    images = directory(working_dir / '{batch_id}/{sample_id}/06.binSegment/images'),
    square_bin = directory(working_dir / '{batch_id}/{sample_id}/06.binSegment/square_bin'),
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    pixel_size = lambda wc: get_parameters(wc.batch_id, "pixel_size"),
    model_file = lambda wc: get_parameters(wc.batch_id, "model_file"),
    segment_method = lambda wc: get_parameters(wc.batch_id, "segment_method"),
    chip_data_dir = lambda wc: Path(sequencing_info[wc.batch_id]["samples"][wc.sample_id]["barcode_file"]).parent,
  resources:
    n_cpus = 10
  priority: 100
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna binSegment \
      --sample {wildcards.sample_id} \
      --thread {resources.n_cpus} \
      --genomeDir {star_index_dir} \
      --model {params.model_file} \
      --segment \
      --pixel-size {params.pixel_size} \
      --input {params.chip_data_dir} \
      --method {params.segment_method} \
      --count \
      --count_detail {input.count_detail_tab} \
      --outdir {params.out_dir}
    '''


rule s07_analysis:
  input:
    unpack(RuleIO(working_dir, "s07_analysis", sequencing_info)),
  output:
    stat_report = working_dir / '{batch_id}/{sample_id}' / f'07.analysis_bin{bin_size}/stat.txt',
    markers = working_dir / '{batch_id}/{sample_id}' / f'07.analysis_bin{bin_size}' / '{sample_id}_markers.tsv',
    markers_raw = working_dir / '{batch_id}/{sample_id}' / f'07.analysis_bin{bin_size}' / '{sample_id}_markers_raw.tsv',
    tsne_coord = working_dir / '{batch_id}/{sample_id}' / f'07.analysis_bin{bin_size}' / '{sample_id}_tsne_coord.tsv',
    bin10_h5ad = working_dir / '{batch_id}/{sample_id}' / f'07.analysis_bin{bin_size}' / '{sample_id}_bin10.h5ad',
    bin20_h5ad = working_dir / '{batch_id}/{sample_id}' / f'07.analysis_bin{bin_size}' / '{sample_id}_bin20.h5ad',
    bin50_h5ad = working_dir / '{batch_id}/{sample_id}' / f'07.analysis_bin{bin_size}' / '{sample_id}_bin50.h5ad',
    bin100_h5ad = working_dir / '{batch_id}/{sample_id}' / f'07.analysis_bin{bin_size}' / '{sample_id}_bin100.h5ad',
    bioinfo_data = directory(working_dir / '{batch_id}/{sample_id}' / f'07.analysis_bin{bin_size}/Bioinfodata'),
    report = working_dir / '{batch_id}/{sample_id}' / f'07.analysis_bin{bin_size}' / '{sample_id}_report.html',
  resources:
    n_cpus = 10
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    pixel_size = lambda wc: get_parameters(wc.batch_id, "pixel_size"),
  priority: 100
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna analysis \
      --sample {wildcards.sample_id} \
      --thread {resources.n_cpus} \
      --genomeDir {star_index_dir} \
      --square_bin_dir {input.square_bin} \
      --pixel-size {params.pixel_size} \
      --bin {bin_size} \
      --outdir {params.out_dir}
    cp {working_dir}/{wildcards.batch_id}/{wildcards.sample_id}/{wildcards.sample_id}_report.html {params.out_dir}
    '''


rule s06_bin_segment_with_image:
  input:
    unpack(RuleIO(working_dir, "s06_bin_segment_with_image", sequencing_info))
  output:
    stat_report = working_dir / '{batch_id}/{sample_id}/06.binSegment_image/stat.txt',
    bin_statistics = working_dir / '{batch_id}/{sample_id}/06.binSegment_image/bin_statistics.html',
    bc_tissue_positions = working_dir / '{batch_id}/{sample_id}/06.binSegment_image/{sample_id}_Barcodes_tissue_positions.csv',
    images = directory(working_dir / '{batch_id}/{sample_id}/06.binSegment_image/images'),
    square_bin = directory(working_dir / '{batch_id}/{sample_id}/06.binSegment_image/square_bin'),
  resources:
    n_cpus = 10
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    pixel_size = lambda wc: get_parameters(wc.batch_id, "pixel_size"),
    model_file = lambda wc: get_parameters(wc.batch_id, "model_file"),
    segment_method = lambda wc: get_parameters(wc.batch_id, "segment_method"),
    chip_data_dir = lambda wc: Path(sequencing_info[wc.batch_id]["samples"][wc.sample_id]["barcode_file"]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna binSegment \
      --sample {wildcards.sample_id} \
      --thread {resources.n_cpus} \
      --genomeDir {star_index_dir} \
      --model {params.model_file} \
      --segment \
      --pixel-size {params.pixel_size} \
      --input {params.chip_data_dir} \
      --method {params.segment_method} \
      --tif {input.image_file} \
      --bs_out {params.out_dir}/bs_out \
      --count \
      --count_detail {input.count_detail_tab} \
      --outdir {params.out_dir} # {{sampledir}}/06.binSegment_image
    '''


rule s07_analysis_with_image:
  input:
    unpack(RuleIO(working_dir, "s07_analysis_with_image", sequencing_info)),
  output:
    stat_report = working_dir / '{batch_id}/{sample_id}' / '07.analysis_bin{bin_size}_with_image/stat.txt',
    markers = working_dir / '{batch_id}/{sample_id}' / '07.analysis_bin{bin_size}_with_image' / '{sample_id}_markers.tsv',
    markers_raw = working_dir / '{batch_id}/{sample_id}' / '07.analysis_bin{bin_size}_with_image' / '{sample_id}_markers_raw.tsv',
    tsne_coord = working_dir / '{batch_id}/{sample_id}' / '07.analysis_bin{bin_size}_with_image' / '{sample_id}_tsne_coord.tsv',
    bin10_h5ad = working_dir / '{batch_id}/{sample_id}' / '07.analysis_bin{bin_size}_with_image' / '{sample_id}_bin10.h5ad',
    bin20_h5ad = working_dir / '{batch_id}/{sample_id}' / '07.analysis_bin{bin_size}_with_image' / '{sample_id}_bin20.h5ad',
    bin50_h5ad = working_dir / '{batch_id}/{sample_id}' / '07.analysis_bin{bin_size}_with_image' / '{sample_id}_bin50.h5ad',
    bin100_h5ad = working_dir / '{batch_id}/{sample_id}' / '07.analysis_bin{bin_size}_with_image' / '{sample_id}_bin100.h5ad',
    bioinfo_data = directory(working_dir / '{batch_id}/{sample_id}' / '07.analysis_bin{bin_size}_with_image/Bioinfodata'),
  resources:
    n_cpus = 10
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    pixel_size = lambda wc: get_parameters(wc.batch_id, "pixel_size"),
  priority: 100
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} celatlas_spatial rna analysis \
      --sample {wildcards.sample_id} \
      --thread {resources.n_cpus} \
      --genomeDir {star_index_dir} \
      --square_bin_dir {input.square_bin} \
      --pixel-size {params.pixel_size} \
      --bin {bin_size} \
      --outdir {params.out_dir}
    cp {working_dir}/{wildcards.batch_id}/{wildcards.sample_id}/{wildcards.sample_id}_report.html {params.out_dir}
    '''


rule s04_estimate_mapping_quality:
  input:
    unpack(RuleIO(working_dir, "s04_estimate_mapping_quality", sequencing_info)),
    genomic_feature_file = parameters["genomic_feature_file"],
    genome_sequence_file = genome_sequence_file,
  output:
    qorts_done = working_dir / '{batch_id}/{sample_id}/04.mapping_reports/{sample_id}_qorts/qorts.done',
    qorts_dir = directory(working_dir / '{batch_id}/{sample_id}/04.mapping_reports/{sample_id}_qorts'),
    qualimap_done = working_dir / '{batch_id}/{sample_id}/04.mapping_reports/{sample_id}_qualimap/qualimap.done',
    qualimap_dir = directory(working_dir / '{batch_id}/{sample_id}/04.mapping_reports/{sample_id}_qualimap'),
  params:
    out_dir = lambda _, output: Path(output[0]).parent
  resources:
    mem = "16G"
  shell:
    '''
    mkdir -p {output.qualimap_dir} {output.qorts_dir}
    apptainer exec -B /mnt/backup {app_image} qualimap --java-mem-size={resources.mem} rnaseq \
      --sorted -bam {input.bam_file} -gtf {input.genomic_feature_file} -outdir {output.qualimap_dir}
    touch {output.qualimap_done}

    apptainer exec -B /mnt/backup {app_image} qorts QC \
      --genomeFA {input.genome_sequence_file} \
      --addFunctions mismatchEngine,cigarMatch,referenceMatch,writeDocs,calcDetailedGeneCounts \
      --maxReadLength 250 --singleEnded --generatePlots \
      {input.bam_file} {input.genomic_feature_file} {output.qorts_dir}
    touch {output.qorts_done}
    '''


rule s04_remove_duplicates:
  input:
    unpack(RuleIO(working_dir, "s04_remove_duplicates", sequencing_info))
  output:
    bam_file = working_dir / '{batch_id}/{sample_id}/04.rmdup/{sample_id}.rmdup.out.bam',
    bai_file = working_dir / '{batch_id}/{sample_id}/04.rmdup/{sample_id}.rmdup.out.bam.bai',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  resources:
    n_cpus = 10
  shell:
    '''
    mkdir -p {params.out_dir}

    awk -F$'\t' -f- <<'EOF' <(apptainer exec -B /mnt/backup {app_image} samtools view -h {input.bam_file}) \
    | apptainer exec -B /mnt/backup {app_image} samtools view -O BAM -o {output.bam_file}.tmp.bam
$0 ~ /^@/ {{ print; next }}
{{ split($1, read_id, "_"); OFS = "\\t"; print $0, "CR:Z:"read_id[1], "UR:Z:"read_id[2] }}
EOF
    apptainer exec -B /mnt/backup {app_image} samtools index -@ {resources.n_cpus} {output.bam_file}.tmp.bam
    apptainer exec -B /mnt/backup {app_image} umi_tools dedup \
      -I {output.bam_file}.tmp.bam --extract-umi-method tag --umi-tag UR --cell-tag CR -S {output.bam_file}
    apptainer exec -B /mnt/backup {app_image} samtools index {output.bam_file}

    rm -f {output.bam_file}.tmp.bam {output.bam_file}.tmp.bam.bai
    '''


rule s05_correct_base_qualities:
  input:
    unpack(RuleIO(working_dir, "s05_correct_base_qualities", sequencing_info))
  output:
    bam_file = working_dir / '{batch_id}/{sample_id}/04.rmdup/{sample_id}.rmdup.base_corrected.out.bam',
    bai_file = working_dir / '{batch_id}/{sample_id}/04.rmdup/{sample_id}.rmdup.base_corrected.out.bam.bai',
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


rule s06_convert_bam_to_cit:
  input:
    unpack(RuleIO(working_dir, "s06_convert_bam_to_cit", sequencing_info)),
  output:
    cit_file = working_dir / '{batch_id}/{sample_id}/06.cit_files/{sample_id}.cit',
    cit_corr_file = working_dir / '{batch_id}/{sample_id}/06.cit_files/{sample_id}.corrected.cit',
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
    '''


rule s07_merge_cits_per_treatment:
  input:
    unpack(RuleIO(working_dir, "s07_merge_cits_per_treatment", sequencing_info))
  output:
    merged_cit = working_dir / '{batch_id}/{sample_id}/06.cit_files/{sample_id}.corrected.with_ctrl.cit',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
    external_controls = lambda wc: get_parameters(wc.batch_id, "external_controls"),
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} \
      gedi -e MergeCIT -p {output.merged_cit} {input.ctrl_cits} {input.sample_cit} {params.external_controls}
    apptainer exec -B /mnt/backup {app_image} gedi -e ReadCount -p {output.merged_cit} 
    '''


rule s08_estimate_nascent_rna_per_treatment:
  input:
    unpack(RuleIO(working_dir, "s08_estimate_nascent_rna_per_treatment", sequencing_info))
  output:
    nascent_rna = working_dir / '{batch_id}/{sample_id}/07.nascent_rna/{sample_id}.tsv.gz',
  params:
    out_prefix = lambda _, output: output.nascent_rna.strip(".tsv.gz"),
    out_dir = lambda _, output: Path(output[0]).parent,
    ctrl_pattern = lambda wc: get_parameters(wc.batch_id, "control_pattern")
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} gedi -e Slam \
      -genomic GEDI_GRCm38 -reads {input.merged_cit} -prefix {params.out_prefix} \
      -minEstimateReads 1000 -trim3p 30 -trim5p 20 -progress -plot -D -full -no4sUpattern {params.ctrl_pattern}
    '''


#
## Estimating nascent RNA per bin.
#
rule s07_collect_barcode_per_bin:
  input:
    unpack(RuleIO(working_dir, "s07_collect_barcode_per_bin", sequencing_info))
  output:
    bin_code = working_dir / '{batch_id}/{sample_id}/07.split_bam_by_bin/bin_code/bin_{bin_code}.txt',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} collect_barcode_per_bin.py \
      {input.raw_bc_file} {input.binned_bc_file}
    '''


rule s08_split_bam_by_bin:
  input:
    unpack(RuleIO(working_dir, "s08_split_bam_by_bin", sequencing_info)),
    bin_code = working_dir / '{batch_id}/{sample_id}/07.split_bam_by_bin/bin_code/bin_{bin_code}.txt',
    # bam_file = working_dir / '{batch_id}/{sample_id}/04.rmdup/{sample_id}.rmdup.out.bam',
    # bai_file = working_dir / '{batch_id}/{sample_id}/04.rmdup/{sample_id}.rmdup.out.bam.bai',
    # square_bin = directory(working_dir / '{batch_id}/{sample_id}/06.binSegment/square_bin'),
  output:
    bam_file = working_dir / '{batch_id}/{sample_id}/07.split_bam_by_bin/{per_sample}.bin_{bin_code}.bam',
    bai_file = working_dir / '{batch_id}/{sample_id}/07.split_bam_by_bin/{per_sample}.bin_{bin_code}.bam.bai',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}

    apptainer exec -B /mnt/backup {app_image} samtools view \
      -d CR:{params.sample_barcode} -O BAM -o {output.bam_file} -@ {resources.n_cpus} {input.bam_file}

    apptainer exec -B /mnt/backup {app_image} samtools view -H {input.bam_file} -D  > {output.bam_file}
    apptainer exec -B /mnt/backup {app_image} gedi -e SplitBAM -p {output.bam_file} {input.bam_file}
    '''


rule s09_convert_bam_to_cit:
  input:
    unpack(RuleIO(working_dir, "s07_split_bam_by_bin", sequencing_info))
    # bam_file = working_dir / '{batch_id}/{sample_id}/07.split_bam_by_bin/{per_sample}.bin_{bin_code}.bam',
    # bai_file = working_dir / '{batch_id}/{sample_id}/07.split_bam_by_bin/{per_sample}.bin_{bin_code}.bam.bai',
  output:
    cit_file = working_dir / '{batch_id}/{sample_id}/09.cit_files/{per_sample}.bin_{bin_code}.cit',
    cit_corr_file = working_dir / '{batch_id}/{sample_id}/09.cit_files/{per_sample}.bin_{bin_code}.corrected.cit',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    mkdir -p {params.out_dir}
    apptainer exec -B /mnt/backup {app_image} gedi -e Bam2CIT -p {output.cit_file} {input.bam_file}
    apptainer exec -B /mnt/backup {app_image} gedi -e CorrectCIT -p {output.cit_corr_file} {output.cit_file}
    '''


rule s10_merge_cits_per_bin:
  input:
    unpack(RuleIO(working_dir, "s09_convert_bam_to_cit", sequencing_info))
  output:
    merged_cit = working_dir / '{batch_id}/{sample_id}/10.cit_files/merged/{per_sample}.all_samples.cit',
  params:
    out_dir = lambda _, output: Path(output[0]).parent,
  shell:
    '''
    '''


# rule s06_merge_cits:
#   input:
#     unpack(RuleIO(working_dir, "s06_merge_cits", sequencing_info))
#   output:
#     merged_cit = working_dir / '{batch_id}/06.cit_files/merged/{batch_id}.all_samples.cit',
#   params:
#     out_dir = lambda _, output: Path(output[0]).parent,
#   shell:
#     '''
#     mkdir -p {params.out_dir}
#     apptainer exec -B /mnt/backup {app_image} gedi -e MergeCIT -p {output.merged_cit} {input}
#     apptainer exec -B /mnt/backup {app_image} gedi -e ReadCount -p {output.merged_cit} 
#     '''


# rule s07_estimate_nascent_rna:
#   input:
#     unpack(RuleIO(working_dir, "s07_estimate_nascent_rna", sequencing_info))
#   output:
#     nascent_rna = working_dir / '{batch_id}/07.nascent_rna/all_samples/{batch_id}.tsv.gz',
#   params:
#     # out_prefix = lambda wc, output: Path(output[0]).parent / f'{wc.batch_id}.nascent_rna',
#     out_prefix = lambda _, output: output.nascent_rna.strip(".tsv.gz"),
#     out_dir = lambda _, output: Path(output[0]).parent,
#     ctrl_pattern = lambda wc: sequencing_info.get(wc.batch_id, {}).get("control_pattern", "4sU")
#   shell:
#     '''
#     mkdir -p {params.out_dir}
#     apptainer exec -B /mnt/backup {app_image} gedi -e Slam \
#       -genomic GEDI_GRChm38 -reads {input.merged_cit} -prefix {params.out_prefix} \
#       -minEstimateReads 5000 -trim3p 30 -trim5p 20 -progress -plot -D -full -allGenes -no4sUpattern {params.ctrl_pattern}
#     '''

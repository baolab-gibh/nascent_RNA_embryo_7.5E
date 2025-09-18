#!/usr/bin/env bash
# File: alternative_splicing.sh
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Oct 25, 2024
# Updated:

project_dir=~/Documents/projects/wp_vasaseq


if [[ -e ${project_dir}/scripts/.env ]]; then source ${project_dir}/scripts/.env/bin/activate; fi

genomic_feature=${project_dir}/outputs/references/genome/mus_musculus.90.gtf
output_dir=${project_dir}/outputs/analysis/alternative_splicing/spladder_out
bam_dir=${project_dir}/outputs/analysis/alternative_splicing/bam_files
bam_files=$(command ls ${bam_dir}/*.bam | xargs -n 10000 | tr ' ' ',')
plot_dir=${project_dir}/outputs/analysis/alternative_splicing/trackplot

# Graph per sample
for per_bam in ${bam_dir}/*.bam; do
  sample_id=$(basename ${per_bam} | sed -e 's/.bam//g')
  if [[ -e ${output_dir}/spladder/genes_graph_conf3.${sample_id}.pickle ]]; then echo $per_bam && continue; fi
  spladder build --parallel 4 -o ${output_dir} -a ${genomic_feature} -b ${per_bam} --merge-strat single --no-extract-ase &
  while [[ $(jobs | wc -l) -ge 10 ]]; do sleep 10; done
done

# Merge graphs
alignments_list=${project_dir}/outputs/analysis/alternative_splicing/alignments_list.txt
spladder build -o ${output_dir} -a ${genomic_feature} -b ${alignments_list} --parallel 4 --merge-strat merge_graphs --no-extract-ase

# Quantifications
for per_bam in ${bam_dir}/*.bam; do
  sample_id=$(basename ${per_bam} | sed -e 's/.bam//g')
  if [[ -e ${output_dir}/spladder/genes_graph_conf3.${sample_id}.pickle ]]; then echo $per_bam && continue; fi
  spladder build -o ${output_dir} -a ${genomic_feature} -b ${per_bam} --merge-strat merge_graphs --no-extract-ase --quantify-graph --qmode single
  while [[ $(jobs | wc -l) -ge 10 ]]; do sleep 10; done
done

# Collect the individual quantifications and aggregate them in a joint database
spladder build -o ${output_dir} -a ${genomic_feature} -b ${alignments_list} --merge-strat merge_graphs --no-extract-ase --quantify-graph --qmode collect

# Call alternative splicing events
spladder build -o ${output_dir} -a ${genomic_feature} -b ${alignments_list} --event-types exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons

# Check example AS
genomic_feature=${project_dir}/outputs/references/genome/mus_musculus.90.sorted.gtf.gz
alignments_list=${project_dir}/outputs/analysis/alternative_splicing/alignments_list.trackplot.txt
event_dir=${project_dir}/outputs/analysis/alternative_splicing/spladder_out


# Prepare trackplot
function def_batches() {
  echo ${3}:$(grep -E "[0-9]+${1}_[ATCG]+$" ${2} | grep -w ${3} | cut -f1 -d$'\t' | xargs -n 1000 | tr ' ' ',')
}

declare -A selected_features
# selected_features=( [Wnt3]=Wnt3 [Lefty2]=Lefty2 [Rhoa]=Rhoa [Bmp4]=Bmp4 [Gas5]=Gas5 [Dll1]=MP )
selected_features=( [Gas5]=P )


# plots
for gene_id in ${!selected_features[*]}; do
  target_region=${selected_features[${gene_id}]}
  # BAM files per batch from given ${target_region}
  batch_1=$(def_batches ${target_region} ${alignments_list} 240409_Lib_embryo)
  batch_2=$(def_batches ${target_region} ${alignments_list} 240612_Lib_28region)
  batch_3=$(def_batches ${target_region} ${alignments_list} 240620_Lib_38region)
  batch_4=$(def_batches ${target_region} ${alignments_list} 240703_Lib_32region)
  batch_5=$(def_batches ${target_region} ${alignments_list} 240710_Lib_37region)
  batch_6=$(def_batches ${target_region} ${alignments_list} 240717_Lib_28region)

  feature_region=$(zgrep \"${gene_id}\"\; ${genomic_feature} | awk '$3 == "gene" {print $1":"$4"-"$5}' | head -n 1)
  echo "[I]: Working on ${gene_id} ${feature_region} ..."

  out_file=${plot_dir}/trackplot.${gene_id}.${target_region}.pdf
  trackplot -e ${feature_region} -r ${genomic_feature} -o ${out_file} --density <(grep -E "[0-9]+${target_region}_[TCGA]+" ${alignments_list} | sort) \
    --dpi 300 --width 10 --height 1 --intron-scale .2 --threshold 3 \
    --normalize-format count \
    --show-junction-num 1>&2 2>/dev/null

  feature_id=$(zgrep \"${gene_id}\"\; ${genomic_feature} | grep -Eom1 "ENSMUSG[0-9]+")
  has_event=$(zcat ${event_dir}/*.confirmed.txt.gz | grep -c ${feature_id})
  if [[ ${has_event} -eq 0 ]]; then
    echo "[W]: No event for ${gene_id}" && continue
  else
    spladder viz -o ${output_dir} -O ${target_region}.${gene_id} -f pdf \
      --track event any \
      --track splicegraph \
      --range gene ${feature_id} \
      --track coverage ${batch_1} \
      --track coverage ${batch_2} \
      --track coverage ${batch_3} \
      --track coverage ${batch_4} \
      --track coverage ${batch_5} \
      --track coverage ${batch_6}
  fi
done


# Call alternative splicing events per layers (L0, Ectoderm, endoderm, mesoderm)

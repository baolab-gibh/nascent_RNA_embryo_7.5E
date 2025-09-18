#!/bin/bash
# File: workflow.sh
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Jun 03, 2024
# Updated: Jun 17, 2024

n_colors=$(tput colors)
if [[ -n ${n_colors} ]]; then
  T_BOLD=$(tput bold)
  T_UNDERLINE=$(tput smul)
  T_RESET=$(tput sgr0)
  T_RED=$(tput setaf 1)
  T_GREEN=$(tput setaf 2)
  T_YELLOW=$(tput setaf 3)
  T_BLACK=$(tput setaf 0)
fi

msg() {
  echo -e "[I]: $@"
}

warn() {
  echo -e "${T_YELLOW}${T_BOLD}[W]: $@${T_RESET}"
}

error() {
  echo -e "${T_RED}${T_BOLD}[E]: $@${T_RESET}"; return -1
}

# trap 'error Errors at $LINENO' ERR


project_dir=~/Documents/projects/wp_vasaseq

if [[ ${USER} == 'zhzhang_gibh' ]]; then
  source ~/tools/miniconda3/bin/activate
  mamba activate sm732 # Snakemake v7.32.4
elif [[ ${USER} =~ 'zzhang' && -d ${project_dir}/scripts/.env ]]; then
  source ${project_dir}/scripts/.env/bin/activate
fi

cd ${project_dir}/temps || return 1

n_jobs=1
prof_dir=${project_dir}/scripts/snakemake/configs


# Configurations
# GEO, VASA, SLAM seq pipeline
work_dir=${project_dir}/outputs/analysis
smk_file=${project_dir}/scripts/snakemake/vasalam_ppl.smk
cfg_file=${project_dir}/scripts/snakemake/vasalam_ppl.yml

# Job graph, not able to show the DAG of modules
snakemake --forceall -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --config myconfigfile=${cfg_file} --rulegraph \
  | dot -Tpdf \
  >| ${work_dir}/preprocessing/geo_vasa_slam/logs/main_ppl.job_graph.pdf

# Touch existing files
snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --config myconfigfile=${cfg_file} --touch 2>/dev/null

# Dry-run
log_file=${work_dir}/preprocessing/geo_vasa_slam/logs/main_ppl.dry_run.txt
snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --config myconfigfile=${cfg_file} --dry-run \
  >| ${log_file}
awk '/^Job stats:$/ {count++} count >= 2 {print}' ${log_file}
less -X -S -N ${log_file}

# Unlock current .snakemake folder
snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --config myconfigfile=${cfg_file} --unlock

# Run
if [[ ${USER} == 'zhzhang_gibh' ]]; then
  snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --config myconfigfile=${cfg_file}
elif [[ ${USER} == 'zzhang' ]]; then
  snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --config myconfigfile=${cfg_file}
fi


# slide-seq, encoded mouse embryo
task_id=slide_encoded
task_id=geo_vasa_slam
task_id=geo_slam
task_id=geo_slam_dm
task_id=slide_decoded
case ${task_id} in
  # slide_decoded)
  #   work_dir=${project_dir}/outputs/analysis/preprocessing/slide_seq_decoded
  #   smk_file=${project_dir}/scripts/snakemake/celatlas_spatial.pipeline.smk
  #   break ;;
  # geo_slam)
  #   work_dir=${project_dir}/outputs/analysis/preprocessing/geo_slam
  #   smk_file=${project_dir}/scripts/snakemake/geo_slam_ppl.smk
  #   break ;;
  slide_decoded)
    work_dir=${project_dir}/outputs/analysis/preprocessing/slide_seq_decoded
    smk_file=${project_dir}/scripts/snakemake/slide_seq.pipeline.smk
    ;;
  slide_encoded)
    work_dir=${project_dir}/outputs/analysis/preprocessing/slide_seq
    smk_file=${project_dir}/scripts/snakemake/slide_seq.smk
    ;;
  peng_2019_nature)
    work_dir=${project_dir}/outputs/analysis/preprocessing/2019_peng_etal_nature
    smk_file=${project_dir}/scripts/snakemake/2019_peng_etal_nature.smk
    ;;
  slam_seq)
    work_dir=${project_dir}/outputs/analysis/preprocessing/slam_seq
    smk_file=${project_dir}/scripts/snakemake/slam_seq.smk
    ;;
  geo_vasa_slam | geo_slam)
    work_dir=${project_dir}/outputs/analysis/preprocessing/geo_vasa_slam
    smk_file=${project_dir}/scripts/snakemake/geo_slam.pipeline.smk
    ;;
  geo_slam_dm)
    work_dir=${project_dir}/outputs/analysis/preprocessing/geo_vasa_slam
    smk_file=${project_dir}/scripts/snakemake/geo_slam.pipeline.v2.smk
    ;;
  *)
    echo "Unknown task_id: ${task_id}"
    return 1
    ;;
esac
msg "Using $(basename ${smk_file})"

if [[ ! -d ${work_dir}/logs ]]; then mkdir -p ${work_dir}/logs/; fi
if [[ $? -eq 0 ]]; then
  # Job graph, not able to show the DAG of modules
  snakemake --forceall -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --rulegraph \
    | dot -Tpdf \
    >| ${work_dir}/logs/${task_id}.job_graph.pdf

  # Touch existing files
  snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --touch

  # Unlock current .snakemake folder
  snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --unlock

  # Dry-run
  dryrun_log=${work_dir}/logs/${task_id}.dryrun.log
  snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --dry-run >| ${dryrun_log}
  if [[ $? == 0 ]]; then awk '/^Job stats:$/ {count++} count >= 2 {print}' ${dryrun_log}; else tail -40 ${dryrun_log}; fi
  # cat ${dryrun_log} | head -50
  # less -X -N ${dryrun_log}

  # Run
  run_log=${work_dir}/logs/${task_id}.run.log
  if [[ ${USER} == 'zhzhang_gibh' ]]; then
    snakemake -kpd ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} |& tee ${run_log}
  elif [[ ${USER} == 'zzhang' ]]; then
    snakemake -kpd ${work_dir} -j ${n_jobs} -s ${smk_file} |& tee ${run_log}
  fi
fi

url=https://document-share.tos-cn-beijing.volces.com/v2/index.html?token=aHR0cHM6Ly9rZWZ1LXN6dG9zLWxpbmsudG9zLWNuLWd1YW5nemhvdS52b2xjZXMuY29tLz9YLVRvcy1BbGdvcml0aG09VE9TNC1ITUFDLVNIQTI1NiZYLVRvcy1DcmVkZW50aWFsPUFLTFRNMlF3T0RZelptVmxOelEyTkRReU9UazBaRFkwWkRSaVlUZGlPRFJsWldFJTJGMjAyNTA2MjElMkZjbi1ndWFuZ3pob3UlMkZ0b3MlMkZyZXF1ZXN0JlgtVG9zLURhdGU9MjAyNTA2MjFUMTExOTI1WiZYLVRvcy1FeHBpcmVzPTEyOTYwMDAmWC1Ub3MtUG9saWN5PWV5SmpiMjVrYVhScGIyNXpJanBiV3lKemRHRnlkSE10ZDJsMGFDSXNJaVJyWlhraUxDSXlNREkxTURZeU1WOWFhRzl1WjB0bFdXRnVSM1ZoYm1kYWFHOTFTbWxoYmt0aGJtZFpZVzVLYVhWWmRXRnVYekV2SWwwc2V5SmlkV05yWlhRaU9pSnJaV1oxTFhONmRHOXpMV3hwYm1zaWZWMTkmWC1Ub3MtU2lnbmF0dXJlPTM3NDc5ODljMDU5ZTNjNTNkODhjNzJjZGI1M2IxZjFhNDM2OTE3NjFjODUzM2VlM2I5Y2FhNTM1YmYyYjc0Y2I=
url=https://document-share.tos-cn-beijing.volces.com/v2/index.html?token=aHR0cHM6Ly9rZWZ1LXN6dG9zLWxpbmsudG9zLWNuLWd1YW5nemhvdS52b2xjZXMuY29tLz9YLVRvcy1BbGdvcml0aG09VE9TNC1ITUFDLVNIQTI1NiZYLVRvcy1DcmVkZW50aWFsPUFLTFRNMlF3T0RZelptVmxOelEyTkRReU9UazBaRFkwWkRSaVlUZGlPRFJsWldFJTJGMjAyNTA3MjMlMkZjbi1ndWFuZ3pob3UlMkZ0b3MlMkZyZXF1ZXN0JlgtVG9zLURhdGU9MjAyNTA3MjNUMTI1NDU5WiZYLVRvcy1FeHBpcmVzPTEyOTYwMDAmWC1Ub3MtUG9saWN5PWV5SmpiMjVrYVhScGIyNXpJanBiV3lKemRHRnlkSE10ZDJsMGFDSXNJaVJyWlhraUxDSXlNREkxTURjeU0xOWFhRzl1WjB0bFdXRnVSM1ZoYm1kYWFHOTFTbWxoYmt0aGJtZFpZVzVLYVhWWmRXRnVYekV2SWwwc2V5SmlkV05yWlhRaU9pSnJaV1oxTFhONmRHOXpMV3hwYm1zaWZWMTkmWC1Ub3MtU2lnbmF0dXJlPTQzZWM5MzFmMWRkNGI5OWZmZWFhNjYzMDFiZDJhMjUxMjVhZDRmZjNkNjMxOGJjNzAxMzA0MDExNWRkNTA5MDY=
url=https://document-share.tos-cn-beijing.volces.com/v2/index.html?token=aHR0cHM6Ly9rZWZ1LXN6dG9zLWxpbmsudG9zLWNuLWd1YW5nemhvdS52b2xjZXMuY29tLz9YLVRvcy1BbGdvcml0aG09VE9TNC1ITUFDLVNIQTI1NiZYLVRvcy1DcmVkZW50aWFsPUFLTFRNMlF3T0RZelptVmxOelEyTkRReU9UazBaRFkwWkRSaVlUZGlPRFJsWldFJTJGMjAyNTA4MTQlMkZjbi1ndWFuZ3pob3UlMkZ0b3MlMkZyZXF1ZXN0JlgtVG9zLURhdGU9MjAyNTA4MTRUMTE0MDM2WiZYLVRvcy1FeHBpcmVzPTEyOTYwMDAmWC1Ub3MtUG9saWN5PWV5SmpiMjVrYVhScGIyNXpJanBiV3lKemRHRnlkSE10ZDJsMGFDSXNJaVJyWlhraUxDSXlNREkxTURneE5GOWFhRzl1WjB0bFdXRnVSM1ZoYm1kYWFHOTFTbWxoYmt0aGJtZFpZVzVLYVhWWmRXRnVYekV2SWwwc2V5SmlkV05yWlhRaU9pSnJaV1oxTFhONmRHOXpMV3hwYm1zaWZWMTkmWC1Ub3MtU2lnbmF0dXJlPTg5ZjQ5Y2FlZjcwNWExZTBkY2U1Y2QyZTY1MGM5ZmY0NmI4NTYxYjU2MDE4MzNjZTc5ODgyNTE1YTMyODRmNTA=
url='https://document-share.tos-cn-beijing.volces.com/v2/index.html?token=aHR0cHM6Ly9rZWZ1LXN6dG9zLWxpbmsudG9zLWNuLWd1YW5nemhvdS52b2xjZXMuY29tLz9YLVRvcy1BbGdvcml0aG09VE9TNC1ITUFDLVNIQTI1NiZYLVRvcy1DcmVkZW50aWFsPUFLTFRNMlF3T0RZelptVmxOelEyTkRReU9UazBaRFkwWkRSaVlUZGlPRFJsWldFJTJGMjAyNTA5MTIlMkZjbi1ndWFuZ3pob3UlMkZ0b3MlMkZyZXF1ZXN0JlgtVG9zLURhdGU9MjAyNTA5MTJUMTEyODEyWiZYLVRvcy1FeHBpcmVzPTEyOTYwMDAmWC1Ub3MtUG9saWN5PWV5SmpiMjVrYVhScGIyNXpJanBiV3lKemRHRnlkSE10ZDJsMGFDSXNJaVJyWlhraUxDSXlNREkxTURreE1sOWFhRzl1WjB0bFdXRnVSM1ZoYm1kYWFHOTFTbWxoYmt0aGJtZFpZVzVLYVhWWmRXRnVYekV2SWwwc2V5SmlkV05yWlhRaU9pSnJaV1oxTFhONmRHOXpMV3hwYm1zaWZWMTkmWC1Ub3MtU2lnbmF0dXJlPWUzMDJlNjhiMjExNDIxNGI3MTJlODU5MTE4ZmNiYzYyM2Q3YWNlZmJiMTU0NzFlNWQ2ODJkMWVlOTE1OTQ0NjA='

tosutil share-cp ${url} ./ -r -u


for per_bam in $(find */*/03.star -name '*_Aligned.sortedByCoord.out.bam'); do
  per_sample=$(cut -f2 -d/ <<<${per_bam})
  bam_dir=$(dirname ${per_bam})

  out_old_bam=${bam_dir}/${per_sample}_Aligned.sortedByCoord.out.old.bam
  out_nascent_bam=${bam_dir}/${per_sample}_Aligned.sortedByCoord.out.nascent.bam
  if [[ -e ${out_nascent_bam} && -e ${out_old_bam} ]]; then continue; fi
  echo "[I]: Working on ${per_sample} ..."
  apptainer exec -B /mnt/backup /home/zzhang/tools/containers/images/celatlas_spatial.sif \
    python /opt/tools/scripts/extract_new_reads.py \
    -@ 10 \
    -n ${out_nascent_bam} \
    -o ${out_old_bam} \
    ${per_bam}
done




#
## Jobs
#

# Plot grand slam results
Rscript ~/Documents/projects/wp_vasaseq/scripts/r/plot_grand_slam_results.r --help
Rscript ~/Documents/projects/wp_vasaseq/scripts/r/plot_grand_slam_results.r \
  -c 20250226LIB1 \
  -o 20250226_encoded_with_polya_embryo.pdf \
  ~/Documents/projects/wp_vasaseq/outputs/analysis/preprocessing/slide_seq/nascent_rna/20250226_encoded_with_polya_embryo/per_sample

Rscript ~/Documents/projects/wp_vasaseq/scripts/r/plot_grand_slam_results.r \
  -c nc_no4suctrl_r1,nc_no4suctrl_r2 \
  -o batch_20250321_byPLY.kideny.pdf \
  -H 20 \
  ~/Documents/projects/wp_vasaseq/outputs/analysis/preprocessing/geo_slam_seq/batch_20250321_byPLY/nascent_rna

Rscript ~/Documents/projects/wp_vasaseq/scripts/r/plot_grand_slam_results.r \
  -c nc_no4suctrl_r3 \
  -o batch_20250321_byLDJ.kideny.pdf \
  -H 12 \
  ~/Documents/projects/wp_vasaseq/outputs/analysis/preprocessing/geo_slam_seq/batch_20250321_byLDJ/nascent_rna

Rscript ~/Documents/projects/wp_vasaseq/scripts/r/plot_grand_slam_results.r \
  -c 20250314Lib1 \
  -o batch_20250314_decoded_embryo.pdf \
  -W 12 \
  ~/Documents/projects/wp_vasaseq/outputs/analysis/preprocessing/slide_seq/nascent_rna/20250314_decoded_embryo

Rscript ~/Documents/projects/wp_vasaseq/scripts/r/plot_grand_slam_results.r \
  -c 20250319Lib1 \
  -o batch_20250319_decoded_embryo.pdf \
  -W 12 \
  ~/Documents/projects/wp_vasaseq/outputs/analysis/preprocessing/slide_seq/nascent_rna/20250319_decoded_embryo

Rscript ~/Documents/projects/wp_vasaseq/scripts/r/plot_grand_slam_results.r \
  -c 20250319Lib1 \
  -o GEO_SLAM_20250420_vasaseq.pdf \
  -W 12 \
  ~/Documents/projects/wp_vasaseq/outputs/analysis/preprocessing/geo_slam/GEO_SLAM_20250420_vasaseq/nascent_rna

Rscript ~/Documents/projects/wp_vasaseq/scripts/r/plot_grand_slam_results.r \
  -c 20250319Lib1 \
  -o GEO_SLAM_20250419_smartseq2.pdf \
  -W 12 \
  ~/Documents/projects/wp_vasaseq/outputs/analysis/preprocessing/geo_slam/GEO_SLAM_20250419_smartseq2/nascent_rna

Rscript ~/Documents/projects/wp_vasaseq/scripts/r/plot_grand_slam_results.r \
  -c 20250319Lib1 \
  -o GEO_SLAM_20250525_Lib1.pdf \
  -W 12 \
  ~/Documents/projects/wp_vasaseq/outputs/analysis/preprocessing/geo_slam/GEO_SLAM_20250525_Lib1/nascent_rna

Rscript ~/Documents/projects/wp_vasaseq/scripts/r/plot_grand_slam_results.r \
  -c 20250319Lib1 \
  -o GEO_SLAM_20250525_Lib2.pdf \
  -W 12 \
  ~/Documents/projects/wp_vasaseq/outputs/analysis/preprocessing/geo_slam/GEO_SLAM_20250525_Lib2/nascent_rna

Rscript ~/Documents/projects/wp_vasaseq/scripts/r/plot_grand_slam_results.r \
  -c 20250319Lib1 \
  -o GEO_SLAM_20250525_Lib3.pdf \
  -W 12 \
  ~/Documents/projects/wp_vasaseq/outputs/analysis/preprocessing/geo_slam/GEO_SLAM_20250525_Lib3/nascent_rna


#
## Estimate development stage 
#
# Overview of available tools
python ${project_dir}/scripts/py3/embryo_development_stage.py --help

# 1. Preprocess data
python ${project_dir}/scripts/py3/embryo_development_stage.py preproc --help
python ${project_dir}/scripts/py3/embryo_development_stage.py preproc \
  -F \
  -o PijuanSala_etal_Nature_2019/Preprocess PijuanSala_etal_Nature_2019/PijuanSala_etal_Nature_2019.raw.h5ad

# 2. Create meta-cells
python ${project_dir}/scripts/py3/embryo_development_stage.py metacell --help
python -m pdb ${project_dir}/scripts/py3/embryo_development_stage.py metacell \
  -F -g louvain preprocessed.h5ad

# 3. Train a model.
python ${project_dir}/scripts/py3/embryo_development_stage.py train --help

# 4. Validate the model.

# 5. Predict new samples using given model.

# 6. Evaluation plots.



# Check FPKM at different thresholds.
echo -e "Region,MinFPKM,NumberOfGenes" > ${project_dir}/outputs/analysis/geo_slam/overview/FPKM_per_threshold.csv
for region in {6,8,10,15,13}EA {10,13,15}MP {6,8,10}EP {2,4,6,8,10,13,15}P {2,4,6,8,10,13,15}A {10,13,15}MA; do
  awk -F$'\t' -v REGION=${region} -f- <<'EOF' \
    $(find ${project_dir}/outputs/analysis/preprocessing/geo_slam/ -name ${region}.genes.read_counts.txt)
$1 ~ /^#/ || $0 ~ /Geneid/ { next }
{ count_tbl[$1] = $8; length_tbl[$1] = $4 - $3; ttl_counts += $8 }
END {
  OFS = ","
  for (i in count_tbl) { fpkm_tbl[i] = (count_tbl[i] * 1000000000) / (ttl_counts * length_tbl[i]) }

  t_tbl[1] = 0; t_tbl[2] = 0.1; t_tbl[3] = 0.5; t_tbl[4] = 1; t_tbl[5] = 5; t_tbl[6] = 10; t_tbl[7] = 50; t_tbl[8] = 100
  for (ii in t_tbl) {
    gene_count = 0
    for (i in fpkm_tbl) { gene_count += fpkm_tbl[i] > t_tbl[ii] ? 1 : 0 }
    print REGION,t_tbl[ii],gene_count
  }
}
EOF
done >> ${project_dir}/outputs/analysis/geo_slam/overview/FPKM_per_threshold.csv
head -10 ${project_dir}/outputs/analysis/geo_slam/overview/FPKM_per_threshold.csv

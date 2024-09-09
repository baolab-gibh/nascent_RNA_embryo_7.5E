#!/bin/bash
# File: workflow.sh
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Jun 03, 2024
# Updated: Jun 17, 2024

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
# SLAM-seq pipeline
work_dir=${project_dir}/outputs/analysis
smk_file=${project_dir}/scripts/snakemake/vasalam_ppl.smk
cfg_file=${project_dir}/scripts/snakemake/vasalam_ppl.yml

# Job graph, not able to show the DAG of modules
snakemake --forceall -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --config myconfigfile=${cfg_file} --rulegraph | dot -Tpdf >| ${work_dir}/preprocessing/logs/main_ppl.job_graph.pdf

# Touch existing files
snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --config myconfigfile=${cfg_file} --touch 2>/dev/null

# Dry-run
snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --config myconfigfile=${cfg_file} --dry-run >| ${work_dir}/preprocessing/logs/main_ppl.dry_run.txt
awk '/^Job stats:$/ {count++} count >= 2 {print}' ${work_dir}/preprocessing/logs/main_ppl.dry_run.txt
cat ${work_dir}/preprocessing/logs/main_ppl.dry_run.txt
less -X -S -N ${work_dir}/preprocessing/logs/main_ppl.dry_run.txt

# Unlock current .snakemake folder
snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --config myconfigfile=${cfg_file} --unlock

# Run
if [[ ${USER} == 'zhzhang_gibh' ]]; then
  snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --config myconfigfile=${cfg_file}
elif [[ ${USER} == 'zzhang' ]]; then
  snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --config myconfigfile=${cfg_file}
fi


# 2019 Peng et al Nature
work_dir=${project_dir}/outputs/analysis
smk_file=${project_dir}/scripts/snakemake/2019_peng_etal_nature.smk

# Job graph, not able to show the DAG of modules
snakemake --forceall -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --rulegraph | dot -Tpdf >| ${work_dir}/preprocessing/logs/2019_peng_etal_nature.job_graph.pdf

# Touch existing files
snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --touch 2>/dev/null

# Dry-run
snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --dry-run >| ${work_dir}/preprocessing/logs/2019_peng_etal_nature.dry_run.txt
awk '/^Job stats:$/ {count++} count >= 2 {print}' ${work_dir}/preprocessing/logs/2019_peng_etal_nature.dry_run.txt
cat ${work_dir}/preprocessing/logs/2019_peng_etal_nature.dry_run.txt
less -X -S -N ${work_dir}/preprocessing/logs/2019_peng_etal_nature.dry_run.txt

# Unlock current .snakemake folder
snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir} --unlock

# Run
if [[ ${USER} == 'zhzhang_gibh' ]]; then
  snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file} --profile ${prof_dir}
elif [[ ${USER} == 'zzhang' ]]; then
  snakemake -d ${work_dir} -j ${n_jobs} -s ${smk_file}
fi



# reference_files
star_genome_index=/public/home/zbmai_gibh/STuDY/vasa-seq-test/mm10

# Split the input fastq file based
awk -f- <<'EOF' <(zcat 240409_Lib_embryo_raw_1.fq.gz | head -6600000) <(zcat 240409_Lib_embryo_raw_2.fq.gz | head -6600000)
BEGION {MEAN_READ_LENGTH = 0; MAX_READ_LENGTH = 0; MIN_READ_LENGTH = 0}

NR % 4 == 1 { read_name = $0 }
NR % 4 == 2 { read_bases = $0 }
NR % 4 == 3 { next }
NR % 4 == 0 {
  read_qual = $0
  if (NR == FNR) {
    bc_bases = substr(read_bases, 1, 6); R1_BC_COUNT[bc_bases]++
    umi_bases = substr(read_bases, 1, 14); R1_UMI_COUNT[umi_bases]++
  } else {
    bc_bases = substr(read_bases, length(read_bases) - 5, 6); R2_BC_COUNT[bc_bases]++
    umi_bases = substr(read_bases, length(read_bases) - 13, 14); R2_UMI_COUNT[umi_bases]++
  }
}

END {
  print "PE,Barcode,ReadCount"
  for (ii in R1_BC_COUNT) {
    if (R1_BC_COUNT[ii] > 10000) {
      print "R1", ii, R1_BC_COUNT[ii]
      r1_n_bc++
    }
  }

  for (ii in R2_BC_COUNT) {
    if (R2_BC_COUNT[ii] > 10000) {
      print "R2", ii, R2_BC_COUNT[ii]
      r2_n_bc++
    }
  }

  print "# Nr. of enriched barcode (R1): " r1_n_bc
  print "# Nr. of enriched barcode (R2): " r2_n_bc
}
EOF


# Count Ti and Tv from the bam files
python -m pdb ${project_dir}/scripts/py3/t2c_counter.py \
  -@ 2 \
  ${project_dir}/outputs/analysis/240703_Lib_32region/read_alignments/11A_TCAGAC/11A_TCAGAC.merged.rmdup.bam \
  ${project_dir}/outputs/analysis/240703_Lib_32region/read_alignments/11EA_TAGTCG/11EA_TAGTCG.merged.rmdup.bam

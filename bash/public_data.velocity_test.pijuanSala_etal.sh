#!/bin/bash

cat <<EOF
A pipeline to analyze the velocity of RNA expression from embryonic single-cell RNA-seq data (E7.0-E7.5)

Data repository:
  - https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6967

Publication:
  - https://dx.doi.org/10.1038%2Fs41586-019-0933-9
EOF


# Meta-variables
project_dir=~/Documents/projects/wp_vasaseq


# Download meta data (aka informations about samples), e.g., sample name, stages, and FASTQ file URLs etc.
dst_file=${project_dir}/inputs/PijuanSala_Nature_2019/E-MTAB-6967.sdrf.txt
src_url=https://www.ebi.ac.uk/biostudies/fire/E-MTAB-/967/E-MTAB-6967/Files/E-MTAB-6967.sdrf.txt
if [[ -e ${dst_file} ]]; then
  echo "[I]: ${dst_file} already exists. Skip downloading."
else
  echo "[I]: Downloading ${src_url} to ${dst_file}"
  wget -c -O ${dst_file} ${src_url}
  if [[ $? -ne 0 ]]; then echo "[E]: Failed to download ${src_url} to ${dst_file}"; exit; fi
fi


# Download processed data.
dst_file=${project_dir}/inputs/PijuanSala_Nature_2019/atlas_data.tar.gz
src_url=https://www.ebi.ac.uk/biostudies/fire/E-MTAB-/967/E-MTAB-6967/Files/atlas_data.tar.gz
if [[ -e ${dst_file} ]]; then
  echo "[I]: ${dst_file} already exists. Skip downloading."
else
  echo "[I]: Downloading ${src_url} to ${dst_file}"
  wget -c -O ${dst_file} ${src_url}
  if [[ $? -ne 0 ]]; then echo "[E]: Failed to download ${src_url} to ${dst_file}"; exit; fi
fi


# Download selected FASTQ files
meta_info_file=${project_dir}/inputs/PijuanSala_Nature_2019/E-MTAB-6967.sdrf.txt
head -1 ${meta_info_file} | tr '\t' '\n' | cat -n
cut -f 1,2,3,6 ${meta_info_file} | csvtk pretty -d$'\t' | head
for sample_id in {10,14,15,30,31,32}; do
  [[ ! -e Sample_${sample_id} ]] && mkdir -p Sample_${sample_id} || echo "[I]: Found dir Sample_${sample_id}"
  cd Sample_${sample_id}
  parallel -j 8 -n 1 curl -L -O -C- {} ::: \
      $(awk -v SAMPLE_ID=${sample_id} -F$'\t' '$1 == "Sample "SAMPLE_ID { OFS="\n"; print $54,$56,$58,$60 }' ${meta_info_file})
  cd -
done

# To download uBAMs
# bam_files=${project_dir}/inputs/PijuanSala_Nature_2019/filereport_read_run_PRJEB28586_tsv.txt
# head -1 ${bam_files} | tr '\t' '\n' | cat -n
# cut -f 10,11 ${bam_files} | csvtk pretty -d$'\t' | head
# $(awk -F$'\t' -v SAMPLE_ID=${sample_id} '$11=="Sample "SAMPLE_ID {split($10,URL,";"); print "https://"URL[1]"\nhttps://"URL[2] }' ${bam_files})


# Cellranger count
# Cellranger: 9.0.0
# Reference: refdata-gex-mm10-2020-A
# Prob sets: Chromium_Mouse_Transcriptome_Probe_Set_v1.1.0_GRCm39-2024-A.csv
# Nucleic acid extraction protocol: 10X Genomics, v2 chemistry

# library_info=${project_dir}/inputs/PijuanSala_Nature_2019/Embryo7.library_info.csv
# (echo fastqs,sample,library_type,; for x in $(ls *.fastq.gz | cut -f1-3 -d_ | sort -u); do echo $(pwd),$x,"Gene Expression,"; done) > ${library_info}
# cellranger count --id ${sample_id} --transcriptome ${reference_dir} --libraries ${library_info} --output-dir ${out_dir} --create-bam true --localcores 4 --localmem 32

library_info_v1=${project_dir}/inputs/PijuanSala_Nature_2019/Embryo7.library_info.csv
library_info_v2=${project_dir}/inputs/PijuanSala_Nature_2019/Embryo7.library_info.v2.txt
awk -F ',' -f- <<'EOF' ${library_info_v1} > ${library_info_v2}
NR == 1 {print "fastq_path\tsample"; next}
{
  split($2, INFO, "_")
  barcode = INFO[3]
  if (whole_dataset[barcode]) { whole_dataset[barcode] = whole_dataset[barcode]","$2 } else { whole_dataset[barcode] = $2 }
}
END {
  base_url = "/home/zzhang/Documents/projects/wp_vasaseq/outputs/analysis/public/PijuanSala_Nature_2019/fastq/outs/fastq_path"
  for (i in whole_dataset) { print base_url"\t"whole_dataset[i] }
}
EOF

# sample_id=22111_2_TCTTAGGC
library_info=${project_dir}/inputs/PijuanSala_Nature_2019/Embryo7.library_info.v2.txt
reference_dir=${HOME}/Documents/projects/resources/10X/Mouse/Reference/refdata-gex-mm10-2020-A/

while read fastq_dir samples; do
  if [[ ${samples} == sample ]]; then continue; fi

  sample_id=$(echo $samples | grep -Eo '[TCGA]+' | uniq)
  out_dir=${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/count/${sample_id}
  cellranger count \
    --id ${sample_id} \
    --transcriptome ${reference_dir} \
    --sample ${samples} \
    --fastqs ${fastq_dir} \
    --output-dir ${out_dir} \
    --create-bam true \
    --localcores 4 \
    --localmem 32

  if [[ $? -ne 0 ]]; then echo "[E]: Failed to count ${samples}"; break; fi
done < ${library_info}


if [[ optional == "yes" ]]; then
  # Aggregate all GEM files
  library_info=${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/count/library_count.csv
  out_dir=${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/aggr
  cellranger aggr --id embryo7 --output-dir ${out_dir} --disable-ui --normalize mapped --csv ${library_info}

  # Merge all BAM files. Here the simply merging does not work due to the duplicated barcodes.
  bam_dir=${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/count/
  out_bamfile=${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/aggr/possorted_genome_bam.bam 
  cd ${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/count
  samtools cat -@ 4 */outs/possorted_genome_bam.bam \
    | samtools sort -@ 4 -O bam -o ../velocity/possorted_genome_bam.bam 
  samtools index -@ 8 ${out_bamfile}
fi


# Obtain velocity results, the time derivative of the gene expression state
in_dir=${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/aggr/
out_dir=${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/velocity/velocyto
masks=${HOME}/Documents/projects/resources/UCSCGenombrowser/mm10_rmsk.gtf
gtf_file=${HOME}/Documents/projects/resources/10X/Mouse/Reference/refdata-gex-mm10-2020-A/genes/genes.gtf
bam_dir=${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/count
# python -m pdb $(which velocyto) run10x -v -@ 4 -m ${masks} ${in_dir} ${gtf_file}
# velocyto run10x -@ 4 -m ${masks} ${in_dir} ${gtf_file}
cd ${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/count
bam_files=$(find ${bam_dir} -type f -name '*.bam')
velocyto run -o ${out_dir} -@ 4 -m ${masks} ${bam_files} ${gtf_file}


# Data analysis
# Author's QC pipepline https://github.com/MarioniLab/EmbryoTimecourse2018/tree/master/analysis_scripts/chimera-wt
in_dir=${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/10X
out_dir=${project_dir}/outputs/analysis/public/PijuanSala_Nature_2019/velocity
Rscript velocity_test.pijuanSala_etal.R ${in_dir} ${out_dir}


# Available samples
# 22109_1_TCGCAATT 22111_1_CAAGTCCA 22111_1_TCGCAATT 22111_2_GTGAGAAG 22089_1_GTCTTTGA 22109_2_TCGCAATT 22108_2_CATGGCAG 
# 22108_2_GTCTTTGA 22109_2_CAAGTCCA 22111_1_GTGAGAAG 22089_2_GTCTTTGA 22111_2_TCGCAATT 22089_2_TCTTAGGC 22089_2_TCGCAATT 
# 22109_1_AGAACGCC 22108_2_GTGAGAAG 22089_2_GTGAGAAG 22108_1_TCGCAATT 22089_1_TCTTAGGC 22109_2_AGAACGCC 22108_1_AGCCCTTT 
# 22109_1_GTGAGAAG 22109_2_AGCCCTTT 22089_1_AGAACGCC 22089_1_GTGAGAAG 22109_2_GTCTTTGA 22109_2_GTGAGAAG 22109_2_CATGGCAG 
# 22111_1_AGAACGCC 22108_2_TCGCAATT 22111_2_TCTTAGGC 22108_1_GTCTTTGA 22108_1_TCTTAGGC 22089_1_CAAGTCCA 22111_1_AGCCCTTT 
# 22108_1_CATGGCAG 22111_2_AGCCCTTT 22108_1_GTGAGAAG 22108_1_AGAACGCC 22109_1_CAAGTCCA 22108_2_AGAACGCC 22089_2_CAAGTCCA 
# 22108_2_AGCCCTTT 22089_1_TCGCAATT 22111_2_CATGGCAG 22089_1_AGCCCTTT 22089_2_CATGGCAG 22108_1_CAAGTCCA 22089_2_AGAACGCC 
# 22111_2_GTCTTTGA 22109_1_TCTTAGGC 22111_1_TCTTAGGC 22108_2_CAAGTCCA 22111_1_CATGGCAG 22111_2_CAAGTCCA 22111_2_AGAACGCC 
# 22109_2_TCTTAGGC 22108_2_TCTTAGGC 22089_2_AGCCCTTT 22109_1_GTCTTTGA 22089_1_CATGGCAG 22109_1_AGCCCTTT 22111_1_GTCTTTGA 
# 22109_1_CATGGCAG

#!/bin/bash
# File: define_erna_regions.sh
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Jun 23, 2024
# Updated:

project_dir=${HOME}/Documents/projects/wp_vasaseq


#
## Define eRNA regions
#
# Reference: https://doi.org/10.1038/s41467-019-12543-5
# Title: Transcriptional landscape and clinical utility of enhancer RNAs for eRNA-targeted therapy in cancer

# Coding or lncRNA regions
gencode_mouse=${HOME}/Documents/projects/resources/Gencode/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M29/GRCm38_mapping/gencode.vM29lift38.annotation.gtf.gz
gencode_sense_regions=${project_dir}/outputs/references/genomic_features/gencode.vM29lift38.gene_features.gtf
awk -F $'\t' -f- <<'EOF' <(zcat ${gencode_mouse}) | sed 's/^chr//g' | sort -k1,1V -k2,2g -k3,3g >| ${gencode_sense_regions}
{
  OFS = "\t"
  split($9, ATTR, "; ")
  for (i in ATTR) { split(ATTR[i], kvp, " "); attr_dict[kvp[1]] = kvp[2] }
  gene_type = attr_dict["gene_type"]; gene_name = attr_dict["gene_name"]; gene_id = attr_dict["gene_id"]
  if ($3 ~ /CDS|start_codon|stop_codon/ && attr_dict["gene_type"] ~ /IG_|TR_|protein_coding/) {
    print $1, $4-1, $5-1, $3";"gene_id";"gene_name";Protein_coding"
  } else if ($3 ~ /transcript/ && attr_dict["gene_type"] ~ /lncRNA|lincRNA|TEC/) {
    print $1, $4-1, $5-1, $3";"gene_id";"gene_name";Non_coding"
  }

  delete attr_dict
}
EOF

# Database: FANTOM (http://fantom.gsc.riken.jp)
# Version: FANTOM5
# URL: https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.bed.gz
# URL: https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.expression.matrix.gz
# URL: https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.expression.usage.matrix.gz
fantom_enhancer=${HOME}/Documents/projects/resources/FANTOM/fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.bed.gz
zcat ${fantom_enhancer} \
  | sed 's/^chr//g' \
  | bedtools intersect -v -a stdin -b ${gencode_sense_regions} \
  | sort -k1,1V -k2,2g -k3,3g \
  | cut -f1-3 \
  | awk '{OFS = "\t"; print $1, $2 - 100 < 0 ? 0 : $2 - 100, $3 + 100}' \
  >| ${project_dir}/outputs/references/genomic_features/fantom5_enhancer.non_overlapped_with_genes.bed


# Database: ENSEMBL (https://www.ensembl.org/)
# Version: release-112
# URL: https://ftp.ensembl.org/pub/release-112/regulation/mus_musculus/mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20221007.gff.gz
ensemnl_enhancer=${HOME}/Documents/projects/resources/Ensembl/ftp.ensembl.org/pub/release-112/regulation/mus_musculus/mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20221007.gff.gz
zcat ${ensemnl_enhancer} \
  | awk '$3 == "enhancer"' \
  | bedtools intersect -v -a stdin -b ${gencode_sense_regions} \
  | sort -k1,1V -k4,4g -k5,5g \
  | cut -f1,4,5 \
  | awk '{OFS = "\t"; print $1, $2 - 100 < 0 ? 0 : $2 - 100, $3 + 100}' \
  >| ${project_dir}/outputs/references/genomic_features/ensembl_enhancer.non_overlapped_with_genes.bed


# Obtain shared enhancers regions between FANTOM5 and ENSEMBL
candidate_erna_regions=${project_dir}/outputs/references/genomic_features/overlapped.fantom5_and_ensembl.gtf
bedtools intersect -wa \
  -a ${project_dir}/outputs/references/genomic_features/fantom5_enhancer.non_overlapped_with_genes.bed \
  -b ${project_dir}/outputs/references/genomic_features/ensembl_enhancer.non_overlapped_with_genes.bed \
  | bedtools merge -i stdin -d 100 \
  | awk -F $'\t' '{OFS = "\t"; print $1, "FANTOM;ENSEMBL", "enhancer", $2+1, $3+1, ".", "*", ".", "enhancer_region "$1":"$2+1"-"$3+1 }' \
  | sort -k1,1V -k4,4g -k5,5g -u \
  >| ${candidate_erna_regions}


# Prepare Salmon index
mm10_fasta=${project_dir}/outputs/references/genome/mus_musculus.90.fasta
candidate_erna_fasta=${project_dir}/outputs/references/genome/overlapped.fantom5_and_ensembl.mus_musculus_90.fasta
bedtools getfasta -fi ${mm10_fasta} -bed ${candidate_erna_regions} > ${candidate_erna_fasta}


#
## Identify potential functions
#
# TOBIAS, for ATAC-seq but potential usable for eRNAs
# HOMMER, findpeaks.pl to define super enhancer

# Transcription unit (TU)
gencode_mouse=${HOME}/Documents/projects/resources/Gencode/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M29/GRCm38_mapping/gencode.vM29lift38.annotation.gtf.gz
gencode_sense_regions=${project_dir}/outputs/references/genomic_features/gencode.vM29lift38.transcription_unit.gtf
awk -F $'\t' -f- <<'EOF' <(zcat ${gencode_mouse}) | sed 's/^chr//g' | sort -k1,1V -k2,2g -k3,3g >| ${gencode_sense_regions}
BEGIN {
  selected_gene_type["\"lincRNA\""] = 1
  selected_gene_type["\"lncRNA\""] = 1
  selected_gene_type["\"miRNA\""] = 1
  selected_gene_type["\"misc_RNA\""] = 1
  selected_gene_type["\"processed_transcript\""] = 1
  selected_gene_type["\"scaRNA\""] = 1
  selected_gene_type["\"scRNA\""] = 1
  selected_gene_type["\"snoRNA\""] = 1
  selected_gene_type["\"snRNA\""] = 1
  selected_gene_type["\"sRNA\""] = 1
  selected_gene_type["\"TEC\""] = 1
} {
  OFS = "\t"
  split($9, ATTR, "; ")
  for (i in ATTR) { split(ATTR[i], kvp, " "); attr_dict[kvp[1]] = kvp[2] }
  gene_type = attr_dict["gene_type"]; gene_name = attr_dict["gene_name"]; gene_id = attr_dict["gene_id"]
  if ($3 ~ /transcript/ && selected_gene_type[gene_type] == 1) {
    print $0
  }

  delete attr_dict
}
EOF

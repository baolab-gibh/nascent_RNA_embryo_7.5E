#!/bin/bash
# File: celatlas_spatial.pipeline.v2.sh
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Apr 09, 2025
# Updated: Apr 09, 2025

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

trap 'error Errors at $LINENO' ERR


project_dir=~/Documents/projects/wp_vasaseq

if [[ ${USER} == 'zhzhang_gibh' ]]; then
  source ~/tools/miniconda3/bin/activate
  mamba activate sm732 # Snakemake v7.32.4
elif [[ ${USER} =~ 'zzhang' && -d ${project_dir}/scripts/.env ]]; then
  source ${project_dir}/scripts/.env/bin/activate
fi

cd ${project_dir}/temps || return 1


# ----- TO be translate into snakemake ------------------
# Slide decoded
# Utility functions
function define_params() {
  local sample_id which_param barcode_map col_idx

  sample_id=$1; which_param=$2
  declare -A barcode_map=(
    # 2025-02-26
    [20250226LIB1]=20250226_decoded_embryo,OST110029,barcodes/OST110029.barcodeToPos.h5,images/32.jpg,raw/20250226LIB1_L04_R1.fq.gz,raw/20250226LIB1_L04_R2.fq.gz
    [20250226LIB2]=20250226_decoded_embryo,OST110030,barcodes/OST110030.barcodeToPos.h5,images/30.jpg,raw/20250226LIB2_L03_R1.fq.gz,raw/20250226LIB2_L03_R2.fq.gz
    [20250226LIB3]=20250226_decoded_embryo,OST110087,barcodes/OST110087.barcodeToPos.h5,images/21.jpg,raw/20250226LIB3_L04_R1.fq.gz,raw/20250226LIB3_L04_R2.fq.gz
    [20250226LIB4]=20250226_decoded_embryo,OST110088,barcodes/OST110088.barcodeToPos.h5,images/19.jpg,raw/20250226LIB4_L03_R1.fq.gz,raw/20250226LIB4_L03_R2.fq.gz

    # 2025-03-14
    [20250314Lib1]=20250314_decoded_embryo,OST110090,OST110090/OST110090.barcodeToPos.h5,images/01.jpg,OST110090/BBV2.4/OST110090_1.fq.gz,OST110090/BBV2.4/OST110090_2.fq.gz
    [20250314Lib2]=20250314_decoded_embryo,OST110091,OST110091/OST110091.barcodeToPos.h5,images/01.jpg,OST110091/BBV2.4/OST110091_1.fq.gz,OST110091/BBV2.4/OST110091_2.fq.gz
    [20250314Lib3]=20250314_decoded_embryo,OST110092,OST110092/OST110092.barcodeToPos.h5,images/01.jpg,OST110092/BBV2.4/OST110092_1.fq.gz,OST110092/BBV2.4/OST110092_2.fq.gz
    [20250314Lib4]=20250314_decoded_embryo,OST110093,OST110093/OST110093.barcodeToPos.h5,images/01.jpg,OST110093/BBV2.4/OST110093_1.fq.gz,OST110093/BBV2.4/OST110093_2.fq.gz
    [20250314Lib5]=20250314_decoded_embryo,OST110094,OST110094/OST110094.barcodeToPos.h5,images/01.jpg,OST110094/BBV2.4/OST110094_1.fq.gz,OST110094/BBV2.4/OST110094_2.fq.gz

    # 2025-03-19
    [20250319Lib1]=20250319_decoded_embryo,OST110083,OST110083/OST110083.barcodeToPos.h5,images/01.jpg,OST110083/BBV2.4/OST110083_1.fq.gz,OST110083/BBV2.4/OST110083_2.fq.gz
    [20250319Lib2]=20250319_decoded_embryo,OST110084,OST110084/OST110084.barcodeToPos.h5,images/01.jpg,OST110084/BBV2.4/OST110084_1.fq.gz,OST110084/BBV2.4/OST110084_2.fq.gz
    [20250319Lib3]=20250319_decoded_embryo,OST110085,OST110085/OST110085.barcodeToPos.h5,images/01.jpg,OST110085/BBV2.4/OST110085_1.fq.gz,OST110085/BBV2.4/OST110085_2.fq.gz
    [20250319Lib4]=20250319_decoded_embryo,OST110086,OST110086/OST110086.barcodeToPos.h5,images/01.jpg,OST110086/BBV2.4/OST110086_1.fq.gz,OST110086/BBV2.4/OST110086_2.fq.gz
    [20250319Lib5]=20250319_decoded_embryo,OST110089,OST110089/OST110089.barcodeToPos.h5,images/01.jpg,OST110089/BBV2.4/OST110089_1.fq.gz,OST110089/BBV2.4/OST110089_2.fq.gz

    # 2025-04-03, kideny
    [20250319Lib1]=20250319_decoded_embryo,OST110083,OST110083/OST110083.barcodeToPos.h5,images/01.jpg,OST110083/BBV2.4/OST110083_1.fq.gz,OST110083/BBV2.4/OST110083_2.fq.gz
    [20250319Lib2]=20250319_decoded_embryo,OST110084,OST110084/OST110084.barcodeToPos.h5,images/01.jpg,OST110084/BBV2.4/OST110084_1.fq.gz,OST110084/BBV2.4/OST110084_2.fq.gz
  )

  case ${which_param} in
    batch) col_idx=1 ;;
    chipid) col_idx=2 ;;
    barcodes) col_idx=3 ;;
    images) col_idx=4 ;;
    R1) col_idx=5 ;;
    R2) col_idx=6 ;;
    *) echo "[E] Unknown parameter: ${which_param}" >&2; return 1 ;;
  esac

  echo $(cut -d',' -f ${col_idx} <<<${barcode_map[${sample_id}]})

  return 0
}

ulimit -n 10240

apptainer_image=~/tools/containers/images/celatlas_spatial.sif
genome_dir=~/Documents/projects/resources/references/GRCm38.celatlas_spatial


n_cpus=10
bin_size=20
bin_size=50
chemistry=customized
chemistry_pattern=C4L15C4L15C4U10T18
insert_size=150
feature_type=gene
expected_cell_num=2000
pixel_size=0.24
segment_method=gene_expr
model_dir=/opt/tools/database/celatlas_spatial/1.5.0/swin_tiny.pth


# sample_id=20250314Lib3
sample_id=20250319Lib4
for sample_id in 2025031{4,9}Lib{1..5}; do
  batch=$(define_params ${sample_id} batch)
  echo "~~~~~~~~~~ Starting ... Work on ${sample_id} ~~~~~~~~~~"
  in_dir=${project_dir}/inputs/slide_seq/${batch}
  fastq_r1=${in_dir}/$(define_params ${sample_id} R1)
  fastq_r2=${in_dir}/$(define_params ${sample_id} R2)
  barcodes=${in_dir}/$(define_params ${sample_id} barcodes)
  image_file=${in_dir}/$(define_params ${sample_id} images)
  chip_id=$(define_params ${sample_id} chipid)
  chip_data_dir=$(dirname ${barcodes})
  sample_out_dir=${project_dir}/outputs/analysis/preprocessing/slide_seq_decoded/${batch}/${sample_id}

  # Normalize sample id
  mkdir -p ${sample_out_dir}/00.sample
  apptainer exec ${apptainer_image} ln -s ${fastq_r1} ${sample_out_dir}/00.sample/${chip_id}_1.fq.gz
  apptainer exec ${apptainer_image} ln -s ${fastq_r2} ${sample_out_dir}/00.sample/${chip_id}_2.fq.gz

  # Create barcode
  mkdir -p ${sample_out_dir}/01.barcode
  apptainer exec ${apptainer_image} celatlas_spatial rna barcode \
    --sample ${chip_id} --chemistry ${chemistry} --pattern ${chemistry_pattern} --thread ${n_cpus} \
    --mode strna --lowNum 2 --output_R1 --resume --gzip \
    --whitelist ${barcodes} \
    --fq1 ${sample_out_dir}/00.sample/${chip_id}_1.fq.gz \
    --fq2 ${sample_out_dir}/00.sample/${chip_id}_2.fq.gz \
    --outdir ${sample_out_dir}/01.barcode

  # Trimming adapter
  mkdir -p ${sample_out_dir}/02.cutadapt
  apptainer exec ${apptainer_image} celatlas_spatial rna cutadapt \
    --cutadapt_param '-g AAGCAGTGGTATCAACGCAGAGATC -g GGAAGAGCGTCGTGTAGGGAAAGAG -a AAAAAAAAAA -a TTTTTTTTTT -n 10 -q 20 --trim-n' \
    --sample ${chip_id} --thread ${n_cpus} \
    --overlap 10 --minimum_length 40 --nextseq_trim 0 --gzip \
    --insert ${insert_size}  --fq ${sample_out_dir}/01.barcode/${chip_id}_2.fq.gz \
    --outdir ${sample_out_dir}/02.cutadapt 

  # Mapping with STAR
  mkdir -p ${sample_out_dir}/03.star
  apptainer exec ${apptainer_image} celatlas_spatial rna star \
    --sample ${chip_id} --thread ${n_cpus} \
    --outFilterMultimapNmax 1 --starMem 30 \
    --genomeDir ${genome_dir} --fq ${sample_out_dir}/02.cutadapt/${chip_id}_clean_2.fq.gz \
    --outdir ${sample_out_dir}/03.star 

  # Count reads by feature counts
  mkdir -p ${sample_out_dir}/04.featureCounts
  apptainer exec ${apptainer_image} celatlas_spatial rna featureCounts \
    --sample ${chip_id} --thread ${n_cpus} --gtf_type ${feature_type} --genomeDir ${genome_dir} \
    --featureCounts_param '-s 1 ' \
    --input ${sample_out_dir}/03.star/${chip_id}_Aligned.sortedByCoord.out.bam \
    --outdir ${sample_out_dir}/04.featureCounts

  # Obtain count details
  mkdir -p ${sample_out_dir}/05.count
  apptainer exec ${apptainer_image} celatlas_spatial rna count \
    --sample ${chip_id} --thread ${n_cpus} --genomeDir ${genome_dir} \
    --expected_cell_num ${expected_cell_num} --cell_calling_method auto --force_cell_num None \
    --bam ${sample_out_dir}/04.featureCounts/${chip_id}_nameSorted.bam \
    --outdir ${sample_out_dir}/05.count

  # Bin segmentation
  mkdir -p ${sample_out_dir}/06.binSegment
  apptainer exec ${apptainer_image} celatlas_spatial rna binSegment \
    --sample ${chip_id} --thread ${n_cpus} --genomeDir ${genome_dir} \
    --model ${model_dir} \
    --segment --pixel-size ${pixel_size} --input ${chip_data_dir} --method ${segment_method} \
    --count --count_detail ${sample_out_dir}/05.count/${chip_id}_count_detail.txt \
    --outdir ${sample_out_dir}/06.binSegment 

  for bin_size in 10 20 50 100; do
    # Analysis, DEG etc.
    mkdir -p ${sample_out_dir}/07.analysis_bin${bin_size}
    apptainer exec ${apptainer_image} celatlas_spatial rna analysis \
      --sample ${chip_id} --thread ${n_cpus} --genomeDir ${genome_dir} \
      --square_bin_dir ${sample_out_dir}/06.binSegment/square_bin \
      --pixel-size ${pixel_size} --bin ${bin_size} \
      --outdir ${sample_out_dir}/07.analysis_bin${bin_size}

    mv -f ${sample_out_dir}/${chip_id}_report.html ${sample_out_dir}/${chip_id}_report_bin${bin_size}.html
  done

  touch ${sample_out_dir}/finished
done

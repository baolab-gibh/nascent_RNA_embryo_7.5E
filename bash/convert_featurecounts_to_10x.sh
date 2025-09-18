project_dir=${HOME}/documents/projects/wp_vasaseq

output_dir=${project_dir}/outputs/analysis/preprocessing/geo_seq/quantification/240409_Lib_embryo/featureCounts/

mkdir -p ${output_dir}

# Barcode list
grep -m 1 -e '^Geneid' ${count_tab} | cut -f9- | tr '\t' '\n' | xargs -I% bash -c 'basename $(dirname %)' > ${barcodes_file}

# Gene list
grep -v -e '^# Program' -e '^Geneid' ${count_tab} | awk -F$'\t' '{OFS="\\t"; print $8, $8}' > ${genes_file}

# Matrix
n_ft=$(wc -l ${genes_file} | cut -f1 -d' ')
n_bc=$(wc -l ${barcodes_file} | cut -f1 -d' ')
awk -v N_FT=${n_ft} -v N_BC=${n_bc} -f- <<'EOF' ${count_tab} > ${matrix_file}
  BEGIN  { OFS = "\\t"; print "%%MatrixMarket matrix coordinate integer general"; print "%"; print N_FT,N_BC,N_FT * N_BC }
  NR > 2 { OFS = "\\t"; i = 9; while (i <= NF) { print NR-2,i-8,$i; i++ } }
EOF

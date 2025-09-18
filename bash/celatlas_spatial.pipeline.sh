sample_name=ST110134_A1
chemistry=BBV2.4
chemistryPattern=C4L15C4L15C4U10T18
sampledir=/mnt/strna/work_project/pipeline/celatlas_spatial/baseST/$sample_name
rawdataDir=/mnt/strna/work_project/rawdata/celescope/$chemistry
segImageDir=/mnt/strna/work_project/rawdata/stomics/binSegment
modelDir=/home/zhoumy/lizt/code/work_space/spatial_bin/swin_tiny.pth
chipdataDir=/mnt/chip_mapping
tifDir=/mnt/strna/work_project/rawdata/celatlas_spatial/images/${sample_name}.tif
whitelistFile=/mnt/chip_mapping/ST_mask/${sample_name}.barcodeToPos.h5
genomeDir=/mnt/strna/work_project/rawdata/celatlas_spatial/reference/Homo_sapiens

ulimit -n 10240
bin=100
pixelSize=0.5
insertR2=150
cell_num=10000
feature_type=gene
mode=strna
method=gene_expr  # gene_expr or image
thread=128
tmp=yes


celatlas_spatial rna sample --outdir ${sampledir}/00.sample --sample ${sample_name} \
--thread ${thread} --chemistry ${chemistry}  --fq1 ${rawdataDir}/${sample_name}_1.fq.gz 


celatlas_spatial rna barcode --outdir ${sampledir}/01.barcode --sample ${sample_name} \
--thread ${thread} --chemistry ${chemistry} --pattern ${chemistryPattern} \
--whitelist ${whitelistFile} --mode ${mode} --lowNum 2 --gzip --output_R1 --resume \
--fq1 ${rawdataDir}/${sample_name}_1.fq.gz \
--fq2 ${rawdataDir}/${sample_name}_2.fq.gz


celatlas_spatial rna cutadapt --outdir ${sampledir}/02.cutadapt --sample ${sample_name} \
--thread ${thread} --minimum_length 20 --nextseq_trim 20 --gzip \
--overlap 10 --insert ${insertR2}  --fq ${sampledir}/01.barcode/${sample_name}_2.fq.gz
 

celatlas_spatial rna star --outdir ${sampledir}/03.star --sample ${sample_name} \
--thread ${thread} --genomeDir ${genomeDir} --outFilterMultimapNmax 1 \
--starMem 30  --fq ${sampledir}/02.cutadapt/${sample_name}_clean_2.fq.gz 


celatlas_spatial rna featureCounts --outdir ${sampledir}/04.featureCounts --sample ${sample_name} \
--thread ${thread} --gtf_type ${feature_type} --genomeDir ${genomeDir} --featureCounts_param '-s 1 ' \
--input ${sampledir}/03.star/${sample_name}_Aligned.sortedByCoord.out.bam


celatlas_spatial rna count --outdir ${sampledir}/05.count --sample ${sample_name} \
--thread ${thread} --genomeDir ${genomeDir} \
--expected_cell_num ${cell_num} --cell_calling_method auto \
--bam ${sampledir}/04.featureCounts/${sample_name}_nameSorted.bam --force_cell_num None


celatlas_spatial rna binSegment --outdir ${sampledir}/06.binSegment --sample ${sample_name} \
--thread ${thread} --genomeDir ${genomeDir} --pixel-size ${pixelSize} --input ${chipdataDir} \
--segment  --model ${modelDir} --method ${method} \
--count --count_detail ${sampledir}/05.count/${sample_name}_count_detail.txt


celatlas_spatial rna analysis --outdir ${sampledir}/07.analysis --sample ${sample_name} \
--thread ${thread} --genomeDir ${genomeDir} --square_bin_dir ${sampledir}/06.binSegment/square_bin \
--pixel-size ${pixelSize} --bin ${bin}


celatlas_spatial rna binSegment --outdir ${sampledir}/06.binSegment --sample ${sample_name} \
--thread ${thread} --genomeDir ${genomeDir} --pixel-size ${pixelSize} --input ${chipdataDir} \
--segment  --tif ${tifDir} --bs_out ${segImageDir} --model ${modelDir} --method ${method} \
--count --count_detail ${sampledir}/05.count/${sample_name}_count_detail.txt


# no-image
celatlas_spatial rna binSegment --outdir ${sampledir}/06.binSegment --sample ${sample_name} \
--thread ${thread} --genomeDir ${genomeDir} --pixel-size ${pixelSize} --input ${chipdataDir} \
--segment  --model ${modelDir} --method ${method} \
--count --count_detail ${sampledir}/05.count/${sample_name}_count_detail.txt


celatlas_spatial rna analysis --outdir ${sampledir}/07.analysis --sample ${sample_name} \
--thread ${thread} --genomeDir ${genomeDir} --square_bin_dir ${sampledir}/06.binSegment/square_bin \
--pixel-size ${pixelSize} --bin ${bin}

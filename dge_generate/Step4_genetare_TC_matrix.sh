#!/bin/bash

## This script excludes the locus which has T->C mutations in control sample

## Pre-install requirment: samtools, perl, Rscript

bam_dir=$1 # bam file folder including bam file which was generated in step3
sample=$2 # sample name with TFEA treatment
#num_barcodes=$3 # number of barcodes shown in the bam file
num_core_barcode=$3 # half of num_barcodes

############ 1. Extract the intermediate file
cd ${bam_dir}
perl /home/songjia/bigdisk/lsc/script/scripts/extract_digital_expression_matrix.pl ${sample}_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected.bam

############ 2. Generate T->C matrix with rds format
Rscript  /home/songjia/bigdisk/lsc/script/scripts/Generate_T_C_matrix.R ${sample}_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_UMI_read.txt ${num_core_barcode} ${sample}_TC_matrix.rds

### In the end you will get a rds files: $(sample)_TC_matrix.rds.

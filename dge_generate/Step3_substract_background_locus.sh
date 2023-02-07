#!/bin/bash

## This script excludes the locus which has T->C mutations in control sample

## Pre-install requirment: sam2tsv, samtools, perl

#tsv_dir=$1 # tsv file folder including tsv files which was generated in step2
sample=$1 # sample name with TFEA treatment
control=$2 # sample name without TFEA treatment
bam_dir=$3 # bam file folder including bam file which was generated in step1
tsv_dir=$4 #sample tsv file folder including bam file which was generated in step2  
tsv_dir2=$5 #control tsv file dir
############ 1. Substract mutations from control samples
#cd ${tsv_dir}
perl /home/songjia/bigdisk/lsc/script/scripts/background_correction.pl -bg ${tsv_dir2}/${control}_both_strand_all_TC.tsv_q27.tsv -in ${tsv_dir}/${sample}_both_strand_all_TC.tsv_q27.tsv

############ 2. Add mutation information back to bam files
perl /home/songjia/bigdisk/lsc/script/scripts/TagIntronicRead_V5.pl -read ${sample}_both_strand_all_TC.tsv_q27.tsv_corrected.tsv -bam ${bam_dir}/${sample}_star_gene_exon_tagged_TagIntronic_clean.bam

### In the end you will get a bam file with "GE:Z:genename--T" or "GE:Z:genename--C" tags. 
### The file format is "${sample}_starAligned.sorted.merged.GeneExonTagged.TagIntronic.clean.${num_barcodes}.TagTC.corrected.bam"




sample=$1
samtools view -S -q 20 ${sample}.sam > ${sample}_MAPQ20.sam
sam2tsv --reference Homo_sapiens.GRCh38.dna.primary_assembly.fa ${sample}_MAPQ20.sam | awk '{if ($9 ~/M|=|X/ && $2==0 && $8=="T") print $0; else if ($9 ~/M|=|X/ && $2==16 && $8=="A") print $0}'  > ${sample}_MAPQ20_both_strand_available_site.tsv
Rscript GetCellBarcode.R ${sample}_TC_matrix.rds
python ExtractMappingInfo.py ${sample}_MAPQ20.sam ${sample}_cell_filter_name.txt
Rscript Count_bycell_cutoff.R ${sample}
Rscript FilterRead.R ${sample}
sort -R ${sample}_filter_count.txt | head -n100000 > ${sample}_filter_sampling.txt
Rscript likehood_theta_for_gene_100times.R ${sample}
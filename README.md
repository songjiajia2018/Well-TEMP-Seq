# Well-TEMP-Seq

### Prerequisites

- Drop-seq 2.5.1
- R 4.1.3
- Pysam
- MASS
- nloptr

----

### Well-TEMP-seq TC calling pipeline includes following steps:

- Step1_Drop_seq.py: utilize the Drop-seq computational pipeline (James Nemesh, McCarroll Lab, version 1.12; Macosko et al., 2015) to map the reads to the genome and tag the reads with cell barcode, UMI barcode and gene annotation in bam files. Next, we extracted intronic reads in bam file because the legacy Drop-seq computational pipeline (version 1.12) only consider exonic reads.

  #### Below steps are referred to scNT pipeline.

- Step2_extract_alignment_info.sh: sam2tsv (https://github.com/lindenb/jvarkit/; version ec2c2364) is used to extract detailed alignment information from bam files and then T-to-C substitutions are identified in both experimental and control samples (without Timelapse chemical conversion reaction, as a control for background mutations).
- Step3_substract_background_locus.sh: exclude the genomic sites with background T-to-C substitutions from the downstream analysis.
- Step4_genetare_TC_matrix.sh: generate labeled and unlabeled gene expression matrix.

### Statistical correction

To address the insufficiency of metabolic RNA labeling, we adopted a statistic model to approximate the real distribution of T-to-C substitution in single cells and estimate the real mutation rate.

### Data

- Raw data files are available at NCBI Gene Expression Omnibus (GEO) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194357).

- Simple demos can be found in `/Well-TEMP-Seq/data`. 2500 entries are extracted from raw fastq files as demo. 400 cells are extracted to demostrate RNA velocity analysis.


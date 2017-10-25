# bulk_RNA_seq_count_base 

A bioinformatic analysis pipeline for bulk RNA-seq differential analysis through count base strtegy

Note: If you already have mapped files, you can start from the Step 4. for sorting

1. Format preprocessing

- Files:
  - bamToFastq.qsub
  - bamToFastq.sh

- Description:
  - These two files change the origin bam format inputs into fastq format. Fastq format can accelerate our following analysis because many tools take fastq format as their default format.

2. FastQC

- Files:
  - fastqc.sh
  - fastqc.qsub

- Description:
  - This shell script help us understand the sequencing qualities of each samples and generate well-organized visualization of qualities

3. HISAT for mapping

- Files:
  - hisat2_mapping.qsub
  - hisat2_mapping.sh

- Description:
  - These files utilize parallel computation to execute  HISAT program. In this script, it will take GRCh38 annotation file to conduct reads mapping with our RNA-seq fastq files.

4. Sorting by name
- Files:
  - sort_name_forCount.sh
  - sort_name_forCount.qsub

- Description:
  - This step help us generate sorted files for the following HTSeq counting

5. HTSeq counting
- Files:
  - htseq_count.sh
  - htseq_count.qsub

- Description:
  - In this step, we use HTSeq program to generate count matrix, which is crucial in this pipeline to understand differential analysis in the following programs

6. Create assemblies.txt
- Files:
  - DESeq2_multiple_groups.R
  
- Description:
  - In this R script, we perform differential analysis and generate multiple figures for data visualization



# Resistome-quantification
__BEFORE YOU START, YOU SHOULD REPLACE THE "***", WHICH LIES AT THE TOP OF THE SCRIPT, WITH THE ACCESSION NUMBER OF THE METAGENOMIC FILE ON THE NCBI SRA DATABASE.__

Your system should also have some software and tools downloaded such as SRA ToolKit, fastq-dump, reformat from the BBmap suite, Bowtie2, Samtools, Bedtools, and Diamond. In addition, the CARD database is a prerequisite, or any other database you prefer to use as your reference.

1. The pipeline downloads raw metagenomic data from NCBI SRA using prefetch from the SRA Toolkit, this gives a raw file in SRA format.
2. It converts the SRA file to FASTQ format using fastq-dump. The output is generated in fastq folder, for paired-end fastq files, R1 and R2 files will be generated.
3. Use reformat from the BBmap suite to interleave the generated forward and reverse fastq files, R_1 and R_2, and generate a single interleaved fastq file for further analysis. (This might make things easier, it is optional, though! You can go without it if you're willing to adjust what follows accordingly).
4. Reads alignment with Bowtie2 to a reference database containing antimicrobial resistance gene sequences.
5. Convert the SAM file to BAM format and sort it for downstream analysis.
6. Convert BAM (the Binary format) to FASTQ format.
7. Quantify the resistome using DIAMOND and the data-driven scripting language, AWK.
   

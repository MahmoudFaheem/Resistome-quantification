# Resistome Quantification Pipeline

## Overview

Before you start, you should replace the "***", which lies at the top of the script, with the accession number of the metagenomic file on the NCBI SRA database.

Your system should also have some software and tools downloaded such as SRA ToolKit, fastq-dump, reformat from the BBmap suite, Bowtie2, Samtools, Bedtools, and Diamond. In addition, the CARD database is a prerequisite, or any other database you prefer to use as your reference.

## Pipeline Steps

1. **Download Raw Metagenomic Data**
    - This step downloads raw metagenomic data from NCBI SRA using prefetch from the SRA Toolkit, which gives a raw file in SRA format.

2. **Convert SRA to FASTQ**
    - Convert the SRA file to FASTQ format using fastq-dump. The output is generated in a fastq folder, with R1 and R2 files for paired-end fastq files.

3. **Interleave FASTQ Files**
    - Use reformat from the BBmap suite to interleave the generated forward and reverse fastq files, R_1 and R_2, and generate a single interleaved fastq file for further analysis. (This step is optional but might make downstream analysis easier.)

4. **Reads Alignment with Bowtie2**
    - Align reads to a reference database containing antimicrobial resistance gene sequences using Bowtie2.

5. **Convert SAM to BAM**
    - Convert the SAM file to BAM format and sort it for downstream analysis.

6. **Convert BAM to FASTQ**
    - Convert BAM (the binary format) to FASTQ format.

7. **Quantify Resistome with DIAMOND**
    - Quantify the resistome using DIAMOND and the data-driven scripting language, AWK.

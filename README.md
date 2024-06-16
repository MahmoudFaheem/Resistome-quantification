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

## How to Use

Follow the steps below to use the Resistome Quantification Pipeline:

### Step 1: Clone Repository

Clone the repository to your local machine and navigate to the project directory:

```bash
git clone https://github.com/mahmoudfaheem/Resistome-quantification.git
cd Resistome-quantification



### Step 2: Update SRA ID

Replace "***" with your actual SRA ID in the metagenomics.yml file.



### Step 3: Commit Changes
Commit the changes to update the SRA ID:
git add .
git commit -m "Update SRA ID"
git push origin main



### Step 4: Trigger Workflow Manually

Run the following command in your terminal to trigger the workflow manually:
curl -X POST -u mahmoudfaheem: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 \
  -H "Accept: application/vnd.github.v3+json" \
  https://api.github.com/repos/mahmoudfaheem/Resistome-quantification/actions/workflows/metagenomics.yml/dispatches \
  -d '{"ref":"main"}'


### Make sure to replace "***" with your actual SRA ID in the metagenomics.yml file.

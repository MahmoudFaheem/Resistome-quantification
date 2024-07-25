#!/bin/bash

# BEFORE YOU START, REPLACE THE "***" WITH THE ACCESSION NUMBER OF THE METAGENOMIC FILE
srr_id="***"
R1="${srr_id}_pass_1.fastq.gz"
R2="${srr_id}_pass_2.fastq.gz"

# 1. Download raw metagenomic data from NCBI SRA using prefetch from the SRA Toolkit
echo "Downloading raw metagenomic data for ${srr_id}"
prefetch ${srr_id}

# 2. Convert the SRA file to FASTQ format using fastq-dump
echo "Converting SRA to FASTQ"
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ~/${srr_id}/${srr_id}.sra

# 4. Interleave the generated forward and reverse fastq files using reformat from the BBmap suite
echo "Interleaving FASTQ files"
reformat.sh in1=~/fastq/${R1} in2=~/fastq/${R2} out=interleaved${srr_id}.fq.gz

# 5. Align reads with Bowtie2 to a reference database
echo "Aligning reads with Bowtie2"
bowtie2 -x reference_index -U ~/interleaved${srr_id}.fq.gz -S ~/aligned${srr_id}.sam

# 6. Convert the SAM file to BAM format and sort it
echo "Converting SAM to BAM and sorting"
samtools view -bS ~/aligned${srr_id}.sam | samtools sort -o ~/sorted${srr_id}.bam
samtools index ~/sorted${srr_id}.bam

# 7. Convert BAM to FASTQ
echo "Converting BAM to FASTQ"
bedtools bamtofastq -i sorted${srr_id}.bam -fq sorted${srr_id}.fastq.gz

# 8. Extract gene names and generate counts using diamond and awk
echo "Extracting gene names and generating counts"
blast_results=~/diamond_${srr_id}
AWK_results=~/awk_diamond_${srr_id}

diamond blastx -d ~/CARD -q ~/sorted${srr_id}.fastq.gz -o ${blast_results} -k 1 -f 6 -b 1
awk 'BEGIN { OFS="\t"; print "Query ID\tSubject ID\tPercentage of identical matches\tAlignment length\tNumber of mismatches\tNumber of gap openings\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject\tExpected value\tBit score\tGene_Count" } { gsub(/.*ARO_Name:/, "", $2); gsub(/\|.*/, "", $2); gene_count[$2]++ } END { for (gene in gene_count) { print gene, gene_count[gene]; total_count += gene_count[gene] } print "Total_Count\t\t\t\t\t\t\t\t\t\t\t\t\t", total_count > "gene_counts.txt"; for (gene in gene_count) { print gene, gene_count[gene] > "unique_gene_sums.txt" } } { print $0 "\t" gene_count[$2] }' "${blast_results}" > "${AWK_results}"

echo "Pipeline completed successfully."

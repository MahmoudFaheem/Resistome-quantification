# BEFORE YOU START, REPLACE THE "***" WITH THE ACCESSION NUMBER OF THE METAGENOMIC FILE

#1. Download raw metagenomic data from NCBI SRA using prefetch from the SRA Toolkit
srr_id="***"
R1="${srr_id}_pass_1.fastq.gz"
R2="${srr_id}_pass_2.fastq.gz"

prefetch ${srr_id}

#2. Convert the SRA file to Fastq format using fastq-dump. The output is generated in a fastq folder, for paired-end fastq files, R1 and R2 files will be generated.
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ~/${srr_id}/${srr_id}.sra

#4. Use reformat from the BBmap suite to interleave the generated forward and reverse fastq files, R_1 and R_2, and generates a single interleaved fastq file for further analysis.
reformat.sh in1=~/fastq/${R1} in2=~/fastq/${R2} out=interleaved${srr_id}.fq.gz

#5. Read Alignment with Bowtie2 to a reference database containing AMR gene sequences
bowtie2 -x reference_index -U ~/interleaved${srr_id}.fq.gz -S ~/aligned${srr_id}.sam

#6. Convert the SAM file to BAM format and sort it for downstream analysis
samtools view -bS ~/aligned${srr_id}.sam | samtools sort -o ~/sorted${srr_id}.bam
samtools index ~/sorted${srr_id}.bam

#7. Convert BAM to FASTQ
bedtools bamtofastq -i sorted${srr_id}.bam -fq sorted${srr_id}.fastq.gz

#8. A single-line script to extract the gene name only and make the header tab-separated, it also creates a new column called "Gene_Count" and assigns a value of +1 to each gene every time it appears. It provides four files, diamond, awk_diamond, gene_sums [showing the number of occurrences for each unique gene], and gene_counts [which is just the total number of all genes combined]:
blast_results=~/diamond_${srr_id}
AWK_results=~/awk_diamond_${srr_id}

diamond blastx -d ~/CARD -q ~/sorted${srr_id}.fastq.gz -o diamond_${srr_id} -k 1 -f 6 -b 1 && awk 'BEGIN { OFS="\t"; print "Query ID\tSubject ID\tPercentage of identical matches\tAlignment length\tNumber of mismatches\tNumber of gap openings\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject\tExpected value\tBit score\tGene_Count" } { gsub(/.*ARO_Name:/, "", $2); gsub(/\|.*/, "", $2); gene_count[$2]++ } END { for (gene in gene_count) { print gene, gene_count[gene]; total_count += gene_count[gene] } print "Total_Count\t\t\t\t\t\t\t\t\t\t\t\t\t", total_count > "gene_counts.txt"; for (gene in gene_count) { print gene, gene_count[gene] > "unique_gene_sums.txt" } } { print $0 "\t" gene_count[$2] }' "$blast_results" > "$AWK_results"




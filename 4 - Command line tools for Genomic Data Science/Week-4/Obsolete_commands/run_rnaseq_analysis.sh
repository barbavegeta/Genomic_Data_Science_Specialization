#!/bin/bash

# =============================================================================
#
# Modern RNA-Seq Analysis Pipeline Script
# Coursera: Johns Hopkins University - Command Line Tools for Genomic Data Science
#
# This script converts the original Module 4 workflow from the outdated
# Tuxedo suite (TopHat/Cufflinks) to the modern HISAT2/StringTie pipeline.
#
# =============================================================================
#
# SCRIPT SETUP
#
# Exit script if any command fails (-e) or if any variable is unset (-u)
set -eu

# ASSUMPTIONS:
# 1. This script is run in a directory containing the necessary input files:
#    - athal_chr.fa (reference genome)
#    - athal_genes.gtf (gene annotation)
#    - Day8.fastq (RNA-seq reads sample 1)
#    - Day16.fastq (RNA-seq reads sample 2)
# 2. All required tools are installed and available in the system's PATH:
#    - hisat2, hisat2-build
#    - samtools
#    - stringtie
#    - gffcompare
# 3. The 'prepDE.py' script (from the StringTie developers) is downloaded and
#    accessible in your PATH or current directory.

echo "Starting the modern RNA-Seq analysis pipeline..."


# =============================================================================
# STAGE 1: GENOME INDEXING
# We use 'hisat2-build' to create a genome index for the HISAT2 aligner.
# =============================================================================
echo -e "\n[STAGE 1] Building HISAT2 genome index..."

mkdir -p athal_index
hisat2-build athal_chr.fa athal_index/athal

echo "Genome indexing complete."


# =============================================================================
# STAGE 2: ALIGNMENT
# We align the RNA-seq reads from both samples to the indexed genome using
# HISAT2. The output is piped to 'samtools' to create a sorted BAM file,
# which is the standard for downstream analysis.
# =============================================================================
echo -e "\n[STAGE 2] Aligning reads to the genome..."

mkdir -p Tophat/Day8
mkdir -p Tophat/Day16

# Align Day 8 sample
hisat2 -x athal_index/athal -U Day8.fastq 2> Tophat/Day8/align_summary.txt | samtools sort -o Tophat/Day8/accepted_hits.bam

# Align Day 16 sample
hisat2 -x athal_index/athal -U Day16.fastq 2> Tophat/Day16/align_summary.txt | samtools sort -o Tophat/Day16/accepted_hits.bam

echo "Alignment complete. Sorted BAM files are in Tophat/Day*/"


# =============================================================================
# STAGE 3: TRANSCRIPT ASSEMBLY
# We use StringTie to assemble transcripts from the aligned reads for each
# sample. This is done in a reference-guided manner using the '-G' flag.
# =============================================================================
echo -e "\n[STAGE 3] Assembling transcripts with StringTie..."

mkdir -p Cufflinks/Day8
mkdir -p Cufflinks/Day16

# Assemble transcripts for Day 8
stringtie Tophat/Day8/accepted_hits.bam -G athal_genes.gtf -o Cufflinks/Day8/transcripts.gtf -l Day8

# Assemble transcripts for Day 16
stringtie Tophat/Day16/accepted_hits.bam -G athal_genes.gtf -o Cufflinks/Day16/transcripts.gtf -l Day16

echo "Transcript assembly complete. Assembled GTFs are in Cufflinks/Day*/"


# =============================================================================
# STAGE 4: COMPARISON WITH REFERENCE ANNOTATION
# We use 'gffcompare' (the successor to cuffcompare) to compare our
# assembled transcripts against the reference annotation.
# =============================================================================
echo -e "\n[STAGE 4] Comparing assembled transcripts with GffCompare..."

# Note: gffcompare is run inside each directory to generate output files there.
(cd Cufflinks/Day8 && gffcompare -r ../../athal_genes.gtf -o cuffcmp transcripts.gtf)
(cd Cufflinks/Day16 && gffcompare -r ../../athal_genes.gtf -o cuffcmp transcripts.gtf)

echo "Comparison complete. GffCompare output is in Cufflinks/Day*/"


# =============================================================================
# STAGE 5: DIFFERENTIAL EXPRESSION ANALYSIS PREPARATION
# The modern replacement for cuffdiff involves multiple steps.
#
# STEP 5.1: Merge Assemblies
# Create a unified, non-redundant set of transcripts from all samples.
# =============================================================================
echo -e "\n[STAGE 5.1] Merging transcript assemblies..."

ls -1 Cufflinks/Day*/transcripts.gtf > GTFs.txt
stringtie --merge -G athal_genes.gtf -o merged.gtf GTFs.txt

echo "Merged GTF created: merged.gtf"


# =============================================================================
# STEP 5.2: Re-quantify Abundances
# Re-run StringTie on each sample, but this time using the merged GTF. The
# '-e' flag restricts quantification to only the transcripts in the merged
# model, ensuring a consistent basis for comparison.
# =============================================================================
echo -e "\n[STAGE 5.2] Re-quantifying abundances against merged model..."

# Create subdirectories for each sample within Cuffdiff/
mkdir -p Cuffdiff/Day8
mkdir -p Cuffdiff/Day16

# stringtie output now goes into these sample-specific subdirectories
stringtie -e -B -G merged.gtf -o Cuffdiff/Day8/transcripts.gtf Tophat/Day8/accepted_hits.bam
stringtie -e -B -G merged.gtf -o Cuffdiff/Day16/transcripts.gtf Tophat/Day16/accepted_hits.bam

echo "Quantification complete. GTF files are now in Cuffdiff/DayX/transcripts.gtf"


# =============================================================================
# STEP 5.3: Generate Count Matrices
# Use the 'prepDE.py' script to convert the StringTie output into gene and
# transcript count matrices. These matrices are the standard input for R
# packages like DESeq2 and edgeR.
# =============================================================================
echo -e "\n[STAGE 5.3] Generating count matrices for R..."

# Remove or comment out these lines as they are no longer needed for this method:
# CURRENT_DIR=$(pwd)
# echo -e "Day8\t${CURRENT_DIR}/Cuffdiff/Day8.gtf" > cuffdiff_gtfs.txt
# echo -e "Day16\t${CURRENT_DIR}/Cuffdiff/Day16.gtf" >> cuffdiff_gtfs.txt

# New prepDE.py command: Use -i with the parent directory.
# This assumes './prepDE.py' is the corrected Python 3 version from earlier steps.
./prepDE.py -i Cuffdiff/
echo "Count matrices 'gene_count_matrix.csv' and 'transcript_count_matrix.csv' created."
echo "These files are ready for differential expression analysis in R."

# =============================================================================
# SCRIPT FINISHED
# =============================================================================
echo -e "\nAnalysis pipeline complete."
echo "You can now use the generated count matrices for statistical analysis in R."
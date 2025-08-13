#!/bin/bash

#########################################
###    Download and Prepare FASTQ     ###
#########################################

# Download .sra file from NCBI using SRA Toolkit
prefetch SRR1107997

# Check that the file was downloaded (it's binary, so head may not look readable)
head SRR1107997.sra

# Convert .sra to .fastq in the background, ignoring hangups (nohup)
nohup fastq-dump SRR1107997.sra &

# Check if FASTQ file was created and view first few lines
head SRR1107997.fastq

# Count total lines in a compressed FASTQ
zcat NA12814_1.fastq.gz | wc -l  # Divide by 4 to get number of reads


########################################
###        SAMTOOLS Essentials       ###
########################################

# Get alignment statistics from a BAM file
samtools flagstat NA12814.bam

# Sort and index BAM
samtools sort NA12814.bam -o NA12814_sorted.bam
samtools index NA12814_sorted.bam

# Merge multiple BAM files
samtools merge merged.bam sample1.bam sample2.bam

# Convert BAM to SAM
samtools view -h NA12814.bam > NA12814.sam

# Convert SAM to BAM
samtools view -Sb NA12814.sam > NA12814_from_sam.bam

# Convert BAM to FASTQ
samtools fastq NA12814.bam -o NA12814.fastq


##########################################
###   GTF/BED to BAM and FASTA Tools   ###
##########################################

# Convert BAM to BED
bedtools bamtobed -i NA12814.bam > NA12814.bed

# Convert BED to BAM
bedtools bedtobam -i RefSeq.bed -g hg38c.hdrs > refseq_from_bed.bam
samtools view refseq_from_bed.bam | more

# Convert GTF to BAM
bedtools bedtobam -i RefSeq.gtf -g hg38c.hdrs > refseq_from_gtf.bam
samtools view refseq_from_gtf.bam | more

# Extract FASTA from GTF
bedtools getfasta -fi hg38c.fa -bed RefSeq.gtf -fo RefSeq.gtf.fasta

# Extract FASTA from BED
bedtools getfasta -fi hg38c.fa -bed RefSeq.bed -fo RefSeq.bed.fasta

# Handle spliced entries using -split
bedtools getfasta -split -fi hg38c.fa -bed RefSeq.bed -fo RefSeq.bed.spliced.fasta


###########################################
###    SAM/BAM Exploration Questions    ###
###########################################

# Q1: Number of alignments in full BAM
samtools flagstat athal_wu_0_A.bam
samtools view athal_wu_0_A.bam | wc -l

# Q2: Alignments with mate unmapped
samtools view athal_wu_0_A.bam | cut -f7 | grep -c '*'

# Q3: Alignments containing a deletion (D)
samtools view athal_wu_0_A.bam | cut -f6 | grep -c 'D'

# Q4: Alignments with mate mapped to same chromosome
samtools view athal_wu_0_A.bam | cut -f7 | grep -c '='

# Q5: Spliced alignments (contain 'N')
samtools view athal_wu_0_A.bam | cut -f6 | grep -c 'N'


##########################################
###   Subset to Region and Re-analyze  ###
##########################################

# Sort and index original BAM
samtools sort -o athal_wu_0_A.sorted.bam athal_wu_0_A.bam
samtools index athal_wu_0_A.sorted.bam

# Extract region of interest
samtools view -b athal_wu_0_A.sorted.bam 'Chr3:11777000-11794000' > athal_wu_0_A.region.bam

# Q6-Q10 repeated for region
samtools flagstat athal_wu_0_A.region.bam
samtools view athal_wu_0_A.region.bam | cut -f7 | grep -c '*'
samtools view athal_wu_0_A.region.bam | cut -f6 | grep -c 'D'
samtools view athal_wu_0_A.region.bam | cut -f7 | grep -c '='
samtools view athal_wu_0_A.region.bam | cut -f6 | grep -c 'N'

Column 1 (QNAME): Read Name
Column 2 (FLAG): Bitwise flag with alignment information
Column 3 (RNAME): Chromosome/Reference name
Column 4 (POS): 1-based start position
Column 5 (MAPQ): Mapping Quality
Column 6 (CIGAR): CIGAR string (describes insertions, deletions, etc.)
  M: Alignment match (can be a sequence match or mismatch)
  I: Insertion into the reference
  D: Deletion from the reference
  S: Soft clipping (bases at the end of the read not aligned)
  H: Hard clipping (bases at the end of the read not present in the sequence)
  N: Skipped region from the reference (used for spliced RNA-Seq alignments)
  =: Sequence match
  X: Sequence mismatch
Column 7 (RNEXT): Mate's chromosome/reference name ("=" if same, "*" if unmapped)

So, you use -f7 to specifically isolate information about the mate's mapping location.
########################################
###      BAM Header and Metadata     ###
########################################

# Q11: Number of genome sequences
samtools view -H athal_wu_0_A.bam | grep -c 'SN:'

# Q12: Length of first sequence
samtools view -H athal_wu_0_A.bam | grep 'SN:' | more

# Q13: Alignment tool used
samtools view -H athal_wu_0_A.bam | grep -i 'PG'

# Q14: Read identifier of first alignment
samtools view athal_wu_0_A.bam | head -1 | cut -f1

# Q15: Start position of read's mate
samtools view athal_wu_0_A.bam | head -n 1 | awk '{print $7":"$8}'


########################################
###        Exon Overlap Analysis     ###
########################################

# Create BED from region BAM
bedtools bamtobed -i athal_wu_0_A.region.bam > region.bed

# Extract exons from GTF to BED format
grep -P "\texon\t" athal_wu_0_A_annot.gtf | awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5}' > exons.bed

# Q16: Number of overlaps
bedtools intersect -abam athal_wu_0_A.region.bam -b athal_wu_0_A_annot.gtf -bed -wo > overlaps.bed
wc -l overlaps.bed

# Q17: Overlaps >=10 bases
cut -f22 overlaps.bed | sort -nrk1 > lengths
cut -f22 overlaps.bed | sort -nrk1 | grep -n '^9' | head -1

# Q18: Number of unique alignments that overlap annotations
cut -f1-12 overlaps.bed | sort -u | wc -l

# Q19: Number of unique exons with mapped reads
cut -f13-21 overlaps.bed | sort -u | wc -l

# Q20: BED records from GTF (count unique transcript entries)
cut -f9 athal_wu_0_A_annot.gtf | cut -d ' ' -f1,2 | sort -u | wc -l


########################################
###       Attribute Field Parsing    ###
########################################

# Intersect RefSeq.gtf and Alus.bed and extract attribute field
tail -n +1 intersect_output.txt | cut -f9 | cut -d ' ' -f2 | sort -u | wc -l

# View extracted field to check correctness
cat intersect_output.txt | cut -f9 | head -n 10

: '
########## Here are the first few standard SAM columns for context: ##########




########## BEDTools Intersect Output (cut -f1-12 and cut -f13-21): ##########

This is different because the file overlaps.bed was created with:
bedtools intersect -abam [alignments.bam] -b [annotations.gtf] -bed -wo > overlaps.bed

- The -wo option (write overlap) appends the number of overlapping bases.
- The -bed option makes sure BAM alignments are converted to BED12 format.

Therefore:
- Columns 1–12: Alignment Record (from the BAM/BED12 input)
- Columns 13–21: Exon Record (from the GTF input)
- Column 22: Overlap length (number of overlapping bases)

To explore the structure of overlaps.bed:


########## Step 1: Look at the first line (will be long!) ##########
head -n 1 overlaps.bed

# Check if columns 1-12 look like BED12 alignments
cut -f1-12 overlaps.bed | head -n 5

# Expected in BED12:
#   Column 1: Chromosome (e.g., Chr3)
#   Columns 2-3: Start/End
#   Column 4: Read name (QNAME)
#   Column 6: Strand (+/-)


########## Step 2: Check columns 13-21 for GTF structure ##########
cut -f13-21 overlaps.bed | head -n 5

# Compare with original GTF file
head athal_wu_0_A_annot.gtf

# In overlaps.bed, Column 15 (original GTF col 3) should be "exon"
# Column 21 should contain attributes like gene_id and transcript_id


########## Step 3: Confirm column 22 shows overlap length ##########
cut -f22 overlaps.bed | head -n 5

# Expect: numeric values showing base-pair overlap between read and exon

Golden Rule:
1. Understand the tool (e.g., bedtools intersect -wo)
2. Hypothesize the structure (e.g., 12 + 9 + 1 = 22 columns)
3. Verify with cut and head
4. Compare to source files

This method gives you confidence that you're selecting the right fields—based not on guesswork, but on data structure evidence.
'

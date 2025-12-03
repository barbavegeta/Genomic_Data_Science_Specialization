Task                    |	Obsolete Tool|	Modern Tool
Spliced read alignment  |	TopHat       |	HISAT2
Transcript assembly     |	Cufflinks    |	StringTie
Transcript comparison   |	Cuffcompare	 | gffcompare
Differential expression |	Cuffdiff     |	Ballgown / DESeq2


TopHat + Cufflinks pipeline into a modern HISAT2 + StringTie version

Modern RNA-seq Pipeline for Genomic Data Science Exam
Key Modern Tools:

Aligner: STAR (Spliced Transcripts Alignment to a Reference)
Transcript Assembler/Quantifier: StringTie
Differential Expression: DESeq2 (or EdgeR)

# Modern RNA-Seq Workflow (for all 35 questions)

### 1. Build HISAT2 index of the reference genome

hisat2-build athal_chr.fa athal_hisat2_index

### 2. Align reads with HISAT2

hisat2 -x athal_hisat2_index -U Day8.fastq -S Day8.sam
hisat2 -x athal_hisat2_index -U Day16.fastq -S Day16.sam

# Convert SAM to sorted BAM:

samtools view -bS Day8.sam | samtools sort -o Day8.sorted.bam
samtools view -bS Day16.sam | samtools sort -o Day16.sorted.bam
samtools index Day8.sorted.bam
samtools index Day16.sorted.bam

### 3. Assemble transcripts with StringTie

stringtie Day8.sorted.bam -G athal_genes.gtf -o Day8.gtf -l Day8
stringtie Day16.sorted.bam -G athal_genes.gtf -o Day16.gtf -l Day16

### 4. Merge assemblies

# Create a `mergelist.txt`:

touch mergelist.txt

Day8.gtf
Day16.gtf

# Then:

stringtie --merge -G athal_genes.gtf -o merged.gtf mergelist.txt

### 5. Quantify transcripts on merged annotation

stringtie Day8.sorted.bam -e -B -G merged.gtf -o Day8_merged.gtf
stringtie Day16.sorted.bam -e -B -G merged.gtf -o Day16_merged.gtf


# This generates expression files for Ballgown or DESeq2 downstream.

### 6. Differential expression analysis

# Using prepDE.py from StringTie or tools like DESeq2 in R on the count tables to find differentially expressed genes and transcripts.
# How to answer the questions with this pipeline
# Example commands for key extractions:

# Count reads in BAM (total alignments)
samtools view -c Day8.sorted.bam
samtools view -c Day16.sorted.bam

# Count mapped reads
samtools view -c -F 4 Day8.sorted.bam
samtools view -c -F 4 Day16.sorted.bam

# Count unmapped reads
samtools view -c -f 4 Day8.sorted.bam
samtools view -c -f 4 Day16.sorted.bam

# Count unique transcripts in GTF
grep -P "\ttranscript\t" Day8.gtf | wc -l
grep -P "\ttranscript\t" Day16.gtf | wc -l

# Count genes in GTF (assuming gene_id attribute present)
grep -Po 'gene_id "[^"]+"' Day8.gtf | sort | uniq | wc -l

# Use gffcompare to classify transcripts
gffcompare -r athal_genes.gtf -o cmp Day8.gtf
# Check cmp.annotated.gtf, cmp.stats, cmp.refmap for detailed stats

# For differential expression, parse DESeq2 R results (pseudo-code):
# summary(res) where res is DESeq2 output dataframe
nrow(subset(res, padj < 0.05))

######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################

# Step 1: Build HISAT2 index for Arabidopsis genome

hisat2-build athal_chr.fa athal_hisat2_index

# Step 2: Align Day8 and Day16 RNA-seq reads with HISAT2

hisat2 -x athal_hisat2_index -U Day8.fastq -S Day8.sam 2> Day8_hisat2.log
hisat2 -x athal_hisat2_index -U Day16.fastq -S Day16.sam 2> Day16_hisat2.log

# Step 3: Convert SAM to sorted BAM and index

samtools view -bS Day8.sam | samtools sort -o Day8.sorted.bam
samtools view -bS Day16.sam | samtools sort -o Day16.sorted.bam
samtools index Day8.sorted.bam
samtools index Day16.sorted.bam

# Step 4: Assemble transcripts with StringTie

stringtie Day8.sorted.bam -G athal_genes.gtf -o Day8.gtf -l Day8
stringtie Day16.sorted.bam -G athal_genes.gtf -o Day16.gtf -l Day16


# Step 5: Merge assemblies for comparison

# Create `mergelist.txt`:

Day8.gtf
Day16.gtf

# Merge:

stringtie --merge -G athal_genes.gtf -o merged.gtf mergelist.txt

# Step 6: Quantify transcripts with merged annotation
stringtie Day8.sorted.bam -e -B -G merged.gtf -o Day8_merged.gtf
stringtie Day16.sorted.bam -e -B -G merged.gtf -o Day16_merged.gtf

# Step 7: Run gffcompare to classify transcripts vs reference
gffcompare -r athal_genes.gtf -o cmp Day8.gtf
gffcompare -r athal_genes.gtf -o cmp Day16.gtf

# Step 8: Extract answers from logs, BAM files, GTF, and gffcompare output

#######################################################################################################################################
#######################################################################################################################################

# Extraction commands for answers

### Q1, Q2: Number of alignments (total reads aligned)

samtools view -c Day8.sorted.bam  # Q1 answer
samtools view -c Day16.sorted.bam # Q2 answer

### Q3, Q4: Number of reads mapped (from HISAT2 logs or BAM)
# Count mapped:

samtools view -c -F 4 Day8.sorted.bam  # Q3
samtools view -c -F 4 Day16.sorted.bam # Q4

### Q5, Q6: Number uniquely aligned reads
# From HISAT2 log files:

grep "aligned exactly 1 time" Day8_hisat2.log  # Q5
grep "aligned exactly 1 time" Day16_hisat2.log # Q6

### Q7, Q8: Number of spliced alignments (junction reads
# In BAM, spliced reads have `N` in the CIGAR string:

samtools view Day8.sorted.bam | grep -c 'N'  # Q7
samtools view Day16.sorted.bam | grep -c 'N' # Q8

### Q9, Q10: Reads left unmapped
samtools view -c -f 4 Day8.sorted.bam  # Q9
samtools view -c -f 4 Day16.sorted.bam # Q10

### Q11, Q12: Number of genes assembled by StringTie
# Count unique gene\_ids in Day8.gtf and Day16.gtf:

grep 'gene_id "' Day8.gtf | sed 's/.*gene_id "\([^"]*\)".*/\1/' | sort | uniq | wc -l  # Q11
grep 'gene_id "' Day16.gtf | sed 's/.*gene_id "\([^"]*\)".*/\1/' | sort | uniq | wc -l # Q12

### Q13, Q14: Number of transcripts assembled
# Count transcripts:

grep 'transcript_id' Day8.gtf | wc -l  # Q13
grep 'transcript_id' Day16.gtf | wc -l # Q14

### Q15, Q16: Number of single transcript genes
# Count genes with exactly one transcript:

# Q15
grep -w 'transcript' Day8.gtf | cut -f9 |  grep -o 'gene_id "[^"]*"' | sort | \
uniq -c | \
awk '$1 == 1' | \
wc -l

# Q16

grep -w 'transcript' Day16.gtf | \
cut -f9 | \
grep -o 'gene_id "[^"]*"' | \
sort | \
uniq -c | \
awk '$1 == 1' | \
wc -l

### Q17, Q18: Number of single exon transcripts
# Count transcripts with one exon:

 # Q17

grep -w 'exon' Day8.gtf | \
cut -f9 | \
grep -o 'transcript_id "[^"]*"' | \
sort | \
uniq -c | \
awk '$1 == 1' | \
wc -l

# Q18

grep -w 'exon' Day16.gtf | \
cut -f9 | \
grep -o 'transcript_id "[^"]*"' | \
sort | \
uniq -c | \
awk '$1 == 1' | \
wc -l

### Q19, Q20: Number of multi-exon transcripts
 # Q19

grep -w 'exon' Day8.gtf | \
cut -f9 | \
grep -o 'transcript_id "[^"]*"' | \
sort | \
uniq -c | \
awk '$1 > 1' | \
wc -l

# Q20

grep -w 'exon' Day16.gtf | \
cut -f9 | \
grep -o 'transcript_id "[^"]*"' | \
sort | \
uniq -c | \
awk '$1 > 1' | \
wc -l

### Q21, Q22: Number of transcripts fully reconstructing annotations

# From gffcompare `.tmap` file, class code `=` means full match:

gffcompare -r athal_genes.gtf -o Day8_cmp Day8.gtf
gffcompare -r athal_genes.gtf -o Day16_cmp Day16.gtf

grep '^=' Day8_cmp.Day8.gtf.tmap | wc -l  # Answer to Q21
grep '^=' Day16_cmp.Day16.gtf.tmap | wc -l  # Answer to Q22

### Q23, Q24: Splice variants of gene AT4G20240
grep 'AT4G20240' Day8.gtf | grep 'transcript_id' | wc -l  # Q23
grep 'AT4G20240' Day16.gtf | grep 'transcript_id' | wc -l # Q24

### Q25, Q26: Number of partial reconstructions (`c` class in gffcompare)
grep '^c' cmp.tmap | wc -l  # Q25 and Q26 for Day8 and Day16

### Q27, Q28: Novel splice variants (`j` class in gffcompare)
grep '^j' cmp.tmap | wc -l  # Q27 and Q28

### Q29, Q30: Transcripts formed in introns (`i` class in gffcompare)
grep '^i' cmp.tmap | wc -l  # Q29 and Q30

### Q31, Q32: Genes and transcripts in merged.gtf

grep 'gene_id "' merged.gtf | sed 's/.*gene_id "\([^"]*\)".*/\1/' | sort | uniq | wc -l  # Q31
grep 'transcript_id' merged.gtf | wc -l  # Q32

### Q33-35: Differential gene and transcript expression

# Load count matrix and metadata
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)

# Q33 total genes tested:
nrow(res)

# Q34 genes differentially expressed (padj < 0.05):
sum(res$padj < 0.05, na.rm=TRUE)

# Q35 transcripts differentially expressed (if separate transcript-level analysis):
# Similar steps on transcript counts

* samtools and grep for read alignment counts and BAM stats.
* StringTie GTF files and awk/grep to parse gene and transcript counts.
* gffcompare to classify transcripts and extract comparison stats.
* DESeq2 or Ballgown in R for differential expression stats.
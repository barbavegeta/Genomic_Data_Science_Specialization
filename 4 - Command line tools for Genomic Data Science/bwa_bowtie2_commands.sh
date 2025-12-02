# 1. Building the Bowtie2 index
# We begin by creating a Bowtie2 index for our reference sequence (IL-2_mRNA.fasta).

bowtie2-build IL-2_mRNA.fasta
# Optionally, we can specify an output directory and index name:

bowtie2-build IL-2_mRNA.fasta hpv/hpv
# This prepares the reference genome so Bowtie2 can efficiently align reads to it.

# 2. Inspecting the read file
# Before mapping, we check how many lines the input FASTQ file contains:

wc -l SRR1107997.fastq
# Since FASTQ files store each read in 4 lines, dividing this number by 4 gives the number of reads.

# 3. Mapping reads with Bowtie2
# We can log the process by redirecting output to a file:

bowtie2 >& bowtie2.log
# For an actual alignment run, we specify the number of threads (-p 4), the index prefix (-x), and the input FASTQ file:

bowtie2 -p 4 -x /. SRR1107997.fastq -S exome.bt2.sam
# The resulting SAM file (exome.bt2.sam) contains the alignment results, which we can inspect:

more exome.bt2.sam | more

# 4. Converting SAM to BAM
# SAM files can be converted to the binary BAM format using samtools:

samtools view -bT exome.bt2.sam > exome.bt2.bam

# 5. Local alignment with Bowtie2
# For more flexible alignments that allow soft clipping, we can run Bowtie2 in local mode:

bowtie2 -p 4 --local -x /. SRR1107997.fastq -S exome.bt2.sam

# 6. Using BWA for alignment
# We can also use the BWA aligner. First, we index the reference:

bwa index IL-2_mRNA.fasta
# Then, log a default bwa mem run:

bwa mem IL-2_mRNA.fasta >& bwa.log
# For an actual mapping run with 4 threads:

bwa mem -t 4 /. SRR1107997.fastq > exome.bwa.sam
# View the SAM:

more exome.bwa.sam
# Convert to BAM:

samtools view -bT /. exome.bwa.sam > exome.bwa.bam
# Inspect the alignment:

samtools view exome.bwa.sam | more

# 7. Pileup generation
# We can generate pileups (base-by-base coverage and variant likelihoods) using samtools mpileup:

samtools mpileup >& mpileup.log
# View a BAM file in SAM format:

samtools view NA12814.bam | more
# Get alignment statistics for a BAM file:

samtools flagstat NA12814.bam
# Generate pileup with a reference genome:

samtools mpileup -f /.Homo.fa sample.bam
# Index a BAM file (required for some downstream tools):

samtools index sample.bam
# Output pileup results to a file:

samtools mpileup -f /.Homo.fasta sample.bam > sample.mpileup
# Inspect pileup file:

more sample.mpileup | more

# 8. Generating VCF or BCF files
# Generate a VCF file directly from BAM:

samtools mpileup -v /.Homo.fasta sample.bam > sample.vcf
more sample.vcf | more
# Generate a BCF binary variant file:

samtools mpileup -g /.Homo.fasta sample.bam > sample.bcf
more sample.bcf | more

# 9. Processing BCF with bcftools
# View BCF contents:

bcftools view sample.bcf
# Log the variant calling process:

bcftools call >& call.log
# Call variants and output as a compressed VCF:

bcftools call -v -m -O z -o sample.vcf.gz sample.bcf
# View compressed VCF:

zcat sample.vcf.gz
# Count variant entries (excluding header lines):

zcat sample.vcf.gz | grep -v "##" | wc -l
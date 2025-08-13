# Question 1
# How many sequences were in the genome?
grep -c "^>" wu_0.v7.fas

# Question 2
# What was the name of the third sequence in the genome file? Give the name only, without the ">" sign.
grep "^>" wu_0.v7.fas | head -3 | tail -1

# Question 3
# What was the name of the last sequence in the genome file? Give the name only, without the ">" sign.
grep "^>" wu_0.v7.fas | tail -1

# Question 4
# How many index files did the operation create?

# Step 1: Suppose the current directory is '/media/sf_gencommand_proj3_data/', which is where supporting data for Exam 3 are stored on the virtual machine for this course. First, we generate the bowtie2 index for the genome file provided. We create a sub-directory 'wu_0' to store the index, then invoke 'bowtie2-build':
mkdir wu_0
bowtie2-build wu_0.v7.fas /media/sf_gencommand_proj3_data/wu_0/wu_0
bowtie2-build >& bowtie2-build.log
# to save and then review the command line usage for the program and decide on the relevant parameters.
# Step 2: To Answer the Question, then simply inspect the content of the index directory and directly observe the number of index files:

ls wu_0/ 

ls wu_0/ 

# Question 5
What is the 3-character extension for the index files created?
bt2
# List the files in the index directory and observe their extension (n.b., all index files have the extension 'bt2'):
ls wu_0/

# Question 6
# How many reads were in the original fastq file?
# Each fastq record is represented on 4 lines in the 'wu_0_A_wgs.fastq' reads file. Simply count the number of lines in the file and divide by 4:
wc -l wu_0_A_wgs.fastq 

# Question 7
# How many matches (alignments) were reported for the original (full-match) setting? Exclude lines in the file containing unmapped reads. 

# Step 1: We first run bowtie2 with two sets of parameters: i) the default parameters, to generate end-to-end read alignments; and ii) the '--local' option, to produce potential partial alignments of a read. For both runs, display the output as SAM alignments (option '-S'):
bowtie2 -x wu_0/wu_0 -U wu_0_A_wgs.fastq -S out.full.sam
bowtie2 -x wu_0/wu_0 -U wu_0_A_wgs.fastq -S out.local.sam --local
# These will create the SAM files 'out.full.sam' and 'out.local.sam'. Upon completion, each run also prints a set of summary statistics on the number of reads unmapped/mapped exactly once/mapped multiple times. Save these to local files, say 'mapped.full.stats' and 'mapped.local.stats'.

# Step 2: To Answer Questions Q7 and Q8, we inspect the SAM files created and determine the number of alignment lines, excluding lines that refer to unmapped reads. A SAM line indicating an unmapped read can be recognized by a '*' in column 3 (chrom). Additionally, we need to exclude the SAM header:
cat out.full.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l
cat out.local.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l

# Question 8
# How many matches (alignments) were reported with the local-match setting? Exclude lines in the file containing unmapped reads. 
# Step 1: We first run bowtie2 with two sets of parameters: i) the default parameters, to generate end-to-end read alignments; and ii) the '--local' option, to produce potential partial alignments of a read. For both runs, display the output as SAM alignments (option '-S'):
bowtie2 -x wu_0/wu_0 -U wu_0_A_wgs.fastq -S out.full.sam
bowtie2 -x wu_0/wu_0 -U wu_0_A_wgs.fastq -S out.local.sam --local
# These will create the SAM files 'out.full.sam' and 'out.local.sam'. Upon completion, each run also prints a set of summary statistics on the number of reads unmapped/mapped exactly once/mapped multiple times. Save these to local files, say 'mapped.full.stats' and 'mapped.local.stats'.
# Step 2: To Answer Questions Q7 and Q8, we inspect the SAM files created and determine the number of alignment lines, excluding lines that refer to unmapped reads. A SAM line indicating an unmapped read can be recognized by a '*' in column 3 (chrom). Additionally, we need to exclude the SAM header:
cat out.full.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l
cat out.local.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l


# Question 9
# How many reads were mapped in the scenario in Question 7?
# This Question explores the difference between reads and alignments. As a reminder, a read may have 0 (unmapped), 1 (unique), or multiple alignments in a SAM file. The information on the number of reads in each category is contained in the 'mapped.full.stats' and 'mapped.local.stats' summary files above: simply add the numbers of reads reported to map exactly once and those reported to match multiple times.
# Note that by default bowtie2 reports only one match per read, so in this case the number of mapped reads at Q9 and the number of alignments at Q7 will be the same.

# Question 10
# How many reads were mapped in the scenario in Question 8?
# This Question explores the difference between reads and alignments. As a reminder, a read may have 0 (unmapped), 1 (unique), or multiple alignments in a SAM file. The information on the number of reads in each category is contained in the 'mapped.full.stats' and 'mapped.local.stats' summary files above: simply add the numbers of reads reported to map exactly once and those reported to match multiple times.
# Note that by default bowtie2 reports only one match per read, so in this case the number of mapped reads at Q10 and the number of alignments at Q8 will be the same.

# Question 11
# How many reads had multiple matches in the scenario in Question 7? You can find this in the bowtie2 summary; note that by default bowtie2 only reports the best match for each read.
# Retrieve this information from the 'mapped.{full,local}.stats' files.

# Question 12
# How many reads had multiple matches in the scenario in Question 8? Use the format above. You can find this in the bowtie2 summary; note that by default bowtie2 only reports the best match for each read.
# Retrieve this information from the 'mapped.{full,local}.stats' files.

# Question 13
# How many alignments contained insertions and/or deletions, in the scenario in Question 7?
# This information is captured in the CIGAR field, marked with 'D' and 'I', respectively:
cut -f6 out.full.sam | grep -c "[I,D]"

# Question 14
# How many alignments contained insertions and/or deletions, in the scenario in Question 8? 
# This information is captured in the CIGAR field, marked with 'D' and 'I', respectively:
cut -f6 out.local.sam | grep -c "[I,D]" 

############# BCFTOOLS VARIANTS

# Question 15
# How many entries were reported for Chr3?
# Step 1: Start by converting the SAM file to BAM format as indicated, then sorting it:
samtools view -bT wu_0.v7.fas out.full.sam > out.full.bam

# then sorting it:
samtools sort out.full.bam out.full.sorted

# This will create the BAM file 'out.full.sorted.bam', which will be used to determine sites of variation.Step 2: Determine candidate sites using 'samtools mpileup', providing the reference fasta genome (option '-f') and using the option '-uv' to report the output in uncompressed VCF format:
###### NEW BCFTOOLS
bcftools mpileup -f wu_0.v7.fas -O b out.full.sorted.bam | bcftools call -mv -Ov > out.full.mpileup.vcf

out.full.mpileup.vcf (generated by bcftools mpileup ... | bcftools call -mv -Ov > out.full.mpileup.vcf)

# Step 3: Count the number of entries in the VCF file located on Chr3. The chromosome information is listed in column 1, once we filter out the header lines (marked with "#"):
cat out.full.mpileup.vcf | grep -v "^#" | cut -f1 | grep -c "^Chr3"

# Potential pitfalls: It is critical to sort the BAM file before analyzing it with samtools.

# Question 16
# How many entries have 'A' as the corresponding genome letter?

# This information is contained in column 4:
cat out.full.mpileup.vcf | grep -v "^#" | cut -f4 | grep -P "^A$"
# where "^A$" tells 'grep' to look for patterns that consist exclusively of 'A' (i.e., between the start '^' and end '$' of the line).

# Question 17
# How many entries have exactly 20 supporting reads (read depth)?

# Read depth is indicated by the 'DP=' field in column 8:
cat out.full.mpileup.vcf | grep -v "^#" | grep -c "DP=20"

# Question 18
# How many entries represent indels?

# This information is marked with the keyword 'INDEL' in the variant line:

cat out.full.mpileup.vcf | grep -v "^#" | grep -c INDEL

# Question 19
# How many entries are reported for position 175672 on Chr1?

# This information is stored in columns 1 and 2 of the VCF file:
cat out.full.mpileup.vcf | grep -v "^#" | cut -f1,2 | grep Chr1 | grep 175672
# then select only the entries corresponding to position 175672.

###### BCFTOOLS VARIANT 
# THIS IS MADE FROM .bcf file
out.final.vcf is about generating the gold-standard, final variant calls through a modular and flexible pipeline that offers more control and refinement.
"""
:
Purpose: This VCF represents the output of the standard, recommended variant calling pipeline. The key advantage is the intermediate .bcf file. This allows for:
Flexibility: You can reuse the out.full.mpileup.bcf (the summarized pileup data) to try different bcftools call parameters, or apply other bcftools sub-commands (like normalization, filtering, or statistics) before generating the final VCF.
Robustness: This two-step process encourages a more controlled and potentially more refined variant calling process, leading to a higher confidence, "final" set of variants after the full calling model has been applied.
"""

# Question 20
How many variants are called on Chr3?

# Step 1: First re-run 'SAMtools mpileup' with the BCF output option '-g':
bcftools mpileup -f wu_0.v7.fas -O b out.full.sorted.bam > out.full.mpileup.bcf

# then call variants using 'BCFtools call' with the multi-allelic caller (option '-m'), showing only variant sites ('-v') and presenting the output in uncompressed VCF format ('-O v'), as instructed:
bcftools call -m -v -O v out.full.mpileup.bcf > out.final.vcf

# Step 2: To Answer the Question, we count all reported variants that show 'Chr3' in column 1:
cat out.final.vcf | grep -v "^#" |  cut -f1 | sort | uniq -c | grep "Chr3"

"""
:
In short: 
The first VCF (out.full.mpileup.vcf) is a direct, all-inclusive variant list from the pileup stage. The second VCF (out.final.vcf) is the result of a standard, multi-step pipeline, often implying a more refined and "final" set of variant calls, particularly when bcftools call applies its full calling model.
"""
# Question 21
# How many variants represent an A->T SNP? If useful, you can use 'grep -P' to allow tabular spaces in the search term.
# This information is stored across columns 2 and 3. An A->T SNP would be represented as an 'A' in column 2 and a 'T' in column 3:

cat out.final.vcf | grep -v "^#" | cut -f4,5 | grep -P "^A\tT$" | wc -l

# Question 22
# How many entries are indels? Answer:
cat out.final.vcf | grep -v "^ #" | grep -c INDEL

# Question 23
# How many entries have precisely 20 supporting reads (read depth)?
cat out.final.vcf | grep -v "^#" | grep -c "DP=20"

# Question 24
# What type of variant (i.e., SNP or INDEL) is called at position 11937923 on Chr3?
# SNP

cat out.final.vcf | grep -v "^#" | cut -f1-5 | grep Chr3 | grep 11937923 
# then inspect the reference and variant sequences in columns 4 and 5 to determine the type of variant, which is a SNP.


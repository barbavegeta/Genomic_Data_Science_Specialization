## Module 2 Exam

#### Q1. How many alignments does the set contain?
samtools view athal_wu_0_A.bam | cut -f3 | grep -v '*' | wc -l
samtools flagstat athal_wu_0_A.bam

#### Q2. How many alignments show the read's mate unmapped?
samtools view athal_wu_0_A.bam | cut -f7 | grep -c '*'

#### Q3. How many alignments contain a deletion (D)?
samtools view athal_wu_0_A.bam | cut -f6 | grep -c 'D'

#### Q4. How many alignments show the read's mate mapped to the same chromosome?
samtools view athal_wu_0_A.bam | cut -f7 | grep -c '='

#### Q5. How many alignments are spliced?
samtools view athal_wu_0_A.bam | cut -f6 | grep -c 'N'

#### Q6. How many alignments does the set contain?
samtools sort -o athal_wu_0_A.sorted.bam athal_wu_0_A.bam
samtools flagstat athal_wu_0_A.region.bam

#### Q7. How many alignments show the read's mate unmapped?
samtools view athal_wu_0_A.region.bam | cut -f7 | grep -c '*'

#### Q8. How many alignments contain a deletion (D)?
samtools view athal_wu_0_A.region.bam | cut -f6 | grep -c 'D'

#### Q9. How many alignments show the read's mate mapped to the same chromosome?
samtools view athal_wu_0_A.bam | cut -f7 | grep -c '='

#### Q10. How many alignments are spliced?
samtools view athal_wu_0_A.bam | cut -f6 | grep -c 'N'

#### Q11. How many sequences are in the genome file?
samtools view -H athal_wu_0_A.bam | grep -c 'SN:'

#### Q12. What is the length of the first file secuence in the genome file?
samtools view -H athal_wu_0_A.bam | grep 'SN:' | more

#### Q13. What alignment tool was used?
samtools view -H athal_wu_0_A.bam | grep '^@PG'

#### Q14. What is the read identifier (name) for the first alignment?
samtools view athal_wu_0_A.bam | head -1 | cut -f1

#### Q15. What is the start position of this read's mate on the genome? Give this as 'chrom:pos' if the read was mapped, or '*' if unmapped.
samtools view athal_wu_0_A.bam | head -n 1 | awk '{print $7":"$8}'

#### Q16. How many overlaps (each overlap is reported on one line) are reported?
bedtools intersect -abam athal_wu_0_A.region.bam -b athal_wu_0_A_annot.gtf -bed -wo > overlaps.bed
wc -l overlaps.bed

#### Q17. How many of these are 10 bases or longer?
cut -f22 overlaps.bed | sort -nrk1 > lengths

#### Q18. How many alignments overlap the annotations?
cut -f1-12 overlaps.bed | sort -u | wc -l

#### Q19. Conversely, how many exons have reads mapped to them?
cut -f13-21 overlaps.bed | sort -u | wc -l

#### Q20. If you were to convert the transcript annotations in the file 'athal_wu_0_A_annot.gtf' into BED format, how many BED records would be generated?
cut -f9 athal_wu_0_A.annot.gtf | cut -d ' ' -f1,2 | sort -u | wc -l
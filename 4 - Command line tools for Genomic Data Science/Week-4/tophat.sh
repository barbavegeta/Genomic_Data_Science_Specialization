mkdir -p /data1/ig=2/florea/Coursera/L4/Tophat/Test1

tophat2 -o /data1/igm2/florea/Coursera/L4/Tophat/Test1 -p 10 -- max-multihits=10 \
-G /data1/igm3/projects/lncRNA/annotation/allmix_nonpseudo.nonredund.gff3 \
-transcriptome-index //data1/igm3/projects/lncRNA/annotation/allmix_nonpseudo_nonredund_bt2index/ \
-r 120 -- mate_std-dev 30 \
/data1/igm3/genomes/hg38/hg38c \
/data1/igm2/florea/Coursera/L4/Test1_1.fastq.gz

# TO THIS

DATADIR=/data1/igm2/florea/Coursera/L4/Data/
WORKDIR=/data1/igm2/florea/Coursera/L4/Tophat/
ANNOT=/data1/igm3/projects/lncRNA/annotation/allmix_nonpseudo.nonredund.gff3
ANNOTIDX=/data1/igm3/projects/lncRNA/annotation/allmix_nonpseudo_nonredund_bt2index
BWT2IDX=/data1/igm3/genomes/hg38/hg38c
mkdir -p $WORKDIR/Test1
tophat2 -o $WORKDIR/Test1 -p 10-max-multihits=10\
-G ANNOT -transcriptome-index $ANNOTIDX \
-r 120-mate_std-dev 30\
$BWT2IDX \
$DATADIR/Test1_1.fastq.gz $DATADIR/Test1_2.fastq.gz

# HOW TO RUN

nohup sh com.tophat >& com.tophat.log

# CUFFLINKS FILE

THDIR=/datal/igm2/florea/Coursera/L4/Tophat/WORKDIR=/data1/igm2/florea/Coursera/Cufflinks
mkdir -p $WORKDIR/Test1
cd $WORKDIR/Test1; cufflinks2 -L BJ-ABZ2 -p 8 $THDIR/Test1/accepted_hits.ban

nohup sh com.cufflinks >& com.cufflinks.log &

#outputs
genes.fpkn_tracking
isoforms.fpkm_tracking
skipped.gtf
transcripts.gtf
more transcripts.gtf

cut -f9 transcripts.gtf | cut -d ' ' -f2 | uniq -u | wc -l

# CUFFDIFF
cuffmerge >& cuffmerge.log

# CUFFDIFF FILE
ANNOT=/data1/igm3/projects/lncRNA/annotation/allmix_nonpseudo.nonredund.gff3
cuffmerge
$ANNOT -p 8-0 $WORKDIR/Cuffmerge $WORKDIR/Cuffmerge/GTFs.txt

# CUFFLINKS
/data1/igm2/florea/Coursera/L4/Cufflinks/Test1/transcripts.gtf
/data1/ign2/florea/Coursera/L4/Cufflinks/Test2/transcripts.gtf
Test3/transcripts. /data1/ign2/florea/Coursera/L4/Cufflinks/Test3/transcripts.gtf
/data1/ign2/florea/Coursera/L4/Cufflinks/Test1/transcripts.gtf oursera/L4/Cufflinks/Test1/transcripts. /data1/igm2/florea/Coursera/L4/Cufflinks/Test1/transcripts.gtf
ra/L4/Cufflin /data1/ign2/florea/Coursera/L4/Cufflinks/Test1/transcripts.gtf


########## AFTER CUFFDIFF FILE

WORKDIR=/data1/igm2/florea/Coursera/L4
THDIR=/data1/igm2/florea/Coursera/L4/Tophat/
ANNOT=/data1/igm3/projects/IncRNA/annotation/allmix_nonpseudo.nonredund.gff3
cuffmerge -g SANNOT -p 80 SWORKDIR/Cuffmerge $WORKDIR/Cuffmerge/GTFs.txt
cuffdiff2 -o SWORKDIR/Cuffdiff -p 10 SWORKDIR/Cuffmerge/merged.gtf \
$THDIR/Test1/accepted_hits.bam, $THDIR/Test2/accepted_hits.bam, $THDIR/Test3/accepted_hits.bam \
$THDIR/Ctrl1/accepted_hits.bam, $THDIR/Ctrl2/accepted_hits.bam, $THDIR/Ctrl3/accepted_hits,bam

nohup ssh com.cuffdiff >& com.cuffdiff.log

#############################

grep yes Cuffdiff/gene_exp.diff | grep chr9


cat Cufflinks/Ctrll/transcripts.gtf | awk '{ if ($1="chr9") print $;}' > Ctrll.gtm
cat Cufflinks/Ctrl2/transcripts.gtf | awk '{ if ($1="chr9") print $1}' > Ctrl2.gtm
cat Cufflinks/Ctrl3/transcripts.gtf | awk '{ if ($1="chr9") print $; } Ctrl3.gtm
cat Cufflinks/Test3/transcripts.gtf | awk '{ if ($1="chr9") print $;}' > Test3.gt
cat Cufflinks/Test2/transcripts.gtf | awk '{ cat Cufflinks/Test1/transcripts.gtf | awk '{ if ($1"chr9") print $;}' > if ($1"chr9") print $;}' > Testi.gtm
Test2.gt

igvtools sort Test1.gtf Test1.sorted.gtf

igvtools sort Ctrl1.gtf Ctrl1.sorted.gtf

igvtools index Ctrl1.sorted.gtf

igvtools index Test1.sorted.gtf
igvtools index Test2.sorted.gtf

samtools view-b Tophat/Test1/accepted_hits.bam "chr9" > chr9.Test1.bam
samtools view chr9.Test1.bam | more

samtols index chr9.Test1.bam # this makes the .bai file

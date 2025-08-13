## Bioconductor for Genomic Data
#20151026
#quiz 3

# Q1  What is the mean expression across all features for sample 5 in the ALL dataset (from the ALL package)?
#load the library
library("ALL")
library("Biobase")
library("genefilter")
#subset the data
sample_5<-ALL[,c(5)]
#get the expression 
sample_5_expression<-exprs(sample_5)
mean(sample_5_expression)
#5.629627

########## ALTERNATIVE ##############
data(ALL)
exprs_data <- exprs(ALL)
mean(exprs_data[, 5])

#############################################################
#############################################################

#Q2:using the Ensembl 75, annotate each feature of the ALL dataset with the Ensembl gene id. how many probesets (features) are annotated with more than one Ensembl gene ID?
library("biomaRt")
library("hgu95av2.db")   # the microarray use in the ALL database

listMarts()
mart<-useMart(host='https://feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
ensembl<-useDataset("hsapiens_gene_ensembl",mart)
class(mart)
listDatasets(mart)
feature_name<-featureNames(ALL)
annotation_ALL<-getBM(attributes=c("ensembl_gene_id","affy_hg_u95av2"),filters="affy_hg_u95av2",values=feature_name,mart=ensembl)
summary(annotation_ALL)
sum(table(annotation_ALL[,2])>1)

########## ALTERNATIVE ##############
mart <- useMart(host='https://feb2014.archive.ensembl.org', biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
all_probes <- featureNames(ALL)

# Get mapping of probes to Ensembl gene IDs for ALL dataset probes
# Suppose all probes are stored in all_probes vector
annot <- getBM(attributes=c("affy_hg_u95av2", "ensembl_gene_id"), filters="affy_hg_u95av2", values = all_probes, mart=mart)

# Count probes with more than one gene
table_probe <- table(annot$affy_hg_u95av2)
sum(table_probe > 1)

#1045

#############################################################
#############################################################

#Q3:How many probesets (Affymetrix IDs) are annotated with one or more genes on the autosomes (chr 1-22)
#參考biomaRt Bioconductor上的說明之task 3
attributes<-listAttributes(ensembl)
list_filter<-listFilters(ensembl)
chrom<-c(1:22)
annotation_ALL_chr<-getBM(attributes=c("ensembl_gene_id","affy_hg_u95av2","chromosome_name"),filters=c("affy_hg_u95av2","chromosome_name"),values=list(feature_name,chrom),mart=ensembl)

table(annotation_ALL_chr[,2])
sum(table(table(annotation_ALL_chr[,2])))
#11016

########## ALTERNATIVE ##############

annot_chr <- getBM(attributes=c("affy_hg_u95av2", "ensembl_gene_id", "chromosome_name"), 
                   filters="affy_hg_u95av2", 
                   values=all_probes, mart=mart)
annot_autosomes <- subset(annot_chr, chromosome_name %in% as.character(1:22))
length(unique(annot_autosomes$affy_hg_u95av2))


#############################################################
#############################################################

# Q4:What is the mean value of the Methylation channel across the features for sample #"5723646052_R04C01" (use the MsetEx dataset from the minfiData package)
#great manual for MsetEx function

#load the library
library(minfiData)
library(minfi)

#explore the data
data(MsetEx)
?MsetEx
head(sampleNames(MsetEx))
head(getMeth(MsetEx))
pData(MsetEx)

#subsetting the data
mean(getMeth(MsetEx)[,2])

#7228.277


########## ALTERNATIVE ##############

data(MsetEx)
meth_values <- getMeth(MsetEx[, "5723646052_R04C01"])
mean(meth_values)

#############################################################
#############################################################

# Q5: Access the processed data from NCBI GEO Accession number GSE788, what is the mean expression level of sample GSM9024?

library(GEOquery)
eList<-getGEO("GSE788")
class(eList)
eData<-eList[[1]]
class(eData)
pData(eData)
names(pData(eData))
mean(exprs(eData)[,2])
#756.432


########## ALTERNATIVE ##############

eList <- getGEO("GSE788")
eData <- eList[[1]]
mean(exprs(eData)[,2])

#############################################################
#############################################################

# Q6:What is the average of the average length across the samples in the expriment?

library(airway)
library(GenomicRanges)

data(airway)
colData(airway)
class(airway)

mean(airway$avgLength)
#113.75



########## ALTERNATIVE ##############
feature_lengths <- width(rowRanges(airway))
mean_length <- mean(feature_lengths)
mean_length

#############################################################
#############################################################

# Q7: What is the number of Ensembl genes which have a count of 1 read or more in sample SRR1039512?
  
table(assay(airway)[,3])
sum(assay(airway,"counts")>=1)

sum(assay(airway)[,3]>=1)
#2599



########## ALTERNATIVE ##############

counts <- assay(airway)
sum(counts[, "SRR1039512"] >= 1)

#############################################################
#############################################################

#Q8:he airway dataset contains more than 64k features. How many of these features overlaps with transcripts on the autosomes (chromosomes 1-22) as represented by the TxDb.Hsapiens.UCSC.hg19.knownGene package?

lapply(c("GenomicFeatures", "TxDb.Hsapiens.UCSC.hg19.knownGene"), library, character.only=TRUE)

#load the packages
# Load packages
library(airway)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

# Load airway data
data(airway)

# Get airway feature ranges
rr <- rowRanges(airway)

# Ensure compatible chromosome names (UCSC style)
seqlevelsStyle(rr) <- "UCSC"

# Get exon ranges grouped by transcript
ex_by_tx <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "tx")

# Collapse each transcript’s exons into one GRanges per transcript (no introns)
collapsed_tx <- reduce(ex_by_tx)

# Flatten the GRangesList into a single GRanges
collapsed_tx <- unlist(collapsed_tx, use.names = FALSE)

# Filter to autosomal chromosomes (chr1–chr22)
collapsed_tx <- keepSeqlevels(collapsed_tx, paste0("chr", 1:22), pruning.mode = "coarse")

# Find overlaps between airway features and exon-only transcript ranges
hits <- findOverlaps(rr, collapsed_tx)

# Count unique overlapping airway features
num_overlap <- length(unique(queryHits(hits)))
cat("Q8 answer:", num_overlap, "\n")



#############################################################
#############################################################

# Question 9 - The expression measures of the airway dataset are the number of reads mapping to each feature. In the previous question we have established that many of these features do not overlap autosomal transcripts from the TxDb.Hsapiens.UCSC.hg19.knownGene. But how many reads map to features which overlaps these transcripts?

overlapping_feature_indices <- unique(queryHits(hits))

# Total reads in this sample
total_reads <- sum(counts[, sample_name])

# Reads in overlapping features only
overlapping_reads <- sum(counts[overlapping_feature_indices, sample_name])

# Fraction
fraction <- overlapping_reads / total_reads
cat("Q9 answer:", fraction, "\n")

# Q9 0.8749633

#############################################################
#############################################################
# Q10 Consider sample SRR1039508 and only consider features which overlap autosomal transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene. We should be able to very roughly divide these transcripts into expressed and non expressed transcript. Expressed transcripts should be marked by H3K4me3 at their promoter. The airway dataset have assayed “airway smooth muscle cells”. In the Roadmap Epigenomics data set, the E096 is supposed to be “lung”. Obtain the H3K4me3 narrowPeaks from the E096 sample using the AnnotationHub package.

library(AnnotationHub)
ah <- AnnotationHub()

query_result <- query(ah, c("E096", "H3K4me3", "narrowPeak"))
narrow_peaks <- query_result[[1]]

# Standardize chromosome style
seqlevelsStyle(narrow_peaks) <- "UCSC"

# Create promoter regions (around transcripts)
tx_promoters <- promoters(transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene), upstream=2200, downstream=200)
tx_promoters <- keepSeqlevels(tx_promoters, paste0("chr", 1:22), pruning.mode="coarse")
seqlevelsStyle(tx_promoters) <- "UCSC"

# Remove chrM from both objects to fix incompatible seqlengths error
rr <- dropSeqlevels(rr, "chrM", pruning.mode = "coarse")
narrow_peaks <- dropSeqlevels(narrow_peaks, "chrM", pruning.mode = "coarse")

# Now run overlaps
peak_hits <- findOverlaps(rr, narrow_peaks)
feature_with_peak <- unique(queryHits(peak_hits))

# Overlap promoters with peaks
peak_hits <- findOverlaps(rr, narrow_peaks)
feature_with_peak <- unique(queryHits(peak_hits))

tx_hits <- findOverlaps(rr, collapsed_tx)
features_autosomal <- unique(queryHits(tx_hits))

# Only keep features overlapping autosomal transcripts
features_with_peak <- intersect(features_autosomal, feature_with_peak)
features_without_peak <- setdiff(features_autosomal, features_with_peak)

# Median counts for SRR1039508
median_with_peak <- median(counts[features_with_peak, sample_name])
median_without_peak <- median(counts[features_without_peak, sample_name])

cat("Q10:\n")
cat("Median counts with H3K4me3:", median_with_peak, "\n")
cat("Median counts without H3K4me3:", median_without_peak, "\n")

# Q10:232


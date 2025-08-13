library(AnnotationHub)
library(GenomicRanges)
library(GenomeInfoDb)

# Question 1
# Use the AnnotationHub package to obtain data on "CpG Islands" in the human genome.
# How many islands exists on the autosomes?

ah <- AnnotationHub()
cpg <- query(ah, c("CpG Islands", "Homo sapiens"))[[1]]

# Use correct seqlevels: chr1 to chr22
autosomal <- keepSeqlevels(cpg, paste0("chr", 1:22), pruning.mode="coarse")

# Now count how many CpG islands are on the autosomes
length(autosomal)

############
#  Question 2
# How many CpG Islands exists on chromosome 4.

cpg_chr4 <- keepSeqlevels(cpg, "chr4", pruning.mode = "coarse")
length(cpg_chr4)

############
#  Question 3

# Obtain the data for the H3K4me3 histone modification for the H1 cell line from Epigenomics Roadmap, using AnnotationHub. Subset these regions to only keep regions mapped to the autosomes (chromosomes 1 to 22).
# How many bases does these regions cover?

query(ah, "H3K4me3")
query(ah, c("H3K4me3", "Homo sapiens"))
query(ah, c("H3K4me3", "H1"))
query(ah, c("H3K4me3", "E003"))  # E003 is Roadmap ID for H1-hESC

# Step 1: Choose the correct dataset

#| ID       | Format        | Description                                                |
#  | -------- | ------------- | ---------------------------------------------------------- |
#  | AH29884  | narrowPeak.gz | Narrow peaks (best for sharp histone marks like H3K4me3) ✅ |
#  | AH28880  | broadPeak.gz  | Broad peaks (better for broader signals like H3K27me3)     |
#  | AH30934  | gappedPeak.gz | Used in some Broad/ENCODE analyses                         |
#  | AH32025+ | bigWig        | Signal tracks (for visualisation, not peak regions)
#|
# This will load the peaks as a GRanges object.
  
h3k4 <- ah[["AH29884"]]

# Step 2: Filter to autosomal chromosomes

h3k4_auto <- keepSeqlevels(h3k4, paste0("chr", 1:22), pruning.mode = "coarse")
sum(width(h3k4_auto))

#################################

# Question 4

# Obtain the data for the H3K27me3 histone modification for the H1 cell line from Epigenomics Roadmap, using the AnnotationHub package. Subset these regions to only keep regions mapped to the autosomes. In the return data, each region has an associated "signalValue". 
# What is the mean signalValue across all regions on the standard chromosomes?

query(ah, "H3K27me3")

query(ah, c("H3K27me3", "E003"))  # E003 is Roadmap ID for H1-hESC

 

h3k27 <- ah[["AH29892"]]
h3k27_auto <- keepSeqlevels(h3k27, paste0("chr", 1:22), pruning.mode="coarse")
mean(mcols(h3k27_auto)$signalValue)

#################################

# Question 5
# Bivalent regions are bound by both H3K4me3 and H3K27me3.
# Using the regions we have obtained above, how many bases on the standard chromosomes are bivalently marked?
  
bivalent <- intersect(h3k4_auto, h3k27_auto)
sum(width(bivalent))

#################################

# Question 6
# We will examine the extent to which bivalent regions overlap CpG Islands.
# how big a fraction (expressed as a number between 0 and 1) of the bivalent regions, overlap one or more CpG Islands?
  
ov <- findOverlaps(bivalent, autosomal)
length(unique(queryHits(ov))) / length(bivalent)

#################################

# Question 7
# How big a fraction (expressed as a number between 0 and 1) of the bases which are part of CpG Islands, are also bivalent marked

# Find overlapping regions between bivalent and autosomal CpG islands
hits <- findOverlaps(bivalent, autosomal)

# Intersect each pair of overlapping regions
overlap_regions <- pintersect(bivalent[queryHits(hits)], autosomal[subjectHits(hits)])

# Proportion of CpG island regions overlapped by bivalent marks
sum(width(overlap_regions)) / sum(width(autosomal))


#################################

# Question 8
# How many bases are bivalently marked within 10kb of CpG Islands?
# Tip: consider using the "resize()"" function.

# Expand CpG island regions by 10kb in both directions (total 20kb)
cpg_expanded <- resize(autosomal, width = 20000, fix = "center")

# Find bivalent marks that fall within these expanded CpG regions
nearby_bivalents <- subsetByOverlaps(bivalent, cpg_expanded)

# Total number of base pairs in these nearby bivalent regions
sum(width(nearby_bivalents))


#################################

# Question 9

# How big a fraction (expressed as a number between 0 and 1) of the human genome is contained in a CpG Island?
# Tip 1: the object returned by AnnotationHub contains "seqlengths".
# Tip 2: you may encounter an integer overflow. As described in the session on R Basic Types, you can address this by converting integers to numeric before summing them, "as.numeric()".

genome_size <- sum(as.numeric(seqlengths(cpg)[paste0("chr", 1:22)]))
cpg_size <- sum(width(autosomal))
tot<-cpg_size / genome_size


#################################

# Question 10

# Compute an odds-ratio for the overlap of bivalent marks with CpG islands.

# Find overlaps
hits <- findOverlaps(bivalent, autosomal)
overlap <- pintersect(bivalent[queryHits(hits)], autosomal[subjectHits(hits)])

# Counts
a <- sum(width(overlap))
b <- sum(width(autosomal)) - a
c <- sum(width(bivalent)) - a
d <- 2.7e9 - (a + b + c)

# Odds ratio
odds_ratio <- (a / b) / (c / d)
odds_ratio


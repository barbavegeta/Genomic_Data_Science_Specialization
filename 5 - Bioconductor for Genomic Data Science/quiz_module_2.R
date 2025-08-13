BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(AnnotationHub)
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Biostrings)

ah <- AnnotationHub()

### Q1: What is the GC content of “chr22” in the “hg19” build of the human genome?
# GC content of chr22 (excluding Ns)

library(BSgenome.Hsapiens.UCSC.hg19)
available.genomes()
seqlengths(BSgenome.Hsapiens.UCSC.hg19)
seqnames(BSgenome.Hsapiens.UCSC.hg19)
chr22 <- BSgenome.Hsapiens.UCSC.hg19$chr22

# Get all base counts
freqs <- alphabetFrequency(chr22, baseOnly = FALSE)

# Extract counts
gc <- freqs["G"] + freqs["C"]
atgc_total <- freqs["A"] + freqs["T"] + freqs["G"] + freqs["C"]
gc_content_chr22 <- gc / atgc_total
gc_content_chr22


################### ALTERNATIVE #########################

# Calculate letter frequencies (probabilities)
A <- letterFrequency(chr22_seq, "A", as.prob = TRUE)
C <- letterFrequency(chr22_seq, "C", as.prob = TRUE)
T <- letterFrequency(chr22_seq, "T", as.prob = TRUE)
G <- letterFrequency(chr22_seq, "G", as.prob = TRUE)
N <- letterFrequency(chr22_seq, "N", as.prob = TRUE)

# Calculate raw counts
nA <- letterFrequency(chr22_seq, "A", as.prob = FALSE)
nC <- letterFrequency(chr22_seq, "C", as.prob = FALSE)
nT <- letterFrequency(chr22_seq, "T", as.prob = FALSE)
nG <- letterFrequency(chr22_seq, "G", as.prob = FALSE)
nN <- letterFrequency(chr22_seq, "N", as.prob = FALSE)

# Calculate frequencies excluding Ns
total_bases <- nA + nC + nT + nG
fA <- nA / total_bases
fC <- nC / total_bases
fT <- nT / total_bases
fG <- nG / total_bases

# Calculate ratio of C*G over total bases
cg_ratio <- (nC * nG) / total_bases

# GC content and adjusted values
gc_content <- (C + G) / (A + C + T + G) 

#######################################################
#######################################################


### Q2: What is the mean GC content of H3K27me3 “narrowPeak” regions of Epigenomics Roadmap from the H1 stem cell line on chr 22.
# Mean GC content in H3K27me3 narrowPeak regions on chr22 (E003)
ah <- AnnotationHub()

query(ah, c("H3K27me3", "E003"))  # E003 is Roadmap ID for H1-hESC
np <- ah[["AH29892"]]

# Load narrowPeak
np_chr22 <- keepSeqlevels(np, "chr22", pruning.mode="coarse")

# Get sequences
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, np_chr22)
gc_percents <- letterFrequency(seqs, letters=c("G","C"), as.prob=TRUE)
mean(rowSums(gc_percents))

################### ALTERNATIVE #########################

# Keep only chr22 data and trim
ah_chr22 <- keepSeqlevels(ah_1, "chr22")
ah_chr22 <- trim(ah_chr22)

# Extract Views on the genome
v1 <- Views(BSgenome.Hsapiens.UCSC.hg19$chr22, ranges(ah_chr22))
mean_gc <- mean(letterFrequency(v1, "GC", as.prob = TRUE))

#######################################################
#######################################################

### Q3: What is the correlation between GC content and “signalValue” of these regions (on chr22)?
# Correlation between GC content and signalValue

gc_vals <- rowSums(letterFrequency(seqs, letters=c("G","C"), as.prob=TRUE))
cor(gc_vals, mcols(np_chr22)$signalValue)

################### ALTERNATIVE #########################

df <- data.frame(GC_content = letterFrequency(v1, "GC", as.prob = TRUE),
                 SignalValue = ah_chr22$signalValue)
correlation <- cor(df$GC_content, df$SignalValue)

#######################################################
#######################################################

### Q4:  What is the correlation between the “signalValue” of the “narrowPeak” regions and the average “fc.signal” across the same regions?
# Correlation of signalValue with average fc.signal


ah_chip <- query(ah, c("fc.signal", "H3K27me3", "E003"))
ah_chip <- ah_chip[["AH32033"]]

# Import signal data for chr22 as Rle
chr22_range <- GRanges("chr22", IRanges(1, 51304566))
ah_chip_bw <- import(ah_chip, which = chr22_range, as = "Rle")

# Create views corresponding to peak ranges
ah_chip_bw_view <- Views(ah_chip_bw$chr22, start = start(ah_chr22), end = end(ah_chr22))

# Calculate means and correlations
mean_signal <- mean(ah_chip_bw_view)
cor_signal <- cor(mean_signal, ah_chr22$signalValue)

#######################################################
#######################################################

### Q5: How many bases on chr22 have an fc.signal greater than or equal to 1?
# Bases on chr22 with fc.signal ≥ 1


count_bases_fc_ge_1 <- sum(ah_chip_bw$chr22 >= 1)


#######################################################
#######################################################

### Q6: Identify the regions of the genome where the signal in E003 is 0.5 or lower and the signal in E055 is 2 or higher.
# Regions where E003 ≤ 0.5 and E055 ≥ 2

ah_E055 <- query(ah, c("H3K27me3", "narrowPeak", "E055"))[["AH30313"]]
ah_E055_fc <- query(ah, c("E055", "fc.signal", "H3K27me3"))[["AH32470"]]

# Import E055 bigwig as Rle for chr22
ah_E055_fc_import <- import(ah_E055_fc, which = chr22_range, as = "Rle")

# Slice by signal thresholds
E055_ir <- slice(ah_E055_fc_import$chr22, lower = 2)
E003_ir <- slice(ah_chip_bw$chr22, upper = 0.5)

# Convert to IRanges and intersect
E055_ir_ranges <- as(E055_ir, "IRanges")
E003_ir_ranges <- as(E003_ir, "IRanges")
comb_ir <- intersect(E055_ir_ranges, E003_ir_ranges)
combined_width <- sum(width(comb_ir))  # ~1869937 bases

combined_width 
#######################################################
#######################################################

### Q7: What is the average observed-to-expected ratio of CpG dinucleotides for CpG Islands on chromosome 22?

# Query CpG islands from AnnotationHub
CpG_islands <- query(ah, c("hg19", "CpG"))
CpG_islands <- subset(CpG_islands, species == "Homo sapiens")
CpG_islands <- CpG_islands[["AH5086"]]

# Extract chr22 CpG islands
CpG_islands_split <- split(CpG_islands, seqnames(CpG_islands))
CpG_islands_chr22 <- CpG_islands_split$chr22

# Create Views of CpG islands on chr22
CpG_island_views <- Views(BSgenome.Hsapiens.UCSC.hg19$chr22,
                          start = start(CpG_islands_chr22),
                          end = end(CpG_islands_chr22))

# Dinucleotide frequency for CpG
CpG_dinuc_freq <- dinucleotideFrequency(CpG_island_views)[, "CG"]

# Frequencies of G and C in CpG islands
CpG_G <- letterFrequency(CpG_island_views, "G")
CpG_C <- letterFrequency(CpG_island_views, "C")

# Calculate observed-to-expected CpG ratio for CpG islands
obs_exp_CpG <- mean((CpG_dinuc_freq * width(CpG_islands_chr22)) / (CpG_G * CpG_C))  # ~0.8341
obs_exp_CpG
#######################################################
#######################################################

### Q8: How many TATA boxes are there on chr 22 of build hg19 of the human genome?
# Count all "TATAAA" motifs on both strands

seq <- BSgenome.Hsapiens.UCSC.hg19$chr22
forward_hits <- countPattern("TATAAA", seq)
reverse_hits <- countPattern("TATAAA", reverseComplement(seq))
total_hits <- forward_hits + reverse_hits
total_hits

################### ALTERNATIVE #########################

TATA <- DNAString("TATAAA")
chr22_seq <- BSgenome.Hsapiens.UCSC.hg19$chr22

# Match pattern on both strands
match_TATA <- matchPattern(TATA, chr22_seq)
match_TATA_rc <- matchPattern(reverseComplement(TATA), chr22_seq)
total_TATA <- length(match_TATA) + length(match_TATA_rc)  # 27263

# Alternative: matchPattern for entire genome or chr22 only
match_TATA <- matchPattern(TATA, chr22_seq)
length(match_TATA)  # Should match total_TATA

#######################################################
#######################################################

### Q9: How many transcript promoters are on chromosome 22, which contain a coding sequence, such as a TATA box on the same strand as the transcript?
# Promoters of coding transcripts with TATA box

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Restrict to chr22
seqlevels(txdb, force = TRUE) <- "chr22"
gr_chr22 <- GRanges("chr22", IRanges(1, 52330658))

# Extract transcripts on chr22
transcripts_chr22 <- subsetByOverlaps(transcripts(txdb), gr_chr22, ignore.strand = TRUE)

# Define promoters (upstream 900bp, downstream 100bp)
promoters_chr22 <- promoters(transcripts_chr22, upstream = 900, downstream = 100)

# Extract coding sequences (CDS)
cds_chr22 <- cds(txdb)

# Find promoters overlapping CDS (coding transcripts)
promoters_with_cds <- promoters_chr22[queryHits(findOverlaps(promoters_chr22, cds_chr22))]

# Extract promoter sequences
promoter_seqs <- Views(chr22_seq, start = start(promoters_with_cds), end = end(promoters_with_cds))

# Count TATA boxes in promoters (same strand matching requires further strand info)
count_TATA_in_promoters <- sum(countPattern(TATA, promoter_seqs)) +
  sum(countPattern(reverseComplement(TATA), promoter_seqs))

count_TATA_in_promoters

#######################################################
#######################################################

### Q10: How many bases on chr22 are part of more than one promoter of a coding sequence?
# Bases in more than one promoter (ignore strand)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# get transcripts
tx <- transcripts(txdb, use.names=TRUE)

# get promoters (900bp upstream, 100bp downstream)
proms_all <- promoters(tx, upstream=900, downstream=100)

# reduce promoters to merge overlaps
reduce_proms <- reduce(proms_all)

# get coverage
coverage_all <- coverage(proms_all)

# for example: sum run lengths where coverage is > 1 for chr22
sum(runLength(coverage_all$chr22)[runValue(coverage_all$chr22) > 1])

################### ALTERNATIVE #########################


## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(GenomicFeatures)

## ----biocLite, eval=FALSE------------------------------------------------
## BiocManager::install("GenomicFeatures")
## BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
## ----txdb----------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

## ----gr------------------------------------------------------------------
gr <- GRanges(seqnames = "chr1", strand = "+", ranges = IRanges(start = 11874, end = 14409))
genes(txdb)
subsetByOverlaps(genes(txdb), gr)
subsetByOverlaps(genes(txdb), gr, ignore.strand = TRUE)

## ----transcripts---------------------------------------------------------
subsetByOverlaps(transcripts(txdb), gr)

## ----exons---------------------------------------------------------------
subsetByOverlaps(exons(txdb), gr)

## ----exonsBy-------------------------------------------------------------
subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)

## ----cds-----------------------------------------------------------------
subsetByOverlaps(cds(txdb), gr)
subsetByOverlaps(cdsBy(txdb, by = "tx"), gr)

## ----transcriptLengths---------------------------------------------------
subset(transcriptLengths(txdb, with.cds_len = TRUE), gene_id == "100287102")

sum(width(subsetByOverlaps(cdsBy(txdb, by "tx"), gr)[["2"]]))
## ----back, child="back.Rmd", echo=FALSE----------------------------------

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()



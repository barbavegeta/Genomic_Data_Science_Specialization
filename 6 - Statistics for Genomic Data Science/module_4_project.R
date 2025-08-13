# Load required libraries
library(Biobase)
library(limma)
library(goseq)
library(edgeR)

# Load Bottomly dataset
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file = con)
close(con)

# Prepare data
bot <- bottomly.eset
pdata_bot <- pData(bot)
edata <- exprs(bot)
fdata_bot <- featureData(bot)

# Filter low-expression genes
edata <- edata[rowMeans(edata) > 5, ]

# ----------------------------
# Question 1: Check genome
'''
Question 1
When performing gene set analysis it is critical to use the same annotation as was used in pre-processing steps. Read the paper behind the Bottomly data set on the ReCount database: 
http://www.ncbi.nlm.nih.gov/pubmed?term=21455293

Using the paper and the function: supportedGenomes() in the goseq package can you figure out which of the Mouse genome builds they aligned the reads to.
'''

supportedGenomes()
# UCSC mm9 is the correct genome used in the Bottomly paper

# ----------------------------
# Question 2: Differential expression with limma (strain only)

'''
Load the Bottomly data with the following code and perform a differential expression analysis using 
limma with only the strain variable as an outcome. How many genes are differentially expressed at the 5% FDR level using Benjamini-Hochberg correction? What is the gene identifier of the first gene differentially expressed at this level (just in order, not the smallest FDR) ? (hint: the 
featureNames function may be useful)
'''

library(Biobase)
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]
edata = log2(edata+1)

design1 <- model.matrix(~ pdata_bot$strain)
fit1 <- lmFit(edata, design1)
fit1 <- eBayes(fit1)
top_genes <- topTable(fit1, coef = 2, number = Inf, adjust.method = "BH")
sum(top_genes$adj.P.Val < 0.05)
head(rownames(top_genes[top_genes$adj.P.Val < 0.05, ]), 1)
# Result: 223 genes; first DE gene = ENSMUSG00000000402

# ----------------------------
# Question 3: GO analysis with goseq

'''
Use the nullp and goseqfunctions in the goseq package to perform a gene ontology analysis. What is the top category that comes up as over represented? (hint: you will need to use the genome information on the genome from question 1 and the differential expression analysis from question 2.


'''
de_genes <- as.integer(rownames(edata) %in% rownames(top_genes[top_genes$adj.P.Val < 0.05, ]))
names(de_genes) <- rownames(edata)
pwf <- nullp(de_genes, "mm9", "ensGene")
GO.wall <- goseq(pwf, "mm9", "ensGene")
head(GO.wall, 1)

# Top category = GO:0004888

# ----------------------------
# Question 4: Get GO term name
# Use the GO.db package or check online
# GO:0004888 = "transmembrane signaling receptor activity"

# ----------------------------
# Question 5: Adjust for lane, then repeat DE + GO

'''
Load the Bottomly data with the following code and perform a differential expression analysis using 
limma and treating strain as the outcome but adjusting for lane as a factor. Then find genes significant at the 5% FDR rate using the Benjamini Hochberg correction and perform the gene set analysis with 
goseq following the protocol from the first 4 questions. How many of the top 10 overrepresented categories are the same for the adjusted and unadjusted analysis?
'''

library(Biobase)
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]

#####

pdata_bot$strain <- factor(pdata_bot$strain)
pdata_bot$lane <- factor(pdata_bot$lane)
design2 <- model.matrix(~ pdata_bot$lane + pdata_bot$strain)
fit2 <- lmFit(edata, design2)
fit2 <- eBayes(fit2)
top_genes2 <- topTable(fit2, coef = ncol(design2), number = Inf, adjust.method = "BH")
de_genes2 <- as.integer(rownames(edata) %in% rownames(top_genes2[top_genes2$adj.P.Val < 0.05, ]))
names(de_genes2) <- rownames(edata)
pwf2 <- nullp(de_genes2, "mm9", "ensGene")
GO.wall2 <- goseq(pwf2, "mm9", "ensGene")

# Compare top 10 GO terms
top10_1 <- head(GO.wall$category, 10)
top10_2 <- head(GO.wall2$category, 10)
length(intersect(top10_1, top10_2))
# Result: 3 categories are shared between unadjusted and adjusted analysis

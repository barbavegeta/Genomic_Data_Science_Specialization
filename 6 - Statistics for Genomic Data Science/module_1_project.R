# Load required libraries
library(Biobase)
library(GenomicRanges)
library(plotrix)


###

######################################################################
######################################################################

# Question 2
# Put the following code chunk at the top of an R markdown document called test.Rmd but set eval=TRUE


```{r setup, eval=FALSE}
knitr::opts_chunk$set(cache=TRUE)
```
Then create the following code chunks


```{r }
x = rnorm(10)
plot(x,pch=19,col="dodgerblue")
```

```{r }
y = rbinom(20,size=1,prob=0.5)
table(y)
```

# Answer
'''
The plot is random the first time you knit the document. It is identical to the first time the second time you knit the document. After removing the folders 
test_cache
test_cachestart verbatim, test_cache, end verbatim and 
test_files
test_filesstart verbatim, test_files, end verbatim they generate new random versions.
'''

########

######################################################################
######################################################################

# Question 3: create a summarizedExperiment object with the following code

library(Biobase)
library(GenomicRanges)
data(sample.ExpressionSet, package = "Biobase")
se = makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet)

'''
Look up the help files for summarizedExperiment with the code?summarizedExperiment. How do you access the genomic data for this object? How do you access the phenotype table? How do you access the feature data? What is the unique additional information provided by rowRanges (se)?

# Get the genomic table with assay (se), get the phenotype table with colData (se), get the feature data with rowData(se). rowRanges (se) gives information on the genomic location and structure of the measured features.

# Get the genomic table with assay (se), get the phenotype table with colData (se), get the feature data with rowRanges (se). rowRanges (se) gives the range of possible values for the expression data.

# Get the genomic table with <code>assay(se)</code>, get the phenotype table with <code>pData(se)</code>, get the feature data with <code>rowData(se)</code>, <code>rowRanges (se)</code> gives information on the genomic location and structure of the measured features.

# Get the genomic table with assay (se), get the phenotype table with colData (se), get the feature data with rowData(se). rowRanges (se) gives the range of possible values for the expression data.
'''
# Access components
genomic_table <- assay(se)
phenotype_table <- colData(se)
feature_data <- rowData(se)
genomic_ranges <- rowRanges(se)

######################################################################
######################################################################

# Question 5: Load the Bottomly and the Bodymap data sets with the following code:
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file = con)
close(con)
bot <- bottomly.eset
pdata_bot <- pData(bot)

con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file = con)
close(con)
bm <- bodymap.eset
pdata_bm <- pData(bm)

# Just considering the phenotype data what are some reasons that the Bottomly data set is likely a better experimental design than the Bodymap data? Imagine the question of interest in the Bottomly data is to compare strains and in the Bodymap data it is to compare tissues.

######################################################################
######################################################################

# Question 6: What are some reasons why this plot is not useful for comparing the number of technical replicates by tissue (you may need to install the plotrix package).

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata_bm=pData(bm)

library(plotrix)
pie3D(pdata_bm$num.tech.reps,labels=pdata_bm$tissue.type)
# Question 7: Heatmap of 500 most highly expressed genes
edata <- exprs(bm)
row_sums <- rowSums(edata)
edata <- edata[order(-row_sums), ]
index <- 1:500
heatmap(edata[index, ], Rowv = NA, Colv = NA)

######################################################################
######################################################################

# Question 7: Load the Bottomly data

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)

'''
Which of the following code chunks will make a heatmap of the 500 most highly expressed genes (as defined by total count), without re-ordering due to clustering? Are the highly expressed samples next to each other in sample order?
'''

row_sums = rowSums(edata)
edata = edata[order(-row_sums),]
index = 1:500
heatmap(edata[index,],Rowv=NA,Colv=NA)

######################################################################
######################################################################

# Question 8: MA plot using log2 and rlog
library(DESeq2)
log2_MA_plot <- function(sample1, sample2) {
  log_counts1 <- log2(sample1 + 1)
  log_counts2 <- log2(sample2 + 1)
  M <- log_counts1 - log_counts2
  A <- (log_counts1 + log_counts2) / 2
  plot(A, M, main = "MA Plot: log2", xlab = "A", ylab = "M", pch = 20, col = "grey")
}

dds <- DESeqDataSetFromMatrix(countData = edata, colData = pdata_bm, design = ~1)
rld <- rlog(dds, blind = TRUE)

rlog_MA_plot <- function(rld_mat) {
  M <- assay(rld_mat)[,1] - assay(rld_mat)[,2]
  A <- (assay(rld_mat)[,1] + assay(rld_mat)[,2]) / 2
  plot(A, M, main = "MA Plot: rlog", xlab = "A", ylab = "M", pch = 20, col = "blue")
}

log2_MA_plot(edata[,1], edata[,2])
rlog_MA_plot(rld)

######################################################################
######################################################################

# Question 9: Load the Montgomery and Pickrell eSet:
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file = con)
close(con)
mp <- montpick.eset
pdata <- pData(mp)
edata <- as.data.frame(exprs(mp))
fdata <- fData(mp)

'''
Cluster the data in three ways:

With no changes to the data

After filtering all genes with rowMeanaless than 100

After taking the log2 transform of the data without filtering

Color the samples by which study they came from (Hint: consider using the function 
myplclust.R in the package rafalib available from CRAN and looking at the argument 
lab.col .)
'''
# Clustering variants
hc_original <- hclust(dist(t(edata)))
hc_filtered <- hclust(dist(t(edata[rowMeans(edata) >= 100, ])))
hc_log2 <- hclust(dist(t(log2(edata + 1))))

# Optional: color by study (lab.col argument requires rafalib)
# library(rafalib)
# myplclust(hc_log2, lab.col = as.numeric(as.factor(pdata$study)))

######################################################################
######################################################################

# Question 10: Load the Montgomery and Pickrell eSet:
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

'''
Cluster the samples using k-means clustering after applying the 
log2transform (be sure to add 1). Set a seed for reproducible results (use 
set.seed(1235). If you choose two clusters, do you get the same two clusters as you get if you use the 
cutree function to cluster the samples into two groups? Which cluster matches most closely to the study labels?

'''
# Compare with cutree
hc <- hclust(dist(t(edata_log)))
cut_clusters <- cutree(hc, k = 2)

# Table comparing kmeans and cutree results
table(kmeans = km$cluster, cutree = cut_clusters, study = pdata$study)

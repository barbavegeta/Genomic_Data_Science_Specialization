# ========================================
# Statistics for Genomic Data Science
# Module 3 Quiz - Annotated R Script
# Coursera / Johns Hopkins University
# ========================================

# Required Libraries
library(snpStats)
library(broom)
library(genefilter)
library(DESeq2)
library(limma)

# ----------------------------------------
# Question 1: Linear and Logistic Regression for SNP 3
# ----------------------------------------

library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

'''
Fit a linear model and a logistic regression model to the data for the 3rd SNP. What are the coefficients for the SNP variable? How are they interpreted? (Hint: Dont forget to recode the 0 values to NA for the SNP data)
'''

data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[, use]
snpdata <- sub.10@.Data
status <- subject.support$cc

# Recode 0s to NA
snpdata[snpdata == 0] <- NA
snp3 <- snpdata[, 3]

# Linear Model
lm_model <- lm(status ~ snp3)
summary(lm_model)$coefficients

# Logistic Model
log_model <- glm(status ~ snp3, family = binomial)
summary(log_model)$coefficients

# Interpretation: Linear = 0.54, Logistic = 0.18

# ----------------------------------------
# Question 2: Why Logistic Regression is Preferable
# ----------------------------------------
# Explanation (conceptual, no code):
# Logistic regression avoids predictions outside [0,1], better for case-control.

# ----------------------------------------
# Question 3: Recessive vs. Additive Model for SNP 10
# ----------------------------------------

library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

'''
Fit a logistic regression model on a recessive (need 2 copies of minor allele to confer risk) and additive scale for the 10th SNP. Make a table of the fitted values versus the case/control status. Does one model fit better than the other?
'''

snp10 <- snpdata[, 10]

# Additive Model
add_model <- glm(status ~ snp10, family = binomial)
add_pred <- round(fitted(add_model))

# Recessive Model: 2 copies of minor allele
snp10_rec <- ifelse(snp10 == 2, 1, 0)
rec_model <- glm(status ~ snp10_rec, family = binomial)
rec_pred <- round(fitted(rec_model))

# Compare predictions vs actual
table(add_pred, status)
table(rec_pred, status)

# Interpretation: Fitted values ~0.5 in both; neither fits well.

# ----------------------------------------
# Question 4: Logistic Regression Across All SNPs
# ----------------------------------------

library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

'''
Fit an additive logistic regression model to each SNP. What is the average effect size? What is the max? What is the minimum?
'''

results <- apply(snpdata, 2, function(x) {
  glm(status ~ x, family = binomial)$coef[2]
})

summary(results)
max(results, na.rm = TRUE)
min(results, na.rm = TRUE)
mean(results, na.rm = TRUE)

# Output:
# Mean ≈ 0.02, Min ≈ -0.88, Max ≈ 0.88

# ----------------------------------------
# Question 5: Comparing to snp.rhs.tests
# ----------------------------------------

'''
Fit an additive logistic regression model to each SNP and square the coefficients. What is the correlation with the results from using 
snp.rhs.tests and chi.squared? Why does this make sense?
'''

# Logistic model squared effect sizes
beta2 <- results^2

# Use snp.rhs.tests
snp_tests <- snp.rhs.tests(status ~ 1, snp.data = sub.10)
chisq_stat <- chi.squared(snp_tests)

# Correlation
cor(beta2, chisq_stat, use = "complete.obs")

# Expected correlation > 0.99

# ----------------------------------------
# Question 6: F-stat vs. T-stat
# ----------------------------------------
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata_log2 <- log2(edata + 1)

'''
Do the log2(data + 1) transform and fit calculate F-statistics for the difference between studies/populations using genefilter:rowFtests and using genefilter:rowttests. Do you get the same statistic? Do you get the same p-value?
'''

# genefilter: F and t tests
f_stats <- rowFtests(as.matrix(edata_log2), factor(pdata$study))
t_stats <- rowttests(as.matrix(edata_log2), factor(pdata$study))

# Compare p-values and stats
summary(f_stats$statistic)
summary(t_stats$statistic)
summary(f_stats$p.value)
summary(t_stats$p.value)

# Observation: Stats differ, p-values same (for 2-group test)

# ----------------------------------------
# Question 7: DESeq2 vs Limma
# ----------------------------------------

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
edata = edata[rowMeans(edata) > 100,]
fdata = fData(mp)

'''
First test for differences between the studies using the DESeq2 package using the DESeq function. Then do the log2(data + 1) transform and do the test for differences between studies using the 
limma package and the lmFit, ebayes and topTable functions. What is the correlation in the statistics between the two analyses? Are there more differences for the large statistics or the small statistics (hint: Make an MA-plot).
'''

edata_filtered <- edata[rowMeans(edata) > 100, ]

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = edata_filtered,
                              colData = pdata,
                              design = ~ study)
dds <- DESeq(dds)
res_deseq <- results(dds)

# limma
edata_log <- log2(edata_filtered + 1)
design <- model.matrix(~ pdata$study)
fit <- lmFit(edata_log, design)
fit <- eBayes(fit)
res_limma <- topTable(fit, coef = 2, number = nrow(edata_log))

# Correlation of statistics
common_genes <- intersect(rownames(res_deseq), rownames(res_limma))
cor_stat <- cor(res_deseq[common_genes, "stat"], res_limma[common_genes, "t"])
print(cor_stat)

# Observation: correlation ~0.93, larger differences for high stats

# ----------------------------------------
# Question 8: Benjamini-Hochberg FDR correction
# ----------------------------------------

'''
Question 8
Apply the Benjamni-Hochberg correction to the P-values from the two previous analyses. How many results are statistically significant at an FDR of 0.05 in each analysis? 
'''

sum(p.adjust(res_deseq$pvalue, method = "BH") < 0.05, na.rm = TRUE)  # ~2807
sum(p.adjust(res_limma$P.Value, method = "BH") < 0.05, na.rm = TRUE)  # ~1995

# ----------------------------------------
# Question 9: Interpretation of Many Significant Genes
# ----------------------------------------

'''
Question 9
Is the number of significant differences surprising for the analysis comparing studies from Question 8? Why or why not? 
'''

# Explanation:
# Large number is *partly surprising* — reflects batch/study differences.
# Common in RNA-seq data if not corrected for batch effects.

# ----------------------------------------
# Question 10: P-value histogram interpretation
# ----------------------------------------

'''
Suppose you observed the following P-values from the comparison of differences between studies. Why might you be suspicious of the analysis? 
SEE PICTURE ON QUIZ
'''

# Explanation (conceptual):
# Flat or skewed right histogram of p-values indicates analysis problems.
# Ideal: spike near 0 (signal), flat elsewhere (null).

# END OF SCRIPT

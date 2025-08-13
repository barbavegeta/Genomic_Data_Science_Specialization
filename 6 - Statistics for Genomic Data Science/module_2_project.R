# ---------------------------------------------
# Module 2 Quiz - Statistics for Genomic Data Science (Coursera)
# Johns Hopkins University
# Complete Code for Questions 1–10
# ---------------------------------------------

# Load required libraries
library(Biobase)
library(limma)
library(sva)

# ---------------------------------------------
# Question 1: PCA Variance Explained
# ---------------------------------------------
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
'''
What percentage of variation is explained by the 1st principal component in the data set if you:

Do no transformations?

log2(data + 1) transform?

log2(data + 1) transform and subtract row means?
'''
# a. No transformation
pc1_var1 <- summary(prcomp(t(edata)))$importance[2,1]

# b. log2(data + 1)
edata_log <- log2(edata + 1)
pc1_var2 <- summary(prcomp(t(edata_log)))$importance[2,1]

# c. log2(data + 1) and subtract row means
edata_log_centered <- edata_log - rowMeans(edata_log)
pc1_var3 <- summary(prcomp(t(edata_log_centered)))$importance[2,1]

cat("Q1 - PC1 variance:\n")
cat("No transform:", pc1_var1, "\nlog2(x+1):", pc1_var2, "\nlog2(x+1) - rowMeans:", pc1_var3, "\n\n")

# ---------------------------------------------
# Question 2: Correlation between SVD and K-means clusters
# ---------------------------------------------
set.seed(333)
km <- kmeans(t(edata_log_centered), centers=2)
svd1 <- svd(edata_log_centered)
correlation_q2 <- cor(svd1$v[,1], km$cluster)
cat("Q2 - Correlation between SVD and K-means:", correlation_q2, "\n\n")

'''
Perform the log2(data + 1) transform and subtract row means from the samples. Set the seed to 
333 and use k-means to cluster the samples into two clusters. Use 
svd to calculate the singular vectors. What is the correlation between the first singular vector and the sample clustering indicator?
'''

# ---------------------------------------------
# Question 3: Gene vs Technical Replicates
# ---------------------------------------------
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

'''
Fit a linear model relating the first gene’s counts to the number of technical replicates, treating the number of replicates as a factor. Plot the data for this gene versus the covariate. Can you think of why this model might not fit well?

'''

gene1 <- edata_bm[1,]
num.reps <- as.factor(pdata_bm$num.tech.reps)
model_q3 <- lm(gene1 ~ num.reps)
plot(num.reps, gene1, main="Q3 - Gene 1 vs Technical Replicates", xlab="Num Tech Reps", ylab="Expression")
abline(model_q3)
summary(model_q3)

# ---------------------------------------------
# Question 4: Linear model with age and sex
# ---------------------------------------------

'''
Fit a linear model relating he first gene’s counts to the age of the person and the sex of the samples. What is the value and interpretation of the coefficient for age?

'''

model_q4 <- lm(gene1 ~ pdata_bm$age + pdata_bm$sex)
summary(model_q4)

# ---------------------------------------------
# Question 5: Regression model to predict population
# ---------------------------------------------
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
X <- model.matrix(~ pdata$population)
fit_q5 <- lm.fit(X, t(edata_log))

cat("Q5 - Dimensions:\n")
cat("Residuals:", dim(fit_q5$residuals), "\n")
cat("Effects:", dim(fit_q5$effects), "\n")
cat("Coefficients:", dim(fit_q5$coefficients), "\n\n")

'''
Perform the log2(data + 1) transform. Then fit a regression model to each sample using population as the outcome. Do this using the 
lm.fit function (hint: dont forget the intercept). What is the dimension of the residual matrix, the effects matrix and the coefficients matrix?
'''

# ---------------------------------------------
# Question 6: What is the effects matrix?
# Explained: The effects matrix contains the projected values on orthogonal basis vectors
# ---------------------------------------------
# Already discussed in code
'''
Perform the log2(data + 1) transform. Then fit a regression model to each sample using population as the outcome. Do this using the 
lm.fit function (hint: dont forget the intercept). What is the effects matrix?
'''

# ---------------------------------------------
# Question 7: lmFit with age as outcome
# ---------------------------------------------

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

'''
Fit many regression models to the expression data where ageis the outcome variable using the 
lmFit from the limma package (hint: you may have to subset the expression data to the samples without missing values of age to get the model to fit). 
What is the coefficient for age for the 1,000th gene? Make a plot of the data and fitted values for this gene. Does the model fit well?
'''

age_idx <- !is.na(pdata_bm$age)
design_q7 <- model.matrix(~ pdata_bm$age[age_idx])
fit_q7 <- lmFit(edata_bm[, age_idx], design_q7)
coef_gene1000 <- fit_q7$coefficients[1000, ]

cat("Q7 - Coefficient for age (gene 1000):", coef_gene1000[2], "\n")
plot(pdata_bm$age[age_idx], edata_bm[1000, age_idx],
     main="Q7 - Gene 1000 vs Age", xlab="Age", ylab="Expression")
abline(lm(edata_bm[1000, age_idx] ~ pdata_bm$age[age_idx]))

# ---------------------------------------------
# Question 8: lmFit with age and tissue.type (overfitting)
# ---------------------------------------------

'''
Fit many regression models to the expression data where age is the outcome variable and 
tissue.type is an adjustment variable using the lmFit function from the limma package 
(hint: you may have to subset the expression data to the samples without missing values of age to get the model to fit). What is wrong with this model?

'''

design_q8 <- model.matrix(~ pdata_bm$age[age_idx] + pdata_bm$tissue.type[age_idx])
tryCatch({
  fit_q8 <- lmFit(edata_bm[, age_idx], design_q8)
}, error = function(e) {
  cat("Q8 - Model fitting error (expected due to overparameterization):", e$message, "\n")
})

# ---------------------------------------------
# Question 9: Why study/population effects are hard to distinguish
# Answer: Each study only measured one population => confounded
# ---------------------------------------------

# ---------------------------------------------
# Question 10: Surrogate variable analysis (SVA)
# ---------------------------------------------

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

'''
Set the seed using the command set.seed(33353) then estimate a single surrogate variable using the 
sva after log2(data + 1) transforming the expression data, removing rows with rowMeans less than 1, 
and treating age as the outcome (hint: you may have to subset the expression data to the samples without 
missing values of age to get the model to fit). What is the correlation between the estimated surrogate for batch and age? Is the surrogate more highly correlated with race or gender?
'''
set.seed(33353)
edata_log <- log2(edata_bm + 1)
edata_filt <- edata_log[rowMeans(edata_log) > 1, ]
age_idx <- !is.na(pdata_bm$age)

mod <- model.matrix(~ pdata_bm$age[age_idx])
mod0 <- model.matrix(~ 1, data = pdata_bm[age_idx, ])
svobj <- sva(edata_filt[, age_idx], mod, mod0)

surr <- svobj$sv[,1]
cor_age <- cor(surr, pdata_bm$age[age_idx])
cor_gender <- cor(surr, as.numeric(as.factor(pdata_bm$gender[age_idx])))
cor_race <- cor(surr, as.numeric(as.factor(pdata_bm$race[age_idx])))

cat("Q10 - Correlation with Age:", cor_age, "\n")
cat("Correlation with Gender:", cor_gender, "\n")
cat("Correlation with Race:", cor_race, "\n")

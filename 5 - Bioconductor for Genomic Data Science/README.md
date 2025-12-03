# Bioconductor for Genomic Data Science (Course 5)

This directory contains my R scripts for the four module quizzes from  
**Course 5 – Bioconductor for Genomic Data Science** in the Johns Hopkins Genomic Data Science Specialization.

These are small, focused scripts – not a full pipeline – and they’re mainly useful as a record of how I worked through the course material in **R/Bioconductor**.

---

## Files

- `quiz_module_1.R`  
 R script with my answers for the **module 1 quiz**.  
 Covers basic Bioconductor setup and introductory operations on common data structures (e.g. importing data, simple manipulations).

- `quiz_module_2.R`  
 Script for the **module 2 quiz**.  
 Implements tasks related to working with genomic ranges and annotations using Bioconductor-style workflows.

- `quiz_module_3.R`  
 Script for the **module 3 quiz**.  
 Contains code for more advanced operations on genomic objects and simple downstream analyses.

- `quiz_module_4.R`  
 Script for the **module 4 quiz**.  
 Implements the final set of exercises, bringing together previous concepts (e.g. operating on expression / genomic data and summarising results).

> Note: the exact content of each script follows the Coursera quiz questions; they are deliberately small and task-specific rather than a general-purpose package.

---

## How to run

1. Open R or RStudio in this folder:

  ```r
  setwd("5 - Bioconductor for Genomic Data Science")

2. Make sure you have R and the relevant Bioconductor packages installed. From R:

if (!require("BiocManager")) {
 install.packages("BiocManager")
}
# Install core Bioconductor packages used in the course
BiocManager::install(c(
 "BiocGenerics",
 "Biobase",
 "GenomicRanges"
 # add others here if a script complains
))  


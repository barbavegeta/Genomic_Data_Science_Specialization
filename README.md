# Genomic Data Science Specialization - Code Exercises

This repository contains my personal solutions, scripts, and notes from the **Johns Hopkins Genomic Data Science Specialization** on Coursera.

The aim of this repo is to document the work I actually did during the courses - from basic sequence processing and algorithms, to command-line workflows, Bioconductor analyses, and introductory statistical genomics - using **Python**, **R**, **Bash**, and **Jupyter Notebooks**.

**Tech stack:** Python · R · Bash · Jupyter Notebooks · Bioconductor · Bowtie2 · BWA · HISAT2 · samtools · bedtools · DESeq2

> This is a learning repository, not a polished production pipeline. It is meant to show my progression through the material and give me a reference for future work in bioinformatics.

---

## Courses covered

This repo includes code and notebooks corresponding to exercises from the following courses in the specialization:

1. **Introduction to Genomic Technologies**
2. **Python for Genomic Data Science**
3. **Algorithms for DNA Sequencing**
4. **Command Line Tools for Genomic Data Science**
5. **Bioconductor for Genomic Data Science**
6. **Statistics for Genomic Data Science**
7. **Genomic Data Science with Galaxy**
8. **Genomic Data Science Capstone**

Not every exercise from every course is represented, but the key coding and command-line components are included.

---

## Repository structure

The exact folder names may evolve, but the structure is organised around courses and shared utilities:

- `3 - Algorithms for DNA Sequencing/`  
  Solutions for the algorithms course: Python and notebook implementations of classic string/sequence algorithms, pattern matching, indexing, and basic read processing.

- `4 - Command line tools for Genomic Data Science/`  
  Shell scripts and command sequences using Unix tools (e.g. `grep`, `awk`, `sed`), plus genomics-specific tools such as `samtools`, `bedtools`, and related utilities for working with FASTQ/BAM/VCF files.

- `5 - Bioconductor for Genomic Data Science/`  
  R scripts and RMarkdown/notebook files using Bioconductor packages for tasks like differential expression, annotation, and basic genomic workflows.

- `6 - Statistics for Genomic Data Science/`  
  R code and notebooks covering statistical concepts applied to genomic data, including basic models, hypothesis testing, and visualisation.

- `Notebooks_commands/`  
  General Jupyter notebooks and command summaries that span multiple courses (e.g. combined notes on tools, small experiments, or scratch work).

- `R/`  
  Shared R helper scripts and utility functions reused across course assignments.

If your local clone has slightly different directory names, the same logic applies: each top-level directory corresponds to a course or shared code.

---

## Getting started

### Prerequisites

You will need:

- A Unix-like environment (Linux, macOS, or WSL on Windows)
- **Python** (≥ 3.10) with Jupyter Notebook or JupyterLab
- **R** (≥ 4.x) with RStudio or another R environment
- For command-line exercises:
  - Standard Unix tools (`grep`, `awk`, `sed`, `sort`, etc.)
  - Genomics tools used in the courses, typically including:
    - `samtools`
    - `bedtools`
    - Other utilities as specified in individual notebooks/scripts

> Note: the exact tool versions are those commonly used during the Coursera course; some commands may need minor adjustment if you are using newer versions.

### Clone the repository

```bash
git clone https://github.com/barbavegeta/Genomic_Data_Science_Specialization.git
cd Genomic_Data_Science_Specialization

# Example - adjust package list as needed
conda create -n genomic-data-science python=3.10 jupyter numpy pandas biopython
conda activate genomic-data-science
jupyter notebook


# Bioconductor
install.packages(c("tidyverse", "data.table"))

if (!require("BiocManager")) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "GenomicRanges",
  "DESeq2",
  "edgeR"
  # add more Bioconductor packages here as used in the scripts
))


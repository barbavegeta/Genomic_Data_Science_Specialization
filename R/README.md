# Shared R / Bioconductor Notes  

This directory contains a set of **R scripts** I wrote while working through the Bioconductor and statistics parts of the Johns Hopkins Genomic Data Science Specialization.  
They are **not** a polished R package. They’re small, focused examples and code snippets showing how to use core Bioconductor classes, packages, and workflows: importing data, working with genomic ranges, handling expression data, and running basic analyses.  
---
## Contents
### Bioconductor setup and general resources
- `Install_Bioconductor.R`  
  Basic Bioconductor installation and update commands.
- `Getting_Data_into_Bioconductor.R`  
  Examples of loading data from files and online resources into Bioconductor objects.
- `Online_Resources.R`  
  Pointers and small examples related to Bioconductor documentation and online help.
### Core R / language features used by Bioconductor
- `R_Base_Types.R`  
  Quick notes on base R data types and structures.
- `R_S4.R`  
  Examples of S4 classes and methods, which underpin many Bioconductor objects.
### Working with sequence data
- `ShortRead.R`  
  Using the **ShortRead** package to read and manipulate sequencing reads (e.g. FASTQ).
- `BSgenome.R`  
  Loading and using whole-genome sequences via **BSgenome**.
- `BSgenome_Views.R`  
  Working with **Views** on BSgenome objects (subsequences, windows).
- `Biostrings.R`  
  Basic operations on DNA/RNA/protein sequences using **Biostrings**.
- `Biostrings_Matching.R`  
  Pattern matching and alignment-style operations with Biostrings.
- `Rsamtools.R`  
  Interacting with BAM/SAM and related formats using **Rsamtools**.
### Genomic ranges and related infrastructure
- `IRanges_Basic.R`  
  Introduction to **IRanges** objects and interval operations.
- `GenomicRanges_GRanges.R`  
  Creating and working with **GRanges** objects.
- `GenomicRanges_GRanges_Usage.R`  
  More applied examples of operations on GRanges (subsetting, overlaps, etc.).
- `GenomicRanges_Lists.R`  
  Using list-like containers of ranges (e.g. **GRangesList**).
- `GenomicRanges_Rle.R`  
  Run-length encoded vectors (**Rle**) and how they fit into genomic workflows.
- `GenomicRanges_seqinfo.R`  
  Managing genome/sequence metadata (`seqinfo`) for GenomicRanges objects.
- `rtracklayer_Import.R`  
  Importing/exporting genome annotation tracks (e.g. BED, GFF, WIG) via **rtracklayer**.
- `GenomicFeatures.R`  
  Building and using TxDb / transcript-level annotation with **GenomicFeatures**.
### Annotation and external resources
- `AnnotationHub.R`  
  Accessing annotation resources via **AnnotationHub**.
- `Usecase_AnnotationHub_GRanges.R`  
  Practical example combining AnnotationHub data with GRanges.
- `GEOquery.R`  
  Using **GEOquery** to pull datasets from GEO and convert them into Bioconductor-friendly objects.
- `biomaRt.R`  
  Annotation and identifier mapping using **biomaRt**.
### Expression data and downstream analyses
- `ExpressionSet.R`  
  Working with the **ExpressionSet** class for microarray / expression data.
- `SummarizedExperiment.R`  
  Using **SummarizedExperiment** for matrix-like assays plus row/column metadata.
- `Count_Based_RNAseq.R`  
  Basic count-based RNA-seq workflow examples (e.g. DE-style input matrices).
- `limma.R`  
  Using **limma** for linear modelling of expression data (microarray-style or voom-transformed counts).
- `oligo.R`  
  Examples using **oligo** for processing Affymetrix and other microarray platforms.
- `minfi.R`  
  Basics of using **minfi** for methylation array data.
---
## How to use these scripts
1. Start R or RStudio in this directory:
   ```r
   setwd("R")
   ```
2. Make sure you have Bioconductor set up:
   ```r
   if (!require("BiocManager")) {
     install.packages("BiocManager")
   }
   ```
3. Install packages as needed. For example:
   ```r
   BiocManager::install(c(
     "Biostrings",
     "ShortRead",
     "GenomicRanges",
     "IRanges",
     "Rsamtools",
     "AnnotationHub",
     "GenomicFeatures",
     "rtracklayer",
     "GEOquery",
     "SummarizedExperiment",
     "limma",
     "oligo",
     "minfi"
     # and others referenced in the scripts if they are missing
   ))
   ```
4. Open or source any script of interest. For example:
   ```r
   source("GenomicRanges_GRanges.R")
   ```
   or open the file in RStudio and run it line by line.
Most scripts are self-contained demonstrations of a specific package or concept and can be run independently once the relevant packages are installed.
---
## Notes
- These files are **reference / practice scripts**, not an integrated analysis pipeline.
- They’re useful as a quick reminder of how to perform common tasks in Bioconductor:

  importing data, working with genomic ranges, accessing annotation resources, and handling expression or methylation data.


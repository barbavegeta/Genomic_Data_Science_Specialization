# Command Line Tools for Genomic Data Science (Course 4)

This directory contains my work for **Course 4 – Command Line Tools for Genomic Data Science** from the Johns Hopkins Genomic Data Science Specialization.

The focus here is on using the Unix command line and common genomics tools to manipulate sequencing data and related file formats. I’ve deliberately kept only:

- My own scripts and command files
- Small example inputs
- A few final result files (e.g. VCFs, count matrices)

Large raw datasets, indexes, and heavy intermediate BAM/FASTQ files are **not** tracked in this repository.

---

## Structure

- `Week-1/`  
  Basic Unix and file-handling exam work:
  - `Module-1-Exam.sh` – script solving the module 1 exam tasks.
  - `apple.*` – small toy genome and condition files used by the script.

- `Week-2/`  
  Command-line work with alignments and genomic intervals:
  - `Module-2-Exam.sh` + `Module-2-Exam.md`, `Module-2-Quiz.md`.
  - Small annotation and interval files (`*.bed`, `*.gtf`, `lengths`), used to practise `samtools` and `bedtools`.
  - Large BAM files are **not** kept here; see comments in the scripts for expected inputs.

- `Week-3/`  
  Variant-calling pipeline:
  - `project_module_3.sh` – bowtie2 + samtools + bcftools pipeline for calling variants on the provided *Arabidopsis thaliana* data.
  - `out.final.vcf`, `full_wu_0_A_variants.vcf`, `local_wu_0_A_variants.vcf` – example VCF outputs.
  - Intermediate BAMs, indexes, and FASTQ files are intentionally omitted from version control.

- `Week-4/`  
  RNA-seq pipeline work:
  - `project_module_4.sh`, `run_rnaseq_analysis.sh`, `tophat.sh` – command-line pipelines using Tophat/Cufflinks/Cuffdiff and/or HISAT2 for RNA-seq.
  - `analysis_r.R` – R script to analyse expression/count matrices (e.g. log2 fold changes, plots).
  - `prepDE.py` – helper script to extract count matrices from Cufflinks output.
  - `gencommand_proj4/commands/` – command files (`com.tophat`, `com.cufflinks`, `com.cuffdiff`, `GTFs.txt`) used in the project.
  - Selected small reference files (e.g. `athal_chr.fa`, `athal_genes.gtf`, `merged.gtf`, count matrices) may be present; large BAM/FASTQ/index files are not kept.

- Top level:
  - `commands.sh` – assorted command snippets used throughout the course.
  - `bwa_bowtie2_commands.sh` – notes and commands for running BWA/Bowtie2 on the IL-2 and NA12814 examples.
  - `IL-2_mRNA.fasta` – small example FASTA used to demonstrate indexing and alignment commands (index artefacts are not tracked).

Modern HISAT2-based pipeline: see 'gencommand_proj4/'

Legacy TopHat/Cufflinks/Cuffdiff example: see 'legacy_tophat_cufflinks_pipeline/' (kept for historical / maintenance context).

---

## Tools used

Across these scripts and command files, I worked with:

- Core Unix tools: `grep`, `awk`, `sed`, `cut`, `sort`, `uniq`, `wc`, `tar`, `gzip`
- Alignment and indexing: `bowtie2`, `bwa`, `tophat`, `hisat2`
- File formats and utilities: `samtools`, `bcftools`, `bedtools`
- RNA-seq tooling: `cufflinks`, `cuffdiff`
- SRA tools: `prefetch` (for downloading from SRA)
- R + Python for downstream analysis (`analysis_r.R`, `prepDE.py`)

The goal here is to demonstrate that I can correctly construct and chain these tools on the command line and interpret the resulting BAM/VCF/GTF/count files.

---

## How to run

1. Use a Unix-like environment (Linux, macOS, or WSL) with the tools above installed and on your `PATH`.

2. Provide or recreate the input data expected by each script:
   - For smaller examples (e.g. `apple.*`, `IL-2_mRNA.fasta`), files are included here.
   - For larger datasets (BAM/FASTQ), refer to the Coursera course resources or public archives and adapt the paths in the scripts.

3. Make scripts executable if needed, e.g.:

   ```bash
   chmod +x commands.sh bwa_bowtie2_commands.sh
   chmod +x Week-1/Module-1-Exam.sh
   chmod +x Week-2/Module-2-Exam.sh
   chmod +x Week-3/project_module_3.sh
   chmod +x Week-4/project_module_4.sh Week-4/run_rnaseq_analysis.sh Week-4/tophat.sh


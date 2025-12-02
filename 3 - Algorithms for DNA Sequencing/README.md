# Algorithms for DNA Sequencing (Course 3)

This folder contains my solutions and notes for **Course 3 – Algorithms for DNA Sequencing** from the Johns Hopkins Genomic Data Science Specialization.

The focus here is on implementing classic algorithms and data structures for analysing DNA sequencing data using **Python** and **Jupyter Notebooks**, with small example read sets.

---

## Contents

### Notebooks

- `Algorithms_Week_1.ipynb`  
- `Algorithms_Week_2.ipynb`  
- `Algorithms_Week_3.ipynb`  
- `Algorithms_Week_4.ipynb`  

Weekly assignment notebooks implementing core algorithms and exercises from the course.

- `Notes_Algorithms_for_DNA_sequencing_week_1.ipynb`  
- `Notes_Algorithms_for_DNA_sequencing_week_2.ipynb`  
- `Notes_Algorithms_for_DNA_sequencing_week_3.ipynb`  
- `Notes_Algorithms_for_DNA_sequencing_week_4.ipynb`  

  Lecture notes and scratch work for each week.

- `hw3_overlap_all.ipynb`  

Additional work on overlap-based assembly and related tasks.

### Python scripts

- `bm_preproc.py` – preprocessing and helper functions for Boyer–Moore or similar string-matching routines used in the assignments.  
- `kmer_index.py` – helpers for building and querying k-mer–based indexes over sequences.

### Example data

Small example datasets provided with the course:

- `ERR037900_1.first1000.fastq`  
- `ERR266411_1.first1000.fastq`  
- `ERR266411_1.for_asm.fastq`  
- `SRR835775_1.first1000.fastq`  
- `ads1_week4_reads.fq`  

Short read sets used in the weekly assignments and assembly exercises.

Reference sequences:

- `chr1.GRCh38.excerpt.fasta`  
- `lambda_virus.fa`  
- `phix.fa`  

FASTA excerpts and small viral genomes used for alignment and algorithm testing.

### De Bruijn graph materials

- `debruijn_graph.sh` – shell script to generate a simple de Bruijn graph or related outputs for one of the assignments.  
- `debruijn_graph.png` – corresponding visualisation of the graph.

---

## How to run

1. Activate the Python environment described in the root `README.md` (e.g. conda env `genomic-data-science`).

2. From the repository root, start Jupyter:

   ```bash
   jupyter notebook

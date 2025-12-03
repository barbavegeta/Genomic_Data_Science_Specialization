# Notebooks & Command Notes (Algorithms / String Processing)  

This directory contains a set of **Jupyter notebooks** I used to practise and internalise the string/sequence algorithms from the Genomic Data Science Specialization (especially the *Algorithms for DNA Sequencing* material).  
They are not graded assignments, but my own working notes and implementations of core ideas: string matching, indexing, approximate matching, alignment, and simple assembly.  

---
## Files
### Section 1 – Basic string and read handling
- `1.01_StringBasics.ipynb`  
Basic Python string operations in the context of DNA sequences.
- `1.02_ManipulatingDNAStrings.ipynb`  
Simple utilities for working with DNA strings (reverse complement, GC content, slicing, etc.).
- `1.03_ParsingRefGenome.ipynb`  
Parsing reference genomes from FASTA and storing them in memory-efficient structures.
- `1.04_WorkingWithSequencingReads.ipynb`  
Reading FASTQ/FASTA, iterating over reads, and basic sanity checks.
- `1.05_AnalyzingReadsByPosition.ipynb`  
Counting and summarising reads by genomic position.
- `1.06_NaiveExactMatching-MatchingArtificialReads.ipynb`  
Naive exact pattern matching on artificial read sets (toy examples).
- `1.07_NaiveExactMatching-MatchingRealReads.ipynb`  
Applying naive exact matching to more realistic read sets to see performance limits.
### Section 2 – Index-based matching
- `2.01_BoyerMoore.ipynb`  
Implementation and experiments with the Boyer–Moore string matching algorithm.
- `2.02_SubstringIndex.ipynb`  
Building and querying substring / k-mer indices for faster pattern matching.
- `2.03_ApproximateMatching.ipynb`  
Approximate matching (with mismatches) using combinations of indexing and local checks.
### Section 3 – Edit distance, alignment, overlaps
- `3.01_EditDistanceDP.ipynb`  
Dynamic programming implementation of edit distance (Levenshtein) and related ideas.
- `3.02_GlobalAlignment.ipynb`  
Global alignment (Needleman–Wunsch-style) via dynamic programming.
- `3.03_FindingOverlaps.ipynb`  
Finding overlaps between reads using simple substring-based methods.
- `3.04_FindingAllOverlaps.ipynb`  
Extending the overlap-finding logic to all read pairs; exploring complexity and optimisations.
### Section 4 – Assembly-style problems
- `4.01_ShortestCommonSuperstring.ipynb`  
Shortest common superstring formulations for toy read sets.
- `4.02_GreedySCS.ipynb`  
Greedy shortest-common-superstring algorithms and experiments.
- `4.03_DeBruijn.ipynb`  
De Bruijn graph construction for k-mers and basic graph-based assembly intuition.
- `README.md`  
This file.
---
## How to run
1. Make sure you have **Python** (≥ 3.10) and Jupyter installed, ideally in the same environment you use for the rest of the Genomic Data Science code.
 Example (using conda):
  ```bash
 conda create -n genomic-data-science python=3.10 jupyter numpy
 conda activate genomic-data-science
  ```
2. From the repository root (or from this folder), start Jupyter:
  ```bash
  jupyter notebook
  ```
3. In the Jupyter UI, navigate to:
  ```text
  Notebooks_commands/
  ```
- and open any `*.ipynb` notebook.
Each notebook is self-contained; when they expect external data, the paths and file names are referenced in comments or in the first cells.
---
## Notes
- These notebooks are **learning notes**, not production-ready code.  
- They are useful as a quick reference for classic string/sequence algorithms that crop up in bioinformatics (pattern matching, indexing, alignment, simple assembly).  
- For polished pipelines, see the course-specific folders (e.g. `3 - Algorithms for DNA Sequencing`).

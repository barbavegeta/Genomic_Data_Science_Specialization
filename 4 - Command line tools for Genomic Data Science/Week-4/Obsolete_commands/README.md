# Legacy TopHat/Cufflinks/Cuffdiff RNA-seq pipeline

This folder contains a legacy RNA-seq workflow using **TopHat**, **Cufflinks**, and **Cuffdiff**.

These tools are now largely superseded by modern aligners and quantification tools
(e.g. **HISAT2**, **STAR**, **StringTie**, **featureCounts**, **Salmon**, **DESeq2**),
but I keep this example to show that I have worked with older-style pipelines and can
understand / maintain them when needed.

Kept here:
- `run_rnaseq_analysis.sh` – main shell pipeline
- `tophat.sh` (if you have it here) – alignment step
- `prepDE.py`, `analysis_r.R` – extraction and analysis of count matrices
- `GTFs.txt`, `cuffdiff_gtfs.txt` – configuration of GTFs
- Small reference files (`athal_chr.fa`, `athal_genes.gtf`, `merged.gtf`)
- Final count matrices and log2 fold-change tables (`gene_count_matrix.csv`,
  `transcript_count_matrix.csv`, `normalized_counts.csv`, `manual_log2_fold_change_results.csv`)
- Minimal Cufflinks and Tophat outputs (e.g. `transcripts.gtf`, `align_summary.txt`)

Large BAMs, index files, and extensive intermediate outputs have been removed to keep
the repository lean and focused on code and key results.
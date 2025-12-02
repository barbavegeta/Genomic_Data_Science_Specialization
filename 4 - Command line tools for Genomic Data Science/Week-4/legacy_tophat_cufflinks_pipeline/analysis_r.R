# Load necessary libraries
# If you don't have them installed, run:
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("readr") # For read_csv
# install.packages("dplyr") # For potential data manipulation, useful later

library(DESeq2)
library(readr) # Recommended for fast and reliable CSV reading
library(dplyr) # Good practice for data manipulation

# --- 1. Set Working Directory ---
# IMPORTANT: Replace this path with the actual directory where your 'gene_count_matrix.csv' is located.
# Use forward slashes (/) or double backslashes (\\) for paths on Windows.
setwd("C:/Users/barba/OneDrive - University College London/Bioinformatics/Coursera/Genomic Data Science Specialization/4 - Command line tools for Genomic Data Science/Week-4/test2")

# --- 2. Load Count Data ---
# It's crucial that 'gene_count_matrix.csv' has gene IDs in the first column
# and then sample counts in subsequent columns, with sample names as headers.
# Example:
# gene_id,Day16,Day8
# gene1,100,120
# gene2,50,45
# ...

# Read the gene count matrix.
# 'col_names = TRUE' (default) ensures the first row is treated as headers.
counts_data <- read_csv("gene_count_matrix.csv")

# Verify the structure after reading (optional, but good for debugging)
print(head(counts_data))
print(colnames(counts_data))

# Set the gene IDs as row names and remove the gene_id column from the data frame
# DESeq2 expects count data with gene IDs as row names.
if ("gene_id" %in% colnames(counts_data)) {
  rownames(counts_data) <- counts_data$gene_id
  counts_data$gene_id <- NULL
} else {
  stop("Error: 'gene_id' column not found. Please check your gene_count_matrix.csv file.")
}

# Convert to a matrix and ensure counts are integers
# DESeq2 requires integer counts. round() is good, but make sure no negative values
# or very small floats might cause issues if not strictly integers.
counts_matrix <- as.matrix(counts_data)
counts_matrix <- round(counts_matrix)

# Important check: Ensure all counts are non-negative integers
if (any(counts_matrix < 0) || any(counts_matrix != floor(counts_matrix))) {
  warning("Count matrix contains non-integer or negative values. This might cause issues with DESeq2.")
}


# --- 3. Create ColData (Sample Information) ---
# This data frame describes your samples and their experimental conditions.
# The row names MUST exactly match the column names of your counts_matrix.
# Define your sample names in the desired order for your coldata
# (e.g., control first, then treatment if you want Day8 as baseline).
sample_names <- c("Day8", "Day16")
# Explicitly define levels for the factor to control the order for comparisons
condition <- factor(c("Day8", "Day16"), levels = c("Day8", "Day16"))

coldata <- data.frame(
  row.names = sample_names,
  condition = condition # This is your primary experimental variable
)

# --- CRITICAL FIX: Reorder counts_matrix columns to match coldata row names ---
# This ensures the sample order is consistent, which is crucial for DESeq2.
# Your 'read_csv' output showed 'Day16, Day8' but we want 'Day8, Day16' in colData.
counts_matrix <- counts_matrix[, rownames(coldata)]

# Now, the crucial check will pass if the above reordering was correct
if (!all(rownames(coldata) == colnames(counts_matrix))) {
  stop("FATAL ERROR: Sample names in colData still do not match column names in counts_matrix after reordering. Review your data and script.")
} else {
  print("Sample order in colData and countData confirmed to match.")
}


# --- 4. Construct DESeqDataSet Object for No-Replicate Analysis ---
# IMPORTANT: With only one sample per condition, DESeq2 CANNOT estimate dispersion
# and therefore CANNOT perform statistical differential expression testing (no p-values).
# We set the design to ~1 to tell DESeq2 there are no experimental groups to compare
# for statistical testing purposes (only for size factor estimation).
print("Constructing DESeqDataSet object for no-replicate analysis...")
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata,
  design = ~ 1 # This design is used when you have no replicates for statistical comparison
)
print("DESeqDataSet object created.")


# --- 5. Estimate Size Factors and Get Normalized Counts ---
# Even without replicates for statistical testing, DESeq2 can still perform
# size factor estimation to normalize for library size differences.
print("Estimating size factors for normalization...")
dds <- estimateSizeFactors(dds)
print("Size factor estimation complete.")

# Get normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
print("Head of normalized counts:")
print(head(normalized_counts))


# --- 6. Calculate Log2 Fold Change Manually (No P-values) ---
# Since DESeq2 cannot perform statistical testing without replicates,
# we manually calculate the log2 fold change based on normalized counts.
# Add a small pseudo-count to avoid log(0) issues for genes with zero counts.
pseudo_count <- 1

# Calculate fold change for Day16 vs Day8
# A positive log2FoldChange means higher expression in Day16 relative to Day8.
log2_fold_change <- log2( (normalized_counts[, "Day16"] + pseudo_count) /
                            (normalized_counts[, "Day8"] + pseudo_count) )

# Create a data frame for these manual results
# This table will NOT have 'stat', 'pvalue', or 'padj' columns.
manual_results <- data.frame(
  baseMean = rowMeans(normalized_counts),
  log2FoldChange = log2_fold_change,
  row.names = rownames(normalized_counts) # Gene IDs as row names
)

# Order results by the absolute log2FoldChange to see the largest changes first
manual_results_ordered <- manual_results[order(abs(manual_results$log2FoldChange), decreasing = TRUE), ]

print("Manual Log2 Fold Change Results (TOP 10 - no p-values/significance):")
print(head(manual_results_ordered, 10)) # Print top 10 for better overview

# --- 7. Save Results ---
# Save the normalized counts
write.csv(normalized_counts, file="normalized_counts.csv", row.names=TRUE)
print("Normalized counts saved to 'normalized_counts.csv'")

# Save the manual log2 fold change results
write.csv(manual_results_ordered, file="manual_log2_fold_change_results.csv", row.names=TRUE)
print("Manual log2 fold change results saved to 'manual_log2_fold_change_results.csv'")


# --- 8. Important Note on Visualization ---
# You CANNOT use plotMA() from DESeq2 directly because the full DESeq() function
# (which calculates statistical results) was not run due to lack of replicates.
# If you want to visualize these non-statistical fold changes, you would need
# to create a custom plot using base R or ggplot2 based on 'manual_results_ordered'.

# Example of a very basic custom MA-like plot (no significance coloring)
plot(manual_results$baseMean, manual_results$log2FoldChange,
     log="x", # Log scale for baseMean (X-axis)
     xlab="Mean of Normalized Counts",
     ylab="Log2 Fold Change (Day16 vs Day8)",
     main="MA-plot (No Replicates - Descriptive Only)",
     pch=20, cex=0.7, col="darkblue")
abline(h=0, col="red", lty=2) # Add a horizontal line at 0 fold change

print("Analysis script finished. Remember results are descriptive only due to lack of replicates.")


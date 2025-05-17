#!/usr/bin/env Rscript

# WGCNA Analysis for Gene Expression Data with Two Genotypes and Three Developmental Stages
# This script performs Weighted Gene Co-expression Network Analysis on gene expression data

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("WGCNA", quietly = TRUE)) BiocManager::install("WGCNA", dependencies = TRUE)
if (!requireNamespace("flashClust", quietly = TRUE)) install.packages("flashClust")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")

library(WGCNA)
library(flashClust)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(DESeq2)
library(edgeR)
library(limma)

# Set WGCNA parameters
options(stringsAsFactors = FALSE)
enableWGCNAThreads() # Enable multi-threading

# Set working directory
# setwd("your/working/directory") # Uncomment and set your working directory

# 1. Data Loading and Preprocessing ------------------------------------------

# File paths - replace with your actual file paths
expression_file <- "expression_data.txt"
metadata_file <- "metadata.csv"

# Load expression data
expr_data <- read.table(expression_file, header = TRUE, sep = "\t", quote = "", row.names = 1)

# Fix column names by removing the "X" prefix that R automatically adds to numeric column names
colnames(expr_data) <- sub("^X", "", colnames(expr_data))

# Load metadata
metadata <- read.csv(metadata_file, header = TRUE)

# Sample ordering check and correction
# Ensure samples in expression data and metadata are in the same order
sample_order_check <- function() {
  # Check if the order of samples matches between metadata and expression data
  if (!identical(as.character(metadata$sample_number), colnames(expr_data))) {
    cat("Warning: Sample order in metadata does not match expression data.\n")
    cat("Reordering expression data columns to match metadata...\n")
    
    # Create a new expression matrix with columns ordered to match metadata
    ordered_expr_data <- expr_data[, match(as.character(metadata$sample_number), colnames(expr_data))]
    
    # Check if any columns couldn't be matched (will be NA)
    if (any(is.na(colnames(ordered_expr_data)))) {
      cat("ERROR: Some samples in metadata couldn't be found in expression data.\n")
      missing_samples <- as.character(metadata$sample_number)[is.na(match(as.character(metadata$sample_number), colnames(expr_data)))]
      cat("Missing samples:", paste(missing_samples, collapse=", "), "\n")
      stop("Cannot proceed with analysis due to missing samples.")
    }
    
    # Replace the original expression data with the ordered one
    expr_data <<- ordered_expr_data
    cat("Expression data columns successfully reordered to match metadata.\n")
  } else {
    cat("Sample order in metadata and expression data match. No reordering needed.\n")
  }
}

# Call the sample order check function
sample_order_check()

# Data preprocessing and normalization --------------------------------------
# Check data type and convert to numeric if needed
cat("Checking data types...\n")

# Function to check if a matrix is numeric
check_and_convert_numeric <- function(data_matrix) {
  # Check if data is already numeric
  if (!is.numeric(as.matrix(data_matrix))) {
    cat("WARNING: Expression data contains non-numeric values.\n")
    
    # Print a sample of the problematic data
    cat("Sample of the data:\n")
    print(data_matrix[1:5, 1:5])
    
    # Try to convert to numeric
    cat("Attempting to convert data to numeric...\n")
    
    # Function to safely convert a single value to numeric
    safe_as_numeric <- function(x) {
      result <- suppressWarnings(as.numeric(as.character(x)))
      if (is.na(result) && !is.na(x)) {
        cat("Conversion failed for value:", x, "\n")
      }
      return(result)
    }
    
    # Apply conversion to all elements
    numeric_matrix <- as.data.frame(lapply(data_matrix, function(col) {
      sapply(col, safe_as_numeric)
    }))
    
    rownames(numeric_matrix) <- rownames(data_matrix)
    
    # Check for NAs introduced by conversion
    na_count <- sum(is.na(numeric_matrix))
    if (na_count > 0) {
      cat("WARNING:", na_count, "values could not be converted to numeric and became NA.\n")
      cat("This may affect analysis results.\n")
      
      # Optionally, remove rows with NAs
      rows_with_na <- rowSums(is.na(numeric_matrix)) > 0
      if (sum(rows_with_na) > 0) {
        cat("Removing", sum(rows_with_na), "rows containing NA values...\n")
        numeric_matrix <- numeric_matrix[!rows_with_na, ]
      }
    }
    
    cat("Conversion complete. Proceeding with numeric data.\n")
    return(numeric_matrix)
  } else {
    cat("Data is already numeric. No conversion needed.\n")
    return(data_matrix)
  }
}

# Apply the conversion function
expr_data <- check_and_convert_numeric(expr_data)

# Detect if we're dealing with count data
tryCatch({
  is_count_data <- all(expr_data == round(expr_data), na.rm = TRUE) && 
                  all(expr_data >= 0, na.rm = TRUE) && 
                  mean(as.matrix(expr_data), na.rm = TRUE) > 10  # Heuristic for count data
  
  cat("Is count data check:", is_count_data, "\n")
  cat("Mean value:", mean(as.matrix(expr_data), na.rm = TRUE), "\n")
}, error = function(e) {
  cat("Error in count data detection:", e$message, "\n")
  cat("Assuming data is count data and proceeding...\n")
  is_count_data <- TRUE
})

# Simple alternative normalization function 
# (used if DESeq2 or edgeR methods fail)
simple_normalization <- function(count_matrix) {
  # First, perform CPM normalization
  lib_sizes <- colSums(count_matrix)
  cpm_matrix <- t(t(count_matrix) / lib_sizes) * 1e6
  
  # Then apply log2 transformation
  log_matrix <- log2(cpm_matrix + 1)
  
  # Return normalized data
  return(log_matrix)
}

# Set a flag to track if we should skip complex normalizations
skip_deseq <- FALSE

if (exists("skip_deseq") && skip_deseq) {
  cat("Using simple CPM + log2 normalization instead of DESeq2/edgeR...\n")
  normalized_data <- simple_normalization(as.matrix(expr_data))
  
  # Create basic visualization
  pdf("simple_normalization.pdf", width = 12, height = 9)
  par(mfrow = c(2, 2))
  
  # Raw data
  boxplot(log2(as.matrix(expr_data) + 1), main = "Raw counts (log2)", 
          xlab = "Samples", ylab = "log2(counts + 1)",
          col = rainbow(ncol(expr_data)), las = 2, cex.axis = 0.7)
  
  # CPM+log2 normalized
  boxplot(normalized_data, main = "CPM + log2 normalized", 
          xlab = "Samples", ylab = "log2(CPM + 1)",
          col = rainbow(ncol(normalized_data)), las = 2, cex.axis = 0.7)
  
  # Density plots
  plot(density(log2(as.matrix(expr_data)[, 1] + 1)), col = 1, main = "Density plots", 
       xlab = "Expression values", ylab = "Density")
  for (i in 2:min(10, ncol(expr_data))) {
    lines(density(log2(as.matrix(expr_data)[, i] + 1)), col = i)
  }
  lines(density(normalized_data[, 1]), col = 1, lty = 2)
  for (i in 2:min(10, ncol(normalized_data))) {
    lines(density(normalized_data[, i]), col = i, lty = 2)
  }
  legend("topright", c("Raw (solid)", "Normalized (dashed)"), lty = 1:2)
  
  # MDS plot
  mds <- cmdscale(dist(t(normalized_data)), k = 2)
  plot(mds, main = "MDS Plot of Samples", 
       xlab = "Dimension 1", ylab = "Dimension 2",
       col = as.numeric(factor(sample_info$condition)), 
       pch = 16)
  legend("topright", levels(factor(sample_info$condition)), 
         col = 1:length(levels(factor(sample_info$condition))), pch = 16)
  
  dev.off()
  
  # Save normalized data
  write.csv(normalized_data, file = "simple_normalized_data.csv")
  
  cat("Simple normalization complete. Using log2(CPM+1) data for WGCNA.\n")
  
  # Use normalized data for subsequent analysis
  expr_data <- normalized_data
}

if (is_count_data) {
  cat("Detected count data, performing RNA-seq normalization...\n")
  
  # Option 1: Variance Stabilizing Transformation (VST) from DESeq2
  # Create DESeq2 object
  sample_info <- data.frame(
    row.names = as.character(colnames(expr_data)),
    condition = paste(metadata$genotype, metadata$development_stage, sep = "_")
  )
  
  # Make sure row.names of sample_info match colnames of expr_data
  if (!identical(rownames(sample_info), colnames(expr_data))) {
    cat("Warning: Sample names in metadata don't match expression data column names.\n")
    cat("Attempting to match sample names...\n")
    
    # Create properly ordered sample_info
    sample_info <- data.frame(
      row.names = colnames(expr_data),
      condition = paste(metadata$genotype[match(colnames(expr_data), as.character(metadata$sample_number))], 
                       metadata$development_stage[match(colnames(expr_data), as.character(metadata$sample_number))], 
                       sep = "_")
    )
  }
  
  cat("Creating DESeq2 dataset...\n")
  
  # Try to create DESeq2 object with error handling
  tryCatch({
    dds <- DESeqDataSetFromMatrix(
      countData = round(as.matrix(expr_data)), # Ensure integers for DESeq2
      colData = sample_info,
      design = ~ condition
    )
    
    # Output dimensions for debugging
    cat("DESeq2 dataset dimensions:", dim(dds)[1], "genes,", dim(dds)[2], "samples\n")
    
    # Filter low count genes
    keep <- rowSums(counts(dds) >= 10) >= 5  # At least 5 samples with at least 10 reads
    
    # Check if keep vector is empty (would happen if conversion issues persist)
    if (sum(keep) == 0) {
      cat("WARNING: No genes passed the low count filtering. This suggests data type issues.\n")
      cat("Proceeding with all genes and basic normalization...\n")
      
      # Skip DESeq2 and just do a simple log transformation
      normalized_data <- log2(as.matrix(expr_data) + 1)
    } else {
      dds <- dds[keep,]
      
      # Variance stabilizing transformation
      cat("Applying variance stabilizing transformation...\n")
      
      # Try VST, but have fallback options
      tryCatch({
        vsd <- vst(dds, blind = TRUE)
        normalized_data <- assay(vsd)
      }, error = function(e) {
        cat("VST failed with error:", e$message, "\n")
        cat("Trying rlog transformation instead...\n")
        
        tryCatch({
          rld <- rlog(dds, blind = TRUE)
          normalized_data <- assay(rld)
        }, error = function(e) {
          cat("rlog also failed with error:", e$message, "\n")
          cat("Falling back to simple log2 transformation...\n")
          
          normalized_data <- log2(as.matrix(counts(dds, normalized=TRUE)) + 1)
        })
      })
      
      # Output pre-normalization stats
      cat("Before normalization:\n")
      cat("  Mean count per gene:", mean(rowMeans(expr_data)), "\n")
      cat("  Median count per gene:", median(rowMeans(expr_data)), "\n")
      cat("  Mean CV per gene:", mean(apply(expr_data, 1, function(x) sd(x)/mean(x))), "\n")
      
      # Output post-normalization stats
      cat("After normalization:\n")
      cat("  Mean value per gene:", mean(rowMeans(normalized_data)), "\n")
      cat("  Median value per gene:", median(rowMeans(normalized_data)), "\n")
      cat("  Mean CV per gene:", mean(apply(normalized_data, 1, function(x) sd(x)/mean(x))), "\n")
      
      # Option 2: TMM normalization with limma-voom (as an alternative)
      # Create DGEList object
      dge <- DGEList(counts = as.matrix(expr_data))
      
      # Filter low expressed genes
      keep <- filterByExpr(dge, group = sample_info$condition)
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      
      # Calculate normalization factors
      dge <- calcNormFactors(dge, method = "TMM")
      
      # Apply voom transformation
      v <- voom(dge, design = model.matrix(~ 0 + sample_info$condition), plot = FALSE)
      
      # Store TMM+voom normalized data as an alternative
      tmm_voom_data <- v$E
      
      # Create comparison plots of normalization methods
      pdf("normalization_comparison.pdf", width = 12, height = 9)
      par(mfrow = c(2, 2))
      
      # Raw data
      boxplot(log2(as.matrix(expr_data) + 1), main = "Raw counts (log2)", 
              xlab = "Samples", ylab = "log2(counts + 1)",
              col = rainbow(ncol(expr_data)), las = 2, cex.axis = 0.7)
      
      # VST normalized
      boxplot(normalized_data, main = "VST normalized", 
              xlab = "Samples", ylab = "VST values",
              col = rainbow(ncol(normalized_data)), las = 2, cex.axis = 0.7)
      
      # TMM+voom normalized
      boxplot(tmm_voom_data, main = "TMM+voom normalized", 
              xlab = "Samples", ylab = "log2 CPM",
              col = rainbow(ncol(tmm_voom_data)), las = 2, cex.axis = 0.7)
      
      # Density plots
      plot(density(log2(as.matrix(expr_data)[, 1] + 1)), col = 1, main = "Density plots", 
           xlab = "Expression values", ylab = "Density")
      for (i in 2:min(10, ncol(expr_data))) {
        lines(density(log2(as.matrix(expr_data)[, i] + 1)), col = i)
      }
      lines(density(normalized_data[, 1]), col = 1, lty = 2)
      for (i in 2:min(10, ncol(normalized_data))) {
        lines(density(normalized_data[, i]), col = i, lty = 2)
      }
      legend("topright", c("Raw (solid)", "Normalized (dashed)"), lty = 1:2)
      
      dev.off()
      
      # Use VST normalized data for WGCNA
      expr_data <- normalized_data
      
      # Save normalized data
      write.csv(normalized_data, file = "vst_normalized_data.csv")
      write.csv(tmm_voom_data, file = "tmm_voom_normalized_data.csv")
      
      cat("Normalization complete. Using variance stabilized data for WGCNA.\n")
    }
  }, error = function(e) {
    cat("Error in DESeq2 processing:", e$message, "\n")
    cat("Falling back to simple normalization...\n")
    
    # Set a flag to skip DESeq2 processing
    assign("skip_deseq", TRUE, envir = .GlobalEnv)
  })
} else {
  cat("Data does not appear to be raw counts. Applying standard scaling...\n")
  # Apply log2 transformation if the data is not already log-transformed
  if (max(expr_data) > 30) {  # Heuristic for detecting non-logged data
    cat("Applying log2 transformation...\n")
    expr_data <- log2(expr_data + 1)
  }
}
  
  # Filter low count genes
  keep <- rowSums(counts(dds) >= 10) >= 5  # At least 5 samples with at least 10 reads
  
  # Check if keep vector is empty (would happen if conversion issues persist)
  if (sum(keep) == 0) {
    cat("WARNING: No genes passed the low count filtering. This suggests data type issues.\n")
    cat("Proceeding with all genes and basic normalization...\n")
    
    # Skip DESeq2 and just do a simple log transformation
    normalized_data <- log2(as.matrix(expr_data) + 1)
  } else {
    dds <- dds[keep,]
    
    # Variance stabilizing transformation
    cat("Applying variance stabilizing transformation...\n")
    
    # Try VST, but have fallback options
    tryCatch({
      vsd <- vst(dds, blind = TRUE)
      normalized_data <- assay(vsd)
    }, error = function(e) {
      cat("VST failed with error:", e$message, "\n")
      cat("Trying rlog transformation instead...\n")
      
      tryCatch({
        rld <- rlog(dds, blind = TRUE)
        normalized_data <- assay(rld)
      }, error = function(e) {
        cat("rlog also failed with error:", e$message, "\n")
        cat("Falling back to simple log2 transformation...\n")
        
        normalized_data <- log2(as.matrix(counts(dds, normalized=TRUE)) + 1)
      })
    })
  
  # Output pre-normalization stats
  cat("Before normalization:\n")
  cat("  Mean count per gene:", mean(rowMeans(expr_data)), "\n")
  cat("  Median count per gene:", median(rowMeans(expr_data)), "\n")
  cat("  Mean CV per gene:", mean(apply(expr_data, 1, function(x) sd(x)/mean(x))), "\n")
  
  # Output post-normalization stats
  cat("After normalization:\n")
  cat("  Mean VST value per gene:", mean(rowMeans(normalized_data)), "\n")
  cat("  Median VST value per gene:", median(rowMeans(normalized_data)), "\n")
  cat("  Mean CV per gene:", mean(apply(normalized_data, 1, function(x) sd(x)/mean(x))), "\n")
  
  # Option 2: TMM normalization with limma-voom (as an alternative)
  # Create DGEList object
  dge <- DGEList(counts = expr_data)
  
  # Filter low expressed genes
  keep <- filterByExpr(dge, group = sample_info$condition)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Calculate normalization factors
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Apply voom transformation
  v <- voom(dge, design = model.matrix(~ 0 + sample_info$condition), plot = FALSE)
  
  # Store TMM+voom normalized data as an alternative
  tmm_voom_data <- v$E
  
  # Create comparison plots of normalization methods
  pdf("normalization_comparison.pdf", width = 12, height = 9)
  par(mfrow = c(2, 2))
  
  # Raw data
  boxplot(log2(expr_data + 1), main = "Raw counts (log2)", 
          xlab = "Samples", ylab = "log2(counts + 1)",
          col = rainbow(ncol(expr_data)), las = 2, cex.axis = 0.7)
  
  # VST normalized
  boxplot(normalized_data, main = "VST normalized", 
          xlab = "Samples", ylab = "VST values",
          col = rainbow(ncol(normalized_data)), las = 2, cex.axis = 0.7)
  
  # TMM+voom normalized
  boxplot(tmm_voom_data, main = "TMM+voom normalized", 
          xlab = "Samples", ylab = "log2 CPM",
          col = rainbow(ncol(tmm_voom_data)), las = 2, cex.axis = 0.7)
  
  # Density plots
  plot(density(log2(expr_data[, 1] + 1)), col = 1, main = "Density plots", 
       xlab = "Expression values", ylab = "Density")
  for (i in 2:min(10, ncol(expr_data))) {
    lines(density(log2(expr_data[, i] + 1)), col = i)
  }
  lines(density(normalized_data[, 1]), col = 1, lty = 2)
  for (i in 2:min(10, ncol(normalized_data))) {
    lines(density(normalized_data[, i]), col = i, lty = 2)
  }
  legend("topright", c("Raw (solid)", "Normalized (dashed)"), lty = 1:2)
  
  dev.off()
  
  # Use VST normalized data for WGCNA
  expr_data <- normalized_data
  
  # Save normalized data
  write.csv(normalized_data, file = "vst_normalized_data.csv")
  write.csv(tmm_voom_data, file = "tmm_voom_normalized_data.csv")
  
  cat("Normalization complete. Using variance stabilized data for WGCNA.\n")
} else {
  cat("Data does not appear to be raw counts. Applying standard scaling...\n")
  # Apply log2 transformation if the data is not already log-transformed
  if (max(expr_data) > 30) {  # Heuristic for detecting non-logged data
    cat("Applying log2 transformation...\n")
    expr_data <- log2(expr_data + 1)
  }
}

# Data preprocessing
# Transpose the expression matrix for WGCNA (samples as rows, genes as columns)
datExpr <- t(expr_data)

# Check for missing values
if (sum(is.na(datExpr)) > 0) {
  cat("Warning: Missing values detected in the expression data!\n")
}

# Check for genes with zero variance
var_genes <- apply(datExpr, 2, var, na.rm = TRUE)
zero_var_genes <- which(var_genes == 0)
if (length(zero_var_genes) > 0) {
  cat("Removing", length(zero_var_genes), "genes with zero variance\n")
  datExpr <- datExpr[, -zero_var_genes]
}

# Create trait data from metadata
# Ensure metadata is ordered the same as samples in datExpr
cat("Adjusting metadata to match expression data...\n")
print(paste("datExpr dimensions:", nrow(datExpr), "samples,", ncol(datExpr), "genes"))
print(paste("metadata dimensions:", nrow(metadata), "rows"))

# Check if sample order in datExpr rownames matches metadata sample_number
if (!identical(rownames(datExpr), as.character(metadata$sample_number))) {
  cat("WARNING: Sample order in expression data doesn't match metadata order.\n")
  cat("Reordering metadata to match expression data samples...\n")
  
  # Try to match metadata rows to datExpr row names
  sample_match <- match(rownames(datExpr), as.character(metadata$sample_number))
  
  # Check if all samples in datExpr were found in metadata
  if (any(is.na(sample_match))) {
    cat("ERROR: Some expression data samples not found in metadata!\n")
    missing_samples <- rownames(datExpr)[is.na(sample_match)]
    cat("Missing samples:", paste(missing_samples, collapse=", "), "\n")
    cat("This could be due to inconsistent naming between files.\n")
    
    # Emergency fix attempt - try to match by position if names don't match
    cat("Attempting emergency fix by position matching...\n")
    if (nrow(datExpr) == nrow(metadata)) {
      cat("Number of samples matches between expression data and metadata.\n")
      cat("Creating metadata mapping by position (assuming same order)...\n")
      
      # Create a mapping between row names and metadata
      mapping <- data.frame(
        expr_id = rownames(datExpr),
        metadata_id = as.character(metadata$sample_number),
        stringsAsFactors = FALSE
      )
      
      # Print the first few mappings to help debug
      cat("Sample mapping (first 10 rows):\n")
      print(head(mapping, 10))
      
      # Use metadata as is, assuming the order matches
      ordered_metadata <- metadata
    } else {
      # If the number of rows doesn't match, try to use the subset that does match
      cat("WARNING: Number of samples doesn't match between files!\n") 
      cat("Will use only the samples that exist in both datasets.\n")
      
      # Find samples that do match
      valid_samples <- rownames(datExpr)[!is.na(sample_match)]
      valid_indices <- sample_match[!is.na(sample_match)]
      
      # Subset metadata to only include matching samples
      ordered_metadata <- metadata[valid_indices, ]
      
      cat("Proceeding with", length(valid_samples), "samples that match between datasets.\n")
    }
  } else {
    # Normal case - we found all samples, just need to reorder
    ordered_metadata <- metadata[sample_match, ]
    cat("Successfully reordered metadata to match expression data.\n")
  }
} else {
  # Order already matches
  ordered_metadata <- metadata
  cat("Sample order in metadata and expression data already match.\n")
}

# Filter data to only include specified developmental stages
cat("Filtering data to only include target developmental stages: 2leaf, 4leaf, and bud...\n")

# Create a data frame to store all sample information together
sample_tracking <- data.frame(
  sample_id = rownames(datExpr),
  original_idx = 1:nrow(datExpr),
  genotype = ordered_metadata$genotype,
  stage = ordered_metadata$development_stage,
  stringsAsFactors = FALSE
)

# Print summary to check data alignment
cat("Pre-filtering sample summary:\n")
print(table(sample_tracking$genotype, sample_tracking$stage))

# Check which samples to keep (2leaf, 4leaf, and bud only)
stages_to_keep <- c("2leaf", "4leaf", "bud")
samples_to_keep <- sample_tracking$stage %in% stages_to_keep

# Report on filtering
cat("Found", sum(samples_to_keep), "samples in target stages out of", length(samples_to_keep), "total samples.\n")
if (sum(!samples_to_keep) > 0) {
  cat("Stages being removed:", 
      paste(unique(as.character(sample_tracking$stage[!samples_to_keep])), collapse=", "), "\n")
}

# Updated sample_tracking
sample_tracking_filtered <- sample_tracking[samples_to_keep, ]

# Verify filtered data has correct distribution
cat("Post-filtering sample summary:\n")
print(table(sample_tracking_filtered$genotype, sample_tracking_filtered$stage))

# Filter expression data
datExpr <- datExpr[samples_to_keep, ]

# Create fresh metadata that perfectly matches the filtered expression data
filtered_metadata <- data.frame(
  sample_number = as.numeric(rownames(datExpr)),
  sample_name = ordered_metadata$sample_name[samples_to_keep],
  Library_Pool = ordered_metadata$Library_Pool[samples_to_keep],
  Primer = ordered_metadata$Primer[samples_to_keep],
  genotype = ordered_metadata$genotype[samples_to_keep],
  development_stage = ordered_metadata$development_stage[samples_to_keep],
  stringsAsFactors = FALSE
)

# Verify the row count matches exactly
cat("After filtering:\n")
cat("  datExpr dimensions:", dim(datExpr)[1], "samples x", dim(datExpr)[2], "genes\n")
cat("  filtered_metadata dimensions:", dim(filtered_metadata)[1], "rows\n")

# Extra verification to ensure complete alignment
if (nrow(datExpr) != nrow(filtered_metadata)) {
  stop("CRITICAL ERROR: Sample count mismatch between expression data and metadata after filtering!")
}

# Double-check that we have all three stages in the filtered data
if (!all(stages_to_keep %in% unique(filtered_metadata$development_stage))) {
  cat("WARNING: Not all target stages present in filtered data!\n")
  cat("Expected stages:", paste(stages_to_keep, collapse=", "), "\n")
  cat("Actual stages:", paste(unique(filtered_metadata$development_stage), collapse=", "), "\n")
}

# Explicitly factor the variables to control levels
filtered_metadata$genotype <- factor(filtered_metadata$genotype, levels = c("Coast", "Inland"))
filtered_metadata$development_stage <- factor(filtered_metadata$development_stage, levels = stages_to_keep)

# Create trait variables in a controlled way
cat("Creating trait variables...\n")

# First verify levels are as expected
cat("Genotype levels:", paste(levels(filtered_metadata$genotype), collapse=", "), "\n")
cat("Stage levels:", paste(levels(filtered_metadata$development_stage), collapse=", "), "\n")

# Method 1: Direct model matrix approach
genotype_dummy <- model.matrix(~ 0 + genotype, data = filtered_metadata)
stage_dummy <- model.matrix(~ 0 + development_stage, data = filtered_metadata)

# Manually rename columns for clarity and consistency
colnames(genotype_dummy) <- levels(filtered_metadata$genotype)  # Should be c("Coast", "Inland")
colnames(stage_dummy) <- paste0("stage_", levels(filtered_metadata$development_stage))  # Using prefix to avoid confusion

# Confirm dimensions explicitly
cat("Dimensions check after filtering:\n")
cat("  datExpr:", dim(datExpr)[1], "rows (samples) x", dim(datExpr)[2], "columns (genes)\n")
cat("  genotype_dummy:", dim(genotype_dummy)[1], "rows x", dim(genotype_dummy)[2], "columns\n")
cat("  stage_dummy:", dim(stage_dummy)[1], "rows x", dim(stage_dummy)[2], "columns\n")

# Combine traits - with explicit checks
if (nrow(genotype_dummy) != nrow(stage_dummy)) {
  stop("CRITICAL ERROR: Row count mismatch between genotype and stage matrices!")
} 

# Combine the trait data
datTraits <- cbind(genotype_dummy, stage_dummy)
rownames(datTraits) <- rownames(datExpr)

# Final verification
cat("Final trait data dimensions:", dim(datTraits)[1], "rows x", dim(datTraits)[2], "columns\n")
cat("Columns in trait data:", paste(colnames(datTraits), collapse=", "), "\n")

cat("Trait data created successfully.\n")

cat("Trait data created successfully.\n")

# 2. Sample Clustering and Outlier Detection ---------------------------------

# Make sure we're using the filtered data for sample clustering
cat("Performing sample clustering...\n")

# Calculate sample distances
tryCatch({
  # Convert to matrix to ensure proper distance calculation
  datExpr_matrix <- as.matrix(datExpr)
  
  # Check for any issues with the data
  if (any(is.na(datExpr_matrix))) {
    cat("WARNING: NAs found in expression data. Removing NAs for distance calculation.\n")
    datExpr_matrix[is.na(datExpr_matrix)] <- 0
  }
  
  if (nrow(datExpr_matrix) < 2) {
    stop("Need at least 2 samples for clustering")
  }
  
  # Calculate distances - using robust correlation measure
  sampleDist <- dist(datExpr_matrix, method = "euclidean")
  
  # Check if distances are valid
  if (any(is.na(sampleDist))) {
    cat("WARNING: NAs in sample distances. Using alternative distance measure.\n")
    sampleDist <- dist(datExpr_matrix, method = "manhattan")
  }
  
  # Create the hierarchical clustering
  sampleTree <- flashClust::hclust(sampleDist, method = "average")
  
  # Plot the sample tree
  pdf("sample_clustering.pdf", width = 12, height = 9)
  par(cex = 0.6)
  par(mar = c(0, 4, 2, 0))
  
  # Use the filtered metadata for plot labels
  plot(sampleTree, main = "Sample Clustering", sub = "", xlab = "", 
       labels = filtered_metadata$sample_name,
       cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  
  # Add colored bands for genotype and developmental stage
  plotDendroAndColors(sampleTree, colors = data.frame(
    Genotype = labels2colors(factor(filtered_metadata$genotype)),
    Stage = labels2colors(factor(filtered_metadata$development_stage))
  ), groupLabels = c("Genotype", "Stage"), 
  main = "Sample dendrogram and trait heatmap")
  
  dev.off()
  
  cat("Sample clustering completed successfully.\n")
}, error = function(e) {
  cat("ERROR in sample clustering:", e$message, "\n")
  cat("Skipping dendrogram creation and continuing with analysis.\n")
})

# Optional: remove outliers if needed
# Define a height cut for outlier detection
# heightCut = 100  # Adjust as needed
# clusterCut = cutreeStatic(sampleTree, cutHeight = heightCut, minSize = 10)
# keepSamples = (clusterCut == 1)
# datExpr = datExpr[keepSamples, ]
# datTraits = datTraits[keepSamples, ]

# 3. Choose Soft-Thresholding Power ------------------------------------------

# Choose a set of soft-thresholding powers
powers <- c(1:30)

# Error handling for power selection
cat("Selecting optimal soft-thresholding power...\n")

tryCatch({
  # First convert datExpr to a numeric matrix if not already
  datExpr_matrix <- as.matrix(datExpr)
  
  # Check for and handle NA values
  if(any(is.na(datExpr_matrix))) {
    cat("WARNING: NAs found in expression data. Imputing with row means for power selection.\n")
    for(i in 1:nrow(datExpr_matrix)) {
      row_mean <- mean(datExpr_matrix[i,], na.rm=TRUE)
      na_idx <- is.na(datExpr_matrix[i,])
      datExpr_matrix[i, na_idx] <- row_mean
    }
  }
  
  # Check for zero variance genes
  var_genes <- apply(datExpr_matrix, 2, var, na.rm = TRUE)
  zero_var_genes <- which(var_genes == 0 | is.na(var_genes))
  if(length(zero_var_genes) > 0) {
    cat("Removing", length(zero_var_genes), "genes with zero variance for power selection.\n")
    datExpr_matrix <- datExpr_matrix[, -zero_var_genes]
  }
  
  cat("Running pickSoftThreshold with", ncol(datExpr_matrix), "genes and", nrow(datExpr_matrix), "samples...\n")
  
  # Try with blockwise approach to reduce memory usage
  cat("Using block-wise approach to reduce memory usage...\n")
  max_block_size <- min(ncol(datExpr_matrix), 5000)
  
  # Call the network topology analysis function with reduced blocksize
  sft <- pickSoftThreshold(
    datExpr_matrix,
    powerVector = powers,
    verbose = 5,
    networkType = "signed",
    corFnc = "bicor", # Use biweight correlation for robustness
    corOptions = list(use = "pairwise.complete.obs", maxPOutliers = 0.1),
    blockSize = max_block_size
  )
  
  # Plot the results
  pdf("soft_thresholding_power.pdf", width = 12, height = 9)
  par(mfrow = c(1, 2))
  cex1 = 0.9
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", 
       ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = paste("Scale independence"))
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = cex1, col = "red")
  abline(h = 0.80, col = "red") # This line corresponds to R^2 = 0.80
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = paste("Mean connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, 
       cex = cex1, col = "red")
  dev.off()
  
  # Select the power based on the plot
  # Find the power where we exceed R^2 cutoff (0.8)
  r2_cutoff <- 0.80
  r2_values <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  best_power_idx <- which(r2_values >= r2_cutoff)[1]
  
  if(!is.na(best_power_idx)) {
    power <- sft$fitIndices[best_power_idx, 1]
    cat("Automatically selected power =", power, "with R^2 =", r2_values[best_power_idx], "\n")
  } else {
    # If no power meets the cutoff, choose based on curve pattern
    # Look for the power where the improvement in R^2 slows down
    r2_diff <- c(0, diff(r2_values))
    elbow_idx <- which(r2_diff < 0.05 & r2_values > 0.1)[1]
    
    if(!is.na(elbow_idx)) {
      power <- sft$fitIndices[elbow_idx, 1]
      cat("Selected power by curve pattern =", power, "with R^2 =", r2_values[elbow_idx], "\n")
    } else {
      power <- 12  # Default value if automatic selection fails
      cat("Could not automatically determine power, using default power of", power, "\n")
    }
  }
}, error = function(e) {
  cat("ERROR in pickSoftThreshold:", e$message, "\n")
  cat("This can happen with small sample sizes or data quality issues.\n")
  
  # Set a default power if the function fails
  power <<- 12
  cat("Using default power of 12 and continuing with analysis.\n")
})

# Make sure power is defined in the global environment
if(!exists("power") || is.na(power)) {
  power <- 12
  cat("Using default power of 12 for network construction.\n")
}

cat("Final selected power for network construction:", power, "\n")

# 4. Construct the Gene Co-expression Network -------------------------------

# Calculate adjacency matrix
adjacency <- adjacency(datExpr, power = power, type = "signed")

# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# 5. Identify Gene Modules --------------------------------------------------

# Hierarchical clustering of genes
geneTree <- flashClust::hclust(as.dist(dissTOM), method = "average")

# Plot the gene clustering tree
pdf("gene_clustering.pdf", width = 12, height = 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# Module identification using dynamic tree cut
# Set the minimum module size
minModuleSize <- 30

# Module identification
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric labels to colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# Plot the module assignment under the gene tree
pdf("module_identification.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                   dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05,
                   main = "Gene dendrogram and module colors")
dev.off()

# 6. Merge Similar Modules -------------------------------------------------

# Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot the result
pdf("module_eigengene_clustering.pdf", width = 12, height = 9)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Choose a height cut of 0.25 (can be adjusted)
MEDissThres <- 0.25
abline(h = MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors <- merge$colors
# Eigengenes of the new merged modules
mergedMEs <- merge$newMEs

# Plot the comparison between original and merged modules
pdf("merged_modules.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                   c("Dynamic Tree Cut", "Merged dynamic"),
                   dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05)
dev.off()

# Replace original module assignments with merged ones
moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs

# 7. Correlate Modules with Traits -----------------------------------------

# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Correlate module eigengenes with traits
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Create a color map for the correlation values
colorMap <- function(corr) {
  colorRampPalette(c("blue", "white", "red"))(100)[as.numeric(cut(corr, breaks = 100))]
}

# Plot the module-trait relationships
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

pdf("module_trait_relationships.pdf", width = 12, height = 12)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
              xLabels = colnames(datTraits),
              yLabels = names(MEs),
              ySymbols = names(MEs),
              colorLabels = FALSE,
              colors = blueWhiteRed(50),
              textMatrix = textMatrix,
              setStdMargins = FALSE,
              cex.text = 0.5,
              zlim = c(-1, 1),
              main = "Module-trait relationships")
dev.off()

# 8. Identify Genes with High Module Membership and Gene Significance -------

# Define the traits of interest
cat("Creating trait correlations for gene significance calculation...\n")

# Verify column names in datTraits
cat("Available traits:", paste(colnames(datTraits), collapse=", "), "\n")

# Genotypes
trait_genotype_Coast <- as.data.frame(datTraits[, "Coast", drop=FALSE])
# No need to rename as column name is already correct

trait_genotype_Inland <- as.data.frame(datTraits[, "Inland", drop=FALSE])
# No need to rename as column name is already correct 

# Developmental stages - using the colnames as they appear in datTraits
stage_cols <- grep("^stage_", colnames(datTraits), value=TRUE)
cat("Stage columns found:", paste(stage_cols, collapse=", "), "\n")

# Create trait dataframes for each stage
trait_stage_2leaf <- as.data.frame(datTraits[, grep("stage_2leaf", colnames(datTraits)), drop=FALSE])
trait_stage_4leaf <- as.data.frame(datTraits[, grep("stage_4leaf", colnames(datTraits)), drop=FALSE])
trait_stage_bud <- as.data.frame(datTraits[, grep("stage_bud", colnames(datTraits)), drop=FALSE])

# Verify the trait dataframes
cat("Trait dataframe dimensions:\n")
cat("  Coast:", dim(trait_genotype_Coast)[1], "rows x", dim(trait_genotype_Coast)[2], "columns\n")
cat("  Inland:", dim(trait_genotype_Inland)[1], "rows x", dim(trait_genotype_Inland)[2], "columns\n")
cat("  2leaf:", dim(trait_stage_2leaf)[1], "rows x", dim(trait_stage_2leaf)[2], "columns\n")
cat("  4leaf:", dim(trait_stage_4leaf)[1], "rows x", dim(trait_stage_4leaf)[2], "columns\n")
cat("  bud:", dim(trait_stage_bud)[1], "rows x", dim(trait_stage_bud)[2], "columns\n")

# Names of all the modules
modNames <- substring(names(MEs), 3)

# Gene significance and module membership for genotype
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

# Gene significance calculation
cat("Calculating gene significance values...\n")

# Calculate number of samples for correlation p-value calculation
nSamples <- nrow(datExpr)
cat("Number of samples for correlation:", nSamples, "\n")

# First calculate the correlation between each gene and each module eigengene
# Ensure the data is properly formatted
if (!is.numeric(as.matrix(datExpr))) {
  cat("WARNING: Converting expression data to numeric for correlation calculations.\n")
  datExpr_matrix <- apply(datExpr, 2, as.numeric)
} else {
  datExpr_matrix <- as.matrix(datExpr)
}

# Gene significance calculations for genotypes using consistent naming
cat("Calculating gene significance for genotypes...\n")
GS_Coast <- as.data.frame(WGCNA::cor(datExpr_matrix, trait_genotype_Coast, use = "p"))
names(GS_Coast) <- "GS_Coast"
GSPvalue_Coast <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(GS_Coast), nSamples))
names(GSPvalue_Coast) <- "p.GS_Coast"

GS_Inland <- as.data.frame(WGCNA::cor(datExpr_matrix, trait_genotype_Inland, use = "p"))
names(GS_Inland) <- "GS_Inland"  
GSPvalue_Inland <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(GS_Inland), nSamples))
names(GSPvalue_Inland) <- "p.GS_Inland"

# Gene significance calculations for developmental stages using consistent naming
cat("Calculating gene significance for developmental stages...\n")
GS_2leaf <- as.data.frame(WGCNA::cor(datExpr_matrix, trait_stage_2leaf, use = "p"))
names(GS_2leaf) <- "GS_2leaf"
GSPvalue_2leaf <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(GS_2leaf), nSamples))
names(GSPvalue_2leaf) <- "p.GS_2leaf"

GS_4leaf <- as.data.frame(WGCNA::cor(datExpr_matrix, trait_stage_4leaf, use = "p"))
names(GS_4leaf) <- "GS_4leaf"
GSPvalue_4leaf <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(GS_4leaf), nSamples))
names(GSPvalue_4leaf) <- "p.GS_4leaf"

GS_bud <- as.data.frame(WGCNA::cor(datExpr_matrix, trait_stage_bud, use = "p"))
names(GS_bud) <- "GS_bud"
GSPvalue_bud <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(GS_bud), nSamples))
names(GSPvalue_bud) <- "p.GS_bud"

cat("Gene significance calculations completed successfully.\n")

# 9. Visualize Module-Trait Relationships -----------------------------------

# Choose a module of interest from heatmap (example: blue module)
module = "blue"  # Change this to your module of interest

# Get genes in the module
module_genes <- colnames(datExpr)[moduleColors == module]
module_gene_indices <- match(module_genes, colnames(datExpr))

# Plot module membership vs. gene significance for genotype
pdf(paste0("MM_vs_GS_", module, "_Coast.pdf"), width = 9, height = 9)
column = match(module, modNames)
moduleGenes = moduleColors == module
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance_Coast[moduleGenes, 1],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Coast genotype",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# Plot module membership vs. gene significance for Inland genotype
pdf(paste0("MM_vs_GS_", module, "_Inland.pdf"), width = 9, height = 9)
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance_Inland[moduleGenes, 1],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Inland genotype",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# Plot for developmental stages
pdf(paste0("MM_vs_GS_", module, "_2leaf.pdf"), width = 9, height = 9)
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance_2leaf[moduleGenes, 1],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for 2-leaf stage",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

pdf(paste0("MM_vs_GS_", module, "_4leaf.pdf"), width = 9, height = 9)
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance_4leaf[moduleGenes, 1],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for 4-leaf stage",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

pdf(paste0("MM_vs_GS_", module, "_bud.pdf"), width = 9, height = 9)
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance_bud[moduleGenes, 1],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for bud stage",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# 10. Export Network Data and Results ---------------------------------------

# Export gene lists for each module
for (module in unique(moduleColors)) {
  if (module != "grey") {  # Skip unassigned genes
    module_genes <- colnames(datExpr)[moduleColors == module]
    write.table(module_genes, file = paste0("module_", module, "_genes.txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

# Create a table with gene module assignments
genes <- colnames(datExpr)
gene_module_table <- data.frame(GeneID = genes, Module = moduleColors)
write.csv(gene_module_table, file = "gene_module_assignments.csv", row.names = FALSE)

# Export module eigengene values
write.csv(MEs, file = "module_eigengenes.csv")

# Export module-trait correlations
module_trait_corr_table <- data.frame(
  Module = rownames(moduleTraitCor),
  moduleTraitCor,
  stringsAsFactors = FALSE
)
write.csv(module_trait_corr_table, file = "module_trait_correlations.csv", row.names = FALSE)

# Export gene significance and module membership values
for (module in unique(moduleColors)) {
  if (module != "grey") {
    module_column <- paste("MM", module, sep = "")
    if (module_column %in% colnames(geneModuleMembership)) {
      module_genes <- colnames(datExpr)[moduleColors == module]
      module_idx <- match(module_genes, colnames(datExpr))
      
      gene_info <- data.frame(
        GeneID = module_genes,
        ModuleMembership = geneModuleMembership[module_idx, module_column],
        MMPvalue = MMPvalue[module_idx, paste("p.MM", module, sep = "")],
        GS_Coast = geneTraitSignificance_Coast[module_idx, 1],
        GSPvalue_Coast = GSPvalue_Coast[module_idx, 1],
        GS_Inland = geneTraitSignificance_Inland[module_idx, 1],
        GSPvalue_Inland = GSPvalue_Inland[module_idx, 1],
        GS_leaf2 = geneTraitSignificance_2leaf[module_idx, 1],
        GSPvalue_leaf2 = GSPvalue_2leaf[module_idx, 1],
        GS_leaf4 = geneTraitSignificance_4leaf[module_idx, 1],
        GSPvalue_leaf4 = GSPvalue_4leaf[module_idx, 1],
        GS_bud = geneTraitSignificance_bud[module_idx, 1],
        GSPvalue_bud = GSPvalue_bud[module_idx, 1]
      )
      
      write.csv(gene_info, file = paste0("module_", module, "_gene_info.csv"), row.names = FALSE)
    }
  }
}

# Export a network for visualization in other tools (e.g., Cytoscape)
# Select a module to export (can be changed to export multiple modules)
module_for_export <- "blue"  # Change to your module of interest
genes_in_module <- colnames(datExpr)[moduleColors == module_for_export]
module_TOM <- TOM[moduleColors == module_for_export, moduleColors == module_for_export]
dimnames(module_TOM) <- list(genes_in_module, genes_in_module)

# Export the network
edge_list <- data.frame(
  Source = rep(rownames(module_TOM), each = ncol(module_TOM)),
  Target = rep(colnames(module_TOM), nrow(module_TOM)),
  Weight = as.vector(module_TOM)
)
# Remove self-loops and keep only stronger connections
edge_list <- edge_list[edge_list$Source != edge_list$Target & edge_list$Weight > 0.1, ]
write.csv(edge_list, file = paste0("module_", module_for_export, "_network.csv"), row.names = FALSE)

# 11. Additional Analysis: Hub Gene Identification --------------------------

# Identify hub genes (genes with highest module membership)
hub_genes <- list()
for (module in unique(moduleColors)) {
  if (module != "grey") {
    module_column <- which(modNames == module)
    if (length(module_column) > 0) {
      module_genes <- colnames(datExpr)[moduleColors == module]
      module_membership <- geneModuleMembership[moduleColors == module, module_column]
      # Sort genes by module membership
      sorted_idx <- order(module_membership, decreasing = TRUE)
      # Select top 10 hub genes
      hub_genes[[module]] <- module_genes[sorted_idx[1:min(10, length(sorted_idx))]]
    }
  }
}

# Export hub genes
hub_gene_table <- do.call(rbind, lapply(names(hub_genes), function(module) {
  data.frame(Module = module, HubGene = hub_genes[[module]], stringsAsFactors = FALSE)
}))
write.csv(hub_gene_table, file = "hub_genes.csv", row.names = FALSE)

# 12. Create Final Report --------------------------------------------------

# Emergency fix for module-trait relationship data
cat("Checking and fixing module-trait data consistency before creating report...\n")

tryCatch({
  # Check if geneTraitSignificance variables are defined with inconsistent names
  gs_vars_exist <- exists("geneTraitSignificance_Coast") || exists("GS_Coast")
  if (gs_vars_exist) {
    cat("Gene significance variables found, checking consistency...\n")
    
    # Standardize names - create consistent GS_ variables regardless of original naming
    if (exists("geneTraitSignificance_Coast") && !exists("GS_Coast")) {
      GS_Coast <- geneTraitSignificance_Coast
      names(GS_Coast) <- "GS_Coast"
      GSPvalue_Coast <- GSPvalue_Coast
      names(GSPvalue_Coast) <- "p.GS_Coast"
      cat("Standardized Coast genotype variables.\n")
    }
    
    if (exists("geneTraitSignificance_Inland") && !exists("GS_Inland")) {
      GS_Inland <- geneTraitSignificance_Inland
      names(GS_Inland) <- "GS_Inland"
      GSPvalue_Inland <- GSPvalue_Inland
      names(GSPvalue_Inland) <- "p.GS_Inland"
      cat("Standardized Inland genotype variables.\n")
    }
    
    if (exists("geneTraitSignificance_2leaf") && !exists("GS_2leaf")) {
      GS_2leaf <- geneTraitSignificance_2leaf
      names(GS_2leaf) <- "GS_2leaf"
      GSPvalue_2leaf <- GSPvalue_2leaf
      names(GSPvalue_2leaf) <- "p.GS_2leaf"
      cat("Standardized 2leaf stage variables.\n")
    }
    
    if (exists("geneTraitSignificance_4leaf") && !exists("GS_4leaf")) {
      GS_4leaf <- geneTraitSignificance_4leaf
      names(GS_4leaf) <- "GS_4leaf"
      GSPvalue_4leaf <- GSPvalue_4leaf
      names(GSPvalue_4leaf) <- "p.GS_4leaf"
      cat("Standardized 4leaf stage variables.\n")
    }
    
    if (exists("geneTraitSignificance_bud") && !exists("GS_bud")) {
      GS_bud <- geneTraitSignificance_bud
      names(GS_bud) <- "GS_bud"
      GSPvalue_bud <- GSPvalue_bud
      names(GSPvalue_bud) <- "p.GS_bud"
      cat("Standardized bud stage variables.\n")
    }
  } else {
    cat("No gene significance variables found. Creating them from scratch...\n")
    
    # Check if datTraits exists and has the expected columns
    if (!exists("datTraits")) {
      cat("WARNING: datTraits not found. Creating minimal trait data for report...\n")
      
      # Create basic trait data from metadata if available
      if (exists("filtered_metadata") || exists("metadata")) {
        use_metadata <- if (exists("filtered_metadata")) filtered_metadata else metadata
        
        # Create traits directly from metadata
        genotype_dummy <- model.matrix(~ 0 + genotype, data = use_metadata)
        colnames(genotype_dummy) <- c("Coast", "Inland")
        
        stage_dummy <- model.matrix(~ 0 + development_stage, data = use_metadata) 
        stage_cols <- c("stage_2leaf", "stage_4leaf", "stage_bud")
        if (ncol(stage_dummy) == length(stage_cols)) {
          colnames(stage_dummy) <- stage_cols
        } else {
          cat("WARNING: Unexpected number of stage columns:", ncol(stage_dummy), "\n")
          colnames(stage_dummy) <- paste0("stage_", 1:ncol(stage_dummy))
        }
        
        # Combine traits
        datTraits <- cbind(genotype_dummy, stage_dummy)
        rownames(datTraits) <- rownames(datExpr)
        
        cat("Created datTraits with columns:", paste(colnames(datTraits), collapse=", "), "\n")
      } else {
        cat("ERROR: No metadata available to create traits. Report may be incomplete.\n")
      }
    }
    
    # Create gene significance values from datTraits
    if (exists("datTraits")) {
      # Ensure datExpr is numeric
      if (!is.numeric(as.matrix(datExpr))) {
        cat("Converting expression data to numeric...\n")
        datExpr_matrix <- apply(datExpr, 2, as.numeric)
      } else {
        datExpr_matrix <- as.matrix(datExpr)
      }
      
      # Calculate number of samples
      nSamples <- nrow(datExpr)
      
      # Check if the expected columns exist in datTraits
      genotype_cols <- c("Coast", "Inland")
      if (all(genotype_cols %in% colnames(datTraits))) {
        cat("Creating gene significance for genotypes...\n")
        # Create trait dataframes
        trait_genotype_Coast <- as.data.frame(datTraits[, "Coast", drop=FALSE])
        trait_genotype_Inland <- as.data.frame(datTraits[, "Inland", drop=FALSE])
        
        # Calculate gene significance
        GS_Coast <- as.data.frame(cor(datExpr_matrix, trait_genotype_Coast, use = "p"))
        names(GS_Coast) <- "GS_Coast"
        GSPvalue_Coast <- as.data.frame(corPvalueStudent(as.matrix(GS_Coast), nSamples))
        names(GSPvalue_Coast) <- "p.GS_Coast"
        
        GS_Inland <- as.data.frame(cor(datExpr_matrix, trait_genotype_Inland, use = "p"))
        names(GS_Inland) <- "GS_Inland"
        GSPvalue_Inland <- as.data.frame(corPvalueStudent(as.matrix(GS_Inland), nSamples))
        names(GSPvalue_Inland) <- "p.GS_Inland"
      } else {
        cat("WARNING: Genotype columns not found in datTraits.\n")
      }
      
      # Look for stage columns
      stage_cols <- grep("stage_|leaf|bud", colnames(datTraits), value=TRUE)
      if (length(stage_cols) >= 3) {
        cat("Found stage columns:", paste(stage_cols, collapse=", "), "\n")
        
        # Identify which columns correspond to each stage
        leaf2_col <- grep("2leaf|leaf2|stage_2", stage_cols, value=TRUE)[1]
        leaf4_col <- grep("4leaf|leaf4|stage_4", stage_cols, value=TRUE)[1]
        bud_col <- grep("bud", stage_cols, value=TRUE)[1]
        
        if (!is.na(leaf2_col)) {
          trait_stage_2leaf <- as.data.frame(datTraits[, leaf2_col, drop=FALSE])
          GS_2leaf <- as.data.frame(cor(datExpr_matrix, trait_stage_2leaf, use = "p"))
          names(GS_2leaf) <- "GS_2leaf"
          GSPvalue_2leaf <- as.data.frame(corPvalueStudent(as.matrix(GS_2leaf), nSamples))
          names(GSPvalue_2leaf) <- "p.GS_2leaf"
        }
        
        if (!is.na(leaf4_col)) {
          trait_stage_4leaf <- as.data.frame(datTraits[, leaf4_col, drop=FALSE])
          GS_4leaf <- as.data.frame(cor(datExpr_matrix, trait_stage_4leaf, use = "p"))
          names(GS_4leaf) <- "GS_4leaf"
          GSPvalue_4leaf <- as.data.frame(corPvalueStudent(as.matrix(GS_4leaf), nSamples))
          names(GSPvalue_4leaf) <- "p.GS_4leaf"
        }
        
        if (!is.na(bud_col)) {
          trait_stage_bud <- as.data.frame(datTraits[, bud_col, drop=FALSE])
          GS_bud <- as.data.frame(cor(datExpr_matrix, trait_stage_bud, use = "p"))
          names(GS_bud) <- "GS_bud"
          GSPvalue_bud <- as.data.frame(corPvalueStudent(as.matrix(GS_bud), nSamples))
          names(GSPvalue_bud) <- "p.GS_bud"
        }
        
        cat("Successfully created gene significance variables for stages.\n")
      } else {
        cat("WARNING: Not enough stage columns found in datTraits.\n")
      }
    }
  }
  
  cat("Module-trait data consistency check completed.\n")
}, error = function(e) {
  cat("WARNING: Error in module-trait data consistency check:", e$message, "\n")
  cat("Some visualizations may be incomplete.\n")
})

# Open PDF file for report
pdf("WGCNA_summary_report.pdf", width = 12, height = 12)

# Memory-efficient MDS plot of modules
cat("Creating MDS visualization with memory constraints...\n")

tryCatch({
  # Instead of using the full TOM matrix, we'll use a subset of genes or work with modules
  
  # Option 1: Use module eigengenes for MDS (very memory efficient)
  cat("Using module eigengenes for MDS plot (memory-efficient approach)...\n")
  
  # MDS on module eigengenes - this uses minimal memory
  ME_dist <- dist(t(MEs))
  ME_mds <- cmdscale(ME_dist, 2)
  
  # Plot MDS of module eigengenes
  plot(ME_mds, type = "n", main = "MDS plot of module eigengenes",
       xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
  text(ME_mds[, 1], ME_mds[, 2], labels = substring(rownames(ME_mds), 3),
       col = substring(rownames(ME_mds), 3), cex = 1.5)
  
  # Option 2: If still needed, try MDS on a sample of genes from each module
  cat("Creating gene-level MDS plot with sampled genes to avoid memory limits...\n")
  
  # Create a random sample of genes, with higher sampling from larger modules
  set.seed(12345) # For reproducibility
  
  # Determine number of genes to sample in total (adjust as needed)
  n_genes_to_sample <- min(5000, ncol(datExpr))
  
  # Initialize a vector to store selected gene indices
  selected_genes <- c()
  
  # Sample genes from each module proportionally to module size
  for (color in unique(moduleColors)) {
    # Get genes in this module
    genes_in_module <- which(moduleColors == color)
    
    # Determine sample size for this module (proportional to module size)
    n_sample <- max(5, min(length(genes_in_module), 
                        round(length(genes_in_module) / ncol(datExpr) * n_genes_to_sample)))
    
    # Sample genes from this module
    if (length(genes_in_module) > n_sample) {
      sampled_genes <- sample(genes_in_module, n_sample)
    } else {
      sampled_genes <- genes_in_module
    }
    
    # Add to selected genes
    selected_genes <- c(selected_genes, sampled_genes)
  }
  
  # Extract TOM subset for selected genes
  dissTOM_subset <- dissTOM[selected_genes, selected_genes]
  
  # Calculate MDS on the subset
  cmd_subset <- cmdscale(as.dist(dissTOM_subset), 2)
  
  # Plot MDS of the sampled genes
  plot(cmd_subset, 
       xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2",
       main = "MDS plot of sampled genes colored by module",
       col = moduleColors[selected_genes], cex = 0.5, pch = 19)
  
  # Add a legend for module colors
  legend("topright", legend = unique(moduleColors), 
         col = unique(moduleColors), pch = 19, cex = 0.6)
  
}, error = function(e) {
  cat("WARNING: MDS visualization failed due to memory constraints:", e$message, "\n")
  cat("Skipping MDS plot and continuing with other visualizations.\n")
  
  # Create a simple text notice in the PDF
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
  text(0.5, 0.5, "MDS plot skipped due to memory constraints.\nConsider running on a higher memory system\nor reducing the number of genes in the analysis.",
       cex = 1.2)
})

# Heatmap of the expression of genes in the most significant modules
# For demonstration, selecting the blue module (change as needed)
significant_module <- "blue"  # Change to your most significant module
genes_in_sig_module <- colnames(datExpr)[moduleColors == significant_module]

if (length(genes_in_sig_module) > 100) {
  # If too many genes, select top 100 by module membership
  module_column <- which(modNames == significant_module)
  module_membership <- geneModuleMembership[moduleColors == significant_module, module_column]
  sorted_idx <- order(module_membership, decreasing = TRUE)
  genes_to_plot <- genes_in_sig_module[sorted_idx[1:100]]
} else {
  genes_to_plot <- genes_in_sig_module
}

# Get expression data for these genes
expr_subset <- t(datExpr[, match(genes_to_plot, colnames(datExpr))])

# Scale the expression data
expr_scaled <- t(scale(t(expr_subset)))

# Create annotation for the heatmap
cat("Creating sample annotation for heatmap...\n")

tryCatch({
  # First check if filtered_metadata exists
  if (exists("filtered_metadata")) {
    cat("Using filtered_metadata for sample annotations...\n")
    sample_annotation <- data.frame(
      Genotype = filtered_metadata$genotype,
      Stage = filtered_metadata$development_stage,
      row.names = rownames(datExpr)
    )
  } else if (exists("datTraits")) {
    # Alternative approach: extract annotations from datTraits
    cat("Extracting sample annotations from datTraits...\n")
    
    # For genotype, determine which column has value 1 for each sample
    genotype_cols <- grep("^Coast$|^Inland$", colnames(datTraits), value=TRUE)
    if (length(genotype_cols) >= 2) {
      genotype_values <- apply(datTraits[, genotype_cols], 1, function(x) {
        genotype_cols[which.max(x)]
      })
    } else {
      genotype_values <- rep("Unknown", nrow(datTraits))
      cat("WARNING: Could not determine genotype from datTraits.\n")
    }
    
    # For stage, look for columns containing stage_ or leaf or bud
    stage_cols <- grep("stage_|leaf|bud", colnames(datTraits), value=TRUE)
    if (length(stage_cols) >= 3) {
      stage_values <- apply(datTraits[, stage_cols], 1, function(x) {
        stage_name <- stage_cols[which.max(x)]
        # Clean up stage name for display
        stage_name <- gsub("stage_", "", stage_name)
        stage_name <- gsub("^leaf", "", stage_name)
        if (grepl("^\\d", stage_name)) {
          stage_name <- paste0(stage_name, "leaf")
        }
        return(stage_name)
      })
    } else {
      stage_values <- rep("Unknown", nrow(datTraits))
      cat("WARNING: Could not determine developmental stage from datTraits.\n")
    }
    
    # Create the annotation dataframe
    sample_annotation <- data.frame(
      Genotype = genotype_values,
      Stage = stage_values,
      row.names = rownames(datExpr)
    )
  } else {
    # Last resort: create minimal annotation
    cat("WARNING: Neither filtered_metadata nor datTraits found.\n")
    cat("Creating minimal sample annotation...\n")
    
    sample_annotation <- data.frame(
      Genotype = rep(c("Coast", "Inland"), each = nrow(datExpr)/2),
      Stage = rep(rep(c("2leaf", "4leaf", "bud"), each = nrow(datExpr)/6), 2),
      row.names = rownames(datExpr)
    )
  }
  
  # Confirm annotation dimensions
  cat("Sample annotation dimensions:", dim(sample_annotation)[1], "rows x", 
      dim(sample_annotation)[2], "columns\n")
  
  # Double-check that we have the right number of rows
  if (nrow(sample_annotation) != nrow(datExpr)) {
    cat("ERROR: Sample annotation row count", nrow(sample_annotation), 
        "doesn't match expression data row count", nrow(datExpr), "\n")
    
    # Emergency fix: create new annotation matching exactly
    cat("Applying emergency fix for annotation...\n")
    sample_annotation <- data.frame(
      Genotype = factor(rep(c("Coast", "Inland"), length.out = nrow(datExpr))),
      Stage = factor(rep(c("2leaf", "4leaf", "bud"), length.out = nrow(datExpr))),
      row.names = rownames(datExpr)
    )
  }
  
  # Print sample counts by group for verification
  cat("Sample counts by group:\n")
  print(table(sample_annotation$Genotype, sample_annotation$Stage))
  
  # Create the heatmap
  pheatmap(expr_scaled, 
           annotation_col = sample_annotation,
           show_rownames = FALSE, 
           cluster_rows = TRUE, 
           cluster_cols = TRUE,
           main = paste("Expression of genes in", significant_module, "module"))
}, error = function(e) {
  cat("ERROR in heatmap creation:", e$message, "\n")
  
  # Simple fallback visualization
  heatmap(expr_scaled, 
          main = paste("Expression of genes in", significant_module, "module"),
          scale = "none", 
          col = colorRampPalette(c("blue", "white", "red"))(100))
  
  cat("Created fallback heatmap without annotations.\n")
})

# Close the PDF file
dev.off()

# 13. Differential network analysis between genotypes -----------------------

# Optional: Perform differential network analysis between genotypes
# Split the data by genotype
coast_samples <- metadata$genotype == "Coast"
inland_samples <- metadata$genotype == "Inland"

# Create separate expression matrices for each genotype
datExpr_coast <- datExpr[coast_samples, ]
datExpr_inland <- datExpr[inland_samples, ]

# Run module preservation analysis
# This checks if modules identified in one genotype are preserved in the other
multiExpr <- list(Coast = list(data = datExpr_coast), Inland = list(data = datExpr_inland))

# Calculate module preservation statistics
mp <- modulePreservation(multiExpr, 
                         networkType = "signed", 
                         referenceNetworks = c(1), 
                         nPermutations = 100, # Increase for publication-quality analysis
                         randomSeed = 12345, 
                         verbose = 3)

# Create a table with the preservation statistics
ref <- 1 # Reference network (Coast)
test <- 2 # Test network (Inland)

# Extract preservation statistics
# Z-summary is a measure of preservation
preservation_stats <- data.frame(
  Module = modNames,
  moduleSize = mp$quality$moduleSize[[ref]],
  Zsummary = mp$preservation$Z$ref.vs.test[[ref]][[test]]$Zsummary$observed,
  medianRank = mp$preservation$rankPres$ref.vs.test[[ref]][[test]]$medianRank$observed
)

# Sort the modules by Z-summary
preservation_stats <- preservation_stats[order(preservation_stats$Zsummary, decreasing = TRUE), ]

# Write the preservation statistics to a file
write.csv(preservation_stats, file = "module_preservation_statistics.csv", row.names = FALSE)

# Plot the preservation statistics
pdf("module_preservation_plot.pdf", width = 10, height = 8)
par(mfrow = c(1, 2))
# Z-summary plot
x <- mp$preservation$Z$ref.vs.test[[ref]][[test]]$Zsummary$observed
y <- mp$preservation$Z$ref.vs.test[[ref]][[test]]$moduleSize
plot(y, x, col = moduleColors, 
     main = "Module preservation - Zsummary",
     xlab = "Module size", ylab = "Zsummary",
     pch = 19, cex = 1.5)
abline(h = 2, col = "blue", lty = 2)
abline(h = 10, col = "red", lty = 2)
text(y, x, labels = modNames, pos = 1, cex = 0.7)
# Median rank plot
y <- mp$preservation$rankPres$ref.vs.test[[ref]][[test]]$medianRank$observed
plot(mp$preservation$rankPres$ref.vs.test[[ref]][[test]]$moduleSize, y,
     col = moduleColors, 
     main = "Module preservation - Median rank",
     xlab = "Module size", ylab = "Median rank",
     pch = 19, cex = 1.5)
text(mp$preservation$rankPres$ref.vs.test[[ref]][[test]]$moduleSize, y,
     labels = modNames, pos = 1, cex = 0.7)
dev.off()

# 14. Stage-specific analysis ---------------------------------------------

# Create combinations of genotype and stage for more detailed analysis
# This allows us to find modules associated with specific combinations
metadata$genotype_stage <- paste(metadata$genotype, metadata$development_stage, sep = "_")

# Create dummy variables for these combinations
genotype_stage_dummy <- model.matrix(~ 0 + factor(metadata$genotype_stage))
colnames(genotype_stage_dummy) <- levels(factor(metadata$genotype_stage))

# Add these to the trait data
datTraits_expanded <- cbind(datTraits, genotype_stage_dummy)

# Correlate modules with these expanded traits
moduleTraitCor_expanded <- cor(MEs, datTraits_expanded, use = "p")
moduleTraitPvalue_expanded <- corPvalueStudent(moduleTraitCor_expanded, nSamples)

# Create a heatmap for these expanded correlations
textMatrix_expanded <- paste(signif(moduleTraitCor_expanded, 2), "\n(", 
                            signif(moduleTraitPvalue_expanded, 1), ")", sep = "")
dim(textMatrix_expanded) <- dim(moduleTraitCor_expanded)

pdf("module_trait_relationships_expanded.pdf", width = 14, height = 12)
labeledHeatmap(Matrix = moduleTraitCor_expanded,
              xLabels = colnames(datTraits_expanded),
              yLabels = names(MEs),
              ySymbols = names(MEs),
              colorLabels = FALSE,
              colors = blueWhiteRed(50),
              textMatrix = textMatrix_expanded,
              setStdMargins = FALSE,
              cex.text = 0.5,
              zlim = c(-1, 1),
              main = "Module-trait relationships (genotype x stage combinations)")
dev.off()

# Write expanded correlation results to file
module_trait_corr_expanded <- data.frame(
  Module = rownames(moduleTraitCor_expanded),
  moduleTraitCor_expanded,
  stringsAsFactors = FALSE
)
write.csv(module_trait_corr_expanded, file = "module_trait_correlations_expanded.csv", row.names = FALSE)

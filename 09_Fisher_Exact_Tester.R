# Simplified Gene Set Enrichment Analysis
# This script tests for enrichment of differentially expressed genes in multiple gene sets
# using only the counts, without requiring actual gene IDs

# Load necessary libraries
library(ggplot2)
library(data.table)

#' Perform enrichment analysis using only gene counts
#' 
#' @param set_de Number of differentially expressed genes in the focused set
#' @param set_total Total number of genes in the focused set
#' @param genome_de Total number of differentially expressed genes genome-wide
#' @param genome_total Total number of genes evaluated genome-wide
#' @param set_name Name of the gene set (for reporting)
#' @return A data.frame with enrichment statistics
analyze_gene_set_counts <- function(set_de, set_total, genome_de, genome_total, set_name) {
  # Check inputs
  if (set_de > set_total) {
    stop("Number of DE genes in set cannot exceed total genes in set")
  }
  if (genome_de > genome_total) {
    stop("Number of DE genes in genome cannot exceed total genes in genome")
  }
  if (set_total > genome_total) {
    stop("Size of gene set cannot exceed total genes in genome")
  }
  
  # Create contingency table
  in_set_de <- set_de
  in_set_not_de <- set_total - set_de
  not_in_set_de <- genome_de - set_de
  not_in_set_not_de <- genome_total - genome_de - in_set_not_de
  
  # Create a contingency matrix
  cont_table <- matrix(c(in_set_de, in_set_not_de, 
                         not_in_set_de, not_in_set_not_de), 
                       nrow = 2, 
                       dimnames = list(
                         c("In Set", "Not in Set"),
                         c("DE", "Not DE")
                       ))
  
  # Perform Fisher's exact test
  fisher_result <- fisher.test(cont_table, alternative = "greater")
  
  # Calculate percentages and fold enrichment
  percent_de_in_set <- (in_set_de / set_total) * 100
  percent_de_genome <- (genome_de / genome_total) * 100
  fold_enrichment <- percent_de_in_set / percent_de_genome
  
  # Return results
  return(data.frame(
    Gene_Set = set_name,
    Total_Genes_In_Set = set_total,
    DE_Genes_In_Set = set_de,
    Percent_DE_In_Set = percent_de_in_set,
    Total_DE_Genes = genome_de,
    Total_Genes = genome_total,
    Percent_DE_Genome = percent_de_genome,
    Fold_Enrichment = fold_enrichment,
    Odds_Ratio = fisher_result$estimate,
    P_Value = fisher_result$p.value,
    FDR = NA,  # To be filled later after analyzing all sets
    Significant = NA,  # To be filled later
    stringsAsFactors = FALSE
  ))
}

#' Run enrichment analysis on multiple gene sets using count data
#' 
#' @param gene_set_data A data.frame with columns for each gene set:
#'        set_name, set_de, set_total, genome_de, genome_total
#' @param fdr_threshold FDR threshold for significance (default: 0.05)
#' @return A data.frame with enrichment statistics for all gene sets
run_multiple_enrichment_counts <- function(gene_set_data, fdr_threshold = 0.05) {
  # Check input data structure
  required_cols <- c("set_name", "set_de", "set_total", "genome_de", "genome_total")
  if (!all(required_cols %in% colnames(gene_set_data))) {
    stop("Input data.frame must contain columns: ", paste(required_cols, collapse=", "))
  }
  
  # Initialize results data frame
  results <- data.frame()
  
  # Process each row of the input data
  for (i in 1:nrow(gene_set_data)) {
    row <- gene_set_data[i, ]
    set_result <- analyze_gene_set_counts(
      row$set_de, 
      row$set_total, 
      row$genome_de, 
      row$genome_total, 
      row$set_name
    )
    results <- rbind(results, set_result)
  }
  
  # Apply FDR correction for multiple testing
  results$FDR <- p.adjust(results$P_Value, method = "BH")
  results$Significant <- results$FDR < fdr_threshold
  
  # Sort by FDR
  results <- results[order(results$FDR), ]
  
  return(results)
}

#' Plot enrichment results for all gene sets
#' 
#' @param results Output from run_multiple_enrichment_counts function
#' @param title Plot title
#' @return ggplot object
plot_enrichment <- function(results, title = "Gene Set Enrichment Analysis") {
  # Prepare data for plotting
  plot_data <- results
  plot_data$Gene_Set <- factor(plot_data$Gene_Set, levels = plot_data$Gene_Set[order(plot_data$Fold_Enrichment)])
  
  # Add significance indicators
  plot_data$sig_symbol <- ifelse(plot_data$FDR < 0.001, "***", 
                            ifelse(plot_data$FDR < 0.01, "**",
                            ifelse(plot_data$FDR < 0.05, "*", "ns")))
  
  # Create plot
  p <- ggplot(plot_data, aes(x = Gene_Set, y = Fold_Enrichment, fill = -log10(FDR))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_text(aes(label = sprintf("%.1f%% vs %.1f%% %s", 
                                Percent_DE_In_Set, 
                                Percent_DE_Genome,
                                sig_symbol)), 
              hjust = -0.1, size = 3) +
    labs(title = title,
         x = "Gene Set",
         y = "Fold Enrichment",
         fill = "-log10(FDR)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))
  
  return(p)
}

#' Export enrichment results to CSV file
#' 
#' @param results Output from run_multiple_enrichment_counts function
#' @param filename Output filename
#' @return NULL
export_results <- function(results, filename = "gene_set_enrichment_results.csv") {
  write.csv(results, file = filename, row.names = FALSE)
  message(paste("Results exported to", filename))
}

# Single set analysis function for quick use
analyze_single_set <- function(set_de, set_total, genome_de, genome_total, set_name = "Gene Set") {
  result <- analyze_gene_set_counts(set_de, set_total, genome_de, genome_total, set_name)
  
  # Print summary
  cat("\n===== Gene Set Enrichment Analysis =====\n")
  cat(sprintf("Gene Set: %s\n", set_name))
  cat(sprintf("DE genes in set: %d out of %d (%.2f%%)\n", 
              result$DE_Genes_In_Set, result$Total_Genes_In_Set, result$Percent_DE_In_Set))
  cat(sprintf("DE genes in genome: %d out of %d (%.2f%%)\n", 
              result$Total_DE_Genes, result$Total_Genes, result$Percent_DE_Genome))
  cat(sprintf("Fold Enrichment: %.2f\n", result$Fold_Enrichment))
  cat(sprintf("Odds Ratio: %.2f\n", result$Odds_Ratio))
  cat(sprintf("P-value: %.2e\n", result$P_Value))
  cat(sprintf("Significant: %s\n", ifelse(result$P_Value < 0.05, "YES", "NO")))
  cat("=======================================\n")
  
  return(result)
}

#=====================================================================
# Example usage:
#=====================================================================

#' Example 1: Quick analysis of a single gene set (your original example)
result <- analyze_single_set(
  set_de = 40,           # DE genes in the focused set
  set_total = 49,        # Total genes in the focused set
  genome_de = 7238,      # Total DE genes genome-wide
  genome_total = 13920,  # Total genes evaluated genome-wide
  set_name = "My Gene Set"
)

#' Example 2: Analyze multiple gene sets at once
# Create a data frame with the counts for multiple gene sets
gene_set_data <- data.frame(
  set_name = c("Gene Set A", "Gene Set B", "Gene Set C", "Gene Set D"),
  set_de = c(40, 25, 60, 15),
  set_total = c(49, 50, 100, 20),
  genome_de = c(7238, 7238, 7238, 7238),
  genome_total = c(13920, 13920, 13920, 13920)
)

# Run the analysis
results <- run_multiple_enrichment_counts(gene_set_data)
print(results)

# Plot results
p <- plot_enrichment(results)
print(p)

# Export results
export_results(results, "multiple_gene_sets_results.csv")

#=====================================================================
# How to use this script with your own data:
#=====================================================================

# Method 1: Edit the gene_set_data data frame above with your own values

# Method 2: Load data from a CSV file
# Your CSV file should have columns: set_name, set_de, set_total, genome_de, genome_total
# Example:
# my_data <- read.csv("my_gene_sets.csv")
# results <- run_multiple_enrichment_counts(my_data)
# p <- plot_enrichment(results)
# print(p)
# export_results(results, "my_results.csv")

# Method 3: Create an R function to interactively analyze sets
analyze_new_set <- function() {
  set_name <- readline(prompt="Enter gene set name: ")
  set_de <- as.numeric(readline(prompt="Enter number of DE genes in set: "))
  set_total <- as.numeric(readline(prompt="Enter total genes in set: "))
  genome_de <- as.numeric(readline(prompt="Enter total DE genes in genome: "))
  genome_total <- as.numeric(readline(prompt="Enter total genes in genome: "))
  
  return(analyze_single_set(set_de, set_total, genome_de, genome_total, set_name))
}

# Uncomment to run interactive mode:
# analyze_new_set()

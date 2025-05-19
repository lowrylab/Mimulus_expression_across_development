# Define output directory (adjust as needed)
output_dir <- "wgcna_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# This script assumes 'moduleTraitCor_expanded' and 'moduleTraitPvalue_expanded' 
# already exist in your workspace (from your previous WGCNA analysis)

# Apply Benjamini-Hochberg correction to adjust for multiple testing
pValueMatrix_adjusted <- matrix(
  p.adjust(as.vector(moduleTraitPvalue_expanded), method = "BH"),
  nrow = nrow(moduleTraitPvalue_expanded),
  ncol = ncol(moduleTraitPvalue_expanded)
)

# Preserve the dimension names
dimnames(pValueMatrix_adjusted) <- dimnames(moduleTraitPvalue_expanded)

# Create a data frame combining correlation values and adjusted p-values
module_trait_results <- data.frame(
  Module = rownames(moduleTraitCor_expanded)
)

# Add correlation values with trait names as column headers
for (trait in colnames(moduleTraitCor_expanded)) {
  module_trait_results[[paste0("cor_", trait)]] <- moduleTraitCor_expanded[, trait]
}

# Add adjusted p-values with trait names as column headers
for (trait in colnames(pValueMatrix_adjusted)) {
  module_trait_results[[paste0("adj.p_", trait)]] <- pValueMatrix_adjusted[, trait]
}

# Write results to a CSV file
write.csv(
  module_trait_results,
  file = file.path(output_dir, "module_trait_correlations_adjusted_pvalues.csv"),
  row.names = FALSE
)

# Create a heatmap text matrix with the adjusted p-values
# Format to 2 decimal places for better readability
textMatrix_adjusted <- paste(
  sprintf("%.2f", moduleTraitCor_expanded), 
  "\n(", 
  sprintf("%.3f", pValueMatrix_adjusted), 
  ")", 
  sep = ""
)
dim(textMatrix_adjusted) <- dim(moduleTraitCor_expanded)

# Create a new heatmap with adjusted p-values and improved visualization
pdf(file.path(output_dir, "module_trait_relationships_adjusted_pvalues.pdf"), 
    width = 16, height = 12)  # Increased width for more space

# Set custom margins to prevent text from being cut off
# Increase left margin (4) significantly
par(mar = c(6, 12, 4, 2) + 0.1)  # Bottom, left, top, right

# Create a heatmap with improved visual settings
WGCNA::labeledHeatmap(
  Matrix = moduleTraitCor_expanded,
  xLabels = colnames(moduleTraitCor_expanded),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = WGCNA::blueWhiteRed(50),
  textMatrix = textMatrix_adjusted,
  setStdMargins = FALSE,  # Don't use standard margins
  cex.text = 0.9,         # Increased text size inside cells (was 0.5)
  cex.lab = 1.1,          # Increased label font size
  cex.axis = 1.1,         # Increased axis font size
  zlim = c(-1, 1),
  main = "Module-trait relationships (with FDR correction)",
  xLabelsAngle = 45,      # Angled x-axis labels for better readability
  yLabelsAdj = 1,         # Right-align y labels
  margins = c(8, 10)      # Additional margin adjustment
)

dev.off()

# Create a second version with different coloring to highlight significant correlations
pdf(file.path(output_dir, "module_trait_relationships_significant.pdf"), 
    width = 16, height = 12)

par(mar = c(6, 12, 4, 2) + 0.1)

# Create a matrix that highlights significant correlations after adjustment
# This will make cells with p < 0.05 after adjustment more visible
significantMatrix <- moduleTraitCor_expanded * (pValueMatrix_adjusted < 0.05)

# Create text that includes asterisks for significant correlations
textMatrix_significant <- paste(
  sprintf("%.2f", moduleTraitCor_expanded), 
  ifelse(pValueMatrix_adjusted < 0.05, "*", ""),  # Add asterisk for significant values
  "\n(", 
  sprintf("%.3f", pValueMatrix_adjusted), 
  ")", 
  sep = ""
)
dim(textMatrix_significant) <- dim(moduleTraitCor_expanded)

WGCNA::labeledHeatmap(
  Matrix = moduleTraitCor_expanded,
  xLabels = colnames(moduleTraitCor_expanded),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = WGCNA::blueWhiteRed(50),
  textMatrix = textMatrix_significant,
  setStdMargins = FALSE,
  cex.text = 0.9,
  cex.lab = 1.1,
  cex.axis = 1.1,
  zlim = c(-1, 1),
  main = "Module-trait relationships (significant correlations marked with *)",
  xLabelsAngle = 45,
  yLabelsAdj = 1,
  margins = c(8, 10)
)

# Add a legend
legend("bottomright", 
       legend = c("* = significant (adjusted p < 0.05)"),
       cex = 0.8, 
       bg = "white")

dev.off()

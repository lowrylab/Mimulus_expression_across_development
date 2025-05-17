# Generalized Script to Visualize Multiple WGCNA Modules

# Load required libraries
library(WGCNA)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(gridExtra)

# Define the module colors to analyze
module_colors <- c("red", "blue", "turquoise", "green", "greenyellow", 
                  "cyan", "black", "lightcyan", "tan")

# Define output directory
output_dir <- "wgcna_results/module_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# This script assumes the following objects exist from your WGCNA analysis:
# - datExpr (transposed expression data with samples as rows)
# - moduleColors (module assignments)
# - filtered_metadata (sample metadata)
# - MEs (module eigengenes, if available)

# Function to analyze and visualize a single module
analyze_module <- function(module_color) {
  
  # Create module-specific output directory
  module_dir <- file.path(output_dir, paste0(module_color, "_module"))
  if (!dir.exists(module_dir)) {
    dir.create(module_dir, recursive = TRUE)
  }
  
  # Extract genes in the current module
  module_genes <- colnames(datExpr)[moduleColors == module_color]
  cat("Number of genes in", module_color, "module:", length(module_genes), "\n")
  
  # Skip if module has no genes
  if (length(module_genes) == 0) {
    cat("No genes found in", module_color, "module. Skipping...\n")
    return(NULL)
  }
  
  # Extract expression data for module genes
  module_expr <- datExpr[, module_genes]
  
  # Get the module eigengene
  ME_name <- paste0("ME", module_color)
  if (exists("MEs") && ME_name %in% colnames(MEs)) {
    module_ME <- MEs[, ME_name]
    cat("Using existing module eigengene from MEs object\n")
  } else {
    # Calculate the first principal component of the expression data
    cat("Calculating module eigengene using principal component analysis\n")
    pca <- prcomp(module_expr, scale. = TRUE)
    module_ME <- pca$x[, 1]
    names(module_ME) <- rownames(module_expr)
    
    # Ensure the eigengene correlates positively with most genes
    mean_cor <- mean(cor(module_expr, module_ME))
    if (mean_cor < 0) {
      module_ME <- -module_ME  # Flip sign if negatively correlated
    }
  }
  
  # Create a data frame for plotting
  expr_data <- data.frame(
    SampleID = rownames(module_expr),
    module_ME = as.numeric(module_ME),
    Genotype = filtered_metadata$genotype,
    Stage = filtered_metadata$development_stage,
    Genotype_Stage = paste(filtered_metadata$genotype, filtered_metadata$development_stage, sep="_")
  )
  
  # Add the expression data for each gene (limited to first 50 genes if there are many)
  genes_to_include <- module_genes
  if (length(module_genes) > 50) {
    # Calculate gene module membership
    gene_MM <- as.data.frame(cor(module_expr, module_ME, use = "p"))
    colnames(gene_MM) <- "ModuleMembership"
    
    # Sort genes by module membership and take top 50
    sorted_genes <- rownames(gene_MM)[order(abs(gene_MM$ModuleMembership), decreasing = TRUE)]
    genes_to_include <- sorted_genes[1:50]
    cat("Too many genes in module. Including only top 50 genes by module membership.\n")
  }
  
  for (gene in genes_to_include) {
    expr_data[[gene]] <- module_expr[, gene]
  }
  
  # Melt the data for easier plotting
  melted_data <- melt(expr_data, 
                     id.vars = c("SampleID", "Genotype", "Stage", "Genotype_Stage", "module_ME"),
                     variable.name = "Gene", 
                     value.name = "Expression")
  
  # 1. Create boxplot of module eigengene by genotype
  p1 <- ggplot(expr_data, aes(x = Genotype, y = module_ME, fill = Genotype)) +
    geom_boxplot() +
    theme_bw() +
    labs(title = paste0(capitalize(module_color), " Module Eigengene by Genotype"),
         y = "Module Eigengene",
         x = "Genotype") +
    scale_fill_manual(values = c("Coast" = "skyblue", "Inland" = "salmon")) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 13))
  
  # 2. Create boxplot of module eigengene by developmental stage
  p2 <- ggplot(expr_data, aes(x = Stage, y = module_ME, fill = Stage)) +
    geom_boxplot() +
    theme_bw() +
    labs(title = paste0(capitalize(module_color), " Module Eigengene by Developmental Stage"),
         y = "Module Eigengene",
         x = "Developmental Stage") +
    scale_fill_brewer(palette = "Set2") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 13))
  
  # 3. Create interaction plot of module eigengene
  p3 <- ggplot(expr_data, aes(x = Stage, y = module_ME, color = Genotype, group = Genotype)) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    theme_bw() +
    labs(title = paste0(capitalize(module_color), " Module Eigengene Interaction Effect"),
         y = "Module Eigengene (Mean ± SE)",
         x = "Developmental Stage") +
    scale_color_manual(values = c("Coast" = "blue", "Inland" = "red")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          legend.position = "right")
  
  # 4. Create boxplot of module eigengene by genotype-stage interaction
  p4 <- ggplot(expr_data, aes(x = Genotype_Stage, y = module_ME, fill = Genotype)) +
    geom_boxplot() +
    theme_bw() +
    labs(title = paste0(capitalize(module_color), " Module Eigengene by Genotype-Stage Interaction"),
         y = "Module Eigengene",
         x = "Genotype-Stage Combination") +
    scale_fill_manual(values = c("Coast" = "skyblue", "Inland" = "salmon")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.title = element_text(size = 13))
  
  # Combine the plots
  combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  # Save the combined plot
  ggsave(file.path(module_dir, paste0(module_color, "_module_eigengene_plots.pdf")), 
         combined_plot, width = 14, height = 10)
  
  # 5. Create a heatmap of all genes in the module
  # Calculate mean expression by group for a cleaner heatmap
  mean_expr_by_group <- matrix(0, 
                             nrow = length(module_genes),
                             ncol = length(unique(expr_data$Genotype_Stage)))
  rownames(mean_expr_by_group) <- module_genes
  colnames(mean_expr_by_group) <- unique(expr_data$Genotype_Stage)
  
  for (gene in module_genes) {
    for (group in unique(expr_data$Genotype_Stage)) {
      samples_in_group <- expr_data$SampleID[expr_data$Genotype_Stage == group]
      mean_expr_by_group[gene, group] <- mean(module_expr[samples_in_group, gene])
    }
  }
  
  # Scale the data by row (gene) for better visualization
  scaled_mean_expr <- t(scale(t(mean_expr_by_group)))
  
  # Calculate gene module membership
  gene_MM <- as.data.frame(cor(module_expr, module_ME, use = "p"))
  colnames(gene_MM) <- "ModuleMembership"
  
  # Create annotation dataframe
  gene_annotation <- data.frame(
    ModuleMembership = gene_MM$ModuleMembership,
    row.names = rownames(gene_MM)
  )
  
  # Create sample annotations
  sample_annotation <- data.frame(
    Genotype = sapply(strsplit(colnames(mean_expr_by_group), "_"), "[", 1),
    Stage = sapply(strsplit(colnames(mean_expr_by_group), "_"), "[", 2),
    row.names = colnames(mean_expr_by_group)
  )
  
  # Define annotation colors
  ann_colors <- list(
    Genotype = c(Coast = "skyblue", Inland = "salmon"),
    Stage = c(`2leaf` = "#66C2A5", `4leaf` = "#FC8D62", bud = "#8DA0CB"),
    ModuleMembership = colorRampPalette(c("white", "darkblue"))(100)
  )
  
  # Create the heatmap
  pdf(file.path(module_dir, paste0(module_color, "_module_genes_heatmap.pdf")), 
      width = 12, height = min(16, max(8, length(module_genes)/15)))
  
  # Determine if we should show row names based on number of genes
  show_names <- length(module_genes) <= 40
  
  pheatmap(scaled_mean_expr,
           annotation_col = sample_annotation,
           annotation_row = gene_annotation,
           annotation_colors = ann_colors,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = show_names,
           main = paste0(capitalize(module_color), " Module Genes Expression Across Genotype-Stage Combinations"),
           fontsize = 12,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  dev.off()
  
  cat("Analysis for", module_color, "module completed and saved to", module_dir, "\n")
  return(TRUE)
}

# Helper function to capitalize first letter of module color name for plot titles
capitalize <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

# Generate a combined summary plot of all modules
create_module_summary <- function(module_colors) {
  plots_list <- list()
  
  # Get all module eigengenes
  all_MEs <- data.frame(row.names = rownames(datExpr))
  
  for (module in module_colors) {
    ME_name <- paste0("ME", module)
    
    # Get the module eigengene
    if (exists("MEs") && ME_name %in% colnames(MEs)) {
      all_MEs[[ME_name]] <- MEs[, ME_name]
    } else {
      # Extract genes in the module
      module_genes <- colnames(datExpr)[moduleColors == module]
      
      # Skip if module has no genes
      if (length(module_genes) == 0) {
        cat("No genes found in", module, "module. Skipping from summary...\n")
        next
      }
      
      # Extract expression data for module genes
      module_expr <- datExpr[, module_genes]
      
      # Calculate the first principal component
      pca <- prcomp(module_expr, scale. = TRUE)
      module_ME <- pca$x[, 1]
      names(module_ME) <- rownames(module_expr)
      
      # Ensure the eigengene correlates positively with most genes
      mean_cor <- mean(cor(module_expr, module_ME))
      if (mean_cor < 0) {
        module_ME <- -module_ME  # Flip sign if negatively correlated
      }
      
      all_MEs[[ME_name]] <- module_ME
    }
    
    # Add to metadata for plotting
    all_MEs$Genotype <- filtered_metadata$genotype
    all_MEs$Stage <- filtered_metadata$development_stage
    all_MEs$Genotype_Stage <- paste(filtered_metadata$genotype, filtered_metadata$development_stage, sep="_")
    
    # Create interaction plot for this module
    plots_list[[module]] <- ggplot(all_MEs, 
                                  aes(x = Stage, y = .data[[ME_name]], 
                                      color = Genotype, group = Genotype)) +
      stat_summary(fun = mean, geom = "point", size = 2) +
      stat_summary(fun = mean, geom = "line", size = 0.8) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
      theme_bw() +
      labs(title = capitalize(module),
           y = "Module Eigengene") +
      scale_color_manual(values = c("Coast" = "blue", "Inland" = "red")) +
      theme(plot.title = element_text(hjust = 0.5, size = 12),
            axis.text = element_text(size = 10),
            axis.title.y = element_text(size = 11),
            axis.title.x = element_blank(),
            legend.position = "none")
  }
  
  # Calculate dimensions for the grid
  n_modules <- length(plots_list)
  n_cols <- min(3, n_modules)
  n_rows <- ceiling(n_modules / n_cols)
  
  # Create a common legend
  legend_plot <- ggplot(all_MEs, aes(x = Stage, y = Genotype, color = Genotype)) +
    geom_point() +
    scale_color_manual(values = c("Coast" = "blue", "Inland" = "red")) +
    theme(legend.position = "bottom")
  
  # Extract the legend
  legend <- get_legend(legend_plot)
  
  # Create a title text grob
  title <- textGrob("Module Eigengenes Across Genotype and Developmental Stage", 
                    gp = gpar(fontsize = 16, fontface = "bold"))
  
  # Combine the plots with a title and shared legend
  combined_summary <- grid.arrange(
    title,
    do.call(arrangeGrob, c(plots_list, ncol = n_cols)),
    legend,
    heights = c(0.1, 0.8, 0.1),
    nrow = 3
  )
  
  # Save the combined summary
  ggsave(file.path(output_dir, "all_modules_eigengene_summary.pdf"), 
         combined_summary, width = min(21, n_cols * 7), height = min(29.7, n_rows * 5 + 2))
  
  cat("Module summary plot created and saved to", output_dir, "\n")
}

# Helper function to extract the legend from a ggplot
get_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) {
    legend <- tmp$grobs[[leg]]
    return(legend)
  } else {
    return(NULL)
  }
}

# Create a function to generate a comprehensive module report
generate_module_report <- function() {
  # Create a PDF with summary information
  pdf(file.path(output_dir, "wgcna_module_report.pdf"), width = 10, height = 12)
  
  # 1. Module sizes bar plot
  module_sizes <- table(moduleColors)
  filtered_sizes <- module_sizes[names(module_sizes) %in% module_colors]
  
  barplot(filtered_sizes, 
          main = "Number of Genes in Each Module",
          col = names(filtered_sizes),
          las = 2,
          cex.names = 0.8,
          ylab = "Number of Genes")
  
  # 2. Create a correlation heatmap between module eigengenes
  all_MEs <- data.frame(row.names = rownames(datExpr))
  
  for (module in module_colors) {
    ME_name <- paste0("ME", module)
    
    # Get the module eigengene
    if (exists("MEs") && ME_name %in% colnames(MEs)) {
      all_MEs[[ME_name]] <- MEs[, ME_name]
    } else {
      # Extract genes in the module
      module_genes <- colnames(datExpr)[moduleColors == module]
      
      # Skip if module has no genes
      if (length(module_genes) == 0) {
        next
      }
      
      # Extract expression data for module genes
      module_expr <- datExpr[, module_genes]
      
      # Calculate the first principal component
      pca <- prcomp(module_expr, scale. = TRUE)
      module_ME <- pca$x[, 1]
      names(module_ME) <- rownames(module_expr)
      
      # Ensure the eigengene correlates positively with most genes
      mean_cor <- mean(cor(module_expr, module_ME))
      if (mean_cor < 0) {
        module_ME <- -module_ME  # Flip sign if negatively correlated
      }
      
      all_MEs[[ME_name]] <- module_ME
    }
  }
  
  # Calculate correlation matrix between module eigengenes
  me_cors <- cor(all_MEs, use = "pairwise.complete.obs")
  
  # Plot heatmap of module eigengene correlations
  par(mar = c(6, 6, 6, 6))
  heatmap(me_cors, 
          main = "Module Eigengene Correlation Heatmap",
          col = colorRampPalette(c("blue", "white", "red"))(100),
          cexRow = 0.9, 
          cexCol = 0.9)
  
  # 3. ANOVA analysis for each module eigengene
  # Create a table to store ANOVA results
  anova_results <- matrix(NA, nrow = length(module_colors), ncol = 4)
  rownames(anova_results) <- module_colors
  colnames(anova_results) <- c("Genotype p-value", "Stage p-value", 
                             "Interaction p-value", "R-squared")
  
  # Perform ANOVA for each module
  for (i in 1:length(module_colors)) {
    module <- module_colors[i]
    ME_name <- paste0("ME", module)
    
    # Skip if module eigengene is not available
    if (!ME_name %in% colnames(all_MEs)) {
      next
    }
    
    # Create a data frame for ANOVA
    anova_data <- data.frame(
      ME = all_MEs[[ME_name]],
      Genotype = filtered_metadata$genotype,
      Stage = filtered_metadata$development_stage
    )
    
    # Fit ANOVA model
    model <- aov(ME ~ Genotype * Stage, data = anova_data)
    anova_table <- summary(model)[[1]]
    
    # Extract p-values
    anova_results[i, 1] <- anova_table[1, 5]  # Genotype p-value
    anova_results[i, 2] <- anova_table[2, 5]  # Stage p-value
    anova_results[i, 3] <- anova_table[3, 5]  # Interaction p-value
    
    # Calculate R-squared
    ss_total <- sum(anova_table[, 2])
    ss_residual <- anova_table[4, 2]
    anova_results[i, 4] <- 1 - (ss_residual / ss_total)
  }
  
  # Format p-values for display
  formatted_results <- anova_results
  for (i in 1:3) {
    formatted_results[, i] <- ifelse(anova_results[, i] < 0.001, "< 0.001",
                                   ifelse(anova_results[, i] < 0.01, "< 0.01",
                                         ifelse(anova_results[, i] < 0.05, "< 0.05",
                                               round(anova_results[, i], 3))))
  }
  
  # Format R-squared for display
  formatted_results[, 4] <- round(as.numeric(anova_results[, 4]), 3)
  
  # Display ANOVA results table
  textplot(formatted_results, 
           show.rownames = TRUE,
           cex = 0.8,
           halign = "center",
           valign = "center",
           mar = c(1, 1, 4, 1))
  title("ANOVA Results for Module Eigengenes", cex.main = 1.2)
  
  dev.off()
  
  cat("Module report generated and saved to", output_dir, "\n")
}

# Function to run the full analysis
run_module_analysis <- function() {
  # Process each module
  for (module in module_colors) {
    cat("\nAnalyzing", module, "module...\n")
    analyze_module(module)
  }
  
  # Create summary visualizations
  create_module_summary(module_colors)
  
  # Generate comprehensive report
  generate_module_report()
  
  cat("\nAll analyses completed!\n")
}

# Run the analysis
run_module_analysis()

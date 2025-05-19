#!/usr/bin/env python3
"""
Gene Enrichment Analysis Tool

This script performs gene enrichment analysis on lists of differentially expressed genes
using various annotation types (GO, Pfam, Panther, KOG, KEGG, etc.) and calculates
statistical significance using the hypergeometric test.

Input:
- Annotation file: Tab-delimited file with gene annotations
- Significance file: CSV file with differentially expressed genes

Output:
- Enrichment results: CSV file with enrichment statistics
- Optional visualization: PDF plot of enrichment results
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict, Counter
from typing import Dict, List, Set, Tuple, Optional, Union

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# Function to check gene overlap between significance and annotation files
def check_gene_overlap(significant_genes: List[str], annotation_genes: Set[str], output_prefix: str = None):
    """
    Check and report the overlap between gene lists from significance and annotation files.
    Optionally writes non-overlapping genes to files for inspection.
    """
    # Convert lists to sets for efficient comparisons
    sig_genes_set = set(significant_genes)
    
    # Calculate overlaps
    overlapping = sig_genes_set.intersection(annotation_genes)
    sig_only = sig_genes_set - annotation_genes
    anno_only = annotation_genes - sig_genes_set
    
    # Calculate percentages
    overlap_percent = len(overlapping) / len(sig_genes_set) * 100 if sig_genes_set else 0
    
    # Print summary
    print("\n--- Gene Overlap Analysis ---")
    print(f"Total significant genes: {len(sig_genes_set)}")
    print(f"Total annotation genes: {len(annotation_genes)}")
    print(f"Overlapping genes: {len(overlapping)} ({overlap_percent:.1f}% of significant genes)")
    print(f"Genes only in significance file: {len(sig_only)}")
    print(f"Genes only in annotation file: {len(annotation_genes) - len(overlapping)}")
    
    # Show sample of non-matching genes
    if sig_only:
        sample_size = min(5, len(sig_only))
        print(f"\nSample of significant genes NOT found in annotation file:")
        for gene in list(sig_only)[:sample_size]:
            print(f"  - {gene}")
    
    # Write non-matching genes to files if output prefix is provided
    if output_prefix and sig_only:
        non_matching_file = f"{output_prefix}_non_matching_genes.txt"
        with open(non_matching_file, 'w') as f:
            f.write("# Significant genes not found in annotation file\n")
            for gene in sorted(sig_only):
                f.write(f"{gene}\n")
        print(f"\nWrote {len(sig_only)} non-matching genes to {non_matching_file}")
    
    return {
        "total_significant": len(sig_genes_set),
        "total_annotation": len(annotation_genes),
        "overlapping": len(overlapping),
        "overlap_percent": overlap_percent,
        "sig_only": len(sig_only),
        "anno_only": len(anno_only),
        "non_matching_genes": list(sig_only)
    }

# Set up command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Gene Enrichment Analysis Tool')
    parser.add_argument('--annotation', '-a', required=True, help='Path to annotation file (tab-delimited)')
    parser.add_argument('--significance', '-s', required=True, help='Path to significance file (CSV)')
    parser.add_argument('--output', '-o', default='enrichment_results', help='Output prefix for result files')
    parser.add_argument('--type', '-t', default='GO', choices=['GO', 'Pfam', 'Panther', 'KOG', 'KO', 'KEGG/ec'],
                        help='Annotation type to use for enrichment analysis')
    parser.add_argument('--pvalue', '-p', type=float, default=0.05, help='Adjusted p-value cutoff for significance')
    parser.add_argument('--log2fc', '-l', type=float, default=0, help='Log2 fold change cutoff (absolute value)')
    parser.add_argument('--visualize', '-v', action='store_true', help='Generate visualization of results')
    parser.add_argument('--top', type=int, default=20, help='Number of top terms to include in visualization')
    parser.add_argument('--correction', '-c', default='fdr_bh', choices=['fdr_bh', 'bonferroni', 'none'],
                        help='Method for multiple testing correction')
    parser.add_argument('--remove-version', action='store_true', default=True, 
                       help='Remove version suffix (e.g., .v2.0) from gene IDs (default: True)')
    parser.add_argument('--force', '-f', action='store_true', 
                       help='Force analysis even with low gene overlap (bypass confirmation prompt)')
    parser.add_argument('--export-gene-lists', action='store_true',
                       help='Export gene lists for each enriched term to separate files')
    return parser.parse_args()

# Function to parse the annotation file
def parse_annotation_file(file_path: str) -> List[Dict]:
    """Parse the gene annotation file in tab-delimited format."""
    annotations = []
    with open(file_path, 'r') as f:
        # Read header line
        header = f.readline().strip().split('\t')
        
        # Remove '#' from first column name if present
        if header[0].startswith('#'):
            header[0] = header[0][1:]
        
        # Read content lines
        for line in f:
            if not line.strip():
                continue
            
            values = line.strip().split('\t')
            entry = {}
            
            # Map values to column names
            for i, column in enumerate(header):
                if i < len(values):
                    # Remove quotes if present
                    value = values[i]
                    if value.startswith('"') and value.endswith('"'):
                        value = value[1:-1]
                    entry[column] = value
                else:
                    entry[column] = ''
            
            annotations.append(entry)
    
    return annotations

# Function to parse the significance file
def parse_significance_file(file_path: str, log2fc_cutoff: float = 0, pvalue_cutoff: float = 0.05, remove_version: bool = True) -> List[str]:
    """Parse the significance file and extract significant gene IDs based on adjusted p-values."""
    # Read CSV file with pandas
    df = pd.read_csv(file_path)
    
    # Handle empty column name for gene IDs
    id_column = df.columns[0]  # First column contains gene IDs
    
    # Prioritize adjusted p-values for filtering
    if 'adj.P.Val' in df.columns:
        # Use adjusted p-values for significance filtering
        if 'logFC' in df.columns:
            significant = df[(df['adj.P.Val'] <= pvalue_cutoff) & (abs(df['logFC']) >= log2fc_cutoff)]
        else:
            significant = df[df['adj.P.Val'] <= pvalue_cutoff]
    # Fall back to raw p-values if adjusted are not available
    elif 'P.Value' in df.columns:
        if 'logFC' in df.columns:
            significant = df[(df['P.Value'] <= pvalue_cutoff) & (abs(df['logFC']) >= log2fc_cutoff)]
        else:
            significant = df[df['P.Value'] <= pvalue_cutoff]
    # If no p-value columns found, filter by fold change only
    elif 'logFC' in df.columns:
        significant = df[abs(df['logFC']) >= log2fc_cutoff]
    else:
        significant = df  # No filtering if required columns are missing
    
    # Extract gene IDs
    gene_ids = significant[id_column].tolist()
    
    # Clean gene IDs (remove quotes and optionally version suffixes like .v2.0)
    clean_gene_ids = []
    for gene_id in gene_ids:
        # Convert to string and remove quotes
        gene_id = str(gene_id).strip('"\'')
        
        # Remove version suffix if requested (like .v2.0)
        if remove_version and ".v" in gene_id:
            # Split at .v and take the first part
            parts = gene_id.split(".v")
            gene_id = parts[0]
        
        clean_gene_ids.append(gene_id)
    
    if gene_ids:
        print(f"Example of original gene ID vs cleaned: {gene_ids[0]} -> {clean_gene_ids[0]}")
        
        # Check if cleaning made a difference
        changed_count = sum(1 for i, g in enumerate(gene_ids) if str(g) != clean_gene_ids[i])
        if changed_count > 0:
            print(f"Cleaned {changed_count} gene IDs ({changed_count/len(gene_ids)*100:.1f}%)")
    
    return clean_gene_ids

# Function to prepare data for enrichment analysis
def prepare_enrichment_data(
        annotations: List[Dict], 
        gene_ids: List[str], 
        term_type: str
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Set[str], Set[str]]:
    """Prepare data structures for enrichment analysis."""
    # Get all unique genes in the annotation data
    all_genes = set()
    for annotation in annotations:
        if 'locusName' in annotation:
            all_genes.add(annotation['locusName'])
    
    # Create a mapping from gene to terms
    gene_to_terms = defaultdict(set)
    for annotation in annotations:
        gene = annotation.get('locusName', '')
        if not gene:
            continue
        
        terms_str = annotation.get(term_type, '')
        if not terms_str:
            continue
        
        # Split terms and clean them
        terms = [term.strip('"\'') for term in terms_str.split(',')]
        gene_to_terms[gene].update(t for t in terms if t)
    
    # Create a mapping from term to genes
    term_to_genes = defaultdict(set)
    for gene, terms in gene_to_terms.items():
        for term in terms:
            term_to_genes[term].add(gene)
    
    # Get set of significant genes that are in the annotation data
    # Note: Genes from significance file have already been cleaned to remove version suffixes
    significant_genes = set(gene_ids).intersection(all_genes)
    
    # Print statistics about gene matching
    print(f"Input genes: {len(gene_ids)}")
    print(f"Total annotation genes: {len(all_genes)}")
    print(f"Matched genes: {len(significant_genes)} ({len(significant_genes)/len(gene_ids)*100:.1f}%)")
    
    # Print sample of genes that weren't found in annotations
    not_found_genes = set(gene_ids) - all_genes
    if not_found_genes:
        sample_size = min(5, len(not_found_genes))
        print(f"Sample of genes not found in annotations: {list(not_found_genes)[:sample_size]}")
    
    return gene_to_terms, term_to_genes, all_genes, significant_genes

# Function to calculate enrichment using hypergeometric test
def calculate_enrichment(
        term_to_genes: Dict[str, Set[str]],
        significant_genes: Set[str],
        all_genes: Set[str],
        pvalue_cutoff: float = 0.05
) -> List[Dict]:
    """Calculate enrichment statistics for each term."""
    results = []
    total_sig_genes = len(significant_genes)
    total_bg_genes = len(all_genes)
    
    for term, genes in term_to_genes.items():
        # Skip terms with no genes
        if not genes:
            continue
        
        # Calculate overlaps
        term_sig_genes = genes.intersection(significant_genes)
        if not term_sig_genes:
            continue
        
        term_bg_genes = genes
        
        # Prepare contingency table for hypergeometric test
        x = len(term_sig_genes)  # Hits in selection
        n = total_sig_genes      # Selection size
        M = len(term_bg_genes)   # Hits in background
        N = total_bg_genes       # Background size
        
        # Calculate p-value using hypergeometric test
        # p-value = P(X >= x), where X follows hypergeometric distribution
        pvalue = stats.hypergeom.sf(x-1, N, M, n)
        
        # Calculate fold enrichment
        expected = (M / N) * n
        fold_enrichment = x / expected if expected > 0 else float('inf')
        
        # Calculate odds ratio
        a = x  # Term in significant
        b = M - x  # Term in background but not in significant
        c = n - x  # Significant without term
        d = N - n - b  # Background without term and not in significant
        
        odds_ratio = (a * d) / (b * c) if (b * c) > 0 else float('inf')
        
        # Store result
        result = {
            'term': term,
            'genes_in_selection': x,
            'genes_in_background': M,
            'percent_in_selection': (x / n * 100) if n > 0 else 0,
            'percent_in_background': (M / N * 100) if N > 0 else 0,
            'fold_enrichment': fold_enrichment,
            'odds_ratio': odds_ratio,
            'pvalue': pvalue,
            'genes': ','.join(term_sig_genes)
        }
        
        results.append(result)
    
    # Sort by p-value
    results.sort(key=lambda r: r['pvalue'])
    
    return results

# Function to apply multiple testing correction - using custom implementation
def apply_multiple_testing_correction(results: List[Dict], method: str = 'fdr_bh') -> List[Dict]:
    """Apply multiple testing correction to the p-values."""
    if not results:
        return results
    
    # Extract p-values
    pvalues = [r['pvalue'] for r in results]
    
    # Apply correction
    if method == 'fdr_bh':
        # Benjamini-Hochberg FDR
        adjusted_pvalues = benjamini_hochberg_correction(pvalues)
        rejected = [p <= 0.05 for p in adjusted_pvalues]
    elif method == 'bonferroni':
        # Bonferroni correction
        adjusted_pvalues = [min(p * len(pvalues), 1.0) for p in pvalues]
        rejected = [p <= 0.05 for p in adjusted_pvalues]
    elif method == 'none':
        # No correction
        adjusted_pvalues = pvalues
        rejected = [p <= 0.05 for p in pvalues]
    else:
        # Default to Benjamini-Hochberg FDR
        adjusted_pvalues = benjamini_hochberg_correction(pvalues)
        rejected = [p <= 0.05 for p in adjusted_pvalues]
    
    # Update results with corrected p-values
    for i, r in enumerate(results):
        r['adjusted_pvalue'] = adjusted_pvalues[i]
        r['significant'] = rejected[i]
    
    return results

# Implementation of the Benjamini-Hochberg correction
def benjamini_hochberg_correction(pvalues):
    """
    Applies Benjamini-Hochberg FDR correction to p-values.
    
    Args:
        pvalues: List of p-values to correct
        
    Returns:
        List of corrected p-values (FDR q-values)
    """
    n = len(pvalues)
    if n == 0:
        return []
    
    # Create index-value pairs
    indexed_pvalues = list(enumerate(pvalues))
    
    # Sort by p-value
    indexed_pvalues.sort(key=lambda x: x[1])
    
    # Calculate ranks (1-based)
    ranks = {idx: rank + 1 for rank, (idx, _) in enumerate(indexed_pvalues)}
    
    # Calculate adjusted p-values
    adjusted_pvalues = [0] * n
    
    # Start with the highest p-value
    prev_adjusted = indexed_pvalues[-1][1]
    adjusted_pvalues[indexed_pvalues[-1][0]] = prev_adjusted
    
    # Process remaining p-values from highest to lowest
    for i in range(n-2, -1, -1):
        idx, pval = indexed_pvalues[i]
        rank = ranks[idx]
        adjusted = min(prev_adjusted, pval * n / rank)
        adjusted_pvalues[idx] = adjusted
        prev_adjusted = adjusted
    
    return adjusted_pvalues
# Function to visualize enrichment results
def visualize_enrichment(results: List[Dict], output_file: str, top_n: int = 20):
    """Create a visualization of enrichment results."""
    if not results:
        print("No results to visualize.")
        return
    
    # Select top N terms by p-value
    top_results = sorted(results, key=lambda r: r['pvalue'])[:top_n]
    
    # Prepare data for plotting
    terms = [r['term'] for r in top_results]
    fold_enrichments = [min(r['fold_enrichment'], 10) for r in top_results]  # Cap at 10 for better visualization
    log10_pvalues = [-math.log10(r['pvalue']) if r['pvalue'] > 0 else 0 for r in top_results]
    gene_counts = [r['genes_in_selection'] for r in top_results]
    
    # Create figure and axis for fold enrichment plot
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Plot fold enrichment as a horizontal bar chart
    bars = ax.barh(terms, fold_enrichments, color='steelblue')
    
    # Add labels for gene counts
    for i, (term, fe, count) in enumerate(zip(terms, fold_enrichments, gene_counts)):
        ax.text(fe + 0.1, i, f"{count} genes", va='center')
    
    # Color the bars according to p-value
    norm = plt.Normalize(0, max(log10_pvalues))
    colors = plt.cm.Reds(norm(log10_pvalues))
    
    for bar, color in zip(bars, colors):
        bar.set_color(color)
    
    # Add p-value information to the bars
    for i, (term, pval) in enumerate(zip(terms, log10_pvalues)):
        if pval > 2:  # Only add text for significant terms (p < 0.01)
            ax.text(0.5, i, f"p={10**-pval:.1e}", va='center', ha='left', color='white', fontweight='bold')
    
    # Create a colorbar legend
    sm = plt.cm.ScalarMappable(cmap=plt.cm.Reds, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label('-log10(p-value)')
    
    # Set title and labels
    ax.set_title('Gene Enrichment Analysis Results', fontsize=16)
    ax.set_xlabel('Fold Enrichment', fontsize=12)
    ax.set_ylabel('Term', fontsize=12)
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close(fig)
    
    # Create a second plot for the p-values
    fig2, ax2 = plt.subplots(figsize=(14, 10))
    ax2.barh(terms, log10_pvalues, color='firebrick')
    ax2.set_title('-log10(p-value) for Enriched Terms', fontsize=16)
    ax2.set_xlabel('-log10(p-value)', fontsize=12)
    ax2.set_ylabel('Term', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_file.replace('.pdf', '_pvalues.pdf'), dpi=300)
    plt.close(fig2)
    
    print(f"Visualizations saved to {output_file} and {output_file.replace('.pdf', '_pvalues.pdf')}")      
    
    # Return to not block execution
    plt.close('all')

# Function to save results to a CSV file
def save_results_to_csv(results: List[Dict], output_file: str):
    """Save enrichment results to a CSV file."""
    if not results:
        print("No results to save.")
        return
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Results saved to {output_file}")
    
    # Count significant terms
    sig_terms = [r for r in results if r.get('adjusted_pvalue', 1) <= 0.05]
    if sig_terms:
        print(f"Found {len(sig_terms)} significantly enriched terms (adj.P.Val ≤ 0.05)")

# Function to export gene lists for each term
def export_gene_lists_by_term(results: List[Dict], output_prefix: str):
    """Export gene lists for each enriched term to separate files."""
    if not results:
        print("No results to export.")
        return
    
    # Create a directory for gene lists if it doesn't exist
    gene_lists_dir = f"{output_prefix}_gene_lists"
    os.makedirs(gene_lists_dir, exist_ok=True)
    
    exported_count = 0
    for result in results:
        # Only export significantly enriched terms
        if result.get('adjusted_pvalue', 1) <= 0.05:
            term = result['term']
            # Create a safe filename
            safe_term = ''.join(c if c.isalnum() else '_' for c in term)
            filename = f"{gene_lists_dir}/{safe_term}.txt"
            
            # Write gene list to file
            with open(filename, 'w') as f:
                f.write(f"# Genes associated with term: {term}\n")
                f.write(f"# P-value: {result['pvalue']}\n")
                f.write(f"# Adjusted P-value: {result['adjusted_pvalue']}\n")
                f.write(f"# Fold Enrichment: {result['fold_enrichment']}\n\n")
                
                genes = result['genes'].split(',')
                for gene in genes:
                    f.write(f"{gene.strip()}\n")
            
            exported_count += 1
    
    if exported_count > 0:
        print(f"Exported gene lists for {exported_count} enriched terms to {gene_lists_dir}/")
    else:
        print("No significantly enriched terms to export gene lists for.")

# Main function
def main():
    args = parse_arguments()
    
    print(f"Loading annotation file: {args.annotation}")
    annotations = parse_annotation_file(args.annotation)
    print(f"Loaded {len(annotations)} annotation entries")
    
    # Extract all unique genes from annotation file
    annotation_genes = set()
    for annotation in annotations:
        if 'locusName' in annotation:
            annotation_genes.add(annotation['locusName'])
    print(f"Found {len(annotation_genes)} unique genes in annotation file")
    
    print(f"Loading significance file: {args.significance}")
    significant_genes = parse_significance_file(
        args.significance, 
        log2fc_cutoff=args.log2fc, 
        pvalue_cutoff=args.pvalue,
        remove_version=args.remove_version
    )
    print(f"Identified {len(significant_genes)} significant genes using adjusted p-value cutoff of {args.pvalue}")
    
    # Check gene overlap between files
    overlap_stats = check_gene_overlap(significant_genes, annotation_genes, args.output)
    
    # If no overlap, exit with warning
    if overlap_stats["overlapping"] == 0:
        print("\n⚠️ WARNING: No overlap between significant genes and annotation genes!")
        print("This suggests a problem with gene ID formatting or incorrect files.")
        print("Please check that gene IDs are in the same format in both files.")
        print("Exiting without performing enrichment analysis.")
        return
    
    # If low overlap, warn but continue if forced or user confirms
    if overlap_stats["overlap_percent"] < 50 and not args.force:
        print(f"\n⚠️ WARNING: Only {overlap_stats['overlap_percent']:.1f}% of significant genes found in annotation file!")
        print("This might indicate a gene ID formatting issue or incorrect files.")
        proceed = input("Do you want to continue with the analysis anyway? (y/n): ")
        if proceed.lower() != 'y':
            print("Exiting at user request.")
            return
    
    print(f"\nPerforming enrichment analysis for {args.type} terms")
    gene_to_terms, term_to_genes, all_genes, significant_genes_in_annotation = prepare_enrichment_data(
        annotations, significant_genes, args.type
    )
    
    print(f"Found {len(significant_genes_in_annotation)} of {len(significant_genes)} significant genes in annotation data")
    print(f"Total background genes: {len(all_genes)}")
    print(f"Total annotation terms: {len(term_to_genes)}")
    
    # Calculate enrichment
    enrichment_results = calculate_enrichment(
        term_to_genes, significant_genes_in_annotation, all_genes, args.pvalue
    )
    print(f"Calculated enrichment for {len(enrichment_results)} terms")
    
    # Apply multiple testing correction
    corrected_results = apply_multiple_testing_correction(enrichment_results, method=args.correction)
    
    # Save results
    output_csv = f"{args.output}_{args.type}_enrichment.csv"
    save_results_to_csv(corrected_results, output_csv)
    
    # Export gene lists by term if requested
    if args.export_gene_lists:
        export_gene_lists_by_term(corrected_results, args.output)
    
    # Visualize if requested
    if args.visualize:
        output_plot = f"{args.output}_{args.type}_enrichment.pdf"
        visualize_enrichment(corrected_results, output_plot, top_n=args.top)
    
    print("Enrichment analysis completed successfully!")

if __name__ == "__main__":
    main()

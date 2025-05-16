#!/usr/bin/env python3
import csv
import argparse
import sys

def clean_gene_id(gene_id):
    """
    Remove quotes and version suffix (.v2.0) from gene IDs
    """
    # Remove quotes if present
    gene_id = gene_id.strip('"')
    
    # Remove version suffix if present
    if ".v" in gene_id:
        gene_id = gene_id.split(".v")[0]
        
    return gene_id

def extract_genes(gene_list_file, data_file, output_file, unmatched_file):
    """
    Extract rows from data_file based on gene IDs in gene_list_file
    """
    # Read gene list
    with open(gene_list_file, 'r') as f:
        gene_list = [line.strip() for line in f if line.strip()]
    
    print(f"Read {len(gene_list)} genes from gene list file")
    
    # Create a set for faster lookups
    gene_set = set(gene_list)
    
    # Track which genes we find matches for
    matched_genes = set()
    
    # Read data file and extract matching rows
    with open(data_file, 'r') as f_in, open(output_file, 'w', newline='') as f_out:
        reader = csv.reader(f_in)
        writer = csv.writer(f_out)
        
        # Write header
        header = next(reader)
        writer.writerow(header)
        
        # Process rows
        rows_written = 0
        total_rows = 0
        
        for row in reader:
            total_rows += 1
            if not row:
                continue
                
            # Extract and clean gene ID from first column
            gene_id = clean_gene_id(row[0])
            
            # Check if this gene is in our list
            if gene_id in gene_set:
                writer.writerow(row)
                matched_genes.add(gene_id)
                rows_written += 1
    
    # Find unmatched genes
    unmatched_genes = gene_set - matched_genes
    
    # Write unmatched genes to file
    with open(unmatched_file, 'w') as f:
        for gene in unmatched_genes:
            f.write(f"{gene}\n")
    
    print(f"Total rows processed: {total_rows}")
    print(f"Matching genes found: {len(matched_genes)}")
    print(f"Unmatched genes: {len(unmatched_genes)}")
    print(f"Wrote {rows_written} rows to {output_file}")
    print(f"Wrote {len(unmatched_genes)} unmatched genes to {unmatched_file}")

def main():
    parser = argparse.ArgumentParser(description='Extract gene data based on gene ID list')
    parser.add_argument('gene_list', help='File containing list of gene IDs, one per line')
    parser.add_argument('data_file', help='CSV file containing gene data')
    parser.add_argument('--output', '-o', default='extracted_genes.csv', 
                        help='Output file for extracted gene data (default: extracted_genes.csv)')
    parser.add_argument('--unmatched', '-u', default='unmatched_genes.txt',
                        help='Output file for unmatched genes (default: unmatched_genes.txt)')
    
    args = parser.parse_args()
    
    try:
        extract_genes(args.gene_list, args.data_file, args.output, args.unmatched)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

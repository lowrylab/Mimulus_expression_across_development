import pandas as pd
import re

def identify_jasmonate_catabolism_genes(annotation_file, output_dir="./", verbose=True):
    """
    Identify genes involved in jasmonate catabolism/inactivation.
    
    Parameters:
    -----------
    annotation_file : str
        Path to the annotation file (TSV format)
    output_dir : str
        Directory to save output files
    verbose : bool
        Print progress information
    """
    # Read the annotation file
    if verbose:
        print(f"Reading annotation file: {annotation_file}")
    
    df = pd.read_csv(annotation_file, sep='\t')
    
    # Define specific keywords for jasmonate catabolism/inactivation
    ja_catabolism_keywords = {
        'CYP94B1': ['CYP94B1', 'cytochrome P450 94B1', 'jasmonoyl-isoleucine 12-hydroxylase'],
        'CYP94B3': ['CYP94B3', 'cytochrome P450 94B3'],
        'CYP94C1': ['CYP94C1', 'cytochrome P450 94C1'],
        'JA-Ile-hydrolase': ['IAR3', 'ILL6', 'amidohydrolase', 'jasmonoyl-isoleucine hydrolase', 'JA-Ile hydrolase'],
        'JA-carboxyl-methyltransferase': ['JMT', 'jasmonic acid carboxyl methyltransferase'],
        'JA-hydroxylases': ['jasmonate hydroxylase', 'jasmonic acid hydroxylase', 'JA hydroxylase'],
        'Sulfotransferase': ['ST2A', 'sulfotransferase', 'hydroxy-jasmonate sulfotransferase'],
        'JA-Ile-oxidase': ['jasmonate oxidase', 'JA-Ile oxidase']
    }
    
    # Dictionary to store identified genes by category
    identified_genes = {category: [] for category in ja_catabolism_keywords}
    
    # CYP94 family is very important, so also look for general CYP94 genes
    general_cyp94_pattern = re.compile(r'\bCYP94[A-Z]?\d*\b', re.IGNORECASE)
    
    # Process the dataframe to find hormone-related genes
    for index, row in df.iterrows():
        if index % 1000 == 0 and verbose:
            print(f"Processing row {index} of {len(df)}")
        
        # Define gene ID
        gene_id = row['locusName']
        if 'transcriptName' in row and pd.notna(row['transcriptName']):
            transcript_id = row['transcriptName']
        else:
            transcript_id = gene_id
            
        # Collect information for reference
        info = {}
        for col in ['arabi-symbol', 'arabi-defline', 'Pfam', 'Panther', 'KOG', 'KEGG/ec', 'KO', 'GO']:
            if col in row and pd.notna(row[col]):
                info[col.replace('-', '_')] = str(row[col])
        
        # Combine all searchable text into a single string
        search_text = ' '.join(str(v).lower() for v in info.values() if pd.notna(v))
        
        # Check for each category of jasmonate catabolism genes
        found_match = False
        for category, keywords in ja_catabolism_keywords.items():
            for keyword in keywords:
                if re.search(r'\b' + re.escape(keyword.lower()) + r'\b', search_text):
                    identified_genes[category].append((transcript_id, keyword, info))
                    found_match = True
                    break
            
            if found_match:
                break
        
        # Also check for general CYP94 pattern for P450s that might not be specifically annotated
        if not found_match:
            for field_value in info.values():
                cyp94_match = general_cyp94_pattern.search(str(field_value))
                if cyp94_match:
                    cyp_id = cyp94_match.group(0)
                    if cyp_id not in [k for k in ja_catabolism_keywords]:  # Avoid double-counting specific CYP94s
                        if 'CYP94_general' not in identified_genes:
                            identified_genes['CYP94_general'] = []
                        identified_genes['CYP94_general'].append((transcript_id, cyp_id, info))
                        break
    
    # Write results to files
    all_catabolism_genes = []
    
    # First, create a detailed annotation file
    with open(f"{output_dir}/jasmonate_catabolism_genes_detailed.tsv", 'w') as f:
        f.write("Category\tGene_ID\tKeyword_Match\tArabidopsis_Symbol\tArabidopsis_Definition\tPfam\tGO\n")
        
        for category, genes in identified_genes.items():
            for gene_id, keyword, info in genes:
                arabi_symbol = info.get('arabi_symbol', '')
                arabi_defline = info.get('arabi_defline', '')
                pfam = info.get('Pfam', '')
                go = info.get('GO', '')
                
                f.write(f"{category}\t{gene_id}\t{keyword}\t{arabi_symbol}\t{arabi_defline}\t{pfam}\t{go}\n")
                all_catabolism_genes.append(gene_id)
    
    # Then create a simple list of all catabolism genes
    with open(f"{output_dir}/jasmonate_catabolism_genes.txt", 'w') as f:
        for gene_id in sorted(set(all_catabolism_genes)):
            f.write(f"{gene_id}\n")
    
    # Count genes by category
    category_counts = {category: len(set([g[0] for g in genes])) for category, genes in identified_genes.items()}
    
    return category_counts

def main():
    # Path to your annotation file
    annotation_file = "Mguttatus_256_v2.0.annotation_info.txt"
    
    # Run the analysis
    stats = identify_jasmonate_catabolism_genes(annotation_file)
    
    print("\nSummary of jasmonate catabolism genes:")
    print("=====================================")
    total_genes = sum(stats.values())
    for category, count in stats.items():
        print(f"{category}: {count} genes")
    print(f"Total: {total_genes} genes")

if __name__ == "__main__":
    main()

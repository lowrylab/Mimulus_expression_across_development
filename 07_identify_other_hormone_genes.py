import pandas as pd
import re

def identify_hormone_pathway_genes(annotation_file, output_dir="./", verbose=True):
    """
    Identify genes involved in major plant hormone pathways.
    
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
    
    # Dictionary to store hormone pathway genes
    hormone_genes = {
        'auxin': [],
        'cytokinin': [],
        'ethylene': [],
        'abscisic_acid': [],
        'brassinosteroid': [],
        'salicylic_acid': [],
        'strigolactone': []
    }
    
    # Define keywords for each hormone pathway
    hormone_keywords = {
        'auxin': [
            'auxin', 'AUX', 'IAA', 'YUCCA', 'YUC', 'TAA', 'PIN', 'ABCB', 'TIR1', 'AFB', 
            'AUX/IAA', 'ARF', 'SAUR', 'GH3', 'auxin efflux', 'auxin influx', 'auxin response',
            'PILS', 'LAX', 'AUX1'
        ],
        'cytokinin': [
            'cytokinin', 'CK', 'IPT', 'isopentenyltransferase', 'CYP735A', 'LOG', 
            'cytokinin oxidase', 'CKX', 'AHK', 'histidine kinase', 'AHP', 'ARR', 'response regulator',
            'histidine phosphotransfer', 'CRF'
        ],
        'ethylene': [
            'ethylene', 'ACC synthase', 'ACS', 'ACC oxidase', 'ACO', 'ETR', 'ERS', 
            'ethylene receptor', 'CTR1', 'EIN2', 'EIN3', 'EIL', 'ERF', 'EBF',
            'ethylene response', 'ethylene insensitive', 'ethylene resistant'
        ],
        'abscisic_acid': [
            'abscisic acid', 'ABA', 'NCED', 'carotenoid cleavage', 'ABA aldehyde oxidase', 
            'AAO', 'ABA receptor', 'PYL', 'PYR', 'RCAR', 'PP2C', 'protein phosphatase 2C', 
            'SnRK', 'ABRE', 'ABF', 'ABI', 'abscisic acid insensitive'
        ],
        'brassinosteroid': [
            'brassinosteroid', 'BRI1', 'BAK1', 'BIN2', 'BSU1', 'BZR', 'BES1', 
            'brassinosteroid insensitive', 'BRZ', 'DWF4', 'ROT3', 'CPD', 'brassinolide',
            'CYP90', 'CYP85A', 'BRL', 'SERK'
        ],
        'salicylic_acid': [
            'salicylic acid', 'SA', 'ICS', 'isochorismate synthase', 'PAL', 'SAGT', 
            'NPR1', 'NPR3', 'NPR4', 'TGA', 'PR1', 'pathogenesis related', 'SABP',
            'EDS', 'PAD', 'SAMT', 'salicylate'
        ],
        'strigolactone': [
            'strigolactone', 'SL', 'CCD', 'carotenoid cleavage dioxygenase', 'CCD7', 'CCD8', 
            'D14', 'D3', 'MAX', 'MORE AXILLARY GROWTH', 'D53', 'SMXL', 'DWARF', 'carlactone'
        ]
    }
    
    # Process the dataframe to find hormone-related genes
    for index, row in df.iterrows():
        if index % 1000 == 0 and verbose:
            print(f"Processing row {index} of {len(df)}")
        
        # Columns to search in
        search_fields = []
        for col in ['arabi-symbol', 'arabi-defline', 'GO', 'Pfam', 'Panther', 'KOG', 'KEGG/ec', 'KO']:
            if col in row and isinstance(row[col], str):
                search_fields.append(row[col].lower())
            
        # Combine all searchable text
        search_text = ' '.join(search_fields)
        
        # Check for each hormone pathway
        for hormone, keywords in hormone_keywords.items():
            for keyword in keywords:
                if re.search(r'\b' + keyword.lower() + r'\b', search_text):
                    gene_id = row['locusName']
                    if 'transcriptName' in row:
                        gene_id = row['transcriptName']
                    
                    # Collect additional information for reference
                    info = {}
                    if 'arabi-symbol' in row and isinstance(row['arabi-symbol'], str):
                        info['arabi_symbol'] = row['arabi-symbol']
                    if 'arabi-defline' in row and isinstance(row['arabi-defline'], str):
                        info['arabi_defline'] = row['arabi-defline']
                    
                    hormone_genes[hormone].append((gene_id, keyword, info))
                    break  # Once a gene is assigned to a pathway, move to the next gene
    
    # Write results to files
    for hormone, genes in hormone_genes.items():
        if genes:
            # Write a simple list file with just gene IDs
            with open(f"{output_dir}/{hormone}_genes.txt", 'w') as f:
                seen_genes = set()
                for gene_id, _, _ in genes:
                    if gene_id not in seen_genes:
                        f.write(f"{gene_id}\n")
                        seen_genes.add(gene_id)
            
            # Write a detailed file with annotations
            with open(f"{output_dir}/{hormone}_genes_annotated.tsv", 'w') as f:
                f.write("Gene_ID\tMatchedKeyword\tArabidopsisSymbol\tArabidopsisDefinition\n")
                seen_genes = set()
                for gene_id, keyword, info in genes:
                    if gene_id not in seen_genes:
                        arabi_symbol = info.get('arabi_symbol', '')
                        arabi_defline = info.get('arabi_defline', '')
                        f.write(f"{gene_id}\t{keyword}\t{arabi_symbol}\t{arabi_defline}\n")
                        seen_genes.add(gene_id)
            
            if verbose:
                print(f"Found {len(set([g[0] for g in genes]))} genes for {hormone} pathway")
    
    # Return statistics
    stats = {hormone: len(set([g[0] for g in genes])) for hormone, genes in hormone_genes.items()}
    return stats

def main():
    # Path to your annotation file
    annotation_file = "Mguttatus_256_v2.0.annotation_info.txt"
    
    # Run the analysis
    stats = identify_hormone_pathway_genes(annotation_file)
    
    print("\nSummary of hormone pathway genes:")
    print("=================================")
    for hormone, count in stats.items():
        print(f"{hormone.replace('_', ' ').title()}: {count} genes")

if __name__ == "__main__":
    main()

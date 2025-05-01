import pandas as pd
import sys

# Load your annotation file with the correct filename
try:
    df = pd.read_csv('Mguttatus_256_v2.0.annotation_info.txt', sep='\t')
    print(f"Successfully loaded file with {len(df)} rows")
    print(f"Columns in the file: {', '.join(df.columns)}")
except Exception as e:
    print(f"Error loading file: {e}")
    sys.exit(1)

# 1. Find genes with relevant GO terms for salt stress and SOS pathway
salt_stress_go_terms = [
    'GO:0042538',  # Hyperosmotic salinity response
    'GO:0009651',  # Response to salt stress
    'GO:0006970',  # Response to osmotic stress
    'GO:0010359',  # Regulation of anion channel activity
    'GO:0006814',  # Sodium ion transport
    'GO:0006885',  # Regulation of pH
    'GO:0005516'   # Calmodulin binding
]

# Safer function to check if a term exists in a string
def safe_contains(value, search_terms):
    if pd.isna(value):
        return False
    value_str = str(value).lower()
    for term in search_terms:
        if str(term).lower() in value_str:
            return True
    return False

# Make sure 'GO' column exists
if 'GO' in df.columns:
    salt_genes_go = df[df['GO'].apply(lambda x: safe_contains(x, salt_stress_go_terms))]
    print(f"Found {len(salt_genes_go)} genes with salt stress GO terms")
else:
    print("Warning: 'GO' column not found. Skipping GO term filtering.")
    salt_genes_go = pd.DataFrame(columns=df.columns)

# 2. Find genes with relevant keywords in Arabidopsis descriptions
sos_keywords = [
    'salt overly sensitive', 'SOS1', 'SOS2', 'SOS3', 'SOS4', 'SOS5',
    'CBL10', 'CBL4', 'CIPK24', 'CIPK8', 'SCABP8',
    'Na+/H+ antiporter', 'NHX1', 'calcium sensor', 'calcineurin B-like',
    'salt tolerance', 'sodium transport', 'sodium/proton exchanger',
    'AtANN4', 'annexin', 'calcium-dependent', 'V-ATPase', 'H+-PPase', 'AVP1',
    'PKS5', 'ABI2', 'GIGANTEA', '14-3-3', 'AKT1'
]

# Make sure 'arabi-defline' column exists
if 'arabi-defline' in df.columns:
    sos_genes_desc = df[df['arabi-defline'].apply(lambda x: safe_contains(x, sos_keywords))]
    print(f"Found {len(sos_genes_desc)} genes with SOS-related keywords in descriptions")
else:
    print("Warning: 'arabi-defline' column not found. Skipping description filtering.")
    sos_genes_desc = pd.DataFrame(columns=df.columns)

# 3. Find genes with relevant Arabidopsis symbols
if 'arabi-symbol' in df.columns:
    sos_genes_symbol = df[df['arabi-symbol'].apply(lambda x: safe_contains(x, sos_keywords))]
    print(f"Found {len(sos_genes_symbol)} genes with SOS-related keywords in symbols")
else:
    print("Warning: 'arabi-symbol' column not found. Skipping symbol filtering.")
    sos_genes_symbol = pd.DataFrame(columns=df.columns)

# 4. Find genes with Pfam domains associated with SOS pathway components
sos_pfam_domains = [
    'PF00069',  # Protein kinase domain (SOS2/CIPK24)
    'PF07887',  # Calcineurin-like phosphoesterase (SOS3/CBL4)
    'PF03020',  # EF-hand calcium binding domains (SOS3/CBL4, CBL10)
    'PF00999',  # Sodium/hydrogen exchanger family (SOS1)
    'PF00571',  # CBS domain (regulatory domain in many transporters)
    'PF00213',  # ATP-dependent proton pump
    'PF00231'   # Vacuolar ATPase
]

if 'Pfam' in df.columns:
    sos_genes_pfam = df[df['Pfam'].apply(lambda x: safe_contains(x, sos_pfam_domains))]
    print(f"Found {len(sos_genes_pfam)} genes with SOS-related Pfam domains")
else:
    print("Warning: 'Pfam' column not found. Skipping Pfam domain filtering.")
    sos_genes_pfam = pd.DataFrame(columns=df.columns)

# 5. Combine all results and remove duplicates
sos_genes = pd.concat([salt_genes_go, sos_genes_desc, sos_genes_symbol, sos_genes_pfam]).drop_duplicates()
print(f"Found {len(sos_genes)} unique genes related to SOS pathway after combining all sources")

# 6. Add SOS pathway component classification
def categorize_sos_gene(row):
    # Default to other salt stress response
    category = 'Other salt stress response'
    
    # Check if necessary columns exist
    desc = str(row.get('arabi-defline', '')) if 'arabi-defline' in row else ''
    symbol = str(row.get('arabi-symbol', '')) if 'arabi-symbol' in row else ''
    
    desc = desc.lower()
    symbol = symbol.lower()
    
    # Categorize based on keywords in description or symbol
    if any(x in desc or x in symbol for x in ['sos1', 'sodium/hydrogen exchanger', 'na+/h+ antiporter']):
        category = 'SOS1 (Na+/H+ antiporter)'
    elif any(x in desc or x in symbol for x in ['sos2', 'cipk24', 'protein kinase']):
        category = 'SOS2/CIPK24 (protein kinase)'
    elif any(x in desc or x in symbol for x in ['sos3', 'cbl4', 'calcium sensor']):
        category = 'SOS3/CBL4 (calcium sensor - root)'
    elif any(x in desc or x in symbol for x in ['cbl10', 'scabp8']):
        category = 'CBL10/SCABP8 (calcium sensor - shoot)'
    elif any(x in desc or x in symbol for x in ['cipk8']):
        category = 'CIPK8 (protein kinase)'
    elif any(x in desc or x in symbol for x in ['nhx1', 'vacuolar', 'na+/h+ exchanger']):
        category = 'NHX1 (vacuolar Na+/H+ exchanger)'
    elif any(x in desc or x in symbol for x in ['avp1', 'h+-ppase', 'pyrophosphatase']):
        category = 'AVP1 (H+-PPase)'
    elif any(x in desc or x in symbol for x in ['v-atpase', 'vacuolar atpase']):
        category = 'V-ATPase'
    elif any(x in desc or x in symbol for x in ['ann4', 'annexin']):
        category = 'AtANN4 (annexin)'
    elif any(x in desc or x in symbol for x in ['akt1', 'potassium channel']):
        category = 'AKT1 (K+ channel)'
    elif any(x in desc or x in symbol for x in ['14-3-3']):
        category = '14-3-3 (inhibitor)'
    elif any(x in desc or x in symbol for x in ['gigantea', 'gi']):
        category = 'GIGANTEA (inhibitor)'
    elif any(x in desc or x in symbol for x in ['abi2']):
        category = 'ABI2 (inhibitor)'
    
    return category

# Add category column
sos_genes['SOS_pathway_component'] = sos_genes.apply(categorize_sos_gene, axis=1)

# 7. Add a column for the evidence/matching criteria
def identify_evidence(row):
    evidence = []
    
    # Check GO terms if column exists
    if 'GO' in row and not pd.isna(row['GO']):
        go_terms = str(row['GO'])
        if any(term in go_terms for term in salt_stress_go_terms):
            evidence.append('GO term match')
    
    # Check description if column exists
    if 'arabi-defline' in row and not pd.isna(row['arabi-defline']):
        desc = str(row['arabi-defline']).lower()
        if any(keyword.lower() in desc for keyword in sos_keywords):
            evidence.append('Description keyword match')
    
    # Check symbol if column exists
    if 'arabi-symbol' in row and not pd.isna(row['arabi-symbol']):
        symbol = str(row['arabi-symbol']).lower()
        if any(keyword.lower() in symbol for keyword in sos_keywords):
            evidence.append('Symbol keyword match')
    
    # Check Pfam if column exists
    if 'Pfam' in row and not pd.isna(row['Pfam']):
        pfam = str(row['Pfam'])
        if any(domain in pfam for domain in sos_pfam_domains):
            evidence.append('Pfam domain match')
    
    return ', '.join(evidence) if evidence else 'Unknown'

sos_genes['Evidence'] = sos_genes.apply(identify_evidence, axis=1)

# 8. Create a single comprehensive output file
# Keep all columns from the original dataframe
sos_genes.to_csv('sos_pathway_and_salt_stress_genes.tsv', sep='\t', index=False)

# 9. Create a summary file with key information
# Use only columns that exist in the dataframe
available_summary_columns = ['SOS_pathway_component', 'Evidence']
for col in df.columns:
    available_summary_columns.append(col)

sos_genes_summary = sos_genes[available_summary_columns]
sos_genes_summary.to_csv('sos_pathway_genes_summary.tsv', sep='\t', index=False)

# Print summary
print(f"\nFound {len(sos_genes)} genes related to the SOS pathway and salt stress response.")
print(f"Complete results saved to 'sos_pathway_and_salt_stress_genes.tsv'")
print(f"Summary results saved to 'sos_pathway_genes_summary.tsv'")

# Print breakdown by category
category_counts = sos_genes['SOS_pathway_component'].value_counts()
print("\nBreakdown by SOS pathway component:")
for category, count in category_counts.items():
    print(f"  {category}: {count} genes")

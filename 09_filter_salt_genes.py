import pandas as pd

# Load your salt genes file
df = pd.read_csv('sos_pathway_and_salt_stress_genes.tsv', sep='\t')

# Create a list of key SOS pathway components to prioritize
key_sos_components = [
    'SOS1 (Na+/H+ antiporter)',
    'SOS2/CIPK24 (protein kinase)',
    'SOS3/CBL4 (calcium sensor - root)',
    'CBL10/SCABP8 (calcium sensor - shoot)',
    'CIPK8 (protein kinase)',
    'NHX1 (vacuolar Na+/H+ exchanger)',
    'AVP1 (H+-PPase)',
    'V-ATPase',
    'AKT1 (K+ channel)',
    'AtANN4 (annexin)',
    '14-3-3 (inhibitor)',
    'GIGANTEA (inhibitor)',
    'ABI2 (inhibitor)'
]

# Filter for high-confidence genes (specific SOS components and/or multiple evidence)
high_confidence_genes = df[(df['SOS_pathway_component'].isin(key_sos_components)) | 
                          (df['Evidence'].str.contains(','))]

# Create a count for evidence types as a confidence score
high_confidence_genes['Evidence_count'] = high_confidence_genes['Evidence'].apply(lambda x: len(x.split(',')))

# Add a column to indicate key Arabidopsis gene symbols related to SOS pathway
sos_keywords = ['SOS1', 'SOS2', 'SOS3', 'SOS4', 'SOS5', 'CBL10', 'CBL4', 'CIPK24', 
                'CIPK8', 'SCABP8', 'NHX1', 'AVP1', 'AKT1', 'AtANN4', 'GI', 'ABI2']

# Check if any SOS-related keyword is in the arabi-symbol or arabi-defline
high_confidence_genes['Direct_SOS_homolog'] = high_confidence_genes.apply(
    lambda row: any(kw.lower() in str(row.get('arabi-symbol', '')).lower() or 
                    kw.lower() in str(row.get('arabi-defline', '')).lower() 
                    for kw in sos_keywords), 
    axis=1
)

# Sort by SOS component type, Direct SOS homolog status, and evidence count
high_confidence_genes = high_confidence_genes.sort_values(
    by=['SOS_pathway_component', 'Direct_SOS_homolog', 'Evidence_count'], 
    ascending=[True, False, False]
)

# Save the high-confidence genes to a new file
high_confidence_genes.to_csv('high_confidence_sos_pathway_genes.tsv', sep='\t', index=False)

# Also create a summary with key columns for quick reference
summary_columns = ['locusName', 'arabi-symbol', 'arabi-defline', 'SOS_pathway_component', 
                  'Evidence', 'Evidence_count', 'Direct_SOS_homolog']
summary_df = high_confidence_genes[summary_columns]
summary_df.to_csv('high_confidence_sos_pathway_genes_summary.tsv', sep='\t', index=False)

# Print summary statistics
print(f"Total salt stress genes: {len(df)}")
print(f"High-confidence SOS pathway genes: {len(high_confidence_genes)}")
print("\nBreakdown by SOS pathway component:")
component_counts = high_confidence_genes['SOS_pathway_component'].value_counts()
for component, count in component_counts.items():
    print(f"  {component}: {count} genes")

print("\nTop 10 highest confidence genes:")
for i, row in high_confidence_genes.sort_values(
    by=['Direct_SOS_homolog', 'Evidence_count'], ascending=[False, False]).head(10).iterrows():
    print(f"  {row['locusName']} - {row.get('arabi-symbol', 'No symbol')} - "
          f"{row['SOS_pathway_component']} - Evidence: {row['Evidence_count']}")

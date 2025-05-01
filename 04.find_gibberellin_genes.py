import pandas as pd

# Load your annotation file with the correct filename
df = pd.read_csv('Mguttatus_256_v2.0.annotation_info.txt', sep='\t')

# 1. Find genes with relevant GO terms
ga_go_terms = ['GO:0009686', 'GO:0009740']  # Gibberellin 
ja_go_terms = ['GO:0009695', 'GO:0009867']  # Jasmonate

ga_genes_go = df[df['GO'].apply(lambda x: any(term in str(x) for term in ga_go_terms))]
ja_genes_go = df[df['GO'].apply(lambda x: any(term in str(x) for term in ja_go_terms))]

# 2. Find genes with relevant keywords in Arabidopsis descriptions
ga_keywords = ['gibberellin', 'DELLA', 'GA20', 'GA3', 'GA2ox', 'GID1']
ja_keywords = ['jasmonate', 'jasmonic', 'LOX', 'AOS', 'AOC', 'OPR', 'JAZ', 'COI1', 'MYC']

ga_genes_desc = df[df['arabi-defline'].apply(lambda x: any(keyword.lower() in str(x).lower() for keyword in ga_keywords))]
ja_genes_desc = df[df['arabi-defline'].apply(lambda x: any(keyword.lower() in str(x).lower() for keyword in ja_keywords))]

# 3. Find genes with relevant Arabidopsis symbols
ga_genes_symbol = df[df['arabi-symbol'].apply(lambda x: any(keyword.lower() in str(x).lower() for keyword in ga_keywords))]
ja_genes_symbol = df[df['arabi-symbol'].apply(lambda x: any(keyword.lower() in str(x).lower() for keyword in ja_keywords))]

# 4. Combine all results and remove duplicates
ga_genes = pd.concat([ga_genes_go, ga_genes_desc, ga_genes_symbol]).drop_duplicates()
ja_genes = pd.concat([ja_genes_go, ja_genes_desc, ja_genes_symbol]).drop_duplicates()

# Output results
ga_genes.to_csv('gibberellin_related_genes.tsv', sep='\t', index=False)
ja_genes.to_csv('jasmonate_related_genes.tsv', sep='\t', index=False)

# Print summary
print(f"Found {len(ga_genes)} gibberellin-related genes and {len(ja_genes)} jasmonate-related genes.")
print("Results saved to 'gibberellin_related_genes.tsv' and 'jasmonate_related_genes.tsv'")

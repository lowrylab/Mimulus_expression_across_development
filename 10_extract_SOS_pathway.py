#!/usr/bin/env python3
"""
extract_SOS_pathway.py - Extract SOS pathway and related salt stress response genes from a plant annotation file

This script:
1. Reads a tab-delimited annotation file
2. Identifies genes involved in the SOS (Salt Overly Sensitive) pathway and related salt stress responses
3. Outputs a filtered list with their annotations

Usage:
python extract_SOS_pathway.py -i <input_annotation_file> -o <output_file> [-v]

Example annotation file format:
pacId	locusName	transcriptName	peptideName	Pfam	Panther	KOG	KEGG/ec	KO	GO	Best-hit-arabi-name	arabi-symbol	arabi-defline
"""

import argparse
import re
import csv
import sys

# Define the reference SOS pathway and related components with their NCBI accessions and Arabidopsis gene IDs
SOS_PATHWAY_REFS = {
    # SOS1 - Na+/H+ antiporter
    'NP_178307': {'name': 'SOS1', 'gene_id': 'AT2G01980', 'alt_names': ['NHX7']},
    
    # SOS2 - Serine/threonine protein kinase
    'NP_172390': {'name': 'SOS2', 'gene_id': 'AT5G35410', 'alt_names': ['CIPK24', 'PKS18', 'SNRK3.24']},
    
    # SOS3 - Calcium sensor
    'NP_198856': {'name': 'SOS3', 'gene_id': 'AT5G24270', 'alt_names': ['CBL4']},
    
    # CBL10 - Calcineurin B-like protein 10
    'NP_174060': {'name': 'CBL10', 'gene_id': 'AT4G33000', 'alt_names': ['SCABP8', 'CALCINEURIN B-LIKE 10']},
    
    # NHX1 - Na+/H+ exchanger 1
    'NP_198067': {'name': 'NHX1', 'gene_id': 'AT5G27150', 'alt_names': ['NA+/H+ EXCHANGER 1']},
    
    # HKT1 - High-affinity K+ transporter 1
    'NP_567354': {'name': 'HKT1', 'gene_id': 'AT4G10310', 'alt_names': ['ATHKT1', 'HKT1;1']},
    
    # CIPK8 - CBL-interacting protein kinase 8
    'NP_172223': {'name': 'CIPK8', 'gene_id': 'AT4G24400', 'alt_names': ['PKS11', 'SNRK3.13']}
}

# Compile all gene IDs, names and alternate names for easy searching
SOS_PATHWAY_NAMES = set()
SOS_PATHWAY_GENE_IDS = set()

for info in SOS_PATHWAY_REFS.values():
    SOS_PATHWAY_NAMES.add(info['name'])
    SOS_PATHWAY_GENE_IDS.add(info['gene_id'])
    for alt_name in info['alt_names']:
        SOS_PATHWAY_NAMES.add(alt_name)

# Define domains and GO terms associated with SOS pathway proteins
SOS_PATHWAY_DOMAINS = {
    # SOS1/NHX1 domains
    'PF00999': 'Na_H_Exchanger - Sodium/hydrogen exchanger family (SOS1, NHX1)',
    'PF07569': 'Na_H_Exchanger_C - Na+/H+ exchanger C-terminus (SOS1)',
    
    # SOS2/CIPK8 domains
    'PF00069': 'Pkinase - Protein kinase domain (SOS2, CIPK8)',
    'PF03822': 'NAF - NAF domain (SOS2, CIPK8)',
    'PF07714': 'Pkinase_Tyr - Protein tyrosine kinase (SOS2, CIPK8)',
    
    # SOS3/CBL10 domains
    'PF07887': 'EF-hand_10 - Calcineurin B-like protein (SOS3, CBL10)',
    'PF13499': 'EF-hand_7 - EF-hand domain pair (SOS3, CBL10)',
    'PF13202': 'EF-hand_5 - EF hand domain (SOS3, CBL10)',
    'PF13405': 'EF-hand_6 - EF-hand domain (SOS3, CBL10)',
    
    # HKT1 domains
    'PF02386': 'TrkH - Cation transport protein (HKT1)',
    'PF02705': 'K_trans - K+ potassium transporter (HKT1)',
    'PF07885': 'Ion_trans_2 - Ion channel (HKT1)'
}

SOS_PATHWAY_GO_TERMS = {
    # General salt stress and ion transport terms
    'GO:0006814': 'sodium ion transport',
    'GO:0006885': 'regulation of pH',
    'GO:0009651': 'response to salt stress',
    'GO:0042542': 'response to hydrogen peroxide',
    'GO:0016036': 'cellular response to phosphate starvation',
    'GO:0042538': 'hyperosmotic salinity response',
    
    # SOS1/NHX1 related
    'GO:0015385': 'sodium:proton antiporter activity',
    'GO:0035725': 'sodium ion transmembrane transport',
    
    # HKT1 related
    'GO:0006813': 'potassium ion transport',
    'GO:0071805': 'potassium ion transmembrane transport',
    'GO:0015079': 'potassium ion transmembrane transporter activity',
    'GO:0022890': 'inorganic cation transmembrane transporter activity',
    
    # SOS2/CIPK8 related
    'GO:0004672': 'protein kinase activity',
    'GO:0006468': 'protein phosphorylation',
    
    # SOS3/CBL10 related
    'GO:0005516': 'calmodulin binding',
    'GO:0005509': 'calcium ion binding',
    
    # Cellular locations
    'GO:0005886': 'plasma membrane',
    'GO:0009705': 'plant-type vacuole membrane',
    'GO:0005774': 'vacuolar membrane'
}

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Extract SOS pathway and related salt stress response genes from annotation file")
    parser.add_argument("-i", "--input", required=True, help="Input annotation file (tab-delimited)")
    parser.add_argument("-o", "--output", required=True, help="Output file for SOS pathway genes")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print detailed progress information")
    parser.add_argument("-b", "--blast", help="Optional BLASTP results file (tabular format) for additional evidence")
    return parser.parse_args()

def read_blast_results(blast_file):
    """
    Read BLAST results in tabular format if available
    Returns a dictionary of query IDs that match SOS pathway components
    """
    if not blast_file:
        return None
    
    blast_hits = {}
    try:
        with open(blast_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                
                query_id = parts[0]
                subject_id = parts[1]
                
                # Check if the subject ID matches any of our reference accessions
                for accession in SOS_PATHWAY_REFS.keys():
                    if accession in subject_id:
                        blast_hits[query_id] = {
                            'subject_id': subject_id,
                            'component': SOS_PATHWAY_REFS[accession]['name']
                        }
                        break
        
        return blast_hits
    except Exception as e:
        print(f"Warning: Could not read BLAST results file: {e}")
        return None

def is_sos_pathway_gene(row, blast_results=None):
    """
    Determine if a gene annotation row potentially represents a SOS pathway or related component
    
    Criteria:
    1. Matches NCBI accessions of reference components
    2. Arabidopsis best hit matches known genes
    3. Symbol or defline contains known keywords
    4. Has domains associated with these proteins
    5. Has GO terms associated with these functions
    """
    # Check if there's a reference to any of the NCBI accessions in any field
    for field in row.values():
        if field:
            for accession in SOS_PATHWAY_REFS.keys():
                if accession in field:
                    return True
    
    # Check Arabidopsis gene IDs
    if row['Best-hit-arabi-name']:
        for gene_id in SOS_PATHWAY_GENE_IDS:
            if gene_id in row['Best-hit-arabi-name']:
                return True
    
    # Check if the symbol contains any pathway component names
    if row['arabi-symbol']:
        for name in SOS_PATHWAY_NAMES:
            if re.search(r'\b' + re.escape(name) + r'\b', row['arabi-symbol'], re.IGNORECASE):
                return True
    
    # Check if the defline contains any pathway component names
    if row['arabi-defline']:
        for name in SOS_PATHWAY_NAMES:
            if re.search(r'\b' + re.escape(name) + r'\b', row['arabi-defline'], re.IGNORECASE):
                return True
    
    # Check for specific terms in defline
    if row['arabi-defline'] and any(re.search(term, row['arabi-defline'], re.IGNORECASE) 
                                  for term in ['salt overly sensitive', 'sodium/hydrogen exchanger',
                                              'Na\+/H\+ antiporter', 'Na\+/H\+ exchanger',
                                              'calcium sensor', 'calcineurin B-like',
                                              'CBL-interacting protein kinase',
                                              'high-affinity potassium transporter']):
        return True
    
    # Check Pfam domains for those associated with pathway components
    if row['Pfam']:
        for domain in SOS_PATHWAY_DOMAINS.keys():
            if domain in row['Pfam']:
                return True
    
    # Check GO terms for relevant functions
    if row['GO']:
        for go_term in SOS_PATHWAY_GO_TERMS.keys():
            if go_term in row['GO']:
                return True
    
    # Check BLAST results if provided
    if blast_results and row['peptideName'] in blast_results:
        return True
    
    return False

def identify_component(row, blast_results=None):
    """
    Identify which specific component this gene represents
    Returns (component, confidence, reasoning)
    """
    component = "Unknown"
    confidence = "Low"
    reasons = []
    score = 0
    
    # Check for direct matches to NCBI accessions (highest confidence)
    for field in row.values():
        if field:
            for accession, info in SOS_PATHWAY_REFS.items():
                if accession in field:
                    component = info['name']
                    score += 5
                    reasons.append(f"Direct match to {info['name']} ({accession})")
                    break
    
    # Check Arabidopsis gene ID
    if row['Best-hit-arabi-name']:
        for accession, info in SOS_PATHWAY_REFS.items():
            if info['gene_id'] in row['Best-hit-arabi-name']:
                component = info['name']
                score += 4
                reasons.append(f"Match to Arabidopsis {info['name']} ({info['gene_id']})")
                break
    
    # Check symbol and defline for component names
    for field in [row['arabi-symbol'], row['arabi-defline']]:
        if field:
            for accession, info in SOS_PATHWAY_REFS.items():
                # Check main name
                if re.search(r'\b' + re.escape(info['name']) + r'\b', field, re.IGNORECASE):
                    component = info['name']
                    score += 3
                    reasons.append(f"Name match to {info['name']}")
                    break
                
                # Check alternative names
                for alt_name in info['alt_names']:
                    if re.search(r'\b' + re.escape(alt_name) + r'\b', field, re.IGNORECASE):
                        component = info['name']
                        score += 3
                        reasons.append(f"Alternative name match to {alt_name} ({info['name']})")
                        break
    
    # Check domain-based evidence
    if row['Pfam']:
        # SOS1/NHX1 specific domains
        if 'PF00999' in row['Pfam']:  # Na_H_Exchanger
            if component == "Unknown":
                if 'PF07569' in row['Pfam']:  # Na_H_Exchanger_C - more likely SOS1
                    component = "SOS1"
                else:
                    component = "NHX1/SOS1"  # Could be either
            score += 2
            reasons.append("Contains Na+/H+ exchanger domain (PF00999)")
        
        # SOS2/CIPK8 specific domains
        if 'PF00069' in row['Pfam'] and 'PF03822' in row['Pfam']:  # Protein kinase + NAF
            if component == "Unknown":
                component = "SOS2/CIPK8"  # Could be either
            score += 2
            reasons.append("Contains protein kinase (PF00069) and NAF (PF03822) domains")
        
        # SOS3/CBL10 specific domains
        ef_hand_domains = ['PF07887', 'PF13499', 'PF13202', 'PF13405']
        ef_hand_count = sum(1 for domain in ef_hand_domains if domain in row['Pfam'])
        if ef_hand_count > 0:
            if component == "Unknown":
                component = "SOS3/CBL10"  # Could be either calcium sensor
            score += ef_hand_count
            reasons.append(f"Contains {ef_hand_count} EF-hand domains associated with calcium binding")
        
        # HKT1 specific domains
        hkt_domains = ['PF02386', 'PF02705', 'PF07885']
        hkt_count = sum(1 for domain in hkt_domains if domain in row['Pfam'])
        if hkt_count > 0:
            if component == "Unknown":
                component = "HKT1"
            score += hkt_count
            reasons.append(f"Contains {hkt_count} domains associated with cation transporters")
    
    # Check GO terms for specific evidence
    if row['GO']:
        # SOS1/NHX1 specific GO terms
        if 'GO:0015385' in row['GO']:  # sodium:proton antiporter activity
            if component == "Unknown" or component == "NHX1/SOS1":
                component = "SOS1/NHX1"  # Could be either
            score += 2
            reasons.append("Has sodium:proton antiporter activity (GO:0015385)")
        
        # HKT1 specific GO terms
        if 'GO:0015079' in row['GO']:  # potassium ion transmembrane transporter activity
            if component == "Unknown":
                component = "HKT1"
            score += 2
            reasons.append("Has potassium ion transmembrane transporter activity (GO:0015079)")
        
        # SOS2/CIPK8 specific terms
        if 'GO:0004672' in row['GO'] and 'GO:0006468' in row['GO']:  # protein kinase activity + phosphorylation
            if component == "Unknown" or component == "SOS2/CIPK8":
                component = "SOS2/CIPK8"
            score += 2
            reasons.append("Has protein kinase activity and involved in phosphorylation")
        
        # SOS3/CBL10 specific GO terms
        if 'GO:0005509' in row['GO']:  # calcium ion binding
            if component == "Unknown" or component == "SOS3/CBL10":
                component = "SOS3/CBL10"
            score += 2
            reasons.append("Has calcium ion binding activity (GO:0005509)")
        
        # General salt stress terms
        salt_stress_terms = ['GO:0009651', 'GO:0042538']
        salt_term_count = sum(1 for term in salt_stress_terms if term in row['GO'])
        if salt_term_count > 0:
            score += salt_term_count
            reasons.append(f"Associated with {salt_term_count} salt stress GO terms")
    
    # Check BLAST results if provided
    if blast_results and row['peptideName'] in blast_results:
        blast_info = blast_results[row['peptideName']]
        component = blast_info['component']
        score += 4
        reasons.append(f"BLAST hit to {component} ({blast_info['subject_id']})")
    
    # Try to resolve ambiguous components
    if component in ["SOS2/CIPK8"]:
        # Try to distinguish SOS2 from CIPK8 by additional evidence
        sos2_indicators = ['CIPK24', 'PKS18', 'AT5G35410']
        cipk8_indicators = ['CIPK8', 'PKS11', 'AT4G24400']
        
        sos2_count = sum(1 for indicator in sos2_indicators if any(
            indicator in field for field in [row['arabi-symbol'], row['arabi-defline'], row['Best-hit-arabi-name']] if field))
        
        cipk8_count = sum(1 for indicator in cipk8_indicators if any(
            indicator in field for field in [row['arabi-symbol'], row['arabi-defline'], row['Best-hit-arabi-name']] if field))
        
        if sos2_count > cipk8_count:
            component = "SOS2"
            reasons.append(f"More evidence for SOS2 than CIPK8 ({sos2_count} vs {cipk8_count})")
        elif cipk8_count > sos2_count:
            component = "CIPK8"
            reasons.append(f"More evidence for CIPK8 than SOS2 ({cipk8_count} vs {sos2_count})")
    
    # Assign confidence level
    if score >= 5:
        confidence = "High"
    elif score >= 3:
        confidence = "Medium"
    else:
        confidence = "Low"
    
    return (component, confidence, "; ".join(reasons))

def read_annotation_file(file_path):
    """Read the tab-delimited annotation file and return a list of dictionaries"""
    try:
        with open(file_path, 'r') as f:
            # Detect if file has header or not
            first_line = f.readline().strip()
            f.seek(0)  # Reset file pointer
            
            if first_line.startswith('#'):
                # File has header starting with #
                reader = csv.DictReader(f, delimiter='\t')
                header = next(reader)  # Skip the header line
                return list(reader)
            else:
                # File doesn't have a header, use default column names
                column_names = [
                    'pacId', 'locusName', 'transcriptName', 'peptideName', 'Pfam', 
                    'Panther', 'KOG', 'KEGG/ec', 'KO', 'GO', 'Best-hit-arabi-name', 
                    'arabi-symbol', 'arabi-defline'
                ]
                reader = csv.DictReader(f, fieldnames=column_names, delimiter='\t')
                return list(reader)
    except Exception as e:
        print(f"Error reading annotation file: {e}")
        sys.exit(1)

def extract_sos_pathway_genes(annotation_data, blast_results=None, verbose=False):
    """Extract SOS pathway and related candidates from annotation data"""
    candidates = []
    
    if verbose:
        print(f"Processing {len(annotation_data)} annotation entries...")
    
    for i, row in enumerate(annotation_data):
        if is_sos_pathway_gene(row, blast_results):
            # Identify which component it is
            component, confidence, reasoning = identify_component(row, blast_results)
            
            # Add candidate to list
            row['component'] = component
            row['confidence'] = confidence
            row['reasoning'] = reasoning
            candidates.append(row)
            
            if verbose and (i+1) % 1000 == 0:
                print(f"Processed {i+1} entries, found {len(candidates)} candidates so far...")
    
    if verbose:
        print(f"Complete! Found {len(candidates)} SOS pathway and related candidates.")
    
    return candidates

def write_output(candidates, output_file):
    """Write SOS pathway candidates to output file"""
    # Define output columns
    output_columns = [
        'locusName', 'transcriptName', 'peptideName', 'component', 'confidence',
        'Pfam', 'Best-hit-arabi-name', 'arabi-symbol', 'arabi-defline', 'reasoning'
    ]
    
    try:
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=output_columns, delimiter='\t', 
                                   extrasaction='ignore')
            writer.writeheader()
            writer.writerows(candidates)
        print(f"SOS pathway and related candidates written to {output_file}")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

def main():
    """Main function to extract SOS pathway and related genes from annotation file"""
    print("\n=== extract_SOS_pathway.py - Extract SOS pathway and related salt stress genes ===")
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Read BLAST results if provided
    blast_results = read_blast_results(args.blast)
    if blast_results:
        print(f"Loaded {len(blast_results)} BLAST hits to reference components.")
    
    # Read annotation file
    print(f"Reading annotation file: {args.input}")
    annotation_data = read_annotation_file(args.input)
    print(f"Read {len(annotation_data)} entries from annotation file.")
    
    # Extract SOS pathway genes
    print("\nExtracting SOS pathway and related candidates...")
    candidates = extract_sos_pathway_genes(annotation_data, blast_results, args.verbose)
    print(f"Found {len(candidates)} gene candidates.")
    
    # Write output
    write_output(candidates, args.output)
    
    # Print summary
    print("\n=== Summary of candidates ===")
    confidence_counts = {}
    component_counts = {}
    
    for candidate in candidates:
        confidence = candidate['confidence']
        component = candidate['component']
        
        confidence_counts[confidence] = confidence_counts.get(confidence, 0) + 1
        component_counts[component] = component_counts.get(component, 0) + 1
    
    print("\nConfidence levels:")
    for confidence, count in sorted(confidence_counts.items(), 
                                   key=lambda x: {"High": 0, "Medium": 1, "Low": 2}.get(x[0], 3)):
        print(f"  {confidence}: {count}")
    
    print("\nComponent types:")
    for component, count in sorted(component_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  {component}: {count}")
    
    print("\nReference components used:")
    for accession, info in SOS_PATHWAY_REFS.items():
        alt_names_str = ", ".join(info['alt_names']) if info['alt_names'] else ""
        print(f"  {info['name']} ({info['gene_id']}): {accession} {alt_names_str}")
    
    print("\nDone! To use these results in downstream analyses, you can:")
    print(f"  1. Filter by confidence level (recommended: High or Medium)")
    print(f"  2. Extract protein/gene sequences using the peptideName/locusName")
    print(f"  3. Perform phylogenetic analysis to confirm relationships with known components")
    print(f"  4. For 'Unknown' or ambiguous components, perform more detailed sequence analysis")

if __name__ == "__main__":
    main()

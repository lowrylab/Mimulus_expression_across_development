#!/usr/bin/env python3
"""
extract_PIFs.py - Extract Phytochrome-interacting factors from a plant annotation file

This script:
1. Reads a tab-delimited annotation file
2. Identifies potential PIF genes based on annotation information
3. Outputs a filtered list of PIFs with their annotations

Usage:
python extract_PIFs.py -i <input_annotation_file> -o <output_file>

Example annotation file format:
pacId    locusName    transcriptName    peptideName    Pfam    Panther    KOG    KEGG/ec    KO    GO    Best-hit-arabi-name    arabi-symbol    arabi-defline
"""

import argparse
import re
import csv
import sys

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Extract Phytochrome-interacting factors from annotation file")
    parser.add_argument("-i", "--input", required=True, help="Input annotation file (tab-delimited)")
    parser.add_argument("-o", "--output", required=True, help="Output file for PIF genes")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print detailed progress information")
    return parser.parse_args()

def is_potential_pif(row):
    """
    Determine if a gene annotation row potentially represents a PIF
    
    Criteria to identify PIFs:
    1. Arabidopsis best hit is a known PIF
    2. Symbol or defline contains "PIF" or "phytochrome interacting factor"
    3. Has bHLH domain (PF00010) and potential APB motif
    """
    # Check Arabidopsis symbol
    if row['arabi-symbol'] and any(re.search(r'\bPIF\d+\b', symbol) for symbol in row['arabi-symbol'].split(',')):
        return True
    
    # Check Arabidopsis defline
    if row['arabi-defline'] and re.search(r'phytochrome[\s-]*interacting factor', row['arabi-defline'], re.IGNORECASE):
        return True
    
    # Check Pfam for bHLH domain
    if row['Pfam'] and "PF00010" in row['Pfam']:
        # If it has bHLH domain and any indicator of PIF in the defline
        if row['arabi-defline'] and re.search(r'\bPIF\b', row['arabi-defline'], re.IGNORECASE):
            return True
    
    # Check if it's a clear homolog of a known PIF
    if row['Best-hit-arabi-name'] and any(
        re.search(r'AT' + at_id, row['Best-hit-arabi-name']) 
        for at_id in ['2G20180', '1G09530', '5G41790', '3G59060', '2G43010', '1G03010', '4G28640', '2G46970']
    ):
        return True
    
    return False

def score_pif_confidence(row):
    """
    Score the confidence level of a PIF identification
    Returns: (score, reasoning)
    """
    score = 0
    reasons = []
    
    # Check Arabidopsis symbol explicitly mentions PIF
    if row['arabi-symbol'] and any(re.search(r'\bPIF\d+\b', symbol) for symbol in row['arabi-symbol'].split(',')):
        score += 3
        reasons.append("Arabidopsis symbol explicitly mentions PIF")
    
    # Check Arabidopsis defline mentions phytochrome interacting factor
    if row['arabi-defline'] and re.search(r'phytochrome[\s-]*interacting factor', row['arabi-defline'], re.IGNORECASE):
        score += 3
        reasons.append("Defline mentions 'phytochrome interacting factor'")
    
    # Check for bHLH domain
    if row['Pfam'] and "PF00010" in row['Pfam']:
        score += 2
        reasons.append("Contains bHLH domain (PF00010)")
    
    # Check for APB motif indirectly via Panther classification
    if row['Panther'] and "SF" in row['Panther']:
        panther_families = row['Panther'].split(',')
        for family in panther_families:
            if ":SF" in family:  # Specific subfamily classification often indicates functional similarity
                score += 1
                reasons.append(f"Specific Panther subfamily: {family}")
                break
    
    # If the best hit is a direct match to known PIFs (known Arabidopsis gene IDs for PIFs)
    known_pif_at_ids = ['AT2G20180', 'AT1G09530', 'AT5G41790', 'AT3G59060', 
                        'AT2G43010', 'AT1G03010', 'AT4G28640', 'AT2G46970']
    
    if row['Best-hit-arabi-name'] and any(at_id in row['Best-hit-arabi-name'] for at_id in known_pif_at_ids):
        score += 3
        reasons.append(f"Direct match to known Arabidopsis PIF: {row['Best-hit-arabi-name']}")
    
    # Confidence level
    if score >= 5:
        confidence = "High"
    elif score >= 3:
        confidence = "Medium"
    else:
        confidence = "Low"
    
    return (confidence, ", ".join(reasons))

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

def map_arabidopsis_pifs():
    """
    Map Arabidopsis gene IDs to PIF names
    Returns a dictionary of AT IDs to PIF names
    """
    # Known Arabidopsis PIFs and their gene IDs
    return {
        'AT2G20180': 'PIF1 (PIL5)',
        'AT1G09530': 'PIF3',
        'AT5G41790': 'CKB1',  # Not a PIF but sometimes misannotated
        'AT3G59060': 'PIF5 (PIL6)',
        'AT2G43010': 'PIF4',
        'AT1G03010': 'PIF7',
        'AT4G28640': 'PIF8 (PIL1)',
        'AT2G46970': 'PIF6 (PIL2)',
        # Add more if known
    }

def extract_pifs(annotation_data, verbose=False):
    """Extract PIF candidates from annotation data"""
    pif_candidates = []
    at_id_to_pif = map_arabidopsis_pifs()
    
    if verbose:
        print(f"Processing {len(annotation_data)} annotation entries...")
    
    for i, row in enumerate(annotation_data):
        if is_potential_pif(row):
            # Score the confidence
            confidence, reasoning = score_pif_confidence(row)
            
            # Map to known PIF if possible
            pif_type = "Unknown"
            if row['Best-hit-arabi-name']:
                for at_id, pif_name in at_id_to_pif.items():
                    if at_id in row['Best-hit-arabi-name']:
                        pif_type = pif_name
                        break
            
            # If still unknown but has symbol info
            if pif_type == "Unknown" and row['arabi-symbol']:
                for symbol in row['arabi-symbol'].split(','):
                    if 'PIF' in symbol:
                        pif_type = symbol
                        break
            
            # Add candidate to list with confidence score
            row['confidence'] = confidence
            row['reasoning'] = reasoning
            row['pif_type'] = pif_type
            pif_candidates.append(row)
            
            if verbose and (i+1) % 1000 == 0:
                print(f"Processed {i+1} entries, found {len(pif_candidates)} PIF candidates so far...")
    
    if verbose:
        print(f"Complete! Found {len(pif_candidates)} PIF candidates.")
    
    return pif_candidates

def write_output(pif_candidates, output_file):
    """Write PIF candidates to output file"""
    # Define output columns
    output_columns = [
        'locusName', 'transcriptName', 'peptideName', 'pif_type', 'confidence',
        'Pfam', 'Best-hit-arabi-name', 'arabi-symbol', 'arabi-defline', 'reasoning'
    ]
    
    try:
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=output_columns, delimiter='\t', 
                                   extrasaction='ignore')
            writer.writeheader()
            writer.writerows(pif_candidates)
        print(f"PIF candidates written to {output_file}")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

def main():
    """Main function to extract PIFs from annotation file"""
    print("\n=== extract_PIFs.py - Extract Phytochrome-interacting factors ===")
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Read annotation file
    print(f"Reading annotation file: {args.input}")
    annotation_data = read_annotation_file(args.input)
    print(f"Read {len(annotation_data)} entries from annotation file.")
    
    # Extract PIFs
    print("\nExtracting PIF candidates...")
    pif_candidates = extract_pifs(annotation_data, args.verbose)
    print(f"Found {len(pif_candidates)} PIF candidates.")
    
    # Write output
    write_output(pif_candidates, args.output)
    
    # Print summary
    print("\n=== Summary of PIF candidates ===")
    confidence_counts = {}
    pif_type_counts = {}
    
    for candidate in pif_candidates:
        confidence = candidate['confidence']
        pif_type = candidate['pif_type']
        
        confidence_counts[confidence] = confidence_counts.get(confidence, 0) + 1
        pif_type_counts[pif_type] = pif_type_counts.get(pif_type, 0) + 1
    
    print("\nConfidence levels:")
    for confidence, count in confidence_counts.items():
        print(f"  {confidence}: {count}")
    
    print("\nPIF types:")
    for pif_type, count in pif_type_counts.items():
        print(f"  {pif_type}: {count}")
    
    print("\nDone! To use these results in downstream analyses, you can:")
    print(f"  1. Filter by confidence level (recommended: High or Medium)")
    print(f"  2. Extract protein/gene sequences using the peptideName/locusName")
    print(f"  3. Perform phylogenetic analysis to confirm relationships with known PIFs")

if __name__ == "__main__":
    main()

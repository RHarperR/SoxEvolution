#!/usr/bin/env python3
"""
Extract specific gene node information from multiple Ranger-DTL output files.
Usage: python extract_ranger_node.py <folder> <gene_node> <output_file>
       python extract_ranger_node.py ranger_results m83 m83_results.tsv
       python extract_ranger_node.py . m83 m83_results.tsv
"""

import sys
import re
import csv
from pathlib import Path


def parse_reconciliation_line(line):
    """
    Parse a reconciliation line and extract event information.
    
    Returns:
        dict with keys: gene_node, event, spe_node_mapping, spe_node_recipient
        or None if not a reconciliation line
    """
    line = line.strip()
    
    # Skip empty lines and headers
    if not line or line.startswith('-') or line.startswith('Species Tree:') or line.startswith('Gene Tree:') or line.startswith('Reconciliation:'):
        return None
    
    # Pattern for leaf nodes: "GCA0017883951: Leaf Node"
    leaf_pattern = r'^(\S+):\s*Leaf Node'
    leaf_match = re.match(leaf_pattern, line)
    if leaf_match:
        return {
            'gene_node': leaf_match.group(1),
            'event': 'Leaf',
            'spe_node_mapping': '',
            'spe_node_recipient': ''
        }
    
    # Pattern for internal nodes with events
    # m83 = LCA[...]: Event, Mapping --> node, [Recipient --> node]
    # Use [^,\s] to exclude comma and whitespace from mapping node capture
    event_pattern = r'^(\S+)\s*=\s*LCA\[.*?\]:\s*(\w+)\s*,\s*Mapping\s*-->\s*([^,\s]+)(?:\s*,\s*Recipient\s*-->\s*(\S+))?'
    event_match = re.match(event_pattern, line)
    
    if event_match:
        gene_node = event_match.group(1)
        event = event_match.group(2)
        mapping = event_match.group(3)
        recipient = event_match.group(4) if event_match.group(4) else ''
        
        return {
            'gene_node': gene_node,
            'event': event,
            'spe_node_mapping': mapping,
            'spe_node_recipient': recipient
        }
    
    return None


def extract_node_from_file(input_file, target_node):
    """
    Extract specific gene node information from a Ranger output file.
    
    Args:
        input_file: Path to the Ranger output file
        target_node: Gene node to search for (e.g., 'm83')
        
    Returns:
        Dictionary with event information, or None if not found
    """
    in_reconciliation = False
    
    with open(input_file, 'r') as f:
        for line in f:
            # Check if we've reached the Reconciliation section
            if line.strip().startswith('Reconciliation:'):
                in_reconciliation = True
                continue
            
            # Stop if we hit the next section or end
            if in_reconciliation and line.strip().startswith('---'):
                break
            
            # Parse lines in the Reconciliation section
            if in_reconciliation:
                # Check if line starts with target node
                if line.strip().startswith(f"{target_node} =") or line.strip().startswith(f"{target_node}:"):
                    event = parse_reconciliation_line(line)
                    if event and event['gene_node'] == target_node:
                        return event
    
    return None


def parse_filename(filename):
    """
    Parse ranger filename to extract parameters.
    Expected format: ranger.D.T.L or ranger.D.T.L.XXX
    
    Returns:
        dict with D, T, L, replicate values
    """
    stem = Path(filename).stem
    parts = stem.split('.')
    
    result = {
        'D': '',
        'T': '',
        'L': '',
        'replicate': ''
    }
    
    if len(parts) >= 4 and parts[0] == 'ranger':
        result['D'] = parts[1]
        result['T'] = parts[2]
        result['L'] = parts[3]
        if len(parts) >= 5:
            result['replicate'] = parts[4]
    
    return result


def process_folder(folder_path, target_node):
    """
    Process all ranger files in a folder and extract specific node information.
    
    Args:
        folder_path: Path to folder containing ranger files
        target_node: Gene node to search for
        
    Returns:
        List of dictionaries with file info and event data
    """
    folder = Path(folder_path)
    results = []
    
    # Find all files starting with 'ranger'
    ranger_files = sorted(folder.glob('ranger*'))
    
    if not ranger_files:
        print(f"Warning: No ranger files found in {folder_path}")
        return results
    
    print(f"Found {len(ranger_files)} ranger files")
    print(f"Searching for gene node: {target_node}")
    print()
    
    processed = 0
    found = 0
    
    for file_path in ranger_files:
        processed += 1
        
        # Parse filename to get parameters
        params = parse_filename(file_path.name)
        
        # Extract node information
        event = extract_node_from_file(file_path, target_node)
        
        if event:
            found += 1
            result = {
                'filename': file_path.name,
                'D': params['D'],
                'T': params['T'],
                'L': params['L'],
                'replicate': params['replicate'],
                'gene_node': event['gene_node'],
                'event': event['event'],
                'spe_node_mapping': event['spe_node_mapping'],
                'spe_node_recipient': event['spe_node_recipient']
            }
            results.append(result)
            
            # Print progress every 100 files
            if processed % 100 == 0:
                print(f"Processed {processed}/{len(ranger_files)} files, found {found} matches")
    
    print(f"\nCompleted: Processed {processed} files, found {found} matches")
    
    return results


def main():
    if len(sys.argv) < 4:
        print("Usage: python extract_ranger_node.py <folder> <gene_node> <output_file>")
        print("\nExamples:")
        print("  python extract_ranger_node.py ranger_results m83 m83_results.tsv")
        print("  python extract_ranger_node.py . m2 m2_results.tsv")
        print("  python extract_ranger_node.py /path/to/results m100 output.tsv")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    target_node = sys.argv[2]
    output_file = sys.argv[3]
    
    # Check if folder exists
    if not Path(folder_path).exists():
        print(f"Error: Folder '{folder_path}' not found!")
        sys.exit(1)
    
    if not Path(folder_path).is_dir():
        print(f"Error: '{folder_path}' is not a directory!")
        sys.exit(1)
    
    # Process all files
    results = process_folder(folder_path, target_node)
    
    if not results:
        print(f"\nWarning: Gene node '{target_node}' not found in any files!")
        return
    
    # Write to output file
    fieldnames = ['filename', 'D', 'T', 'L', 'replicate', 'gene_node', 'event', 'spe_node_mapping', 'spe_node_recipient']
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nResults written to: {output_file}")
    print(f"Total records: {len(results)}")
    
    # Print summary statistics
    event_counts = {}
    for result in results:
        event_type = result['event']
        event_counts[event_type] = event_counts.get(event_type, 0) + 1
    
    print("\nEvent Type Summary:")
    for event_type, count in sorted(event_counts.items()):
        print(f"  {event_type}: {count}")
    
    # Print parameter summary if multiple files
    if len(results) > 1:
        d_values = set(r['D'] for r in results if r['D'])
        l_values = set(r['L'] for r in results if r['L'])
        replicates = set(r['replicate'] for r in results if r['replicate'])
        
        print("\nParameter Summary:")
        if d_values:
            print(f"  D values: {sorted(d_values)}")
        if l_values:
            print(f"  L values: {sorted(l_values)}")
        if replicates:
            print(f"  Replicates: {len(replicates)} unique values")


if __name__ == "__main__":
    main()
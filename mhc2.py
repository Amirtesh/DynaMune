#!/usr/bin/env python3
# version1.0
"""
MHC-II Binding Prediction using IEDB Next-Generation API
Reads sequences from test.fasta and alleles from alleles.txt
Predicts binding for peptide lengths 8-11
"""

import requests
import json
import time
import csv
from pathlib import Path

# API Configuration
API_BASE = "https://api-nextgen-tools.iedb.org/api/v1"
PIPELINE_ENDPOINT = f"{API_BASE}/pipeline"
RESULTS_ENDPOINT = f"{API_BASE}/results"

def read_fasta(fasta_file):
    """Read sequences from FASTA file"""
    sequences = []
    current_seq = ""
    current_header = ""
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append((current_header, current_seq))
                current_header = line
                current_seq = ""
            else:
                current_seq += line
        
        if current_seq:
            sequences.append((current_header, current_seq))
    
    # Format as FASTA text
    fasta_text = "\n".join([f"{header}\n{seq}" for header, seq in sequences])
    return fasta_text

def read_alleles(alleles_file):
    """Read alleles from text file (one per line or comma-separated)"""
    with open(alleles_file, 'r') as f:
        content = f.read().strip()
    
    # If already comma-separated, return as is
    if ',' in content and '\n' not in content:
        return content
    
    # If one per line, convert to comma-separated
    alleles = [line.strip() for line in content.split('\n') if line.strip()]
    return ','.join(alleles)

def submit_mhc2_prediction(fasta_text, alleles, peptide_lengths=[9, 15]):
    """
    Submit MHC-II prediction job to IEDB API
    
    Args:
        fasta_text: FASTA formatted sequences
        alleles: Comma-separated allele list (e.g., "HLA-DRB1*01:01,HLA-DRB1*03:01")
        peptide_lengths: List [min, max] for peptide length range
    
    Returns:
        result_id: Job result ID for polling
        pipeline_id: Pipeline ID
    """
    
    payload = {
        "pipeline_title": "MHC-II Binding Prediction",
        "run_stage_range": [1, 1],
        "stages": [
            {
                "stage_number": 1,
                "tool_group": "mhcii",  # Changed from "mhci" to "mhcii"
                "input_sequence_text": fasta_text,
                "input_parameters": {
                    "alleles": alleles,
                    "peptide_length_range": peptide_lengths,
                    "predictors": [
                        {
                            "type": "binding",
                            "method": "netmhciipan_el"  # IEDB recommended method for MHC-II
                        }
                    ]
                }
            }
        ]
    }
    
    print("Submitting MHC-II prediction job...")
    print(f"Alleles: {alleles}")
    print(f"Peptide lengths: {peptide_lengths}")
    
    response = requests.post(
        PIPELINE_ENDPOINT,
        headers={
            "accept": "application/json",
            "Content-Type": "application/json"
        },
        json=payload
    )
    
    if response.status_code != 200:
        raise Exception(f"API Error: {response.status_code}\n{response.text}")
    
    result = response.json()
    print(f"Job submitted successfully!")
    print(f"Result ID: {result['result_id']}")
    print(f"Pipeline ID: {result['pipeline_id']}")
    
    return result['result_id'], result['pipeline_id']

def poll_results(result_id, max_wait=600, poll_interval=10):
    """
    Poll for job results
    
    Args:
        result_id: Job result ID
        max_wait: Maximum wait time in seconds
        poll_interval: Time between polls in seconds
    
    Returns:
        results: JSON results from API
    """
    
    results_url = f"{RESULTS_ENDPOINT}/{result_id}"
    elapsed = 0
    
    print(f"\nPolling for results (max wait: {max_wait}s)...")
    
    while elapsed < max_wait:
        response = requests.get(results_url)
        
        if response.status_code != 200:
            print(f"Warning: Got status code {response.status_code}, retrying...")
            time.sleep(poll_interval)
            elapsed += poll_interval
            continue
            
        results = response.json()
        
        # Check if results are ready
        status = results.get('status', 'pending')
        if status in ['complete', 'done']:
            print("Results ready!")
            return results
        elif status == 'error':
            raise Exception(f"Job failed: {results.get('message', 'Unknown error')}")
        else:
            print(f"Job status: {status} (waited {elapsed}s)")
        
        time.sleep(poll_interval)
        elapsed += poll_interval
    
    raise TimeoutError(f"Results not ready after {max_wait}s")

def parse_and_save_results(results, output_file="mhc2_prediction.csv"):
    """
    Parse results and save to CSV
    
    Args:
        results: JSON results from API
        output_file: Output CSV filename
    """
    
    print(f"\nParsing results and saving to {output_file}...")
    
    # Navigate to the data tables
    # Structure: results -> data -> results (list of tables)
    if 'data' not in results:
        raise Exception("No data found in results")
    
    data = results['data']
    
    if 'results' not in data:
        raise Exception("No results found in data")
    
    results_list = data['results']
    
    # results_list is a list of data tables
    if not isinstance(results_list, list) or len(results_list) == 0:
        raise Exception("Results is not a list or is empty")
    
    # Find the predictions table (look for one with allele/peptide columns)
    predictions_table = None
    for table in results_list:
        if 'table_columns' in table:
            col_names = [col['name'] for col in table['table_columns']]
            if 'allele' in col_names or 'peptide' in col_names:
                predictions_table = table
                break
    
    if predictions_table is None:
        print(f"Available tables: {[t.get('type', 'unknown') for t in results_list]}")
        raise Exception("No predictions table found")
    
    # Extract column names
    columns = [col['name'] for col in predictions_table['table_columns']]
    
    # Get the data - it might be in 'table_data' as a dict with column names as keys
    # containing lists of values
    if 'table_data' not in predictions_table:
        raise Exception("No table_data found in predictions table")
    
    table_data = predictions_table['table_data']
    
    # If table_data is a dict with column names as keys
    if isinstance(table_data, dict):
        # Convert dict format to rows
        num_rows = len(table_data[columns[0]])
        rows = []
        for i in range(num_rows):
            row = [table_data[col][i] if col in table_data else None for col in columns]
            rows.append(row)
    else:
        # table_data is already a list of rows
        rows = table_data
    
    print(f"Found {len(rows)} predictions with columns: {columns}")
    
    # Write to CSV
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(columns)  # Write header
        writer.writerows(rows)    # Write data rows
    
    print(f"Successfully saved {len(rows)} predictions to {output_file}")
    
    # Print summary statistics
    print("\n=== Prediction Summary ===")
    print(f"Total predictions: {len(rows)}")
    
    # Try to find allele column index for summary
    try:
        allele_idx = columns.index('allele')
        allele_counts = {}
        for row in rows:
            allele = row[allele_idx]
            allele_counts[allele] = allele_counts.get(allele, 0) + 1
        
        print("\nPredictions by allele:")
        for allele, count in sorted(allele_counts.items()):
            print(f"  {allele}: {count}")
    except (ValueError, IndexError):
        pass  # Allele column not found, skip summary

def main():
    """Main execution function"""
    
    import sys
    
    # Check command line arguments
    if len(sys.argv) != 3:
        print("Usage: python3 mhc2_prediction.py <fasta_file> <alleles_file>")
        print("Example: python3 mhc2_prediction.py test.fasta alleles.txt")
        sys.exit(1)
    
    # Get file paths from command line
    fasta_file = sys.argv[1]
    alleles_file = sys.argv[2]
    output_file = "mhc2_prediction.csv"
    
    # Check if input files exist
    if not Path(fasta_file).exists():
        print(f"Error: {fasta_file} not found!")
        sys.exit(1)
    
    if not Path(alleles_file).exists():
        print(f"Error: {alleles_file} not found!")
        sys.exit(1)
    
    try:
        # Read input files
        print("Reading input files...")
        fasta_text = read_fasta(fasta_file)
        alleles = read_alleles(alleles_file)
        
        print(f"Loaded sequences from {fasta_file}")
        print(f"Loaded alleles: {alleles}")
        
        # Submit prediction job
        result_id, pipeline_id = submit_mhc2_prediction(
            fasta_text=fasta_text,
            alleles=alleles,
            peptide_lengths=[9, 15]  # Recommended range for MHC-II
        )
        
        # Poll for results
        results = poll_results(result_id, max_wait=600, poll_interval=10)
        
        # Parse and save results
        parse_and_save_results(results, output_file)
        
        print(f"\n✓ MHC-II prediction completed successfully!")
        print(f"✓ Results saved to: {output_file}")
        
    except Exception as e:
        print(f"\n✗ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

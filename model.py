#!/usr/bin/env python3
# version1.0
"""
ESMFold API Protein Structure Prediction
Predict protein 3D structures using Meta's ESMFold public API.
NO authentication required - completely FREE!

Usage: python3 model.py sequence.txt

The ESM Atlas API provides fast structure prediction (up to 60x faster than AlphaFold2)
without requiring multiple sequence alignments or any login credentials.
"""

import sys
import requests
from pathlib import Path


def read_sequence_file(filepath):
    """
    Read protein sequence from a text file.
    
    Args:
        filepath: Path to the sequence file
        
    Returns:
        Clean protein sequence string
    """
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Remove FASTA header if present
        lines = content.strip().split('\n')
        sequence = ''
        for line in lines:
            if not line.startswith('>') and not line.startswith(';'):
                # Remove whitespace and common formatting
                sequence += line.strip().replace(' ', '').replace('\n', '').replace('\r', '')
        
        # Validate sequence (basic check for amino acid letters)
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        sequence = sequence.upper()
        
        if not sequence:
            raise ValueError("No sequence found in file!")
        
        # Check for invalid characters
        invalid_chars = set(sequence) - valid_aa
        if invalid_chars:
            print(f"Warning: Found non-standard amino acid characters: {invalid_chars}")
            print("These will be removed from the sequence.")
            # Remove invalid characters
            sequence = ''.join(c for c in sequence if c in valid_aa)
        
        if not sequence:
            raise ValueError("No valid amino acids found after cleaning!")
        
        return sequence
    
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found!")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)


def predict_structure_esmfold(sequence, output_file="predicted_structure.pdb"):
    """
    Predict protein structure using ESMFold public API.
    
    Args:
        sequence: Protein amino acid sequence (single letter code)
        output_file: Path to save the PDB file
        
    Returns:
        True if successful, False otherwise
    """
    # ESM Atlas API endpoint
    api_url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    
    # Check sequence length (API has a maximum of 400 amino acids)
    if len(sequence) > 400:
        print(f"\nWarning: Sequence length ({len(sequence)}) exceeds ESMFold API limit (400 AA).")
        print("Truncating to first 400 amino acids...")
        sequence = sequence[:400]
    
    print("\n" + "="*60)
    print("ESMFold Structure Prediction")
    print("="*60)
    print(f"Sequence length: {len(sequence)} amino acids")
    print(f"First 50 AA: {sequence[:50]}...")
    print("\nSubmitting to ESM Atlas API...")
    print("This typically takes 10-30 seconds depending on sequence length.")
    print("-"*60)
    
    try:
        # Make POST request to ESMFold API
        response = requests.post(
            api_url,
            data=sequence,
            headers={'Content-Type': 'text/plain'},
            timeout=120  # 2 minute timeout
        )
        
        # Check if request was successful
        response.raise_for_status()
        
        # Get PDB content
        pdb_content = response.text
        
        # Validate PDB content
        if not pdb_content or not pdb_content.startswith('HEADER'):
            raise ValueError("Invalid PDB response from API")
        
        # Save to file
        with open(output_file, 'w') as f:
            f.write(pdb_content)
        
        print("-"*60)
        print("✓ Prediction completed successfully!")
        print("\n" + "="*60)
        print("SUCCESS!")
        print("="*60)
        print(f"\nPredicted structure saved to: {Path(output_file).absolute()}")
        print("\nYou can visualize the structure using:")
        print("  • PyMOL:     pymol predicted_structure.pdb")
        print("  • ChimeraX:  chimerax predicted_structure.pdb")
        print("  • Online:    https://www.rcsb.org/3d-view")
        print("  • Mol*:      https://molstar.org/viewer/")
        
        # Extract and display confidence information if available
        try:
            confidence_lines = [line for line in pdb_content.split('\n') if 'pLDDT' in line or 'B-factor' in line]
            if confidence_lines:
                print("\nNote: Confidence scores are stored in the B-factor column")
                print("      Higher values (closer to 100) indicate higher confidence")
        except:
            pass
        
        print("="*60 + "\n")
        return True
        
    except requests.exceptions.Timeout:
        print("\n✗ Error: Request timed out. The server may be busy.")
        print("   Please try again in a few moments.")
        return False
    
    except requests.exceptions.HTTPError as e:
        print(f"\n✗ HTTP Error: {e}")
        if e.response.status_code == 400:
            print("   The sequence may contain invalid characters.")
        elif e.response.status_code == 413:
            print("   The sequence is too long for the API.")
        elif e.response.status_code == 503:
            print("   The ESMFold server is currently unavailable.")
            print("   Please try again later.")
        return False
    
    except requests.exceptions.ConnectionError:
        print("\n✗ Connection Error: Could not connect to ESM Atlas API.")
        print("   Please check your internet connection.")
        return False
    
    except ValueError as e:
        print(f"\n✗ Error: {e}")
        return False
    
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        return False


def main():
    """Main execution function."""
    
    # ASCII art banner
    print("""
╔══════════════════════════════════════════════════════════════╗
║                                                              ║
║   ESMFold Protein Structure Prediction (FREE API)           ║
║   No authentication required | No installation needed       ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
""")
    
    # Check command line arguments
    if len(sys.argv) != 2:
        print("Usage: python3 model.py sequence.txt")
        print("\nThe sequence.txt file should contain:")
        print("  • Plain protein sequence (single letter amino acid code)")
        print("  • FASTA format is also supported")
        print("  • Maximum 400 amino acids")
        print("\nExample sequence.txt:")
        print("  MKTIIALSYIFCLVFADYKDDDDK...")
        print("\nOr in FASTA format:")
        print("  >MyProtein")
        print("  MKTIIALSYIFCLVFADYKDDDDK...")
        sys.exit(1)
    
    sequence_file = sys.argv[1]
    
    # Read sequence from file
    print(f"Reading sequence from: {sequence_file}")
    sequence = read_sequence_file(sequence_file)
    
    # Predict structure
    success = predict_structure_esmfold(
        sequence=sequence,
        output_file="predicted_structure.pdb"
    )
    
    if not success:
        print("\nPrediction failed. Please check the error messages above.")
        sys.exit(1)


if __name__ == "__main__":
    main()

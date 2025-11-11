#!/usr/bin/env python3
"""
ProDy Deformation Analysis Script

This script compares the intrinsic dynamics (NMA) of a reference structure 
(e.g., unbound vaccine) to an observed conformational change (e.g., 
when bound in a complex).

Usage:
  python3 prody_deformation_analysis.py \
    --reference protein_A.pdb \
    --target protein_A_complex.pdb \
    --ref_chain A \
    --target_chain C

This compares the 'unbound' protein_A (chain A) with the 'bound'
conformation (chain C) from the complex.
"""

import sys
import os
import argparse
import numpy as np  # <-- THIS IS THE FIX
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform

try:
    from prody import *
except ImportError:
    print("Error: ProDy is not installed. Install with: pip install prody")
    sys.exit(1)


def setup_output_directory(pdb_name):
    """Create output directory for results"""
    output_dir = f"{pdb_name}_Deformation_Analysis"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir


def plot_deformation_overlap(overlaps_sq, cumulative_overlap, output_dir):
    """Plot the results of the deformation analysis"""
    print("Generating output plot...")
    
    n_modes = len(overlaps_sq)
    mode_indices = np.arange(1, n_modes + 1)
    
    fig, ax1 = plt.subplots(figsize=(12, 7))
    
    # Bar plot for individual overlap
    color = 'tab:blue'
    ax1.set_xlabel('ANM Mode Index', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Squared Overlap (Contribution)', color=color, fontsize=12, fontweight='bold')
    ax1.bar(mode_indices, overlaps_sq, color=color, alpha=0.7, label='Individual Mode Contribution')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_ylim([0, 1])
    ax1.set_xlim([0.5, n_modes + 0.5])

    # Line plot for cumulative overlap
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('Cumulative Overlap', color=color, fontsize=12, fontweight='bold')
    ax2.plot(mode_indices, cumulative_overlap, color=color, marker='o', linestyle='--', linewidth=2.5, label='Cumulative Overlap')
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylim([0, 1.05])

    # Add threshold lines
    ax2.axhline(y=0.5, color='gray', linestyle=':', linewidth=1.5, label='50% Cumulative')
    ax2.axhline(y=0.8, color='black', linestyle=':', linewidth=1.5, label='80% Cumulative')
    
    plt.title('NMA Mode Contribution to Conformational Change', fontsize=16, fontweight='bold')
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    
    # Add a combined legend
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='center right')
    
    plt.savefig(f'{output_dir}/deformation_overlap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/deformation_overlap.png")


def save_results_to_file(overlaps_sq, cumulative_overlap, output_dir, n_modes):
    """Save numerical overlap data to a text file"""
    print("Saving numerical results to file...")
    
    with open(f'{output_dir}/deformation_results.txt', 'w') as f:
        f.write("NMA Mode Deformation Analysis Results\n")
        f.write("="*60 + "\n\n")
        f.write("This file shows how much each NMA mode (intrinsic motion) \n")
        f.write("contributes to the observed conformational change (deformation).\n\n")
        
        f.write(f"{'Mode':<5} {'Squared_Overlap':<18} {'Cumulative_Overlap':<20}\n")
        f.write(f"{'-'*4:<5} {'-'*17:<18} {'-'*19:<20}\n")
        
        for i in range(n_modes):
            sq_overlap = overlaps_sq[i]
            cum_overlap = cumulative_overlap[i]
            f.write(f"{i+1:<5} {sq_overlap:<18.4f} {cum_overlap:<20.4f}\n")
            
    print(f"Saved: {output_dir}/deformation_results.txt")


def check_modes(value):
    """Custom argparse type to check for mode range"""
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid mode value: {value}. Must be an integer.")
    
    if ivalue < 1 or ivalue > 100:
        raise argparse.ArgumentTypeError(f"Number of modes ({value}) must be between 1 and 100.")
    return ivalue


def main():
    parser = argparse.ArgumentParser(
        description='ProDy Deformation Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python3 prody_deformation_analysis.py \
    --reference unbound_protein.pdb \
    --target bound_complex.pdb \
    --ref_chain A \
    --target_chain C
        """)
    
    parser.add_argument('--reference', required=True, help='Reference PDB file (e.g., unbound state)')
    parser.add_argument('--target', required=True, help='Target PDB file (e.g., bound state)')
    parser.add_argument('--ref_chain', required=True, help='Chain ID of the protein in the reference file')
    parser.add_argument('--target_chain', required=True, help='Chain ID of the protein in the target file')
    parser.add_argument('--nmodes', type=check_modes, default=20, help='Number of NMA modes to calculate (default: 20, max: 100)')
    parser.add_argument('--cutoff', type=float, default=15.0, help='ANM cutoff distance in Å (default: 15.0)')
    parser.add_argument('--selection', default='calpha', help='Atom selection (default: calpha)')
    
    args = parser.parse_args()

    # 1. Setup
    pdb_name = os.path.splitext(os.path.basename(args.reference))[0]
    output_dir = setup_output_directory(pdb_name)

    # 2. Load Structures
    print(f"Loading reference: {args.reference} (Chain: {args.ref_chain})")
    ref_struct = parsePDB(args.reference)
    ref_atoms = ref_struct.select(f'chain {args.ref_chain} and {args.selection}')
    if ref_atoms is None:
        print(f"Error: No atoms found for chain {args.ref_chain} and selection '{args.selection}' in {args.reference}")
        sys.exit(1)

    print(f"Loading target: {args.target} (Chain: {args.target_chain})")
    target_struct = parsePDB(args.target)
    target_atoms = target_struct.select(f'chain {args.target_chain} and {args.selection}')
    if target_atoms is None:
        print(f"Error: No atoms found for chain {args.target_chain} and selection '{args.selection}' in {args.target}")
        sys.exit(1)

    print(f"Found {ref_atoms.numAtoms()} atoms in reference, {target_atoms.numAtoms()} atoms in target.")

    # 3. Match and Align
    print("Matching and aligning structures...")
    try:
        matches = matchChains(ref_atoms, target_atoms)
        if not matches:
            raise ValueError("No matching chains found.")
        
        # Robustly get the matched AtomMap objects
        match_tuple = matches[0]
        ref_atoms_matched = match_tuple[0]
        target_atoms_matched = match_tuple[1]
        
        print(f"Successfully matched {ref_atoms_matched.numAtoms()} atoms.")

        # Superpose the structures
        # Use the function prody.superpose(), which works on AtomMap objects
        superpose(target_atoms_matched, ref_atoms_matched)

    except Exception as e:
        print(f"Error during atom matching: {e}")
        print("Please ensure the specified chains correspond to the same protein.")
        sys.exit(1)

    # 4. Perform ANM on Reference
    print("\nPerforming ANM on reference structure...")
    anm = ANM(ref_atoms_matched.getTitle())
    anm.buildHessian(ref_atoms_matched, cutoff=args.cutoff)
    anm.calcModes(n_modes=args.nmodes)
    
    model = anm[:args.nmodes]
    print(f"Calculated {len(model)} non-zero modes.")

    # 5. Calculate Deformation Vector
    print("\nCalculating deformation vector...")
    # The deformation vector is the 3D vector field showing the movement
    # from the reference to the target conformation.
    deform_vec = calcDeformVector(ref_atoms_matched, target_atoms_matched)

    # 6. Calculate Overlap
    print("Calculating overlap between NMA modes and deformation...")
    
    # This is the fix, following the ProDy tutorial's method
    # This calculates the squared overlap for all modes in the 'model'
    overlaps_sq = calcOverlap(model, deform_vec)

    # This is the key fix for the "IndexError: invalid index to scalar variable."
    # It ensures overlaps_sq is *always* an array, even if calcOverlap
    # returns a single scalar value (which can happen).
    overlaps_sq = np.atleast_1d(overlaps_sq)
    
    cumulative_overlap = np.cumsum(overlaps_sq)

    # 7. Report Results
    print("\n" + "="*60)
    print("DEFORMATION ANALYSIS RESULTS")
    print("="*60)
    
    # Find how many modes to reach 80% overlap
    try:
        modes_80 = np.where(cumulative_overlap > 0.8)[0][0] + 1
        print(f"\n{modes_80} modes are required to capture 80% of the deformation.")
    except IndexError:
        print("\nThe selected modes do not capture 80% of the deformation.")

    print(f"\nOverlap of deformation with ANM Mode 1: {overlaps_sq[0]*100:.2f}%")
    if len(overlaps_sq) > 1:
        print(f"Overlap of deformation with ANM Mode 2: {overlaps_sq[1]*100:.2f}%")
    if len(overlaps_sq) > 2:
        print(f"Overlap of deformation with ANM Mode 3: {overlaps_sq[2]*100:.2f}%")
    
    if len(cumulative_overlap) >= 3:
        print(f"\nCumulative overlap of first 3 modes: {cumulative_overlap[2]*100:.2f}%")
    else:
        print(f"\nCumulative overlap of all {len(cumulative_overlap)} modes: {cumulative_overlap[-1]*100:.2f}%")
    print("\n" + "="*60)

    # 8. Plot and Save
    plot_deformation_overlap(overlaps_sq, cumulative_overlap, output_dir)
    save_results_to_file(overlaps_sq, cumulative_overlap, output_dir, len(overlaps_sq))
    
    print("\nDeformation analysis complete.")
    print(f"All results saved to: {output_dir}/")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# version1.0
"""
ProDy Perturbation Response Scanning (PRS) Analysis Script

This script performs PRS analysis to identify residues that regulate
protein dynamics and allosteric communication pathways.

Usage:
  # Whole structure (default)
  python3 prody_prs_analysis.py --pdb structure.pdb

  # Specific chain only
  python3 prody_prs_analysis.py --pdb structure.pdb --chain A

  # Specific residue range
  python3 prody_prs_analysis.py --pdb structure.pdb --chain A --residue-range 20-150

  # Exclude regions (e.g., flexible tails)
  python3 prody_prs_analysis.py --pdb structure.pdb --chain A --exclude-residues 1-15,200-220

  # Custom selection (epitope, interface, etc.)
  python3 prody_prs_analysis.py --pdb structure.pdb --selection "chain A and resnum 25 to 120"

  # With specific effector-sensor analysis
  python3 prody_prs_analysis.py --pdb structure.pdb --chain A \
    --effector-residues 45,67,89 --sensor-residues 120,145,167
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

try:
    from prody import *  # type: ignore
except ImportError:
    print("Error: ProDy is not installed. Install with: pip install prody")
    sys.exit(1)


def setup_output_directory(pdb_name):
    """Create output directory for results"""
    output_dir = f"{pdb_name}_PRS_Analysis"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir


def analyze_chain_structure(atoms):
    """Analyze chain structure and create mapping"""
    chains = atoms.getChids()
    resnums = atoms.getResnums()
    resnames = atoms.getResnames()
    
    residue_indices = np.arange(1, len(chains) + 1)
    
    chain_info = []
    current_chain = chains[0]
    start_idx = 0
    
    for i in range(1, len(chains)):
        if chains[i] != current_chain:
            chain_info.append({
                'chain': current_chain,
                'start_idx': start_idx,
                'end_idx': i,
                'start_resnum': resnums[start_idx],
                'end_resnum': resnums[i-1],
                'length': i - start_idx
            })
            current_chain = chains[i]
            start_idx = i
    
    chain_info.append({
        'chain': current_chain,
        'start_idx': start_idx,
        'end_idx': len(chains),
        'start_resnum': resnums[start_idx],
        'end_resnum': resnums[-1],
        'length': len(chains) - start_idx
    })
    
    return residue_indices, chain_info


def print_chain_summary(chain_info):
    """Print summary of chain structure"""
    print(f"\n{'='*60}")
    print("CHAIN STRUCTURE ANALYSIS")
    print(f"{'='*60}")
    print(f"Number of chains: {len(chain_info)}")
    print(f"\nChain details:")
    for info in chain_info:
        print(f"  Chain {info['chain']}: {info['length']} residues "
              f"(PDB: {info['start_resnum']}-{info['end_resnum']}, "
              f"Index: {info['start_idx']+1}-{info['end_idx']})")
    print(f"{'='*60}")


def identify_interface_residues(atoms, chain_info, distance_cutoff=10.0):
    """Identify residues at the interface between chains"""
    if len(chain_info) < 2:
        return [], {}
    
    print(f"\nIdentifying interface residues (cutoff: {distance_cutoff} Å)...")
    
    coords = atoms.getCoords()
    chains = atoms.getChids()
    interface_residues = []
    interface_pairs = {}
    
    for i, info1 in enumerate(chain_info[:-1]):
        for info2 in chain_info[i+1:]:
            mask1 = (chains == info1['chain'])
            mask2 = (chains == info2['chain'])
            coords1 = coords[mask1]
            coords2 = coords[mask2]
            
            for idx1 in range(len(coords1)):
                global_idx1 = info1['start_idx'] + idx1
                for idx2 in range(len(coords2)):
                    global_idx2 = info2['start_idx'] + idx2
                    dist = np.linalg.norm(coords1[idx1] - coords2[idx2])
                    if dist < distance_cutoff:
                        if global_idx1 not in interface_residues:
                            interface_residues.append(global_idx1)
                            interface_pairs[global_idx1] = []
                        if global_idx2 not in interface_residues:
                            interface_residues.append(global_idx2)
                            interface_pairs[global_idx2] = []
                        interface_pairs[global_idx1].append((global_idx2, dist))
                        interface_pairs[global_idx2].append((global_idx1, dist))
    
    interface_residues.sort()
    print(f"Found {len(interface_residues)} interface residues")
    return interface_residues, interface_pairs


def add_chain_boundaries(ax, chain_info, y_fraction=0.95):
    """Add chain boundary markers and labels to plot"""
    ylim = ax.get_ylim()
    y_pos = ylim[0] + (ylim[1] - ylim[0]) * y_fraction
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(chain_info)))
    
    for i, info in enumerate(chain_info):
        if i > 0:
            ax.axvline(x=info['start_idx']+1, color='gray', linestyle='--', 
                      alpha=0.4, linewidth=1.5, zorder=1)
        
        ax.add_patch(plt.Rectangle((info['start_idx']+0.5, ylim[0]), 
                               info['length'], ylim[1]-ylim[0],
                               facecolor=colors[i], alpha=0.1, zorder=0))
        
        mid_point = info['start_idx'] + info['length'] / 2
        ax.text(mid_point, y_pos, f"Chain {info['chain']}\n({info['length']} res)", 
               ha='center', va='top', fontsize=9, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.4', facecolor=colors[i], 
                        alpha=0.7, edgecolor='black', linewidth=1))


def parse_residue_range(range_str):
    """Parse residue range string like '20-150' into (start, end)"""
    if not range_str:
        return None
    try:
        start, end = map(int, range_str.split('-'))
        return (start, end)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid residue range: {range_str}. Use format: START-END")


def parse_exclude_residues(exclude_str):
    """Parse exclude residues string like '1-15,200-220' into list of ranges"""
    if not exclude_str:
        return []
    ranges = []
    for part in exclude_str.split(','):
        try:
            start, end = map(int, part.strip().split('-'))
            ranges.append((start, end))
        except ValueError:
            raise argparse.ArgumentTypeError(f"Invalid exclude range: {part}. Use format: START1-END1,START2-END2")
    return ranges


def parse_residue_list(res_string):
    """Parse comma-separated residue numbers"""
    if not res_string:
        return None
    try:
        return [int(r.strip()) for r in res_string.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid residue list: {res_string}. Use comma-separated integers.")


def build_selection_string(chain=None, residue_range=None, exclude_ranges=None, base_selection='calpha'):
    """Build ProDy selection string from parameters"""
    parts = [base_selection]
    
    if chain:
        parts.append(f"chain {chain}")
    
    if residue_range:
        start, end = residue_range
        parts.append(f"resnum {start} to {end}")
    
    if exclude_ranges:
        for start, end in exclude_ranges:
            parts.append(f"not (resnum {start} to {end})")
    
    selection = ' and '.join(parts)
    return selection


def calculate_prs_matrix(anm):
    print("\nCalculating PRS matrix...")
    
    result = calcPerturbResponse(anm)

    # Handle different ProDy versions:
    if not isinstance(result, (tuple, list)):
        prs_matrix = result
        print("Detected old ProDy: calcPerturbResponse returned a single matrix.")
    else:
        prs_matrix = result[0]
        print(f"Detected ProDy version returning {len(result)} outputs.")
    
    print(f"PRS matrix shape: {prs_matrix.shape}")
    return prs_matrix


def calculate_effectiveness_sensitivity(prs_matrix):
    """
    Calculate effectiveness (ability to perturb others) and 
    sensitivity (ability to be perturbed by others)
    """
    effectiveness = np.sum(np.abs(prs_matrix), axis=1)
    sensitivity = np.sum(np.abs(prs_matrix), axis=0)
    return effectiveness, sensitivity


def identify_key_residues(effectiveness, sensitivity, atoms, top_n=10):
    """Identify top regulatory residues"""
    resnums = atoms.getResnums()
    resnames = atoms.getResnames()
    
    # Find top effective residues
    top_eff_idx = np.argsort(effectiveness)[-top_n:][::-1]
    
    # Find top sensitive residues
    top_sens_idx = np.argsort(sensitivity)[-top_n:][::-1]
    
    # Find residues that are both effective and sensitive (hubs)
    norm_eff = (effectiveness - effectiveness.min()) / (effectiveness.max() - effectiveness.min() + 1e-10)
    norm_sens = (sensitivity - sensitivity.min()) / (sensitivity.max() - sensitivity.min() + 1e-10)
    hub_score = norm_eff * norm_sens
    top_hub_idx = np.argsort(hub_score)[-top_n:][::-1]
    
    results = {
        'top_effective': [(resnums[i], resnames[i], effectiveness[i]) for i in top_eff_idx],
        'top_sensitive': [(resnums[i], resnames[i], sensitivity[i]) for i in top_sens_idx],
        'top_hubs': [(resnums[i], resnames[i], hub_score[i]) for i in top_hub_idx],
        'hub_score': hub_score
    }
    
    return results


def analyze_specific_pathways(prs_matrix, effector_res, sensor_res, atoms):
    """Analyze communication between specific effector and sensor residues"""
    if effector_res is None or sensor_res is None:
        return None
    
    resnums = atoms.getResnums()
    
    # Map residue numbers to indices
    effector_idx = []
    for r in effector_res:
        idx = np.where(resnums == r)[0]
        if len(idx) > 0:
            effector_idx.append(idx[0])
        else:
            print(f"Warning: Effector residue {r} not found in selection")
    
    sensor_idx = []
    for r in sensor_res:
        idx = np.where(resnums == r)[0]
        if len(idx) > 0:
            sensor_idx.append(idx[0])
        else:
            print(f"Warning: Sensor residue {r} not found in selection")
    
    if not effector_idx or not sensor_idx:
        print("Warning: Could not find all specified effector/sensor residues")
        return None
    
    # Extract submatrix
    pathway_matrix = prs_matrix[np.ix_(effector_idx, sensor_idx)]
    pathway_strength = np.mean(np.abs(pathway_matrix))
    
    return {
        'effector_idx': effector_idx,
        'sensor_idx': sensor_idx,
        'pathway_matrix': pathway_matrix,
        'pathway_strength': pathway_strength,
        'effector_res': [effector_res[i] for i in range(len(effector_idx))],
        'sensor_res': [sensor_res[i] for i in range(len(sensor_idx))]
    }


def plot_prs_heatmap(prs_matrix, atoms, output_dir, vmax_percentile=95):
    """Plot PRS matrix as heatmap"""
    print("\nGenerating PRS heatmap...")
    
    vmax = np.percentile(np.abs(prs_matrix), vmax_percentile)
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    im = ax.imshow(prs_matrix, cmap='RdBu_r', aspect='auto', 
                   vmin=-vmax, vmax=vmax, interpolation='nearest')
    
    ax.set_xlabel('Sensor Residue Index', fontsize=12, fontweight='bold')
    ax.set_ylabel('Effector Residue Index', fontsize=12, fontweight='bold')
    ax.set_title('PRS Matrix: Perturbation Response', fontsize=14, fontweight='bold')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Response Magnitude', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/prs_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/prs_heatmap.png")


def plot_effectiveness_sensitivity(effectiveness, sensitivity, atoms, output_dir, key_residues, residue_indices, chain_info):
    """Plot effectiveness and sensitivity profiles"""
    print("\nGenerating effectiveness/sensitivity plots...")
    
    resnums = atoms.getResnums()
    hub_score = key_residues['hub_score']
    
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))
    
    # Plot 1: Effectiveness
    axes[0].bar(resnums, effectiveness, color='steelblue', alpha=0.7)
    axes[0].set_xlabel('Residue Number', fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Effectiveness', fontsize=11, fontweight='bold')
    axes[0].set_title('Effectiveness (Ability to Perturb Others)', fontsize=12, fontweight='bold')
    axes[0].grid(axis='y', alpha=0.3)
    
    # Highlight top effective residues
    top_eff = key_residues['top_effective']
    for resnum, resname, score in top_eff[:5]:
        idx = np.where(resnums == resnum)[0]
        if len(idx) > 0:
            axes[0].bar(resnum, effectiveness[idx[0]], color='red', alpha=0.8)
    
    # Plot 2: Sensitivity
    axes[1].bar(resnums, sensitivity, color='darkgreen', alpha=0.7)
    axes[1].set_xlabel('Residue Number', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('Sensitivity', fontsize=11, fontweight='bold')
    axes[1].set_title('Sensitivity (Ability to be Perturbed)', fontsize=12, fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)
    
    # Highlight top sensitive residues
    top_sens = key_residues['top_sensitive']
    for resnum, resname, score in top_sens[:5]:
        idx = np.where(resnums == resnum)[0]
        if len(idx) > 0:
            axes[1].bar(resnum, sensitivity[idx[0]], color='orange', alpha=0.8)
    
    # Plot 3: Hub score
    axes[2].bar(resnums, hub_score, color='purple', alpha=0.7)
    axes[2].set_xlabel('Residue Number', fontsize=11, fontweight='bold')
    axes[2].set_ylabel('Hub Score', fontsize=11, fontweight='bold')
    axes[2].set_title('Hub Score (Both Effective and Sensitive)', fontsize=12, fontweight='bold')
    axes[2].grid(axis='y', alpha=0.3)
    
    # Highlight top hubs
    top_hubs = key_residues['top_hubs']
    for resnum, resname, score in top_hubs[:5]:
        idx = np.where(resnums == resnum)[0]
        if len(idx) > 0:
            axes[2].bar(resnum, hub_score[idx[0]], color='darkred', alpha=0.8)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/effectiveness_sensitivity.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/effectiveness_sensitivity.png")


def plot_pathway_analysis(pathway_data, output_dir):
    """Plot specific pathway analysis if provided"""
    if pathway_data is None:
        return
    
    print("\nGenerating pathway analysis plot...")
    
    pathway_matrix = pathway_data['pathway_matrix']
    effector_res = pathway_data['effector_res']
    sensor_res = pathway_data['sensor_res']
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    im = ax.imshow(pathway_matrix, cmap='RdBu_r', aspect='auto', interpolation='nearest')
    
    ax.set_xlabel('Sensor Residues', fontsize=12, fontweight='bold')
    ax.set_ylabel('Effector Residues', fontsize=12, fontweight='bold')
    ax.set_title('Effector-Sensor Pathway Strength', fontsize=14, fontweight='bold')
    
    ax.set_xticks(range(len(sensor_res)))
    ax.set_xticklabels(sensor_res)
    ax.set_yticks(range(len(effector_res)))
    ax.set_yticklabels(effector_res)
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Response Magnitude', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pathway_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/pathway_analysis.png")


def save_results_to_file(effectiveness, sensitivity, atoms, key_residues, 
                         pathway_data, output_dir, selection_info):
    """Save comprehensive results to text file"""
    print("\nSaving numerical results to file...")
    
    resnums = atoms.getResnums()
    resnames = atoms.getResnames()
    
    with open(f'{output_dir}/prs_results.txt', 'w') as f:
        f.write("PERTURBATION RESPONSE SCANNING (PRS) ANALYSIS RESULTS\n")
        f.write("="*70 + "\n\n")
        
        f.write("SELECTION DETAILS\n")
        f.write("-"*70 + "\n")
        f.write(f"Selection string: {selection_info}\n")
        f.write(f"Number of residues analyzed: {len(resnums)}\n")
        f.write(f"Residue range: {resnums[0]} to {resnums[-1]}\n\n")
        
        f.write("PRS identifies residues that regulate protein dynamics and\n")
        f.write("allosteric communication pathways.\n\n")
        
        # Top effective residues
        f.write("\nTOP 10 MOST EFFECTIVE RESIDUES (Best Perturbators)\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Rank':<6} {'Residue':<12} {'Effectiveness':<15}\n")
        f.write("-"*70 + "\n")
        for rank, (resnum, resname, score) in enumerate(key_residues['top_effective'], 1):
            f.write(f"{rank:<6} {resname}{resnum:<8} {score:<15.4f}\n")
        
        # Top sensitive residues
        f.write("\n\nTOP 10 MOST SENSITIVE RESIDUES (Best Responders)\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Rank':<6} {'Residue':<12} {'Sensitivity':<15}\n")
        f.write("-"*70 + "\n")
        for rank, (resnum, resname, score) in enumerate(key_residues['top_sensitive'], 1):
            f.write(f"{rank:<6} {resname}{resnum:<8} {score:<15.4f}\n")
        
        # Top hub residues
        f.write("\n\nTOP 10 HUB RESIDUES (Both Effective and Sensitive)\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Rank':<6} {'Residue':<12} {'Hub Score':<15}\n")
        f.write("-"*70 + "\n")
        for rank, (resnum, resname, score) in enumerate(key_residues['top_hubs'], 1):
            f.write(f"{rank:<6} {resname}{resnum:<8} {score:<15.4f}\n")
        
        # Pathway analysis if available
        if pathway_data is not None:
            f.write("\n\nSPECIFIC PATHWAY ANALYSIS\n")
            f.write("-"*70 + "\n")
            f.write(f"Effector residues: {pathway_data['effector_res']}\n")
            f.write(f"Sensor residues: {pathway_data['sensor_res']}\n")
            f.write(f"Average pathway strength: {pathway_data['pathway_strength']:.4f}\n\n")
            
            f.write("Detailed pathway matrix:\n")
            f.write(f"{'Effector':<12}")
            for sensor in pathway_data['sensor_res']:
                f.write(f"{sensor:<12}")
            f.write("\n")
            
            for i, eff in enumerate(pathway_data['effector_res']):
                f.write(f"{eff:<12}")
                for j in range(len(pathway_data['sensor_res'])):
                    f.write(f"{pathway_data['pathway_matrix'][i, j]:<12.4f}")
                f.write("\n")
        
        # Full profile
        f.write("\n\nFULL RESIDUE PROFILES\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Residue':<12} {'Effectiveness':<18} {'Sensitivity':<18}\n")
        f.write("-"*70 + "\n")
        for i, (resnum, resname) in enumerate(zip(resnums, resnames)):
            f.write(f"{resname}{resnum:<8} {effectiveness[i]:<18.4f} {sensitivity[i]:<18.4f}\n")
    
    print(f"Saved: {output_dir}/prs_results.txt")


def main():
    parser = argparse.ArgumentParser(
        description='ProDy Perturbation Response Scanning (PRS) Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Whole structure (default Calpha network)
  python3 prody_prs_analysis.py --pdb structure.pdb

  # Specific chain only
  python3 prody_prs_analysis.py --pdb structure.pdb --chain A

  # Specific residue range (e.g., construct region)
  python3 prody_prs_analysis.py --pdb structure.pdb --chain A --residue-range 20-150

  # Exclude flexible regions
  python3 prody_prs_analysis.py --pdb structure.pdb --chain A --exclude-residues 1-15,200-220

  # Custom selection (epitope, interface, etc.)
  python3 prody_prs_analysis.py --pdb structure.pdb --selection "chain A and resnum 25 to 120"

  # With pathway analysis
  python3 prody_prs_analysis.py --pdb structure.pdb --chain A \
    --effector-residues 45,67,89 --sensor-residues 120,145,167
        """)
    
    parser.add_argument('--pdb', required=True, help='Input PDB file')
    parser.add_argument('--chain', default=None, help='Chain ID to analyze (optional)')
    parser.add_argument('--residue-range', type=str, default=None,
                       help='Residue range to analyze, e.g., 20-150 (optional)')
    parser.add_argument('--exclude-residues', type=str, default=None,
                       help='Residue ranges to exclude, e.g., 1-15,200-220 (optional)')
    parser.add_argument('--selection', type=str, default=None,
                       help='Custom ProDy selection string (overrides other selection options)')
    parser.add_argument('--cutoff', type=float, default=15.0, help='ANM cutoff distance in Å (default: 15.0)')
    parser.add_argument('--nmodes', type=int, default=20, help='Number of modes for ANM (default: 20)')
    parser.add_argument('--effector-residues', type=str, default=None,
                       help='Comma-separated effector residue numbers (optional)')
    parser.add_argument('--sensor-residues', type=str, default=None,
                       help='Comma-separated sensor residue numbers (optional)')
    parser.add_argument('--top-n', type=int, default=10, help='Number of top residues to report (default: 10)')
    
    args = parser.parse_args()
    
    # Parse parameters
    residue_range = parse_residue_range(args.residue_range)
    exclude_ranges = parse_exclude_residues(args.exclude_residues)
    effector_res = parse_residue_list(args.effector_residues)
    sensor_res = parse_residue_list(args.sensor_residues)
    
    # Setup
    pdb_name = os.path.splitext(os.path.basename(args.pdb))[0]
    output_dir = setup_output_directory(pdb_name)
    
    # Build selection
    if args.selection:
        selection_string = args.selection
        print(f"Using custom selection: {selection_string}")
    else:
        selection_string = build_selection_string(
            chain=args.chain,
            residue_range=residue_range,
            exclude_ranges=exclude_ranges,
            base_selection='calpha'
        )
        print(f"Built selection: {selection_string}")
    
    # Load structure
    print(f"\nLoading PDB: {args.pdb}")
    structure = parsePDB(args.pdb)
    atoms = structure.select(selection_string)
    
    if atoms is None:
        print(f"Error: No atoms found with selection '{selection_string}'")
        sys.exit(1)
    
    print(f"Found {atoms.numAtoms()} atoms for analysis")
    
    # Analyze chain structure
    residue_indices, chain_info = analyze_chain_structure(atoms)
    print_chain_summary(chain_info)
    
    # Identify interface residues (if multi-chain)
    interface_residues, interface_pairs = identify_interface_residues(
        atoms, chain_info, distance_cutoff=10.0
    )
    
    resnums = atoms.getResnums()
    print(f"Residue range in selection: {resnums[0]} to {resnums[-1]}")
    
    # Build ANM
    print(f"\nBuilding ANM with {args.nmodes} modes (cutoff: {args.cutoff} Å)...")
    anm = ANM(pdb_name)
    anm.buildHessian(atoms, cutoff=args.cutoff)
    anm.calcModes(n_modes=args.nmodes)
    
    # Calculate PRS matrix
    prs_matrix = calculate_prs_matrix(anm)
    
    # Calculate effectiveness and sensitivity
    effectiveness, sensitivity = calculate_effectiveness_sensitivity(prs_matrix)
    
    # Identify key residues
    key_residues = identify_key_residues(effectiveness, sensitivity, atoms, args.top_n)
    
    # Analyze specific pathways if provided
    pathway_data = analyze_specific_pathways(prs_matrix, effector_res, sensor_res, atoms)
    
    # Report results
    print("\n" + "="*70)
    print("PRS ANALYSIS RESULTS")
    print("="*70)
    
    print("\nTop 5 Most Effective Residues (Best Perturbators):")
    for rank, (resnum, resname, score) in enumerate(key_residues['top_effective'][:5], 1):
        print(f"  {rank}. {resname}{resnum}: {score:.4f}")
    
    print("\nTop 5 Most Sensitive Residues (Best Responders):")
    for rank, (resnum, resname, score) in enumerate(key_residues['top_sensitive'][:5], 1):
        print(f"  {rank}. {resname}{resnum}: {score:.4f}")
    
    print("\nTop 5 Hub Residues (Both Effective and Sensitive):")
    for rank, (resnum, resname, score) in enumerate(key_residues['top_hubs'][:5], 1):
        print(f"  {rank}. {resname}{resnum}: {score:.4f}")
    
    if pathway_data is not None:
        print(f"\nEffector-Sensor Pathway Strength: {pathway_data['pathway_strength']:.4f}")
    
    print("\n" + "="*70)
    
    # Generate plots
    plot_prs_heatmap(prs_matrix, atoms, output_dir)
    plot_effectiveness_sensitivity(effectiveness, sensitivity, atoms, output_dir, key_residues, residue_indices, chain_info)
    plot_pathway_analysis(pathway_data, output_dir)
    
    # Save results
    save_results_to_file(effectiveness, sensitivity, atoms, key_residues, 
                        pathway_data, output_dir, selection_string)
    
    print("\nPRS analysis complete.")
    print(f"All results saved to: {output_dir}/")


if __name__ == '__main__':
    main()

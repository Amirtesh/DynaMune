#!/usr/bin/env python3
# version1.0
"""
ProDy Pocket Dynamics Analysis Script

This script analyzes the dynamic behavior of binding pockets using
conformational sampling from Normal Mode Analysis.

Usage:
  # Whole structure (default - will detect pockets automatically)
  python3 prody_pocket_dynamics.py --pdb structure.pdb

  # Specific chain (recommended for multi-chain structures)
  python3 prody_pocket_dynamics.py --pdb structure.pdb --chain A

  # Specific residue range (e.g., analyze pocket in one domain)
  python3 prody_pocket_dynamics.py --pdb structure.pdb --chain A --residue-range 50-250

  # Skip flexible tails to reduce noise
  python3 prody_pocket_dynamics.py --pdb structure.pdb --chain A \
    --residue-range 20-300 --pocket-threshold 3.0

  # Analyze specific pocket residues
  python3 prody_pocket_dynamics.py --pdb structure.pdb --chain A \
    --pocket-residues 45,67,89,120,145
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cdist, pdist, squareform

try:
    from prody import *  # type: ignore
except ImportError:
    print("Error: ProDy is not installed. Install with: pip install prody")
    sys.exit(1)


def setup_output_directory(pdb_name):
    """Create output directory for results"""
    output_dir = f"{pdb_name}_Pocket_Dynamics"
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


def parse_residue_range(range_str):
    """Parse residue range string like '50-250' into (start, end)"""
    if not range_str:
        return None
    try:
        start, end = map(int, range_str.split('-'))
        return (start, end)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid residue range: {range_str}. Use format: START-END")


def parse_residue_list(res_string):
    """Parse comma-separated residue numbers"""
    if not res_string:
        return None
    try:
        return [int(r.strip()) for r in res_string.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid residue list: {res_string}. Use comma-separated integers.")


def build_selection_string(chain=None, residue_range=None, base_selection='calpha'):
    """Build ProDy selection string from parameters"""
    parts = [base_selection]
    
    if chain:
        parts.append(f"chain {chain}")
    
    if residue_range:
        start, end = residue_range
        parts.append(f"resnum {start} to {end}")
    
    selection = ' and '.join(parts)
    return selection


def detect_pockets_simple(atoms, threshold=10.0):
    """
    Simple pocket detection based on spatial clustering
    Identifies regions with high local density (potential pockets)
    """
    print(f"\nDetecting potential pockets (threshold: {threshold} Å)...")
    
    coords = atoms.getCoords()
    n_atoms = len(coords)
    
    # Calculate distance matrix
    dist_matrix = squareform(pdist(coords))
    
    # Count neighbors within threshold
    neighbor_counts = np.sum(dist_matrix < threshold, axis=1)
    
    # Identify residues with high local density (potential pocket-forming)
    # Use percentile to identify top regions
    pocket_threshold = np.percentile(neighbor_counts, 75)
    pocket_residues = np.where(neighbor_counts >= pocket_threshold)[0]
    
    print(f"Identified {len(pocket_residues)} residues as potential pocket-forming")
    
    return pocket_residues, neighbor_counts


def get_pocket_residues(atoms, pocket_res_list):
    """Get indices of specific pocket residues"""
    resnums = atoms.getResnums()
    
    pocket_indices = []
    for res in pocket_res_list:
        idx = np.where(resnums == res)[0]
        if len(idx) > 0:
            pocket_indices.append(idx[0])
        else:
            print(f"Warning: Residue {res} not found in selection")
    
    if not pocket_indices:
        print("Error: No valid pocket residues found")
        return None
    
    print(f"Using {len(pocket_indices)} specified pocket residues")
    return np.array(pocket_indices)


def generate_conformers(anm, atoms, n_conformers=20, rmsd=1.0):
    """
    Generate conformational ensemble from ANM modes
    """
    print(f"\nGenerating {n_conformers} conformers (RMSD: {rmsd} Å)...")
    
    # Sample conformations along the most significant modes
    ensemble = sampleModes(anm[:10], atoms, n_confs=n_conformers, rmsd=rmsd)
    
    print(f"Generated ensemble with {ensemble.numConfs()} conformations")
    return ensemble


def analyze_pocket_volume_changes(ensemble, pocket_indices):
    """
    Analyze pocket volume changes across conformational ensemble
    Approximate volume by convex hull or distance-based measures
    """
    print("\nAnalyzing pocket volume changes...")
    
    n_confs = ensemble.numConfs()
    volumes = np.zeros(n_confs)
    
    for i in range(n_confs):
        coords = ensemble.getCoordsets(i)[pocket_indices]
        
        # Simple volume approximation: bounding box volume
        # More sophisticated: use ConvexHull, but this is faster
        min_coords = coords.min(axis=0)
        max_coords = coords.max(axis=0)
        dimensions = max_coords - min_coords
        volume = np.prod(dimensions)
        
        volumes[i] = volume
    
    volume_range = volumes.max() - volumes.min()
    volume_std = volumes.std()
    mean_volume = volumes.mean()
    
    print(f"Mean pocket volume: {mean_volume:.2f} ų")
    print(f"Volume range: {volume_range:.2f} ų ({volume_range/mean_volume*100:.1f}% of mean)")
    print(f"Volume std: {volume_std:.2f} ų")
    
    return volumes


def analyze_pocket_flexibility(ensemble, pocket_indices):
    print("\nAnalyzing pocket flexibility...")

    # Get raw coordinate sets: (n_confs, n_atoms, 3)
    coords = ensemble.getCoordsets()

    # === CRITICAL ===
    # Subset only pocket atoms
    pocket_coords = coords[:, pocket_indices, :]
    # pocket_coords shape: (n_confs, 101, 3)

    # Compute RMSF manually on pocket residue coordinates
    mean_coords = pocket_coords.mean(axis=0)
    rmsf = np.sqrt(((pocket_coords - mean_coords)**2).sum(axis=2).mean(axis=0))

    return rmsf



def analyze_pocket_shape_changes(ensemble, pocket_indices):
    """
    Analyze shape changes using principal component analysis of distances
    """
    print("\nAnalyzing pocket shape changes...")
    
    n_confs = ensemble.numConfs()
    n_pocket = len(pocket_indices)
    
    # Calculate distance matrices for each conformer
    distance_matrices = []
    for i in range(n_confs):
        coords = ensemble.getCoordsets(i)[pocket_indices]
        dist_matrix = squareform(pdist(coords))
        distance_matrices.append(dist_matrix.flatten())
    
    distance_matrices = np.array(distance_matrices)
    
    # Calculate variance in distances (shape flexibility)
    distance_variance = distance_matrices.var(axis=0)
    mean_shape_variance = distance_variance.mean()
    
    print(f"Mean shape variance: {mean_shape_variance:.4f} ų")
    
    return distance_variance, distance_matrices


def calculate_pocket_accessibility(atoms, pocket_indices, probe_radius=1.4):
    """
    Estimate pocket accessibility/openness
    Based on coordination number and local geometry
    """
    print("\nCalculating pocket accessibility...")
    
    coords = atoms.getCoords()
    pocket_coords = coords[pocket_indices]
    
    # Calculate center of pocket
    pocket_center = pocket_coords.mean(axis=0)
    
    # Calculate distances from center to pocket residues
    distances_to_center = np.linalg.norm(pocket_coords - pocket_center, axis=1)
    
    mean_radius = distances_to_center.mean()
    
    # Accessibility score based on pocket openness
    accessibility = mean_radius / len(pocket_indices)**0.5
    
    print(f"Pocket radius: {mean_radius:.2f} Å")
    print(f"Accessibility score: {accessibility:.3f}")
    
    return accessibility, mean_radius, pocket_center


def plot_pocket_volume_dynamics(volumes, output_dir):
    """Plot pocket volume changes across ensemble"""
    print("\nGenerating pocket volume plot...")
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot 1: Volume trajectory
    axes[0].plot(range(len(volumes)), volumes, 'o-', color='steelblue', linewidth=2, markersize=6)
    axes[0].axhline(y=volumes.mean(), color='red', linestyle='--', linewidth=2, label='Mean')
    axes[0].fill_between(range(len(volumes)), 
                         volumes.mean() - volumes.std(), 
                         volumes.mean() + volumes.std(),
                         alpha=0.3, color='red', label='±1 SD')
    axes[0].set_xlabel('Conformer Index', fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Pocket Volume (ų)', fontsize=11, fontweight='bold')
    axes[0].set_title('Pocket Volume Dynamics', fontsize=12, fontweight='bold')
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    
    # Plot 2: Volume distribution
    axes[1].hist(volumes, bins=20, color='darkgreen', alpha=0.7, edgecolor='black')
    axes[1].axvline(x=volumes.mean(), color='red', linestyle='--', linewidth=2, label='Mean')
    axes[1].set_xlabel('Pocket Volume (ų)', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('Frequency', fontsize=11, fontweight='bold')
    axes[1].set_title('Pocket Volume Distribution', fontsize=12, fontweight='bold')
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pocket_volume_dynamics.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/pocket_volume_dynamics.png")


def plot_pocket_flexibility(rmsf, atoms, pocket_indices, output_dir):
    """Plot RMSF for pocket residues"""
    print("\nGenerating pocket flexibility plot...")
    
    # Global residue numbers, subset to pocket residues only
    resnums = atoms.getResnums()[pocket_indices]
    resnames = atoms.getResnames()[pocket_indices]

    # Safety: ensure matching lengths
    if len(resnums) != len(rmsf):
        raise ValueError(f"Length mismatch: resnums={len(resnums)}, rmsf={len(rmsf)}")

    fig, ax = plt.subplots(figsize=(14, 6))

    # Color scale based on RMSF values
    colors = plt.cm.RdYlBu_r(rmsf / (rmsf.max() if rmsf.max() != 0 else 1))

    # Bar plot
    ax.bar(resnums, rmsf, color=colors, edgecolor='black', linewidth=0.5)

    # Labels and formatting
    ax.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
    ax.set_ylabel('RMSF (Å)', fontsize=12, fontweight='bold')
    ax.set_title('Pocket Residue Flexibility', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

    # Add mean line
    ax.axhline(y=rmsf.mean(), color='red', linestyle='--', linewidth=2, label='Mean RMSF')
    ax.legend()

    # Save
    plt.tight_layout()
    outfile = os.path.join(output_dir, 'pocket_flexibility.png')
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {outfile}")



def plot_pocket_shape_variance(distance_variance, pocket_indices, atoms, output_dir):
    """Plot shape variance matrix"""
    print("\nGenerating shape variance plot...")
    
    n_pocket = len(pocket_indices)
    resnums = atoms.getResnums()[pocket_indices]
    
    # Reshape to matrix
    variance_matrix = squareform(distance_variance[:n_pocket*(n_pocket-1)//2])
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    im = ax.imshow(variance_matrix, cmap='YlOrRd', aspect='auto', interpolation='nearest')
    
    ax.set_xlabel('Pocket Residue Index', fontsize=12, fontweight='bold')
    ax.set_ylabel('Pocket Residue Index', fontsize=12, fontweight='bold')
    ax.set_title('Pocket Shape Variance (Distance Fluctuations)', fontsize=14, fontweight='bold')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Distance Variance (ų)', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pocket_shape_variance.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/pocket_shape_variance.png")


def save_results_to_file(volumes, rmsf, accessibility, mean_radius, 
                         pocket_indices, atoms, output_dir, selection_info):
    """Save comprehensive results to text file"""
    print("\nSaving numerical results to file...")
    
    resnums = atoms.getResnums()[pocket_indices]
    resnames = atoms.getResnames()[pocket_indices]
    
    with open(f'{output_dir}/pocket_dynamics_results.txt', 'w') as f:
        f.write("POCKET DYNAMICS ANALYSIS RESULTS\n")
        f.write("="*70 + "\n\n")
        
        f.write("SELECTION DETAILS\n")
        f.write("-"*70 + "\n")
        f.write(f"Selection string: {selection_info}\n")
        f.write(f"Number of pocket residues: {len(pocket_indices)}\n\n")
        
        f.write("POCKET PROPERTIES\n")
        f.write("-"*70 + "\n")
        f.write(f"Mean pocket volume: {volumes.mean():.2f} ų\n")
        f.write(f"Volume range: {volumes.max() - volumes.min():.2f} ų\n")
        f.write(f"Volume std deviation: {volumes.std():.2f} ų\n")
        f.write(f"Volume coefficient of variation: {volumes.std()/volumes.mean()*100:.1f}%\n\n")
        
        f.write(f"Mean pocket RMSF: {rmsf.mean():.3f} Å\n")
        f.write(f"Max pocket RMSF: {rmsf.max():.3f} Å\n")
        f.write(f"Min pocket RMSF: {rmsf.min():.3f} Å\n\n")
        
        f.write(f"Pocket radius: {mean_radius:.2f} Å\n")
        f.write(f"Accessibility score: {accessibility:.3f}\n\n")
        
        f.write("\nPOCKET RESIDUES AND FLEXIBILITY\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Residue':<12} {'RMSF (Å)':<15} {'Flexibility Rank':<20}\n")
        f.write("-"*70 + "\n")
        
        # Sort by RMSF
        sorted_indices = np.argsort(rmsf)[::-1]
        for rank, idx in enumerate(sorted_indices, 1):
            resnum = resnums[idx]
            resname = resnames[idx]
            f.write(f"{resname}{resnum:<8} {rmsf[idx]:<15.3f} {rank:<20}\n")
        
        f.write("\n\nVOLUME TRAJECTORY\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Conformer':<12} {'Volume (ų)':<15}\n")
        f.write("-"*70 + "\n")
        for i, vol in enumerate(volumes):
            f.write(f"{i:<12} {vol:<15.2f}\n")
    
    print(f"Saved: {output_dir}/pocket_dynamics_results.txt")


def save_pymol_script(pocket_indices, rmsf, atoms, output_dir):
    """Generate PyMOL script to visualize pocket and flexibility"""
    resnums = atoms.getResnums()[pocket_indices]
    
    # Normalize RMSF for coloring
    norm_rmsf = (rmsf - rmsf.min()) / (rmsf.max() - rmsf.min())
    
    with open(f'{output_dir}/visualize_pocket.pml', 'w') as f:
        f.write("# PyMOL script to visualize pocket dynamics\n\n")
        
        f.write("# Select pocket residues\n")
        residue_list = '+'.join([str(r) for r in resnums])
        f.write(f"select pocket, resi {residue_list}\n\n")
        
        f.write("# Show pocket\n")
        f.write("show surface, pocket\n")
        f.write("set surface_color, marine, pocket\n")
        f.write("set transparency, 0.3, pocket\n\n")
        
        f.write("# Color by flexibility (RMSF)\n")
        f.write("# Blue = rigid, Red = flexible\n")
        for i, (resnum, flex) in enumerate(zip(resnums, norm_rmsf)):
            if flex < 0.33:
                color = "blue"
            elif flex < 0.67:
                color = "yellow"
            else:
                color = "red"
            f.write(f"color {color}, resi {resnum}\n")
        
        f.write("\n# Display settings\n")
        f.write("hide everything\n")
        f.write("show cartoon\n")
        f.write("show sticks, pocket\n")
        f.write("set cartoon_fancy_helices, 1\n")
        f.write("set stick_radius, 0.2\n")
        f.write("bg_color white\n")
        f.write("zoom pocket\n")
    
    print(f"Saved: {output_dir}/visualize_pocket.pml")


def main():
    parser = argparse.ArgumentParser(
        description='ProDy Pocket Dynamics Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect pockets (whole structure)
  python3 prody_pocket_dynamics.py --pdb structure.pdb

  # Specific chain (recommended for multi-chain)
  python3 prody_pocket_dynamics.py --pdb structure.pdb --chain A

  # Analyze one domain (skip flexible tails)
  python3 prody_pocket_dynamics.py --pdb structure.pdb --chain A --residue-range 50-250

  # Analyze specific pocket residues
  python3 prody_pocket_dynamics.py --pdb structure.pdb --chain A \
    --pocket-residues 45,67,89,120,145
        """)
    
    parser.add_argument('--pdb', required=True, help='Input PDB file')
    parser.add_argument('--chain', default=None, help='Chain ID to analyze (optional)')
    parser.add_argument('--residue-range', type=str, default=None,
                       help='Residue range to analyze, e.g., 50-250 (optional)')
    parser.add_argument('--selection', type=str, default=None,
                       help='Custom ProDy selection string (overrides other options)')
    parser.add_argument('--pocket-residues', type=str, default=None,
                       help='Comma-separated pocket residue numbers (optional, auto-detect if not provided)')
    parser.add_argument('--pocket-threshold', type=float, default=10.0,
                       help='Threshold for pocket detection in Å (default: 10.0)')
    parser.add_argument('--cutoff', type=float, default=15.0,
                       help='ANM cutoff distance in Å (default: 15.0)')
    parser.add_argument('--nmodes', type=int, default=20,
                       help='Number of modes for analysis (default: 20)')
    parser.add_argument('--conformers', type=int, default=20,
                       help='Number of conformers to generate (default: 20)')
    parser.add_argument('--rmsd', type=float, default=1.0,
                       help='RMSD for conformer generation in Å (default: 1.0)')
    
    args = parser.parse_args()
    
    # Parse parameters
    residue_range = parse_residue_range(args.residue_range)
    pocket_res_list = parse_residue_list(args.pocket_residues)
    
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
    
    resnums = atoms.getResnums()
    print(f"Residue range in selection: {resnums[0]} to {resnums[-1]}")
    
    # Determine pocket residues
    if pocket_res_list:
        pocket_indices = get_pocket_residues(atoms, pocket_res_list)
        if pocket_indices is None:
            sys.exit(1)
    else:
        pocket_indices, neighbor_counts = detect_pockets_simple(atoms, args.pocket_threshold)
    
    # Build ANM
    print(f"\nBuilding ANM with {args.nmodes} modes (cutoff: {args.cutoff} Å)...")
    anm = ANM(pdb_name)
    anm.buildHessian(atoms, cutoff=args.cutoff)
    anm.calcModes(n_modes=args.nmodes)
    
    # Generate conformational ensemble
    ensemble = generate_conformers(anm, atoms, n_conformers=args.conformers, rmsd=args.rmsd)
    ensemble.setAtoms(atoms)

    
    # Analyze pocket properties
    volumes = analyze_pocket_volume_changes(ensemble, pocket_indices)
    
    # Reset ensemble for RMSF calculation
    ensemble = generate_conformers(anm, atoms, n_conformers=args.conformers, rmsd=args.rmsd)
    ensemble.setAtoms(atoms)

    rmsf = analyze_pocket_flexibility(ensemble, pocket_indices)
    
    # Analyze shape changes
    distance_variance, distance_matrices = analyze_pocket_shape_changes(ensemble, pocket_indices)
    
    # Calculate accessibility
    accessibility, mean_radius, pocket_center = calculate_pocket_accessibility(
        atoms, pocket_indices, probe_radius=1.4
    )
    
    # Report results
    print("\n" + "="*70)
    print("POCKET DYNAMICS ANALYSIS RESULTS")
    print("="*70)
    
    print(f"\nPocket contains {len(pocket_indices)} residues")
    print(f"Mean volume: {volumes.mean():.2f} ų (range: {volumes.max()-volumes.min():.2f} ų)")
    print(f"Volume variation: {volumes.std()/volumes.mean()*100:.1f}%")
    print(f"\nMean pocket RMSF: {rmsf.mean():.3f} Å")
    print(f"Pocket radius: {mean_radius:.2f} Å")
    print(f"Accessibility score: {accessibility:.3f}")
    
    resnums_global = atoms.getResnums()        # full protein (342)
    resnames_global = atoms.getResnames()      # full protein (342)

    print("\nTop 5 most flexible pocket residues:")

# Local pocket indices sorted by RMSF (largest first)
    top_indices = rmsf.argsort()[::-1][:5]

    for rank, local_idx in enumerate(top_indices, 1):
         global_idx = pocket_indices[local_idx]   # map pocket → global index
         print(f"  {rank}. {resnames_global[global_idx]}{resnums_global[global_idx]}: RMSF = {rmsf[local_idx]:.3f} Å")
    
    print("\n" + "="*70)
    
    # Generate plots
    plot_pocket_volume_dynamics(volumes, output_dir)
    plot_pocket_flexibility(rmsf, atoms, pocket_indices, output_dir)
    plot_pocket_shape_variance(distance_variance, pocket_indices, atoms, output_dir)
    
    # Save results
    save_results_to_file(volumes, rmsf, accessibility, mean_radius,
                        pocket_indices, atoms, output_dir, selection_string)
    save_pymol_script(pocket_indices, rmsf, atoms, output_dir)
    
    print("\nPocket dynamics analysis complete.")
    print(f"All results saved to: {output_dir}/")
    print(f"\nTo visualize in PyMOL, run: pymol {args.pdb} {output_dir}/visualize_pocket.pml")


if __name__ == '__main__':
    main()

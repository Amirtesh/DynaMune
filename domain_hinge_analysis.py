#!/usr/bin/env python3
# version1.0
"""
ProDy Domain Detection and Hinge Analysis Script

This script identifies protein domains based on correlated dynamics and
detects hinge regions that facilitate domain movements.

Usage:
  # Whole structure (default)
  python3 prody_domain_hinge_analysis.py --pdb structure.pdb

  # Specific chain (e.g., multi-chain antigen)
  python3 prody_domain_hinge_analysis.py --pdb structure.pdb --chain A

  # Specific residue range
  python3 prody_domain_hinge_analysis.py --pdb structure.pdb --chain A --residue-range 1-300

  # Custom domain threshold
  python3 prody_domain_hinge_analysis.py --pdb structure.pdb --chain A --domain-threshold 0.6
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform

try:
    from prody import *  # type: ignore
except ImportError:
    print("Error: ProDy is not installed. Install with: pip install prody")
    sys.exit(1)


def setup_output_directory(pdb_name):
    """Create output directory for results"""
    output_dir = f"{pdb_name}_Domain_Hinge_Analysis"
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
    """Parse residue range string like '1-300' into (start, end)"""
    if not range_str:
        return None
    try:
        start, end = map(int, range_str.split('-'))
        return (start, end)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid residue range: {range_str}. Use format: START-END")


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


def calculate_cross_correlations(anm, n_modes=20):
    """Calculate cross-correlation matrix from ANM modes"""
    print("\nCalculating cross-correlations...")
    
    # Get the covariance matrix from the ANM modes
    cov_matrix = calcCovariance(anm[:n_modes])
    
    # Calculate cross-correlations
    cross_corr = calcCrossCorr(anm[:n_modes])
    
    print(f"Cross-correlation matrix shape: {cross_corr.shape}")
    return cross_corr


def identify_domains_from_correlations(cross_corr, threshold=0.6, method='average'):
    """
    Identify domains using hierarchical clustering on cross-correlations
    
    Parameters:
    - cross_corr: Cross-correlation matrix
    - threshold: Correlation threshold for domain separation (default: 0.6)
    - method: Linkage method for hierarchical clustering
    """
    print(f"\nIdentifying domains (threshold: {threshold})...")
    
    # Convert correlation to distance (1 - correlation)
    distance_matrix = 1 - np.abs(cross_corr)
    
    # Ensure diagonal is zero
    np.fill_diagonal(distance_matrix, 0)
    
    # Convert to condensed distance matrix for linkage
    condensed_dist = squareform(distance_matrix, checks=False)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_dist, method=method)
    
    # Cut the dendrogram to form clusters
    # Use threshold as height cutoff
    clusters = fcluster(linkage_matrix, t=1-threshold, criterion='distance')
    
    n_domains = len(np.unique(clusters))
    print(f"Found {n_domains} domains")
    
    return clusters, linkage_matrix, n_domains


def identify_hinge_residues(cross_corr, clusters, atoms, percentile=10):
    """
    Identify hinge residues based on low correlation with their domain
    and high correlation with other domains
    """
    print("\nIdentifying hinge residues...")
    
    n_residues = len(clusters)
    resnums = atoms.getResnums()
    
    hinge_scores = np.zeros(n_residues)
    
    for i in range(n_residues):
        domain_i = clusters[i]
        
        # Get residues in same domain
        same_domain = clusters == domain_i
        
        # Get residues in different domains
        diff_domain = clusters != domain_i
        
        # Calculate average correlation with same domain (excluding self)
        same_domain_corr = np.mean(np.abs(cross_corr[i, same_domain & (np.arange(n_residues) != i)]))
        
        # Calculate average correlation with different domains
        if np.any(diff_domain):
            diff_domain_corr = np.mean(np.abs(cross_corr[i, diff_domain]))
        else:
            diff_domain_corr = 0
        
        # Hinge score: low correlation with own domain, moderate with others
        # Lower same_domain_corr and moderate diff_domain_corr indicate hinge
        hinge_scores[i] = diff_domain_corr / (same_domain_corr + 0.01)
    
    # Identify top hinge residues
    hinge_threshold = np.percentile(hinge_scores, 100 - percentile)
    hinge_residues = np.where(hinge_scores >= hinge_threshold)[0]
    
    print(f"Identified {len(hinge_residues)} potential hinge residues")
    
    return hinge_scores, hinge_residues


def calculate_mechanical_stiffness(anm, atoms):
    """
    Calculate mechanical stiffness (inverse of flexibility) for each residue
    Hinge residues typically have lower stiffness
    """
    print("\nCalculating mechanical stiffness...")
    
    # Calculate square fluctuations (flexibility)
    sqflucts = calcSqFlucts(anm[:20])
    
    # Stiffness is inverse of flexibility
    stiffness = 1.0 / (sqflucts + 0.001)
    
    # Normalize
    stiffness = (stiffness - stiffness.min()) / (stiffness.max() - stiffness.min())
    
    return stiffness, sqflucts


def plot_cross_correlation_matrix(cross_corr, clusters, atoms, output_dir):
    """Plot cross-correlation matrix with domain boundaries"""
    print("\nGenerating cross-correlation heatmap...")
    
    resnums = atoms.getResnums()
    n_domains = len(np.unique(clusters))
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    im = ax.imshow(cross_corr, cmap='RdBu_r', aspect='auto', 
                   vmin=-1, vmax=1, interpolation='nearest')
    
    # Add domain boundaries
    domain_changes = np.where(np.diff(clusters) != 0)[0] + 1
    for change in domain_changes:
        ax.axhline(y=change-0.5, color='black', linewidth=2)
        ax.axvline(x=change-0.5, color='black', linewidth=2)
    
    ax.set_xlabel('Residue Index', fontsize=12, fontweight='bold')
    ax.set_ylabel('Residue Index', fontsize=12, fontweight='bold')
    ax.set_title(f'Cross-Correlation Matrix ({n_domains} Domains)', fontsize=14, fontweight='bold')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Cross-Correlation', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/cross_correlation_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/cross_correlation_matrix.png")


def plot_domain_assignment(clusters, atoms, output_dir):
    """Plot domain assignment along sequence"""
    print("\nGenerating domain assignment plot...")
    
    resnums = atoms.getResnums()
    n_domains = len(np.unique(clusters))
    
    # Create color map for domains
    colors = plt.cm.tab10(np.linspace(0, 1, n_domains))
    domain_colors = [colors[c-1] for c in clusters]
    
    fig, ax = plt.subplots(figsize=(14, 4))
    
    ax.bar(resnums, np.ones(len(resnums)), color=domain_colors, width=1.0, edgecolor='none')
    
    ax.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
    ax.set_ylabel('Domain', fontsize=12, fontweight='bold')
    ax.set_title(f'Domain Assignment ({n_domains} Domains)', fontsize=14, fontweight='bold')
    ax.set_ylim([0, 1.5])
    ax.set_yticks([])
    
    # Add domain labels
    for domain_id in np.unique(clusters):
        domain_residues = np.where(clusters == domain_id)[0]
        center = np.mean(resnums[domain_residues])
        ax.text(center, 0.5, f'D{domain_id}', ha='center', va='center', 
               fontsize=14, fontweight='bold', color='white',
               bbox=dict(boxstyle='round', facecolor='black', alpha=0.6))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/domain_assignment.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/domain_assignment.png")


def plot_hinge_analysis(hinge_scores, hinge_residues, stiffness, clusters, atoms, output_dir, residue_indices, chain_info):
    """Plot hinge analysis results"""
    print("\nGenerating hinge analysis plots...")
    
    resnums = atoms.getResnums()
    
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))
    
    # Plot 1: Hinge scores
    axes[0].plot(residue_indices, hinge_scores, color='steelblue', linewidth=2, zorder=2)
    axes[0].scatter(residue_indices[hinge_residues], hinge_scores[hinge_residues], 
                   color='red', s=50, zorder=5, label='Potential Hinges')
    axes[0].set_xlabel('Residue Index', fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Hinge Score', fontsize=11, fontweight='bold')
    axes[0].set_title('Hinge Score (Inter-domain / Intra-domain Correlation)', 
                     fontsize=12, fontweight='bold')
    axes[0].legend()
    axes[0].grid(alpha=0.3, zorder=1)
    add_chain_boundaries(axes[0], chain_info, y_fraction=0.92)
    
    # Plot 2: Mechanical stiffness
    axes[1].plot(residue_indices, stiffness, color='darkgreen', linewidth=2, zorder=2)
    axes[1].scatter(residue_indices[hinge_residues], stiffness[hinge_residues], 
                   color='red', s=50, zorder=5, label='Potential Hinges')
    axes[1].set_xlabel('Residue Index', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('Mechanical Stiffness', fontsize=11, fontweight='bold')
    axes[1].set_title('Mechanical Stiffness (Normalized)', fontsize=12, fontweight='bold')
    axes[1].legend()
    axes[1].grid(alpha=0.3, zorder=1)
    add_chain_boundaries(axes[1], chain_info, y_fraction=0.92)
    
    # Plot 3: Combined view with domains
    n_domains = len(np.unique(clusters))
    colors = plt.cm.tab10(np.linspace(0, 1, n_domains))
    
    for domain_id in np.unique(clusters):
        domain_mask = clusters == domain_id
        axes[2].fill_between(residue_indices, 0, 1, where=domain_mask, 
                           alpha=0.3, color=colors[domain_id-1], 
                           label=f'Domain {domain_id}')
    
    # Normalize hinge scores for visualization
    norm_hinge = (hinge_scores - hinge_scores.min()) / (hinge_scores.max() - hinge_scores.min())
    axes[2].plot(residue_indices, norm_hinge, color='black', linewidth=2, label='Hinge Score', zorder=2)
    axes[2].scatter(residue_indices[hinge_residues], norm_hinge[hinge_residues], 
                   color='red', s=100, zorder=5, marker='^', label='Hinges')
    
    axes[2].set_xlabel('Residue Index', fontsize=11, fontweight='bold')
    axes[2].set_ylabel('Normalized Score', fontsize=11, fontweight='bold')
    axes[2].set_title('Domains and Hinge Regions', fontsize=12, fontweight='bold')
    axes[2].legend(loc='upper right', ncol=2)
    axes[2].set_ylim([0, 1.2])
    add_chain_boundaries(axes[2], chain_info, y_fraction=1.15)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/hinge_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/hinge_analysis.png")


def plot_dendrogram(linkage_matrix, atoms, output_dir):
    """Plot hierarchical clustering dendrogram"""
    print("\nGenerating dendrogram...")
    
    fig, ax = plt.subplots(figsize=(16, 8))
    
    dendrogram(linkage_matrix, ax=ax, color_threshold=0.5)
    
    ax.set_xlabel('Residue Index', fontsize=12, fontweight='bold')
    ax.set_ylabel('Distance (1 - |Correlation|)', fontsize=12, fontweight='bold')
    ax.set_title('Hierarchical Clustering Dendrogram', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/dendrogram.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/dendrogram.png")


def save_results_to_file(clusters, hinge_scores, hinge_residues, stiffness, 
                         atoms, output_dir, selection_info, n_domains):
    """Save comprehensive results to text file"""
    print("\nSaving numerical results to file...")
    
    resnums = atoms.getResnums()
    resnames = atoms.getResnames()
    
    with open(f'{output_dir}/domain_hinge_results.txt', 'w') as f:
        f.write("DOMAIN DETECTION AND HINGE ANALYSIS RESULTS\n")
        f.write("="*70 + "\n\n")
        
        f.write("SELECTION DETAILS\n")
        f.write("-"*70 + "\n")
        f.write(f"Selection string: {selection_info}\n")
        f.write(f"Number of residues analyzed: {len(resnums)}\n")
        f.write(f"Residue range: {resnums[0]} to {resnums[-1]}\n")
        f.write(f"Number of domains identified: {n_domains}\n\n")
        
        # Domain assignments
        f.write("\nDOMAIN ASSIGNMENTS\n")
        f.write("-"*70 + "\n")
        for domain_id in np.unique(clusters):
            domain_residues = np.where(clusters == domain_id)[0]
            start_res = resnums[domain_residues[0]]
            end_res = resnums[domain_residues[-1]]
            f.write(f"\nDomain {domain_id}:\n")
            f.write(f"  Residue range: {start_res} - {end_res}\n")
            f.write(f"  Number of residues: {len(domain_residues)}\n")
            f.write(f"  Residues: {', '.join([f'{resnames[i]}{resnums[i]}' for i in domain_residues[:10]])}")
            if len(domain_residues) > 10:
                f.write(f" ... ({len(domain_residues)-10} more)")
            f.write("\n")
        
        # Hinge residues
        f.write("\n\nHINGE RESIDUES\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Residue':<12} {'Hinge Score':<15} {'Stiffness':<15} {'Domain':<10}\n")
        f.write("-"*70 + "\n")
        
        # Sort by hinge score
        sorted_hinges = sorted(hinge_residues, key=lambda x: hinge_scores[x], reverse=True)
        for idx in sorted_hinges:
            resnum = resnums[idx]
            resname = resnames[idx]
            f.write(f"{resname}{resnum:<8} {hinge_scores[idx]:<15.4f} "
                   f"{stiffness[idx]:<15.4f} {clusters[idx]:<10}\n")
        
        # Full profile
        f.write("\n\nFULL RESIDUE PROFILES\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Residue':<12} {'Domain':<10} {'Hinge Score':<15} {'Stiffness':<15}\n")
        f.write("-"*70 + "\n")
        for i, (resnum, resname) in enumerate(zip(resnums, resnames)):
            f.write(f"{resname}{resnum:<8} {clusters[i]:<10} "
                   f"{hinge_scores[i]:<15.4f} {stiffness[i]:<15.4f}\n")
    
    print(f"Saved: {output_dir}/domain_hinge_results.txt")
    
    # Save PyMOL script for visualization
    save_pymol_script(clusters, hinge_residues, atoms, output_dir)


def save_pymol_script(clusters, hinge_residues, atoms, output_dir):
    """Generate PyMOL script to visualize domains and hinges"""
    resnums = atoms.getResnums()
    n_domains = len(np.unique(clusters))
    
    with open(f'{output_dir}/visualize_domains.pml', 'w') as f:
        f.write("# PyMOL script to visualize domains and hinges\n\n")
        f.write("# Color domains\n")
        
        colors = ['red', 'green', 'blue', 'yellow', 'magenta', 'cyan', 
                 'orange', 'purple', 'lime', 'pink']
        
        for domain_id in np.unique(clusters):
            domain_residues = np.where(clusters == domain_id)[0]
            residue_list = '+'.join([str(resnums[i]) for i in domain_residues])
            color = colors[(domain_id - 1) % len(colors)]
            f.write(f"select domain_{domain_id}, resi {residue_list}\n")
            f.write(f"color {color}, domain_{domain_id}\n")
        
        f.write("\n# Highlight hinge residues\n")
        hinge_list = '+'.join([str(resnums[i]) for i in hinge_residues])
        f.write(f"select hinges, resi {hinge_list}\n")
        f.write("show spheres, hinges\n")
        f.write("color white, hinges\n")
        f.write("set sphere_scale, 0.5, hinges\n")
        
        f.write("\n# Display settings\n")
        f.write("hide everything\n")
        f.write("show cartoon\n")
        f.write("set cartoon_fancy_helices, 1\n")
        f.write("set ray_shadows, 0\n")
        f.write("bg_color white\n")
    
    print(f"Saved: {output_dir}/visualize_domains.pml")


def main():
    parser = argparse.ArgumentParser(
        description='ProDy Domain Detection and Hinge Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Whole structure (default)
  python3 prody_domain_hinge_analysis.py --pdb structure.pdb

  # Specific chain (for multi-chain antigens)
  python3 prody_domain_hinge_analysis.py --pdb structure.pdb --chain A

  # Specific residue range
  python3 prody_domain_hinge_analysis.py --pdb structure.pdb --chain A --residue-range 1-300

  # Custom domain threshold (higher = fewer domains)
  python3 prody_domain_hinge_analysis.py --pdb structure.pdb --chain A --domain-threshold 0.7
        """)
    
    parser.add_argument('--pdb', required=True, help='Input PDB file')
    parser.add_argument('--chain', default=None, help='Chain ID to analyze (optional)')
    parser.add_argument('--residue-range', type=str, default=None,
                       help='Residue range to analyze, e.g., 1-300 (optional)')
    parser.add_argument('--selection', type=str, default=None,
                       help='Custom ProDy selection string (overrides other options)')
    parser.add_argument('--cutoff', type=float, default=15.0, 
                       help='ANM cutoff distance in Å (default: 15.0)')
    parser.add_argument('--nmodes', type=int, default=20, 
                       help='Number of modes for analysis (default: 20)')
    parser.add_argument('--domain-threshold', type=float, default=0.6,
                       help='Correlation threshold for domain separation (default: 0.6, range: 0-1)')
    parser.add_argument('--hinge-percentile', type=float, default=10,
                       help='Percentile for hinge identification (default: 10)')
    
    args = parser.parse_args()
    
    # Parse parameters
    residue_range = parse_residue_range(args.residue_range)
    
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
    
    # Build ANM
    print(f"\nBuilding ANM with {args.nmodes} modes (cutoff: {args.cutoff} Å)...")
    anm = ANM(pdb_name)
    anm.buildHessian(atoms, cutoff=args.cutoff)
    anm.calcModes(n_modes=args.nmodes)
    
    # Calculate cross-correlations
    cross_corr = calculate_cross_correlations(anm, args.nmodes)
    
    # Identify domains
    clusters, linkage_matrix, n_domains = identify_domains_from_correlations(
        cross_corr, threshold=args.domain_threshold
    )
    
    # Identify hinge residues
    hinge_scores, hinge_residues = identify_hinge_residues(
        cross_corr, clusters, atoms, percentile=args.hinge_percentile
    )
    
    # Calculate mechanical stiffness
    stiffness, sqflucts = calculate_mechanical_stiffness(anm, atoms)
    
    # Report results
    print("\n" + "="*70)
    print("DOMAIN AND HINGE ANALYSIS RESULTS")
    print("="*70)
    print(f"\nNumber of domains identified: {n_domains}")
    
    resnums = atoms.getResnums()
    resnames = atoms.getResnames()
    
    print("\nDomain assignments:")
    for domain_id in np.unique(clusters):
        domain_residues = np.where(clusters == domain_id)[0]
        start_res = resnums[domain_residues[0]]
        end_res = resnums[domain_residues[-1]]
        print(f"  Domain {domain_id}: residues {start_res}-{end_res} ({len(domain_residues)} residues)")
    
    print(f"\nIdentified {len(hinge_residues)} hinge residues:")
    sorted_hinges = sorted(hinge_residues, key=lambda x: hinge_scores[x], reverse=True)
    for idx in sorted_hinges[:10]:
        resnum = resnums[idx]
        resname = resnames[idx]
        print(f"  {resname}{resnum}: hinge score = {hinge_scores[idx]:.4f}")
    
    if len(hinge_residues) > 10:
        print(f"  ... and {len(hinge_residues)-10} more (see output file)")
    
    print("\n" + "="*70)
    
    # Generate plots
    plot_cross_correlation_matrix(cross_corr, clusters, atoms, output_dir)
    plot_domain_assignment(clusters, atoms, output_dir)
    plot_hinge_analysis(hinge_scores, hinge_residues, stiffness, clusters, atoms, output_dir, residue_indices, chain_info)
    plot_dendrogram(linkage_matrix, atoms, output_dir)
    
    # Save results
    save_results_to_file(clusters, hinge_scores, hinge_residues, stiffness, 
                        atoms, output_dir, selection_string, n_domains)
    
    print("\nDomain and hinge analysis complete.")
    print(f"All results saved to: {output_dir}/")
    print(f"\nTo visualize in PyMOL, run: pymol {args.pdb} {output_dir}/visualize_domains.pml")


if __name__ == '__main__':
    main()

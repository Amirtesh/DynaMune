#!/usr/bin/env python3
# version1.0
"""
Comprehensive ProDy NMA & Ensemble PCA Analysis Script for Multi-Chain Complexes

Usage:
  # Basic NMA (like prody_nma.py)
  python3 prody_nma_pca_combined.py protein.pdb --method ANM

  # NMA with PCA on generated ensemble (the new feature)
  python3 prody_nma_pca_combined.py protein.pdb --method ANM --ensemble-pca

  # Full analysis with all features
  python3 prody_nma_pca_combined.py complex.pdb --method ANM --conformers 20 --ensemble-pca --allosteric --clustenmd
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform

try:
    from prody import *  # type: ignore
except ImportError:
    print("Error: ProDy is not installed. Install with: pip install prody")
    sys.exit(1)


def setup_output_directory(pdb_name):
    """Create output directory for results"""
    output_dir = f"{pdb_name}_NMA_PCA_results"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir


def analyze_chain_structure(calphas):
    """Analyze chain structure and create mapping"""
    chains = calphas.getChids()
    resnums = calphas.getResnums()
    resnames = calphas.getResnames()
    
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
    
    return residue_indices, chain_info, chains, resnums, resnames


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


def identify_interface_residues(calphas, chain_info, distance_cutoff=10.0):
    """Identify residues at the interface between chains"""
    if len(chain_info) < 2:
        return [], {}
    
    print(f"\nIdentifying interface residues (cutoff: {distance_cutoff} Å)...")
    
    coords = calphas.getCoords()
    chains = calphas.getChids()
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


def perform_anm_analysis(structure, calphas, n_modes=20, cutoff=15.0):
    """Perform ANM analysis"""
    print(f"\n{'='*60}")
    print("PERFORMING ANM ANALYSIS")
    print(f"{'='*60}")
    print(f"Number of Cα atoms: {calphas.numAtoms()}")
    print(f"Cutoff distance: {cutoff} Å")
    print(f"Number of modes to calculate: {n_modes}")
    
    anm = ANM(structure.getTitle())
    anm.buildHessian(calphas, cutoff=cutoff)
    anm.calcModes(n_modes=n_modes)
    
    print(f"\nANM model built successfully")
    print(f"Calculated {len(anm)} non-zero modes")
    print(f"First 5 eigenvalues: {anm.getEigvals()[:5].round(3)}")
    
    return anm


def perform_gnm_analysis(structure, calphas, n_modes=20, cutoff=10.0):
    """Perform GNM analysis"""
    print(f"\n{'='*60}")
    print("PERFORMING GNM ANALYSIS")
    print(f"{'='*60}")
    print(f"Number of Cα atoms: {calphas.numAtoms()}")
    print(f"Cutoff distance: {cutoff} Å")
    print(f"Number of modes to calculate: {n_modes}")
    
    gnm = GNM(structure.getTitle())
    gnm.buildKirchhoff(calphas, cutoff=cutoff)
    gnm.calcModes(n_modes=n_modes)
    
    print(f"\nGNM model built successfully")
    print(f"Calculated {len(gnm)} non-zero modes")
    print(f"First 5 eigenvalues: {gnm.getEigvals()[:5].round(3)}")
    
    return gnm

def perform_pca_analysis(ensemble, n_components=20):
    """Perform PCA using ProDy on centered ensemble"""
    print(f"\n{'='*60}")
    print("PERFORMING PCA ANALYSIS")
    print(f"{'='*60}")
    print(f"Number of conformations: {ensemble.numCoordsets()}")
    print(f"Number of atoms: {ensemble.numAtoms()}")
    print(f"Number of components: {n_components}")

    from prody import PCA

    # Apply PCA directly using ProDy
    pca = PCA(ensemble.getTitle() + '_pca')
    pca.buildCovariance(ensemble)
    pca.calcModes(n_modes=n_components)

    # Extract variance
    eigvals = pca.getEigvals()
    total_variance = eigvals.sum()
    variance_explained = eigvals / total_variance * 100
    cumulative_variance = np.cumsum(variance_explained)

    print(f"\nPCA completed successfully")
    print(f"First 5 eigenvalues: {eigvals[:5].round(3)}")

    print(f"\nVariance explained by first 5 PCs:")
    for i in range(min(5, len(variance_explained))):
        print(f"  PC{i+1}: {variance_explained[i]:.2f}% (cumulative: {cumulative_variance[i]:.2f}%)")

    return pca, variance_explained, cumulative_variance



def calculate_mode_overlap(model1, model2, n_modes=None):
    """Calculate overlap between two sets of normal modes"""
    if n_modes is None:
        n_modes = min(len(model1), len(model2), 20) # Limit to first 20 for speed
    
    overlap_matrix = np.zeros((n_modes, n_modes))
    
    for i in range(n_modes):
        mode1 = model1[i].getEigvec()
        for j in range(n_modes):
            mode2 = model2[j].getEigvec()
            overlap = np.abs(np.dot(mode1, mode2))
            overlap_matrix[i, j] = overlap
    
    return overlap_matrix


def calculate_cumulative_overlap(model1, model2, n_modes=None):
    """Calculate cumulative overlap between mode subspaces"""
    if n_modes is None:
        n_modes = min(len(model1), len(model2), 20)
    
    cumulative_overlap = []
    
    for n in range(1, n_modes + 1):
        # Using subspace overlap (trace of dot product of covariance matrices)
        # This is more stable and standard than summing squared overlaps
        subspace1 = model1[:n].getEigvecs()
        subspace2 = model2[:n].getEigvecs()
        
        # Calculate P.T * Q * Q.T * P
        dot_product = np.dot(subspace1.T, subspace2)
        overlap_matrix = np.dot(dot_product, dot_product.T)
        
        # Cumulative overlap is trace(overlap_matrix) / n
        overlap = np.trace(overlap_matrix) / n
        cumulative_overlap.append(overlap)
    
    return np.array(cumulative_overlap)


def perform_clustenmd_analysis(ensemble, model, output_dir):
    """Perform ClustENMD-style analysis"""
    print(f"\n{'='*60}")
    print("PERFORMING CLUSTENMD-STYLE ANALYSIS")
    print(f"{'='*60}")
    
    n_confs = ensemble.numCoordsets()
    n_modes = min(len(model), 20)
    
    # Project ensemble onto normal modes
    projections = np.zeros((n_confs, n_modes))
    ref_coords = ensemble.getCoords()
    
    for i in range(n_confs):
        coords = ensemble.getCoordsets(i)
        displacement = (coords - ref_coords).flatten()
        
        for j in range(n_modes):
            mode_vector = model[j].getEigvec()
            projections[i, j] = np.dot(displacement, mode_vector)
    
    # Cluster based on projections (first 3 modes)
    proj_subset = projections[:, :3]
    distances = pdist(proj_subset)
    linkage_matrix = linkage(distances, method='ward')
    
    # Get clusters
    n_clusters = min(5, n_confs // 2) # Ensure at least 2 confs per cluster
    if n_clusters < 2:
        n_clusters = 2
    clusters = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
    
    print(f"Projected {n_confs} conformations onto {n_modes} modes")
    print(f"Identified {n_clusters} clusters")
    
    # Plot projections
    fig = plt.figure(figsize=(16, 10))
    
    # 2D projections
    ax1 = plt.subplot(2, 3, 1)
    scatter = ax1.scatter(projections[:, 0], projections[:, 1], 
                         c=clusters, cmap='viridis', s=50, alpha=0.7)
    ax1.set_xlabel('Mode 1 Projection', fontweight='bold')
    ax1.set_ylabel('Mode 2 Projection', fontweight='bold')
    ax1.set_title('Mode 1 vs Mode 2', fontweight='bold')
    plt.colorbar(scatter, ax=ax1, label='Cluster')
    ax1.grid(True, alpha=0.3)
    
    ax2 = plt.subplot(2, 3, 2)
    scatter = ax2.scatter(projections[:, 1], projections[:, 2], 
                         c=clusters, cmap='viridis', s=50, alpha=0.7)
    ax2.set_xlabel('Mode 2 Projection', fontweight='bold')
    ax2.set_ylabel('Mode 3 Projection', fontweight='bold')
    ax2.set_title('Mode 2 vs Mode 3', fontweight='bold')
    plt.colorbar(scatter, ax=ax2, label='Cluster')
    ax2.grid(True, alpha=0.3)
    
    ax3 = plt.subplot(2, 3, 3)
    scatter = ax3.scatter(projections[:, 0], projections[:, 2], 
                         c=clusters, cmap='viridis', s=50, alpha=0.7)
    ax3.set_xlabel('Mode 1 Projection', fontweight='bold')
    ax3.set_ylabel('Mode 3 Projection', fontweight='bold')
    ax3.set_title('Mode 1 vs Mode 3', fontweight='bold')
    plt.colorbar(scatter, ax=ax3, label='Cluster')
    ax3.grid(True, alpha=0.3)
    
    # 3D projection
    try:
        from mpl_toolkits.mplot3d import Axes3D
        ax4 = plt.subplot(2, 3, 4, projection='3d')
        scatter = ax4.scatter(projections[:, 0], projections[:, 1], projections[:, 2],
                             c=clusters, cmap='viridis', s=50, alpha=0.7)
        ax4.set_xlabel('Mode 1', fontweight='bold')
        ax4.set_ylabel('Mode 2', fontweight='bold')
        ax4.set_zlabel('Mode 3', fontweight='bold')
        ax4.set_title('3D Projection Space', fontweight='bold')
    except ImportError:
        print("Note: mpl_toolkits.mplot3d not found. Skipping 3D plot.")
        ax4 = plt.subplot(2, 3, 4)
        ax4.text(0.5, 0.5, '3D plot skipped', ha='center', va='center')
    
    # Dendrogram
    ax5 = plt.subplot(2, 3, 5)
    dendrogram(linkage_matrix, ax=ax5, truncate_mode='lastp', p=12)
    ax5.set_xlabel('Conformation Index / Cluster', fontweight='bold')
    ax5.set_ylabel('Distance', fontweight='bold')
    ax5.set_title('Hierarchical Clustering', fontweight='bold')
    
    # Projection distribution
    ax6 = plt.subplot(2, 3, 6)
    for i in range(min(5, n_modes)):
        ax6.hist(projections[:, i], bins=20, alpha=0.5, label=f'Mode {i+1}')
    ax6.set_xlabel('Projection Value', fontweight='bold')
    ax6.set_ylabel('Frequency', fontweight='bold')
    ax6.set_title('Projection Distributions', fontweight='bold')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/clustenmd_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save cluster assignments
    with open(f'{output_dir}/cluster_assignments.txt', 'w') as f:
        f.write("Conformation\tCluster\tMode1_Proj\tMode2_Proj\tMode3_Proj\n")
        for i in range(n_confs):
            f.write(f"{i+1}\t{clusters[i]}\t{projections[i,0]:.4f}\t"
                   f"{projections[i,1]:.4f}\t{projections[i,2]:.4f}\n")
    
    print(f"Saved: {output_dir}/clustenmd_analysis.png")
    print(f"Saved: {output_dir}/cluster_assignments.txt")
    
    return projections, clusters


def analyze_allosteric_pathways(model, interface_residues, chain_info, output_dir):
    """Identify potential allosteric communication pathways"""
    print(f"\n{'='*60}")
    print("ANALYZING ALLOSTERIC PATHWAYS")
    print(f"{'='*60}")
    
    if not interface_residues:
        print("No interface residues - skipping allosteric analysis")
        return
    
    # Calculate cross-correlation
    try:
        cross_corr = calcCrossCorr(model)
    except:
        print("Could not calculate cross-correlation - skipping")
        return
    
    # Identify highly correlated residue pairs
    threshold = 0.7
    n_residues = cross_corr.shape[0]
    
    allosteric_pairs = []
    for i in range(n_residues):
        for j in range(i+1, n_residues):
            if abs(cross_corr[i, j]) > threshold:
                # Check if one is interface and one is not
                i_is_interface = i in interface_residues
                j_is_interface = j in interface_residues
                
                if i_is_interface != j_is_interface:
                    allosteric_pairs.append({
                        'res1': i,
                        'res2': j,
                        'correlation': cross_corr[i, j],
                        'distance': abs(i - j)
                    })
    
    # Sort by correlation strength
    allosteric_pairs.sort(key=lambda x: abs(x['correlation']), reverse=True)
    
    print(f"Found {len(allosteric_pairs)} potential allosteric connections")
    print(f"Top 5 allosteric pairs:")
    for pair in allosteric_pairs[:5]:
        print(f"  Residue {pair['res1']+1} <-> {pair['res2']+1}: "
              f"correlation = {pair['correlation']:.3f}")
    
    # Save results
    with open(f'{output_dir}/allosteric_pathways.txt', 'w') as f:
        f.write("ALLOSTERIC PATHWAY ANALYSIS\n")
        f.write("="*60 + "\n\n")
        f.write(f"Correlation threshold: {threshold}\n")
        f.write(f"Number of potential allosteric connections: {len(allosteric_pairs)}\n\n")
        f.write("Res1\tRes2\tCorrelation\tSequence_Distance\n")
        for pair in allosteric_pairs:
            f.write(f"{pair['res1']+1}\t{pair['res2']+1}\t"
                   f"{pair['correlation']:.4f}\t{pair['distance']}\n")
    
    print(f"Saved: {output_dir}/allosteric_pathways.txt")


def calculate_bfactors(model, calphas):
    """Calculate theoretical B-factors from normal modes"""
    print("\nCalculating theoretical B-factors...")
    msf = calcSqFlucts(model)
    bfactors = 8 * np.pi**2 / 3 * msf
    return msf, bfactors


def compare_with_experimental_bfactors(calphas, theoretical_bfactors, output_dir):
    """Compare theoretical B-factors with experimental values"""
    print("\nComparing with experimental B-factors...")
    
    experimental_bfactors = calphas.getBetas()
    
    if experimental_bfactors is None or np.all(experimental_bfactors == 0):
        print("No experimental B-factors available")
        return
    
    # Calculate correlation
    mask = ~np.isnan(experimental_bfactors) & ~np.isnan(theoretical_bfactors)
    if np.sum(mask) < 10:
        print("Insufficient data for comparison")
        return
    
    correlation = np.corrcoef(experimental_bfactors[mask], 
                              theoretical_bfactors[mask])[0, 1]
    
    print(f"Correlation with experimental B-factors: {correlation:.3f}")
    
    # Plot comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Scatter plot
    axes[0].scatter(experimental_bfactors[mask], theoretical_bfactors[mask], 
                   alpha=0.5, s=20)
    axes[0].plot([experimental_bfactors[mask].min(), experimental_bfactors[mask].max()],
                [experimental_bfactors[mask].min(), experimental_bfactors[mask].max()],
                'r--', linewidth=2, label='y=x')
    axes[0].set_xlabel('Experimental B-factors (Å²)', fontweight='bold')
    axes[0].set_ylabel('Theoretical B-factors (Å²)', fontweight='bold')
    axes[0].set_title(f'B-factor Correlation (r={correlation:.3f})', fontweight='bold')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Profile plot
    residue_indices = np.arange(1, len(calphas) + 1)
    axes[1].plot(residue_indices, experimental_bfactors, 'b-', 
                label='Experimental', linewidth=2, alpha=0.7)
    axes[1].plot(residue_indices, theoretical_bfactors, 'r-', 
                label='Theoretical', linewidth=2, alpha=0.7)
    axes[1].set_xlabel('Residue Index', fontweight='bold')
    axes[1].set_ylabel('B-factor (Å²)', fontweight='bold')
    axes[1].set_title('B-factor Profiles', fontweight='bold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/bfactor_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_dir}/bfactor_comparison.png")


def add_chain_boundaries(ax, chain_info, y_fraction=0.95):
    """Add chain boundary markers and labels to plot"""
    ylim = ax.get_ylim()
    y_pos = ylim[0] + (ylim[1] - ylim[0]) * y_fraction
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(chain_info)))
    
    for i, info in enumerate(chain_info):
        # Add vertical line at chain boundary (except first)
        if i > 0:
            ax.axvline(x=info['start_idx']+1, color='gray', linestyle='--', 
                      alpha=0.4, linewidth=1.5, zorder=1)
        
        # Add colored background for each chain
        ax.add_patch(Rectangle((info['start_idx']+0.5, ylim[0]), 
                               info['length'], ylim[1]-ylim[0],
                               facecolor=colors[i], alpha=0.1, zorder=0))
        
        # Add chain label
        mid_point = info['start_idx'] + info['length'] / 2
        ax.text(mid_point, y_pos, f"Chain {info['chain']}\n({info['length']} res)", 
               ha='center', va='top', fontsize=9, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.4', facecolor=colors[i], 
                        alpha=0.7, edgecolor='black', linewidth=1))


def plot_pca_analysis(pca, ensemble, variance_explained, cumulative_variance, 
                     output_dir, residue_indices, chain_info):
    """Create comprehensive PCA analysis plots"""
    print("\nGenerating PCA analysis plots...")
    
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)
    
    # 1. Variance explained
    ax1 = fig.add_subplot(gs[0, 0])
    n_components = len(variance_explained)
    ax1.bar(range(1, min(21, n_components+1)), 
           variance_explained[:20], alpha=0.7, color='steelblue')
    ax1.set_xlabel('Principal Component', fontweight='bold')
    ax1.set_ylabel('Variance Explained (%)', fontweight='bold')
    ax1.set_title('Scree Plot', fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # 2. Cumulative variance
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(range(1, min(21, n_components+1)), 
            cumulative_variance[:20], 'o-', linewidth=2, markersize=6)
    ax2.axhline(y=90, color='r', linestyle='--', label='90% threshold')
    ax2.set_xlabel('Number of Components', fontweight='bold')
    ax2.set_ylabel('Cumulative Variance (%)', fontweight='bold')
    ax2.set_title('Cumulative Variance Explained', fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. PC1 vs PC2 scatter
    ax3 = fig.add_subplot(gs[0, 2])
    projection = calcProjection(ensemble, pca[:2])
    ax3.scatter(projection[:, 0], projection[:, 1], alpha=0.6, s=30)
    ax3.set_xlabel(f'PC1 ({variance_explained[0]:.1f}%)', fontweight='bold')
    ax3.set_ylabel(f'PC2 ({variance_explained[1]:.1f}%)', fontweight='bold')
    ax3.set_title('PC1 vs PC2 Projection', fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # 4-6. PC shapes
    for i in range(3):
        ax = fig.add_subplot(gs[1, i])
        if i < len(pca):
            mode = pca[i]
            mode_array = mode.getArrayNx3()
            mode_magnitude = np.sqrt((mode_array**2).sum(axis=1))
            
            ax.plot(residue_indices, mode_magnitude, linewidth=2, zorder=2)
            add_chain_boundaries(ax, chain_info, y_fraction=0.92)
            ax.set_xlabel('Residue Index', fontweight='bold', fontsize=9)
            ax.set_ylabel('Magnitude', fontweight='bold', fontsize=9)
            ax.set_title(f'PC{i+1} ({variance_explained[i]:.1f}%)', 
                        fontweight='bold')
            ax.grid(True, alpha=0.3, zorder=1)
    
    # 7. Eigenvalue spectrum
    ax7 = fig.add_subplot(gs[2, 0])
    eigenvalues = pca.getEigvals()
    ax7.plot(range(1, min(21, len(eigenvalues)+1)), eigenvalues[:20], 
            'o-', linewidth=2, markersize=5)
    ax7.set_xlabel('PC Index', fontweight='bold')
    ax7.set_ylabel('Eigenvalue', fontweight='bold')
    ax7.set_title('PCA Eigenvalue Spectrum', fontweight='bold')
    ax7.set_yscale('log')
    ax7.grid(True, alpha=0.3)
    
    # 8. RMSF from PCA
    ax8 = fig.add_subplot(gs[2, 1])
    pca_msf = calcSqFlucts(pca)
    ax8.plot(residue_indices, pca_msf, linewidth=2, color='green', zorder=2)
    add_chain_boundaries(ax8, chain_info, y_fraction=0.92)
    ax8.set_xlabel('Residue Index', fontweight='bold')
    ax8.set_ylabel('Mean Square Fluctuation (Å²)', fontweight='bold')
    ax8.set_title('PCA-derived Fluctuations', fontweight='bold')
    ax8.grid(True, alpha=0.3, zorder=1)
    
    # 9. Cross-correlation from PCA
    ax9 = fig.add_subplot(gs[2, 2])
    try:
        pca_corr = calcCrossCorr(pca)
        im = ax9.imshow(pca_corr, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
        for info in chain_info[1:]:
            ax9.axhline(y=info['start_idx'], color='black', linewidth=1, alpha=0.7)
            ax9.axvline(x=info['start_idx'], color='black', linewidth=1, alpha=0.7)
        ax9.set_xlabel('Residue Index', fontweight='bold')
        ax9.set_ylabel('Residue Index', fontweight='bold')
        ax9.set_title('PCA Cross-Correlation', fontweight='bold')
        plt.colorbar(im, ax=ax9, label='Correlation')
    except:
        ax9.text(0.5, 0.5, 'Cross-correlation\ncalculation error',
                ha='center', va='center', transform=ax9.transAxes)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pca_comprehensive_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_dir}/pca_comprehensive_analysis.png")


def plot_nma_pca_comparison(anm, pca, output_dir, residue_indices, chain_info):
    """Compare ANM and PCA results"""
    print("\nGenerating ANM vs PCA comparison plots...")
    
    # Calculate overlaps
    n_modes = min(len(anm), len(pca), 10)
    overlap_matrix = calculate_mode_overlap(anm, pca, n_modes)
    cumulative_overlap = calculate_cumulative_overlap(anm, pca, n_modes)
    
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    # 1. Mode overlap matrix
    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.imshow(overlap_matrix, cmap='YlOrRd', vmin=0, vmax=1, aspect='auto')
    ax1.set_xlabel('PCA Mode', fontweight='bold')
    ax1.set_ylabel('ANM Mode', fontweight='bold')
    ax1.set_title('Mode Overlap Matrix', fontweight='bold')
    plt.colorbar(im1, ax=ax1, label='|Overlap|')
    
    # 2. Cumulative overlap
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(range(1, len(cumulative_overlap)+1), cumulative_overlap, 
            'o-', linewidth=2, markersize=6)
    ax2.axhline(y=0.5, color='r', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Number of Modes', fontweight='bold')
    ax2.set_ylabel('Cumulative Overlap', fontweight='bold')
    ax2.set_title('Subspace Overlap', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 1])
    
    # 3. Eigenvalue comparison
    ax3 = fig.add_subplot(gs[0, 2])
    anm_eigvals = anm.getEigvals()[:n_modes]
    pca_eigvals = pca.getEigvals()[:n_modes]
    x = range(1, n_modes+1)
    ax3.plot(x, anm_eigvals, 'o-', label='ANM', linewidth=2, markersize=6)
    ax3.plot(x, pca_eigvals, 's-', label='PCA', linewidth=2, markersize=6)
    ax3.set_xlabel('Mode/PC Index', fontweight='bold')
    ax3.set_ylabel('Eigenvalue', fontweight='bold')
    ax3.set_title('Eigenvalue Comparison', fontweight='bold')
    ax3.set_yscale('log')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Fluctuation comparison
    ax4 = fig.add_subplot(gs[1, :])
    anm_msf = calcSqFlucts(anm)
    pca_msf = calcSqFlucts(pca)
    ax4.plot(residue_indices, anm_msf, linewidth=2, label='ANM', alpha=0.8, zorder=2)
    ax4.plot(residue_indices, pca_msf, linewidth=2, label='PCA', alpha=0.8, zorder=2)
    add_chain_boundaries(ax4, chain_info, y_fraction=0.92)
    ax4.set_xlabel('Residue Index', fontweight='bold')
    ax4.set_ylabel('Mean Square Fluctuation (Å²)', fontweight='bold')
    ax4.set_title('ANM vs PCA Fluctuations', fontweight='bold')
    ax4.legend(fontsize=11)
    ax4.grid(True, alpha=0.3, zorder=1)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/anm_pca_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_dir}/anm_pca_comparison.png")
    print(f"Average mode overlap (first {n_modes} modes): {np.mean(overlap_matrix):.3f}")
    print(f"Cumulative overlap ({n_modes} modes): {cumulative_overlap[-1]:.3f}")
    
    # Save overlap data
    with open(f'{output_dir}/anm_pca_overlap.txt', 'w') as f:
        f.write("ANM vs PCA Mode Overlap Analysis\n")
        f.write("="*60 + "\n\n")
        f.write(f"Number of modes compared: {n_modes}\n")
        f.write(f"Average overlap: {np.mean(overlap_matrix):.4f}\n")
        f.write(f"Maximum overlap: {np.max(overlap_matrix):.4f}\n")
        f.write(f"Cumulative overlap: {cumulative_overlap[-1]:.4f}\n\n")
        f.write("Mode\tANM_Eigenvalue\tPCA_Eigenvalue\tMax_Overlap_with_PC\n")
        for i in range(n_modes):
            max_overlap = np.max(overlap_matrix[i, :])
            f.write(f"{i+1}\t{anm_eigvals[i]:.6f}\t{pca_eigvals[i]:.6f}\t{max_overlap:.4f}\n")
    
    print(f"Saved: {output_dir}/anm_pca_overlap.txt")


def plot_comprehensive_analysis(model, calphas, msf, bfactors, output_dir, 
                                method_name, residue_indices, chain_info, 
                                interface_residues):
    """Create comprehensive analysis plots with chain information
       (This is the detailed version from prody_nma.py)"""
    print("\nGenerating comprehensive analysis plots...")
    
    chains = calphas.getChids()
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(20, 14))
    gs = GridSpec(4, 3, figure=fig, hspace=0.35, wspace=0.3)
    
    # 1. Square Fluctuations with chain info
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(residue_indices, msf, linewidth=2, color='blue', zorder=2)
    
    # Highlight interface residues
    if interface_residues:
        ax1.scatter(residue_indices[interface_residues], msf[interface_residues],
                   color='red', s=30, alpha=0.6, label='Interface', zorder=3)
    
    add_chain_boundaries(ax1, chain_info)
    ax1.set_xlabel('Residue Index (continuous)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Mean Square Fluctuation (Å²)', fontsize=12, fontweight='bold')
    ax1.set_title(f'{method_name} - Square Fluctuations (Multi-Chain)', 
                 fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, zorder=1)
    ax1.set_xlim(residue_indices[0], residue_indices[-1])
    if interface_residues:
        ax1.legend(loc='upper right')
    
    # Save individual MSF plot
    fig_msf = plt.figure(figsize=(14, 7))
    ax_msf = plt.gca()
    ax_msf.plot(residue_indices, msf, linewidth=2, color='blue', zorder=2)
    if interface_residues:
        ax_msf.scatter(residue_indices[interface_residues], msf[interface_residues],
                      color='red', s=40, alpha=0.6, label='Interface residues', zorder=3)
    add_chain_boundaries(ax_msf, chain_info)
    plt.xlabel('Residue Index (continuous)', fontsize=12, fontweight='bold')
    plt.ylabel('Mean Square Fluctuation (Å²)', fontsize=12, fontweight='bold')
    plt.title(f'{method_name} - Square Fluctuations', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3, zorder=1)
    plt.xlim(residue_indices[0], residue_indices[-1])
    if interface_residues:
        plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/01_square_fluctuations.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Theoretical B-factors
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.plot(residue_indices, bfactors, 's-', alpha=0.7, markersize=2, color='red')
    ax2.set_xlabel('Residue Index', fontsize=10, fontweight='bold')
    ax2.set_ylabel('B-factor (Å²)', fontsize=10, fontweight='bold')
    ax2.set_title('Theoretical B-factors', fontsize=11, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Save individual B-factor plot
    fig_bf = plt.figure(figsize=(14, 7))
    ax_bf = plt.gca()
    ax_bf.plot(residue_indices, bfactors, 's-', alpha=0.7, markersize=3, color='red', zorder=2)
    add_chain_boundaries(ax_bf, chain_info, y_fraction=0.92)
    plt.xlabel('Residue Index', fontsize=12, fontweight='bold')
    plt.ylabel('B-factor (Å²)', fontsize=12, fontweight='bold')
    plt.title('Theoretical B-factors', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3, zorder=1)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/02_bfactors.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Eigenvalue spectrum
    ax3 = fig.add_subplot(gs[1, 1])
    eigenvalues = model.getEigvals()
    ax3.plot(range(1, len(eigenvalues)+1), eigenvalues, 'o-', linewidth=2, markersize=5)
    ax3.set_xlabel('Mode Index', fontsize=10, fontweight='bold')
    ax3.set_ylabel('Eigenvalue', fontsize=10, fontweight='bold')
    ax3.set_title('Eigenvalue Spectrum', fontsize=11, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale('log')
    
    # Save individual eigenvalue plot
    fig_eig = plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(eigenvalues)+1), eigenvalues, 'o-', linewidth=2, markersize=6)
    plt.xlabel('Mode Index', fontsize=12, fontweight='bold')
    plt.ylabel('Eigenvalue', fontsize=12, fontweight='bold')
    plt.title('Eigenvalue Spectrum', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/03_eigenvalues.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Mode shape for slowest mode
    ax4 = fig.add_subplot(gs[1, 2])
    if method_name == 'ANM':
        mode_shape = model[0].getArrayNx3()
        mode_magnitude = np.sqrt((mode_shape**2).sum(axis=1))
    else:
        mode_magnitude = np.abs(model[0].getArray())
    ax4.plot(residue_indices, mode_magnitude, linewidth=2, color='red')
    ax4.set_xlabel('Residue Index', fontsize=10, fontweight='bold')
    ax4.set_ylabel('Mode Magnitude', fontsize=10, fontweight='bold')
    ax4.set_title(f'Mode 1 (λ={model[0].getEigval():.3f})', fontsize=11, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    # Save individual mode 1 plot
    fig_mode1 = plt.figure(figsize=(14, 7))
    ax_m1 = plt.gca()
    ax_m1.plot(residue_indices, mode_magnitude, linewidth=2, color='red', zorder=2)
    add_chain_boundaries(ax_m1, chain_info)
    plt.xlabel('Residue Index', fontsize=12, fontweight='bold')
    plt.ylabel('Mode Magnitude', fontsize=12, fontweight='bold')
    plt.title(f'Mode 1 Shape (λ={model[0].getEigval():.3f})', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3, zorder=1)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/04_mode1_shape.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. Cross-correlation matrix
    ax5 = fig.add_subplot(gs[2, 0])
    if method_name == 'ANM':
        try:
            cross_corr = calcCrossCorr(model)
            im1 = ax5.imshow(cross_corr, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
            
            # Add chain boundary lines
            for info in chain_info[1:]:
                ax5.axhline(y=info['start_idx'], color='black', linewidth=1, alpha=0.7)
                ax5.axvline(x=info['start_idx'], color='black', linewidth=1, alpha=0.7)
            
            ax5.set_xlabel('Residue Index', fontsize=10, fontweight='bold')
            ax5.set_ylabel('Residue Index', fontsize=10, fontweight='bold')
            ax5.set_title('Cross-Correlation', fontsize=11, fontweight='bold')
            plt.colorbar(im1, ax=ax5, label='Correlation')
        except Exception as e:
            ax5.text(0.5, 0.5, f'Cross-correlation error',
                    ha='center', va='center', transform=ax5.transAxes)
            ax5.set_title('Cross-Correlation', fontsize=11, fontweight='bold')
    else:
        ax5.text(0.5, 0.5, 'Cross-correlation\nonly for ANM',
                ha='center', va='center', transform=ax5.transAxes)
        ax5.set_title('Cross-Correlation', fontsize=11, fontweight='bold')
    
    # Save individual cross-correlation
    if method_name == 'ANM':
        try:
            fig_cc = plt.figure(figsize=(12, 10))
            cross_corr_ind = calcCrossCorr(model)
            im = plt.imshow(cross_corr_ind, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
            
            # Add chain boundaries
            for info in chain_info[1:]:
                plt.axhline(y=info['start_idx'], color='black', linewidth=2, alpha=0.8)
                plt.axvline(x=info['start_idx'], color='black', linewidth=2, alpha=0.8)
            
            # Add chain labels
            for info in chain_info:
                mid = info['start_idx'] + info['length'] / 2
                plt.text(mid, -20, f"Chain {info['chain']}", ha='center', fontsize=10, fontweight='bold')
                plt.text(-20, mid, f"Chain {info['chain']}", ha='right', va='center', fontsize=10, fontweight='bold', rotation=90)
            
            plt.xlabel('Residue Index', fontsize=13, fontweight='bold')
            plt.ylabel('Residue Index', fontsize=13, fontweight='bold')
            plt.title('Cross-Correlation Matrix (Multi-Chain)', fontsize=15, fontweight='bold')
            cbar = plt.colorbar(im, label='Correlation')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/05_cross_correlation.png', dpi=300, bbox_inches='tight')
            plt.close()
        except:
            pass
    
    # 6. Covariance matrix
    ax6 = fig.add_subplot(gs[2, 1])
    if method_name == 'ANM':
        try:
            covariance = calcCovariance(model)
            n_residues = calphas.numAtoms()
            cov_matrix = np.zeros((n_residues, n_residues))
            for i in range(n_residues):
                for j in range(n_residues):
                    cov_matrix[i,j] = covariance[3*i:3*i+3, 3*j:3*j+3].trace()
            
            vmin = np.percentile(cov_matrix, 5)
            vmax = np.percentile(cov_matrix, 95)
            
            im2 = ax6.imshow(cov_matrix, cmap='plasma', aspect='auto', vmin=vmin, vmax=vmax)
            
            # Add chain boundaries
            for info in chain_info[1:]:
                ax6.axhline(y=info['start_idx'], color='white', linewidth=1, alpha=0.7)
                ax6.axvline(x=info['start_idx'], color='white', linewidth=1, alpha=0.7)
            
            ax6.set_xlabel('Residue Index', fontsize=10, fontweight='bold')
            ax6.set_ylabel('Residue Index', fontsize=10, fontweight='bold')
            ax6.set_title('Covariance Matrix', fontsize=11, fontweight='bold')
            plt.colorbar(im2, ax=ax6, label='Covariance')
        except Exception as e:
            ax6.text(0.5, 0.5, f'Covariance error',
                    ha='center', va='center', transform=ax6.transAxes)
            ax6.set_title('Covariance', fontsize=11, fontweight='bold')
    else:
        ax6.text(0.5, 0.5, 'Covariance\nonly for ANM',
                ha='center', va='center', transform=ax6.transAxes)
        ax6.set_title('Covariance', fontsize=11, fontweight='bold')
    
    # Save individual covariance
    if method_name == 'ANM':
        try:
            fig_cov = plt.figure(figsize=(12, 10))
            covariance_ind = calcCovariance(model)
            n_residues_ind = calphas.numAtoms()
            cov_matrix_ind = np.zeros((n_residues_ind, n_residues_ind))
            for i in range(n_residues_ind):
                for j in range(n_residues_ind):
                    cov_matrix_ind[i,j] = covariance_ind[3*i:3*i+3, 3*j:3*j+3].trace()
            
            vmin = np.percentile(cov_matrix_ind, 5)
            vmax = np.percentile(cov_matrix_ind, 95)
            
            im = plt.imshow(cov_matrix_ind, cmap='plasma', aspect='auto', vmin=vmin, vmax=vmax)
            
            # Add chain boundaries
            for info in chain_info[1:]:
                plt.axhline(y=info['start_idx'], color='white', linewidth=2, alpha=0.8)
                plt.axvline(x=info['start_idx'], color='white', linewidth=2, alpha=0.8)
            
            # Add chain labels
            for info in chain_info:
                mid = info['start_idx'] + info['length'] / 2
                plt.text(mid, -20, f"Chain {info['chain']}", ha='center', fontsize=10, fontweight='bold')
                plt.text(-20, mid, f"Chain {info['chain']}", ha='right', va='center', fontsize=10, fontweight='bold', rotation=90)
            
            plt.xlabel('Residue Index', fontsize=13, fontweight='bold')
            plt.ylabel('Residue Index', fontsize=13, fontweight='bold')
            plt.title('Covariance Matrix (Multi-Chain)', fontsize=15, fontweight='bold')
            cbar = plt.colorbar(im, label='Covariance')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/06_covariance.png', dpi=300, bbox_inches='tight')
            plt.close()
        except:
            pass
    
    # 7. Collectivity
    ax7 = fig.add_subplot(gs[2, 2])
    collectivity = calcCollectivity(model)
    ax7.plot(range(1, len(collectivity)+1), collectivity, 'o-', linewidth=2, markersize=5)
    ax7.set_xlabel('Mode Index', fontsize=10, fontweight='bold')
    ax7.set_ylabel('Collectivity', fontsize=10, fontweight='bold')
    ax7.set_title('Mode Collectivity', fontsize=11, fontweight='bold')
    ax7.grid(True, alpha=0.3)
    ax7.set_ylim([0, 1])
    
    # Save individual collectivity
    fig_coll = plt.figure(figsize=(10, 6))
    collectivity_ind = calcCollectivity(model)
    plt.plot(range(1, len(collectivity_ind)+1), collectivity_ind, 'o-', linewidth=2, markersize=6)
    plt.xlabel('Mode Index', fontsize=12, fontweight='bold')
    plt.ylabel('Collectivity', fontsize=12, fontweight='bold')
    plt.title('Mode Collectivity', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.ylim([0, 1])
    plt.tight_layout()
    plt.savefig(f'{output_dir}/07_collectivity.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 8. Per-chain fluctuation comparison
    ax8 = fig.add_subplot(gs[3, :])
    for info in chain_info:
        start = info['start_idx']
        end = info['end_idx']
        chain_indices = residue_indices[start:end] - start
        ax8.plot(chain_indices, msf[start:end], linewidth=2, 
                label=f"Chain {info['chain']}", alpha=0.8)
    ax8.set_xlabel('Residue Index (within chain)', fontsize=11, fontweight='bold')
    ax8.set_ylabel('MSF (Å²)', fontsize=11, fontweight='bold')
    ax8.set_title('Per-Chain Fluctuation Comparison', fontsize=12, fontweight='bold')
    ax8.legend()
    ax8.grid(True, alpha=0.3)
    
    # Save individual per-chain comparison
    fig_pc = plt.figure(figsize=(14, 7))
    for info in chain_info:
        start = info['start_idx']
        end = info['end_idx']
        chain_indices = residue_indices[start:end] - start
        plt.plot(chain_indices, msf[start:end], linewidth=2.5, 
                label=f"Chain {info['chain']} ({info['length']} residues)", alpha=0.8)
    plt.xlabel('Residue Index (within chain)', fontsize=12, fontweight='bold')
    plt.ylabel('Mean Square Fluctuation (Å²)', fontsize=12, fontweight='bold')
    plt.title('Per-Chain Fluctuation Comparison', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/08_per_chain_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save comprehensive plot
    plt.savefig(f'{output_dir}/comprehensive_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/comprehensive_analysis.png")
    print(f"Saved: Individual plots (01-08) in {output_dir}/")


def plot_mode_shapes(model, calphas, output_dir, method_name, residue_indices, 
                    chain_info, n_modes=None):
    """Plot shapes of the slowest modes with chain information
       (This is the detailed version from prody_nma.py)"""
    print("\nGenerating mode shape plots...")
    
    if n_modes is None:
        n_modes = len(model)
    n_plots = min(n_modes, len(model))
    
    # Combined plot (first 10 modes)
    n_combined = min(10, n_plots)
    fig, axes = plt.subplots(5, 2, figsize=(16, 18))
    axes = axes.flatten()
    
    for i in range(n_combined):
        mode = model[i]
        if method_name == 'ANM':
            mode_array = mode.getArrayNx3()
            mode_magnitude = np.sqrt((mode_array**2).sum(axis=1))
        else:
            mode_magnitude = np.abs(mode.getArray())
        
        axes[i].plot(residue_indices, mode_magnitude, linewidth=2, zorder=2)
        add_chain_boundaries(axes[i], chain_info, y_fraction=0.92)
        axes[i].set_xlabel('Residue Index', fontsize=9)
        axes[i].set_ylabel('Magnitude', fontsize=9)
        axes[i].set_title(f'Mode {i+1} (λ={mode.getEigval():.3f})', 
                         fontsize=10, fontweight='bold')
        axes[i].grid(True, alpha=0.3, zorder=1)
        
        if method_name == 'GNM':
            try:
                hinges = calcHinges(mode)
                if len(hinges) > 0:
                    axes[i].plot(residue_indices[hinges], mode_magnitude[hinges], 
                               'r*', markersize=8, label='Hinges', zorder=3)
                    axes[i].legend(fontsize=7)
            except:
                pass
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/mode_shapes_combined.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/mode_shapes_combined.png")
    
    # Individual mode plots
    print(f"Saving individual mode shape plots for all {n_plots} modes...")
    for i in range(n_plots):
        fig_single = plt.figure(figsize=(14, 7))
        ax = plt.gca()
        mode = model[i]
        
        if method_name == 'ANM':
            mode_array = mode.getArrayNx3()
            mode_magnitude = np.sqrt((mode_array**2).sum(axis=1))
        else:
            mode_magnitude = np.abs(mode.getArray())
        
        ax.plot(residue_indices, mode_magnitude, linewidth=2.5, zorder=2)
        add_chain_boundaries(ax, chain_info)
        plt.xlabel('Residue Index (continuous)', fontsize=12, fontweight='bold')
        plt.ylabel('Mode Magnitude', fontsize=12, fontweight='bold')
        plt.title(f'Mode {i+1} Shape (λ={mode.getEigval():.3f})', 
                 fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3, zorder=1)
        
        if method_name == 'GNM':
            try:
                hinges = calcHinges(mode)
                if len(hinges) > 0:
                    plt.plot(residue_indices[hinges], mode_magnitude[hinges], 
                           'r*', markersize=12, label='Hinge sites', zorder=3)
                    plt.legend(fontsize=10)
            except:
                pass
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/mode_{i+1:02d}_shape.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"Saved: Individual mode shapes (mode_01 to mode_{n_plots:02d}) in {output_dir}/")


def analyze_interface_dynamics(msf, interface_residues, interface_pairs, chain_info, output_dir):
    """Analyze dynamics at the interface with statistical testing
       (This is the detailed version from enhanced_prody.py)"""
    if not interface_residues:
        return
    
    print("\nAnalyzing interface dynamics...")
    
    interface_msf = msf[interface_residues]
    non_interface_msf = np.delete(msf, interface_residues)
    
    # Perform statistical test
    t_stat, p_value = stats.ttest_ind(interface_msf, non_interface_msf)
    
    stats_dict = {
        'interface_mean': np.mean(interface_msf),
        'interface_std': np.std(interface_msf),
        'interface_median': np.median(interface_msf),
        'non_interface_mean': np.mean(non_interface_msf),
        'non_interface_std': np.std(non_interface_msf),
        'non_interface_median': np.median(non_interface_msf),
        'ratio': np.mean(interface_msf) / np.mean(non_interface_msf),
        't_statistic': t_stat,
        'p_value': p_value
    }
    
    with open(f'{output_dir}/interface_analysis.txt', 'w') as f:
        f.write("INTERFACE DYNAMICS ANALYSIS\n")
        f.write("="*60 + "\n\n")
        f.write(f"Number of interface residues: {len(interface_residues)}\n")
        f.write(f"Total residues: {len(msf)}\n")
        f.write(f"Interface percentage: {100*len(interface_residues)/len(msf):.1f}%\n\n")
        f.write(f"Interface MSF: {stats_dict['interface_mean']:.4f} ± {stats_dict['interface_std']:.4f} Å²\n")
        f.write(f"Interface median: {stats_dict['interface_median']:.4f} Å²\n")
        f.write(f"Non-interface MSF: {stats_dict['non_interface_mean']:.4f} ± {stats_dict['non_interface_std']:.4f} Å²\n")
        f.write(f"Non-interface median: {stats_dict['non_interface_median']:.4f} Å²\n")
        f.write(f"Interface/Non-interface ratio: {stats_dict['ratio']:.2f}\n\n")
        f.write(f"Statistical test (t-test):\n")
        f.write(f"  t-statistic: {stats_dict['t_statistic']:.4f}\n")
        f.write(f"  p-value: {stats_dict['p_value']:.4e}\n")
        f.write(f"  Significant: {'Yes' if stats_dict['p_value'] < 0.05 else 'No'} (α=0.05)\n\n")
        
        if stats_dict['ratio'] < 0.8:
            f.write("Interpretation: Interface is MORE RIGID than bulk protein\n")
            f.write("  → Suggests stable binding interface\n")
        elif stats_dict['ratio'] > 1.2:
            f.write("Interpretation: Interface is MORE FLEXIBLE than bulk protein\n")
            f.write("  → Suggests dynamic/adaptable binding\n")
        else:
            f.write("Interpretation: Interface has SIMILAR flexibility to bulk protein\n")
            f.write("  → Suggests balanced rigidity/flexibility\n")
    
    print(f"Saved: {output_dir}/interface_analysis.txt")
    print(f"  Interface MSF: {stats_dict['interface_mean']:.4f} Å² (p={stats_dict['p_value']:.4e})")


def save_output_files(model, calphas, msf, bfactors, output_dir, method_name, 
                     residue_indices, chains, resnames, resnums, interface_residues):
    """Save analysis results to files with chain information
       (This is the detailed version from prody_nma.py)"""
    print("\nSaving output files...")
    
    # Save model
    saveModel(model, f'{output_dir}/{method_name}_model')
    print(f"Saved: {output_dir}/{method_name}_model.{method_name.lower()}.npz")
    
    # Save NMD file
    writeNMD(f'{output_dir}/{method_name}_modes.nmd', model, calphas)
    print(f"Saved: {output_dir}/{method_name}_modes.nmd")
    
    # Save data with chain information
    with open(f'{output_dir}/square_fluctuations.txt', 'w') as f:
        f.write("Index\tChain\tResName\tResNum\tMSF\tBfactor_theory\tInterface\n")
        for i, (chain, resname, resnum) in enumerate(zip(chains, resnames, resnums)):
            is_interface = "Yes" if i in interface_residues else "No"
            f.write(f"{residue_indices[i]}\t{chain}\t{resname}\t{resnum}\t"
                   f"{msf[i]:.4f}\t{bfactors[i]:.4f}\t{is_interface}\n")
    print(f"Saved: {output_dir}/square_fluctuations.txt")
    
    # Eigenvalues
    with open(f'{output_dir}/eigenvalues.txt', 'w') as f:
        f.write("Mode\tEigenvalue\tCollectivity\n")
        collectivity = calcCollectivity(model)
        for i, (eigval, coll) in enumerate(zip(model.getEigvals(), collectivity)):
            f.write(f"{i+1}\t{eigval:.6f}\t{coll:.6f}\n")
    print(f"Saved: {output_dir}/eigenvalues.txt")
    
    # Cross-correlation
    if method_name == 'ANM':
        try:
            cross_corr = calcCrossCorr(model)
            np.savetxt(f'{output_dir}/cross_correlation.txt', cross_corr, fmt='%.4f')
            print(f"Saved: {output_dir}/cross_correlation.txt")
        except:
            print("Cross-correlation matrix not saved (calculation error)")
    
    # Save PDB with theoretical B-factors
    writePDB(f'{output_dir}/{method_name}_bfactors.pdb', calphas, beta=bfactors)
    print(f"Saved: {output_dir}/{method_name}_bfactors.pdb")
    
    # Save interface residues
    if interface_residues:
        with open(f'{output_dir}/interface_residues.txt', 'w') as f:
            f.write("Index\tChain\tResName\tResNum\tMSF\n")
            for idx in interface_residues:
                f.write(f"{residue_indices[idx]}\t{chains[idx]}\t{resnames[idx]}\t"
                       f"{resnums[idx]}\t{msf[idx]:.4f}\n")
        print(f"Saved: {output_dir}/interface_residues.txt")


def generate_conformers(model, calphas, output_dir, method_name, n_conformers=10):
    """Generate sample conformers along slowest modes
       (This is the detailed version from prody_nma.py)"""
    if method_name != 'ANM':
        print(f"\nConformer generation only available for ANM")
        return None
    
    print(f"\nGenerating {n_conformers} sample conformers...")
    
    try:
        ensemble = sampleModes(model[:3], calphas, n_confs=n_conformers, rmsd=2.0)
        
        if ensemble is None:
            print("Warning: No conformers generated")
            return None
        
        print(f"Generated {ensemble.numCoordsets()} conformers")
        
        # Save individual conformers
        for i in range(ensemble.numCoordsets()):
            coords = ensemble.getCoordsets(i)
            conf_atoms = calphas.copy()
            conf_atoms.setCoords(coords)
            writePDB(f'{output_dir}/conformer_{i+1:02d}.pdb', conf_atoms)
        
        print(f"Saved: Individual conformers (conformer_01 to conformer_{ensemble.numCoordsets():02d})")
        
        # Save ensemble file
        try:
            with open(f'{output_dir}/conformers_ensemble.pdb', 'w') as f:
                for i in range(ensemble.numCoordsets()):
                    coords = ensemble.getCoordsets(i)
                    conf_atoms = calphas.copy()
                    conf_atoms.setCoords(coords)
                    
                    f.write(f"MODEL     {i+1:4d}\n")
                    for j, atom in enumerate(conf_atoms):
                        f.write(f"ATOM  {j+1:5d}  CA  {atom.getResname():3s} "
                               f"{atom.getChid()}{atom.getResnum():4d}    "
                               f"{coords[j][0]:8.3f}{coords[j][1]:8.3f}{coords[j][2]:8.3f}"
                               f"  1.00  0.00           C  \n")
                    f.write("ENDMDL\n")
            print(f"Saved: {output_dir}/conformers_ensemble.pdb")
        except Exception as e:
            print(f"Note: Ensemble PDB not created: {str(e)}")
        
        return ensemble
        
    except Exception as e:
        print(f"Warning: Could not generate conformers: {str(e)}")
        return None


def print_summary(model, calphas, msf, bfactors, method_name, chain_info, interface_residues):
    """Print analysis summary"""
    print(f"\n{'='*60}")
    print(f"{method_name} ANALYSIS SUMMARY")
    print(f"{'='*60}")
    print(f"Structure: {model.getTitle()}")
    print(f"Total Cα atoms: {calphas.numAtoms()}")
    print(f"Number of chains: {len(chain_info)}")
    print(f"Number of modes: {len(model)}")
    print(f"Cutoff distance: {model.getCutoff()} Å")
    
    print(f"\nChain composition:")
    for info in chain_info:
        print(f"  Chain {info['chain']}: {info['length']} residues")
    
    if interface_residues:
        print(f"\nInterface residues: {len(interface_residues)} "
              f"({100*len(interface_residues)/len(msf):.1f}%)")
    
    print(f"\nEigenvalues (first 5 modes):")
    for i, eigval in enumerate(model.getEigvals()[:5]):
        print(f"  Mode {i+1}: {eigval:.6f}")
    
    print(f"\nFluctuation statistics:")
    print(f"  Mean MSF: {np.mean(msf):.4f} Å²")
    print(f"  Std MSF: {np.std(msf):.4f} Å²")
    print(f"  Max MSF: {np.max(msf):.4f} Å² (index {np.argmax(msf)+1})")
    print(f"  Min MSF: {np.min(msf):.4f} Å² (index {np.argmin(msf)+1})")
    
    print(f"\n{'='*60}")


def main():
    parser = argparse.ArgumentParser(
        description='Enhanced ProDy NMA with Ensemble PCA Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 prody_nma.py protein.pdb
  python3 prody_nma.py protein.pdb --method ANM --nmodes 20 --conformers 20
  python3 prody_nma.py protein.pdb --method ANM --ensemble-pca --allosteric --cutoff 18 --gamma 0.5
        """)

    # Core inputs
    parser.add_argument('pdb_file', help='Input PDB file')
    parser.add_argument('--method', choices=['ANM', 'GNM'], default='ANM')
    parser.add_argument('--nmodes', type=int, default=20)
    parser.add_argument('--selection', default='calpha')
    parser.add_argument('--conformers', type=int, default=10)
    parser.add_argument('--interface-cutoff', type=float, default=10.0)
    parser.add_argument('--ensemble-pca', action='store_true')
    parser.add_argument('--allosteric', action='store_true')
    parser.add_argument('--clustenmd', action='store_true')
    parser.add_argument('--cutoff', type=float, default=None)
    parser.add_argument('--gamma', type=float, default=1.0)

    args = parser.parse_args()

    # Validate PDB file
    if not os.path.exists(args.pdb_file):
        print(f"Error: File '{args.pdb_file}' not found!")
        sys.exit(1)

    # Resolve cutoff
    cutoff = args.cutoff if args.cutoff is not None else (15.0 if args.method == 'ANM' else 10.0)
    gamma = args.gamma

    print(f"\n{'='*60}")
    print("Enhanced ProDy NMA with Ensemble PCA Analysis")
    print(f"{'='*60}")
    print(f"Input File: {args.pdb_file}")
    print(f"NMA Method: {args.method}")
    print(f"Cutoff: {cutoff} Å | Gamma: {gamma}")
    print(f"Modes: {args.nmodes} | Conformers: {args.conformers}")
    print(f"Selection: {args.selection}")
    print(f"Ensemble PCA: {args.ensemble_pca} | Allosteric: {args.allosteric} | ClustENMD: {args.clustenmd}")
    print(f"{'='*60}")

    # Load structure
    structure = parsePDB(args.pdb_file)
    calphas = structure.select(args.selection)
    if calphas is None or calphas.numAtoms() == 0:
        print(f"Error: No atoms found with selection '{args.selection}'!")
        sys.exit(1)

    coords = calphas.getCoords()

    pdb_name = os.path.splitext(os.path.basename(args.pdb_file))[0]
    output_dir = setup_output_directory(pdb_name)

    residue_indices, chain_info, chains, resnums, resnames = analyze_chain_structure(calphas)
    print_chain_summary(chain_info)

    interface_residues, interface_pairs = identify_interface_residues(
        calphas, chain_info, args.interface_cutoff)

    # NMA computation
    if args.method == 'ANM':
        model = ANM('anm_model')
        model.buildHessian(coords, cutoff=cutoff, gamma=gamma)
    else:
        model = GNM('gnm_model')
        model.buildKirchhoff(coords, cutoff=cutoff)

    model.calcModes(n_modes=args.nmodes)

    # Fluctuation and plot generation
    msf, bfactors = calculate_bfactors(model, calphas)
    plot_comprehensive_analysis(model, calphas, msf, bfactors, output_dir,
                                args.method, residue_indices, chain_info, interface_residues)
    plot_mode_shapes(model, calphas, output_dir, args.method, residue_indices, chain_info)
    compare_with_experimental_bfactors(calphas, bfactors, output_dir)
    analyze_interface_dynamics(msf, interface_residues, interface_pairs, chain_info, output_dir)

    # Allosteric pathway analysis
    if args.allosteric:
        if args.method != 'ANM':
            print("\nWARNING: Allosteric analysis only available for ANM. Skipping.")
        else:
            analyze_allosteric_pathways(model, interface_residues, chain_info, output_dir)

    # Generate conformers
    ensemble = generate_conformers(model, calphas, output_dir, args.method, args.conformers) \
               if args.method == 'ANM' else None

    # Ensemble PCA analysis
    if args.ensemble_pca and ensemble is not None:
        print(f"\n{'=' * 60}")
        print("ENSEMBLE PCA ANALYSIS (CENTERED)")
        print(f"{'=' * 60}")

        import numpy as np
        from prody import Ensemble, PCA

        # Convert to NumPy coordinates
        ensemble_coords = np.array([conf.getCoords() for conf in ensemble])

        # Center the ensemble
        ensemble_coords_centered = ensemble_coords - ensemble_coords.mean(axis=0)

        # Rebuild as ProDy ensemble
        centered_ensemble = Ensemble('centered_ensemble')
        centered_ensemble.setCoords(ensemble_coords_centered[0])
        for coords in ensemble_coords_centered:
            centered_ensemble.addCoordset(coords)

        # Perform PCA using ProDy
        pca, variance_explained, cumulative_variance = perform_pca_analysis(
            centered_ensemble, args.nmodes
        )

        # PCA visualization and comparison
        plot_pca_analysis(
            pca, centered_ensemble, variance_explained, cumulative_variance,
            output_dir, residue_indices, chain_info
        )
        plot_nma_pca_comparison(
            model, pca, output_dir, residue_indices, chain_info
        )

        # ClustENMD
        if args.clustenmd:
            perform_clustenmd_analysis(centered_ensemble, model, output_dir)

    # Final summary
    print_summary(model, calphas, msf, bfactors, args.method, chain_info, interface_residues)

    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print(f"Results saved to: {output_dir}/")
    print(f"\nKey output files:")
    print(f"  - comprehensive_analysis.png")
    print(f"  - 01_square_fluctuations.png")
    print(f"  - eigenvalues.txt")
    print(f"  - mode_shapes_combined.png")
    if args.ensemble_pca and ensemble is not None:
        print(f"  - pca_comprehensive_analysis.png")
        print(f"  - anm_pca_comparison.png")
    if args.clustenmd and args.ensemble_pca and ensemble is not None:
        print(f"  - clustenmd_analysis.png")
    if interface_residues:
        print(f"  - interface_analysis.txt")
    if args.allosteric:
        print(f"  - allosteric_pathways.txt")
    print(f"\n{'='*60}\n")




if __name__ == '__main__':
    main()

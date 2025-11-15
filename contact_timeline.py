#!/usr/bin/env python3
# version1.0
"""
ProDy Inter-Chain Contact Timeline Analysis Script

This script tracks residue-residue contacts between two chains across
an ensemble of conformers generated from ANM modes.

Usage:
  # Basic usage with default parameters
  python3 inter_chain_contact_timeline.py --pdb structure.pdb

  # Specify chains explicitly
  python3 inter_chain_contact_timeline.py --pdb structure.pdb --chain-a A --chain-b B

  # Custom ANM parameters
  python3 inter_chain_contact_timeline.py --pdb structure.pdb --anm-modes 10 --num-conformers 50

  # Custom contact cutoff
  python3 inter_chain_contact_timeline.py --pdb structure.pdb --contact-cutoff 5.0
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cdist

try:
    from prody import *  # type: ignore
except ImportError:
    print("Error: ProDy is not installed. Install with: pip install prody")
    sys.exit(1)


def setup_output_directory(pdb_name, output_dir=None):
    """Create output directory for results"""
    if output_dir is None:
        output_dir = f"{pdb_name}_Contact_Timeline"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir


def analyze_structure_chains(structure):
    """Analyze available chains in the structure"""
    all_chains = structure.select('calpha')
    chain_ids = np.unique(all_chains.getChids())
    
    chain_info = {}
    for chain_id in chain_ids:
        chain_atoms = structure.select(f'calpha and chain {chain_id}')
        chain_info[chain_id] = {
            'num_residues': chain_atoms.numAtoms(),
            'residue_range': (chain_atoms.getResnums()[0], chain_atoms.getResnums()[-1])
        }
    
    return chain_ids, chain_info


def print_structure_summary(chain_ids, chain_info):
    """Print summary of structure chains"""
    print(f"\n{'='*60}")
    print("STRUCTURE ANALYSIS")
    print(f"{'='*60}")
    print(f"Number of chains: {len(chain_ids)}")
    print(f"\nChain details:")
    for chain_id in sorted(chain_ids):
        info = chain_info[chain_id]
        print(f"  Chain {chain_id}: {info['num_residues']} residues "
              f"(residues {info['residue_range'][0]}-{info['residue_range'][1]})")
    print(f"{'='*60}")


def select_chains(structure, chain_a=None, chain_b=None):
    """Select two chains for analysis"""
    all_chains = structure.select('calpha')
    available_chains = sorted(np.unique(all_chains.getChids()))
    
    if len(available_chains) < 2:
        print(f"Error: Structure must have at least 2 chains. Found: {len(available_chains)}")
        sys.exit(1)
    
    # Auto-select chains if not specified
    if chain_a is None:
        chain_a = available_chains[0]
        print(f"Auto-selected chain A: {chain_a}")
    
    if chain_b is None:
        if len(available_chains) > 1:
            chain_b = available_chains[1]
        else:
            print("Error: Need at least 2 chains for analysis")
            sys.exit(1)
        print(f"Auto-selected chain B: {chain_b}")
    
    # Verify chains exist
    if chain_a not in available_chains:
        print(f"Error: Chain {chain_a} not found. Available chains: {available_chains}")
        sys.exit(1)
    
    if chain_b not in available_chains:
        print(f"Error: Chain {chain_b} not found. Available chains: {available_chains}")
        sys.exit(1)
    
    # Select chain atoms
    atoms_a = structure.select(f'calpha and chain {chain_a}')
    atoms_b = structure.select(f'calpha and chain {chain_b}')
    
    print(f"\nSelected chains for analysis:")
    print(f"  Chain {chain_a}: {atoms_a.numAtoms()} residues")
    print(f"  Chain {chain_b}: {atoms_b.numAtoms()} residues")
    
    return atoms_a, atoms_b, chain_a, chain_b


def build_anm_and_sample(structure, n_modes=5, n_confs=20, amplitude=1.0, cutoff=15.0):
    """Build ANM model and generate conformer ensemble"""
    print(f"\n{'='*60}")
    print("ANM CONSTRUCTION AND CONFORMER SAMPLING")
    print(f"{'='*60}")
    
    # Select all CA atoms for ANM
    all_ca = structure.select('calpha')
    
    print(f"Building ANM with {n_modes} modes (cutoff: {cutoff} Å)...")
    anm = ANM('ANM_Model')
    anm.buildHessian(all_ca, cutoff=cutoff)
    anm.calcModes(n_modes=n_modes)
    
    print(f"ANM calculation complete")
    print(f"  Total modes calculated: {anm.numModes()}")
    print(f"  Modes used for sampling: {n_modes}")
    
    # Sample conformers
    print(f"\nGenerating {n_confs} conformers per mode...")
    print(f"  Mode amplitude: {amplitude} Å")
    
    ensemble = sampleModes(modes=anm[:n_modes], atoms=all_ca, 
                          n_confs=n_confs, rmsd=amplitude)
    
    total_conformers = ensemble.numConfs()
    print(f"Generated ensemble with {total_conformers} conformers")
    print(f"{'='*60}")
    
    return anm, ensemble, all_ca


def calculate_contact_matrix(coords_a, coords_b, cutoff=4.5):
    """
    Calculate contact matrix between two sets of coordinates
    
    Returns:
        Binary matrix where 1 = contact, 0 = no contact
    """
    distances = cdist(coords_a, coords_b)
    contacts = (distances <= cutoff).astype(int)
    return contacts, distances


def build_contact_timeline(ensemble, atoms_a, atoms_b, all_ca_ref, chain_a_id, chain_b_id, contact_cutoff):
    """
    Build contact timeline across all conformers
    
    Args:
        ensemble: ProDy ensemble of conformers
        atoms_a: Atoms object for chain A
        atoms_b: Atoms object for chain B
        all_ca_ref: Reference CA atoms (for chain mapping)
        chain_a_id: Chain ID for chain A
        chain_b_id: Chain ID for chain B
        contact_cutoff: Distance cutoff for contacts
    
    Returns:
        contact_timeline: Binary matrix (num_residue_pairs, num_conformers)
        contact_pairs: List of tuples (res_a_idx, res_b_idx)
        residue_info: Dict with residue information
    """
    print(f"\n{'='*60}")
    print("CONTACT TIMELINE CONSTRUCTION")
    print(f"{'='*60}")
    print(f"Contact cutoff: {contact_cutoff} Å")
    
    n_conformers = ensemble.numConfs()
    n_res_a = atoms_a.numAtoms()
    n_res_b = atoms_b.numAtoms()
    
    print(f"Analyzing {n_conformers} conformers...")
    print(f"  Chain A: {n_res_a} residues")
    print(f"  Chain B: {n_res_b} residues")
    print(f"  Maximum possible contacts: {n_res_a * n_res_b}")
    
    # Get residue information
    resnums_a = atoms_a.getResnums()
    resnames_a = atoms_a.getResnames()
    resnums_b = atoms_b.getResnums()
    resnames_b = atoms_b.getResnames()
    
    # Get chain masks from reference CA atoms
    all_chids = all_ca_ref.getChids()
    chain_a_mask = all_chids == chain_a_id
    chain_b_mask = all_chids == chain_b_id
    
    # Initialize contact storage
    all_contacts = []
    
    # Process each conformer
    for conf_idx in range(n_conformers):
        if (conf_idx + 1) % 10 == 0:
            print(f"  Processing conformer {conf_idx + 1}/{n_conformers}...")
        
        # Get coordinates for this conformer
        coords = ensemble.getCoordsets(conf_idx)
        
        # Extract chain-specific coordinates using masks
        coords_a = coords[chain_a_mask]
        coords_b = coords[chain_b_mask]
        
        # Calculate contacts
        contacts, _ = calculate_contact_matrix(coords_a, coords_b, contact_cutoff)
        all_contacts.append(contacts)
    
    # Stack into 3D array: (conformers, res_a, res_b)
    contact_array = np.array(all_contacts)
    
    # Find all unique contact pairs that occur at least once
    ever_contacted = contact_array.sum(axis=0) > 0
    contact_pairs = np.argwhere(ever_contacted)
    
    print(f"\nContact analysis complete:")
    print(f"  Unique residue pairs that contact: {len(contact_pairs)}")
    
    # Build timeline matrix: (contact_pairs, conformers)
    contact_timeline = np.zeros((len(contact_pairs), n_conformers), dtype=int)
    
    for pair_idx, (res_a_idx, res_b_idx) in enumerate(contact_pairs):
        contact_timeline[pair_idx, :] = contact_array[:, res_a_idx, res_b_idx]
    
    # Calculate contact statistics
    contact_persistence = contact_timeline.sum(axis=1) / n_conformers
    
    print(f"\nContact persistence statistics:")
    print(f"  Always present (100%): {np.sum(contact_persistence == 1.0)} contacts")
    print(f"  Frequently present (>75%): {np.sum(contact_persistence > 0.75)} contacts")
    print(f"  Sometimes present (25-75%): {np.sum((contact_persistence >= 0.25) & (contact_persistence <= 0.75))} contacts")
    print(f"  Rarely present (<25%): {np.sum(contact_persistence < 0.25)} contacts")
    print(f"{'='*60}")
    
    # Store residue information
    residue_info = {
        'chain_a_id': chain_a_id,
        'chain_b_id': chain_b_id,
        'resnums_a': resnums_a,
        'resnames_a': resnames_a,
        'resnums_b': resnums_b,
        'resnames_b': resnames_b,
        'contact_pairs': contact_pairs,
        'contact_persistence': contact_persistence
    }
    
    return contact_timeline, contact_pairs, residue_info


def plot_contact_timeline_heatmap(contact_timeline, residue_info, output_dir):
    """Plot contact timeline as heatmap"""
    print("\nGenerating contact timeline heatmap...")
    
    n_contacts, n_conformers = contact_timeline.shape
    
    if n_contacts == 0:
        print("Warning: No contacts found. Skipping heatmap.")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, max(8, n_contacts * 0.15)))
    
    # Plot heatmap
    im = ax.imshow(contact_timeline, cmap='RdYlGn', aspect='auto', 
                   interpolation='nearest', vmin=0, vmax=1)
    
    # Set labels
    ax.set_xlabel('Conformer Index', fontsize=12, fontweight='bold')
    ax.set_ylabel('Contact Pair', fontsize=12, fontweight='bold')
    
    chain_a_id = residue_info['chain_a_id']
    chain_b_id = residue_info['chain_b_id']
    ax.set_title(f'Inter-Chain Contact Timeline: Chain {chain_a_id} - Chain {chain_b_id}\n'
                f'({n_contacts} contact pairs, {n_conformers} conformers)',
                fontsize=14, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Contact Presence', fontsize=11, fontweight='bold')
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(['No Contact', 'Contact'])
    
    # Add y-axis labels for contact pairs (if not too many)
    if n_contacts <= 50:
        contact_pairs = residue_info['contact_pairs']
        resnums_a = residue_info['resnums_a']
        resnames_a = residue_info['resnames_a']
        resnums_b = residue_info['resnums_b']
        resnames_b = residue_info['resnames_b']
        
        ytick_labels = []
        for res_a_idx, res_b_idx in contact_pairs:
            label = f"{resnames_a[res_a_idx]}{resnums_a[res_a_idx]}-{resnames_b[res_b_idx]}{resnums_b[res_b_idx]}"
            ytick_labels.append(label)
        
        ax.set_yticks(range(n_contacts))
        ax.set_yticklabels(ytick_labels, fontsize=8)
    else:
        ax.set_ylabel(f'Contact Pair Index (1-{n_contacts})', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/contact_timeline_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/contact_timeline_heatmap.png")


def plot_contact_persistence(residue_info, output_dir):
    """Plot contact persistence distribution"""
    print("\nGenerating contact persistence plot...")
    
    contact_persistence = residue_info['contact_persistence']
    contact_pairs = residue_info['contact_pairs']
    
    if len(contact_persistence) == 0:
        print("Warning: No contacts found. Skipping persistence plot.")
        return
    
    # Sort by persistence
    sorted_indices = np.argsort(contact_persistence)[::-1]
    sorted_persistence = contact_persistence[sorted_indices]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: Persistence histogram
    ax1.hist(contact_persistence, bins=20, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.axvline(x=0.5, color='red', linestyle='--', linewidth=2, label='50% threshold')
    ax1.set_xlabel('Contact Persistence (fraction of conformers)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Number of Contact Pairs', fontsize=11, fontweight='bold')
    ax1.set_title('Distribution of Contact Persistence', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Plot 2: Sorted persistence
    x_indices = np.arange(len(sorted_persistence))
    colors = plt.cm.RdYlGn(sorted_persistence)
    
    ax2.bar(x_indices, sorted_persistence, color=colors, edgecolor='black', linewidth=0.5)
    ax2.axhline(y=0.75, color='green', linestyle='--', linewidth=1.5, alpha=0.7, label='75% threshold')
    ax2.axhline(y=0.5, color='orange', linestyle='--', linewidth=1.5, alpha=0.7, label='50% threshold')
    ax2.axhline(y=0.25, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='25% threshold')
    
    ax2.set_xlabel('Contact Pair (sorted by persistence)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Contact Persistence', fontsize=11, fontweight='bold')
    ax2.set_title('Contact Pairs Sorted by Persistence', fontsize=12, fontweight='bold')
    ax2.set_ylim([0, 1.05])
    ax2.legend()
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/contact_persistence.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/contact_persistence.png")


def plot_contact_network(residue_info, contact_timeline, output_dir, persistence_threshold=0.5):
    """Plot contact network showing stable vs. dynamic contacts"""
    print("\nGenerating contact network visualization...")
    
    contact_persistence = residue_info['contact_persistence']
    contact_pairs = residue_info['contact_pairs']
    resnums_a = residue_info['resnums_a']
    resnums_b = residue_info['resnums_b']
    
    if len(contact_persistence) == 0:
        print("Warning: No contacts found. Skipping network plot.")
        return
    
    # Categorize contacts
    stable_contacts = contact_persistence >= persistence_threshold
    dynamic_contacts = contact_persistence < persistence_threshold
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    # Plot for each chain showing contact residues
    for ax_idx, (ax, chain_res) in enumerate([(axes[0], resnums_a), (axes[1], resnums_b)]):
        chain_id = residue_info['chain_a_id'] if ax_idx == 0 else residue_info['chain_b_id']
        
        # Count contacts per residue
        contact_counts = np.zeros(len(chain_res))
        stable_counts = np.zeros(len(chain_res))
        
        for pair_idx, (res_a_idx, res_b_idx) in enumerate(contact_pairs):
            res_idx = res_a_idx if ax_idx == 0 else res_b_idx
            contact_counts[res_idx] += 1
            if stable_contacts[pair_idx]:
                stable_counts[res_idx] += 1
        
        # Plot
        x = np.arange(len(chain_res))
        ax.bar(x, contact_counts, color='lightblue', edgecolor='black', 
               linewidth=0.5, label='All contacts', alpha=0.7)
        ax.bar(x, stable_counts, color='darkgreen', edgecolor='black', 
               linewidth=0.5, label=f'Stable (≥{persistence_threshold*100:.0f}%)', alpha=0.9)
        
        ax.set_xlabel(f'Chain {chain_id} Residue Index', fontsize=11, fontweight='bold')
        ax.set_ylabel('Number of Inter-chain Contacts', fontsize=11, fontweight='bold')
        ax.set_title(f'Chain {chain_id} Contact Profile', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3, axis='y')
        
        # Add residue numbers on x-axis (sample if too many)
        if len(chain_res) <= 50:
            ax.set_xticks(x)
            ax.set_xticklabels(chain_res, rotation=90, fontsize=8)
        else:
            step = len(chain_res) // 20
            ax.set_xticks(x[::step])
            ax.set_xticklabels(chain_res[::step], rotation=90, fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/contact_network.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/contact_network.png")


def save_contact_timeline_matrix(contact_timeline, residue_info, output_dir):
    """Save contact timeline as CSV and text files"""
    print("\nSaving contact timeline matrix...")
    
    contact_pairs = residue_info['contact_pairs']
    resnums_a = residue_info['resnums_a']
    resnames_a = residue_info['resnames_a']
    resnums_b = residue_info['resnums_b']
    resnames_b = residue_info['resnames_b']
    contact_persistence = residue_info['contact_persistence']
    chain_a_id = residue_info['chain_a_id']
    chain_b_id = residue_info['chain_b_id']
    
    # Save as CSV
    csv_path = f'{output_dir}/contact_timeline.csv'
    with open(csv_path, 'w') as f:
        # Header
        f.write(f"# Inter-Chain Contact Timeline: Chain {chain_a_id} - Chain {chain_b_id}\n")
        f.write(f"# Rows: Contact pairs (residue pairs)\n")
        f.write(f"# Columns: Conformer index\n")
        f.write(f"# Values: 1 = contact present, 0 = no contact\n#\n")
        
        # Column headers
        f.write("ResA,ResB,Persistence," + ",".join([f"Conf_{i}" for i in range(contact_timeline.shape[1])]) + "\n")
        
        # Data rows
        for pair_idx, (res_a_idx, res_b_idx) in enumerate(contact_pairs):
            res_a_label = f"{resnames_a[res_a_idx]}{resnums_a[res_a_idx]}"
            res_b_label = f"{resnames_b[res_b_idx]}{resnums_b[res_b_idx]}"
            persistence = contact_persistence[pair_idx]
            
            row_data = ",".join(map(str, contact_timeline[pair_idx, :]))
            f.write(f"{res_a_label},{res_b_label},{persistence:.4f},{row_data}\n")
    
    print(f"Saved: {csv_path}")


def save_detailed_report(contact_timeline, residue_info, anm, output_dir, args):
    """Save detailed text report"""
    print("\nSaving detailed report...")
    
    contact_pairs = residue_info['contact_pairs']
    resnums_a = residue_info['resnums_a']
    resnames_a = residue_info['resnames_a']
    resnums_b = residue_info['resnums_b']
    resnames_b = residue_info['resnames_b']
    contact_persistence = residue_info['contact_persistence']
    chain_a_id = residue_info['chain_a_id']
    chain_b_id = residue_info['chain_b_id']
    
    report_path = f'{output_dir}/contact_analysis_report.txt'
    
    with open(report_path, 'w') as f:
        f.write("="*70 + "\n")
        f.write("INTER-CHAIN CONTACT TIMELINE ANALYSIS REPORT\n")
        f.write("="*70 + "\n\n")
        
        # Analysis parameters
        f.write("ANALYSIS PARAMETERS\n")
        f.write("-"*70 + "\n")
        f.write(f"PDB file: {args.pdb}\n")
        f.write(f"Chain A: {chain_a_id} ({len(resnums_a)} residues)\n")
        f.write(f"Chain B: {chain_b_id} ({len(resnums_b)} residues)\n")
        f.write(f"Contact cutoff: {args.contact_cutoff} Å\n")
        f.write(f"ANM modes: {args.anm_modes}\n")
        f.write(f"Conformers per mode: {args.num_conformers}\n")
        f.write(f"Total conformers: {contact_timeline.shape[1]}\n")
        f.write(f"Mode amplitude: {args.mode_amplitude} Å\n")
        f.write(f"ANM cutoff: {args.anm_cutoff} Å\n\n")
        
        # Summary statistics
        f.write("CONTACT SUMMARY\n")
        f.write("-"*70 + "\n")
        f.write(f"Total unique contact pairs identified: {len(contact_pairs)}\n")
        f.write(f"Stable contacts (>75% persistence): {np.sum(contact_persistence > 0.75)}\n")
        f.write(f"Frequent contacts (50-75% persistence): {np.sum((contact_persistence >= 0.5) & (contact_persistence <= 0.75))}\n")
        f.write(f"Variable contacts (25-50% persistence): {np.sum((contact_persistence >= 0.25) & (contact_persistence < 0.5))}\n")
        f.write(f"Rare contacts (<25% persistence): {np.sum(contact_persistence < 0.25)}\n\n")
        
        # Contact pairs sorted by persistence
        f.write("CONTACT PAIRS SORTED BY PERSISTENCE\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Rank':<6} {'Chain A Residue':<18} {'Chain B Residue':<18} {'Persistence':<12} {'Category'}\n")
        f.write("-"*70 + "\n")
        
        sorted_indices = np.argsort(contact_persistence)[::-1]
        
        for rank, idx in enumerate(sorted_indices, 1):
            res_a_idx, res_b_idx = contact_pairs[idx]
            res_a_label = f"{resnames_a[res_a_idx]}{resnums_a[res_a_idx]}"
            res_b_label = f"{resnames_b[res_b_idx]}{resnums_b[res_b_idx]}"
            persistence = contact_persistence[idx]
            
            if persistence > 0.75:
                category = "Stable"
            elif persistence >= 0.5:
                category = "Frequent"
            elif persistence >= 0.25:
                category = "Variable"
            else:
                category = "Rare"
            
            f.write(f"{rank:<6} {res_a_label:<18} {res_b_label:<18} {persistence:<12.4f} {category}\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*70 + "\n")
    
    print(f"Saved: {report_path}")


def save_pymol_script(residue_info, contact_timeline, output_dir, persistence_threshold=0.75):
    """Generate PyMOL script to visualize stable contacts"""
    print("\nGenerating PyMOL visualization script...")
    
    contact_pairs = residue_info['contact_pairs']
    contact_persistence = residue_info['contact_persistence']
    resnums_a = residue_info['resnums_a']
    resnums_b = residue_info['resnums_b']
    chain_a_id = residue_info['chain_a_id']
    chain_b_id = residue_info['chain_b_id']
    
    # Identify stable contacts
    stable_indices = np.where(contact_persistence >= persistence_threshold)[0]
    
    script_path = f'{output_dir}/visualize_contacts.pml'
    
    with open(script_path, 'w') as f:
        f.write("# PyMOL script to visualize inter-chain contacts\n\n")
        
        f.write("# Color chains\n")
        f.write(f"color blue, chain {chain_a_id}\n")
        f.write(f"color red, chain {chain_b_id}\n\n")
        
        f.write("# Show cartoon representation\n")
        f.write("hide everything\n")
        f.write("show cartoon\n")
        f.write("set cartoon_fancy_helices, 1\n\n")
        
        f.write(f"# Highlight stable contact residues (persistence ≥ {persistence_threshold*100:.0f}%)\n")
        
        # Collect all stable contact residues
        stable_res_a = set()
        stable_res_b = set()
        
        for idx in stable_indices:
            res_a_idx, res_b_idx = contact_pairs[idx]
            stable_res_a.add(resnums_a[res_a_idx])
            stable_res_b.add(resnums_b[res_b_idx])
        
        if stable_res_a:
            res_list_a = '+'.join(map(str, sorted(stable_res_a)))
            f.write(f"select stable_contacts_A, chain {chain_a_id} and resi {res_list_a}\n")
            f.write("show sticks, stable_contacts_A\n")
            f.write("color cyan, stable_contacts_A\n\n")
        
        if stable_res_b:
            res_list_b = '+'.join(map(str, sorted(stable_res_b)))
            f.write(f"select stable_contacts_B, chain {chain_b_id} and resi {res_list_b}\n")
            f.write("show sticks, stable_contacts_B\n")
            f.write("color magenta, stable_contacts_B\n\n")
        
        f.write("# Draw distance lines for stable contacts\n")
        for idx in stable_indices[:20]:  # Limit to first 20 for clarity
            res_a_idx, res_b_idx = contact_pairs[idx]
            resnum_a = resnums_a[res_a_idx]
            resnum_b = resnums_b[res_b_idx]
            f.write(f"distance stable_contact_{idx}, chain {chain_a_id} and resi {resnum_a} and name CA, "
                   f"chain {chain_b_id} and resi {resnum_b} and name CA\n")
        
        f.write("\nhide labels\n")
        f.write("color yellow, stable_contact_*\n")
        f.write("set dash_width, 2\n\n")
        
        f.write("# Display settings\n")
        f.write("bg_color white\n")
        f.write("set ray_shadows, 0\n")
        f.write("set antialias, 2\n")
        f.write("zoom\n\n")
        
        f.write(f"# Summary: {len(stable_indices)} stable contacts visualized\n")
    
    print(f"Saved: {script_path}")


def main():
    parser = argparse.ArgumentParser(
        description='ProDy Inter-Chain Contact Timeline Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (auto-detect chains)
  python3 inter_chain_contact_timeline.py --pdb structure.pdb

  # Specify chains explicitly
  python3 inter_chain_contact_timeline.py --pdb structure.pdb --chain-a A --chain-b B

  # Custom ANM parameters
  python3 inter_chain_contact_timeline.py --pdb structure.pdb --anm-modes 10 --num-conformers 50

  # Custom contact cutoff (looser contacts)
  python3 inter_chain_contact_timeline.py --pdb structure.pdb --contact-cutoff 5.5

  # Full customization
  python3 inter_chain_contact_timeline.py --pdb structure.pdb --chain-a H --chain-b L \\
      --anm-modes 8 --num-conformers 30 --contact-cutoff 4.0 --mode-amplitude 1.5
        """)
    
    parser.add_argument('--pdb', required=True, 
                       help='Input PDB file (must contain at least 2 chains)')
    parser.add_argument('--chain-a', default=None, 
                       help='Chain ID for first chain (default: auto-select first chain)')
    parser.add_argument('--chain-b', default=None, 
                       help='Chain ID for second chain (default: auto-select second chain)')
    parser.add_argument('--anm-modes', type=int, default=5, 
                       help='Number of ANM modes to use (default: 5)')
    parser.add_argument('--num-conformers', type=int, default=20, 
                       help='Number of conformers per mode (default: 20)')
    parser.add_argument('--contact-cutoff', type=float, default=4.5, 
                       help='Distance cutoff for contacts in Å (default: 4.5)')
    parser.add_argument('--mode-amplitude', type=float, default=1.0, 
                       help='Amplitude for mode deformation in Å (default: 1.0)')
    parser.add_argument('--anm-cutoff', type=float, default=15.0, 
                       help='ANM cutoff distance in Å (default: 15.0)')
    parser.add_argument('--output', default=None, 
                       help='Output directory (default: <pdb_name>_Contact_Timeline)')
    parser.add_argument('--persistence-threshold', type=float, default=0.75,
                       help='Persistence threshold for stable contacts (default: 0.75)')
    
    args = parser.parse_args()
    
    # Validate parameters
    if args.anm_modes < 1:
        print("Error: --anm-modes must be at least 1")
        sys.exit(1)
    
    if args.num_conformers < 1:
        print("Error: --num-conformers must be at least 1")
        sys.exit(1)
    
    if args.contact_cutoff <= 0:
        print("Error: --contact-cutoff must be positive")
        sys.exit(1)
    
    if not (0 < args.persistence_threshold < 1):
        print("Error: --persistence-threshold must be between 0 and 1")
        sys.exit(1)
    
    # Setup
    pdb_name = os.path.splitext(os.path.basename(args.pdb))[0]
    output_dir = setup_output_directory(pdb_name, args.output)
    
    print("="*70)
    print("INTER-CHAIN CONTACT TIMELINE ANALYSIS")
    print("="*70)
    print(f"PDB file: {args.pdb}")
    print(f"Output directory: {output_dir}")
    
    # Load structure
    print(f"\nLoading structure...")
    structure = parsePDB(args.pdb)
    
    if structure is None:
        print(f"Error: Could not load PDB file: {args.pdb}")
        sys.exit(1)
    
    # Analyze chains
    chain_ids, chain_info = analyze_structure_chains(structure)
    print_structure_summary(chain_ids, chain_info)
    
    # Select chains for analysis
    atoms_a, atoms_b, chain_a_id, chain_b_id = select_chains(
        structure, args.chain_a, args.chain_b
    )
    
    # Build ANM and generate ensemble
    anm, ensemble, all_ca = build_anm_and_sample(
        structure, 
        n_modes=args.anm_modes,
        n_confs=args.num_conformers,
        amplitude=args.mode_amplitude,
        cutoff=args.anm_cutoff
    )
    
    # Build contact timeline
    contact_timeline, contact_pairs, residue_info = build_contact_timeline(
        ensemble, atoms_a, atoms_b, all_ca, chain_a_id, chain_b_id, args.contact_cutoff
    )
    
    # Generate visualizations
    print(f"\n{'='*70}")
    print("GENERATING VISUALIZATIONS")
    print(f"{'='*70}")
    
    plot_contact_timeline_heatmap(contact_timeline, residue_info, output_dir)
    plot_contact_persistence(residue_info, output_dir)
    plot_contact_network(residue_info, contact_timeline, output_dir, 
                        persistence_threshold=args.persistence_threshold)
    
    # Save results
    print(f"\n{'='*70}")
    print("SAVING RESULTS")
    print(f"{'='*70}")
    
    save_contact_timeline_matrix(contact_timeline, residue_info, output_dir)
    save_detailed_report(contact_timeline, residue_info, anm, output_dir, args)
    save_pymol_script(residue_info, contact_timeline, output_dir, 
                     persistence_threshold=args.persistence_threshold)
    
    # Final summary
    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*70}")
    print(f"\nResults saved to: {output_dir}/")
    print(f"\nGenerated files:")
    print(f"  - contact_timeline_heatmap.png : Heatmap of all contacts across conformers")
    print(f"  - contact_persistence.png : Distribution and ranking of contact persistence")
    print(f"  - contact_network.png : Per-chain contact profiles")
    print(f"  - contact_timeline.csv : Binary contact matrix (importable to Excel/Python)")
    print(f"  - contact_analysis_report.txt : Detailed text report")
    print(f"  - visualize_contacts.pml : PyMOL script for 3D visualization")
    
    print(f"\nTo visualize in PyMOL:")
    print(f"  pymol {args.pdb} {output_dir}/visualize_contacts.pml")
    
    print(f"\nContact summary:")
    print(f"  Total contact pairs: {len(contact_pairs)}")
    print(f"  Stable contacts (≥{args.persistence_threshold*100:.0f}%): "
          f"{np.sum(residue_info['contact_persistence'] >= args.persistence_threshold)}")
    print(f"  Conformers analyzed: {contact_timeline.shape[1]}")
    
    print(f"\n{'='*70}")


if __name__ == '__main__':
    main()

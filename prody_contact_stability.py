#!/usr/bin/env python3
# version1.0
"""
ProDy Contact, SASA, and Stability Analysis Script

This script performs comprehensive structural analysis including:
1. Intermolecular contact mapping (between chains)
2. Intramolecular contact mapping (within chains)
3. Hydrophobicity and solvent accessibility (SASA)
4. Noncovalent interaction mapping (H-bonds, salt bridges)

Usage:
  # Analyze all chains
  python3 prody_contact_stability_analysis.py --pdb complex.pdb

  # Analyze interface between specific chains
  python3 prody_contact_stability_analysis.py --pdb complex.pdb --chain-a A --chain-b B

  # Full analysis for single chain
  python3 prody_contact_stability_analysis.py --pdb protein.pdb --chain A
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from itertools import combinations

try:
    from prody import *  # type: ignore
    # Try to import SASA calculation
    try:
        from prody.proteins import calcSASA  # type: ignore
        SASA_AVAILABLE = True
    except ImportError:
        SASA_AVAILABLE = False
        print("Warning: calcSASA not available in this ProDy version. SASA analysis will use approximation.")
except ImportError:
    print("Error: ProDy is not installed. Install with: pip install prody")
    sys.exit(1)

# Hydrophobicity scales (Kyte-Doolittle)
HYDROPHOBICITY = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
    'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
    'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}

# Aromatic residues for pi interactions
AROMATIC = {'PHE', 'TYR', 'TRP', 'HIS'}

# Charged residues
POSITIVE = {'ARG', 'LYS', 'HIS'}
NEGATIVE = {'ASP', 'GLU'}


def setup_output_directory(pdb_name):
    """Create output directory for results"""
    output_dir = f"{pdb_name}_Contact_Stability_Analysis"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir


def analyze_chain_structure(structure):
    """Analyze and report chain structure"""
    chain_info = {}
    
    for chain in structure.iterChains():
        chain_id = chain.getChid()
        ca_atoms = chain.select('calpha')
        if ca_atoms:
            chain_info[chain_id] = {
                'n_residues': ca_atoms.numAtoms(),
                'start_resnum': ca_atoms.getResnums()[0],
                'end_resnum': ca_atoms.getResnums()[-1],
                'residues': ca_atoms
            }
    
    return chain_info


def print_chain_summary(chain_info):
    """Print summary of chain structure"""
    print(f"\n{'='*60}")
    print("CHAIN STRUCTURE ANALYSIS")
    print(f"{'='*60}")
    print(f"Number of chains: {len(chain_info)}")
    print(f"\nChain details:")
    for chain_id, info in sorted(chain_info.items()):
        print(f"  Chain {chain_id}: {info['n_residues']} residues "
              f"(PDB: {info['start_resnum']}-{info['end_resnum']})")
    print(f"{'='*60}")


# ============================================================================
# 1. INTERMOLECULAR CONTACT MAPPING
# ============================================================================

def calculate_intermolecular_contacts(structure, chain_a, chain_b, cutoff=5.0):
    """
    Calculate contacts between two chains
    
    Parameters:
    - structure: ProDy structure
    - chain_a, chain_b: Chain IDs
    - cutoff: Distance cutoff in Angstroms (default: 5.0)
    """
    print(f"\n{'='*60}")
    print(f"INTERMOLECULAR CONTACTS: Chain {chain_a} ↔ Chain {chain_b}")
    print(f"{'='*60}")
    
    # Select chains
    atoms_a = structure.select(f'chain {chain_a}')
    atoms_b = structure.select(f'chain {chain_b}')
    
    if atoms_a is None or atoms_b is None:
        print(f"Error: Could not find chain {chain_a} or {chain_b}")
        return None
    
    # Get CA atoms for residue-level analysis
    ca_a = atoms_a.select('calpha')
    ca_b = atoms_b.select('calpha')
    
    # Get all heavy atoms for precise distance calculation
    heavy_a = atoms_a.select('not hydrogen')
    heavy_b = atoms_b.select('not hydrogen')
    
    contacts = []
    contact_counts_a = defaultdict(int)
    contact_counts_b = defaultdict(int)
    
    print(f"Analyzing contacts (cutoff: {cutoff} Å)...")
    
    # For each residue in chain A
    for i, ca_res_a in enumerate(ca_a):
        resnum_a = ca_res_a.getResnum()
        resname_a = ca_res_a.getResname()
        
        # Get all heavy atoms for this residue
        res_atoms_a = heavy_a.select(f'resnum {resnum_a}')
        if res_atoms_a is None:
            continue
        
        # Check against all residues in chain B
        for j, ca_res_b in enumerate(ca_b):
            resnum_b = ca_res_b.getResnum()
            resname_b = ca_res_b.getResname()
            
            # Quick check: if CA distance > cutoff + 10, skip
            ca_dist = np.linalg.norm(ca_res_a.getCoords() - ca_res_b.getCoords())
            if ca_dist > cutoff + 10:
                continue
            
            # Get all heavy atoms for this residue
            res_atoms_b = heavy_b.select(f'resnum {resnum_b}')
            if res_atoms_b is None:
                continue
            
            # Calculate minimum distance between any heavy atoms
            coords_a = res_atoms_a.getCoords()
            coords_b = res_atoms_b.getCoords()
            
            # Calculate all pairwise distances
            distances = np.sqrt(np.sum((coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :])**2, axis=2))
            min_dist = np.min(distances)
            
            if min_dist <= cutoff:
                contacts.append({
                    'chain_a': chain_a,
                    'resnum_a': resnum_a,
                    'resname_a': resname_a,
                    'chain_b': chain_b,
                    'resnum_b': resnum_b,
                    'resname_b': resname_b,
                    'distance': min_dist
                })
                contact_counts_a[resnum_a] += 1
                contact_counts_b[resnum_b] += 1
    
    print(f"Found {len(contacts)} residue-residue contacts")
    print(f"Interface residues: {len(contact_counts_a)} (chain {chain_a}), "
          f"{len(contact_counts_b)} (chain {chain_b})")
    
    return {
        'contacts': contacts,
        'counts_a': contact_counts_a,
        'counts_b': contact_counts_b,
        'ca_a': ca_a,
        'ca_b': ca_b
    }


# ============================================================================
# 2. INTRAMOLECULAR CONTACT MAPPING
# ============================================================================

def calculate_intramolecular_contacts(chain_atoms, chain_id, contact_cutoff=5.0, neighbor_cutoff=8.0):
    """
    Calculate contacts and neighbors within a single chain
    
    Parameters:
    - chain_atoms: ProDy atoms for the chain
    - chain_id: Chain identifier
    - contact_cutoff: Distance for contact (Å, default: 5.0)
    - neighbor_cutoff: Distance for CA neighbors (Å, default: 8.0)
    """
    print(f"\n{'='*60}")
    print(f"INTRAMOLECULAR CONTACTS: Chain {chain_id}")
    print(f"{'='*60}")
    
    ca_atoms = chain_atoms.select('calpha')
    heavy_atoms = chain_atoms.select('not hydrogen')
    
    n_residues = ca_atoms.numAtoms()
    resnums = ca_atoms.getResnums()
    resnames = ca_atoms.getResnames()
    
    # Contact matrix
    contact_matrix = np.zeros((n_residues, n_residues))
    neighbor_matrix = np.zeros((n_residues, n_residues))
    contact_list = []
    
    print(f"Calculating contacts for {n_residues} residues...")
    
    for i in range(n_residues):
        resnum_i = resnums[i]
        resname_i = resnames[i]
        
        # Get heavy atoms for residue i
        res_atoms_i = heavy_atoms.select(f'resnum {resnum_i}')
        if res_atoms_i is None:
            continue
        
        for j in range(i + 1, n_residues):
            resnum_j = resnums[j]
            resname_j = resnames[j]
            
            # Skip sequential neighbors (i, i+1, i+2)
            if abs(j - i) <= 2:
                continue
            
            # CA distance for neighbor check
            ca_i_coords = ca_atoms[i].getCoords()
            ca_j_coords = ca_atoms[j].getCoords()
            ca_dist = np.linalg.norm(ca_i_coords - ca_j_coords)
            
            if ca_dist <= neighbor_cutoff:
                neighbor_matrix[i, j] = 1
                neighbor_matrix[j, i] = 1
                
                # Get heavy atoms for residue j
                res_atoms_j = heavy_atoms.select(f'resnum {resnum_j}')
                if res_atoms_j is None:
                    continue
                
                # Check for contact - calculate pairwise distances
                coords_i = res_atoms_i.getCoords()
                coords_j = res_atoms_j.getCoords()
                
                # Calculate all pairwise distances
                distances = np.sqrt(np.sum((coords_i[:, np.newaxis, :] - coords_j[np.newaxis, :, :])**2, axis=2))
                min_dist = np.min(distances)
                
                if min_dist <= contact_cutoff:
                    contact_matrix[i, j] = 1
                    contact_matrix[j, i] = 1
                    contact_list.append({
                        'res_i': f"{resname_i}{resnum_i}",
                        'res_j': f"{resname_j}{resnum_j}",
                        'distance': min_dist
                    })
    
    contact_counts = np.sum(contact_matrix, axis=1)
    neighbor_counts = np.sum(neighbor_matrix, axis=1)
    
    print(f"Total contacts: {len(contact_list)}")
    print(f"Average contacts per residue: {np.mean(contact_counts):.2f}")
    print(f"Average neighbors per residue: {np.mean(neighbor_counts):.2f}")
    
    return {
        'contact_matrix': contact_matrix,
        'neighbor_matrix': neighbor_matrix,
        'contact_counts': contact_counts,
        'neighbor_counts': neighbor_counts,
        'contact_list': contact_list,
        'resnums': resnums,
        'resnames': resnames
    }


# ============================================================================
# 3. SASA AND HYDROPHOBICITY ANALYSIS
# ============================================================================

def calculate_sasa_hydrophobicity(structure, chain_id=None):
    """
    Calculate solvent accessible surface area and hydrophobicity
    
    Parameters:
    - structure: ProDy structure
    - chain_id: Chain to analyze (None for all)
    """
    global SASA_AVAILABLE
    
    print(f"\n{'='*60}")
    if chain_id:
        print(f"SASA & HYDROPHOBICITY: Chain {chain_id}")
        selection = structure.select(f'chain {chain_id}')
    else:
        print("SASA & HYDROPHOBICITY: All chains")
        selection = structure
    print(f"{'='*60}")
    
    ca_atoms = selection.select('calpha')
    resnums = ca_atoms.getResnums()
    resnames = ca_atoms.getResnames()
    
    # Calculate SASA
    print("Calculating SASA...")
    
    use_sasa = SASA_AVAILABLE
    if use_sasa:
        try:
            sasa_array = calcSASA(selection)
        except Exception as e:
            print(f"Warning: calcSASA failed ({e}), using approximation method")
            use_sasa = False
    
    if not use_sasa:
        # Use contact-based approximation
        print("Using contact-based SASA approximation...")
        all_atoms = selection.select('not hydrogen')
        residue_sasa = []
        
        for i, (resnum, resname) in enumerate(zip(resnums, resnames)):
            res_atoms = all_atoms.select(f'resnum {resnum}')
            if res_atoms is None:
                residue_sasa.append(0)
                continue
            
            # Count exposed atoms (those with few neighbors)
            res_coords = res_atoms.getCoords()
            n_atoms = len(res_coords)
            
            # Count neighbors within 5Å for each atom
            exposed_count = 0
            for atom_coord in res_coords:
                # Calculate distances to all other atoms
                all_coords = all_atoms.getCoords()
                distances = np.sqrt(np.sum((all_coords - atom_coord)**2, axis=1))
                n_neighbors = np.sum((distances > 0.1) & (distances < 5.0))
                
                # If fewer than 15 neighbors, consider exposed
                if n_neighbors < 15:
                    exposed_count += 1
            
            # Approximate SASA based on exposed atoms
            approx_sasa = (exposed_count / n_atoms) * 100 if n_atoms > 0 else 0
            residue_sasa.append(approx_sasa)
        
        residue_sasa = np.array(residue_sasa)
    else:
        # Get per-residue SASA from calcSASA
        residue_sasa = []
        
        for i, (resnum, resname) in enumerate(zip(resnums, resnames)):
            # Get SASA for this residue
            res_atoms = selection.select(f'resnum {resnum}')
            if res_atoms is None:
                residue_sasa.append(0)
                continue
            
            # Sum SASA for all atoms in residue
            res_indices = res_atoms.getIndices()
            local_indices = [np.where(selection.getIndices() == idx)[0][0] for idx in res_indices]
            res_sasa = np.sum(sasa_array[local_indices])
            residue_sasa.append(res_sasa)
        
        residue_sasa = np.array(residue_sasa)
    
    # Get hydrophobicity scores
    residue_hydro = []
    for resname in resnames:
        hydro = HYDROPHOBICITY.get(resname, 0)
        residue_hydro.append(hydro)
    
    residue_hydro = np.array(residue_hydro)
    
    # Classify exposure
    exposure_class = np.zeros(len(residue_sasa), dtype=int)
    exposure_class[residue_sasa > 40] = 2  # Surface
    exposure_class[(residue_sasa > 10) & (residue_sasa <= 40)] = 1  # Partially exposed
    # 0 = Core (buried)
    
    print(f"Core residues (SASA < 10): {np.sum(exposure_class == 0)}")
    print(f"Partially exposed (10-40): {np.sum(exposure_class == 1)}")
    print(f"Surface residues (> 40): {np.sum(exposure_class == 2)}")
    
    return {
        'sasa': residue_sasa,
        'hydrophobicity': residue_hydro,
        'exposure_class': exposure_class,
        'resnums': resnums,
        'resnames': resnames,
        'ca_atoms': ca_atoms
    }


# ============================================================================
# 4. NONCOVALENT INTERACTION MAPPING
# ============================================================================

def find_hydrogen_bonds(structure, chain_a=None, chain_b=None, dist_cutoff=3.5, angle_cutoff=120):
    """
    Identify hydrogen bonds
    
    Parameters:
    - dist_cutoff: H-bond distance cutoff (Å)
    - angle_cutoff: D-H...A angle cutoff (degrees)
    """
    print(f"\n{'='*60}")
    print("HYDROGEN BOND ANALYSIS")
    print(f"{'='*60}")
    
    if chain_a and chain_b:
        donors_sel = structure.select(f'chain {chain_a} and (name N or name O)')
        acceptors_sel = structure.select(f'chain {chain_b} and (name O or name N)')
        print(f"Analyzing inter-chain H-bonds: {chain_a} → {chain_b}")
    else:
        donors_sel = structure.select('name N or name O')
        acceptors_sel = structure.select('name O or name N')
        print("Analyzing all H-bonds")
    
    hbonds = []
    
    # Simple geometric criterion
    for donor in donors_sel:
        donor_resnum = donor.getResnum()
        donor_resname = donor.getResname()
        donor_chain = donor.getChid()
        
        for acceptor in acceptors_sel:
            acceptor_resnum = acceptor.getResnum()
            acceptor_resname = acceptor.getResname()
            acceptor_chain = acceptor.getChid()
            
            # Skip same residue
            if donor_resnum == acceptor_resnum and donor_chain == acceptor_chain:
                continue
            
            dist = np.linalg.norm(donor.getCoords() - acceptor.getCoords())
            
            if dist <= dist_cutoff:
                hbonds.append({
                    'donor': f"{donor_resname}{donor_resnum}",
                    'donor_chain': donor_chain,
                    'acceptor': f"{acceptor_resname}{acceptor_resnum}",
                    'acceptor_chain': acceptor_chain,
                    'distance': dist
                })
    
    print(f"Found {len(hbonds)} potential hydrogen bonds")
    return hbonds


def find_salt_bridges(structure, chain_a=None, chain_b=None, cutoff=4.0):
    """
    Identify salt bridges between charged residues
    """
    print(f"\n{'='*60}")
    print("SALT BRIDGE ANALYSIS")
    print(f"{'='*60}")
    
    if chain_a and chain_b:
        positive_sel = structure.select(f'chain {chain_a} and calpha and (resname ARG LYS HIS)')
        negative_sel = structure.select(f'chain {chain_b} and calpha and (resname ASP GLU)')
        print(f"Analyzing inter-chain salt bridges: {chain_a} ↔ {chain_b}")
    else:
        positive_sel = structure.select('calpha and (resname ARG LYS HIS)')
        negative_sel = structure.select('calpha and (resname ASP GLU)')
        print("Analyzing all salt bridges")
    
    salt_bridges = []
    
    if positive_sel is None or negative_sel is None:
        print("No charged residues found")
        return salt_bridges
    
    for pos_res in positive_sel:
        pos_resnum = pos_res.getResnum()
        pos_resname = pos_res.getResname()
        pos_chain = pos_res.getChid()
        
        for neg_res in negative_sel:
            neg_resnum = neg_res.getResnum()
            neg_resname = neg_res.getResname()
            neg_chain = neg_res.getChid()
            
            # Skip same residue
            if pos_resnum == neg_resnum and pos_chain == neg_chain:
                continue
            
            dist = np.linalg.norm(pos_res.getCoords() - neg_res.getCoords())
            
            if dist <= cutoff:
                salt_bridges.append({
                    'positive': f"{pos_resname}{pos_resnum}",
                    'pos_chain': pos_chain,
                    'negative': f"{neg_resname}{neg_resnum}",
                    'neg_chain': neg_chain,
                    'distance': dist
                })
    
    print(f"Found {len(salt_bridges)} potential salt bridges")
    return salt_bridges


def find_aromatic_interactions(structure, chain_a=None, chain_b=None, cutoff=7.0):
    """
    Identify pi-pi and cation-pi interactions
    """
    print(f"\n{'='*60}")
    print("AROMATIC INTERACTION ANALYSIS")
    print(f"{'='*60}")
    
    if chain_a and chain_b:
        aromatic_a = structure.select(f'chain {chain_a} and calpha and (resname PHE TYR TRP HIS)')
        aromatic_b = structure.select(f'chain {chain_b} and calpha and (resname PHE TYR TRP HIS)')
        cation_a = structure.select(f'chain {chain_a} and calpha and (resname ARG LYS)')
        cation_b = structure.select(f'chain {chain_b} and calpha and (resname ARG LYS)')
        print(f"Analyzing aromatic interactions: {chain_a} ↔ {chain_b}")
    else:
        aromatic_a = structure.select('calpha and (resname PHE TYR TRP HIS)')
        aromatic_b = aromatic_a
        cation_a = structure.select('calpha and (resname ARG LYS)')
        cation_b = cation_a
        print("Analyzing all aromatic interactions")
    
    interactions = []
    
    # Pi-Pi interactions
    if aromatic_a is not None and aromatic_b is not None:
        for arom1 in aromatic_a:
            res1_num = arom1.getResnum()
            res1_name = arom1.getResname()
            chain1 = arom1.getChid()
            
            for arom2 in aromatic_b:
                res2_num = arom2.getResnum()
                res2_name = arom2.getResname()
                chain2 = arom2.getChid()
                
                # Skip same residue
                if res1_num == res2_num and chain1 == chain2:
                    continue
                
                # If analyzing within same chain, avoid duplicates
                if chain_a is None and chain1 == chain2 and res1_num >= res2_num:
                    continue
                
                dist = np.linalg.norm(arom1.getCoords() - arom2.getCoords())
                
                if dist <= cutoff:
                    interactions.append({
                        'type': 'pi-pi',
                        'res1': f"{res1_name}{res1_num}",
                        'chain1': chain1,
                        'res2': f"{res2_name}{res2_num}",
                        'chain2': chain2,
                        'distance': dist
                    })
    
    # Cation-Pi interactions
    if aromatic_a is not None and cation_b is not None:
        for arom in aromatic_a:
            arom_num = arom.getResnum()
            arom_name = arom.getResname()
            arom_chain = arom.getChid()
            
            for cation in cation_b:
                cat_num = cation.getResnum()
                cat_name = cation.getResname()
                cat_chain = cation.getChid()
                
                # Skip same residue
                if arom_num == cat_num and arom_chain == cat_chain:
                    continue
                
                dist = np.linalg.norm(arom.getCoords() - cation.getCoords())
                
                if dist <= cutoff:
                    interactions.append({
                        'type': 'cation-pi',
                        'res1': f"{arom_name}{arom_num}",
                        'chain1': arom_chain,
                        'res2': f"{cat_name}{cat_num}",
                        'chain2': cat_chain,
                        'distance': dist
                    })
    
    # Reverse direction for cation-pi if analyzing between chains
    if chain_a and chain_b and aromatic_b is not None and cation_a is not None:
        for arom in aromatic_b:
            arom_num = arom.getResnum()
            arom_name = arom.getResname()
            arom_chain = arom.getChid()
            
            for cation in cation_a:
                cat_num = cation.getResnum()
                cat_name = cation.getResname()
                cat_chain = cation.getChid()
                
                dist = np.linalg.norm(arom.getCoords() - cation.getCoords())
                
                if dist <= cutoff:
                    interactions.append({
                        'type': 'cation-pi',
                        'res1': f"{cat_name}{cat_num}",
                        'chain1': cat_chain,
                        'res2': f"{arom_name}{arom_num}",
                        'chain2': arom_chain,
                        'distance': dist
                    })
    
    print(f"Found {len(interactions)} aromatic interactions")
    return interactions


# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

def plot_intermolecular_contact_heatmap(inter_results, output_dir, chain_a, chain_b):
    """Plot heatmap of intermolecular contacts"""
    if inter_results is None:
        return
    
    print("\nGenerating intermolecular contact heatmap...")
    
    ca_a = inter_results['ca_a']
    ca_b = inter_results['ca_b']
    contacts = inter_results['contacts']
    
    resnums_a = ca_a.getResnums()
    resnums_b = ca_b.getResnums()
    
    # Create contact matrix
    contact_matrix = np.zeros((len(resnums_a), len(resnums_b)))
    
    for contact in contacts:
        idx_a = np.where(resnums_a == contact['resnum_a'])[0]
        idx_b = np.where(resnums_b == contact['resnum_b'])[0]
        if len(idx_a) > 0 and len(idx_b) > 0:
            contact_matrix[idx_a[0], idx_b[0]] = 1
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    im = ax.imshow(contact_matrix, cmap='YlOrRd', aspect='auto', interpolation='nearest')
    
    ax.set_xlabel(f'Chain {chain_b} Residue Number', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'Chain {chain_a} Residue Number', fontsize=12, fontweight='bold')
    ax.set_title(f'Intermolecular Contacts: {chain_a} ↔ {chain_b}', fontsize=14, fontweight='bold')
    
    # Set ticks
    step_a = max(1, len(resnums_a) // 20)
    step_b = max(1, len(resnums_b) // 20)
    ax.set_yticks(np.arange(0, len(resnums_a), step_a))
    ax.set_yticklabels(resnums_a[::step_a])
    ax.set_xticks(np.arange(0, len(resnums_b), step_b))
    ax.set_xticklabels(resnums_b[::step_b])
    
    plt.colorbar(im, ax=ax, label='Contact')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/intermolecular_contacts_{chain_a}_{chain_b}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/intermolecular_contacts_{chain_a}_{chain_b}.png")


def plot_contact_counts(inter_results, output_dir, chain_a, chain_b):
    """Plot contact counts per residue"""
    if inter_results is None:
        return
    
    print("\nGenerating contact count plots...")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))
    
    # Chain A contacts
    ca_a = inter_results['ca_a']
    resnums_a = ca_a.getResnums()
    counts_a = [inter_results['counts_a'].get(rn, 0) for rn in resnums_a]
    
    ax1.bar(resnums_a, counts_a, color='steelblue', edgecolor='black', linewidth=0.5)
    ax1.set_xlabel(f'Chain {chain_a} Residue Number', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Number of Contacts', fontsize=11, fontweight='bold')
    ax1.set_title(f'Interface Residues: Chain {chain_a}', fontsize=12, fontweight='bold')
    ax1.grid(alpha=0.3)
    
    # Chain B contacts
    ca_b = inter_results['ca_b']
    resnums_b = ca_b.getResnums()
    counts_b = [inter_results['counts_b'].get(rn, 0) for rn in resnums_b]
    
    ax2.bar(resnums_b, counts_b, color='coral', edgecolor='black', linewidth=0.5)
    ax2.set_xlabel(f'Chain {chain_b} Residue Number', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Number of Contacts', fontsize=11, fontweight='bold')
    ax2.set_title(f'Interface Residues: Chain {chain_b}', fontsize=12, fontweight='bold')
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/interface_contact_counts_{chain_a}_{chain_b}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/interface_contact_counts_{chain_a}_{chain_b}.png")


def plot_intramolecular_contacts(intra_results, output_dir, chain_id):
    """Plot intramolecular contact matrix and counts"""
    print(f"\nGenerating intramolecular contact plots for chain {chain_id}...")
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Contact matrix
    ax = axes[0]
    im = ax.imshow(intra_results['contact_matrix'], cmap='Blues', aspect='auto', interpolation='nearest')
    ax.set_xlabel('Residue Index', fontsize=11, fontweight='bold')
    ax.set_ylabel('Residue Index', fontsize=11, fontweight='bold')
    ax.set_title(f'Contact Matrix: Chain {chain_id}', fontsize=12, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Contact')
    
    # Contact counts
    ax = axes[1]
    resnums = intra_results['resnums']
    counts = intra_results['contact_counts']
    
    ax.bar(resnums, counts, color='steelblue', edgecolor='black', linewidth=0.5)
    ax.set_xlabel('Residue Number', fontsize=11, fontweight='bold')
    ax.set_ylabel('Number of Contacts', fontsize=11, fontweight='bold')
    ax.set_title(f'Contact Count per Residue: Chain {chain_id}', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/intramolecular_contacts_chain_{chain_id}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/intramolecular_contacts_chain_{chain_id}.png")


def plot_sasa_hydrophobicity(sasa_results, output_dir, chain_id):
    """Plot SASA and hydrophobicity profiles"""
    print(f"\nGenerating SASA and hydrophobicity plots...")
    
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))
    
    resnums = sasa_results['resnums']
    sasa = sasa_results['sasa']
    hydro = sasa_results['hydrophobicity']
    exposure = sasa_results['exposure_class']
    
    # SASA profile
    ax = axes[0]
    colors = ['blue' if e == 0 else 'orange' if e == 1 else 'red' for e in exposure]
    ax.bar(resnums, sasa, color=colors, edgecolor='black', linewidth=0.3)
    ax.set_xlabel('Residue Number', fontsize=11, fontweight='bold')
    ax.set_ylabel('SASA (Ų)', fontsize=11, fontweight='bold')
    ax.set_title(f'Solvent Accessible Surface Area: Chain {chain_id}', fontsize=12, fontweight='bold')
    ax.axhline(y=10, color='gray', linestyle='--', alpha=0.5, label='Core threshold')
    ax.axhline(y=40, color='gray', linestyle='--', alpha=0.5, label='Surface threshold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Hydrophobicity
    ax = axes[1]
    colors = ['red' if h < 0 else 'blue' for h in hydro]
    ax.bar(resnums, hydro, color=colors, edgecolor='black', linewidth=0.3)
    ax.set_xlabel('Residue Number', fontsize=11, fontweight='bold')
    ax.set_ylabel('Hydrophobicity Score', fontsize=11, fontweight='bold')
    ax.set_title(f'Hydrophobicity Profile (Kyte-Doolittle): Chain {chain_id}', fontsize=12, fontweight='bold')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.grid(alpha=0.3)
    
    # Combined view
    ax = axes[2]
    ax2 = ax.twinx()
    
    # SASA on left axis
    ax.fill_between(resnums, sasa, alpha=0.3, color='blue', label='SASA')
    ax.plot(resnums, sasa, color='blue', linewidth=2)
    ax.set_xlabel('Residue Number', fontsize=11, fontweight='bold')
    ax.set_ylabel('SASA (Ų)', fontsize=11, fontweight='bold', color='blue')
    ax.tick_params(axis='y', labelcolor='blue')
    
    # Hydrophobicity on right axis
    ax2.plot(resnums, hydro, color='red', linewidth=2, label='Hydrophobicity')
    ax2.set_ylabel('Hydrophobicity Score', fontsize=11, fontweight='bold', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    ax.set_title(f'SASA vs Hydrophobicity: Chain {chain_id}', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/sasa_hydrophobicity_chain_{chain_id}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/sasa_hydrophobicity_chain_{chain_id}.png")


def plot_noncovalent_interactions(hbonds, salt_bridges, aromatic, output_dir, chain_a=None, chain_b=None):
    """Plot summary of noncovalent interactions"""
    print("\nGenerating noncovalent interaction summary...")
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Count interactions by type
    interaction_types = ['H-bonds', 'Salt Bridges', 'Aromatic']
    counts = [len(hbonds), len(salt_bridges), len(aromatic)]
    colors = ['skyblue', 'coral', 'lightgreen']
    
    ax = axes[0]
    ax.bar(interaction_types, counts, color=colors, edgecolor='black', linewidth=2)
    ax.set_ylabel('Number of Interactions', fontsize=11, fontweight='bold')
    ax.set_title('Noncovalent Interaction Summary', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    # H-bond distance distribution
    ax = axes[1]
    if len(hbonds) > 0:
        distances = [hb['distance'] for hb in hbonds]
        ax.hist(distances, bins=20, color='skyblue', edgecolor='black', alpha=0.7)
        ax.axvline(x=np.mean(distances), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(distances):.2f} Å')
        ax.legend()
    ax.set_xlabel('Distance (Å)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax.set_title('H-bond Distance Distribution', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    # Salt bridge distance distribution
    ax = axes[2]
    if len(salt_bridges) > 0:
        distances = [sb['distance'] for sb in salt_bridges]
        ax.hist(distances, bins=20, color='coral', edgecolor='black', alpha=0.7)
        ax.axvline(x=np.mean(distances), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(distances):.2f} Å')
        ax.legend()
    ax.set_xlabel('Distance (Å)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax.set_title('Salt Bridge Distance Distribution', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    
    suffix = f"_{chain_a}_{chain_b}" if chain_a and chain_b else "_all"
    plt.savefig(f'{output_dir}/noncovalent_interactions{suffix}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/noncovalent_interactions{suffix}.png")


# ============================================================================
# OUTPUT FUNCTIONS
# ============================================================================

def save_intermolecular_results(inter_results, output_dir, chain_a, chain_b):
    """Save intermolecular contact results to file"""
    if inter_results is None:
        return
    
    print(f"\nSaving intermolecular contact results...")
    
    filename = f"{output_dir}/intermolecular_contacts_{chain_a}_{chain_b}.txt"
    
    with open(filename, 'w') as f:
        f.write("INTERMOLECULAR CONTACT ANALYSIS\n")
        f.write("="*70 + "\n\n")
        f.write(f"Chain A: {chain_a}\n")
        f.write(f"Chain B: {chain_b}\n")
        f.write(f"Total contacts: {len(inter_results['contacts'])}\n")
        f.write(f"Interface residues (chain {chain_a}): {len(inter_results['counts_a'])}\n")
        f.write(f"Interface residues (chain {chain_b}): {len(inter_results['counts_b'])}\n\n")
        
        f.write("CONTACT LIST\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Chain A Residue':<20} {'Chain B Residue':<20} {'Distance (Å)':<15}\n")
        f.write("-"*70 + "\n")
        
        for contact in sorted(inter_results['contacts'], key=lambda x: x['distance']):
            res_a = f"{contact['resname_a']}{contact['resnum_a']}"
            res_b = f"{contact['resname_b']}{contact['resnum_b']}"
            f.write(f"{res_a:<20} {res_b:<20} {contact['distance']:<15.3f}\n")
        
        f.write("\n\nINTERFACE RESIDUES (Chain A)\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Residue':<15} {'Contact Count':<15}\n")
        f.write("-"*70 + "\n")
        for resnum, count in sorted(inter_results['counts_a'].items(), key=lambda x: x[1], reverse=True):
            f.write(f"{resnum:<15} {count:<15}\n")
        
        f.write("\n\nINTERFACE RESIDUES (Chain B)\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Residue':<15} {'Contact Count':<15}\n")
        f.write("-"*70 + "\n")
        for resnum, count in sorted(inter_results['counts_b'].items(), key=lambda x: x[1], reverse=True):
            f.write(f"{resnum:<15} {count:<15}\n")
    
    print(f"Saved: {filename}")


def save_intramolecular_results(intra_results, output_dir, chain_id):
    """Save intramolecular contact results to file"""
    print(f"\nSaving intramolecular contact results for chain {chain_id}...")
    
    filename = f"{output_dir}/intramolecular_contacts_chain_{chain_id}.txt"
    
    with open(filename, 'w') as f:
        f.write(f"INTRAMOLECULAR CONTACT ANALYSIS: Chain {chain_id}\n")
        f.write("="*70 + "\n\n")
        f.write(f"Total contacts: {len(intra_results['contact_list'])}\n")
        f.write(f"Average contacts per residue: {np.mean(intra_results['contact_counts']):.2f}\n")
        f.write(f"Average neighbors per residue: {np.mean(intra_results['neighbor_counts']):.2f}\n\n")
        
        f.write("CONTACT LIST\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Residue i':<15} {'Residue j':<15} {'Distance (Å)':<15}\n")
        f.write("-"*70 + "\n")
        
        for contact in sorted(intra_results['contact_list'], key=lambda x: x['distance']):
            f.write(f"{contact['res_i']:<15} {contact['res_j']:<15} {contact['distance']:<15.3f}\n")
        
        f.write("\n\nPER-RESIDUE CONTACT PROFILE\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Residue':<15} {'Contacts':<15} {'Neighbors':<15}\n")
        f.write("-"*70 + "\n")
        
        for i, (resnum, resname) in enumerate(zip(intra_results['resnums'], intra_results['resnames'])):
            res_label = f"{resname}{resnum}"
            contacts = int(intra_results['contact_counts'][i])
            neighbors = int(intra_results['neighbor_counts'][i])
            f.write(f"{res_label:<15} {contacts:<15} {neighbors:<15}\n")
    
    print(f"Saved: {filename}")


def save_sasa_results(sasa_results, output_dir, chain_id):
    """Save SASA and hydrophobicity results to file"""
    print(f"\nSaving SASA and hydrophobicity results...")
    
    filename = f"{output_dir}/sasa_hydrophobicity_chain_{chain_id}.txt"
    
    exposure_names = {0: 'Core', 1: 'Partial', 2: 'Surface'}
    
    with open(filename, 'w') as f:
        f.write(f"SASA AND HYDROPHOBICITY ANALYSIS: Chain {chain_id}\n")
        f.write("="*70 + "\n\n")
        
        n_core = np.sum(sasa_results['exposure_class'] == 0)
        n_partial = np.sum(sasa_results['exposure_class'] == 1)
        n_surface = np.sum(sasa_results['exposure_class'] == 2)
        
        f.write(f"Core residues (SASA < 10): {n_core}\n")
        f.write(f"Partially exposed (10-40): {n_partial}\n")
        f.write(f"Surface residues (> 40): {n_surface}\n\n")
        
        f.write("PER-RESIDUE PROFILE\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Residue':<15} {'SASA (Ų)':<15} {'Hydrophobicity':<15} {'Exposure':<15}\n")
        f.write("-"*70 + "\n")
        
        for i, (resnum, resname) in enumerate(zip(sasa_results['resnums'], sasa_results['resnames'])):
            res_label = f"{resname}{resnum}"
            sasa = sasa_results['sasa'][i]
            hydro = sasa_results['hydrophobicity'][i]
            exposure = exposure_names[sasa_results['exposure_class'][i]]
            f.write(f"{res_label:<15} {sasa:<15.2f} {hydro:<15.2f} {exposure:<15}\n")
    
    print(f"Saved: {filename}")


def save_noncovalent_results(hbonds, salt_bridges, aromatic, output_dir, chain_a=None, chain_b=None):
    """Save noncovalent interaction results to file"""
    print("\nSaving noncovalent interaction results...")
    
    suffix = f"_{chain_a}_{chain_b}" if chain_a and chain_b else "_all"
    filename = f"{output_dir}/noncovalent_interactions{suffix}.txt"
    
    with open(filename, 'w') as f:
        f.write("NONCOVALENT INTERACTION ANALYSIS\n")
        f.write("="*70 + "\n\n")
        
        if chain_a and chain_b:
            f.write(f"Analyzing interactions between Chain {chain_a} and Chain {chain_b}\n\n")
        else:
            f.write("Analyzing all interactions\n\n")
        
        f.write(f"Total hydrogen bonds: {len(hbonds)}\n")
        f.write(f"Total salt bridges: {len(salt_bridges)}\n")
        f.write(f"Total aromatic interactions: {len(aromatic)}\n\n")
        
        # Hydrogen bonds
        f.write("HYDROGEN BONDS\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Donor':<20} {'Chain':<10} {'Acceptor':<20} {'Chain':<10} {'Distance (Å)':<15}\n")
        f.write("-"*70 + "\n")
        
        for hb in sorted(hbonds, key=lambda x: x['distance']):
            f.write(f"{hb['donor']:<20} {hb['donor_chain']:<10} {hb['acceptor']:<20} "
                   f"{hb['acceptor_chain']:<10} {hb['distance']:<15.3f}\n")
        
        # Salt bridges
        f.write("\n\nSALT BRIDGES\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Positive':<20} {'Chain':<10} {'Negative':<20} {'Chain':<10} {'Distance (Å)':<15}\n")
        f.write("-"*70 + "\n")
        
        for sb in sorted(salt_bridges, key=lambda x: x['distance']):
            f.write(f"{sb['positive']:<20} {sb['pos_chain']:<10} {sb['negative']:<20} "
                   f"{sb['neg_chain']:<10} {sb['distance']:<15.3f}\n")
        
        # Aromatic interactions
        f.write("\n\nAROMATIC INTERACTIONS\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Type':<15} {'Residue 1':<20} {'Chain':<10} {'Residue 2':<20} {'Chain':<10} {'Distance (Å)':<15}\n")
        f.write("-"*70 + "\n")
        
        for ar in sorted(aromatic, key=lambda x: x['distance']):
            f.write(f"{ar['type']:<15} {ar['res1']:<20} {ar['chain1']:<10} {ar['res2']:<20} "
                   f"{ar['chain2']:<10} {ar['distance']:<15.3f}\n")
    
    print(f"Saved: {filename}")


def save_pymol_visualization_scripts(inter_results, intra_results, sasa_results, 
                                     hbonds, salt_bridges, output_dir, chain_a=None, chain_b=None):
    """Generate PyMOL scripts for visualization"""
    print("\nGenerating PyMOL visualization scripts...")
    
    # Intermolecular contacts script
    if inter_results and chain_a and chain_b:
        filename = f"{output_dir}/visualize_interface_{chain_a}_{chain_b}.pml"
        with open(filename, 'w') as f:
            f.write("# PyMOL script to visualize interface\n\n")
            f.write("hide everything\n")
            f.write("show cartoon\n")
            f.write("set cartoon_fancy_helices, 1\n\n")
            
            f.write(f"# Color chains\n")
            f.write(f"color lightblue, chain {chain_a}\n")
            f.write(f"color lightpink, chain {chain_b}\n\n")
            
            f.write("# Highlight interface residues\n")
            if len(inter_results['counts_a']) > 0:
                res_list_a = '+'.join([str(rn) for rn in inter_results['counts_a'].keys()])
                f.write(f"select interface_A, chain {chain_a} and resi {res_list_a}\n")
                f.write("color blue, interface_A\n")
                f.write("show sticks, interface_A\n\n")
            
            if len(inter_results['counts_b']) > 0:
                res_list_b = '+'.join([str(rn) for rn in inter_results['counts_b'].keys()])
                f.write(f"select interface_B, chain {chain_b} and resi {res_list_b}\n")
                f.write("color red, interface_B\n")
                f.write("show sticks, interface_B\n\n")
            
            f.write("bg_color white\n")
            f.write("set ray_shadows, 0\n")
        
        print(f"Saved: {filename}")
    
    # SASA visualization script
    if sasa_results:
        chain_id = chain_a if chain_a else 'A'
        filename = f"{output_dir}/visualize_sasa_chain_{chain_id}.pml"
        with open(filename, 'w') as f:
            f.write("# PyMOL script to visualize SASA and exposure\n\n")
            f.write("hide everything\n")
            f.write("show cartoon\n")
            f.write("set cartoon_fancy_helices, 1\n\n")
            
            f.write("# Color by exposure\n")
            exposure = sasa_results['exposure_class']
            resnums = sasa_results['resnums']
            
            core_res = [str(rn) for i, rn in enumerate(resnums) if exposure[i] == 0]
            partial_res = [str(rn) for i, rn in enumerate(resnums) if exposure[i] == 1]
            surface_res = [str(rn) for i, rn in enumerate(resnums) if exposure[i] == 2]
            
            if core_res:
                f.write(f"select core, resi {'+'.join(core_res)}\n")
                f.write("color blue, core\n\n")
            
            if partial_res:
                f.write(f"select partial, resi {'+'.join(partial_res)}\n")
                f.write("color orange, partial\n\n")
            
            if surface_res:
                f.write(f"select surface, resi {'+'.join(surface_res)}\n")
                f.write("color red, surface\n\n")
            
            f.write("bg_color white\n")
            f.write("set ray_shadows, 0\n")
        
        print(f"Saved: {filename}")
    
    # Noncovalent interactions script
    if hbonds or salt_bridges:
        suffix = f"_{chain_a}_{chain_b}" if chain_a and chain_b else "_all"
        filename = f"{output_dir}/visualize_interactions{suffix}.pml"
        with open(filename, 'w') as f:
            f.write("# PyMOL script to visualize noncovalent interactions\n\n")
            f.write("hide everything\n")
            f.write("show cartoon\n")
            f.write("set cartoon_fancy_helices, 1\n\n")
            
            # Highlight residues involved in H-bonds
            if hbonds:
                hb_residues = set()
                for hb in hbonds:
                    # Extract residue numbers from strings like "ALA123"
                    donor_num = ''.join(filter(str.isdigit, hb['donor']))
                    acceptor_num = ''.join(filter(str.isdigit, hb['acceptor']))
                    hb_residues.add(donor_num)
                    hb_residues.add(acceptor_num)
                
                if hb_residues:
                    f.write(f"select hbond_residues, resi {'+'.join(hb_residues)}\n")
                    f.write("color cyan, hbond_residues\n")
                    f.write("show sticks, hbond_residues\n\n")
            
            # Highlight residues involved in salt bridges
            if salt_bridges:
                sb_residues = set()
                for sb in salt_bridges:
                    pos_num = ''.join(filter(str.isdigit, sb['positive']))
                    neg_num = ''.join(filter(str.isdigit, sb['negative']))
                    sb_residues.add(pos_num)
                    sb_residues.add(neg_num)
                
                if sb_residues:
                    f.write(f"select saltbridge_residues, resi {'+'.join(sb_residues)}\n")
                    f.write("color yellow, saltbridge_residues\n")
                    f.write("show sticks, saltbridge_residues\n\n")
            
            f.write("bg_color white\n")
            f.write("set ray_shadows, 0\n")
        
        print(f"Saved: {filename}")


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='ProDy Contact, SASA, and Stability Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze all chains
  python3 prody_contact_stability_analysis.py --pdb complex.pdb

  # Analyze interface between specific chains
  python3 prody_contact_stability_analysis.py --pdb complex.pdb --chain-a A --chain-b B

  # Full analysis for single chain
  python3 prody_contact_stability_analysis.py --pdb protein.pdb --chain A

  # Custom cutoffs
  python3 prody_contact_stability_analysis.py --pdb complex.pdb --chain-a A --chain-b B --contact-cutoff 4.5
        """)
    
    parser.add_argument('--pdb', required=True, help='Input PDB file')
    parser.add_argument('--chain-a', '--chain', dest='chain_a', default=None,
                       help='First chain or single chain to analyze')
    parser.add_argument('--chain-b', default=None,
                       help='Second chain for interface analysis (optional)')
    parser.add_argument('--contact-cutoff', type=float, default=5.0,
                       help='Contact distance cutoff in Å (default: 5.0)')
    parser.add_argument('--neighbor-cutoff', type=float, default=8.0,
                       help='Neighbor distance cutoff in Å (default: 8.0)')
    parser.add_argument('--hbond-cutoff', type=float, default=3.5,
                       help='H-bond distance cutoff in Å (default: 3.5)')
    parser.add_argument('--saltbridge-cutoff', type=float, default=4.0,
                       help='Salt bridge distance cutoff in Å (default: 4.0)')
    parser.add_argument('--skip-sasa', action='store_true',
                       help='Skip SASA calculation (faster)')
    
    args = parser.parse_args()
    
    # Setup
    pdb_name = os.path.splitext(os.path.basename(args.pdb))[0]
    output_dir = setup_output_directory(pdb_name)
    
    # Load structure
    print(f"\nLoading PDB: {args.pdb}")
    structure = parsePDB(args.pdb)
    
    if structure is None:
        print(f"Error: Could not parse PDB file: {args.pdb}")
        sys.exit(1)
    
    # Analyze chain structure
    chain_info = analyze_chain_structure(structure)
    print_chain_summary(chain_info)
    
    # Determine analysis mode
    if args.chain_a and args.chain_b:
        # Interface analysis mode
        print(f"\nMode: Interface analysis between chains {args.chain_a} and {args.chain_b}")
        
        # 1. Intermolecular contacts
        inter_results = calculate_intermolecular_contacts(
            structure, args.chain_a, args.chain_b, cutoff=args.contact_cutoff
        )
        
        # 2. Noncovalent interactions (interface)
        hbonds = find_hydrogen_bonds(structure, args.chain_a, args.chain_b, 
                                     dist_cutoff=args.hbond_cutoff)
        salt_bridges = find_salt_bridges(structure, args.chain_a, args.chain_b,
                                        cutoff=args.saltbridge_cutoff)
        aromatic = find_aromatic_interactions(structure, args.chain_a, args.chain_b)
        
        # Generate plots
        plot_intermolecular_contact_heatmap(inter_results, output_dir, args.chain_a, args.chain_b)
        plot_contact_counts(inter_results, output_dir, args.chain_a, args.chain_b)
        plot_noncovalent_interactions(hbonds, salt_bridges, aromatic, output_dir, 
                                     args.chain_a, args.chain_b)
        
        # Save results
        save_intermolecular_results(inter_results, output_dir, args.chain_a, args.chain_b)
        save_noncovalent_results(hbonds, salt_bridges, aromatic, output_dir,
                               args.chain_a, args.chain_b)
        
        # Optional: analyze individual chains
        for chain_id in [args.chain_a, args.chain_b]:
            if chain_id in chain_info:
                chain_atoms = structure.select(f'chain {chain_id}')
                
                # Intramolecular contacts
                intra_results = calculate_intramolecular_contacts(
                    chain_atoms, chain_id, 
                    contact_cutoff=args.contact_cutoff,
                    neighbor_cutoff=args.neighbor_cutoff
                )
                plot_intramolecular_contacts(intra_results, output_dir, chain_id)
                save_intramolecular_results(intra_results, output_dir, chain_id)
                
                # SASA analysis
                if not args.skip_sasa:
                    sasa_results = calculate_sasa_hydrophobicity(structure, chain_id)
                    plot_sasa_hydrophobicity(sasa_results, output_dir, chain_id)
                    save_sasa_results(sasa_results, output_dir, chain_id)
        
        # PyMOL scripts
        sasa_results = None if args.skip_sasa else calculate_sasa_hydrophobicity(structure, args.chain_a)
        save_pymol_visualization_scripts(inter_results, None, sasa_results,
                                        hbonds, salt_bridges, output_dir,
                                        args.chain_a, args.chain_b)
        
    elif args.chain_a:
        # Single chain analysis mode
        print(f"\nMode: Single chain analysis (chain {args.chain_a})")
        
        if args.chain_a not in chain_info:
            print(f"Error: Chain {args.chain_a} not found in structure")
            sys.exit(1)
        
        chain_atoms = structure.select(f'chain {args.chain_a}')
        
        # 1. Intramolecular contacts
        intra_results = calculate_intramolecular_contacts(
            chain_atoms, args.chain_a,
            contact_cutoff=args.contact_cutoff,
            neighbor_cutoff=args.neighbor_cutoff
        )
        
        # 2. SASA analysis
        if not args.skip_sasa:
            sasa_results = calculate_sasa_hydrophobicity(structure, args.chain_a)
        else:
            sasa_results = None
        
        # 3. Noncovalent interactions (within chain)
        chain_struct = structure.select(f'chain {args.chain_a}')
        hbonds = find_hydrogen_bonds(chain_struct)
        salt_bridges = find_salt_bridges(chain_struct)
        aromatic = find_aromatic_interactions(chain_struct)
        
        # Generate plots
        plot_intramolecular_contacts(intra_results, output_dir, args.chain_a)
        if sasa_results:
            plot_sasa_hydrophobicity(sasa_results, output_dir, args.chain_a)
        plot_noncovalent_interactions(hbonds, salt_bridges, aromatic, output_dir)
        
        # Save results
        save_intramolecular_results(intra_results, output_dir, args.chain_a)
        if sasa_results:
            save_sasa_results(sasa_results, output_dir, args.chain_a)
        save_noncovalent_results(hbonds, salt_bridges, aromatic, output_dir)
        
        # PyMOL scripts
        save_pymol_visualization_scripts(None, intra_results, sasa_results,
                                        hbonds, salt_bridges, output_dir,
                                        args.chain_a, None)
        
    else:
        # Analyze all chains
        print("\nMode: Multi-chain analysis (all chains)")
        
        for chain_id in sorted(chain_info.keys()):
            print(f"\n{'='*60}")
            print(f"Analyzing Chain {chain_id}")
            print(f"{'='*60}")
            
            chain_atoms = structure.select(f'chain {chain_id}')
            
            # Intramolecular contacts
            intra_results = calculate_intramolecular_contacts(
                chain_atoms, chain_id,
                contact_cutoff=args.contact_cutoff,
                neighbor_cutoff=args.neighbor_cutoff
            )
            plot_intramolecular_contacts(intra_results, output_dir, chain_id)
            save_intramolecular_results(intra_results, output_dir, chain_id)
            
            # SASA analysis
            if not args.skip_sasa:
                sasa_results = calculate_sasa_hydrophobicity(structure, chain_id)
                plot_sasa_hydrophobicity(sasa_results, output_dir, chain_id)
                save_sasa_results(sasa_results, output_dir, chain_id)
        
        # Analyze all pairwise interfaces
        print(f"\n{'='*60}")
        print("PAIRWISE INTERFACE ANALYSIS")
        print(f"{'='*60}")
        
        chain_ids = sorted(chain_info.keys())
        for i, chain_a in enumerate(chain_ids):
            for chain_b in chain_ids[i+1:]:
                print(f"\nAnalyzing interface: {chain_a} ↔ {chain_b}")
                
                # Intermolecular contacts
                inter_results = calculate_intermolecular_contacts(
                    structure, chain_a, chain_b, cutoff=args.contact_cutoff
                )
                
                if inter_results and len(inter_results['contacts']) > 0:
                    # Noncovalent interactions
                    hbonds = find_hydrogen_bonds(structure, chain_a, chain_b,
                                                dist_cutoff=args.hbond_cutoff)
                    salt_bridges = find_salt_bridges(structure, chain_a, chain_b,
                                                    cutoff=args.saltbridge_cutoff)
                    aromatic = find_aromatic_interactions(structure, chain_a, chain_b)
                    
                    # Generate plots
                    plot_intermolecular_contact_heatmap(inter_results, output_dir, chain_a, chain_b)
                    plot_contact_counts(inter_results, output_dir, chain_a, chain_b)
                    plot_noncovalent_interactions(hbonds, salt_bridges, aromatic, output_dir,
                                                 chain_a, chain_b)
                    
                    # Save results
                    save_intermolecular_results(inter_results, output_dir, chain_a, chain_b)
                    save_noncovalent_results(hbonds, salt_bridges, aromatic, output_dir,
                                           chain_a, chain_b)
                    
                    # PyMOL scripts
                    save_pymol_visualization_scripts(inter_results, None, None,
                                                    hbonds, salt_bridges, output_dir,
                                                    chain_a, chain_b)
                else:
                    print(f"  No significant interface found between {chain_a} and {chain_b}")
    
    # Summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nAll results saved to: {output_dir}/")
    print("\nGenerated files:")
    print("  - Contact matrices and heatmaps")
    print("  - SASA and hydrophobicity profiles")
    print("  - Noncovalent interaction summaries")
    print("  - PyMOL visualization scripts")
    print("\nTo visualize in PyMOL:")
    print(f"  pymol {args.pdb} {output_dir}/visualize_*.pml")
    print("="*70)


if __name__ == '__main__':
    main()

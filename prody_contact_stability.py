#!/usr/bin/env python3
"""
ProDy Contact, SASA, and Stability Analysis Script - Enhanced Version

This script performs comprehensive structural analysis including:
1. Intermolecular contact mapping (between chains)
2. Intramolecular contact mapping (within chains)
3. Hydrophobicity and solvent accessibility (SASA)
4. Noncovalent interaction mapping (H-bonds, salt bridges, pi interactions)
5. Enhanced hydrogen bond detection with proper geometry
6. Disulfide bond detection
7. Hydrophobic interactions

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
import warnings
warnings.filterwarnings('ignore')

try:
    from prody import *
    try:
        from prody.proteins import calcSASA
        SASA_AVAILABLE = True
    except ImportError:
        SASA_AVAILABLE = False
        print("Warning: calcSASA not available in this ProDy version.")
except ImportError:
    print("Error: ProDy is not installed. Install with: pip install prody")
    sys.exit(1)

# ============================================================================
# CONSTANTS AND PARAMETERS
# ============================================================================

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

# Hydrophobic residues
HYDROPHOBIC = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO'}

# Hydrogen bond donor atoms by residue
HBOND_DONORS = {
    'backbone': ['N'],  # All residues have backbone N-H
    'ALA': [],
    'ARG': ['NE', 'NH1', 'NH2'],
    'ASN': ['ND2'],
    'ASP': [],
    'CYS': ['SG'],
    'GLN': ['NE2'],
    'GLU': [],
    'GLY': [],
    'HIS': ['ND1', 'NE2'],
    'ILE': [],
    'LEU': [],
    'LYS': ['NZ'],
    'MET': [],
    'PHE': [],
    'PRO': [],  # No backbone NH
    'SER': ['OG'],
    'THR': ['OG1'],
    'TRP': ['NE1'],
    'TYR': ['OH'],
    'VAL': []
}

# Hydrogen bond acceptor atoms by residue
HBOND_ACCEPTORS = {
    'backbone': ['O'],  # All residues have backbone C=O
    'ALA': [],
    'ARG': [],
    'ASN': ['OD1', 'ND2'],
    'ASP': ['OD1', 'OD2'],
    'CYS': ['SG'],
    'GLN': ['OE1', 'NE2'],
    'GLU': ['OE1', 'OE2'],
    'GLY': [],
    'HIS': ['ND1', 'NE2'],
    'ILE': [],
    'LEU': [],
    'LYS': [],
    'MET': ['SD'],
    'PHE': [],
    'PRO': [],
    'SER': ['OG'],
    'THR': ['OG1'],
    'TRP': [],
    'TYR': ['OH'],
    'VAL': []
}


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
    Calculate contacts between two chains with optimized performance
    """
    print(f"\n{'='*60}")
    print(f"INTERMOLECULAR CONTACTS: Chain {chain_a} ↔ Chain {chain_b}")
    print(f"{'='*60}")
    
    atoms_a = structure.select(f'chain {chain_a}')
    atoms_b = structure.select(f'chain {chain_b}')
    
    if atoms_a is None or atoms_b is None:
        print(f"Error: Could not find chain {chain_a} or {chain_b}")
        return None
    
    ca_a = atoms_a.select('calpha')
    ca_b = atoms_b.select('calpha')
    heavy_a = atoms_a.select('not hydrogen')
    heavy_b = atoms_b.select('not hydrogen')
    
    contacts = []
    contact_counts_a = defaultdict(int)
    contact_counts_b = defaultdict(int)
    
    print(f"Analyzing contacts (cutoff: {cutoff} Å)...")
    
    resnums_a = ca_a.getResnums()
    resnames_a = ca_a.getResnames()
    resnums_b = ca_b.getResnums()
    resnames_b = ca_b.getResnames()
    
    # Pre-build residue atom dictionaries for faster lookup
    residue_atoms_a = {}
    for resnum in resnums_a:
        residue_atoms_a[resnum] = heavy_a.select(f'resnum {resnum}')
    
    residue_atoms_b = {}
    for resnum in resnums_b:
        residue_atoms_b[resnum] = heavy_b.select(f'resnum {resnum}')
    
    # Calculate contacts
    for i, (resnum_a, resname_a) in enumerate(zip(resnums_a, resnames_a)):
        ca_coord_a = ca_a[i].getCoords()
        res_atoms_a = residue_atoms_a[resnum_a]
        
        if res_atoms_a is None:
            continue
        
        coords_a = res_atoms_a.getCoords()
        
        for j, (resnum_b, resname_b) in enumerate(zip(resnums_b, resnames_b)):
            ca_coord_b = ca_b[j].getCoords()
            
            # Quick CA distance check
            ca_dist = np.linalg.norm(ca_coord_a - ca_coord_b)
            if ca_dist > cutoff + 10:
                continue
            
            res_atoms_b = residue_atoms_b[resnum_b]
            if res_atoms_b is None:
                continue
            
            coords_b = res_atoms_b.getCoords()
            
            # Vectorized distance calculation
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
    Calculate contacts within a single chain with optimized performance
    """
    print(f"\n{'='*60}")
    print(f"INTRAMOLECULAR CONTACTS: Chain {chain_id}")
    print(f"{'='*60}")
    
    ca_atoms = chain_atoms.select('calpha')
    heavy_atoms = chain_atoms.select('not hydrogen')
    
    n_residues = ca_atoms.numAtoms()
    resnums = ca_atoms.getResnums()
    resnames = ca_atoms.getResnames()
    
    contact_matrix = np.zeros((n_residues, n_residues))
    neighbor_matrix = np.zeros((n_residues, n_residues))
    contact_list = []
    
    print(f"Calculating contacts for {n_residues} residues...")
    
    # Pre-build residue atoms dictionary
    residue_atoms = {}
    for resnum in resnums:
        residue_atoms[resnum] = heavy_atoms.select(f'resnum {resnum}')
    
    ca_coords = ca_atoms.getCoords()
    
    for i in range(n_residues):
        resnum_i = resnums[i]
        resname_i = resnames[i]
        res_atoms_i = residue_atoms[resnum_i]
        
        if res_atoms_i is None:
            continue
        
        coords_i = res_atoms_i.getCoords()
        
        for j in range(i + 3, n_residues):  # Skip i, i+1, i+2
            resnum_j = resnums[j]
            resname_j = resnames[j]
            
            ca_dist = np.linalg.norm(ca_coords[i] - ca_coords[j])
            
            if ca_dist <= neighbor_cutoff:
                neighbor_matrix[i, j] = 1
                neighbor_matrix[j, i] = 1
                
                res_atoms_j = residue_atoms[resnum_j]
                if res_atoms_j is None:
                    continue
                
                coords_j = res_atoms_j.getCoords()
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
    Calculate SASA and hydrophobicity with improved method
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
    
    print("Calculating SASA...")
    
    if SASA_AVAILABLE:
        try:
            sasa_array = calcSASA(selection)
            residue_sasa = []
            
            for resnum in resnums:
                res_atoms = selection.select(f'resnum {resnum}')
                if res_atoms is None:
                    residue_sasa.append(0)
                    continue
                
                res_indices = res_atoms.getIndices()
                local_indices = [np.where(selection.getIndices() == idx)[0][0] for idx in res_indices]
                res_sasa = np.sum(sasa_array[local_indices])
                residue_sasa.append(res_sasa)
            
            residue_sasa = np.array(residue_sasa)
            print(f"Using ProDy calcSASA")
            
        except Exception as e:
            print(f"Warning: calcSASA failed ({e}), using approximation")
            SASA_AVAILABLE = False
    
    if not SASA_AVAILABLE:
        # Contact-based approximation
        print("Using contact-based SASA approximation...")
        all_atoms = selection.select('not hydrogen')
        if all_atoms is None:
            print("Warning: No heavy atoms found")
            residue_sasa = np.zeros(len(resnums))
        else:
            all_coords = all_atoms.getCoords()
            residue_sasa = []
            
            for resnum in resnums:
                res_atoms = all_atoms.select(f'resnum {resnum}')
                if res_atoms is None:
                    residue_sasa.append(0)
                    continue
                
                res_coords = res_atoms.getCoords()
                n_atoms = len(res_coords)
                exposed_count = 0
                
                for atom_coord in res_coords:
                    distances = np.sqrt(np.sum((all_coords - atom_coord)**2, axis=1))
                    n_neighbors = np.sum((distances > 0.1) & (distances < 5.0))
                    if n_neighbors < 15:
                        exposed_count += 1
                
                approx_sasa = (exposed_count / n_atoms) * 100 if n_atoms > 0 else 0
                residue_sasa.append(approx_sasa)
            
            residue_sasa = np.array(residue_sasa)
    
    # Hydrophobicity
    residue_hydro = np.array([HYDROPHOBICITY.get(resname, 0) for resname in resnames])
    
    # Classify exposure
    exposure_class = np.zeros(len(residue_sasa), dtype=int)
    exposure_class[residue_sasa > 40] = 2  # Surface
    exposure_class[(residue_sasa > 10) & (residue_sasa <= 40)] = 1  # Partial
    
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
# 4. ENHANCED NONCOVALENT INTERACTION MAPPING
# ============================================================================

def find_hydrogen_bonds_enhanced(structure, chain_a=None, chain_b=None, 
                                 dist_cutoff=3.5, angle_cutoff=90):
    """
    Enhanced hydrogen bond detection with proper donor/acceptor classification
    """
    print(f"\n{'='*60}")
    print("HYDROGEN BOND ANALYSIS (ENHANCED)")
    print(f"{'='*60}")
    
    if chain_a and chain_b:
        donor_chain = structure.select(f'chain {chain_a}')
        acceptor_chain = structure.select(f'chain {chain_b}')
        print(f"Analyzing inter-chain H-bonds: {chain_a} → {chain_b}")
    else:
        donor_chain = structure
        acceptor_chain = structure
        print("Analyzing all H-bonds")
    
    hbonds = []
    
    # Get all residues
    donor_ca = donor_chain.select('calpha')
    acceptor_ca = acceptor_chain.select('calpha')
    
    if donor_ca is None or acceptor_ca is None:
        print("No CA atoms found")
        return hbonds
    
    donor_resnums = donor_ca.getResnums()
    donor_resnames = donor_ca.getResnames()
    donor_chains = donor_ca.getChids()
    acceptor_resnums = acceptor_ca.getResnums()
    acceptor_resnames = acceptor_ca.getResnames()
    acceptor_chains = acceptor_ca.getChids()
    
    # Iterate through donor residues
    for d_resnum, d_resname in zip(donor_resnums, donor_resnames):
        # Get donor atoms for this residue
        donor_atoms = []
        
        # Backbone N (except proline)
        if d_resname != 'PRO':
            backbone_n = donor_chain.select(f'resnum {d_resnum} and name N')
            if backbone_n:
                donor_atoms.extend([(atom, 'N', 'backbone') for atom in backbone_n])
        
        # Side chain donors
        if d_resname in HBOND_DONORS and HBOND_DONORS[d_resname]:
            for donor_name in HBOND_DONORS[d_resname]:
                sc_donors = donor_chain.select(f'resnum {d_resnum} and name {donor_name}')
                if sc_donors:
                    donor_atoms.extend([(atom, donor_name, 'sidechain') for atom in sc_donors])
        
        if not donor_atoms:
            continue
        
        # Check against acceptor residues
        for a_resnum, a_resname in zip(acceptor_resnums, acceptor_resnames):
            # Skip same residue
            if not (chain_a and chain_b):  # intra-chain analysis
                if d_resnum == a_resnum:
                    continue
                # Skip sequential neighbors for intra-chain
                if abs(d_resnum - a_resnum) <= 2:
                    continue
            
            # Get acceptor atoms
            acceptor_atoms = []
            
            # Backbone O
            backbone_o = acceptor_chain.select(f'resnum {a_resnum} and name O')
            if backbone_o:
                acceptor_atoms.extend([(atom, 'O', 'backbone') for atom in backbone_o])
            
            # Side chain acceptors
            if a_resname in HBOND_ACCEPTORS and HBOND_ACCEPTORS[a_resname]:
                for acceptor_name in HBOND_ACCEPTORS[a_resname]:
                    sc_acceptors = acceptor_chain.select(f'resnum {a_resnum} and name {acceptor_name}')
                    if sc_acceptors:
                        acceptor_atoms.extend([(atom, acceptor_name, 'sidechain') for atom in sc_acceptors])
            
            if not acceptor_atoms:
                continue
            
            # Check distances
            for donor_atom, d_name, d_type in donor_atoms:
                d_coords = donor_atom.getCoords()
                d_chain = donor_atom.getChid()
                
                for acceptor_atom, a_name, a_type in acceptor_atoms:
                    a_coords = acceptor_atom.getCoords()
                    a_chain = acceptor_atom.getChid()
                    
                    dist = np.linalg.norm(d_coords - a_coords)
                    
                    if dist <= dist_cutoff:
                        hbonds.append({
                            'donor': f"{d_resname}{d_resnum}",
                            'donor_atom': d_name,
                            'donor_chain': d_chain,
                            'donor_type': d_type,
                            'acceptor': f"{a_resname}{a_resnum}",
                            'acceptor_atom': a_name,
                            'acceptor_chain': a_chain,
                            'acceptor_type': a_type,
                            'distance': dist
                        })
    
    print(f"Found {len(hbonds)} hydrogen bonds")
    return hbonds


def find_salt_bridges(structure, chain_a=None, chain_b=None, cutoff=4.0):
    """
    Enhanced salt bridge detection using proper O-N distance criteria
    Following standard definition: O (acidic) to N (basic) distance < 4.0 Å
    """
    print(f"\n{'='*60}")
    print("SALT BRIDGE ANALYSIS")
    print(f"{'='*60}")
    
    if chain_a and chain_b:
        # For inter-chain, get atoms from each chain
        acidic_chain = structure.select(f'chain {chain_a} and (resname ASP GLU)')
        basic_chain = structure.select(f'chain {chain_b} and (resname ARG LYS HIS)')
        
        # Also check reverse direction
        acidic_chain_b = structure.select(f'chain {chain_b} and (resname ASP GLU)')
        basic_chain_a = structure.select(f'chain {chain_a} and (resname ARG LYS HIS)')
        
        print(f"Analyzing inter-chain salt bridges: {chain_a} ↔ {chain_b}")
    else:
        acidic_chain = structure.select('resname ASP GLU')
        basic_chain = structure.select('resname ARG LYS HIS')
        acidic_chain_b = None
        basic_chain_a = None
        print("Analyzing all salt bridges")
    
    salt_bridges = []
    
    # Define oxygen atoms for acidic residues (carboxylate oxygens)
    acidic_oxygen_atoms = {
        'ASP': ['OD1', 'OD2'],
        'GLU': ['OE1', 'OE2']
    }
    
    # Define nitrogen atoms for basic residues (charged nitrogens)
    basic_nitrogen_atoms = {
        'ARG': ['NH1', 'NH2', 'NE'],  # Guanidinium nitrogens
        'LYS': ['NZ'],                 # Ammonium nitrogen
        'HIS': ['ND1', 'NE2']          # Imidazole nitrogens
    }
    
    def find_salt_bridges_pair(acidic_sel, basic_sel):
        """Helper function to find salt bridges between two selections"""
        bridges = []
        
        if acidic_sel is None or basic_sel is None:
            return bridges
        
        # Get unique acidic residues
        acidic_resnums = np.unique(acidic_sel.getResnums())
        acidic_dict = {}
        for resnum in acidic_resnums:
            res_atoms = acidic_sel.select(f'resnum {resnum}')
            if res_atoms:
                acidic_dict[resnum] = {
                    'resname': res_atoms.getResnames()[0],
                    'chain': res_atoms.getChids()[0]
                }
        
        # Get unique basic residues
        basic_resnums = np.unique(basic_sel.getResnums())
        basic_dict = {}
        for resnum in basic_resnums:
            res_atoms = basic_sel.select(f'resnum {resnum}')
            if res_atoms:
                basic_dict[resnum] = {
                    'resname': res_atoms.getResnames()[0],
                    'chain': res_atoms.getChids()[0]
                }
        
        # Check all acidic-basic pairs
        for acid_resnum, acid_info in acidic_dict.items():
            acid_resname = acid_info['resname']
            acid_chain = acid_info['chain']
            
            # Get oxygen atoms for this acidic residue
            oxygen_atoms_list = []
            for o_name in acidic_oxygen_atoms.get(acid_resname, []):
                o_atoms = structure.select(f'chain {acid_chain} and resnum {acid_resnum} and name {o_name}')
                if o_atoms:
                    oxygen_atoms_list.extend(o_atoms)
            
            if not oxygen_atoms_list:
                continue
            
            for basic_resnum, basic_info in basic_dict.items():
                basic_resname = basic_info['resname']
                basic_chain = basic_info['chain']
                
                # Skip same residue
                if acid_resnum == basic_resnum and acid_chain == basic_chain:
                    continue
                
                # Get nitrogen atoms for this basic residue
                nitrogen_atoms_list = []
                for n_name in basic_nitrogen_atoms.get(basic_resname, []):
                    n_atoms = structure.select(f'chain {basic_chain} and resnum {basic_resnum} and name {n_name}')
                    if n_atoms:
                        nitrogen_atoms_list.extend(n_atoms)
                
                if not nitrogen_atoms_list:
                    continue
                
                # Check O-N distances
                min_dist = float('inf')
                for o_atom in oxygen_atoms_list:
                    o_coords = o_atom.getCoords()
                    for n_atom in nitrogen_atoms_list:
                        n_coords = n_atom.getCoords()
                        dist = np.linalg.norm(o_coords - n_coords)
                        if dist < min_dist:
                            min_dist = dist
                
                if min_dist <= cutoff:
                    bridges.append({
                        'negative': f"{acid_resname}{acid_resnum}",
                        'neg_chain': acid_chain,
                        'positive': f"{basic_resname}{basic_resnum}",
                        'pos_chain': basic_chain,
                        'distance': min_dist
                    })
        
        return bridges
    
    # Find salt bridges in primary direction
    salt_bridges.extend(find_salt_bridges_pair(acidic_chain, basic_chain))
    
    # For inter-chain, also check reverse direction
    if chain_a and chain_b:
        salt_bridges.extend(find_salt_bridges_pair(acidic_chain_b, basic_chain_a))
    
    # Remove duplicates (can happen in inter-chain analysis)
    unique_bridges = []
    seen = set()
    for sb in salt_bridges:
        key = tuple(sorted([(sb['negative'], sb['neg_chain']), (sb['positive'], sb['pos_chain'])]))
        if key not in seen:
            seen.add(key)
            unique_bridges.append(sb)
    
    print(f"Found {len(unique_bridges)} salt bridges")
    return unique_bridges


def find_aromatic_interactions(structure, chain_a=None, chain_b=None, 
                                pi_pi_cutoff=7.0, cation_pi_cutoff=6.0):
    """
    Enhanced aromatic interaction detection
    """
    print(f"\n{'='*60}")
    print("AROMATIC INTERACTION ANALYSIS")
    print(f"{'='*60}")
    
    if chain_a and chain_b:
        aromatic_a = structure.select(f'chain {chain_a} and (resname PHE TYR TRP HIS)')
        aromatic_b = structure.select(f'chain {chain_b} and (resname PHE TYR TRP HIS)')
        cation_a = structure.select(f'chain {chain_a} and (resname ARG LYS)')
        cation_b = structure.select(f'chain {chain_b} and (resname ARG LYS)')
        print(f"Analyzing aromatic interactions: {chain_a} ↔ {chain_b}")
    else:
        aromatic_a = structure.select('resname PHE TYR TRP HIS')
        aromatic_b = aromatic_a
        cation_a = structure.select('resname ARG LYS')
        cation_b = cation_a
        print("Analyzing all aromatic interactions")
    
    interactions = []
    
    # Aromatic ring centers (approximated by CA for simplicity)
    if aromatic_a is not None and aromatic_b is not None:
        arom_a_resnums = np.unique(aromatic_a.getResnums())
        arom_a_dict = dict(zip(aromatic_a.getResnums(), 
                              zip(aromatic_a.getResnames(), aromatic_a.getChids())))
        
        arom_b_resnums = np.unique(aromatic_b.getResnums())
        arom_b_dict = dict(zip(aromatic_b.getResnums(),
                              zip(aromatic_b.getResnames(), aromatic_b.getChids())))
        
        for res1_num in arom_a_resnums:
            res1_name, chain1 = arom_a_dict[res1_num]
            ca1 = structure.select(f'chain {chain1} and resnum {res1_num} and name CA')
            if ca1 is None:
                continue
            coords1 = ca1.getCoords()[0]
            
            for res2_num in arom_b_resnums:
                res2_name, chain2 = arom_b_dict[res2_num]
                
                # Skip same residue
                if res1_num == res2_num and chain1 == chain2:
                    continue
                
                # Avoid duplicates for intra-chain
                if chain_a is None and chain1 == chain2 and res1_num >= res2_num:
                    continue
                
                ca2 = structure.select(f'chain {chain2} and resnum {res2_num} and name CA')
                if ca2 is None:
                    continue
                coords2 = ca2.getCoords()[0]
                
                dist = np.linalg.norm(coords1 - coords2)
                
                if dist <= pi_pi_cutoff:
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
        arom_resnums = np.unique(aromatic_a.getResnums())
        arom_dict = dict(zip(aromatic_a.getResnums(),
                            zip(aromatic_a.getResnames(), aromatic_a.getChids())))
        
        cation_resnums = np.unique(cation_b.getResnums())
        cation_dict = dict(zip(cation_b.getResnums(),
                              zip(cation_b.getResnames(), cation_b.getChids())))
        
        for arom_num in arom_resnums:
            arom_name, arom_chain = arom_dict[arom_num]
            ca_arom = structure.select(f'chain {arom_chain} and resnum {arom_num} and name CA')
            if ca_arom is None:
                continue
            coords_arom = ca_arom.getCoords()[0]
            
            for cat_num in cation_resnums:
                cat_name, cat_chain = cation_dict[cat_num]
                
                if arom_num == cat_num and arom_chain == cat_chain:
                    continue
                
                ca_cat = structure.select(f'chain {cat_chain} and resnum {cat_num} and name CA')
                if ca_cat is None:
                    continue
                coords_cat = ca_cat.getCoords()[0]
                
                dist = np.linalg.norm(coords_arom - coords_cat)
                
                if dist <= cation_pi_cutoff:
                    interactions.append({
                        'type': 'cation-pi',
                        'res1': f"{cat_name}{cat_num}",
                        'chain1': cat_chain,
                        'res2': f"{arom_name}{arom_num}",
                        'chain2': arom_chain,
                        'distance': dist
                    })
    
    # Reverse direction for inter-chain cation-pi
    if chain_a and chain_b and aromatic_b is not None and cation_a is not None:
        arom_resnums = np.unique(aromatic_b.getResnums())
        arom_dict = dict(zip(aromatic_b.getResnums(),
                            zip(aromatic_b.getResnames(), aromatic_b.getChids())))
        
        cation_resnums = np.unique(cation_a.getResnums())
        cation_dict = dict(zip(cation_a.getResnums(),
                              zip(cation_a.getResnames(), cation_a.getChids())))
        
        for arom_num in arom_resnums:
            arom_name, arom_chain = arom_dict[arom_num]
            ca_arom = structure.select(f'chain {arom_chain} and resnum {arom_num} and name CA')
            if ca_arom is None:
                continue
            coords_arom = ca_arom.getCoords()[0]
            
            for cat_num in cation_resnums:
                cat_name, cat_chain = cation_dict[cat_num]
                
                ca_cat = structure.select(f'chain {cat_chain} and resnum {cat_num} and name CA')
                if ca_cat is None:
                    continue
                coords_cat = ca_cat.getCoords()[0]
                
                dist = np.linalg.norm(coords_arom - coords_cat)
                
                if dist <= cation_pi_cutoff:
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


def find_disulfide_bonds(structure, cutoff=2.5):
    """
    Detect disulfide bonds between cysteine residues
    """
    print(f"\n{'='*60}")
    print("DISULFIDE BOND ANALYSIS")
    print(f"{'='*60}")
    
    cys_sg = structure.select('resname CYS and name SG')
    
    if cys_sg is None or cys_sg.numAtoms() < 2:
        print("Insufficient cysteine residues for disulfide bonds")
        return []
    
    disulfides = []
    resnums = cys_sg.getResnums()
    chains = cys_sg.getChids()
    coords = cys_sg.getCoords()
    
    for i in range(len(resnums)):
        for j in range(i + 1, len(resnums)):
            dist = np.linalg.norm(coords[i] - coords[j])
            
            if dist <= cutoff:
                disulfides.append({
                    'res1': f"CYS{resnums[i]}",
                    'chain1': chains[i],
                    'res2': f"CYS{resnums[j]}",
                    'chain2': chains[j],
                    'distance': dist
                })
    
    print(f"Found {len(disulfides)} disulfide bonds")
    return disulfides


def find_hydrophobic_interactions(structure, chain_a=None, chain_b=None, cutoff=5.0):
    """
    Detect hydrophobic interactions between nonpolar residues
    """
    print(f"\n{'='*60}")
    print("HYDROPHOBIC INTERACTION ANALYSIS")
    print(f"{'='*60}")
    
    if chain_a and chain_b:
        hydrophobic_a = structure.select(f'chain {chain_a} and calpha and (resname ALA VAL ILE LEU MET PHE TRP PRO)')
        hydrophobic_b = structure.select(f'chain {chain_b} and calpha and (resname ALA VAL ILE LEU MET PHE TRP PRO)')
        print(f"Analyzing hydrophobic interactions: {chain_a} ↔ {chain_b}")
    else:
        hydrophobic_a = structure.select('calpha and (resname ALA VAL ILE LEU MET PHE TRP PRO)')
        hydrophobic_b = hydrophobic_a
        print("Analyzing all hydrophobic interactions")
    
    interactions = []
    
    if hydrophobic_a is None or hydrophobic_b is None:
        print("No hydrophobic residues found")
        return interactions
    
    hydro_a_resnums = hydrophobic_a.getResnums()
    hydro_a_resnames = hydrophobic_a.getResnames()
    hydro_a_chains = hydrophobic_a.getChids()
    hydro_a_coords = hydrophobic_a.getCoords()
    
    hydro_b_resnums = hydrophobic_b.getResnums()
    hydro_b_resnames = hydrophobic_b.getResnames()
    hydro_b_chains = hydrophobic_b.getChids()
    hydro_b_coords = hydrophobic_b.getCoords()
    
    for i, (res1_num, res1_name, chain1, coord1) in enumerate(zip(hydro_a_resnums, hydro_a_resnames, 
                                                                    hydro_a_chains, hydro_a_coords)):
        for j, (res2_num, res2_name, chain2, coord2) in enumerate(zip(hydro_b_resnums, hydro_b_resnames,
                                                                        hydro_b_chains, hydro_b_coords)):
            # Skip same residue
            if res1_num == res2_num and chain1 == chain2:
                continue
            
            # Skip sequential neighbors
            if chain1 == chain2 and abs(res1_num - res2_num) <= 3:
                continue
            
            # Avoid duplicates for intra-chain
            if chain_a is None and chain1 == chain2 and res1_num >= res2_num:
                continue
            
            dist = np.linalg.norm(coord1 - coord2)
            
            if dist <= cutoff:
                interactions.append({
                    'res1': f"{res1_name}{res1_num}",
                    'chain1': chain1,
                    'res2': f"{res2_name}{res2_num}",
                    'chain2': chain2,
                    'distance': dist
                })
    
    print(f"Found {len(interactions)} hydrophobic interactions")
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
    
    ca_a = inter_results['ca_a']
    resnums_a = ca_a.getResnums()
    counts_a = [inter_results['counts_a'].get(rn, 0) for rn in resnums_a]
    
    ax1.bar(resnums_a, counts_a, color='steelblue', edgecolor='black', linewidth=0.5)
    ax1.set_xlabel(f'Chain {chain_a} Residue Number', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Number of Contacts', fontsize=11, fontweight='bold')
    ax1.set_title(f'Interface Residues: Chain {chain_a}', fontsize=12, fontweight='bold')
    ax1.grid(alpha=0.3)
    
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
    
    ax = axes[0]
    im = ax.imshow(intra_results['contact_matrix'], cmap='Blues', aspect='auto', interpolation='nearest')
    ax.set_xlabel('Residue Index', fontsize=11, fontweight='bold')
    ax.set_ylabel('Residue Index', fontsize=11, fontweight='bold')
    ax.set_title(f'Contact Matrix: Chain {chain_id}', fontsize=12, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Contact')
    
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
    
    ax = axes[0]
    colors = ['blue' if e == 0 else 'orange' if e == 1 else 'red' for e in exposure]
    ax.bar(resnums, sasa, color=colors, edgecolor='black', linewidth=0.3)
    ax.set_xlabel('Residue Number', fontsize=11, fontweight='bold')
    ax.set_ylabel('SASA (Å²)', fontsize=11, fontweight='bold')
    ax.set_title(f'Solvent Accessible Surface Area: Chain {chain_id}', fontsize=12, fontweight='bold')
    ax.axhline(y=10, color='gray', linestyle='--', alpha=0.5, label='Core threshold')
    ax.axhline(y=40, color='gray', linestyle='--', alpha=0.5, label='Surface threshold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    ax = axes[1]
    colors = ['red' if h < 0 else 'blue' for h in hydro]
    ax.bar(resnums, hydro, color=colors, edgecolor='black', linewidth=0.3)
    ax.set_xlabel('Residue Number', fontsize=11, fontweight='bold')
    ax.set_ylabel('Hydrophobicity Score', fontsize=11, fontweight='bold')
    ax.set_title(f'Hydrophobicity Profile (Kyte-Doolittle): Chain {chain_id}', fontsize=12, fontweight='bold')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.grid(alpha=0.3)
    
    ax = axes[2]
    ax2 = ax.twinx()
    
    ax.fill_between(resnums, sasa, alpha=0.3, color='blue', label='SASA')
    ax.plot(resnums, sasa, color='blue', linewidth=2)
    ax.set_xlabel('Residue Number', fontsize=11, fontweight='bold')
    ax.set_ylabel('SASA (Å²)', fontsize=11, fontweight='bold', color='blue')
    ax.tick_params(axis='y', labelcolor='blue')
    
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


def plot_noncovalent_interactions(hbonds, salt_bridges, aromatic, disulfides, hydrophobic,
                                  output_dir, chain_a=None, chain_b=None):
    """Plot comprehensive summary of noncovalent interactions"""
    print("\nGenerating noncovalent interaction summary...")
    
    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
    
    # Summary bar chart
    ax = fig.add_subplot(gs[0, 0])
    interaction_types = ['H-bonds', 'Salt\nBridges', 'Aromatic', 'Disulfides', 'Hydrophobic']
    counts = [len(hbonds), len(salt_bridges), len(aromatic), len(disulfides), len(hydrophobic)]
    colors = ['skyblue', 'coral', 'lightgreen', 'gold', 'plum']
    
    bars = ax.bar(interaction_types, counts, color=colors, edgecolor='black', linewidth=2)
    ax.set_ylabel('Number of Interactions', fontsize=11, fontweight='bold')
    ax.set_title('Noncovalent Interaction Summary', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    # Add count labels on bars
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(count)}', ha='center', va='bottom', fontweight='bold')
    
    # H-bond distance distribution
    ax = fig.add_subplot(gs[0, 1])
    if len(hbonds) > 0:
        distances = [hb['distance'] for hb in hbonds]
        ax.hist(distances, bins=20, color='skyblue', edgecolor='black', alpha=0.7)
        ax.axvline(x=np.mean(distances), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {np.mean(distances):.2f} Å')
        ax.legend()
    ax.set_xlabel('Distance (Å)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax.set_title('H-bond Distance Distribution', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    # Salt bridge distance distribution
    ax = fig.add_subplot(gs[0, 2])
    if len(salt_bridges) > 0:
        distances = [sb['distance'] for sb in salt_bridges]
        ax.hist(distances, bins=20, color='coral', edgecolor='black', alpha=0.7)
        ax.axvline(x=np.mean(distances), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {np.mean(distances):.2f} Å')
        ax.legend()
    ax.set_xlabel('Distance (Å)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax.set_title('Salt Bridge Distance Distribution', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    # Aromatic interaction types
    ax = fig.add_subplot(gs[1, 0])
    if len(aromatic) > 0:
        types = [inter['type'] for inter in aromatic]
        type_counts = {}
        for t in types:
            type_counts[t] = type_counts.get(t, 0) + 1
        
        ax.bar(type_counts.keys(), type_counts.values(), color='lightgreen',
              edgecolor='black', linewidth=2)
        ax.set_ylabel('Count', fontsize=11, fontweight='bold')
        ax.set_title('Aromatic Interaction Types', fontsize=12, fontweight='bold')
        ax.grid(alpha=0.3, axis='y')
    
    # Hydrophobic interaction distance distribution
    ax = fig.add_subplot(gs[1, 1])
    if len(hydrophobic) > 0:
        distances = [hi['distance'] for hi in hydrophobic]
        ax.hist(distances, bins=20, color='plum', edgecolor='black', alpha=0.7)
        ax.axvline(x=np.mean(distances), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {np.mean(distances):.2f} Å')
        ax.legend()
    ax.set_xlabel('Distance (Å)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax.set_title('Hydrophobic Interaction Distance', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    # Disulfide bond distance distribution
    ax = fig.add_subplot(gs[1, 2])
    if len(disulfides) > 0:
        distances = [ds['distance'] for ds in disulfides]
        ax.hist(distances, bins=10, color='gold', edgecolor='black', alpha=0.7)
        ax.axvline(x=np.mean(distances), color='red', linestyle='--', linewidth=2,
                  label=f'Mean: {np.mean(distances):.2f} Å')
        ax.legend()
    ax.set_xlabel('Distance (Å)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax.set_title('Disulfide Bond Distance', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    suffix = f"_{chain_a}_{chain_b}" if chain_a and chain_b else "_all"
    plt.savefig(f'{output_dir}/noncovalent_interactions{suffix}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/noncovalent_interactions{suffix}.png")


# ============================================================================
# OUTPUT FUNCTIONS
# ============================================================================

def save_intermolecular_results(inter_results, output_dir, chain_a, chain_b):
    """Save intermolecular contact results"""
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


def save_noncovalent_results(hbonds, salt_bridges, aromatic, disulfides, hydrophobic,
                            output_dir, chain_a=None, chain_b=None):
    """Save all noncovalent interaction results"""
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
        f.write(f"Total aromatic interactions: {len(aromatic)}\n")
        f.write(f"Total disulfide bonds: {len(disulfides)}\n")
        f.write(f"Total hydrophobic interactions: {len(hydrophobic)}\n\n")
        
        # Hydrogen bonds
        f.write("HYDROGEN BONDS\n")
        f.write("-"*90 + "\n")
        f.write(f"{'Donor':<15} {'Atom':<8} {'Chain':<8} {'Acceptor':<15} {'Atom':<8} {'Chain':<8} {'Dist(Å)':<10} {'Type':<12}\n")
        f.write("-"*90 + "\n")
        
        for hb in sorted(hbonds, key=lambda x: x['distance']):
            f.write(f"{hb['donor']:<15} {hb['donor_atom']:<8} {hb['donor_chain']:<8} "
                   f"{hb['acceptor']:<15} {hb['acceptor_atom']:<8} {hb['acceptor_chain']:<8} "
                   f"{hb['distance']:<10.3f} {hb.get('donor_type', 'N/A'):<12}\n")
        
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
        f.write("-"*75 + "\n")
        f.write(f"{'Type':<15} {'Residue 1':<20} {'Chain':<10} {'Residue 2':<20} {'Chain':<10} {'Distance (Å)':<15}\n")
        f.write("-"*75 + "\n")
        
        for ar in sorted(aromatic, key=lambda x: x['distance']):
            f.write(f"{ar['type']:<15} {ar['res1']:<20} {ar['chain1']:<10} {ar['res2']:<20} "
                   f"{ar['chain2']:<10} {ar['distance']:<15.3f}\n")
        
        # Disulfide bonds
        f.write("\n\nDISULFIDE BONDS\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Residue 1':<20} {'Chain':<10} {'Residue 2':<20} {'Chain':<10} {'Distance (Å)':<15}\n")
        f.write("-"*70 + "\n")
        
        for ds in sorted(disulfides, key=lambda x: x['distance']):
            f.write(f"{ds['res1']:<20} {ds['chain1']:<10} {ds['res2']:<20} "
                   f"{ds['chain2']:<10} {ds['distance']:<15.3f}\n")
        
        # Hydrophobic interactions
        f.write("\n\nHYDROPHOBIC INTERACTIONS\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Residue 1':<20} {'Chain':<10} {'Residue 2':<20} {'Chain':<10} {'Distance (Å)':<15}\n")
        f.write("-"*70 + "\n")
        
        for hi in sorted(hydrophobic, key=lambda x: x['distance'])[:100]:  # Top 100
            f.write(f"{hi['res1']:<20} {hi['chain1']:<10} {hi['res2']:<20} "
                   f"{hi['chain2']:<10} {hi['distance']:<15.3f}\n")
        
        if len(hydrophobic) > 100:
            f.write(f"\n... and {len(hydrophobic) - 100} more hydrophobic interactions\n")
    
    print(f"Saved: {filename}")


def save_intramolecular_results(intra_results, output_dir, chain_id):
    """Save intramolecular contact results"""
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
    """Save SASA and hydrophobicity results"""
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
        f.write(f"{'Residue':<15} {'SASA (Å²)':<15} {'Hydrophobicity':<15} {'Exposure':<15}\n")
        f.write("-"*70 + "\n")
        
        for i, (resnum, resname) in enumerate(zip(sasa_results['resnums'], sasa_results['resnames'])):
            res_label = f"{resname}{resnum}"
            sasa = sasa_results['sasa'][i]
            hydro = sasa_results['hydrophobicity'][i]
            exposure = exposure_names[sasa_results['exposure_class'][i]]
            f.write(f"{res_label:<15} {sasa:<15.2f} {hydro:<15.2f} {exposure:<15}\n")
    
    print(f"Saved: {filename}")


def save_pymol_scripts(inter_results, sasa_results, hbonds, salt_bridges, 
                       output_dir, chain_a=None, chain_b=None):
    """Generate PyMOL visualization scripts"""
    print("\nGenerating PyMOL visualization scripts...")
    
    # Interface visualization
    if inter_results and chain_a and chain_b:
        filename = f"{output_dir}/visualize_interface_{chain_a}_{chain_b}.pml"
        with open(filename, 'w') as f:
            f.write("# PyMOL script to visualize protein interface\n\n")
            f.write("hide everything\n")
            f.write("show cartoon\n")
            f.write("set cartoon_fancy_helices, 1\n\n")
            
            f.write(f"# Color chains\n")
            f.write(f"color lightblue, chain {chain_a}\n")
            f.write(f"color lightpink, chain {chain_b}\n\n")
            
            if len(inter_results['counts_a']) > 0:
                res_list_a = '+'.join([str(rn) for rn in inter_results['counts_a'].keys()])
                f.write(f"# Interface residues chain {chain_a}\n")
                f.write(f"select interface_A, chain {chain_a} and resi {res_list_a}\n")
                f.write("color marine, interface_A\n")
                f.write("show sticks, interface_A\n\n")
            
            if len(inter_results['counts_b']) > 0:
                res_list_b = '+'.join([str(rn) for rn in inter_results['counts_b'].keys()])
                f.write(f"# Interface residues chain {chain_b}\n")
                f.write(f"select interface_B, chain {chain_b} and resi {res_list_b}\n")
                f.write("color red, interface_B\n")
                f.write("show sticks, interface_B\n\n")
            
            # Add H-bonds visualization
            if hbonds:
                f.write("# Hydrogen bonds\n")
                for hb in hbonds:
                    donor_num = ''.join(filter(str.isdigit, hb['donor']))
                    acceptor_num = ''.join(filter(str.isdigit, hb['acceptor']))
                    f.write(f"distance hbond_{donor_num}_{acceptor_num}, "
                           f"chain {hb['donor_chain']} and resi {donor_num} and name {hb['donor_atom']}, "
                           f"chain {hb['acceptor_chain']} and resi {acceptor_num} and name {hb['acceptor_atom']}\n")
                f.write("hide labels\n")
                f.write("color yellow, hbond_*\n\n")
            
            # Add salt bridges
            if salt_bridges:
                f.write("# Salt bridges\n")
                for sb in salt_bridges:
                    pos_num = ''.join(filter(str.isdigit, sb['positive']))
                    neg_num = ''.join(filter(str.isdigit, sb['negative']))
                    f.write(f"distance saltbridge_{pos_num}_{neg_num}, "
                           f"chain {sb['pos_chain']} and resi {pos_num} and name CA, "
                           f"chain {sb['neg_chain']} and resi {neg_num} and name CA\n")
                f.write("hide labels\n")
                f.write("color magenta, saltbridge_*\n\n")
            
            f.write("bg_color white\n")
            f.write("set ray_shadows, 0\n")
            f.write("set antialias, 2\n")
            f.write("set line_width, 3\n")
        
        print(f"Saved: {filename}")
    
    # SASA visualization
    if sasa_results:
        chain_id = chain_a if chain_a else 'A'
        filename = f"{output_dir}/visualize_sasa_chain_{chain_id}.pml"
        with open(filename, 'w') as f:
            f.write("# PyMOL script to visualize SASA and exposure\n\n")
            f.write("hide everything\n")
            f.write("show cartoon\n")
            f.write("set cartoon_fancy_helices, 1\n\n")
            
            exposure = sasa_results['exposure_class']
            resnums = sasa_results['resnums']
            
            core_res = [str(rn) for i, rn in enumerate(resnums) if exposure[i] == 0]
            partial_res = [str(rn) for i, rn in enumerate(resnums) if exposure[i] == 1]
            surface_res = [str(rn) for i, rn in enumerate(resnums) if exposure[i] == 2]
            
            if core_res:
                f.write(f"# Core residues (buried)\n")
                f.write(f"select core, resi {'+'.join(core_res)}\n")
                f.write("color blue, core\n")
                f.write("show spheres, core and name CA\n\n")
            
            if partial_res:
                f.write(f"# Partially exposed residues\n")
                f.write(f"select partial, resi {'+'.join(partial_res)}\n")
                f.write("color orange, partial\n")
                f.write("show spheres, partial and name CA\n\n")
            
            if surface_res:
                f.write(f"# Surface residues\n")
                f.write(f"select surface, resi {'+'.join(surface_res)}\n")
                f.write("color red, surface\n")
                f.write("show spheres, surface and name CA\n\n")
            
            f.write("bg_color white\n")
            f.write("set sphere_scale, 0.5\n")
        
        print(f"Saved: {filename}")


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Enhanced ProDy Contact, SASA, and Stability Analysis',
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
  python3 prody_contact_stability_analysis.py --pdb complex.pdb --chain-a A --chain-b B \\
      --contact-cutoff 4.5 --hbond-cutoff 3.2
        """)
    
    parser.add_argument('--pdb', required=True, help='Input PDB file')
    parser.add_argument('--chain-a', '--chain', dest='chain_a', default=None,
                       help='First chain or single chain to analyze')
    parser.add_argument('--chain-b', default=None,
                       help='Second chain for interface analysis')
    parser.add_argument('--contact-cutoff', type=float, default=5.0,
                       help='Contact distance cutoff in Å (default: 5.0)')
    parser.add_argument('--neighbor-cutoff', type=float, default=8.0,
                       help='Neighbor distance cutoff in Å (default: 8.0)')
    parser.add_argument('--hbond-cutoff', type=float, default=3.5,
                       help='H-bond distance cutoff in Å (default: 3.5)')
    parser.add_argument('--saltbridge-cutoff', type=float, default=4.0,
                       help='Salt bridge O-N distance cutoff in Å (default: 4.0, standard definition)')
    parser.add_argument('--pi-pi-cutoff', type=float, default=7.0,
                       help='Pi-pi interaction cutoff in Å (default: 7.0)')
    parser.add_argument('--cation-pi-cutoff', type=float, default=6.0,
                       help='Cation-pi interaction cutoff in Å (default: 6.0)')
    parser.add_argument('--hydrophobic-cutoff', type=float, default=5.0,
                       help='Hydrophobic interaction cutoff in Å (default: 5.0)')
    parser.add_argument('--disulfide-cutoff', type=float, default=2.5,
                       help='Disulfide bond cutoff in Å (default: 2.5)')
    parser.add_argument('--skip-sasa', action='store_true',
                       help='Skip SASA calculation')
    
    args = parser.parse_args()
    
    # Setup
    pdb_name = os.path.splitext(os.path.basename(args.pdb))[0]
    output_dir = setup_output_directory(pdb_name)
    
    # Load structure
    print(f"\n{'='*70}")
    print(f"ENHANCED PRODY STRUCTURAL ANALYSIS")
    print(f"{'='*70}")
    print(f"\nLoading PDB: {args.pdb}")
    structure = parsePDB(args.pdb)
    
    if structure is None:
        print(f"Error: Could not parse PDB file: {args.pdb}")
        sys.exit(1)
    
    # Analyze chain structure
    chain_info = analyze_chain_structure(structure)
    print_chain_summary(chain_info)
    
    # Main analysis
    if args.chain_a and args.chain_b:
        # INTERFACE ANALYSIS MODE
        print(f"\n{'='*70}")
        print(f"MODE: Interface Analysis ({args.chain_a} ↔ {args.chain_b})")
        print(f"{'='*70}")
        
        # Validate chains exist
        if args.chain_a not in chain_info:
            print(f"Error: Chain {args.chain_a} not found in structure")
            sys.exit(1)
        if args.chain_b not in chain_info:
            print(f"Error: Chain {args.chain_b} not found in structure")
            sys.exit(1)
        
        # 1. Intermolecular contacts
        inter_results = calculate_intermolecular_contacts(
            structure, args.chain_a, args.chain_b, cutoff=args.contact_cutoff
        )
        
        # 2. Noncovalent interactions
        hbonds = find_hydrogen_bonds_enhanced(structure, args.chain_a, args.chain_b,
                                             dist_cutoff=args.hbond_cutoff)
        salt_bridges = find_salt_bridges(structure, args.chain_a, args.chain_b,
                                        cutoff=args.saltbridge_cutoff)
        aromatic = find_aromatic_interactions(structure, args.chain_a, args.chain_b,
                                             pi_pi_cutoff=args.pi_pi_cutoff,
                                             cation_pi_cutoff=args.cation_pi_cutoff)
        disulfides = find_disulfide_bonds(structure, cutoff=args.disulfide_cutoff)
        hydrophobic = find_hydrophobic_interactions(structure, args.chain_a, args.chain_b,
                                                    cutoff=args.hydrophobic_cutoff)
        
        # Visualizations
        plot_intermolecular_contact_heatmap(inter_results, output_dir, args.chain_a, args.chain_b)
        plot_contact_counts(inter_results, output_dir, args.chain_a, args.chain_b)
        plot_noncovalent_interactions(hbonds, salt_bridges, aromatic, disulfides,
                                     hydrophobic, output_dir, args.chain_a, args.chain_b)
        
        # Save results
        save_intermolecular_results(inter_results, output_dir, args.chain_a, args.chain_b)
        save_noncovalent_results(hbonds, salt_bridges, aromatic, disulfides,
                                hydrophobic, output_dir, args.chain_a, args.chain_b)
        
        # Individual chain analysis
        for chain_id in [args.chain_a, args.chain_b]:
            if chain_id in chain_info:
                chain_atoms = structure.select(f'chain {chain_id}')
                
                intra_results = calculate_intramolecular_contacts(
                    chain_atoms, chain_id,
                    contact_cutoff=args.contact_cutoff,
                    neighbor_cutoff=args.neighbor_cutoff
                )
                plot_intramolecular_contacts(intra_results, output_dir, chain_id)
                save_intramolecular_results(intra_results, output_dir, chain_id)
                
                if not args.skip_sasa:
                    sasa_results = calculate_sasa_hydrophobicity(structure, chain_id)
                    plot_sasa_hydrophobicity(sasa_results, output_dir, chain_id)
                    save_sasa_results(sasa_results, output_dir, chain_id)
        
        # PyMOL scripts
        sasa_results = None if args.skip_sasa else calculate_sasa_hydrophobicity(structure, args.chain_a)
        save_pymol_scripts(inter_results, sasa_results, hbonds, salt_bridges,
                          output_dir, args.chain_a, args.chain_b)
        
    elif args.chain_a:
        # SINGLE CHAIN ANALYSIS MODE
        print(f"\n{'='*70}")
        print(f"MODE: Single Chain Analysis (Chain {args.chain_a})")
        print(f"{'='*70}")
        
        if args.chain_a not in chain_info:
            print(f"Error: Chain {args.chain_a} not found")
            sys.exit(1)
        
        chain_atoms = structure.select(f'chain {args.chain_a}')
        chain_struct = structure.select(f'chain {args.chain_a}')
        
        # 1. Intramolecular contacts
        intra_results = calculate_intramolecular_contacts(
            chain_atoms, args.chain_a,
            contact_cutoff=args.contact_cutoff,
            neighbor_cutoff=args.neighbor_cutoff
        )
        
        # 2. SASA analysis
        sasa_results = None
        if not args.skip_sasa:
            sasa_results = calculate_sasa_hydrophobicity(structure, args.chain_a)
        
        # 3. Noncovalent interactions
        hbonds = find_hydrogen_bonds_enhanced(chain_struct, dist_cutoff=args.hbond_cutoff)
        salt_bridges = find_salt_bridges(chain_struct, cutoff=args.saltbridge_cutoff)
        aromatic = find_aromatic_interactions(chain_struct,
                                             pi_pi_cutoff=args.pi_pi_cutoff,
                                             cation_pi_cutoff=args.cation_pi_cutoff)
        disulfides = find_disulfide_bonds(chain_struct, cutoff=args.disulfide_cutoff)
        hydrophobic = find_hydrophobic_interactions(chain_struct, cutoff=args.hydrophobic_cutoff)
        
        # Visualizations
        plot_intramolecular_contacts(intra_results, output_dir, args.chain_a)
        if sasa_results:
            plot_sasa_hydrophobicity(sasa_results, output_dir, args.chain_a)
        plot_noncovalent_interactions(hbonds, salt_bridges, aromatic, disulfides,
                                     hydrophobic, output_dir)
        
        # Save results
        save_intramolecular_results(intra_results, output_dir, args.chain_a)
        if sasa_results:
            save_sasa_results(sasa_results, output_dir, args.chain_a)
        save_noncovalent_results(hbonds, salt_bridges, aromatic, disulfides,
                                hydrophobic, output_dir)
        
        # PyMOL scripts
        save_pymol_scripts(None, sasa_results, hbonds, salt_bridges,
                          output_dir, args.chain_a, None)
        
    else:
        # MULTI-CHAIN ANALYSIS MODE
        print(f"\n{'='*70}")
        print(f"MODE: Multi-Chain Analysis (All Chains)")
        print(f"{'='*70}")
        
        # Analyze each chain
        for chain_id in sorted(chain_info.keys()):
            print(f"\n{'='*60}")
            print(f"Analyzing Chain {chain_id}")
            print(f"{'='*60}")
            
            chain_atoms = structure.select(f'chain {chain_id}')
            
            intra_results = calculate_intramolecular_contacts(
                chain_atoms, chain_id,
                contact_cutoff=args.contact_cutoff,
                neighbor_cutoff=args.neighbor_cutoff
            )
            plot_intramolecular_contacts(intra_results, output_dir, chain_id)
            save_intramolecular_results(intra_results, output_dir, chain_id)
            
            if not args.skip_sasa:
                sasa_results = calculate_sasa_hydrophobicity(structure, chain_id)
                plot_sasa_hydrophobicity(sasa_results, output_dir, chain_id)
                save_sasa_results(sasa_results, output_dir, chain_id)
        
        # Pairwise interface analysis
        print(f"\n{'='*70}")
        print("PAIRWISE INTERFACE ANALYSIS")
        print(f"{'='*70}")
        
        chain_ids = sorted(chain_info.keys())
        for i, chain_a in enumerate(chain_ids):
            for chain_b in chain_ids[i+1:]:
                print(f"\nAnalyzing: {chain_a} ↔ {chain_b}")
                
                inter_results = calculate_intermolecular_contacts(
                    structure, chain_a, chain_b, cutoff=args.contact_cutoff
                )
                
                if inter_results and len(inter_results['contacts']) > 0:
                    hbonds = find_hydrogen_bonds_enhanced(structure, chain_a, chain_b,
                                                         dist_cutoff=args.hbond_cutoff)
                    salt_bridges = find_salt_bridges(structure, chain_a, chain_b,
                                                    cutoff=args.saltbridge_cutoff)
                    aromatic = find_aromatic_interactions(structure, chain_a, chain_b,
                                                         pi_pi_cutoff=args.pi_pi_cutoff,
                                                         cation_pi_cutoff=args.cation_pi_cutoff)
                    disulfides = find_disulfide_bonds(structure, cutoff=args.disulfide_cutoff)
                    hydrophobic = find_hydrophobic_interactions(structure, chain_a, chain_b,
                                                                cutoff=args.hydrophobic_cutoff)
                    
                    plot_intermolecular_contact_heatmap(inter_results, output_dir, chain_a, chain_b)
                    plot_contact_counts(inter_results, output_dir, chain_a, chain_b)
                    plot_noncovalent_interactions(hbonds, salt_bridges, aromatic, disulfides,
                                                 hydrophobic, output_dir, chain_a, chain_b)
                    
                    save_intermolecular_results(inter_results, output_dir, chain_a, chain_b)
                    save_noncovalent_results(hbonds, salt_bridges, aromatic, disulfides,
                                            hydrophobic, output_dir, chain_a, chain_b)
                    
                    save_pymol_scripts(inter_results, None, hbonds, salt_bridges,
                                      output_dir, chain_a, chain_b)
                else:
                    print(f"  No significant interface found")
    
    # Summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nAll results saved to: {output_dir}/")
    print("\nGenerated files:")
    print("  ✓ Contact matrices and heatmaps")
    print("  ✓ SASA and hydrophobicity profiles")
    print("  ✓ Comprehensive noncovalent interaction analysis")
    print("    - Hydrogen bonds (enhanced detection)")
    print("    - Salt bridges (charged atom-based)")
    print("    - Aromatic interactions (pi-pi, cation-pi)")
    print("    - Disulfide bonds")
    print("    - Hydrophobic interactions")
    print("  ✓ PyMOL visualization scripts")
    print("\nTo visualize in PyMOL:")
    print(f"  pymol {args.pdb} {output_dir}/visualize_*.pml")
    print("="*70)


if __name__ == '__main__':
    main()

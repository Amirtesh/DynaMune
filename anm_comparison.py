#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANM comparison with sequence-based chain alignment & mapping.

Computes:
  • Matched CA-only structures for apo and complex
  • ANM Mode Overlap (diagonal + best-match + full heatmap)
  • ΔRMSF (complex – apo)

Outputs:
  matched_apo_ca.pdb
  matched_complex_ca.pdb
  anm_overlap_matrix.txt
  anm_overlap_matrix.npy
  anm_overlap_diag.txt
  anm_overlap_max_per_mode.txt
  anm_overlap_bestmatch_idx.txt
  anm_overlap_summary.png
  anm_overlap_matrix_heatmap.png
  anm_overlap_summary.txt
  rmsf_apo.txt
  rmsf_complex.txt
  delta_rmsf.txt
  delta_rmsf_plot.png
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from prody import parsePDB, ANM, calcOverlap
from Bio.PDB.Polypeptide import three_to_one
from Bio import pairwise2


# -----------------------------------------------------------
# Utility
# -----------------------------------------------------------
def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def three_to_one_safe(res):
    try:
        return three_to_one(res.capitalize())
    except:
        return "X"


def atomgroup_ca_info(struct, chain):
    sel = struct.select(f"chain {chain} and calpha")
    if sel is None:
        return [], [], []
    return sel.getCoords(), sel.getResnames(), sel.getResnums()


def seq_from_resnames(res_list):
    return "".join([three_to_one_safe(r) for r in res_list])


# -----------------------------------------------------------
# Sequence alignment-based mapping
# -----------------------------------------------------------
def align_and_map(seqA, seqB):
    aln = pairwise2.align.globalxx(seqA, seqB, one_alignment_only=True)
    if len(aln) == 0:
        raise ValueError("Alignment failed")

    aA, aB, score, start, end = aln[0]

    idxA = 0
    idxB = 0
    matchesA = []
    matchesB = []

    for x, y in zip(aA, aB):
        if x != "-":
            curA = idxA
            idxA += 1
        else:
            curA = None

        if y != "-":
            curB = idxB
            idxB += 1
        else:
            curB = None

        if curA is not None and curB is not None:
            matchesA.append(curA)
            matchesB.append(curB)

    return matchesA, matchesB


# -----------------------------------------------------------
# Write matched CA-only PDB
# -----------------------------------------------------------
def write_ca_pdb(coords, resnames, outpath, chain="A"):
    with open(outpath, "w") as fh:
        serial = 1
        for i, (c, r) in enumerate(zip(coords, resnames)):
            fh.write(
                (
                    "ATOM  {serial:5d}  CA  {res:>3s} {chain:1s}{rnum:4d}    "
                    "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
                ).format(
                    serial=serial,
                    res=r,
                    chain=chain,
                    rnum=i + 1,
                    x=c[0],
                    y=c[1],
                    z=c[2],
                )
            )
            serial += 1


# -----------------------------------------------------------
# Build ANM
# -----------------------------------------------------------
def build_anm(pdb_path, modes=20):
    ag = parsePDB(pdb_path)
    sel = ag.select("calpha")
    if sel is None:
        raise ValueError(f"No CA atoms found in {pdb_path}")

    anm = ANM("ANM Model")
    anm.buildHessian(sel)
    anm.calcModes(modes)
    return anm


# -----------------------------------------------------------
# ΔRMSF - FIXED VERSION
# -----------------------------------------------------------
def compute_and_plot_delta_rmsf(anm_apo, anm_complex, outdir):
    """
    Calculate RMSF from ANM modes using manual calculation.
    """
    print("  > Calculating RMSF...")
    
    # Use first 10 modes for RMSF calculation
    modes_apo = anm_apo[:anm_apo.numModes()]
    modes_complex = anm_complex[:anm_complex.numModes()]
    
    # Calculate squared fluctuations manually from modes
    sqfluct_apo = np.zeros(modes_apo.numAtoms())
    sqfluct_complex = np.zeros(modes_complex.numAtoms())
    
    for i in range(modes_apo.numModes()):
        eigvec_apo = modes_apo[i].getEigvec()
        eigval_apo = modes_apo[i].getEigval()
        # Sum x, y, z components for each atom
        sqfluct_apo += eigval_apo * (eigvec_apo[::3]**2 + eigvec_apo[1::3]**2 + eigvec_apo[2::3]**2)
        
        eigvec_complex = modes_complex[i].getEigvec()
        eigval_complex = modes_complex[i].getEigval()
        sqfluct_complex += eigval_complex * (eigvec_complex[::3]**2 + eigvec_complex[1::3]**2 + eigvec_complex[2::3]**2)
    
    # RMSF is the square root of squared fluctuations
    rmsf_apo = np.sqrt(sqfluct_apo)
    rmsf_complex = np.sqrt(sqfluct_complex)

    np.savetxt(os.path.join(outdir, "rmsf_apo.txt"), rmsf_apo, fmt="%.6f")
    np.savetxt(os.path.join(outdir, "rmsf_complex.txt"), rmsf_complex, fmt="%.6f")

    delta = rmsf_complex - rmsf_apo
    np.savetxt(os.path.join(outdir, "delta_rmsf.txt"), delta, fmt="%.6f")

    plt.figure(figsize=(12, 5))
    plt.axhline(0, color="black", linewidth=0.8)
    plt.plot(delta, color="firebrick", linewidth=1.2)
    plt.xlabel("Residue Index (aligned)")
    plt.ylabel("ΔRMSF (Complex - Apo)")
    plt.title("Change in Residue Flexibility After Binding")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "delta_rmsf_plot.png"), dpi=300)
    plt.close()
    
    # Print summary statistics
    print(f"\n  ΔRMSF Summary:")
    print(f"    Mean: {np.mean(delta):.4f}")
    print(f"    Std:  {np.std(delta):.4f}")
    print(f"    Min:  {np.min(delta):.4f} (residue {np.argmin(delta)+1})")
    print(f"    Max:  {np.max(delta):.4f} (residue {np.argmax(delta)+1})")

    return delta


# -----------------------------------------------------------
# FIXED Overlap Calculation
# -----------------------------------------------------------
def compute_and_plot_overlap(anm_apo, anm_complex, nmodes, outdir):
    """
    Compute full overlap matrix (nmodes × nmodes) and extract:
    - Diagonal overlaps (mode i vs mode i)
    - Best-match overlaps (max overlap for each apo mode)
    """
    print(f"  > Computing full {nmodes}×{nmodes} overlap matrix...")

    # This returns an nmodes × nmodes matrix
    overlap_mat = calcOverlap(anm_apo[:nmodes], anm_complex[:nmodes])
    overlap_mat = np.asarray(overlap_mat)
    
    # Verify shape
    print(f"  > Overlap matrix shape: {overlap_mat.shape}")
    
    # Save full matrix
    np.save(os.path.join(outdir, "anm_overlap_matrix.npy"), overlap_mat)
    np.savetxt(os.path.join(outdir, "anm_overlap_matrix.txt"), overlap_mat, fmt="%.6f")

    # Extract diagonal (apo mode i vs complex mode i)
    diag = np.diag(overlap_mat)
    np.savetxt(os.path.join(outdir, "anm_overlap_diag.txt"), diag, fmt="%.6f")

    # Find best match for each apo mode (max absolute value across complex modes)
    max_per_mode = np.max(np.abs(overlap_mat), axis=1)
    best_idx = np.argmax(np.abs(overlap_mat), axis=1)

    np.savetxt(os.path.join(outdir, "anm_overlap_max_per_mode.txt"), max_per_mode, fmt="%.6f")
    np.savetxt(os.path.join(outdir, "anm_overlap_bestmatch_idx.txt"), best_idx + 1, fmt="%d")

    # Summary line plot (20 points, not 400!)
    modes_idx = np.arange(1, nmodes + 1)
    plt.figure(figsize=(10, 5))
    plt.plot(modes_idx, diag, marker="o", linewidth=1.5, markersize=5, 
             label="Diagonal Overlap (mode i vs i)")
    plt.plot(modes_idx, max_per_mode, marker="s", linestyle="--", linewidth=1.5, 
             markersize=5, label="Best-Match Overlap (max |overlap|)")
    plt.axhline(0, color="gray", linewidth=0.8, linestyle=":")
    plt.xlabel("Apo Mode Index", fontsize=11)
    plt.ylabel("Overlap", fontsize=11)
    plt.title(f"ANM Mode Overlap Summary ({nmodes} modes)", fontsize=12)
    plt.legend(frameon=True, loc="best")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "anm_overlap_summary.png"), dpi=300)
    plt.close()

    # Heatmap of full matrix
    print("  > Creating overlap heatmap...")
    try:
        import seaborn as sns
        plt.figure(figsize=(8, 7))
        sns.heatmap(
            overlap_mat,
            cmap="coolwarm",
            center=0,
            square=True,
            cbar_kws={"label": "Overlap"},
            xticklabels=5,  # Show every 5th label
            yticklabels=5,
        )
        plt.xlabel("Complex Mode Index", fontsize=11)
        plt.ylabel("Apo Mode Index", fontsize=11)
        plt.title(f"ANM Overlap Matrix ({nmodes}×{nmodes})", fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "anm_overlap_matrix_heatmap.png"), dpi=300)
        plt.close()
    except ImportError:
        plt.figure(figsize=(8, 7))
        im = plt.imshow(
            overlap_mat,
            cmap="coolwarm",
            aspect="auto",
            origin="lower",
            vmin=-1,
            vmax=1,
        )
        plt.colorbar(im, label="Overlap")
        plt.xlabel("Complex Mode Index", fontsize=11)
        plt.ylabel("Apo Mode Index", fontsize=11)
        plt.title(f"ANM Overlap Matrix ({nmodes}×{nmodes})", fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "anm_overlap_matrix_heatmap.png"), dpi=300)
        plt.close()

    # Summary text file
    with open(os.path.join(outdir, "anm_overlap_summary.txt"), "w") as f:
        f.write("ANM Overlap Summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Number of modes analyzed: {nmodes}\n")
        f.write(f"Overlap matrix shape: {overlap_mat.shape}\n\n")
        f.write(f"Mean diagonal overlap: {np.mean(diag):.6f}\n")
        f.write(f"Mean best-match overlap: {np.mean(max_per_mode):.6f}\n\n")
        f.write("Mode\tDiagonal\tBest-Match\tBest-Match-Index\n")
        f.write("-" * 60 + "\n")
        for i in range(nmodes):
            f.write(
                f"{i+1}\t{diag[i]:.6f}\t{max_per_mode[i]:.6f}\t{best_idx[i]+1}\n"
            )

    print(f"  > Mean diagonal overlap: {np.mean(diag):.4f}")
    print(f"  > Mean best-match overlap: {np.mean(max_per_mode):.4f}")

    return overlap_mat, diag, max_per_mode, best_idx


# -----------------------------------------------------------
# MAIN
# -----------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="ANM Apo–Complex Comparison with Chain Alignment"
    )
    parser.add_argument("--apo", required=True, help="Apo PDB file")
    parser.add_argument("--complex", required=True, help="Complex PDB file")
    parser.add_argument("--c_apo", required=True, help="Chain ID in apo")
    parser.add_argument("--c_complex", required=True, help="Chain ID in complex")
    parser.add_argument("--modes", type=int, default=20, help="Number of modes (default: 20)")
    parser.add_argument("--outdir", default="anm_comparison_aligned_results", 
                       help="Output directory")

    args = parser.parse_args()
    ensure_dir(args.outdir)

    print("\n" + "="*60)
    print("ANM APO-COMPLEX COMPARISON")
    print("="*60)
    
    print("\n[1/6] Parsing PDB files...")
    struct_apo = parsePDB(args.apo)
    struct_complex = parsePDB(args.complex)

    coordsA, resA, _ = atomgroup_ca_info(struct_apo, args.c_apo)
    coordsB, resB, _ = atomgroup_ca_info(struct_complex, args.c_complex)

    if len(coordsA) == 0 or len(coordsB) == 0:
        raise ValueError("Chain selections returned no CA atoms.")

    print(f"  Apo chain {args.c_apo}: {len(coordsA)} CA atoms")
    print(f"  Complex chain {args.c_complex}: {len(coordsB)} CA atoms")

    seqA = seq_from_resnames(resA)
    seqB = seq_from_resnames(resB)

    print("\n[2/6] Aligning sequences...")
    matchesA, matchesB = align_and_map(seqA, seqB)
    print(f"  Matched residues: {len(matchesA)}")

    # Build matched lists
    mA_coords = [coordsA[i] for i in matchesA]
    mA_resn   = [resA[i] for i in matchesA]

    mB_coords = [coordsB[i] for i in matchesB]
    mB_resn   = [resB[i] for i in matchesB]

    # Write matched CA PDBs
    print("\n[3/6] Writing matched CA-only PDB files...")
    apo_pdb = os.path.join(args.outdir, "matched_apo_ca.pdb")
    complex_pdb = os.path.join(args.outdir, "matched_complex_ca.pdb")

    write_ca_pdb(mA_coords, mA_resn, apo_pdb, chain="A")
    write_ca_pdb(mB_coords, mB_resn, complex_pdb, chain="B")
    print(f"  Saved: matched_apo_ca.pdb")
    print(f"  Saved: matched_complex_ca.pdb")

    print(f"\n[4/6] Building ANM models ({args.modes} modes)...")
    anm_apo = build_anm(apo_pdb, modes=args.modes)
    anm_complex = build_anm(complex_pdb, modes=args.modes)
    print(f"  Apo ANM: {anm_apo.numModes()} modes")
    print(f"  Complex ANM: {anm_complex.numModes()} modes")

    print(f"\n[5/6] Computing mode overlap...")
    compute_and_plot_overlap(anm_apo, anm_complex, args.modes, args.outdir)

    print(f"\n[6/6] Computing ΔRMSF...")
    compute_and_plot_delta_rmsf(anm_apo, anm_complex, args.outdir)

    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"Results saved in: {args.outdir}/")
    print("\nGenerated files:")
    files = [
        "matched_apo_ca.pdb",
        "matched_complex_ca.pdb",
        "anm_overlap_matrix.txt",
        "anm_overlap_matrix.npy",
        "anm_overlap_diag.txt",
        "anm_overlap_max_per_mode.txt",
        "anm_overlap_bestmatch_idx.txt",
        "anm_overlap_summary.png",
        "anm_overlap_matrix_heatmap.png",
        "anm_overlap_summary.txt",
        "rmsf_apo.txt",
        "rmsf_complex.txt",
        "delta_rmsf.txt",
        "delta_rmsf_plot.png",
    ]
    for fname in files:
        fpath = os.path.join(args.outdir, fname)
        if os.path.exists(fpath):
            print(f"  ✓ {fname}")


if __name__ == "__main__":
    main()

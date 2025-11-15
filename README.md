# DynaMune

## Unified Platform for Protein Dynamics with Integrated Immunoinformatics Tools

This platform provides an integrated, end-to-end environment for multi-epitope immunoinformatics, construct assembly, structural modeling, and advanced protein dynamics analysis. It combines epitope prediction, safety screening, construct design, 3D structure prediction, and ProDy-based Normal Mode and deformability analyses into a unified workflow. The tool supports both vaccine-oriented and general protein structural studies, including ligand-binding and allostery exploration.

---

## Features

### 1. Protein Structure and Dynamics Analysis (ProDy-Powered)

A comprehensive toolkit for exploring protein motions, flexibility, allostery, and interface behavior using Elastic Network Models (GNM/ANM) and deformation analysis.

#### A. Normal Mode Analysis (NMA)

* Supports both GNM and ANM with adjustable cutoffs and mode counts
* Outputs include:

  * Mean-square fluctuations (MSF), theoretical B-factors
  * Covariance and cross-correlation matrices
  * Eigenvalues and visualizable mode vectors
* Mode animations and conformer ensemble generation
* PCA projection of conformer space

#### B. Pocket Dynamics

* Detect pocket residues automatically
* Analyze pocket volume fluctuation across ANM modes
* Compute RMSF of pocket residues, pocket radius, and accessibility score
* Export ranked flexibility tables and pocket fluctuation plots

#### C. Perturbation Response Scanning (PRS)

* Apply systematic perturbations to individual residues
* Identify allosteric sensor and effector residues
* Visualize PRS sensitivity and effectiveness matrices (heatmaps)
* Export numeric matrices and top-ranked residue lists

#### D. Apo vs Holo Comparative Dynamics

* Upload apo (unbound) and holo (bound) structures
* Compute deformation vectors and project them onto ANM modes
* Identify active modes contributing to functional conformational changes

#### E. Interface Dynamics

* Automatically detect inter-chain or protein–ligand interface residues
* Quantify their local flexibility and dynamics
* Export interface mobility statistics

#### F. Domain & Hinge Detection

* Extract domain-wise motions from mode vectors
* Detect hinge regions driving major conformational transitions
* Visualize per-residue mobility and deformation arrows

#### G. Contact and Stability Analysis

* Analyze inter-chain residue contacts, hydrogen bonds, salt bridges, and interface stability
* Customizable cutoffs for:

  * Contact distance
  * Hydrogen bond angle and distance
  * Salt bridge electrostatic distance
* Outputs:

  * Plots and visualizations of interface metrics
  * Hydrophobic/hydrophilic interaction scores
  * Residue-wise contact heatmaps and PNG exports
* Supports dual-chain systems (A/B or auto-detect)
* Interactive form with pre-filled cutoff defaults and editable fields
* Full output (plots and data) downloadable as individual files or a zip archive

### H. Inter-Chain Contact Timeline Analysis

* Simulate contact evolution across an ANM-driven ensemble
* Track the persistence of inter-chain residue-residue contacts through conformers
* Outputs include:

  * **Contact Timeline Heatmap** — visualizing when each contact forms/breaks
  * **Contact Persistence Plot** — highlighting residue pairs with highest sustained interaction
  * **Contact Network Graph** — showing global inter-residue interaction network
* Parameters (user-editable or defaulted):

  * ANM cutoff distance
  * Number of conformers to generate
  * Contact distance threshold
  * Mode amplitude
  * Mode count
  * Persistence threshold
* Supports up to 2 chains (user-specified or auto-detected)
* Results can be downloaded individually or as a complete archive

---


### 2. Deformability Analysis (ProDy-Powered)

* Compare any two conformations (e.g., mutant vs wild-type, apo vs complex)
* Compute deformation vectors between structures
* Assess overlap with NMA modes to identify dominant contributors
* Export overlap plots and numeric reports

---

### 3. Integrated Immunoinformatics Suite

This module complements the dynamics analysis with a streamlined vaccine design workflow from epitope prediction to structure modeling.

#### A. Epitope Prediction and Filtering

* Predict B-cell epitopes using bundled NetBCE
* Predict MHC-I and MHC-II binding using IEDB (NetMHCpan & NetMHCIIpan EL)
* Filter using allergenicity (AlgPred2) and toxicity (ToxinPred3)
* Upload known epitopes directly if desired

#### B. Construct Assembly

* Combine epitopes into multi-epitope vaccine constructs
* Add adjuvants and standard linkers (AAY, GPGPG, EAAAK)
* Export as a clean FASTA sequence

#### C. Structure Prediction (for Constructs)

* Predict 3D structure of final construct using ESMFold
* PDB automatically stored and passed to NMA modules for structural analysis

---

### 4. Session-Based Execution

* Every analysis creates a unique `results/` folder containing:

  * Structure and epitope files
  * Dynamics plots, matrices, and conformers
  * Pocket and PRS results
  * Construct FASTA and predicted structure (if applicable)
* Ensures reproducibility and traceability across workflows

---

### 📁 Project Structure 

```
VaccineBuilder/
├── app.py
├── mhc1.py
├── mhc2.py
├── model.py
├── prody_nma.py
├── deformation_analysis.py
├── pocket_dynamics.py
├── prs_analysis.py
├── comparative_anm.py
├── contact_stability.py
├── contact_timeline.py          
├── requirements.txt
├── templates/
├── static/
├── NetBCE/
└── results/
```

---

## Workflows

### A. Sequence-Based Prediction

1. Upload sequence.
2. Predict epitopes.
3. Filter and annotate.
4. Design construct.
5. Predict structure.
6. Run NMA or Pocket/PRS as needed.

### B. Start from Epitopes

1. Paste known epitopes.
2. Skip prediction.
3. Continue with construct → structure → dynamics.

### C. Standalone NMA & Pocket Analysis

1. Upload any protein PDB.
2. Select ANM/GNM.
3. Run dynamics, pocket, PRS, etc.
4. Download heatmaps, PDBs, tables.

### D. Deformability Analysis (Comparative)

1. Upload unbound and bound PDBs.
2. Run deformation vector + mode projection.
3. Analyze contributing modes.

### E. Contact and Stability Analysis

1. Upload multi-chain PDB (max 2 chains).
2. Adjust contact and interaction cutoffs (or keep defaults).
3. View and download:

   * PNG interaction plots
   * Hydrophobicity/charge-based metrics
   * Contact matrices and summaries
4. Download all results or zip archive

### F. Inter-Chain Contact Timeline

1. Upload multi-chain PDB structure (2 chains max)
2. Choose ANM parameters (cutoff, modes, amplitude)
3. Run simulation to generate conformer ensemble
4. View timeline heatmap, contact persistence graph, and contact network
5. Export figures or all results together

---

## Requirements

* Python 3.7+
* Flask, ProDy, Matplotlib, Biopython, Pandas
* External access for IEDB APIs
* NetBCE locally available

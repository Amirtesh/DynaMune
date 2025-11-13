# Unified Multi-Epitope Immunoinformatics and Structural Analysis Platform

This platform provides an integrated, end-to-end environment for multi-epitope immunoinformatics, construct assembly, structural modeling, and protein dynamics analysis. It combines epitope prediction, safety screening, construct design, 3D structure prediction, and ProDy-based Normal Mode and deformability analyses into a single workflow. The tool supports both vaccine-oriented and general protein analysis use cases.

## Features

### 1. Epitope Prediction

* B-cell epitope prediction using bundled NetBCE (local execution).
* MHC-I binding prediction using IEDB Next-Generation API (NetMHCpan 4.1 EL model).
* MHC-II binding prediction using IEDB Next-Generation API (NetMHCIIpan 4.1 EL model).
* Configurable allele selection for MHC-I and MHC-II.
* Automatic parsing and export of epitope prediction tables.

### 2. Safety and Filtering

* Allergenicity prediction using AlgPred2.
* Toxicity prediction using ToxinPred3.
* Annotation of all predicted epitopes with safety metadata.
* Workflow for starting directly from user-provided epitopes.

### 3. Construct Assembly

* Flexible construct builder supporting multiple design patterns.
* Library of adjuvants.
* Automated linker insertion (AAY, GPGPG, EAAAK, and flexible linkers).
* Final construct exported in FASTA format.

### 4. Physicochemical Characterization

* ProtParam-like computation using Biopython.
* Outputs molecular weight, theoretical pI, instability index, aliphatic index, GRAVY, and amino acid composition.

### 5. Structure Prediction

* Integrated ESMFold-based prediction backend.
* Automated PDB file generation.
* Structures saved to session-specific results directory.

### 6. ProDy-Based Structural Dynamics and NMA

The platform integrates a fully featured ProDy pipeline including:

#### Normal Mode Analysis (NMA)

* ANM or GNM with configurable cutoffs.
* Calculation of mean-square fluctuations and theoretical B-factors.
* Mode shape generation and visualization.
* Covariance and cross-correlation matrices.
* Chain-aware analysis.
* Interface residue detection with user-defined distance thresholds.

#### Conformer Ensemble Generation and PCA

* ANM-driven conformer ensemble generation.
* Principal Component Analysis (PCA) of ensemble.
* Comparison of ANM modes with PCA eigenvectors.

#### Allosteric and Interface Dynamics

* Allosteric pathway analysis for ANM models.
* Interface dynamics quantification.
* Export of all results, plots, and numerical matrices.

### 7. Deformability Analysis

* Reference vs target structure comparison.
* Computation of deformation vectors.
* Squared overlap and cumulative overlap against NMA modes.
* Identification of key contributing modes.
* Plots and numerical data exported.

### 8. Session-Based Execution and Reproducibility

* Each analysis session receives a dedicated directory under `results/` containing:

  * Epitope prediction tables
  * Final construct
  * Structure prediction output (PDB)
  * NMA and PCA analysis results
  * Deformability analysis output
* Ensures reproducibility and organized storage.

## Project Structure

```
VaccineBuilder/
├── app.py                      # Flask backend
├── mhc1.py                     # MHC-I prediction wrapper for IEDB API
├── mhc2.py                     # MHC-II prediction wrapper for IEDB API
├── model.py                    # ESMFold structure prediction
├── prody_nma.py                # ProDy Normal Mode Analysis workflow
├── deformation_analysis.py     # ProDy-based deformability comparison
├── requirements.txt            # Python dependencies
├── templates/                  # HTML templates
├── static/                     # CSS/JS assets
├── NetBCE/                     # Bundled NetBCE model and scripts
└── results/                    # Session-specific analysis outputs
```

## Workflows

### Sequence-Based Prediction Workflow

1. Input FASTA sequence.
2. Select alleles for MHC-I and MHC-II.
3. Run epitope predictions.
4. Review and filter epitopes.
5. Perform safety screening.
6. Assemble construct.
7. Run ProtParam analysis.
8. Predict structure with ESMFold.
9. Run NMA and/or deformability analysis.
10. Export all results.

### Start from Predicted Epitopes

1. Provide precomputed B-cell, CTL, and HTL epitopes.
2. Perform safety filtering.
3. Assemble construct and run downstream analysis.

### Standalone Normal Mode Analysis

1. Upload any PDB structure (single or multi-chain).
2. Configure ANM/GNM parameters.
3. Perform full NMA, PCA, allosteric, interface, and ensemble analysis.
4. Download outputs.

### Standalone Deformability Analysis

1. Upload reference (unbound) and target (bound) structures.
2. Specify chain IDs.
3. Compute deformation vectors.
4. Calculate squared overlaps and cumulative contributions.
5. Export results.

## Dependencies

* Python 3.7+
* Flask
* ProDy
* Matplotlib
* Biopython
* Requests
* Pandas
* NetBCE (bundled)

## Notes

* IEDB Next-Generation API access requires an active internet connection.
* ESMFold prediction uses the platform-defined model backend.
* The platform is intended for computational research and does not provide clinical recommendations.

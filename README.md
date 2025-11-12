# 🧬 Multi-Epitope Vaccine Builder

A professional **Flask-based** web platform for designing multi-epitope vaccines using state-of-the-art **B-cell** and **T-cell** epitope prediction methods and modern machine-learning models (**NetBCE**, **NetMHCpan 4.1**, **NetMHCIIpan 4.1**).

---

## 🎯 Overview

The Multi-Epitope Vaccine Builder enables researchers to:

* Predict B-cell and T-cell epitopes from protein sequences
* Select and combine promising epitopes
* Design vaccine constructs with multiple architectural patterns
* Analyze physicochemical properties (ProtParam)
* Export predictions (CSV) and final constructs (FASTA)

---

## 🌟 Key Features

### 1) Dual Workflow Support

* **Prediction Workflow:** Start from a protein FASTA sequence and run integrated predictors.
* **Direct Input Workflow:** Paste pre-predicted epitopes and proceed directly to construct design.

### 2) Advanced Epitope Prediction

* **B-Cell (Linear):** **NetBCE** deep learning model.
* **MHC-I (CTL):** **NetMHCpan 4.1** (8–11mers; wide HLA coverage).
* **MHC-II (HTL):** **NetMHCIIpan 4.1** (typically 15mers; DR/DP/DQ; multi-allele support).

### 3) Interactive Results

* Paginated DataTables with sort/filter/search.
* Cross-page selection with persistent state and real-time selection counters.
* One-click CSV export.

### 4) Flexible Construct Design

Choose from **7 construct patterns** (sequential and alternating):

**Sequential**

1. CTL → HTL → B-cell (default)
2. HTL → CTL → B-cell
3. B-cell → HTL → CTL

**Alternating** 

4. CTL-HTL-B (repeating)
5. HTL-B-CTL (repeating)
6. CTL-B-HTL (repeating)
7. B-CTL-HTL (repeating)

**Validated linker set**

* **EAAAK:** Adjuvant separator (rigid)
* **AAY:** CTL (MHC-I) epitope linker
* **GPGPG:** HTL (MHC-II) epitope linker (flexible)
* **KK:** B-cell epitope linker (cleavage-prone)

### 5) Adjuvant Integration (9 options)

* β-defensin (human defensin-3)
* PADRE (Pan HLA-DR epitope)
* Cholera Toxin B subunit (CTB)
* 50S Ribosomal Protein L7/L12 (M. tuberculosis)
* RS09 (TLR4 agonist peptide)
* HBHA fragment
* Flagellin (FliC fragment)
* Hsp70 fragment
* HABA peptide

### 6) Protein Analysis (ProtParam)

* Molecular weight, AA composition, theoretical pI
* Instability index, aliphatic index, GRAVY
* Secondary structure fractions (helix/turn/sheet)
* Extinction coefficients

### 7) 3D Structure Prediction & Analysis

* **Interactive 3D Visualization**: Predict protein structure using ESMFold API
* **3Dmol.js Integration**: View predicted structures in interactive 3D viewer with cartoon representation
* **PDB Download**: Export predicted structures in PDB format
* **ProDy Coarse Grained Normal Mode Analysis (CG-NMA)**: 
  * Analyze protein flexibility and dynamics using ANM (Anisotropic Network Model) or GNM (Gaussian Network Model)
  * Support for custom PDB structure uploads
  * Adjustable number of modes (1-100, default: 20)
  * Configurable number of conformers (1-50, default: 10)
  * Configurable interface cutoff (5.0-30.0 Å, default: 10.0)
  * **Advanced Analysis Features** (always enabled):
    * **Ensemble PCA**: Principal Component Analysis for ensemble dynamics
    * **Allosteric Analysis**: Identify allosteric communication pathways
    * **ClustENMD**: Cluster-based Elastic Network Model analysis
  * Comprehensive analysis including MSF, B-factors, eigenvalues, cross-correlation, covariance, and collectivity
  * PCA comprehensive analysis with ANM-PCA comparison plots
  * Mode shape visualization and conformer generation
  * Downloadable results (PNG plots and complete ZIP package)

### 8) Deformability Analysis

* **Conformational Change Analysis**: Compare intrinsic dynamics of unbound vs. bound protein structures
* **ProDy-based Analysis**: Calculate overlap between NMA modes and observed conformational changes
* **Dual Structure Upload**: Upload reference (unbound) and target (bound complex) PDB structures
* **Chain-specific Analysis**: Specify chain IDs for precise protein comparison
* **Deformation Vector Calculation**: Measure conformational changes upon binding
* **Mode Overlap Quantification**: Determine which normal modes contribute to conformational changes
* **Comprehensive Visualization**: 
  * Individual mode contributions (squared overlap)
  * Cumulative overlap analysis
  * 80% threshold identification for key modes
* **Downloadable Results**: PNG plots and detailed numerical results

---

## 📁 Project Structure

```
VaccineBuilder/
├── app.py                      # Flask app (routes & logic)
├── mhc1.py                     # MHC-I predictions (IEDB API wrapper)
├── mhc2.py                     # MHC-II predictions (IEDB API wrapper)
├── model.py                    # ESMFold 3D structure prediction
├── prody_nma.py                # ProDy Coarse Grained Normal Mode Analysis
├── deformation_analysis.py     # ProDy Deformability Analysis
├── requirements.txt            # Python dependencies
├── templates/                  # Jinja2 HTML templates
│   ├── index.html
│   ├── results.html
│   ├── construct.html
│   ├── prody_nma.html
│   ├── standalone_nma.html
│   └── deformability.html
├── static/                     # CSS/JS assets
│   ├── css/
│   │   └── style.css
│   └── js/
│       ├── main.js
│       ├── results.js
│       └── construct.js
├── NetBCE/                     # NetBCE tool & models
│   ├── prediction/
│   │   └── NetBCE_prediction.py
│   ├── models/
│   └── result/
└── results/                    # Generated outputs (session-based folders)
    └── {YYYYMMDD_HHMMSS_xxxxxxxx}/
        ├── predicted_structure.pdb
        ├── *_NMA_PCA_results/       # ProDy CG-NMA analysis outputs (with PCA)
        └── *_Deformation_Analysis/  # Deformability analysis outputs
```

---

## 📋 Prerequisites

* **Python**: 3.7+
* **pip** (package installer)
* **External**: Internet access for IEDB API (mhc1.py/mhc2.py) and ESMFold API (model.py)
* **Bundled**: NetBCE tool under `NetBCE/`
* **Dependencies**: ProDy, Matplotlib (installed via requirements.txt)

---

### Workflow A — Start from Prediction

1. **Input**: Paste protein sequence in FASTA (header line begins with `>`).
2. **Select alleles**:

   * MHC-I: e.g., `HLA-A*02:01,HLA-A*24:02`
   * MHC-II: e.g., `HLA-DRB1*01:01,HLA-DRB1*03:01`
3. **Run**: Start predictions (NetBCE, NetMHCpan, NetMHCIIpan).
4. **Review**: Explore tables, filter/sort, select epitopes.
5. **Proceed**: Move to construct screen.

### Workflow B — Start from Predicted Epitopes

1. Choose “Start from Predicted Epitopes.”
2. Paste epitopes (one per line) for **B-cell**, **CTL**, **HTL**.
3. Proceed directly to construct design.

### Build the Construct

1. **Pick an adjuvant** (from the 9 options).
2. **Choose a pattern** (one of the 7 designs).
3. **Build** and inspect the final sequence with linkers.
4. **Analyze** (optional): run ProtParam.
5. **Predict 3D Structure** (optional): Generate and visualize protein structure interactively.
6. **Run ProDy CG-NMA Analysis** (optional): Perform Coarse Grained Normal Mode Analysis for flexibility and dynamics insights.
7. **Export**: copy to clipboard, **Download FASTA**, **Download CSV**, **Download PDB**, **Download Analysis Plots**.

### Standalone CG-NMA Workflow

1. **Upload Custom PDB**: Navigate to the "Coarse Grained Normal Mode Analysis" workflow from the homepage.
2. **Upload Structure**: Upload your own PDB file (single or multi-chain complexes supported).
3. **Configure Analysis**: 
   * Choose ANM or GNM method
   * Set number of modes (1-100, default: 20)
   * Set number of conformers (1-50, default: 10)
   * Adjust interface cutoff (5.0-30.0 Å, default: 10.0)
   * PCA, Allosteric, and ClustENMD analyses are automatically enabled
4. **Run Analysis**: Generate comprehensive flexibility analysis including:
   * MSF, B-factors, eigenvalues, cross-correlation, covariance, and mode shapes
   * PCA comprehensive analysis with ANM-PCA comparison
   * ClustENMD cluster-based analysis
5. **Download Results**: Export individual plots (Comprehensive Analysis, Mode Shapes, PCA Comprehensive, ANM-PCA Comparison, ClustENMD) or complete ZIP package.

### Deformability Analysis Workflow

1. **Access Workflow**: Navigate to "Deformability Analysis" from the homepage.
2. **Upload Structures**: 
   * Reference PDB: Unbound vaccine structure (standalone protein)
   * Target PDB: Vaccine in complex (bound state, e.g., with MHC)
3. **Specify Chains**: Enter single-letter chain IDs for both structures.
4. **Configure Parameters**: Set number of NMA modes (1-100, default: 20).
5. **Run Analysis**: System performs:
   * Structure matching and alignment
   * ANM calculation on reference structure
   * Deformation vector calculation
   * Mode overlap quantification
6. **View Results**: Interactive plot showing individual and cumulative mode contributions.
7. **Download**: Export deformation overlap plot and numerical results.

---

## 📊 Outputs

**Generated files (in `results/{session_id}/`):**

* `B_Cell_Epitopes.csv` — NetBCE predictions
* `mhc1_prediction.csv` — MHC-I binding predictions
* `mhc2_prediction.csv` — MHC-II binding predictions
* `vaccine_construct.fasta` — Final vaccine construct
* `predicted_structure.pdb` — 3D structure prediction (ESMFold)
* `*_NMA_PCA_results/` — ProDy CG-NMA outputs with PCA
  * `comprehensive_analysis.png` — Multi-panel analysis (MSF, B-factors, eigenvalues, cross-correlation, covariance, collectivity)
  * `mode_shapes_combined.png` — Grid of first 6-12 mode shapes
  * `pca_comprehensive_analysis.png` — PCA analysis with variance and PC projections
  * `anm_pca_comparison.png` — Comparison between ANM and PCA modes
  * `clustenmd_analysis.png` — Cluster-based Elastic Network Model analysis
  * `individual_mode_*.png` — Individual mode shape plots
  * `conformers.pdb` — Generated conformers
  * Text files with detailed analysis data
* `*_Deformation_Analysis/` — Deformability analysis outputs
  * `deformation_overlap.png` — Mode contribution and cumulative overlap plot
  * `deformation_results.txt` — Numerical overlap data for each mode

**Formats**

* **CSV:** Comma-separated with headers
* **FASTA:** Standard, wrapped at 80 chars/line
* **PDB:** Protein Data Bank format for 3D structures
* **PNG/SVG:** High-quality plots for analysis and publication

---

## 🔬 Scientific Background

### Epitope Predictors

* **NetBCE**: Deep learning on experimental B-cell epitope data; outputs linear epitope scores.
* **NetMHCpan 4.1**: Pan-specific MHC-I binding prediction (HLA-A/B/C; 8–11mers).
* **NetMHCIIpan 4.1**: Pan-specific MHC-II binding prediction (HLA-DR/DP/DQ; ~15mers).

### Design Rationale

* **Sequential patterns** may amplify specific arm(s) of the response.
* **Alternating patterns** can promote a more balanced immune response.
* **Linkers** prevent unintended fusion and help preserve structural/antigenic independence:

  * EAAAK (rigid adjuvant spacer), AAY (CTL), GPGPG (HTL), KK (B-cell).

---

## 📝 Best Practices

### Epitope Selection

1. Prefer high prediction scores and binder ranks.
2. Cover diverse HLA alleles for population breadth.
3. Balance B-cell, CTL, HTL epitope counts (often 3–5 each).
4. Exclude sequences with close human homology.

### Construct Design

1. Start with **sequential** patterns for simpler analysis.
2. Explore **alternating** patterns for balance.
3. Keep length practical (≈200–500 aa).
4. Validate stability/solubility (ProtParam).

### Adjuvant Choice

* **PADRE**: broadly effective for T-cell help.
* **β-defensin**: antimicrobial and immunostimulatory properties.
* **CTB**: strong for mucosal immunity.

---

## ⚠️ Error Handling & Troubleshooting

**Validation**

* FASTA headers/characters, allele formats, and inputs are validated.
* Subprocess and timeout handling (≈5 minutes/prediction); clear user-level errors.

**Common Issues**

* *Slow predictions*: long sequences increase runtime (especially NetMHCPan and NetMHCIIPan since they use an API).
* *No results*: ensure all three predictors complete; check browser console.
* *Selection not persisting*: refresh; selection is managed globally.
* *ProtParam error*: ensure only valid amino acids are present.
* *IEDB timeout*: verify internet connectivity; service downtime happens.
* *3D structure prediction slow*: ESMFold API can take 1-3 minutes depending on sequence length and server load.
* *ProDy CG-NMA timeout*: Large proteins or high mode counts may take longer; reduce number of modes or conformers for faster analysis.
* *Deformability analysis chain matching*: Ensure correct chain IDs are specified; chains must have sufficient sequence similarity (>25%) and overlap (>70%) for successful matching.
* *3D viewer lag*: Disable surface rendering (already optimized); ensure modern browser with WebGL support.

---

> *Design multi-epitope vaccines using advanced neural-network algorithms (NetBCE, NetMHCpan 4.1, NetMHCIIpan 4.1).*

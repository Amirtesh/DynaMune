# Quick Start Guide

## Get DynaMune Running in 5 Minutes

### Step 1: Install

```bash
# Clone the repository
git clone https://github.com/Amirtesh/DynaMune.git
cd DynaMune

# Create and activate environment
conda create -n dynamune python=3.10
conda activate dynamune

# Install DynaMune
pip install -e .
```

### Step 2: Verify Installation

```bash
# Check that commands are available
dynamune-web --help
dynamune-nma --help
```

### Step 3: Start Using DynaMune

#### Option A: Web Interface (Recommended for Beginners)

```bash
# Start the local server
dynamune-web

# Open your browser to: http://127.0.0.1:5000
```

Navigate through the interface:
- **Normal Mode Analysis**: Upload a PDB file, configure parameters, run analysis
- **Pocket Dynamics**: Analyze binding pocket flexibility
- **PRS Analysis**: Identify allosteric residues
- **Vaccine Design**: Start with antigen sequence, design multi-epitope constructs

#### Option B: Command-Line (for Batch Processing)

```bash
# Download a test structure
wget https://files.rcsb.org/download/1AKE.pdb

# Run Normal Mode Analysis
dynamune-nma 1AKE.pdb --method ANM --conformers 50

# Results will be in: 1AKE_NMA_PCA_results/
```

### Step 4: Explore Example Workflows

#### Example 1: Comprehensive Protein Dynamics

```bash
# Analyze a protein structure
PDB="protein.pdb"

# Normal modes and ensemble
dynamune-nma $PDB --method ANM --ensemble-pca --conformers 100

# Pocket dynamics
dynamune-pocket $PDB --auto-detect

# Allosteric analysis
dynamune-prs $PDB --cutoff 15

# Domain decomposition
dynamune-domain-hinge $PDB
```

#### Example 2: Comparative Apo vs Holo

```bash
# Compare two conformational states
dynamune-anm-compare apo.pdb holo.pdb --modes 20
dynamune-deformation apo.pdb holo.pdb
```

#### Example 3: Interface Stability

```bash
# Analyze protein-protein complex
dynamune-contact-stability complex.pdb --chains A,B
dynamune-contact-timeline complex.pdb --chains A,B --conformers 100
```

### Step 5: Access Your Results

Results are organized by analysis type:

```
results/
└── [session_id]/
    ├── [protein]_NMA_PCA_results/
    │   ├── plots/
    │   ├── data/
    │   └── conformers.pdb
    ├── [protein]_PRS_results/
    └── ...
```

All outputs include:
- **PNG plots**: Publication-ready visualizations
- **CSV/TXT data**: Numerical results for further analysis
- **PDB files**: Conformer ensembles, modified structures

---

## Common Use Cases

### 1. I want to analyze protein flexibility

**Web Interface:**
1. Go to "Normal Mode Analysis"
2. Upload PDB
3. Select ANM method
4. Click "Analyze"
5. Download RMSF plots and B-factors

**Command Line:**
```bash
dynamune-nma protein.pdb --method ANM
```

### 2. I want to identify allosteric sites

**Web Interface:**
1. Go to "PRS Analysis"
2. Upload PDB
3. Set cutoff (default: 15 Å)
4. View sensor/effector residues

**Command Line:**
```bash
dynamune-prs protein.pdb --cutoff 15
```

### 3. I want to design a multi-epitope vaccine

**Web Interface:**
1. Go to "Epitope Prediction"
2. Upload antigen sequence
3. Predict B-cell and T-cell epitopes
4. Filter by allergenicity/toxicity
5. Assemble construct
6. Predict structure
7. Analyze dynamics

### 4. I need to compare apo and holo structures

**Command Line:**
```bash
# Mode overlap analysis
dynamune-anm-compare apo.pdb holo.pdb

# Deformation analysis
dynamune-deformation apo.pdb holo.pdb

# Check both results directories
```

---

## Tips and Tricks

### Performance Optimization

```bash
# For large proteins (>1000 residues):
# - Reduce conformer count
dynamune-nma large_protein.pdb --conformers 20

# - Use GNM instead of ANM (faster)
dynamune-nma large_protein.pdb --method GNM

# - Reduce mode count
dynamune-nma large_protein.pdb --modes 10
```

### Customizing Output

```bash
# Specify output directory
dynamune-nma protein.pdb --output-dir my_results/

# Adjust parameters
dynamune-nma protein.pdb --cutoff 12 --gamma 0.8 --modes 20
```

### Batch Processing

```bash
# Process all PDB files in a directory
for pdb in structures/*.pdb; do
    basename=$(basename "$pdb" .pdb)
    dynamune-nma "$pdb" --output-dir "results/$basename"
done
```

---

## Getting Help

```bash
# General help
dynamune-web --help

# Tool-specific help
dynamune-nma --help
dynamune-prs --help
dynamune-pocket --help
```

**Documentation:**
- Full README: `README.md`
- GitHub Issues: https://github.com/Amirtesh/DynaMune/issues

**Common Issues:**
- Command not found → Ensure virtual environment is activated
- Import errors → Reinstall with `pip install -e .`
- Port in use → Try different port: `dynamune-web --port 8080`

---

## What's Next?

- 📖 Read the full [README.md](README.md) for comprehensive documentation
- 🔬 Explore advanced features in each analysis module
- 🧬 Try the vaccine design workflow
- 📊 Integrate results with your visualization pipeline (PyMOL, Chimera)
- 🤝 Contribute improvements via GitHub

---

**Happy analyzing! 🚀**

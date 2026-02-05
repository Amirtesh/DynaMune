# DynaMune

**Comprehensive Protein Dynamics Analysis Platform**

DynaMune is a local-first scientific software tool for analyzing protein dynamics, flexibility, and allosteric communication using Normal Mode Analysis (NMA) and related computational methods. Built on ProDy, it provides both an intuitive web interface and powerful command-line tools for structural bioinformatics research.

## Key Features

- **Normal Mode Analysis (NMA)** - Predict protein flexibility and conformational dynamics
- **Perturbation Response Scanning (PRS)** - Identify allosteric hotspots and signaling pathways  
- **Domain & Hinge Detection** - Automatically identify rigid domains and flexible hinges
- **Pocket Dynamics** - Analyze binding site flexibility and druggability
- **Contact Analysis** - Track inter-residue interactions and stability
- **Deformation Analysis** - Quantify local and global structural changes
- **ANM Comparison** - Compare dynamics between different protein states

## Installation

### Prerequisites
- Python 3.9 or higher
- 4GB+ RAM (16GB+ recommended for large proteins)

### Quick Install

```bash
# Clone the repository
git clone https://github.com/Amirtesh/DynaMune.git
cd DynaMune

# Create conda environment (recommended)
conda create -n dynamune python=3.10
conda activate dynamune

# Install DynaMune
pip install -e .

# Verify installation
python3 verify_install.py
```

### Alternative: Using the automated script

```bash
bash install.sh
```

## Usage

### Web Interface (Recommended for Most Users)

Start the interactive web application:

```bash
dynamune-web
```

Then open your browser to `http://localhost:5000`

**Web Interface Options:**
```bash
dynamune-web --host 0.0.0.0 --port 8080      # Custom host/port
dynamune-web --debug                          # Enable debug mode
```

### Command-Line Tools

For batch processing, automation, or large-scale analyses:

#### 1. Normal Mode Analysis
```bash
dynamune-nma structure.pdb --method ANM --modes 20
```

#### 2. Perturbation Response Scanning
```bash
dynamune-prs structure.pdb --residue A:100 --output prs_results/
```

#### 3. Pocket Dynamics
```bash
dynamune-pocket structure.pdb --pocket-residues "A:50,A:51,A:52"
```

#### 4. Domain & Hinge Analysis
```bash
dynamune-domain-hinge structure.pdb --modes 10
```

#### 5. Deformation Analysis
```bash
dynamune-deformation structure.pdb --modes 5
```

#### 6. Contact Stability
```bash
dynamune-contact-stability structure.pdb --cutoff 4.5
```

#### 7. Contact Timeline (Multi-structure)
```bash
dynamune-contact-timeline trajectory.dcd --topology structure.pdb
```

#### 8. ANM Comparison
```bash
dynamune-anm-compare apo.pdb holo.pdb --output comparison/
```

### Getting Help

**Open comprehensive documentation:**
```bash
dynamune-help
```

This opens detailed use cases, scientific validation, FAQs, and workflow examples in your browser.

**Command-specific help:**
```bash
dynamune-nma --help
dynamune-prs --help
# etc.
```

## Quick Start Example

**Analyze a protein structure in 3 steps:**

```bash
# 1. Download a structure (e.g., PDB ID: 1ubq)
wget https://files.rcsb.org/download/1UBQ.pdb

# 2. Run Normal Mode Analysis
dynamune-nma 1UBQ.pdb --method ANM --modes 20

# 3. View results in results/ directory
ls results/
```

Or use the web interface for interactive visualization and analysis.

## Project Structure

```
DynaMune/
├── dynamune/           # Main package
│   ├── cli/           # Command-line tools
│   └── web.py         # Web interface launcher
├── templates/         # HTML templates for web UI
├── static/           # CSS, JavaScript, assets
├── NetBCE/           # B-cell epitope prediction module
├── results/          # Analysis outputs (auto-created)
├── setup.py          # Package configuration
└── requirements.txt  # Dependencies
```

## System Requirements

- **Minimum:** Python 3.9, 4GB RAM, 2 CPU cores
- **Recommended:** Python 3.10+, 16GB RAM, 4+ CPU cores
- **Large proteins (>2000 residues):** 32GB+ RAM recommended

## Output Files

All analyses generate results in the `results/` directory:

- **PDB files** - Modified structures, conformational ensembles
- **CSV/TXT** - Numerical data, residue scores, matrices
- **PNG/SVG** - Plots, heatmaps, visualizations
- **JSON** - Structured metadata for programmatic access

## Scientific Foundation

DynaMune implements methods based on:

- **Elastic Network Models (ENM)** - Anisotropic Network Model (ANM), Gaussian Network Model (GNM)
- **ProDy framework** - Industry-standard protein dynamics toolkit
- **Published algorithms** - PRS (Atilgan et al.), domain detection (Emekli et al.), and more

Validated against experimental data including:
- X-ray crystallography B-factors
- NMR dynamics measurements  
- Hydrogen-deuterium exchange (HDX-MS)
- Literature-reported allosteric sites

## Citation

If you use DynaMune in your research, please cite:

```bibtex
@software{dynamune2025,
  title={DynaMune: Comprehensive Protein Dynamics Analysis Platform},
  author={Amirtesh},
  year={2025},
  url={https://github.com/Amirtesh/DynaMune}
}
```

## Troubleshooting

**Installation issues:**
```bash
# Update pip first
pip install --upgrade pip

# Install with verbose output
pip install -e . -v
```

**"Command not found" errors:**
```bash
# Reinstall package
pip install -e . --force-reinstall
```

**Import errors (ProDy, Flask, etc.):**
```bash
# Reinstall dependencies
pip install -r requirements.txt --force-reinstall
```

**Large protein timeout:**
- Use CLI tools instead of web interface for >2000 residue proteins
- Reduce number of modes calculated (`--modes 10` instead of default 20)
- Consider domain-based analysis for very large complexes

## Support

For comprehensive documentation, use cases, and FAQs:
```bash
dynamune-help
```

Report issues: https://github.com/Amirtesh/DynaMune/issues

## License

DynaMune is released under the MIT License.


---

**Made for structural biologists, by structural biologists.**  
Built on ProDy | Runs locally | No cloud dependencies | Open source

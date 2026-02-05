# DynaMune

**DynaMune** is an integrated, ensemble-based protein dynamics analysis platform
built on elastic network models (ENM) and normal mode analysis (NMA). It provides
a unified, reproducible workflow for analyzing collective motions, allosteric
coupling, domain and hinge behavior, deformation patterns, pocket dynamics, and
protein–protein interface interactions without the computational cost of
molecular dynamics simulations.

The platform is implemented in Python, uses ProDy as its computational backbone,
and supports both command-line and optional web-based interfaces.

---

## Scope of the JOSS submission

The JOSS paper associated with this repository describes the **core DynaMune
protein dynamics framework**, implemented under the `dynamune/` package and its
associated command-line tools.

The Flask-based web interface and immunoinformatics utilities (e.g., NetBCE, MHC
modules) included in this repository are **optional extensions** provided for
interactive exploration and demonstration purposes. These components are **not
required** for core DynaMune functionality and are **not part of the primary
software scope evaluated in the JOSS submission**.

---

## Key features

### Core protein dynamics analysis
- Elastic network models (GNM and ANM)
- Normal mode analysis (NMA)
- Principal component analysis (PCA)
- Ensemble generation for statistical interpretation
- Perturbation response scanning (PRS)
- Domain and hinge residue identification
- Deformation and flexibility mapping
- Pocket dynamics analysis
- Apo–complex comparative dynamics

### Interface analysis
- Automated identification of protein–protein interface contacts
- Non-covalent interaction extraction (hydrogen bonds, salt bridges, disulfides)
- Interface contact stability and summary statistics

### Interfaces
- Command-line interface for reproducible batch workflows
- Optional Flask-based web interface for interactive exploration

---

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/Amirtesh/DynaMune.git
cd DynaMune
pip install -r requirements.txt
```

Alternatively, install using the provided setup script:

```bash
bash install.sh
```

## Usage

### Primary interface: command-line interface (CLI)

The primary and recommended interface for DynaMune is the command-line interface,
which supports reproducible and scriptable workflows.

Example:

```bash
dynamune nma --pdb input.pdb --chain A
```

Available modules include:

- `nma`
- `prs`
- `domain_hinge`
- `deformation`
- `pocket_dynamics`
- `contact_stability`
- `contact_timeline`

Use the built-in help for details:

```bash
dynamune help
```

### Optional web interface (interactive frontend)

An optional Flask-based web interface is provided for interactive analysis and
visual exploration. This interface is not required for batch execution or
reproducible workflows.

To launch the web application:

```bash
python app.py
```

The web interface exposes selected DynaMune modules through a browser-based UI
and is intended for exploratory use and demonstrations.

---

## Verification

Basic installation and functionality can be verified using:

```bash
python verify_install.py
```

This script checks core dependencies and validates that key analysis modules can
be executed successfully.

Representative validation case studies are documented in the JOSS paper and
accompanying figures.

---

## Repository structure (simplified)

```
DynaMune/
├── dynamune/           # Core DynaMune package (CLI + analysis modules)
├── figures/            # Figures used in the JOSS paper
├── paper.md            # JOSS manuscript
├── paper.bib           # Bibliography
├── app.py              # Optional Flask web interface
├── requirements.txt
├── README.md
└── LICENSE
```

Additional directories (e.g., NetBCE, MHC utilities) provide optional extensions
and are not required for core protein dynamics analysis.

---

## Validation overview

DynaMune has been evaluated on representative benchmark systems commonly used in
protein dynamics studies, including adenylate kinase and the ACE2–SARS-CoV-2
spike complex. These case studies demonstrate recovery of canonical domain
motions, hinge-localized flexibility patterns, and automated extraction of
interface interaction features consistent with prior literature.

Detailed validation results are presented in the JOSS manuscript.

---

## License

This project is released under the MIT License. See the LICENSE file for
details.

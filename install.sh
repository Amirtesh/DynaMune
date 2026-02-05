#!/bin/bash
#
# DynaMune Installation Script
# 
# This script automates the installation of DynaMune and its dependencies.
# Run with: bash install.sh
#

set -e  # Exit on error

echo "======================================================================"
echo "  DynaMune Installation Script"
echo "======================================================================"
echo ""

# Check Python version
echo "Checking Python version..."
python3 --version || {
    echo "Error: Python 3 not found. Please install Python 3.9 or higher."
    exit 1
}

# Check if we're in a virtual environment
if [[ -z "${VIRTUAL_ENV}" ]] && [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo ""
    echo "WARNING: No virtual environment detected!"
    echo ""
    echo "It's recommended to install DynaMune in a virtual environment."
    echo ""
    echo "Options:"
    echo "  1. Create conda environment:  conda create -n dynamune python=3.10"
    echo "                                conda activate dynamune"
    echo ""
    echo "  2. Create venv environment:   python3 -m venv venv"
    echo "                                source venv/bin/activate"
    echo ""
    read -p "Continue anyway? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Installation cancelled."
        exit 1
    fi
fi

# Upgrade pip
echo ""
echo "Upgrading pip..."
python3 -m pip install --upgrade pip

# Install package in editable mode
echo ""
echo "Installing DynaMune package..."
pip install -e .

# Verify installation
echo ""
echo "Verifying installation..."
python3 verify_install.py

echo ""
echo "======================================================================"
echo "  Installation Complete!"
echo "======================================================================"
echo ""
echo "Quick Start:"
echo "  - Start web interface:  dynamune-web"
echo "  - Run NMA analysis:     dynamune-nma --help"
echo "  - See documentation:    cat README.md"
echo ""
echo "For detailed instructions, see:"
echo "  - README.md       (comprehensive documentation)"
echo "  - QUICKSTART.md   (5-minute getting started guide)"
echo ""

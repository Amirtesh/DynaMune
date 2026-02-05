#!/usr/bin/env python3
"""
CLI wrapper for ProDy Normal Mode Analysis

This module provides a command-line entry point for the NMA analysis tool.
It wraps the original prody_nma.py script with proper path handling.
"""

import sys
import os
from pathlib import Path


def get_script_path():
    """Get the path to the original prody_nma.py script"""
    # When installed, the script is in the package root
    package_root = Path(__file__).parent.parent.parent
    script_path = package_root / "prody_nma.py"
    
    if not script_path.exists():
        raise FileNotFoundError(
            f"prody_nma.py not found at {script_path}. "
            "Please ensure DynaMune is properly installed."
        )
    
    return script_path


def main():
    """Entry point for dynamune-nma command"""
    script_path = get_script_path()
    
    # Execute the script with the original arguments
    import runpy
    sys.argv[0] = str(script_path)
    
    try:
        runpy.run_path(str(script_path), run_name="__main__")
    except SystemExit:
        # Allow the script to exit normally
        pass


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
CLI wrapper for PRS Analysis

This module provides a command-line entry point for the PRS analysis tool.
"""

import sys
from pathlib import Path


def get_script_path():
    """Get the path to the original prs_analysis.py script"""
    package_root = Path(__file__).parent.parent.parent
    script_path = package_root / "prs_analysis.py"
    
    if not script_path.exists():
        raise FileNotFoundError(
            f"prs_analysis.py not found at {script_path}. "
            "Please ensure DynaMune is properly installed."
        )
    
    return script_path


def main():
    """Entry point for dynamune-prs command"""
    script_path = get_script_path()
    
    import runpy
    sys.argv[0] = str(script_path)
    
    try:
        runpy.run_path(str(script_path), run_name="__main__")
    except SystemExit:
        pass


if __name__ == "__main__":
    main()

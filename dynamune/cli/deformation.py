#!/usr/bin/env python3
"""CLI wrapper for Deformation Analysis"""
import sys
from pathlib import Path


def main():
    script_path = Path(__file__).parent.parent.parent / "deformation_analysis.py"
    if not script_path.exists():
        raise FileNotFoundError(f"deformation_analysis.py not found at {script_path}")
    
    import runpy
    sys.argv[0] = str(script_path)
    try:
        runpy.run_path(str(script_path), run_name="__main__")
    except SystemExit:
        pass


if __name__ == "__main__":
    main()

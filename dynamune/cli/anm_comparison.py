#!/usr/bin/env python3
"""CLI wrapper for ANM Comparison Analysis"""
import sys
from pathlib import Path


def main():
    script_path = Path(__file__).parent.parent.parent / "anm_comparison.py"
    if not script_path.exists():
        raise FileNotFoundError(f"anm_comparison.py not found at {script_path}")
    
    import runpy
    sys.argv[0] = str(script_path)
    try:
        runpy.run_path(str(script_path), run_name="__main__")
    except SystemExit:
        pass


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Installation verification script for DynaMune

Run this after installation to verify everything is set up correctly.
"""

import sys
import importlib
from pathlib import Path


def check_python_version():
    """Verify Python version is 3.9+"""
    version = sys.version_info
    print(f"✓ Python {version.major}.{version.minor}.{version.micro}")
    
    if version < (3, 9):
        print(f"✗ ERROR: Python 3.9+ required, found {version.major}.{version.minor}")
        return False
    return True


def check_package_installed():
    """Check if dynamune package is importable"""
    try:
        import dynamune
        print(f"✓ DynaMune package installed (v{dynamune.__version__})")
        return True
    except ImportError as e:
        print(f"✗ ERROR: Could not import dynamune package: {e}")
        print("  Run: pip install -e .")
        return False


def check_dependencies():
    """Check if key dependencies are available"""
    required = {
        'prody': 'ProDy',
        'flask': 'Flask',
        'Bio': 'Biopython',
        'numpy': 'NumPy',
        'scipy': 'SciPy',
        'pandas': 'Pandas',
        'matplotlib': 'Matplotlib',
        'sklearn': 'scikit-learn',
    }
    
    all_ok = True
    for module, name in required.items():
        try:
            importlib.import_module(module)
            print(f"✓ {name}")
        except ImportError:
            print(f"✗ {name} not found")
            all_ok = False
    
    return all_ok


def check_cli_commands():
    """Check if CLI entry points are available"""
    import shutil
    
    commands = [
        'dynamune-web',
        'dynamune-nma',
        'dynamune-prs',
        'dynamune-pocket',
        'dynamune-deformation',
        'dynamune-domain-hinge',
        'dynamune-anm-compare',
        'dynamune-contact-stability',
        'dynamune-contact-timeline',
    ]
    
    all_ok = True
    for cmd in commands:
        if shutil.which(cmd):
            print(f"✓ {cmd}")
        else:
            print(f"✗ {cmd} not found in PATH")
            all_ok = False
    
    return all_ok


def check_files():
    """Verify essential files exist"""
    package_root = Path(__file__).parent
    
    essential_files = [
        'app.py',
        'prody_nma.py',
        'prs_analysis.py',
        'pocket_dynamics.py',
        'requirements.txt',
        'templates',
        'static',
    ]
    
    all_ok = True
    for item in essential_files:
        path = package_root / item
        if path.exists():
            print(f"✓ {item}")
        else:
            print(f"✗ {item} not found")
            all_ok = False
    
    return all_ok


def main():
    print("="*70)
    print("  DynaMune Installation Verification")
    print("="*70)
    print()
    
    print("Checking Python version...")
    python_ok = check_python_version()
    print()
    
    print("Checking DynaMune package...")
    package_ok = check_package_installed()
    print()
    
    print("Checking dependencies...")
    deps_ok = check_dependencies()
    print()
    
    print("Checking CLI commands...")
    cli_ok = check_cli_commands()
    print()
    
    print("Checking essential files...")
    files_ok = check_files()
    print()
    
    print("="*70)
    if all([python_ok, package_ok, deps_ok, cli_ok, files_ok]):
        print("✓ All checks passed! DynaMune is ready to use.")
        print()
        print("Quick start:")
        print("  Web interface:  dynamune-web")
        print("  CLI help:       dynamune-nma --help")
        print("  Documentation:  See README.md and QUICKSTART.md")
        return 0
    else:
        print("✗ Some checks failed. Please review the errors above.")
        print()
        print("Common fixes:")
        print("  1. Ensure virtual environment is activated")
        print("  2. Reinstall: pip install -e .")
        print("  3. Check dependencies: pip install -r requirements.txt")
        return 1


if __name__ == "__main__":
    sys.exit(main())

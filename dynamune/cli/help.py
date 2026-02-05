"""
DynaMune Help Command - Opens comprehensive use cases and documentation
"""
import webbrowser
from pathlib import Path

def main():
    """Open the DynaMune use cases documentation in the default web browser."""
    # Get the package root directory
    package_root = Path(__file__).parent.parent.parent
    use_cases_file = package_root / "templates" / "use_cases.html"
    
    if not use_cases_file.exists():
        print(f"Error: Documentation file not found at {use_cases_file}")
        print("Please ensure DynaMune is properly installed.")
        return 1
    
    # Open in default web browser
    file_url = f"file://{use_cases_file.absolute()}"
    print(f"Opening DynaMune documentation in your browser...")
    print(f"Location: {use_cases_file}")
    
    try:
        webbrowser.open(file_url)
        print("\n✓ Documentation opened successfully!")
        print("\nFor interactive analysis, run: dynamune-web")
    except Exception as e:
        print(f"\nError opening browser: {e}")
        print(f"\nYou can manually open: {use_cases_file}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())

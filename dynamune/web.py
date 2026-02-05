#!/usr/bin/env python3
"""
Web interface launcher for DynaMune

Starts the Flask-based local web server for interactive protein dynamics analysis.
"""

import os
import sys
from pathlib import Path
import argparse


def main():
    """Launch the DynaMune web interface"""
    parser = argparse.ArgumentParser(
        description="Launch DynaMune web interface",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Start on default port 5000
  dynamune-web

  # Start on custom port
  dynamune-web --port 8080

  # Enable debug mode (development only)
  dynamune-web --debug

  # Bind to specific host (e.g., for network access)
  dynamune-web --host 0.0.0.0 --port 8080
        """
    )
    
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host to bind to (default: 127.0.0.1 for local access only)"
    )
    
    parser.add_argument(
        "--port",
        type=int,
        default=5000,
        help="Port to run the server on (default: 5000)"
    )
    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run in debug mode (auto-reload on code changes, detailed errors)"
    )
    
    args = parser.parse_args()
    
    # Find the app.py file
    package_root = Path(__file__).parent.parent
    app_path = package_root / "app.py"
    
    if not app_path.exists():
        print(f"Error: app.py not found at {app_path}", file=sys.stderr)
        print("Please ensure DynaMune is properly installed.", file=sys.stderr)
        sys.exit(1)
    
    # Change to the package root directory so Flask can find templates and static files
    original_dir = os.getcwd()
    os.chdir(package_root)
    
    try:
        # Import and run the Flask app
        sys.path.insert(0, str(package_root))
        
        # Import the app module
        import app as app_module
        
        print(f"\n{'='*70}")
        print(f"  DynaMune Web Interface")
        print(f"{'='*70}")
        print(f"  Server starting...")
        print(f"  Host: {args.host}")
        print(f"  Port: {args.port}")
        print(f"  URL:  http://{args.host}:{args.port}")
        print(f"  Debug mode: {'ON' if args.debug else 'OFF'}")
        print(f"{'='*70}")
        print(f"\n  Press Ctrl+C to stop the server\n")
        
        # Run the Flask app
        app_module.app.run(
            host=args.host,
            port=args.port,
            debug=args.debug
        )
    
    except KeyboardInterrupt:
        print("\n\nShutting down DynaMune web server...")
        sys.exit(0)
    
    except Exception as e:
        print(f"\nError starting web server: {e}", file=sys.stderr)
        sys.exit(1)
    
    finally:
        os.chdir(original_dir)


if __name__ == "__main__":
    main()

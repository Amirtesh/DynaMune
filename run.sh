#!/bin/bash

# Automated Multi Epitope Vaccine Builder - Run Script

echo "=========================================="
echo "  Automated Multi Epitope Vaccine Builder"
echo "=========================================="
echo ""

# Check if Python is installed
if ! command -v python3 &> /dev/null
then
    echo "Error: Python 3 is not installed"
    exit 1
fi

echo "✓ Python 3 found"

# Check if pip is installed
if ! command -v pip3 &> /dev/null
then
    echo "Error: pip3 is not installed"
    exit 1
fi

echo "✓ pip3 found"

# Install dependencies
echo ""
echo "Installing dependencies..."
pip3 install -r requirements.txt

# Check if installation was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to install dependencies"
    exit 1
fi

echo ""
echo "✓ All dependencies installed"
echo ""
echo "=========================================="
echo "  Starting Flask Application"
echo "=========================================="
echo ""
echo "Server will be available at:"
echo "  Local:   http://localhost:5000"
echo "  Network: http://0.0.0.0:5000"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

# Run the Flask app
python3 app.py

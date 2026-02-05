#!/usr/bin/env python3
"""
DynaMune: Unified Platform for Protein Dynamics and Immunoinformatics
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README for long description
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

# Read requirements
requirements_file = Path(__file__).parent / "requirements.txt"
requirements = []
if requirements_file.exists():
    requirements = [
        line.strip() 
        for line in requirements_file.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]

setup(
    name="dynamune",
    version="1.0.0",
    author="DynaMune Development Team",
    author_email="",  # Add your email here
    description="Unified platform for protein dynamics analysis with integrated immunoinformatics tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Amirtesh/DynaMune",  # Update with actual repo URL
    packages=find_packages(exclude=["tests", "docs", "examples"]),
    include_package_data=True,
    package_data={
        "dynamune": [
            "templates/*.html",
            "templates/**/*.html",
            "static/**/*",
            "static/**/**/*",
            "NetBCE/**/*",
            "NetBCE/**/**/*",
        ],
    },
    install_requires=requirements,
    python_requires=">=3.9",
    entry_points={
        "console_scripts": [
            # Main web interface
            "dynamune-web=dynamune.web:main",
            "dynamune-help=dynamune.cli.help:main",
            
            # Core analysis tools
            "dynamune-nma=dynamune.cli.nma:main",
            "dynamune-prs=dynamune.cli.prs:main",
            "dynamune-pocket=dynamune.cli.pocket_dynamics:main",
            "dynamune-deformation=dynamune.cli.deformation:main",
            "dynamune-domain-hinge=dynamune.cli.domain_hinge:main",
            "dynamune-anm-compare=dynamune.cli.anm_comparison:main",
            "dynamune-contact-stability=dynamune.cli.contact_stability:main",
            "dynamune-contact-timeline=dynamune.cli.contact_timeline:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",  # Update if different
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    keywords="protein dynamics structural biology normal mode analysis immunoinformatics vaccine design",
    project_urls={
        "Bug Reports": "https://github.com/Amirtesh/DynaMune/issues",
        "Source": "https://github.com/Amirtesh/DynaMune",
    },
)

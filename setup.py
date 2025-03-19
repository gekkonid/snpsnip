#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name="snpsnip",
    version="0.1.0",
    description="Interactive VCF filtering tool",
    author="SNPSnip Team",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "snpsnip": ["templates/*", "static/*"],
    },
    install_requires=[
        "flask>=2.0.0",
        "waitress>=2.0.0",
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "scikit-learn>=1.0.0",
    ],
    entry_points={
        "console_scripts": [
            "snpsnip=snpsnip:main",
        ],
    },
    python_requires=">=3.7",
)

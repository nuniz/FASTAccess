"""
Setup script for fastaccess library.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

setup(
    name="fastaccess",
    version="0.1.0",
    author="Bioinformatics Developer",
    description="Efficient random access to subsequences in large FASTA files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/fastaccess",
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*"]),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
        "Typing :: Typed",
    ],
    python_requires=">=3.7",
    install_requires=[
        # No dependencies - stdlib only!
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
            "mypy>=1.0",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords="bioinformatics fasta genomics sequence-analysis random-access",
)

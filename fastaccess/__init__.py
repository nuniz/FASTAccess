"""
fastaccess - Efficient random access to subsequences in FASTA files.

Provides indexed, random-access retrieval of subsequences from large multi-record
FASTA files using 1-based inclusive coordinates.
"""

from .api import FastaStore

__version__ = "0.1.0"
__all__ = ["FastaStore"]

#!/usr/bin/env python3
"""
Simple benchmark for fastaccess library.
"""

import time
import tempfile
from pathlib import Path
from fastaccess.api import FastaStore


def create_test_fasta(path, num_sequences=10, seq_length=1_000_000):
    """Create a test FASTA file with multiple large sequences."""
    print(f"Creating test FASTA: {num_sequences} sequences × {seq_length:,} bp")
    with open(path, 'w') as f:
        for i in range(num_sequences):
            f.write(f'>seq{i+1} Test sequence {i+1}\n')
            # Write in 60-character lines
            seq = ('ACGT' * (seq_length // 4))[:seq_length]
            for j in range(0, len(seq), 60):
                f.write(seq[j:j+60] + '\n')
    
    file_size = Path(path).stat().st_size
    print(f"Created file: {file_size / 1024 / 1024:.1f} MB")
    return file_size


def benchmark_index_building(fasta_path):
    """Benchmark index building time."""
    print("\n--- Index Building ---")
    start = time.time()
    fa = FastaStore(fasta_path)
    elapsed = time.time() - start
    
    num_seqs = len(fa.list_sequences())
    total_bases = sum(fa.get_length(name) for name in fa.list_sequences())
    
    print(f"Time: {elapsed:.3f} seconds")
    print(f"Sequences indexed: {num_seqs}")
    print(f"Total bases: {total_bases:,}")
    print(f"Throughput: {total_bases / elapsed / 1_000_000:.1f} Mbp/s")
    
    return fa


def benchmark_fetching(fa, num_fetches=1000):
    """Benchmark subsequence fetching."""
    print("\n--- Subsequence Fetching ---")
    
    sequences = fa.list_sequences()
    
    # Small fetches (100 bp)
    print(f"\nSmall fetches ({num_fetches} × 100 bp):")
    start = time.time()
    for i in range(num_fetches):
        seq_name = sequences[i % len(sequences)]
        fa.fetch(seq_name, 1000, 1099)
    elapsed = time.time() - start
    print(f"Time: {elapsed:.3f} seconds")
    print(f"Average: {elapsed / num_fetches * 1000:.3f} ms per fetch")
    
    # Medium fetches (10 KB)
    print(f"\nMedium fetches ({num_fetches // 10} × 10 KB):")
    start = time.time()
    for i in range(num_fetches // 10):
        seq_name = sequences[i % len(sequences)]
        fa.fetch(seq_name, 5000, 14999)
    elapsed = time.time() - start
    print(f"Time: {elapsed:.3f} seconds")
    print(f"Average: {elapsed / (num_fetches // 10) * 1000:.3f} ms per fetch")
    
    # Large fetches (100 KB)
    print(f"\nLarge fetches ({num_fetches // 100} × 100 KB):")
    start = time.time()
    for i in range(num_fetches // 100):
        seq_name = sequences[i % len(sequences)]
        fa.fetch(seq_name, 10000, 109999)
    elapsed = time.time() - start
    print(f"Time: {elapsed:.3f} seconds")
    print(f"Average: {elapsed / (num_fetches // 100) * 1000:.3f} ms per fetch")


def benchmark_batch_fetching(fa, batch_size=100):
    """Benchmark batch fetching."""
    print(f"\n--- Batch Fetching ---")
    print(f"Batch size: {batch_size} queries")
    
    sequences = fa.list_sequences()
    queries = [
        (sequences[i % len(sequences)], 1000 + i * 100, 1099 + i * 100)
        for i in range(batch_size)
    ]
    
    start = time.time()
    results = fa.fetch_many(queries)
    elapsed = time.time() - start
    
    print(f"Time: {elapsed:.3f} seconds")
    print(f"Average: {elapsed / batch_size * 1000:.3f} ms per query")
    print(f"Total bases retrieved: {sum(len(s) for s in results):,}")


def main():
    """Run benchmarks."""
    print("=" * 60)
    print("fastaccess Performance Benchmark")
    print("=" * 60)
    
    # Create temporary test file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        temp_path = f.name
    
    try:
        # Create test data: 10 sequences of 1 MB each = ~10 MB file
        file_size = create_test_fasta(temp_path, num_sequences=10, seq_length=1_000_000)
        
        # Benchmark index building
        fa = benchmark_index_building(temp_path)
        
        # Benchmark fetching
        benchmark_fetching(fa, num_fetches=1000)
        
        # Benchmark batch fetching
        benchmark_batch_fetching(fa, batch_size=100)
        
        print("\n" + "=" * 60)
        print("Benchmark complete!")
        print("=" * 60)
        
    finally:
        # Clean up
        Path(temp_path).unlink(missing_ok=True)


if __name__ == "__main__":
    main()

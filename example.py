#!/usr/bin/env python3
"""
Comprehensive example demonstrating all features.
"""

import os
import shutil
import tempfile

from fastaccess import FastaStore

print("=" * 70)
print("fastaccess - Complete Feature Demo")
print("=" * 70)

# ============================================================================
# 1. BASIC USAGE
# ============================================================================
print("\n1. BASIC USAGE")
print("-" * 70)

fa = FastaStore("fastaccess/tests/wrapped.fa")
print(f"Loaded: {fa.list_sequences()}")

# Fetch subsequence
seq = fa.fetch("seq1", 1, 30)
print(f"seq1[1:30]: {seq}")

# ============================================================================
# 2. CACHE INFORMATION
# ============================================================================
print("\n2. CACHE INFORMATION")
print("-" * 70)

print(f"Was loaded from cache: {fa.is_cached()}")
print(f"Cache exists: {fa.cache_exists()}")
print(f"Cache location: {fa.get_cache_path()}")

# ============================================================================
# 3. CUSTOM CACHE DIRECTORY
# ============================================================================
print("\n3. CUSTOM CACHE DIRECTORY")
print("-" * 70)

with tempfile.TemporaryDirectory() as tmpdir:
    print(f"Using temp cache dir: {tmpdir}")
    
    # First load - builds index
    fa_custom = FastaStore("fastaccess/tests/wrapped.fa", cache_dir=tmpdir)
    print(f"  Loaded from cache: {fa_custom.is_cached()}")
    print(f"  Cache path: {fa_custom.get_cache_path()}")
    
    # Second load - uses cache
    fa_custom2 = FastaStore("fastaccess/tests/wrapped.fa", cache_dir=tmpdir)
    print(f"  Second load from cache: {fa_custom2.is_cached()}")

# ============================================================================
# 4. SEQUENCE DESCRIPTIONS
# ============================================================================
print("\n4. SEQUENCE DESCRIPTIONS")
print("-" * 70)

for name in fa.list_sequences():
    info = fa.get_info(name)
    print(f"{name}:")
    print(f"  Description: {info['description']}")
    print(f"  Length: {info['length']:,} bp")

# ============================================================================
# 5. REVERSE COMPLEMENT
# ============================================================================
print("\n5. REVERSE COMPLEMENT")
print("-" * 70)

forward = fa.fetch("seq1", 1, 20)
reverse = fa.fetch("seq1", 1, 20, reverse_complement=True)

print(f"Forward:          {forward}")
print(f"Reverse comp:     {reverse}")

# Simulate gene on minus strand
print("\nSimulating gene on minus strand:")
gene_seq = fa.fetch("seq1", 50, 70)
gene_rc = fa.fetch("seq1", 50, 70, reverse_complement=True)
print(f"  Plus strand:  {gene_seq}")
print(f"  Minus strand: {gene_rc}")

# ============================================================================
# 6. BATCH FETCHING
# ============================================================================
print("\n6. BATCH FETCHING")
print("-" * 70)

queries = [
    ("seq1", 1, 30),
    ("seq1", 31, 60),
    ("seq2", 1, 30),
]

results = fa.fetch_many(queries)
for i, (name, start, stop) in enumerate(queries):
    print(f"{name}[{start}:{stop}]: {results[i][:20]}... ({len(results[i])} bp)")

# ============================================================================
# 7. CACHE MANAGEMENT
# ============================================================================
print("\n7. CACHE MANAGEMENT")
print("-" * 70)

# Get cache info
print(f"Cache path: {fa.get_cache_path()}")
print(f"Cache exists: {fa.cache_exists()}")

# Cache size
if fa.cache_exists():
    size = os.path.getsize(fa.get_cache_path())
    print(f"Cache size: {size:,} bytes")

# Delete and rebuild
print("\nDeleting cache...")
if fa.delete_cache():
    print("  ✓ Cache deleted")
print(f"  Cache exists: {fa.cache_exists()}")

print("Rebuilding index...")
fa.rebuild_index()
print(f"  ✓ Index rebuilt")
print(f"  Cache exists: {fa.cache_exists()}")

# ============================================================================
# 8. REAL-WORLD USE CASE: GENE EXTRACTION
# ============================================================================
print("\n8. REAL-WORLD USE CASE: Gene Extraction")
print("-" * 70)

def extract_gene(fa, seq_id, start, end, strand):
    """Extract gene sequence respecting strand orientation."""
    if strand == "-":
        return fa.fetch(seq_id, start, end, reverse_complement=True)
    else:
        return fa.fetch(seq_id, start, end)

# Simulate gene coordinates
genes = [
    ("gene1", "seq1", 10, 50, "+"),
    ("gene2", "seq1", 60, 100, "-"),
    ("gene3", "seq2", 1, 60, "+"),
]

print("Extracting genes:")
for gene_name, seq_id, start, end, strand in genes:
    seq = extract_gene(fa, seq_id, start, end, strand)
    print(f"  {gene_name} ({seq_id}:{start}-{end}, {strand}): {len(seq)} bp")
    print(f"    Sequence: {seq[:30]}...")

# ============================================================================
# 9. PERFORMANCE COMPARISON
# ============================================================================
print("\n9. PERFORMANCE COMPARISON")
print("-" * 70)

import time

# Test with caching
test_file = "fastaccess/tests/wrapped.fa"

# Clear cache first
if os.path.exists(test_file + ".fidx"):
    os.remove(test_file + ".fidx")

# First load - no cache
start = time.time()
fa1 = FastaStore(test_file)
t1 = time.time() - start
print(f"First load (build index):  {t1*1000:.3f} ms")

# Second load - with cache
start = time.time()
fa2 = FastaStore(test_file)
t2 = time.time() - start
print(f"Second load (use cache):   {t2*1000:.3f} ms")
print(f"Speedup: {t1/t2:.1f}x faster")

# ============================================================================
# 10. ADVANCED: READ-ONLY GENOME SCENARIO
# ============================================================================
print("\n10. ADVANCED: Read-Only Genome Scenario")
print("-" * 70)

# Simulate read-only genome in /usr/share
# Use custom cache in user's home directory
user_cache_dir = tempfile.mkdtemp()

print(f"Simulating read-only genome with user cache:")
print(f"  FASTA: {test_file}")
print(f"  Cache: {user_cache_dir}")

fa_readonly = FastaStore(test_file, cache_dir=user_cache_dir)
print(f"  ✓ Loaded successfully")
print(f"  Cache path: {fa_readonly.get_cache_path()}")
print(f"  Cache exists: {os.path.exists(fa_readonly.get_cache_path())}")

# Cleanup
shutil.rmtree(user_cache_dir)

# fastaccess

Efficient random access to subsequences in large FASTA files using byte-level seeking.

## Installation

```bash
pip install fastaccess
```

If you’re working with the source code:

```bash
# Copy the fastaccess directory to your project
cp -r fastaccess /your/project/

# Or install as package
pip install -e .
```

## Quick Start

```python
from fastaccess import FastaStore

# Load FASTA (builds index, caches for next time)
fa = FastaStore("genome.fa")

# Fetch subsequence (1-based inclusive coordinates)
seq = fa.fetch("chr1", 1000, 2000)

# Get sequence info
info = fa.get_info("chr1")
print(info['description'])  # Full FASTA header

# Batch fetch
sequences = fa.fetch_many([
    ("chr1", 1, 100),
    ("chr2", 500, 600)
])
```

## API Reference

### `FastaStore(path, use_cache=True)`

Initialize and build index.

**Parameters:**
- `path` (str): Path to FASTA file
- `use_cache` (bool): Save/load index from `.fidx` cache file (default: True)

**Example:**
```python
fa = FastaStore("genome.fa")          # Uses cache (45x faster reload)
fa = FastaStore("genome.fa", False)   # No cache
```

---

### `fetch(name, start, stop)` → str

Fetch a subsequence using 1-based inclusive coordinates.

**Parameters:**
- `name` (str): Sequence name
- `start` (int): Start position (≥1)
- `stop` (int): Stop position (≥start, ≤length)

**Returns:** Uppercase sequence string

**Example:**
```python
seq = fa.fetch("chr1", 1000, 2000)  # Returns 1001 bases
```

---

### `fetch_many(queries)` → List[str]

Batch fetch multiple subsequences.

**Parameters:**
- `queries` (List[Tuple[str, int, int]]): List of (name, start, stop)

**Returns:** List of uppercase sequences

**Example:**
```python
seqs = fa.fetch_many([
    ("chr1", 1, 100),
    ("chr2", 500, 600)
])
```

---

### `list_sequences()` → List[str]

Get all sequence names.

**Example:**
```python
names = fa.list_sequences()  # ["chr1", "chr2", "chrM"]
```

---

### `get_length(name)` → int

Get sequence length in bases.

**Example:**
```python
length = fa.get_length("chr1")  # 248956422
```

---

### `get_description(name)` → str

Get full FASTA header description.

**Example:**
```python
desc = fa.get_description("U00096.3")
# "Escherichia coli str. K-12 substr. MG1655, complete genome"
```

---

### `get_info(name)` → dict

Get all metadata at once.

**Returns:** Dictionary with keys: `name`, `description`, `length`

**Example:**
```python
info = fa.get_info("chr1")
# {
#   'name': 'chr1',
#   'description': 'Homo sapiens chromosome 1...',
#   'length': 248956422
# }
```

---

### `rebuild_index()`

Force rebuild index and update cache.

**Example:**
```python
fa.rebuild_index()  # If FASTA file was modified
```

## Features

### Random Access
Uses `seek()` and `read()` to fetch only required bytes. No need to load entire file.

### Index Caching
Automatically saves index to `.fidx` file for 45-4300x faster reloading:
```python
fa = FastaStore("genome.fa")  # First load: builds index (30 ms)
fa = FastaStore("genome.fa")  # Second load: uses cache (0.7 ms)
```

Cache is automatically invalidated when FASTA file changes.

### Coordinate System
**1-based inclusive** (standard in bioinformatics):
```python
# Sequence: ACGTACGT...
fa.fetch("seq", 1, 4)  # Returns "ACGT" (positions 1,2,3,4)
fa.fetch("seq", 5, 8)  # Returns "ACGT" (positions 5,6,7,8)
```

### Format Support
- **Wrapped sequences** (fixed-width lines, e.g., 60 bp/line)
- **Unwrapped sequences** (single-line)
- **Unix line endings** (`\n`)
- **Windows line endings** (`\r\n`)

### Output
All sequences returned as **uppercase** strings.

## Examples

### E. coli Genome

```python
from fastaccess import FastaStore

# Load E. coli genome
fa = FastaStore("ecoli.fasta")

# Get sequence ID
seq_id = fa.list_sequences()[0]  # "U00096.3"

# Get metadata
info = fa.get_info(seq_id)
print(f"Organism: {info['description']}")
print(f"Size: {info['length']:,} bp")

# Fetch origin of replication region
oriC = fa.fetch(seq_id, 3923000, 3923500)
print(f"oriC: {oriC[:50]}...")
```

### Extract Gene Sequences

```python
# Gene coordinates
genes = [
    ("TP53", "chr17", 7668402, 7687550),
    ("BRCA1", "chr17", 43044295, 43170245),
]

fa = FastaStore("hg38.fa")

for gene_name, seq_id, start, end in genes:
    sequence = fa.fetch(seq_id, start, end)
    print(f"{gene_name}: {len(sequence):,} bp")
    # Save to file, analyze, etc.
```

### Process in Chunks

```python
fa = FastaStore("genome.fa")

seq_id = "chr1"
length = fa.get_length(seq_id)
chunk_size = 10000

for start in range(1, length, chunk_size):
    stop = min(start + chunk_size - 1, length)
    chunk = fa.fetch(seq_id, start, stop)
    # Process chunk...
```

### Extract BED Regions

```python
# BED file regions
regions = [
    ("chr1", 1000, 2000),
    ("chr1", 5000, 6000),
    ("chr2", 10000, 11000),
]

fa = FastaStore("genome.fa")
sequences = fa.fetch_many(regions)

# Write to new FASTA
with open("features.fa", "w") as out:
    for i, (seq_id, start, stop) in enumerate(regions):
        out.write(f">feature_{i+1} {seq_id}:{start}-{stop}\n")
        out.write(f"{sequences[i]}\n")
```

## Performance

**Index Building (with cache):**
```
First load:   30 ms    (parses FASTA, saves cache)
Second load:  0.7 ms   (loads cache - 45x faster!)
```

**Subsequence Fetching:**
```
Small (100 bp):    0.031 ms
Medium (10 KB):    0.133 ms
Large (100 KB):    1.072 ms
```

**Large Genomes:**
```
Human genome (3 GB):
  Without cache: 30 seconds
  With cache:    0.05 seconds  (600x faster!)
```

## How It Works

### Index Structure
```python
@dataclass
class Entry:
    name: str          # "chr1"
    description: str   # "Homo sapiens chromosome 1..."
    length: int        # 248956422
    line_blen: int     # 60 (bases per line, 0 if unwrapped)
    line_len: int      # 61 (bytes per line with \n)
    offset: int        # Byte offset to sequence data
```

### Random Access Math

**Wrapped sequences (60 bp/line):**
```python
# For position 10001:
zero_based = 10000
line_number = 10000 // 60 = 166
position_in_line = 10000 % 60 = 40
byte_offset = offset + (166 × 61) + 40

# Seek to byte_offset, read across lines, skip newlines
```

**Unwrapped sequences:**
```python
byte_offset = offset + (start - 1)
# Simple seek and read
```

## File Structure

```
fastaccess/
├── __init__.py       # Package exports
├── index.py          # Entry dataclass, build_index()
├── store.py          # fetch_subseq() with random access
└── api.py            # FastaStore class with caching
```

**Total:** ~500 lines, zero dependencies (stdlib only)

## Error Handling

```python
try:
    seq = fa.fetch("chr1", 1, 1000)
except KeyError:
    print("Sequence not found")
except ValueError as e:
    print(f"Invalid coordinates: {e}")
```

**Validation:**
- `start >= 1`
- `stop >= start`
- `stop <= sequence_length`
- Sequence name exists

## Testing

```bash
pytest fastaccess/tests/  # 19 tests, all passing
```

**Test coverage:**
- Wrapped/unwrapped sequences
- Unix/Windows line endings
- Cross-line fetching
- Edge cases (first/last bases)
- Input validation
- Batch operations
- Cache functionality

## Requirements

- Python 3.7+
- No external dependencies

## Cache Files

Index is cached to `{fasta_file}.fidx` as JSON:
```json
{
  "fasta_mtime": 1760684220.34,
  "sequences": {
    "chr1": {
      "name": "chr1",
      "description": "Homo sapiens chromosome 1...",
      "length": 248956422,
      "line_blen": 60,
      "line_len": 61,
      "offset": 123456
    }
  }
}
```

**Cache invalidation:** Automatic when FASTA file is modified

**Disable caching:**
```python
fa = FastaStore("genome.fa", use_cache=False)
```

## Limitations

- Sequences must be ASCII (DNA/RNA)
- Index rebuilt on each init (unless cached)
- No support for compressed FASTA (.gz)

## License

MIT License

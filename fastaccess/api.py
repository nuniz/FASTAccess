"""
High-level API for FASTA random access.
"""

from typing import List, Tuple
from .index import build_index, Index
from .store import fetch_subseq


class FastaStore:
    """
    High-level interface for indexed FASTA file access.
    
    Builds an in-memory index on initialization and provides methods
    for efficient random access to subsequences. Index is cached to disk
    for faster reloading.
    
    Example:
        >>> fa = FastaStore("genome.fa")
        >>> seq = fa.fetch("chr1", 10001, 20000)
        >>> print(len(seq), seq[:30])
        10000 ACGTACGTACGTACGTACGTACGTACGTAC
        
        >>> batch = fa.fetch_many([("chr1", 1, 60), ("chr2", 500, 560)])
        >>> for s in batch:
        ...     print(len(s))
        60
        61
    """
    
    def __init__(self, path: str, use_cache: bool = True, cache_dir: str = None):
        """
        Initialize the FastaStore and build or load the index.
        
        Args:
            path: Path to the FASTA file
            use_cache: If True, save/load index from cache file (.fidx)
            cache_dir: Optional directory for cache file. If None, uses same dir as FASTA.
                      Useful when FASTA directory is read-only.
        """
        import os
        
        self.path = path
        self.use_cache = use_cache
        
        # Determine cache path
        if cache_dir is not None:
            # Custom cache directory
            os.makedirs(cache_dir, exist_ok=True)
            fasta_basename = os.path.basename(path)
            self.cache_path = os.path.join(cache_dir, fasta_basename + '.fidx')
        else:
            # Default: same directory as FASTA
            self.cache_path = path + '.fidx'
        
        self._loaded_from_cache = False
        
        # Try to load from cache if available
        if use_cache and self._load_cache():
            self._loaded_from_cache = True
        else:
            # Build new index
            self.index: Index = build_index(path)
            
            # Save to cache
            if use_cache:
                self._save_cache()
    
    def _get_fasta_mtime(self) -> float:
        """Get modification time of FASTA file."""
        import os
        return os.path.getmtime(self.path)
    
    def _save_cache(self) -> None:
        """Save index to cache file."""
        import json
        import os
        
        cache_data = {
            'fasta_mtime': self._get_fasta_mtime(),
            'sequences': {}
        }
        
        for name, entry in self.index.items():
            cache_data['sequences'][name] = {
                'name': entry.name,
                'description': entry.description,
                'length': entry.length,
                'line_blen': entry.line_blen,
                'line_len': entry.line_len,
                'offset': entry.offset
            }
        
        try:
            with open(self.cache_path, 'w') as f:
                json.dump(cache_data, f, indent=2)
        except Exception:
            # Silently fail if we can't write cache
            pass
    
    def _load_cache(self) -> bool:
        """Load index from cache file if valid. Returns True if successful."""
        import json
        import os
        from .index import Entry
        
        # Check if cache file exists
        if not os.path.exists(self.cache_path):
            return False
        
        try:
            with open(self.cache_path, 'r') as f:
                cache_data = json.load(f)
            
            # Verify FASTA file hasn't been modified
            if cache_data.get('fasta_mtime') != self._get_fasta_mtime():
                return False
            
            # Rebuild index from cache
            self.index = {}
            for name, data in cache_data['sequences'].items():
                self.index[name] = Entry(
                    name=data['name'],
                    description=data['description'],
                    length=data['length'],
                    line_blen=data['line_blen'],
                    line_len=data['line_len'],
                    offset=data['offset']
                )
            
            return True
            
        except Exception:
            # If cache is corrupted or incompatible, rebuild
            return False
    
    def rebuild_index(self) -> None:
        """Force rebuild of the index and update cache."""
        self.index = build_index(self.path)
        if self.use_cache:
            self._save_cache()
    
    def is_cached(self) -> bool:
        """Check if this instance was loaded from cache."""
        return self._loaded_from_cache
    
    def cache_exists(self) -> bool:
        """Check if a cache file exists for this FASTA."""
        import os
        return os.path.exists(self.cache_path)
    
    def get_cache_path(self) -> str:
        """Get the path to the cache file."""
        return self.cache_path
    
    def delete_cache(self) -> bool:
        """
        Delete the cache file if it exists.
        
        Returns:
            True if cache was deleted, False if it didn't exist
        """
        import os
        if os.path.exists(self.cache_path):
            try:
                os.remove(self.cache_path)
                return True
            except Exception:
                return False
        return False
    
    def fetch(self, name: str, start: int, stop: int, reverse_complement: bool = False) -> str:
        """
        Fetch a single subsequence using 1-based inclusive coordinates.
        
        Args:
            name: Sequence name
            start: Start position (1-based, inclusive)
            stop: Stop position (1-based, inclusive)
            reverse_complement: If True, return reverse complement of the sequence
            
        Returns:
            Uppercase string containing the requested subsequence
            
        Raises:
            KeyError: If sequence name not found
            ValueError: If coordinates are invalid
        """
        seq = fetch_subseq(self.path, self.index, name, start, stop)
        
        if reverse_complement:
            seq = self._reverse_complement(seq)
        
        return seq
    
    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement of a DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                     'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S',
                     'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
                     'V': 'B', 'D': 'H', 'H': 'D'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    
    def fetch_many(self, queries: List[Tuple[str, int, int]]) -> List[str]:
        """
        Fetch multiple subsequences in batch.
        
        Args:
            queries: List of (name, start, stop) tuples
            
        Returns:
            List of uppercase strings, one for each query
            
        Raises:
            KeyError: If any sequence name not found
            ValueError: If any coordinates are invalid
        """
        return [self.fetch(n, s, e) for (n, s, e) in queries]
    
    def list_sequences(self) -> List[str]:
        """
        Get a list of all sequence names in the FASTA file.
        
        Returns:
            List of sequence names
        """
        return list(self.index.keys())
    
    def get_length(self, name: str) -> int:
        """
        Get the length of a sequence.
        
        Args:
            name: Sequence name
            
        Returns:
            Length of the sequence in bases
            
        Raises:
            KeyError: If sequence name not found
        """
        if name not in self.index:
            raise KeyError(f"Sequence '{name}' not found in index")
        return self.index[name].length
    
    def get_description(self, name: str) -> str:
        """
        Get the description of a sequence.
        
        Args:
            name: Sequence name
            
        Returns:
            Description text (everything after the name in the header)
            
        Raises:
            KeyError: If sequence name not found
        """
        if name not in self.index:
            raise KeyError(f"Sequence '{name}' not found in index")
        return self.index[name].description
    
    def get_info(self, name: str) -> dict:
        """
        Get all metadata for a sequence.
        
        Args:
            name: Sequence name
            
        Returns:
            Dictionary with name, description, and length
            
        Raises:
            KeyError: If sequence name not found
        """
        if name not in self.index:
            raise KeyError(f"Sequence '{name}' not found in index")
        entry = self.index[name]
        return {
            'name': entry.name,
            'description': entry.description,
            'length': entry.length
        }

"""
Random access subsequence retrieval from indexed FASTA files.
"""

from typing import Dict
from .index import Index, Entry


def fetch_subseq(path: str, index: Index, name: str, start: int, stop: int) -> str:
    """
    Fetch a subsequence using 1-based inclusive coordinates.
    
    For wrapped sequences:
        - Calculate which line the start position is on
        - Calculate byte offset within that line
        - Read across lines, skipping newline bytes
    
    For unwrapped sequences:
        - Simple seek and read operation
    
    Args:
        path: Path to the FASTA file
        index: Pre-built index dictionary
        name: Sequence name
        start: Start position (1-based, inclusive)
        stop: Stop position (1-based, inclusive)
        
    Returns:
        Uppercase string containing the requested subsequence
        
    Raises:
        KeyError: If sequence name not found
        ValueError: If coordinates are invalid
    """
    # Validate sequence name exists
    if name not in index:
        raise KeyError(f"Sequence '{name}' not found in index")
    
    entry = index[name]
    
    # Validate coordinates
    if start < 1:
        raise ValueError(f"Start position must be >= 1, got {start}")
    if stop < start:
        raise ValueError(f"Stop position must be >= start, got start={start}, stop={stop}")
    if stop > entry.length:
        raise ValueError(
            f"Stop position {stop} exceeds sequence length {entry.length} for '{name}'"
        )
    
    # Calculate number of bases to read
    num_bases = stop - start + 1
    
    with open(path, 'rb') as f:
        if entry.line_blen == 0:
            # Unwrapped sequence - simple seek and read
            byte_pos = entry.offset + (start - 1)
            f.seek(byte_pos)
            data = f.read(num_bases)
            return data.decode('ascii', errors='ignore').upper()
        else:
            # Wrapped sequence - need to skip newlines
            # Calculate starting position
            # start-1 because we're converting from 1-based to 0-based
            zero_based_start = start - 1
            
            # Which line (0-indexed) does our start position fall on?
            blocks_before = zero_based_start // entry.line_blen
            
            # Position within that line
            within_line = zero_based_start % entry.line_blen
            
            # Byte offset to seek to
            byte_pos = entry.offset + blocks_before * entry.line_len + within_line
            
            f.seek(byte_pos)
            
            # Read bases, skipping newlines
            result = []
            bases_read = 0
            bases_remaining_in_line = entry.line_blen - within_line
            
            while bases_read < num_bases:
                # How many bases to read from current line?
                to_read = min(bases_remaining_in_line, num_bases - bases_read)
                
                # Read the bases
                chunk = f.read(to_read)
                if not chunk:
                    break
                
                result.append(chunk)
                bases_read += to_read
                
                # If we need more bases, skip the newline and continue
                if bases_read < num_bases:
                    # Skip newline bytes
                    newline_size = entry.line_len - entry.line_blen
                    f.read(newline_size)
                    
                    # Next line has full line_blen bases available
                    bases_remaining_in_line = entry.line_blen
            
            # Combine all chunks and return uppercase
            combined = b''.join(result)
            return combined.decode('ascii', errors='ignore').upper()

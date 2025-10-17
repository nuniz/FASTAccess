"""
Index building for FASTA files.

Parses FASTA files in a single pass to build an in-memory index that enables
efficient random access to subsequences.
"""

from dataclasses import dataclass
from typing import Dict


@dataclass
class Entry:
    """Index entry for a single FASTA record."""
    name: str          # Sequence name (from header)
    description: str   # Full description text after name
    length: int        # Total number of bases in the sequence
    line_blen: int     # Bases per full line; 0 if unwrapped (single line)
    line_len: int      # Bytes per line including newline(s); 0 if unwrapped
    offset: int        # Byte offset where sequence data starts


Index = Dict[str, Entry]


def build_index(path: str) -> Index:
    """
    Build an in-memory index of all FASTA records in the file.
    
    Single pass over FASTA:
      - Parse headers (>name ...)
      - Detect newline style (\n or \r\n)
      - Detect first sequence line width (line_blen) and full line bytes (line_len)
      - Accumulate total base length (strip newline bytes)
      - Record byte offset where sequence starts
    
    Args:
        path: Path to the FASTA file
        
    Returns:
        Dictionary mapping sequence name to Entry with index information
    """
    index: Index = {}
    
    with open(path, 'rb') as f:
        current_name = None
        current_description = ""
        current_offset = 0
        current_length = 0
        current_line_blen = 0
        current_line_len = 0
        newline_size = 0
        first_seq_line = True
        
        while True:
            line_start = f.tell()
            line = f.readline()
            
            if not line:
                # End of file - save last entry if exists
                if current_name is not None:
                    index[current_name] = Entry(
                        name=current_name,
                        description=current_description,
                        length=current_length,
                        line_blen=current_line_blen,
                        line_len=current_line_len,
                        offset=current_offset
                    )
                break
            
            # Check if this is a header line
            if line.startswith(b'>'):
                # Save previous entry if exists
                if current_name is not None:
                    index[current_name] = Entry(
                        name=current_name,
                        description=current_description,
                        length=current_length,
                        line_blen=current_line_blen,
                        line_len=current_line_len,
                        offset=current_offset
                    )
                
                # Parse new header - extract name and description
                header = line[1:].decode('ascii', errors='ignore').strip()
                parts = header.split(maxsplit=1)  # Split into name and rest
                current_name = parts[0] if parts else header
                current_description = parts[1] if len(parts) > 1 else ""
                
                # Next byte is where sequence starts
                current_offset = f.tell()
                current_length = 0
                current_line_blen = 0
                current_line_len = 0
                first_seq_line = True
                
            else:
                # This is a sequence line
                # Detect newline style on first sequence line
                if first_seq_line:
                    # Determine newline size
                    if line.endswith(b'\r\n'):
                        newline_size = 2
                    elif line.endswith(b'\n'):
                        newline_size = 1
                    else:
                        # Last line with no newline, or unwrapped
                        newline_size = 0
                    
                    # Strip newline to get bases
                    bases = line.rstrip(b'\r\n')
                    bases_count = len(bases)
                    
                    # Check if wrapped or unwrapped
                    # Peek at next byte to see if there's another line
                    next_pos = f.tell()
                    next_byte = f.read(1)
                    
                    if next_byte and not next_byte.startswith(b'>'):
                        # There's another sequence line, so this is wrapped
                        current_line_blen = bases_count
                        current_line_len = bases_count + newline_size
                    else:
                        # Unwrapped (single line) or last line
                        current_line_blen = 0
                        current_line_len = 0
                    
                    # Seek back
                    f.seek(next_pos)
                    
                    current_length += bases_count
                    first_seq_line = False
                else:
                    # Subsequent sequence lines
                    bases = line.rstrip(b'\r\n')
                    current_length += len(bases)
    
    return index

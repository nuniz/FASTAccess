"""
Comprehensive test suite for fastaccess library.
"""

import os
# Add parent directory to path for imports
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from fastaccess.api import FastaStore

# Test fixtures paths
TEST_DIR = Path(__file__).parent
WRAPPED_FA = TEST_DIR / "wrapped.fa"
UNWRAPPED_FA = TEST_DIR / "unwrapped.fa"
WINDOWS_FA = TEST_DIR / "windows.fa"


class TestWrappedFasta:
    """Tests for wrapped FASTA files (60 bp per line)."""
    
    def setup_method(self):
        """Initialize FastaStore for each test."""
        self.store = FastaStore(str(WRAPPED_FA))
    
    def test_index_building(self):
        """Test that index is built correctly."""
        assert "seq1" in self.store.index
        assert "seq2" in self.store.index
        
        entry1 = self.store.index["seq1"]
        assert entry1.length == 180  # 3 lines * 60 bp
        assert entry1.line_blen == 60
        assert entry1.line_len == 61  # 60 bp + 1 newline
        
        entry2 = self.store.index["seq2"]
        assert entry2.length == 240  # 4 lines * 60 bp
    
    def test_subsequence_within_one_line(self):
        """Test fetching subsequence within a single line."""
        # First 30 bases of seq1 (within first line)
        seq = self.store.fetch("seq1", 1, 30)
        assert len(seq) == 30
        assert seq == "ACGTACGTACGTACGTACGTACGTACGTAC"
        
        # Middle 20 bases of first line
        seq = self.store.fetch("seq1", 21, 40)
        assert len(seq) == 20
        assert seq == "ACGTACGTACGTACGTACGT"
    
    def test_subsequence_across_line_breaks(self):
        """Test fetching subsequence that spans multiple lines."""
        # Fetch across first and second line
        seq = self.store.fetch("seq1", 50, 70)
        assert len(seq) == 21
        # Should concatenate without newlines
        # Line 1 is ACGT*15, pos 50-60 is "CGTACGTACGT" (11 bases from line 1)
        # Line 2 is TGCA*15, pos 61-70 is "TGCATGCATG" (10 bases from line 2)
        expected = "CGTACGTACGT" + "TGCATGCATG"
        assert seq == expected.upper()
        
        # Fetch across all three lines
        seq = self.store.fetch("seq1", 55, 125)
        assert len(seq) == 71
        assert "\n" not in seq
        assert "\r" not in seq
    
    def test_first_and_last_base(self):
        """Test fetching first and last bases of a record."""
        # First base
        seq = self.store.fetch("seq1", 1, 1)
        assert seq == "A"
        
        # Last base (position 180)
        seq = self.store.fetch("seq1", 180, 180)
        assert seq == "C"
        
        # Full sequence
        seq = self.store.fetch("seq1", 1, 180)
        assert len(seq) == 180
    
    def test_second_sequence(self):
        """Test fetching from second sequence in file."""
        seq = self.store.fetch("seq2", 1, 60)
        assert len(seq) == 60
        assert seq == "A" * 60
        
        seq = self.store.fetch("seq2", 61, 120)
        assert seq == "T" * 60
        
        seq = self.store.fetch("seq2", 121, 180)
        assert seq == "G" * 60
        
        seq = self.store.fetch("seq2", 181, 240)
        assert seq == "C" * 60


class TestUnwrappedFasta:
    """Tests for unwrapped FASTA files (single line sequences)."""
    
    def setup_method(self):
        """Initialize FastaStore for each test."""
        self.store = FastaStore(str(UNWRAPPED_FA))
    
    def test_index_building(self):
        """Test that unwrapped sequences are indexed correctly."""
        entry = self.store.index["single1"]
        assert entry.line_blen == 0  # Unwrapped
        assert entry.line_len == 0
        assert entry.length == 180
    
    def test_fetch_from_unwrapped(self):
        """Test fetching from unwrapped sequences."""
        seq = self.store.fetch("single1", 1, 30)
        assert len(seq) == 30
        assert seq == "ACGTACGTACGTACGTACGTACGTACGTAC"
        
        # Fetch from middle
        seq = self.store.fetch("single1", 100, 120)
        assert len(seq) == 21
        
        # Full sequence
        seq = self.store.fetch("single1", 1, 180)
        assert len(seq) == 180


class TestWindowsFasta:
    """Tests for FASTA files with Windows-style CRLF newlines."""
    
    def setup_method(self):
        """Initialize FastaStore for each test."""
        self.store = FastaStore(str(WINDOWS_FA))
    
    def test_crlf_handling(self):
        """Test that CRLF newlines are handled correctly."""
        entry = self.store.index["winseq"]
        # Should detect 2-byte newlines
        assert entry.line_len == 62  # 60 bp + 2 bytes (\r\n)
        assert entry.length == 120  # 2 lines * 60 bp
    
    def test_fetch_with_crlf(self):
        """Test that fetching works correctly with CRLF."""
        # Fetch within first line
        seq = self.store.fetch("winseq", 1, 30)
        assert len(seq) == 30
        assert seq == "ACGTACGTACGTACGTACGTACGTACGTAC"
        
        # Fetch across CRLF boundary
        seq = self.store.fetch("winseq", 50, 70)
        assert len(seq) == 21
        assert "\r" not in seq
        assert "\n" not in seq


class TestInputValidation:
    """Tests for input validation and error handling."""
    
    def setup_method(self):
        """Initialize FastaStore for each test."""
        self.store = FastaStore(str(WRAPPED_FA))
    
    def test_sequence_not_found(self):
        """Test that KeyError is raised for non-existent sequence."""
        with pytest.raises(KeyError, match="not found"):
            self.store.fetch("nonexistent", 1, 10)
    
    def test_start_less_than_one(self):
        """Test that ValueError is raised for start < 1."""
        with pytest.raises(ValueError, match="Start position must be >= 1"):
            self.store.fetch("seq1", 0, 10)
        
        with pytest.raises(ValueError, match="Start position must be >= 1"):
            self.store.fetch("seq1", -5, 10)
    
    def test_stop_less_than_start(self):
        """Test that ValueError is raised for stop < start."""
        with pytest.raises(ValueError, match="Stop position must be >= start"):
            self.store.fetch("seq1", 10, 5)
    
    def test_stop_exceeds_length(self):
        """Test that ValueError is raised for stop > length."""
        with pytest.raises(ValueError, match="exceeds sequence length"):
            self.store.fetch("seq1", 1, 200)  # seq1 length is 180
    
    def test_valid_edge_cases(self):
        """Test valid edge cases."""
        # start == stop (single base)
        seq = self.store.fetch("seq1", 50, 50)
        assert len(seq) == 1
        
        # Full sequence
        seq = self.store.fetch("seq1", 1, 180)
        assert len(seq) == 180


class TestBatchFetching:
    """Tests for batch fetching operations."""
    
    def setup_method(self):
        """Initialize FastaStore for each test."""
        self.store = FastaStore(str(WRAPPED_FA))
    
    def test_fetch_many(self):
        """Test fetching multiple subsequences."""
        queries = [
            ("seq1", 1, 60),
            ("seq1", 61, 120),
            ("seq2", 1, 60),
            ("seq2", 61, 120),
        ]
        
        results = self.store.fetch_many(queries)
        assert len(results) == 4
        assert all(len(r) == 60 for r in results)
        
        # Verify content
        assert results[0] == "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        assert results[2] == "A" * 60
        assert results[3] == "T" * 60
    
    def test_fetch_many_mixed_sizes(self):
        """Test fetching subsequences of different sizes."""
        queries = [
            ("seq1", 1, 10),
            ("seq2", 50, 100),
            ("seq1", 170, 180),
        ]
        
        results = self.store.fetch_many(queries)
        assert len(results) == 3
        assert len(results[0]) == 10
        assert len(results[1]) == 51
        assert len(results[2]) == 11


class TestAPIHelpers:
    """Tests for additional API helper methods."""
    
    def setup_method(self):
        """Initialize FastaStore for each test."""
        self.store = FastaStore(str(WRAPPED_FA))
    
    def test_list_sequences(self):
        """Test listing all sequence names."""
        names = self.store.list_sequences()
        assert "seq1" in names
        assert "seq2" in names
        assert len(names) == 2
    
    def test_get_length(self):
        """Test getting sequence lengths."""
        assert self.store.get_length("seq1") == 180
        assert self.store.get_length("seq2") == 240
        
        with pytest.raises(KeyError):
            self.store.get_length("nonexistent")


class TestUppercaseOutput:
    """Tests to ensure output is always uppercase."""
    
    def test_lowercase_input(self):
        """Test that lowercase input is converted to uppercase."""
        # Create a temp file with lowercase bases
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(">lowercase\n")
            f.write("acgtacgt" * 10 + "\n")
            temp_path = f.name
        
        try:
            store = FastaStore(temp_path)
            seq = store.fetch("lowercase", 1, 10)
            assert seq == "ACGTACGTAC"
            assert seq.isupper()
        finally:
            os.unlink(temp_path)


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])


import pytest
from src.smith_waterman import SmithWaterman


def test_identical_sequences():
    """Perfect match should give maximum score."""
    sw = SmithWaterman(match_score=2)
    _, _, score = sw.align("ACGT", "ACGT")
    assert score == 8  # 4 matches × 2

def test_no_overlap():
    """Completely different sequences should give score of 0."""
    sw = SmithWaterman()
    _, _, score = sw.align("AAAA", "CCCC")
    assert score == 0

def test_partial_match():
    """Should find the best local match."""
    sw = SmithWaterman()
    a1, a2, score = sw.align("TTTACGT", "ACGT")
    assert "ACGT" in a1 or score > 0

def test_gap_penalty():
    """Gap penalties should reduce score."""
    sw_tight = SmithWaterman(gap_open=-5)
    sw_loose = SmithWaterman(gap_open=-1)
    seq1, seq2 = "ACGT", "A-GT"
    # Tight gaps should score lower
    _, _, s1 = sw_tight.align("ACGT", "AGT")
    _, _, s2 = sw_loose.align("ACGT", "AGT")
    assert s2 >= s1



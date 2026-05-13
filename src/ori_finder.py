
"""
Ori-Finder: Replication Origin Detection
Uses GC-skew analysis and frequent k-mer detection (DnaA boxes)
"""

import numpy as np
from collections import Counter
from typing import List, Tuple, Dict


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))


def gc_skew(genome: str, window: int = 1000, step: int = 100) -> np.ndarray:
    """
    Calculate GC skew = (G - C) / (G + C) in sliding windows.
    The minimum of cumulative GC skew approximates oriC.

    Args:
        genome: DNA sequence string
        window: Window size in bp
        step: Step size in bp

    Returns:
        Array of (position, skew) tuples
    """
    positions, skews = [], []
    for i in range(0, len(genome) - window, step):
        chunk = genome[i:i+window]
        g = chunk.count('G')
        c = chunk.count('C')
        skew = (g - c) / (g + c) if (g + c) > 0 else 0
        positions.append(i + window // 2)
        skews.append(skew)
    return np.array(positions), np.array(skews)


def cumulative_skew(genome: str) -> np.ndarray:
    """
    Calculate position-by-position cumulative GC skew.
    The minimum point predicts oriC.
    """
    skew = [0]
    for base in genome.upper():
        if base == 'G':   skew.append(skew[-1] + 1)
        elif base == 'C': skew.append(skew[-1] - 1)
        else:             skew.append(skew[-1])
    return np.array(skew)


def predict_oric(genome: str) -> int:
    """Return the predicted position of oriC based on min cumulative GC skew."""
    skew = cumulative_skew(genome)
    return int(np.argmin(skew))


def kmer_frequencies(seq: str, k: int) -> Dict[str, int]:
    """Count all k-mers and their reverse complements in a sequence."""
    counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k].upper()
        rc   = reverse_complement(kmer)
        canonical = min(kmer, rc)  # canonical form
        counts[canonical] += 1
    return dict(counts)


def find_dnaa_boxes(genome: str, oric_pos: int,
                    window: int = 500, k: int = 9) -> List[Tuple[str, int]]:
    """
    Search for frequent 9-mers near predicted oriC — these are candidate DnaA boxes.
    E. coli DnaA box consensus: TTATCCACA

    Returns top 10 most frequent k-mers near oriC.
    """
    start = max(0, oric_pos - window)
    end   = min(len(genome), oric_pos + window)
    region = genome[start:end]

    freqs = kmer_frequencies(region, k)
    sorted_kmers = sorted(freqs.items(), key=lambda x: -x[1])
    return sorted_kmers[:10]


def hamming_distance(s1: str, s2: str) -> int:
    """Count mismatches between two equal-length strings."""
    assert len(s1) == len(s2), 'Sequences must be same length'
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def find_approximate_pattern(genome: str, pattern: str, d: int) -> List[int]:
    """
    Find all positions where pattern occurs with at most d mismatches.
    Used to find DnaA binding sites.
    """
    k = len(pattern)
    positions = []
    for i in range(len(genome) - k + 1):
        if hamming_distance(genome[i:i+k].upper(), pattern.upper()) <= d:
            positions.append(i)
    return positions

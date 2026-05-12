"""
Smith-Waterman Local Sequence Alignment
Implementation from scratch — O(mn) time and space
"""

import numpy as np
from typing import Tuple


class SmithWaterman:
    """
    Smith-Waterman local alignment algorithm.
    
    Parameters:
        match_score (int): Score for matching characters. Default: 2
        mismatch_score (int): Penalty for mismatches. Default: -1
        gap_open (int): Gap opening penalty. Default: -2
    """

    def __init__(self, match_score=2, mismatch_score=-1, gap_open=-2):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_open = gap_open

    def score(self, a: str, b: str) -> int:
        """Return match or mismatch score for two characters."""
        return self.match_score if a == b else self.mismatch_score

    def align(self, seq1: str, seq2: str) -> Tuple[str, str, float]:
        """
        Perform local alignment between seq1 and seq2.

        Returns:
            aligned_seq1: Aligned version of seq1 (with gaps as '-')
            aligned_seq2: Aligned version of seq2 (with gaps as '-')
            score: Alignment score
        """
        m, n = len(seq1), len(seq2)

        # Step 1: Initialize scoring matrix
        H = np.zeros((m + 1, n + 1), dtype=float)
        traceback = np.zeros((m + 1, n + 1), dtype=int)
        # traceback: 0=stop, 1=diagonal, 2=up, 3=left

        best_score = 0
        best_pos = (0, 0)

        # Step 2: Fill the matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                diag = H[i-1][j-1] + self.score(seq1[i-1], seq2[j-1])
                up   = H[i-1][j]   + self.gap_open
                left = H[i][j-1]   + self.gap_open

                H[i][j] = max(0, diag, up, left)

                if H[i][j] == 0:          traceback[i][j] = 0
                elif H[i][j] == diag:     traceback[i][j] = 1
                elif H[i][j] == up:       traceback[i][j] = 2
                else:                     traceback[i][j] = 3

                if H[i][j] >= best_score:
                    best_score = H[i][j]
                    best_pos = (i, j)

        # Step 3: Traceback
        aligned1, aligned2 = [], []
        i, j = best_pos

        while traceback[i][j] != 0:
            if traceback[i][j] == 1:
                aligned1.append(seq1[i-1])
                aligned2.append(seq2[j-1])
                i -= 1; j -= 1
            elif traceback[i][j] == 2:
                aligned1.append(seq1[i-1])
                aligned2.append('-')
                i -= 1
            else:
                aligned1.append('-')
                aligned2.append(seq2[j-1])
                j -= 1

        return ''.join(reversed(aligned1)), ''.join(reversed(aligned2)), best_score

    def format_alignment(self, seq1: str, seq2: str) -> str:
        """Pretty-print alignment result."""
        a1, a2, score = self.align(seq1, seq2)
        match_line = ''.join('|' if c1 == c2 else ' '
                             for c1, c2 in zip(a1, a2))
        return (
            f'Score: {score}\n'
            f'Seq1: {a1}\n'
            f'      {match_line}\n'
            f'Seq2: {a2}'
        )

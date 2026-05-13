
import sys
import os
sys.path.insert(0, os.path.abspath('.'))
#!/usr/bin/env python3
"""Run ori-finder on E. profundum genome."""

import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from Bio import SeqIO
from src.ori_finder import cumulative_skew, predict_oric, find_dnaa_boxes


def load_genome(fasta_path: str) -> str:
    """Load first sequence from FASTA file."""
    record = next(SeqIO.parse(fasta_path, 'fasta'))
    return str(record.seq)


def main(fasta_path: str):
    print(f'Loading genome from {fasta_path}...')
    genome = load_genome(fasta_path)
    print(f'Genome length: {len(genome):,} bp')

    # 1. Predict oriC
    oric = predict_oric(genome)
    print(f'\nPredicted oriC position: {oric:,} bp')

    # 2. Find DnaA boxes near oriC
    dnaa = find_dnaa_boxes(genome, oric)
    print('\nTop candidate DnaA box sequences near oriC:')
    for kmer, count in dnaa[:5]:
        print(f'  {kmer}  (count: {count})')

    # 3. Plot cumulative GC skew
    skew = cumulative_skew(genome)
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(skew, color='steelblue', linewidth=0.8, alpha=0.9)
    ax.axvline(oric, color='red', linestyle='--', linewidth=2,
               label=f'Predicted oriC: {oric:,} bp')
    ax.set_xlabel('Genomic Position (bp)', fontsize=12)
    ax.set_ylabel('Cumulative GC Skew', fontsize=12)
    ax.set_title('E. profundum SS9 — Cumulative GC Skew & Predicted oriC',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('data/results/oric_gcskew.png', dpi=150, bbox_inches='tight')
    print('\nPlot saved to data/results/oric_gcskew.png')
    plt.show()


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'data/sample/e_profundum_genome.fasta'
    main(path)



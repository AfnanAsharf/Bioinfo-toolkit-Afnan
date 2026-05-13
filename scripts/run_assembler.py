
import sys
import os
sys.path.insert(0, os.path.abspath('.'))

from Bio import SeqIO
from src.debruijn_assembler import DeBruijnAssembler
from src.smith_waterman import SmithWaterman

def load_reads(path):
    reads = []
    for rec in SeqIO.parse(path, 'fasta'):
        reads.append(str(rec.seq))
    return reads

path = sys.argv[1] if len(sys.argv) > 1 else 'data/sample/sequence_16S.fasta'
print(f'Loading reads from {path}...')
reads = load_reads(path)
print(f'Loaded {len(reads)} reads')

assembler = DeBruijnAssembler(k=21)
contigs = assembler.get_contigs(reads)
contigs.sort(key=len, reverse=True)
print(f'Generated {len(contigs)} contigs')
if contigs:
    print(f'Longest contig: {len(contigs[0])} bp')
    print(f'Sequence: {contigs[0]}')
    # Compare assembled contig to original using Smith-Waterman
print('\nComparing assembled contig to original sequence using Smith-Waterman...')
original = reads[0]
assembled = contigs[0]
sw = SmithWaterman()
a1, a2, score = sw.align(assembled[:200], original[:200])
print(f'Alignment score: {score}')
print(f'Assembled: {a1[:60]}')
print(f'Original:  {a2[:60]}')
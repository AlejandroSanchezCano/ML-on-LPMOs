import os
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from Bio import Align
from parse_structures import AlphaFoldStructure
from variables import FASTA, PLOTS, SEQUENCE_CLUSTERING, CAZY_EXPANDED

# List all sequences
with open(f'{FASTA}/His1.fasta', 'r') as file:
    sequences = file.read().split('\n\n')[:-1]
    sequences = [fasta.split('\n')[1] for fasta in sequences]
    ids = [fasta[1:].split('\n')[0] for fasta in sequences]

# Configure pairwise alignemnt
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = Align.substitution_matrices.load("PAM70")
aligner.open_gap_score = -4
aligner.extend_gap_score = -0.5

# Compute pairwise alignment -> O(nlog(n))
n_seqs = len(sequences)
matrix = np.zeros([n_seqs, n_seqs])
for i in tqdm(range(n_seqs)):
    for j in range(i, n_seqs):
        score = aligner.align(sequences[i], sequences[j]).score
        matrix[i][j] = int(score)

# Triangular matrix -> full matrix
X = matrix + matrix.T - np.diag(np.diag(matrix))
np.save(f'{SEQUENCE_CLUSTERING}/his1_matrix', X)


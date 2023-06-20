import os
import numpy as np
import seaborn as sns
from tqdm import tqdm
from variables import FASTA, PLOTS
from Bio import Align
from parse_structures import AlphaFoldStructure


# List all sequences
with open(f'{FASTA}/His1.fasta', 'r') as file:
    sequences = file.read().split('\n\n')[:-1]
    sequences = [fasta.split('\n')[1] for fasta in sequences]

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

# Plot
heatmap = sns.clustermap(X)
heatmap.savefig(f'{PLOTS}/sequence_heatmap.png', transparent=True)

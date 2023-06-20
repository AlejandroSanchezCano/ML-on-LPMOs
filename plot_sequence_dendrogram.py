import os
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from Bio import Align
from parse_structures import AlphaFoldStructure
from variables import FASTA, PLOTS, SEQUENCE_CLUSTERING, CAZY_EXPANDED

print('started')

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

# Assign color to family
df = pd.read_pickle(f'{CAZY_EXPANDED}/AA_filtered')
families = ['AA0', 'AA9', 'AA10', 'AA11', 'AA13', 'AA14', 'AA15', 'AA16', 'AA17']
families_colors = dict(zip(families, sns.color_palette()))
joined_df = pd.merge(
    left = df,
    right = pd.DataFrame({'UniProt' : ids}),
    on = 'UniProt'
).drop_duplicates('UniProt')

colors = joined_df['Family'].map(families_colors).to_numpy()

# Plot
heatmap = sns.clustermap(X, row_cluster = False, column_colors = colors)
heatmap.savefig(f'{PLOTS}/sequence_heatmap.png', transparent=True)

import pickle
import numpy as np
import pandas as pd
import seaborn as sns
from variables import SEQUENCE_CLUSTERING, CAZY_DATA, FASTA, PLOTS

# Import matrix
matrix = np.load(f'{SEQUENCE_CLUSTERING}/his1_matrix.npy')

# Get AlphaFold IDs 
with open(f'{FASTA}/His1.fasta', 'r') as file:
    sequences = file.read().split('\n\n')[:-1]
    ids = [fasta[1:].split('\n')[0] for fasta in sequences]

# AlphaFold IDs -> corresponding families
with open(f'{CAZY_DATA}/uniprot_family.pkl', 'rb') as dic:
    dic = pickle.load(dic)


# Assign color to each family
families = ['AA0', 'AA9', 'AA10', 'AA11', 'AA13', 'AA14', 'AA15', 'AA16', 'AA17']
families_colors = dict(zip(families, sns.color_palette()))
colors = [families_colors[dic[id]] for id in ids]

# Plot
heatmap = sns.clustermap(matrix, col_colors = colors, row_colors = colors)
heatmap.savefig(f'{PLOTS}/sequence_heatmap.png', transparent=True)

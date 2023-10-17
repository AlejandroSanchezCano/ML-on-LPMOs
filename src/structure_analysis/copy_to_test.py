import os
import pickle
import argparse
import numpy as np
import pandas as pd
import geometricus as gm
from sklearn.decomposition import PCA
#from variables import CAZY_DATA, STRUCTURE_CLUSTERING
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
import umap as UMAP

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

'''
# Argument handling
parser = argparse.ArgumentParser(
                    prog = 'geometricus_embedder',
                    description = 'Runs Geometricus on a set of protein structures',
                    epilog = 'See ya'
                    )

parser.add_argument('--proteins', '-p', 
                    type = str, 
                    default = None,
                    help = 'Set of proteins')
parser.add_argument('--resolution', '-r', 
                    type = float, 
                    default = 1.,
                    help = 'Coarse-grainness of the shapemers')
parser.add_argument('--threads', '-t', 
                    type = int, 
                    default = 4,
                    help = 'Number of threads use in multithreading')
parser.add_argument('--kmer', '-k',
                    type = int,
                    default = 8,
                    help = 'k-mer-based fragmentation'
                    )
parser.add_argument('--radius', '-o',
                    type = int,
                    default = 5,
                    help = 'radius-based fragmentation'
                    )
'''

# Get argument values
#args = parser.parse_args()
structures_folder = '/home/sanch046/lustre/ML-on-LPMOs/data/AF_core' #args.proteins
resolution = 2 #args.resolution
n_threads = 4 #args.threads
kmer = 6 #args.kmer
radius = 0 #args.radius

# Handle kmer/radius fragmentation logic
if kmer and radius:
    split = [
        gm.SplitInfo(gm.SplitType.KMER, kmer), 
        gm.SplitInfo(gm.SplitType.RADIUS, radius)
        ]
elif kmer and not radius:
    split = [gm.SplitInfo(gm.SplitType.KMER, kmer)]
elif not kmer and radius:
    split = [gm.SplitInfo(gm.SplitType.RADIUS, radius)]
else:
    raise ValueError('Fragmentation types cannot be both null')

# Input proteins
proteins = [structure for structure in os.listdir(structures_folder)]

# Run Geometricus
invariants, _ = gm.get_invariants_for_structures(
    structures_folder, 
    n_threads = n_threads,
    split_infos = split,
    moment_types = ["O_3", "O_4", "O_5", "F"]
    )

shapemer_class = gm.Geometricus.from_invariants(
    invariants, 
    protein_keys = proteins, 
    resolution = resolution
    )

shapemer_count_matrix = shapemer_class.get_count_matrix()

# Normalizationfor protein length
shapemer_sum = np.sum(shapemer_count_matrix, axis = 1)
normalized_matrix = shapemer_count_matrix/shapemer_sum[:, None]
normalized_matrix *= 100

proteins = [protein.replace('.pdb', '') for protein in shapemer_class.protein_keys]
shapemers = shapemer_class.shapemer_keys
matrix = pd.DataFrame(normalized_matrix, index = proteins, columns = shapemers)
matrix.to_pickle('matrix.pkl')

print(matrix)
# Save normalized embedding
#file_name = structures_folder.split('/')[-1]
#np.save(f'{STRUCTURE_CLUSTERING}/{file_name}_resolution{resolution}.npy', normalized_matrix)

'''
# UMAP
umap = UMAP.UMAP(metric = "cosine").fit_transform(normalized_matrix)

# PCA
pca = PCA().fit_transform(normalized_matrix)

# PCA scree plot
PC_values = np.arange(pca.n_components_) + 1
PC_explained_variance_ratio = pca.explained_variance_ratio_

#PLOT FOR NOW BUT THIS WILL BE REMOVED SINCE IT WILL BE IN THE PLOTLY   
fig = plt.plot(PC_values, PC_explained_variance_ratio, 'ro-', linewidth=2)
plt.title('Scree Plot')
plt.xlabel('Principal Component')
plt.ylabel('Proportion of Variance Explained')
plt.show()


# Add family
with open(f'{CAZY_DATA}/uniprot_family.pkl', 'rb') as dic:
    dic = pickle.load(dic)

# Add annotations to embeddings
protein_names = [name.split('.')[0] for name in shapemer_class.protein_keys]
umap = pd.DataFrame(umap[:, 0:2], index = protein_names)
pca = pd.DataFrame(pca[:, 0:2], index = protein_names)
families = [dic[protein] for protein in protein_names]
umap['Family'] = families
pca['Family'] = families

# Disregard AA0
umap = umap[umap['Family'] != 'AA0']
pca = pca[pca['Family'] != 'AA0']

##########ADDITION
import seaborn as sns

'''

'''
# One Hot Encode families into colors
dummies = pd.get_dummies(umap['Family'])

# Translate family to color
family_color = {
    'AA9' : '#1f77b4',
    'AA10' : '#ff7f0e',
    'AA11' : '#2ca02c',
    'AA13' : '#d62728',
    'AA14' : '#9467bd',
    'AA15' : '#8c564b',
    'AA16' : '#e377c2',
    'AA17' : '#bcbd22'
}
for column in dummies:
    dummies[column] = dummies[column].replace(1, family_color[column])
    dummies[column] = dummies[column].replace(0, 'grey')

import plotly.graph_objects as go

# Scatter plots
scatter_umap = go.Scatter(
    x = umap[0],
    y = umap[1],
    name = "Default UMAP",
    mode = 'markers',
    visible = True,
    hovertemplate = umap.index.to_list(),
    marker = dict(color = ['grey']*len(dummies))
)

scatter_pca = go.Scatter(
    x = pca[0],
    y = pca[1],
    name = "Default PCA",
    mode = 'markers',
    visible = False,
    hovertemplate = pca.index.to_list(),
    marker=dict(color = ['grey']*len(dummies))
)

fig = go.FigureWidget([scatter_umap, scatter_pca])

button_layer_1_height = 1.3
fig.update_layout(
    updatemenus=[
        dict(
            buttons=list([
                
                dict(
                    label = "None",
                    args = [{"marker": {"color":['grey']*len(dummies)}}],
                    method = "restyle"
                ),  
                
                dict(
                    label = "AA9",
                    args = [{"marker": {"color": dummies['AA9']}}],
                    method = "restyle"
                ),    

                dict(
                    label = "AA10",
                    args = [{"marker": {"color": dummies['AA10']}}],
                    method = "restyle"
                ),   
                
                dict(
                    label = "AA11",
                    args = [{"marker": {"color": dummies['AA11']}}],
                    method = "restyle"
                ),   

                dict(
                    label = "AA13",
                    args = [{"marker": {"color": dummies['AA13']}}],
                    method = "restyle"
                ), 
                
                dict(
                    label = "AA14",
                    args = [{"marker": {"color": dummies['AA14']}}],
                    method = "restyle"
                ),   

                dict(
                    label = "AA15",
                    args = [{"marker": {"color": dummies['AA15']}}],
                    method = "restyle"
                ),   

                dict(
                    label = "AA16",
                    args = [{"marker": {"color": dummies['AA16']}}],
                    method = "restyle"
                ),                  
                
                dict(
                    label = "AA17",
                    args = [{"marker": {"color": dummies['AA17']}}],
                    method = "restyle"
                ),   

            ]), 
            direction="down",
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.11,
            xanchor="left",
            y=button_layer_1_height,
            yanchor="top"
        ),
        dict(
            buttons=list([
                dict(label="UMAP",
                     method="update",
                     args=[{"visible": [True, False]},
                           {"title": "UMAP"}]),
                dict(label="PCA",
                     method="update",
                     args=[{"visible": [False, True]},
                           {"title": "PCA"}]
                ),
            ]),
            direction="down",
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0,
            xanchor="left",
            y=button_layer_1_height,
            yanchor="top"
        )
    ]
)

# create our callback function
def selection_fn(trace, points, selector):
    # Give shapemer and protein info to normalized_matrix
    df = pd.DataFrame(
        data = normalized_matrix, 
        columns = shapemer_class.shapemer_keys,
        index = [key.split('.')[0] for key in shapemer_class.protein_keys]
        )
    
    # Calculate most frequent shapemer/protein
    point_indeces = points.point_inds
    if point_indeces:
        print(point_indeces)
        points = df.iloc[point_indeces].to_numpy()
        point_names = set(df.iloc[point_indeces].index)
        index = np.unravel_index(np.argmax(points, axis=None), points.shape)
        most_frequent_shapemer = df.columns[index[1]]
        print(most_frequent_shapemer)

        # Map shapemer to residues
        selectedname_residue = dict()
        name_residue = dict(shapemer_class.map_shapemer_to_residues(most_frequent_shapemer))
        for name in point_names:
            if name + '.pdb' in name_residue:
                selectedname_residue[name] = name_residue[name + '.pdb']
        print(selectedname_residue)

fig.data[0].on_selection(selection_fn)
fig.data[1].on_selection(selection_fn)

fig
'''
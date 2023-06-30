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
family_color = {
    'AA0' : '#e377c2',
    'AA9' : '#1f77b4',
    'AA10' : '#ff7f0e',
    'AA11' : '#2ca02c',
    'AA13' : '#d62728',
    'AA14' : '#9467bd',
    'AA15' : '#8c564b',
    'AA16' : '#e377c2',
    'AA17' : '#bcbd22'
}
colors = [family_color[dic[id]] for id in ids]

# Plot
heatmap = sns.clustermap(matrix, col_colors = colors, row_colors = colors)
heatmap.savefig(f'{PLOTS}/sequence_heatmap.png', transparent=True)

# Average hierarchical clustering
from scipy.cluster import hierarchy
import fastcluster
Z = fastcluster.linkage(matrix, 'average')

def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick

# Get Newick format
tree = hierarchy.to_tree(Z, False)
newick = get_newick(tree, tree.dist, leaf_names = ids)
print(newick)
with open(f'{SEQUENCE_CLUSTERING}/his1_tree.txt', 'w') as tree_file:
    tree_file.write(newick)

#  Write annotation file
with open(f'{SEQUENCE_CLUSTERING}/his1_annotations.txt', 'w') as annotations_file:
    settings = '''
DATASET_COLORSTRIP
#lines starting with a hash are comments and ignored during parsing
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throught this file (except in the SEPARATOR line, which uses space).

#SEPARATOR TAB
SEPARATOR SPACE
#SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL color_strip1

#dataset color (can be changed later)
COLOR #000

#optional settings

#all other optional settings can be set or changed later in the web interface (under 'Datasets' tab)
COLOR_BRANCHES 0
#maximum width
STRIP_WIDTH 25

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
MARGIN 0

#border width; if set above 0, a black border of specified width (in pixels) will be drawn around the color strip 
BORDER_WIDTH 1
BORDER_COLOR #000

#show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
SHOW_INTERNAL 0

#In colored strip charts, each ID is associated to a color. Color can be specified in hexadecimal, RGB or RGBA notation
#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#Actual data follows after the "DATA" keyword
DATA
'''
    annotations_file.write(settings)

    for id in ids:
        annotations_file.write(f'{id} {family_color[dic[id]]}\n')
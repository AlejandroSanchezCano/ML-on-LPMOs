import pickle
import numpy as np
import pandas as pd
import seaborn as sns
from variables import SEQUENCE_CLUSTERING, CAZY_DATA, FASTA, PLOTS

# Import matrix
matrix = np.load(f'{SEQUENCE_CLUSTERING}/core_matrix.npy')

# Get AlphaFold IDs 
with open(f'{FASTA}/core.fasta', 'r') as file:
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
heatmap.savefig(f'{PLOTS}/core_sequence_heatmap.png', transparent=True)

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
with open(f'{SEQUENCE_CLUSTERING}/core_tree.txt', 'w') as tree_file:
    tree_file.write(newick)

#  Write annotation color file
with open(f'{SEQUENCE_CLUSTERING}/core_annotations_color.txt', 'w') as annotations_file:
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


#  Write annotation binary file
from variables import CAZY_EXPANDED, SUPERVISED

df_labels = pd.read_pickle(f'{SUPERVISED}/regioselectivity_substratespecificity_db.pkl')

def f(regio, substrate):
    if regio == 'C1':
        c1 = '1'
        c4 = '0'
    elif regio == 'C4':
        c1 = '0'
        c4 = '1'
    else:
        c1, c4 = '1', '1'
    
    if 'chitin' in substrate:
        s = '0'
    elif 'cellulose' in substrate:
        s = '1'
    else:
        s = '-1'
    return c1, c4, s

l = df_labels.apply(lambda x: f(x.Regioselectivity, x['Substrate specificity']), axis = 1)
d = dict(zip(df_labels['UniProt'], l))

with open(f'{SEQUENCE_CLUSTERING}/core_annotations_binary.txt', 'w') as annotations_file:
    settings = '''DATASET_BINARY
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throught this file (except in the SEPARATOR line, which uses space).

#SEPARATOR TAB
#SEPARATOR SPACE
SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL,binary_data

#dataset color (can be changed later)
COLOR,#ff0000

#Binary datasets can contain one or more values for each node. Each value will be represented by a symbol (defined in FIELD_SHAPES) with corresponding color and label (from FIELD_COLORS and FIELD_LABELS). Possible values (defined under DATA below) for each node are 1 (filled shapes), 0 (empty shapes) and -1 (completely ommited).

#define colors for each individual field column (if not defined all symbols will use the main dataset color, defined in COLOR)
#shapes for each field column; possible choices are
#1: rectangle 
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
FIELD_LABELS,rl0,rl1,rl2
FIELD_COLORS,#000,#000,#66f
FIELD_SHAPES,4,5,3

#all other optional settings can be set or changed later in the web interface (under 'Datasets' tab)

#show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
SHOW_INTERNAL,1

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
MARGIN,0

#symbol height factor; Default symbol height will be slightly less than the available space between leaves, but you can set a multiplication factor here to increase/decrease it (values from 0 to 1 will decrease it, values above 1 will increase it)
HEIGHT_FACTOR,1

#increase/decrease the spacing between individual levels, when there is more than one binary level defined 
SYMBOL_SPACING,10

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#Actual data follows after the "DATA" keyword
DATA

#Example dataset with 4 columns (circle, left triangle, right triangle and rectangle):
#FIELD_SHAPES,2,4,5,1
#FIELD_LABELS,f1,f2,f3,f4
#FIELD_COLORS,#ff0000,#00ff00,#ffff00,#0000ff
#DATA
#node 9606 will have a filled circle, empty left triangle, nothing in the 3rd column and an empty rectangle
#9606,1,0,-1,0
'''


    annotations_file.write(settings)


    for id in ids:
        if id  in d:
            labels = ','.join(d[id])
            annotations_file.write(f'{id},{labels}\n')
        else:
            annotations_file.write(f'{id},-1,-1,-1\n')
        

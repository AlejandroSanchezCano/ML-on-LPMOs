import pickle
import argparse
import pandas as pd
import fastcluster
from annotations import Annotation
from scipy.cluster import hierarchy
from variables import CAZY_DATA, SUPERVISED

# Argument handling
parser = argparse.ArgumentParser(
                    prog = 'plot_dendrogram',
                    description = 'Generates Newick tree anda annotations file from a matrix',
                    epilog = 'See ya'
                    )

parser.add_argument('--matrix', '-m', 
                    type = str, 
                    default = None,
                    help = 'Input matrix path')
parser.add_argument('--output', '-o', 
                    type = str, 
                    default = None,
                    help = 'Output path')

# Get argument values
args = parser.parse_args()
matrix = args.matrix
output = args.output

# Import matrix
matrix = pd.read_pickle(matrix)

# Get AlphaFold IDs 
ids = matrix.index.to_list()

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
protein_color = {id:family_color[dic[id]] for id in ids}

# Supervised labels
def format_to_itol(regio, substrate):
    if regio == 'C1':
        c1 = '1'
        c4 = '0'
    elif regio == 'C4':
        c1 = '0'
        c4 = '1'
    else:
        c1, c4 = '1', '1'
    
    if 'cellulose' in substrate:
        s = '1'
    elif 'chitin' in substrate:
        s = '0'
    else:
        s = '-1'
    return c1, c4, s

labels = pd.read_pickle(f'{SUPERVISED}/regioselectivity_substratespecificity_db.pkl')


labels_itol = labels.apply(lambda x: format_to_itol(x['Regioselectivity'], x['Substrate specificity']), axis = 1)
binary_labels = dict(zip(labels['UniProt'], labels_itol))

# Average hierarchical clustering
Z = fastcluster.linkage(matrix, 'average')
tree = hierarchy.to_tree(Z, False)

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
    

# Tree Newick format
newick = get_newick(tree, tree.dist, leaf_names = ids)
with open(f'{output}_tree.txt', 'w') as tree_file:
    tree_file.write(newick)

# Annotations
annotation = Annotation()
colorstrip = annotation.colorstrip(protein_color)
print(colorstrip)
binary = annotation.binary(binary_labels)
print(binary)
with open(f'{output}_colorstrip.txt', 'w') as color_file:
    color_file.write(colorstrip)
with open(f'{output}_binary.txt', 'w') as binary_file:
    binary_file.write(binary)
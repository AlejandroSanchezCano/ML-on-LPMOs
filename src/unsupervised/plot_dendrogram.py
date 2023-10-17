'''
Plot dendrogram

Takes a FASTA file containing LPMO entries. Generates the tree file in newick 
format. Additionally, generates the colostrip and binary annotation files to be 
plotted in iTOL.

Functions
---------
handle_arguments
get_newick
make_tree
make_colorstrip
to_binary
make_binary
'''


import scipy
import pickle
import argparse
import numpy as np
import fastcluster
import pandas as pd
from .annotations import Annotation
from ..config.config_parser import config

def handle_arguments() -> argparse.Namespace:
    '''
    Handles the arguments passed via command line.

    Return
    ------
    args : argparse.Namespace
        Arguments
    '''
    
    parser = argparse.ArgumentParser(
        prog = 'plot_dendrogram',
        description = 'Generates tree from similarity matrix and'\
            ' additional annotation files for iTOL',
        epilog = 'See ya'
        )
    
    # Add arguments
    parser.add_argument(
        '-i', '--input_matrix',
        type = str,
        help = 'Input similarity matrix path'
        )
    
    parser.add_argument(
        '-o', '--output_dir',
        type = str,
        help = 'Output path where to store the generated files'
        )
    
    parser.add_argument(
        '-c', '--cut_tree',
        action = 'store_true',
        help = 'Whether to cut the hierarchical tree.'
        )

    # Arguments from Namespace object
    return parser.parse_args()

def get_newick( 
        node : scipy.cluster.hierarchy.ClusterNode, 
        parent_dist : float, 
        leaf_names : list, 
        newick : None = ''
        ) -> str:
    """
    Uses recursion to convert scipy.cluster.hierarchy.to_tree() output to 
    Newick format.

    Parameters
    ----------
    node : scipy.cluster.hierarchy.ClusterNode
        Output of scipy.cluster.hierarchy.to_tree()
    
    parent_dist : float
        Output of scipy.cluster.hierarchy.to_tree().dist
    
    leaf_names : list
        List of leaf names
    
    newick : None
        Leave empty, this variable is used in recursion

    Returns
    -------
    newick : str
        Tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (
            leaf_names[node.id], 
            parent_dist - node.dist, 
            newick
            )
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(
            node.get_left(), 
            node.dist, 
            leaf_names, 
            newick = newick
            )
        newick = get_newick(
            node.get_right(), 
            node.dist, 
            leaf_names, 
            newick = ",%s" % (newick)
            )
        newick = "(%s" % (newick)
        
        return newick

def make_linkage(args : argparse.Namespace) -> np.ndarray:
    '''
    Converts a similarity matrix to a scipy linkage matrix.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments.
    
    Returns
    -------
    Z : np.ndarray
        Scipy linkage matrix
    '''

    # Unpickle similarity matrix
    similarity_matrix = pd.read_pickle(args.input_matrix)

    # Average hierarchical clustering
    Z = fastcluster.linkage(similarity_matrix, 'average')

    return Z

def make_tree(args : argparse.Namespace, Z : np.ndarray):
    '''
    Converts a linkage matrix to a scipy linkage matrix.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments.

    Z : np.ndarray
        Scipy linkage matrix
    
    Returns
    -------
    None
    '''

    # Unpickle similarity matrix
    similarity_matrix = pd.read_pickle(args.input_matrix)
    
    # Make tree
    tree = scipy.cluster.hierarchy.to_tree(Z, False)

    # Convert tree to newick format
    newick = get_newick(
        node = tree, 
        parent_dist = tree.dist, 
        leaf_names = similarity_matrix.index.to_list()
        )

    # Store newick tree file
    with open(f'{args.output_dir}/tree.txt', 'w') as tree_file:
        tree_file.write(newick)

def cut_tree(args : argparse.Namespace, Z : np.ndarray):
    '''
    Cut the hierarchical tree to retrieve clusters.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments.

    Z : np.ndarray
        Scipy linkage matrix.
    
    Returns
    -------
    None.
    '''

    # Unpickle similarity matrix to get UniProt IDs
    similarity_matrix = pd.read_pickle(args.input_matrix)
    uniprots = similarity_matrix.columns

    # Cut tree
    cutree = scipy.cluster.hierarchy.cut_tree(Z)

    # Initialize matrix with UniProt IDs per height (column) 
    # and custer number (row)
    matrix = [[[] for _ in range(len(cutree))] for _ in range(len(cutree))]

    for height, agglomeration in enumerate(cutree.T):
        for uniprot, cluster in zip(uniprots, agglomeration):
            matrix[cluster][height].append(uniprot) 
    
    # Store matrix
    pd.to_pickle(pd.DataFrame(matrix), f'{args.output_dir}/agglomeration.pkl')


def make_colorstrip(args : argparse.Namespace):
    '''
    Makes colorstrip annotation file for iTOL plotting of a tree.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments.
    
    Returns
    -------
    None.
    '''

    # Unpickle UniProt ID : family fictionary
    with open(f'{config["data"]}/uniprot2family.pkl', 'rb') as handle:
        uniprot2family = pickle.load(handle)

    # Import family : color dictioonary
    family2color = config['family2color']

    # Build UniProt ID : color dictionary
    uniprot2colors = {uniprot:family2color[family] \
              for uniprot, family in uniprot2family.items()}
    
    # Make colorstrip annotation file
    colorstrip = Annotation().colorstrip(uniprot2colors)

    # Store colostrip annotation file
    with open(f'{args.output_dir}/colorstrip.txt', 'w') as annotation_file:
        annotation_file.write(colorstrip)

def to_binary(
        regio : str, 
        substrate : tuple[str]
        ) -> tuple[str, str, str]:

    '''
    Translate the labels of the entries used for supervised learning to binary
    labels. With regioselectivity, the label 1 is assigned to C1 and C1/C4 
    LPMOs and the label 0 is assigned to C4 and C1/C4 LPMOs. With substrate 
    specificity, the label 1 is assigned to cellulose-degrading LPMOs, the
    label 0 to chitin-degrading LPMOs and the label -1 to those which degrade
    neither cellulose not chitin.

    Parameters
    ----------
    regio : str
        Regioselectivity.

    substrate : Tuple[str]
        Substrate specificity

    Returns
    -------
    c1 : str
        '1' if the entry cleaves at C1 or C1/C4. '0' if the entry cleaves at C4
        or C1/C4.

    c4 : str
        '0' if the entry cleaves at C1 or C1/C4. '1' if the entry cleaves at C4
        or C1/C4.

    ss : str
        '1' if the entry cleaves cellulose. '0' if the entry cleaves chitin.
        '-1' if neither.

    '''
    # Regioselectivity
    if regio == 'C1':
        c1 = '1'
        c4 = '0'
    elif regio == 'C4':
        c1 = '0'
        c4 = '1'
    else:
        c1, c4 = '1', '1'
    
    # Substrate specificity
    if 'chitin' in substrate:
        ss = '0'
    elif 'cellulose' in substrate:
        ss = '1'
    else:
        ss = '-1'
    return c1, c4, ss

def make_binary(args):
    '''
    Makes binary annotation file for iTOL plotting of a tree.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments.
    
    Returns
    -------
    None.
    '''
    
    # Read supervised data frame´
    supervised = pd.read_pickle(f'{config["supervised"]}/labels.pkl')

    # Make binary labels
    transform_labels_to_binary = lambda df: to_binary(
        df['Regioselectivity'], 
        df['Substrate specificity']
        )
    binary_labels = supervised.apply(transform_labels_to_binary, axis = 1)
    uniprot2binary = dict(zip(supervised['UniProt'], binary_labels))

    # Make binary annotation file
    binary = Annotation().binary(uniprot2binary)

    # Store binary annotation file
    with open(f'{args.output_dir}/binary.txt', 'w') as annotation_file:
        annotation_file.write(binary)

def main():
    '''Program flow.'''
    
    # Arguments
    args = handle_arguments()

    # Make tree file
    Z = make_linkage(args)
    make_tree(args, Z)

    # Cut tree
    if args.cut_tree:
        cut_tree(args, Z)
    

    # Make annotation files
    make_colorstrip(args)
    make_binary(args)

if __name__ == '__main__':
    main()
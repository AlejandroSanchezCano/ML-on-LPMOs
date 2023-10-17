import pickle
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
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
    
    # Arguments from Namespace object
    return parser.parse_args()

def pca(args):

    # Import matrix
    matrix = pd.read_pickle(args.input_matrix)

    # Import uniprot2family
    with open(f'{config["data"]}/uniprot2family.pkl', 'rb') as handle:
        uniprot2family = pickle.load(handle)

    # Family index
    new_index = [uniprot2family[uniprot] for uniprot in matrix.index]

    # PCA
    model = PCA(n_components = 2, random_state = 1999)
    space = model.fit_transform(matrix)
    scores = pd.DataFrame(
        data = space, 
        columns = [0, 1], 
        index = new_index
        )

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
    print(scores)
    
    # Plot PCA
    for family in family_color: 
        scores_family = scores.loc[family]
        plt.scatter(
            scores_family[0], 
            scores_family[1], 
            color = family_color[family], 
            label = family,
            s = 2 
            )

    plt.xlabel(f'PC1 {round(model.explained_variance_ratio_[0]*100, 2)}%')
    plt.ylabel(f'PC2 {round(model.explained_variance_ratio_[1]*100, 2)}%')
    plt.legend()
    plt.show()

def main():
    '''Program flow.'''
    
    # Arguments
    args = handle_arguments()

    # PCA
    pca(args)


if __name__ == '__main__':
    main()
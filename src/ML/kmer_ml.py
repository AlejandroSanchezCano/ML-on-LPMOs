'''
Sequence-based ML 

Run machine learning algorithms and perform feature selection with kmers as
features.

Functions
---------
handle_arguments
kmer
'''

import re
import argparse
import numpy as np
import pandas as pd
from collections import Counter
from ..config.config_parser import config
from sklearn.linear_model import LogisticRegression
from ..filtering.parse_pdb import AlphaFoldStructure
from sklearn.model_selection import LeaveOneOut, cross_val_score

def handle_arguments() -> argparse.Namespace:
    '''
    Handles the arguments passed via command line

    Return
    ------
    args : argparse.Namespace
        Arguments
    '''

    # Define parser 
    parser = argparse.ArgumentParser(prog = 'Machine Learning')

    # Add arguments
    parser.add_argument(
        '--problem', '-p',
        type = str,
        default = 'C1C4',
        choices = ['C1', 'C4', 'C1C4', 'SS'],
        help = 'Problem'
    )

    parser.add_argument(
        '--kmer_length', '-k',
        type = int,
        default = 3,
        help = 'K-mer length'
    )

    return parser.parse_args()

def kmer_ML(args : argparse.Namespace, k : int):

    # Import supervised database of characterized LPMOs
    db = pd.read_pickle(f'{config["supervised"]}/{args.problem}.pkl')

    # Initialize proto-embedding
    uniprot2kmercounts = {}

    # Get sequences
    for uniprot in db['UniProt']:
        structure = AlphaFoldStructure(
            f"{config['AF_labels_core']}/{args.problem}/{uniprot}.pdb"
            )
        seq = structure.seq
        
        # Count kmers
        kmer_counts = Counter(seq[i : i + k] for i in range(len(seq)+1-k))
        
        # Add to dict
        uniprot2kmercounts[structure.id] = kmer_counts

    # Make embedding
    embedding = pd.DataFrame.from_dict(
        data = uniprot2kmercounts, 
        orient = 'index'
        ).fillna(0)
    
    # Use same order in database and in emebdding data frames
    sorter = embedding.index.to_list()
    db['UniProt'] = db['UniProt'].astype('category')
    db['UniProt'] = db['UniProt'].cat.set_categories(sorter)
    db = db.sort_values(by = 'UniProt')

    # Logistic regression with LOOCV
    cv = LeaveOneOut()
    logreg = LogisticRegression(
        max_iter = 1000,
        random_state = 1999
        )
    scores = cross_val_score(
        estimator = logreg, 
        X = embedding * 100, 
        y = db['y'], 
        #scoring = 'neg_root_mean_squared_error',
        cv = cv, 
        n_jobs = -1
        )
    accuracy = np.mean(scores)

    # Run ML algorithm with all the data
    X_train = embedding * 100
    y_train = db['y']
    logreg = LogisticRegression(
        max_iter = 1000,
        random_state = 1999,
        )
    logreg.fit(X_train, y_train)

    # Feature importance logistic regression
    feature_importance = pd.DataFrame({
        'Feature': X_train.columns, 
        'Importance': logreg.coef_[0]
        })
    feature_importance = feature_importance.sort_values(
        by = 'Importance',
        ascending = False,
        key = lambda series : np.abs(series),
        ignore_index = True
        )
    report_name = f'{config["ML"]}/kmer/{args.kmer_length}_{args.problem}.txt'
    feature_importance.iloc[0:10].to_csv(report_name, sep = ' ')
    
    # Plot feature importance
    plot = feature_importance.plot(
        x = 'Feature', 
        y = 'Importance', 
        kind = 'barh',
        )
    fig = plot.get_figure()
    fig.savefig(
        f"{config['plots']}/{args.kmer_length}_{args.problem}.png", 
        transparent = True
        )

    # Map kmer to sequence
    for rank in range(0, 4):
        kmer = feature_importance.at[rank, 'Feature']
        importance = feature_importance.at[rank, 'Importance']
        proteins = [uniprot for uniprot, kmers in uniprot2kmercounts.items() if kmer in kmers]

        with open(report_name, 'a') as report:
            report.write(f'\nAccuracy = {accuracy}\n')
            report.write(f"\nkmer: {kmer}, rank: {rank + 1}, importance: {importance}, Number of proteins with k-mer: {len(proteins)}\n\n")

            for pdb_name in proteins:
                start = db[db['UniProt'] == pdb_name]['Core'].values
                start = start[0][0]
                counts = uniprot2kmercounts[pdb_name][kmer]
                seq = AlphaFoldStructure(f"{config['AF_labels_core']}/{args.problem}/{pdb_name}.pdb").seq
                indeces = [m.start() + start for m in re.finditer(kmer, seq)]
                report.write(f'{pdb_name} --> {indeces}\n')

            report.write(f'\n{"-"*80}\n')

    return accuracy

def main():
    '''Program flow.'''
    # Argument handling
    args = handle_arguments()

    # Perform ML on kmer counts
    embedding = kmer_ML(args = args, k = args.kmer_length) * 100
    print('Accuracy: ', embedding)

if __name__ == '__main__':
    main()

import shap
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn import metrics
import matplotlib.pyplot as plt
from ..config.config_parser import config
from .embedder import GeometricusExperiment
from sklearn.tree import DecisionTreeClassifier
from lazypredict.Supervised import LazyClassifier
from sklearn.linear_model import LogisticRegression
from ..filtering.parse_pdb import AlphaFoldStructure
from sklearn.model_selection import train_test_split
from sklearn.model_selection import LeaveOneOut, cross_val_score


def handle_arguments() -> argparse.Namespace:
    '''
    Handles the arguments passed via command line

    Return
    ------
    args : argparse.Namespace
        Arguments
    '''

    # Define parser and subparser
    parser = argparse.ArgumentParser(prog='Machine Learning')
    subparsers = parser.add_subparsers(
        help = 'sub-command help', 
        required = True,
        dest = 'subparser'
        )

    # Subparser optimization
    parser_optimization = subparsers.add_parser(
        'optimization', 
        prog = 'Hyperparameter optimization',
        description = 'Runs Geometricus on a set of protein structures with\
            different parameters to determine the optimal set',
        help='optimization help'
        )

    # Arguments subparser optimization
    parser_optimization.add_argument(
        '--kmers', '-k',
        type = np.ndarray,
        default = np.array([6]),#np.arange(0,12),
        help = 'k-mer-based fragmentation'
        )
    parser_optimization.add_argument(
        '--radiuses', '-o',
        type = np.ndarray,
        default = np.array([0]),#np.arange(0,12),
        help = 'Radius-based fragmentation'
        )

    parser_optimization.add_argument(
        '--resolutions', '-r', 
        type = np.ndarray, 
        default = np.array([2]),#np.array([0.001, 0.01, 0.1, 0.5, 1, 2, 4, 6, 8, 10]),
        help = 'Coarse-grainness of the shapemers'
        )
    
    parser_optimization.add_argument(
        '--problem', '-p',
        type = str,
        default = 'C1C4',
        choices = ['C1', 'C4', 'C1C4', 'SS'],
        help = 'Problem'
    )

    # Subparser feature selection
    parser_feature_selection = subparsers.add_parser(
        'selection', 
        prog = 'Hyperparameter optimization',
        description = 'Runs Geometricus on a set of protein structures with\
            the best parameter set and run feature selection',
        help='Selection help'
        )
    
    # Arguments feature selection
    parser_feature_selection.add_argument(
        '--kmer', '-k',
        type = int,
        default = 0,
        help = 'k-mer-based fragmentation'
        )
    parser_feature_selection.add_argument(
        '--radius', '-o',
        type = int,
        default = 6,
        help = 'Radius-based fragmentation'
        )

    parser_feature_selection.add_argument(
        '--resolution', '-r', 
        type = float, 
        default = 2.,
        help = 'Coarse-grainness of the shapemers'
        )
    
    parser_feature_selection.add_argument(
        '--problem', '-p',
        type = str,
        default = 'C1C4',
        choices = ['C1', 'C4', 'C1C4', 'SS'],
        help = 'Problem'
    )
    
    # Arguments from Namespace object
    return parser.parse_args()

def optimize(args):
    # Import supervised database of characterized LPMOs
    db = pd.read_pickle(f'{config["supervised"]}/{args.problem}.pkl')
    
    # Initialize output variable
    accuracies = {}

    # Iterate over the hyperparameters
    pbar = tqdm(total = args.kmers.size * args.radiuses.size * args.resolutions.size)
    for kmer in args.kmers:
        for radius in args.radiuses:
            
            # Deal with case when both splits are 0
            if kmer == 0 and radius == 0:
                pbar.update(args.resolutions.size)
                continue

            # Calculate invariants from GeometricusExperiment
            experiment = GeometricusExperiment(
                input_dir = f"{config['AF_labels_core']}/{args.problem}"
                )
            invariants = experiment.run_invariants(kmer, radius)

            # The bottleneck in a GeometricusExperiment run is calculating the 
            # invariants, so the resultion gets applied after that
            for resolution in args.resolutions:
                print(kmer, radius, resolution)
                
                # Calculate matrix
                experiment.run_embedding(invariants, resolution)
                #print(experiment.embedding(1)[(6, 0, 5, 12)].loc['Q9S296'])
                print(experiment.embedding(100))

                # Accuracy = 0 when only one shapemer is detected
                if experiment.embedding(100).shape[1] == 1:
                    accuracies[(kmer, radius, resolution)] = 0
                
                # Use same order in database and in emebdding data frames
                sorter = experiment.embedding(100).index.to_list()
                db['UniProt'] = db['UniProt'].astype('category')
                db['UniProt'] = db['UniProt'].cat.set_categories(sorter)
                db = db.sort_values(by = 'UniProt')

                # Train test
                #split = train_test_split(
                #    experiment.matrix(100),
                #    db['y'], 
                #    test_size = 0.25, 
                #    random_state = 1999
                #    )
                #X_train, X_test, y_train, y_test = split

                # Determine best ML algorithm with lazy predict
                #clf = LazyClassifier(
                #    verbose = 0,
                #    ignore_warnings = True, 
                #    custom_metric = None
                #    )
                #models, predictions = clf.fit(X_train, X_test, y_train, y_test)
                #print(models)
                #print(predictions)

                # Logistic regression with LOOCV
                cv = LeaveOneOut()
                logreg = LogisticRegression(
                    max_iter = 1000,
                    random_state = 1999
                    )
                scores = cross_val_score(
                    estimator = logreg, 
                    X = experiment.embedding(100), 
                    y = db['y'], 
                    #scoring = 'neg_root_mean_squared_error',
                    cv = cv, 
                    n_jobs = -1
                    )
                
                b_accuracy = metrics.balanced_accuracy_score(db['y'], scores)
                accuracy = metrics.accuracy_score(db['y'], scores)
                accuracies[(kmer, radius, resolution)] = (accuracy * 100, b_accuracy * 100)

                print(db['y'])
                print(scores)
                print(accuracy)

                # Update progress bar
                pbar.update()

    # Close progress bar
    pbar.close()
    
    # Sort accuracies dict
    parameter_options = dict(sorted(
        accuracies.items(), 
        key = lambda item: item[1],
        reverse = True
        ))
    
    return parameter_options

def feature_importance(args):
    # Import supervised database of characterized LPMOs
    db = pd.read_pickle(f'{config["supervised"]}/{args.problem}.pkl')

    # Run Geometricus
    experiment = GeometricusExperiment(
                input_dir = f"{config['AF_labels_core']}/{args.problem}"
                )
    invariants = experiment.run_invariants(args.kmer, args.radius)
    experiment.run_embedding(invariants, args.resolution)
    X = experiment.embedding(100)
    shapemer_class = experiment.shapemers

    # Use same order in database and in emebdding data frames
    sorter = experiment.embedding(100).index.to_list()
    db['UniProt'] = db['UniProt'].astype('category')
    db['UniProt'] = db['UniProt'].cat.set_categories(sorter)
    db = db.sort_values(by = 'UniProt')

    # Run ML algorithm with all the data
    X_train = X
    y_train = db['y']
    logreg = LogisticRegression(
        max_iter = 1000,
        random_state = 1999,
        )
    logreg.fit(X_train, y_train)
    
    # Feature importance logistic regression
    feature_importance = pd.DataFrame({'Feature': X.columns, 'Importance': logreg.coef_[0]})
    feature_importance = feature_importance.sort_values(
        by = 'Importance',
        ascending = False,
        key = lambda series : np.abs(series),
        ignore_index = True
        )
    report_name = f'{config["ML"]}/shapemer/{args.problem}_{args.kmer}_{args.radius}_{args.resolution}.txt'
    feature_importance.iloc[0:10].to_csv(report_name, sep = ' ')

    # Get protein sequences
    seqs = {}
    for file in X.index.to_list():
        structure = AlphaFoldStructure(f"{config['AF_labels_core']}/{args.problem}/{file}.pdb")
        seqs[file + '.pdb'] = structure.seq


    protein_popularity = {}
    for rank in range(0, 5):

        # Map shapemer to residues
        shapemer = feature_importance.at[rank, 'Feature']
        importance = feature_importance.at[rank, 'Importance']
        residue_indices = shapemer_class.map_shapemer_to_residues(shapemer)

        # Residue indeces per protein
        with open(report_name, 'a') as report:
            report.write(f"\nshapemer: {shapemer}, rank: {rank + 1}, importance: {importance}, Number of proteins with k-mer: {len(residue_indices)}\n\n")

            for pdb_name in residue_indices:
                start = db[db['UniProt'] == pdb_name.replace('.pdb', '')]['Core'].values
                start = start[0][0]
                indeces = [i + start for i in residue_indices[pdb_name]]
                residues = [seqs[pdb_name][i - start] for i in indeces]
                report.write(f"{pdb_name} --> {'+'.join(map(str, indeces))} --> {residues}\n")

            report.write(f'\n{"-"*80}\n')

        # See which protein is gets more repeated
        for pdb_name in residue_indices:
            if pdb_name not in protein_popularity:
                protein_popularity[pdb_name] = 1
            else:
                protein_popularity[pdb_name] += 1
    print(protein_popularity)

    # Feature importance SHAP
    #explainer_lpmo = shap.Explainer(logreg.predict, X_test, max_evals = 1000)
    #shap_values = explainer_lpmo(X_test)
    #shap.plots.bar(shap_values)
    #plt.savefig("shap_bar.png", transparent = False)

    #shap.summary_plot(shap_values, plot_type='violin', show = False)
    #plt.savefig("shap_violin.png", transparent = False)



def main():
    '''Program flow.'''

    # Argument handling
    args = handle_arguments()
    
    # Optimize hyperparameters
    if args.subparser == 'optimization':
        parameter_set = optimize(args)
        print(parameter_set)

    # Feature selection
    elif args.subparser == 'selection':
        feature_importance(args)

if __name__ == '__main__':
    main()
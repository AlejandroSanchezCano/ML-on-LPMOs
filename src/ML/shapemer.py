'''
shape-mer-based supervised machine learning

This file holds the Shapemer class, which allows to perform structure-based 
supervised machine learning with protein sequence shape-mer counts as the
embedding.

Class
-----
Shapemer
'''

import os
import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
from ..config.config_parser import config
from sklearn.model_selection import  KFold
from sklearn.preprocessing import MaxAbsScaler
from sklearn.linear_model import LogisticRegression
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning


warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
import umap as UMAP
import geometricus as gm

class Shapemer():
    
    def __init__(self, problem, input_dir):
        self.problem = problem
        self.db = pd.read_pickle(
            f'{config["supervised"]}/{problem}.pkl'
            ).sort_values('UniProt')
        self.input_dir = input_dir
        self.y = pd.DataFrame(self.db['y']).set_index(self.db['UniProt'])['y']

        self.X = None
        self.proteins = None

###############################################################################
####################          CREATE EMBEDDING         ########################
###############################################################################

    def make_invariants(self, kmer : int, radius : int) -> list:
        '''
        Run first part of a Geometricus experiment: invariants calculation.

        Parameters
        ----------
        kmer : int
            K-mer fragmentation.
        
        radius : int
            Radius fragmentation.

        Returns
        -------
        invariants : list
            Moment invariants.
        '''
        
        # Determine split_info
        if not kmer:
            split_infos = [gm.SplitInfo(gm.SplitType.RADIUS, int(radius))]
        elif not radius:
            split_infos = [gm.SplitInfo(gm.SplitType.KMER, int(kmer))]
        else:
            split_infos = [
                gm.SplitInfo(gm.SplitType.KMER, int(kmer)), 
                gm.SplitInfo(gm.SplitType.RADIUS, int(radius))
                ]
        
        # Invariants
        invariants, _ = gm.get_invariants_for_structures(
            self.input_dir, 
            n_threads = 20,
            split_infos = split_infos,
            moment_types = ["O_3", "O_4", "O_5", "F"]
            )
        
        return invariants
    
    def make_embedding(self, invariants : list, resolution : float):
        '''
        Use Geometricus to calculate the shapemers from the invariants.
        Additionally, normalize the embedding and add the column and index
        names (shapemers and UniProt IDs, respectively).

        Parameters
        ----------
        invariants : list
            Moment invariants.

        resolution : float
            Shapemer resolution.

        Returns
        -------
        X : pd.DataFrame
            Embedding
        '''
    
        # Shapemer class
        self.shapemers = gm.Geometricus.from_invariants(
                invariants, 
                protein_keys = os.listdir(self.input_dir),
                resolution = resolution
                )

        # Shapemer count matrix
        shapemer_count_matrix = self.shapemers.get_count_matrix()

        # Normalize by length
        shapemer_sum = np.sum(shapemer_count_matrix, axis = 1)
        embedding = shapemer_count_matrix/shapemer_sum[:, None]

        # Shapemer count data frame
        proteins = [protein.replace('.pdb', '') for protein in self.shapemers.protein_keys]
        self.proteins = proteins
        shapemers = self.shapemers.shapemer_keys
        self.X = pd.DataFrame(embedding, index = proteins, columns = shapemers)

        # Sort embedding and response alphabetically
        self.X = self.X.sort_index()

        return self.X
    
    def normalize_features(self) -> pd.DataFrame:
        # The columns and indeces are lost during normalization
        columns = self.X.columns
        index = self.X.index

        # Normalize features
        transformer = MaxAbsScaler().fit(self.X)
        scaled_embedding = transformer.transform(self.X)

        # Get columns and indeces back
        self.X = pd.DataFrame(
            data = scaled_embedding, 
            columns = columns, 
            index = index
            )
        
        # Shuffle to avoid homology bias
        #self.X = self.X.sample(frac = 1, axis = 0, random_state = 1999)
        #self.y = self.y.sample(frac = 1, random_state = 1999)

        return self.X
    
    def scalar_embedding(self, scalar : int) -> pd.DataFrame:
        self.X *= scalar
        return self.X

###############################################################################
#################         DIMENSIONALITY REDUCTION         ####################
###############################################################################   

    def umap(self):
        umap = UMAP.UMAP(metric = "cosine").fit_transform(self.X)[:, :2]
        umap = pd.DataFrame(umap, index = self.proteins, columns = [0, 1])
        plt.scatter(umap[0], umap[1])
        plt.show()

        return umap

###############################################################################
####################         KFOLD PERFORMANCE         ########################
###############################################################################   
    
    def __measure_performance(self, y_train_pred, y_test_pred, y_probas):

        metric = {
            'train_accuracy' : metrics.accuracy_score(self.y, y_train_pred),
            'test_accuracy' : metrics.accuracy_score(self.y, y_test_pred),
            'b_accuracy' : metrics.balanced_accuracy_score(self.y, y_test_pred),
            'f1' : metrics.f1_score(self.y, y_test_pred),
            'auc' : metrics.roc_auc_score(self.y, y_probas[:, 1])
        }          

        return metric
    
    def __plot_roc(self, y_probas):
        
        y_pred_proba = y_probas[:, 1]
        fpr, tpr, _ = metrics.roc_curve(self.y,  y_pred_proba)     
        auc = metrics.roc_auc_score(self.y, y_pred_proba)  
        plt.plot(fpr, tpr, label = "auc=" + str(auc))
        plt.legend(loc = 4)
        #plt.savefig( "a.png", transparent = False) 
        plt.show()

    def __plot_confusion_matrix(self, y_pred):
        cm = metrics.confusion_matrix(self.y, y_pred)
        plot = metrics.ConfusionMatrixDisplay(confusion_matrix = cm)
        plot.plot()
        plt.show()
    
    def k_fold(self, X, y, k, plot = True):
        
        # Set CV extremes
        if k == 'loocv':
            k = len(X)
        elif k == 'validation_set':
            k = 2
        
        # K-fold cross-validation
        kf = KFold(n_splits = k, shuffle = False)

        # Initialize values that will be used for calculating metrics
        y_train_pred = [[] for _ in range(len(X))]
        y_test_pred = []
        y_probas = []

        # Iterate over 
        for train_index, test_index in kf.split(X):

            # Test and train sets
            X_train = X.iloc[train_index]
            y_train = y.iloc[train_index]
            X_test = X.iloc[test_index]
            y_test = y.iloc[test_index]

            # Train logistic regression model
            model = LogisticRegression(
                max_iter = 10_000,
                random_state = 1999
                )
            fit = model.fit(X_train, y_train)

            # Calculate output values
            y_test_pred += fit.predict(X_test).tolist()
            y_probas += fit.predict_proba(X_test).tolist() 

            # Organize y_train_pred
            for index, y_pred in zip(train_index, fit.predict(X_train).tolist()):
                y_train_pred[index].append(y_pred)

        # Calculate mode of y_train_pred
        y_train_pred = [max(i, key = i.count) for i in y_train_pred]

        # Create numpy arrays
        y_train_pred = np.array(y_train_pred)
        y_test_pred = np.array(y_test_pred)
        y_probas = np.array(y_probas)
        
        # Calculate performance
        metrics = self.__measure_performance(
            y_train_pred = y_train_pred,
            y_test_pred = y_test_pred,
            y_probas = y_probas
            )
        
        # Manage plotting
        if plot:
            self.__plot_roc(y_probas = y_probas)
            self.__plot_confusion_matrix(y_pred = y_pred)

        return metrics

    def range_kfold(self, X, y):
        
        # List with accumulated metrics as tuple per split
        accumulated_metrics = {
            'train_accuracy' : [],
            'test_accuracy' : [],
            'b_accuracy' : [],
            'f1' : [],
            'auc' : []
        }

        # Iterate over all splits
        split_range = range(2, len(X) + 1)
        for n_splits in tqdm(split_range):
            metrics = self.k_fold(X = X, y = y, k = n_splits, plot = False)
            for metric in metrics:
                accumulated_metrics[metric].append(metrics[metric])
        
        # Plot
        for metric, values in accumulated_metrics.items():
            plt.plot(split_range, values, label = metric)

        plt.legend()
        plt.show()

        return accumulated_metrics

###############################################################################
################          HYPERPARAMETER OPTIMIZATION         #################
###############################################################################

    def optimize_geometricus(self, kmers, radiuses, resolutions):
        
        # Initialize output variable
        accuracies = {
            'train_accuracy' : [],
            'test_accuracy' : [],
            'b_accuracy' : [],
            'f1' : [],
            'auc' : []
        }

        combinations = []

        # Iterate over the hyperparameters
        pbar = tqdm(total = kmers.size * radiuses.size * resolutions.size)
        for kmer in kmers:
            for radius in radiuses:
                
                # Deal with case when both splits are 0
                if kmer == 0 and radius == 0:
                    # Add combinations
                    for resolution in resolutions:
                        combinations.append((kmer, radius, resolution))

                        #Add accuracy = 0
                        for metric in accuracies:
                            accuracies[metric].append(0)
                    
                    pbar.update(resolutions.size)
                    continue

                # Calculate invariants from GeometricusExperiment
                invariants = self.make_invariants(kmer = kmer, radius = radius)

                # The bottleneck in a GeometricusExperiment run is calculating the 
                # invariants, so the resultion gets applied after that
                for resolution in resolutions:
                    print(kmer, radius, resolution)

                    # Add comnination
                    combinations.append((kmer, radius, resolution))
                    
                    # Calculate embedding
                    X = self.make_embedding(
                        invariants = invariants, 
                        resolution = resolution
                        )
                    X = self.normalize_features()
                    X = self.scalar_embedding(100)
                
                    # Accuracy = 0 when only one shapemer is detected
                    if self.X.shape[1] == 1:
                        for metric in accuracies:
                            accuracies[metric].append(0)
                        continue
                    
                    # Logistic regression with k-fold
                    result = self.k_fold(X = self.X, y = self.y, k = 20, plot = False)

                    # Append results
                    for metric in result:
                        accuracies[metric].append(result[metric])

                    # Update progress bar
                    pbar.update()

        # Close progress bar
        pbar.close()

        # Sort accuracies dict
        parameter_options = pd.DataFrame(accuracies, index = combinations)
        
        return parameter_options

# Initialize and embedding
shapemer = Shapemer(problem = 'SS', input_dir = f"{config['AF_labels_core']}/SS")
#shapemer = Shapemer(problem = 'C1C4', input_dir = f"{config['AF_core']}")
#invariants = shapemer.make_invariants(kmer = 6, radius = 0)
#X = shapemer.make_embedding(invariants = invariants, resolution = 2)
#X = shapemer.normalize_features()
#X = shapemer.scalar_embedding(100)


# K-fold CV
# result = shapemer.k_fold(X = X, y = shapemer.y, k = 'loocv', plot = False)
# result2 = shapemer.range_kfold(X = X, y = shapemer.y)


# Hyper parameter optimization
kmers = np.arange(0, 12)
radiuses = np.arange(0, 12)
resolutions = np.array([0.001, 0.01, 0.1, 0.5, 1, 1.5, 2, 2.5, 4, 6, 8, 10])
combinations = shapemer.optimize_geometricus(kmers, radiuses, resolutions)
scores = combinations.sort_values('b_accuracy', ascending = False)
print(scores.iloc[:100])
scores.to_csv('optimization_ss.csv')


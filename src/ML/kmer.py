'''
K-mer-based supervised machine learning

This file holds the Kmer class, which allows to perform sequence-based 
supervised machine learning with protein sequence k-mer counts as the
embedding.

Class
-----
Kmer
'''

import os
import re
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from collections import Counter
import sklearn.metrics as metrics
from collections import defaultdict
from ..config.config_parser import config
from sklearn.model_selection import  KFold
from sklearn.preprocessing import MaxAbsScaler
from sklearn.linear_model import LogisticRegression
from ..filtering.parse_pdb import AlphaFoldStructure

class KMer():

    def __init__(self, k, problem):
        self.k = k
        self.problem = problem
        self.title = 'Substrate specificity' if problem == 'SS' else 'Regioselectivity'
        self.db = pd.read_pickle(f'{config["supervised"]}/{problem}.pkl')
        self.y = pd.DataFrame(self.db['y']).set_index(self.db['UniProt'])['y']

        self.X = None
        self.important_features = None
        self.mapped_features = None
        self.accuracy2n_features = {}
        self.top = {}

###############################################################################
####################          CREATE EMBEDDING         ########################
###############################################################################

    def make_embedding(self):

        # Initialize proto-embedding
        uniprot2kmercounts = {}

        # Get sequences
        for uniprot in self.db['UniProt']:
            structure = AlphaFoldStructure(
                f"{config['AF_labels_core']}/{self.problem}/{uniprot}.pdb"
                )
            seq = structure.seq
            
            # Count kmers
            kmer_counts = Counter(
                seq[i : i + self.k] for i in range(len(seq) + 1 - self.k)
                )
            
            # Add to dict
            uniprot2kmercounts[structure.id] = kmer_counts

        # Make embedding
        self.X = pd.DataFrame(uniprot2kmercounts).T.fillna(0)
        self.X_pure = self.X
        return self.X
    
    def normalize_features(self):
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
        self.X = self.X.sample(frac = 1, axis = 0, random_state = 1999)
        self.y = self.y.sample(frac = 1, random_state = 1999)

        return self.X

    def scalar_embedding(self, scalar):
        self.X *= scalar
        return self.X

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
        plt.ylabel('True positive rate')
        plt.xlabel('False positive rate')
        plt.legend(loc = 4)
        plt.title(self.title)
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

            #print(len(train_index), train_index)
            #print(len(y_train_pred), y_train_pred)
            #print('-'*80)

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
            self.__plot_confusion_matrix(y_pred = y_test_pred)

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

        plt.ylabel('Metric value')
        plt.xlabel('Number of splits')
        plt.title(self.title)
        plt.legend()
        plt.show()

        return accumulated_metrics

###############################################################################
####################         CALCULATE FEATURES        ########################
###############################################################################

    def feature_importance(self):
        # Run ML algorithm with all the data
        X_train = self.X
        y_train = self.y
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
        self.important_features = feature_importance


        # Plot feature importance
        #plot = feature_importance.plot(
        #    x = 'Feature', 
        #     y = 'Importance', 
        #    kind = 'barh',
        #    )
        #plt.show()
        #fig = plot.get_figure()
        #fig.savefig(
        #    f"{config['plots']}/{self.k}_{self.problem}.png", 
        #    transparent = True
        #    )

        return self.important_features

###############################################################################
#########################        MAP FEATURES        ##########################
###############################################################################

    def map_features(self, n_features):

        # Initialize return variable
        mapped_features = {
            'k-mer' : [],
            'Importance' : [],
            'Number of proteins' : [],
            'Proteins' : [],
            'Indeces' : []
        }
        
        # Iterate over the top features
        for rank in range(n_features):
            # Get k-mer, coefficient and proteins that have said k-mer
            kmer = self.important_features.at[rank, 'Feature']
            importance = self.important_features.at[rank, 'Importance']
            proteins = self.X[self.X[kmer] != 0].index.values

            # Iterate over the proteins that have the k-mer to get the index
            indeces = []
            for uniprot in proteins:

                # Find protein start
                uniprot2core = dict(zip(self.db['UniProt'], self.db['Core']))
                protein_start = uniprot2core[uniprot][0]

                # Get proteinsequence
                seq = AlphaFoldStructure(
                    f"{config['AF_labels_core']}/{self.problem}/{uniprot}.pdb"
                    ).seq
                
                # Get indeces
                indeces += [[re_match.start() + protein_start \
                           for re_match in re.finditer(kmer, seq)]]
                
            # Add to proto data frame
            mapped_features['k-mer'].append(kmer)
            mapped_features['Importance'].append(importance)
            mapped_features['Number of proteins'].append(len(proteins))
            mapped_features['Proteins'].append(proteins)
            mapped_features['Indeces'].append(indeces)

        # Build data frame
        self.mapped_features = pd.DataFrame(mapped_features)

        return self.mapped_features
    
    def map_features_report(self, n_features):
        

        # Make template of dataframe with kmers as rows and proteins as columns
        important_kmers = self.important_features.head(n_features)['Feature']
        template = self.X_pure[important_kmers].T
        template['Importance'] = self.important_features['Importance'].iloc[:n_features].to_list()
        
        # Fill dataframe with index of features
        for uniprot in tqdm(template.drop('Importance', axis = 1)):
            
            # Find protein start
            uniprot2core = dict(zip(self.db['UniProt'], self.db['Core']))
            protein_start = uniprot2core[uniprot][0]

            # Get proteinsequence
            seq = AlphaFoldStructure(
                f"{config['AF_labels_core']}/{self.problem}/{uniprot}.pdb"
                ).seq
                
            # Iterate over the kmers
            for kmer in important_kmers:

                # Get indeces
                indeces = [re_match.start() + protein_start + 1\
                           for re_match in re.finditer(kmer, seq)]
                indeces = ', '.join(map(str, indeces))
                
                # Replace by the index
                template.at[kmer, uniprot] = indeces
        
        return template.T

    def features_report(self, n_features):
        
        path = f"{config['ML']}/kmer/{self.k}_{self.problem}.txt"
        with open(path, 'w') as report:
            # Write important features
            self.important_features.iloc[0:10].to_csv(report, sep = ' ')

            # Write scores
            scores = self.k_fold(X = self.X, y = self.y, k = 'loocv')
            for name, value in scores.items():
               report.write(f'\n{name} = {value}\n')

            # Write report per k-mer
            for rank in range(n_features):
                row = self.mapped_features.iloc[rank]
                kmer, importance, n_proteins, proteins, indeces = row
                report.write(f"\nkmer: {kmer}, rank: {rank + 1}, importance: {importance}, Number of proteins with k-mer: {n_proteins}\n\n")

                # Proteins that have a k-mer
                for protein, index in zip(proteins, indeces):
                    consecutive = '+'.join(map(str, [index[0] + i for i in range(self.k)]))
                    report.write(f'{protein} --> {consecutive}\n')
                
                # Write separation
                report.write(f'\n{"-"*80}\n')

###############################################################################
######################        PYMOL VISUALIZATION        ######################
############################################################################### 

    def _hex_to_RGB(self, hex_str):
        """ #FFFFFF -> [255,255,255]"""
        #Pass 16 to the integer function for change of base
        return [int(hex_str[i:i+2], 16) for i in range(1,6,2)]

    def _get_color_gradient(self, c1 = "#2BBFA7", c2 = "#FFC600", n = 33):
        """
        Given two hex colors, returns a color gradient
        with n colors.
        """
        assert n > 1
        c1_rgb = np.array(self._hex_to_RGB(c1))/255
        c2_rgb = np.array(self._hex_to_RGB(c2))/255
        mix_pcts = [x/(n-1) for x in range(n)]
        rgb_colors = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
        return ["0x" + "".join([format(int(round(val*255)), "02x") for val in item]) for item in rgb_colors]


    def report2pymol(self, n_features):

        # Initialize export variable
        # Initiate at 0 cause that does nto affect Pymol
        uniprot2indeces = {uniprot:['0' for _ in range(n_features)] for uniprot in self.db['UniProt']}

        # Iterate over the mapped features rows
        for rank in range(len(self.mapped_features)):
            row = self.mapped_features.iloc[rank]
            *_, proteins, indeces = row

            # Make protein to indeces dictionary (e.g. UniProt ID : [10])
            for protein, index in zip(proteins, indeces):

                # Add k consecutive residues
                new_index = []
                for i in index:
                    new_index += [i + k for k in range(1, self.k + 1)]
                
                # Add notation pymol (34+35)
                new_index = '+'.join(map(str, new_index))

                # Add to dict
                uniprot2indeces[protein][rank] = new_index
            
        # Generate color gradient
        gradient = self._get_color_gradient(
            c1 = "#2BBFA7",
            c2 = "#FFC600",
            n = len(self.mapped_features)
            )
        
        # Make pml file
        code = f'cd {config["AF_all"]}\n\n'
        for uniprot in uniprot2indeces:
            code += f'# {uniprot}\n'
            code += f'load {uniprot}.pdb\n'
            code += 'bg white\n'
            code += 'color white\n'
            for n, rank_indeces in enumerate(uniprot2indeces[uniprot]):
                kmer = self.mapped_features["k-mer"].iloc[n]
                sign = '+' if self.mapped_features["Importance"].iloc[n] > 0 else '-'
                code += f'select selection, resi {rank_indeces}\n'
                code += f'color {gradient[n]}, selection\n'
                code += 'show sticks, selection\n'
                code += f'label first selection, "{sign}{kmer}_{n + 1}"\n'
                code += 'set label_color, black\n'
                code += 'set label_size, -2\n'
            code += f'save {config["pymol"]}/{self.k}_{self.problem}/{uniprot}.pse\n'
            code += 'delete all\n'
            code += '\n\n'

        # Save pml file
        file_name = 'run_pymol.pml'
        with open(file_name, 'w') as run_file:
            run_file.write(code)

###############################################################################
#####################        PORTION OF FEATURES        #######################
###############################################################################    

    def kfold_per_feature(self, k):
        # Initialize variables
        top = {
            'train_accuracy' : [],
            'test_accuracy' : [],
            'b_accuracy' : [],
            'f1' : [],
            'auc' : []
        }
        bad = {
            'train_accuracy' : [],
            'test_accuracy' : [],
            'b_accuracy' : [],
            'f1' : [],
            'auc' : []
        }

        # Iterate over the feature index (of importance)
        for n_features in tqdm(range(1, len(self.important_features))):

            # Get top and non_top features
            top_features = self.important_features['Feature'][0:n_features]
            
            # Make new embeddings
            embedding_top_features = self.X[top_features]
            embedding_bad_features = self.X.drop(top_features, axis = 1)

            # Prediction only with top features
            performance_top = self.k_fold(
                X = embedding_top_features, 
                y = self.y,
                k = k,
                plot = False
                )
            performance_bad= self.k_fold(
                X = embedding_bad_features, 
                y = self.y,
                k = k,
                plot = False
                )

            # Append results
            for metric in top:
                top[metric].append(performance_top[metric])
                bad[metric].append(performance_bad[metric])
        
        # Make attribute (of numpy arrays)
        self.top = {key : np.array(value) for key, value in top.items()}

        # Plot
        x = list(range(1, len(self.important_features)))
        plt.plot(x, top['b_accuracy'])
        plt.plot(x, bad['b_accuracy'][::-1])
        plt.xlabel('Feature importance index')
        plt.ylabel('Balanced accuracy')
        plt.title(self.title)
        plt.legend(['Top features', 'Bad features'])
        plt.show()

        # Save picture
        #plt.savefig(
        #    f'{config["plots"]}/portion_{self.k}_{self.problem}.png', 
        #    transparent = False
        #    )

        return top, bad
        
    def n_features_for_metric(self, metric, value):

        # Index of the element
        index = np.argmax(self.top[metric] >= value)
        
        if index == 0 and value > np.max(self.top[metric]):
            return -1
        else:
            return index + 1
    
###############################################################################
#########        EXCHANGE RESPONSE WITH MOST SIMILAR PROTEIN        ###########
############################################################################### 

    def most_similar_cluster(self):

        # Import agglomeration
        agglomeration = pd.read_pickle(f'{config["Dendrogram"]}/labels/agglomeration.pkl')
 
        # Initialize UniProt ID : most similar UniProt IDs
        uniprot2mostsimilars = {}

        # UniProt ID : binary response
        uniprot2response = dict(zip(self.db['UniProt'], self.db['y']))

        # Import sequence-based similarity matrix
        matrix = pd.read_pickle(f'{config["Dendrogram"]}/labels/similarity_matrix.pkl')

        # Remove entries that don't apply to the specific problem
        to_drop = set(matrix.columns) - set(self.db['UniProt'])
        matrix = matrix.drop(to_drop, axis = 0)
        matrix = matrix.drop(to_drop, axis = 1)

        # Identify most similar neighbour
        for uniprot in matrix.columns:
            similarities = matrix.loc[uniprot].drop(index = uniprot)
            ordered_similarities = similarities.sort_values(ascending = False)
            uniprot2mostsimilars[uniprot] = ordered_similarities.index.to_list()

        # Initiate plotting variable
        accumulated_metrics = {
            'train_accuracy' : [],
            'test_accuracy' : [],
            'b_accuracy' : [],
            'f1' : [],
            'auc' : []
        }

        # Iterate over clustering steps
        for step in tqdm(range(len(agglomeration) - 1)):
            # Get clusters per set
            clusters = agglomeration[step]
            # Initialize UniProt ID : response of the most similar UniProt ID
            uniprot2similarresponse = {}
            # Iterate over the UniProt IDs
            for uniprot in matrix.columns:
                # Calculate most similar UniProt outside of cluster
                for cluster in clusters:
                    if uniprot in cluster:

                        other_cluster = [uniprot \
                                for uniprot in uniprot2mostsimilars[uniprot] \
                                if uniprot not in cluster]
                        
                        uniprot2similarresponse[uniprot] = uniprot2response[other_cluster[0]]
                        
            # Create similar response
            y_similar = self.db['UniProt'].apply(
                lambda uniprot : uniprot2similarresponse[uniprot]
            )
            sampled_y_similar = y_similar.sample(frac = 1, random_state = 1999)

            # Count mismatches
            mismatch = np.sum(self.db['y'] != y_similar)

            # Perform k-fold 
            # Error may arise when all of the one class (e.g. C4) are converted
            # to the other class (e.g. C1) so the class imbalance is that strong
            # that you cannot train the model.
            try:
                metrics = self.k_fold(
                    X = self.X, 
                    y = sampled_y_similar,
                    k = 20, 
                    plot = False
                )
            except ValueError:
                metrics = {metric : None for metric in metrics}

            # Append to plotting variable
            for metric in metrics:
                accumulated_metrics[metric].append(metrics[metric])

        # Plot
        plt.plot(
            range(len(agglomeration) - 1), 
            accumulated_metrics['b_accuracy'],
            'g'
            )
        plt.xlabel('Clustering steps')
        plt.ylabel('Balanced accuracy')
        plt.title(self.title)
        plt.show()

        
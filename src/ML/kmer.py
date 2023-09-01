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
import math
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from collections import Counter
from collections import defaultdict
from ..config.config_parser import config
from sklearn.preprocessing import MaxAbsScaler
from sklearn.linear_model import LogisticRegression
from ..filtering.parse_pdb import AlphaFoldStructure
from sklearn.model_selection import LeaveOneOut, cross_validate

class KMer():

    def __init__(self, k, problem):
        self.k = k
        self.problem = problem
        self.db = pd.read_pickle(f'{config["supervised"]}/{problem}.pkl')


        self.embedding = None
        self.scores = None
        self.important_features = None
        self.mapped_features = None
        self.accuracy2n_features = {}

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
        self.embedding = pd.DataFrame.from_dict(
            data = uniprot2kmercounts, 
            orient = 'index'
            ).fillna(0) 
            
        # Use same order in database and in emebdding data frames
        sorter = self.embedding.index.to_list()
        self.db['UniProt'] = self.db['UniProt'].astype('category')
        self.db['UniProt'] = self.db['UniProt'].cat.set_categories(sorter)
        self.db = self.db.sort_values(by = 'UniProt')

        return self.embedding
    
    def normalize_features(self):
        # The columns and indeces are lost during normalization
        columns = self.embedding.columns
        index = self.embedding.index

        # Normalize features
        transformer = MaxAbsScaler().fit(self.embedding)
        scaled_embedding = transformer.transform(self.embedding)

        # Get columns and indeces back
        self.embedding = pd.DataFrame(
            data = scaled_embedding, 
            columns = columns, 
            index = index
            )

        return self.embedding

    def scalar_embedding(self, scalar):
        self.embedding *= scalar
        return self.embedding
    
    def loocv(self, X, y):
        cv = LeaveOneOut()
        logreg = LogisticRegression(
            max_iter = 10000,
            random_state = 1999
            )
        scores = cross_validate(
            estimator = logreg, 
            X = X, 
            y = y, 
            scoring = ['accuracy'], #'balanced_accuracy'],# 'f1'],
            cv = cv, 
            n_jobs = -1
            )
        self.scores = {key:np.mean(value) for key, value in scores.items()\
                    if key.startswith('test')}
        
        return self.scores
    
    def feature_importance(self):
        # Run ML algorithm with all the data
        X_train = self.embedding
        y_train = self.db['y']
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

        #report_name = f'{config["ML"]}/kmer/{self.k}_{self.problem}.txt'
        #feature_importance.iloc[0:10].to_csv(report_name, sep = ' ')
        '''
        # Plot feature importance
        plot = feature_importance.plot(
            x = 'Feature', 
            y = 'Importance', 
            kind = 'barh',
            )
        fig = plot.get_figure()
        fig.savefig(
            f"{config['plots']}/{self.k}_{self.problem}.png", 
            transparent = True
            )
        '''
        
        return self.important_features
        
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
            proteins = self.embedding[self.embedding[kmer] != 0].index.values

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
    
    def features_report(self, n_features):
        
        path = f"{config['ML']}/kmer/{self.k}_{self.problem}.txt"
        with open(path, 'w') as report:
            # Write important features
            self.important_features.iloc[0:10].to_csv(report, sep = ' ')

            # Write scores
            for name, value in self.scores.items():
               report.write(f'\n{name} = {value}\n')

            # Write report per k-mer
            for rank in range(n_features):
                row = self.mapped_features.iloc[rank]
                kmer, importance, n_proteins, proteins, indeces = row
                report.write(f"\nkmer: {kmer}, rank: {rank + 1}, importance: {importance}, Number of proteins with k-mer: {n_proteins}\n\n")

                # Proteins that have a k-mer
                for protein, index in zip(proteins, indeces):
                    report.write(f'{protein} --> {index}\n')
                
                # Write separation
                report.write(f'\n{"-"*80}\n')

    def report2pymol(self):
        # Initialize export variable
        uniprot2indeces = defaultdict(list)
        
        # Iterate over the mapped features rows
        for rank in range(len(self.mapped_features)):
            row = self.mapped_features.iloc[rank]
            *_, proteins, indeces = row
            # Make protein to indeces dictionary (e.g. UniProt ID : [10])
            for protein, index in zip(proteins, indeces):
                uniprot2indeces[protein] += index
        
        # Make directory with AF files and code
        dir_path = f'{config["pymol"]}/{self.k}_{self.problem}/UniProt'
        os.system(f'mkdir -p {dir_path}')

        # Export .pdb files to local
        for uniprot in uniprot2indeces:
            uniprot_path = f'{config["AF_all"]}/{uniprot}.pdb'
            os.system(f'cp {uniprot_path} {dir_path}')
        
        # Code for .pml file
        code = 'cd C:/Users/alexs/OneDrive/Escritorio\n'
        for uniprot, indeces in uniprot2indeces.items():
            # Load file
            local_path = f'M2/MSc_thesis/pymol/{self.k}_{self.problem}'
            code += f'load {local_path}/UniProt/{uniprot}.pdb\n'

            # Color protein
            code += 'color white\n'

            # Color selection
            for index in indeces:
                code += f'select selection, resi {index + 1}-{index + self.k}\n'
                code += 'color red, selection\n' 
                code += 'show sticks, selection\n'
    
            # Save session
            code += f'save {local_path}/Sessions/{uniprot}.pse\n\n'

            break

        # Make .pml file
        file_name = f'{config["pymol"]}/{self.k}_{self.problem}/run_pymol.pml'
        with open(file_name, 'w') as run_file:
            run_file.write(code)

    def loocv_portion(self):
        # Initialize return variables
        top_features_accuracies, bad_features_accuracies = [], []

        # Iterate over the feature index (of importance)
        for n_features in tqdm(range(1, len(self.important_features))):

            # Get top and non_top features
            features = self.important_features['Feature']
            top_features = self.important_features['Feature'][0:n_features]
            
            # Make new embeddings
            embedding_top_features = self.embedding[top_features]
            embedding_bad_features = self.embedding.drop(
                top_features, 
                axis = 1
                )

            # Prediction only with top features
            acc_top = self.loocv(X = embedding_top_features, y = self.db['y'])
            acc_bad = self.loocv(X = embedding_bad_features, y = self.db['y'])

            # Append results
            top_features_accuracies.append(acc_top['test_accuracy'])
            bad_features_accuracies.append(acc_bad['test_accuracy'])

        # Plot
        x = list(range(1, len(self.important_features)))
        plt.plot(x, top_features_accuracies)
        plt.plot(x, bad_features_accuracies[::-1])
        plt.xlabel('Feature importance index')
        plt.ylabel('Accuracy')
        plt.legend(['Top features', 'Bad features'])

        # Save picture
        plt.savefig(
            f'{config["plots"]}/portion_{self.k}_{self.problem}.png', 
            transparent = False
            )
        
        # Calculate how many features are needed for a given accuracy
        for n, accuracy in enumerate(top_features_accuracies):
            # Start with the accuracy given by the best feature
            if n == 0:
                self.accuracy2n_features[accuracy] = n + 1
                milestone = 0.1 + math.floor(accuracy * 100)/100

            # Stop if we reach 100% accuracy
            elif accuracy == 1:
                self.accuracy2n_features[accuracy] =  n + 1
                break

            # Add number of features every time we find a milestone
            elif accuracy > milestone:
                self.accuracy2n_features[accuracy] = n + 1
                milestone = 0.1 + math.floor(accuracy * 100)/100
        print(self.accuracy2n_features)
            

kmer = KMer(k = 2, problem = 'C1C4')
embedding = kmer.make_embedding()
embedding = kmer.normalize_features()
embedding = kmer.scalar_embedding(100)
print('Calculate accuracy')
accuracy_loocv = kmer.loocv(X = kmer.embedding, y = kmer.db['y'])
print('Important features')
important_features = kmer.feature_importance()
print('Map features')
mapped_features = kmer.map_features(3)
print('Generate report')
kmer.features_report(3)
#print('Portion')
#kmer.loocv_portion()
print('To PyMol')
kmer.report2pymol()


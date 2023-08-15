'''
Geometricus embedder

This file holds the GeometricusExperiment class, which serves as an interface 
with Geometricus and the following experiments one could do with a Geometricus
embedding such as dimensionality reduction. 

Class
-----
GeometricusExperiment
'''

import os
import warnings
import numpy as np
import pandas as pd
import geometricus as gm
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from typing import Any, Union, Callable

# Import umap without warnings 
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
import umap

class GeometricusExperiment():
    '''
    Class that represents a Geometricus experiment.
    This class can is able to take a directory with AlphaFold structure .pdb
    files and apply Geometricus. The output will be an embedding where one
    can run machine learning algorithms. This class also provides methods to 
    run dimensionality reduction techniques and plot the results.

    Attributes
    ----------
    input_dir : str
        Directory containing AlphaFold structure .pdb files to run a
        Geometricus experiment on.

    shapemers : list
        List of shapemers

    _embedding : pd.DataFrame
        Untransformed embedding
    '''

    def __init__(self, input_dir : str):
        '''
        GeometricusExperiment constructor

        Parameters
        ----------
        input_dir : str
            Directory containing AlphaFold structure .pdb files to run a
            Geometricus experiment on.
        '''

        # Input
        self.input_dir = input_dir 

        # Attributes
        self.shapemers = None
        self._embedding = None     

    def __repr__(self):
        return 'GeometricusExperiment' + str(self.__dict__)
    
    def run_invariants(self, kmer : int, radius : int) -> list:
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
            n_threads = 4,
            split_infos = split_infos,
            moment_types = ["O_3", "O_4", "O_5", "F"]
            )
        
        return invariants

    def embedding(self, multiplier : int = 100) -> pd.DataFrame:
        '''
        The normalization of the embedding matrix leaves a sparse matrix with 
        non-zero values quite close to 0. So much that the machine learning and
        dimensionality reduction algorithms have difficulties differnciating
        them. Therefore, it's necessary to modify the embedding data. This is 
        done by simply multypling the embedding matrix by a factor.

        Parameters
        ----------
        multiplier : int
            Correction factor of the embedding to make the non-zero values
            outstand from the widespread zero values. Teh default is 100.

        Returns
        -------
        corrected_embedding : pd.DataFrame
            Corrected embedding.

        '''

        corrected_embedding = self._embedding * multiplier
        return corrected_embedding
    
    def run_embedding(self, invariants : list, resolution : float):
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
        None.
        '''

        # Shapemer class
        self.shapemers = gm.Geometricus.from_invariants(
                invariants, 
                protein_keys = os.listdir(self.input_dir),
                resolution = resolution
                )

        # Shapemer count matrix
        shapemer_count_matrix = self.shapemers.get_count_matrix()
        shapemer_sum = np.sum(shapemer_count_matrix, axis = 1)
        embedding = shapemer_count_matrix/shapemer_sum[:, None]

        # Shapemer count data frame
        proteins = [protein.replace('.pdb', '') for protein in self.shapemers.protein_keys]
        shapemers = self.shapemers.shapemer_keys
        embedding = pd.DataFrame(embedding, index = proteins, columns = shapemers)
        self._embedding = embedding

    def pca(self, *args : Any, **kwargs : Any) -> pd.DataFrame:
        '''
        Calculate PCA scores from embedding.

        Parameters
        ----------
        *args : Any
            PCA arguments.
        **kwargs : Any
            PCA key-word arguments.

        Returns
        -------
        scores : pd.DataFrame
            PCA score data frame.
        '''
        space = PCA(*args, **kwargs).fit_transform(self.embedding())
        scores = pd.DataFrame(
            data = space, 
            columns = [0, 1], 
            index = self._embedding.index
            )
        return scores
    
    def umap(self, *args : Any, **kwargs : Any) -> pd.DataFrame:
        '''
        Calculate UMAP scores from embedding.

        Parameters
        ----------
        *args : Any
            UMAP arguments.
        **kwargs : Any
            UMAP key-word arguments.

        Returns
        -------
        scores : pd.DataFrame
            UMAP score data frame.
        '''

        space = umap.UMAP(*args, **kwargs).fit_transform(self.embedding())
        scores = pd.DataFrame(
            data = space, 
            columns = [0, 1], 
            index = self._embedding.index
            )
        return scores

    def plot_per_category(
            self, 
            uniprot2category : dict[str, str], 
            reduction : str = 'umap',
            sort_key : Callable = None,
            specific_categories : Union[bool, list] = None
         ): 
        '''
        Plot the dimensionality reduction scores per category such as family,
        regioselectivity or substrate specificity.

        Parameters
        ----------
        uniprot2category : dict[str, str]
            Dictionary with UniProt IDs as keys and the category the belong to 
            as values.

        reduction : str
            Dimensionality reduction technique. Choices are 'PCA' for Principal
            Component Analysis and 'UMAP' for Uniform Manifold Approximation 
            and Projection. The default option is UMAP.

        sort_key : function
            Function used to sort the entries in a specific way. This is
            specifically useful for sorting families, where the string 'AA9' 
            should go before 'AA10'.

        specific_categories : Union[bool, list]
            By default, the method generates a plot of all the categories. But
            if one wishes to focus one or a set of specific category or 
            highlight them, they should be specified here.
        '''

        # Handle dimensionality reduction technique
        if reduction == 'umap':
            scores = self.umap(metric = "cosine", random_state = 1999)
        elif reduction == 'pca':
            scores = self.pca(n_components = 2, random_state = 1999)
        
        # Sort keys if necessary (e.g. families) 
        categories = np.unique(list(uniprot2category.values()))
        if sort_key:
            categories = sorted(categories, key = sort_key)
        
        # Plot
        fig, ax = plt.subplots()
        for category in categories:

            # Get UniProt IDs per family
            per_category = [uniprot for uniprot in uniprot2category\
                if uniprot2category[uniprot] == category]

            filter_per_category = scores.loc[per_category]

            # Plot normal
            if not specific_categories or category not in specific_categories:
                ax.scatter(
                    filter_per_category[0], 
                    filter_per_category[1], 
                    label = category
                    )
            
            # Plot when selecting specific classes
            else:
                ax.scatter(
                    filter_per_category[0], 
                    filter_per_category[1],  
                    label = category,
                    alpha = 0.2
                    )
                
        ax.legend()
        plt.savefig('foo.png')
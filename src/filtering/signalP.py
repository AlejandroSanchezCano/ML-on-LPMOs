'''
SignalP results parser

This file holds the SignalP class, which parses the results file from SignalP 
experiments and retrieves relevant information

Class
-----
SignalP
'''

import os
import re
import pandas as pd
from typing import Generator
from .parse_pdb import AlphaFoldStructure
from ..config.config_parser import config

class SignalP():
    '''
    Class that represents a SignalP result.
    This class is able to parse results in both cases when --organism is set 
    to 'other' or 'eukarya'.

    Attributes
    ----------
    df : pd.DataFrame
        Pandas data frame containing the UniProt ID, the SignalP prediction(s),
        where the cut happens if a signal peptide is predicted and different
        probabilities associated with the predictions

    Methods
    -------
    __init__
    _open
    _parse
    '''

    def __init__(self, path : str):
        '''
        SignalP constructor

        Parameters
        ----------
        path : str
            Path to the result file of SignalP
        '''
        self.df = None
        self.other = None
        self.signalp = []
        self._parse(SignalP._open(path))   
        #self._his_around_cut()
    
    def _open(path):
        '''
        Opens txt file and returns a per-line generator

        Parameters
        ---------
        path : str
            Path to the result file of SignalP

        Yields
        ------
        line : str
            Lines of results file in generator format
        '''
        with open(path, 'r') as file:
            for line in file.readlines():
                yield line

    def _parse(self, lines : Generator[str, None, None]):
        '''
        Retrieves relevant information from each line and puts it in a data
        frame.

        Parameters
        ----------
        lines : Generator[str, None, None]
            File lines
        '''

        # Parse lines
        matrix = []
        for line in lines:

            # Get headers
            if line.startswith('#'):
                header = line.strip('# ').rstrip('\n').split('\t')
            
            # Split line
            else:
                row = line.rstrip('\n').split('\t')
                matrix.append(row)
        
        # Make data frame
        self.df = pd.DataFrame(matrix, columns = header)
        self.df = self.df.set_index('ID')

        # Add cut column
        regex = lambda x : int(re.search(r'pos: ([0-9]+)-', x).group(1))\
              if x else 0
        self.df['CUT'] = self.df['CS Position'].apply(regex)
        
    def _his_around_cut(self):
        '''
        Identifies His1 residue in a window of +3 and -3 residues of the 
        cutting point.
        '''

        # Iterate over non-other proteins
        detected = self.df[self.df['Prediction'] != 'OTHER'][['ID', 'CS Position']]
        for uniprot, cut in detected.itertuples(index = False):
            # Get sequence
            structure = AlphaFoldStructure(f"{config['AF_all']}/{uniprot}.pdb")
            seq = structure.seq

            # Identify His
            window_length = 3
            cut = int(cut.split('-')[0])
            window = seq[cut - window_length : cut + window_length]

            if seq[cut] == 'H' or window.find('H') != -1:
                self.signalp.append(uniprot)
            else:
                self.other.append(uniprot)

def main():
    '''Program flow.'''
    
    # Run SignalP
    os.system(f'signalp6 --fastafile {config["FASTA"]}/all.fasta --output_dir {config["SignalP"]} --format none')

if __name__ == '__main__':
    main()
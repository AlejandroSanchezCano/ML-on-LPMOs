'''
Filter and disregard AlphaFold .pdb structures that have a general bad
structure. For that,  the protein will be diveded into peptides of length 
peptide_length. Then, the average pLDDT per peptide is calculated. The protein 
that do not have a peptide with an average pLDDT above a cutoff will not 
subject of future analysis.
It takes 6 min for 4000 proteins.

Functions
---------
handle_arguments
bad_structures
main
'''

import os
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import Union
from .parse_pdb import AlphaFoldStructure
from ..config.config_parser import config

def handle_arguments() -> argparse.Namespace:
    '''
    Handles the arguments passed via command line

    Return
    ------
    args : argparse.Namespace
        Arguments
    '''
    
    parser = argparse.ArgumentParser(
        prog = 'general_filtering',
        description = 'Filters out structure of an overall bad quality',
        epilog = 'See ya'
        )
    
    # Add arguments
    parser.add_argument(
        '-l', '--length', 
        type = int, 
        default = 20,
        help = 'Peptide length resulting from in-silico fragmentation'
        )
    
    parser.add_argument(
        '-c', '--cutoff',
        type = int,
        default = 90,
        help = 'Average pLDDT cutoff of the fragmented peptides' 
    )

    # Arguments from Namespace object
    return parser.parse_args()

def bad_quality(
        peptide_length: int,
        cutoff: Union[int, float] 
        ) -> dict[str, list[AlphaFoldStructure]]:    
    '''
    Impose a filtering based on overall poor quality structure. To do so, the
    protein will be diveded into peptides of length peptide_length. Then,
    the average pLDDT per peptide is calculated. The protein that do not have
    a peptide with an average pLDDT above a cutoff will not subject of future
    analysis.
    Takes 6 min for 3984 proteins.

    Parameters
    ----------
    peptide_length : int, optional
        Length of the peptides that the protein is divided into. The default is
        20.
        
    cutoff : Union[int, float], optional
        pLDDT cutoff used to assess peptide quality. The default is 90.
    
    Returns
    -------
    result : dict[str, list[AlphaFoldStructure]]
        Dictionary containing the AlphaFoldStructure objects that pass the
        filter (values of the key 'pass') and don't (values of the key 
        'nopass').
    '''

    # Instanciate return variable
    result = {
        'pass' : [],
        'nopass' : []
    }

    for structure in tqdm(os.listdir(config['AF_all'])):
        # Get per-residue pLDDT
        af_structure = AlphaFoldStructure(f"{config['AF_all']}/{structure}")
        quality = af_structure.pLDDT
        # Divide protein in peptides -> matrix of pLDDT values
        peptides = [
            quality[position : position + peptide_length] \
            for position in range(0, len(quality), peptide_length)
            ]
        # Calculate mean pLDDT per peptide
        mean_pLDDT_per_peptide = np.array([peptide.mean()\
                                           for peptide in peptides])
        
        # Assess quality > cutoff
        if any(mean_pLDDT_per_peptide >= cutoff):
            result['pass'].append(af_structure.id)
        else:
            result['nopass'].append(af_structure.id)
        
    return result

def update_db(result : dict) -> None:
    '''
    Add a boolean column to the data frame database reflecting whether the 
    UniProt entries have passed the general filtering.

    Parameters
    ----------
    result : dict[str, list[AlphaFoldStructure]]
        Dictionary containing the AlphaFoldStructure objects that pass the
        filter (values of the key 'pass') and don't (values of the key 
        'nopass').

    Returns
    -------
    None
    '''
    
    # Import database
    df = pd.read_pickle(f'{config["CAZy_expanded"]}/all.pkl')

    # Lambda function -> True (pass), False (nopass), None (not applicable)
    general_filter = lambda  alphafold : None if not alphafold else\
          False if alphafold in set(result['nopass']) else True
    
    # Create new column
    df['General filter'] = df['AlphaFold'].apply(general_filter)

    # Store database
    df.to_pickle(f'{config["CAZy_expanded"]}/all.pkl')

def main():
    '''Program flow.'''

    # Arguments
    args = handle_arguments()

    # Filter
    result = bad_quality(args.length, args.cutoff)

    # Update database
    update_db(result)

if __name__ == '__main__':
    main()
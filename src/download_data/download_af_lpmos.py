'''
Download AlphaFold files

This file uses the curl command to download AlphaFold structures from the
UniProt ID. No a priori information is used to know whether each UniProt ID
corrspond to a deleted entry or simply does not have a structure. Therefore, 
running curl will sometimes result in a file with an error message. The low
size of those files (127 bytes) to remove them.

Functions
---------
download_af_files
main
'''

import os
import pandas as pd
from tqdm import tqdm
from ..config.config_parser import config

def download_af_files(df: pd.DataFrame) -> None:
    '''
    Download .pdb structure from the web using curl command calls.
    UniProt IDs are the key that identifies each structure. The logic is:
        1) Get whatever response is received from the curl command.
        2) If there is no AlphaFold structure associated with the unique
           UniProt ID, the response received results in a text file much 
           shorter than a regular .pdb file, so all files < 300 bytes are 
           removed.
           
    Parameters
    ----------
    df : pd.DataFrame
        Data frame containing a the UniProt IDs of the protein whose
        .pdb structure will be obtained

    Returns
    -------
    None.

    '''

    # UniProt IDs
    uniprot_ids = df['UniProt'].dropna().drop_duplicates()

    # Download AF files
    for uniprot_id in tqdm(uniprot_ids):
        version = 'v4'
        model_url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_{version}.pdb'
        os.system(f"curl {model_url} -o {config['AF_all']}/{uniprot_id}.pdb")

    # Remove failed AF files (< 200 bytes)
    for file in tqdm(os.listdir(config['AF_all'])):
        size = os.path.getsize(f"{config['AF_all']}/{file}")
        if size < 300:
            os.remove(f"{config['AF_all']}/{file}")

def main():
    '''Program flow.'''

    # Input supreme frame
    all = pd.read_pickle(f"{config['CAZy_expanded']}/all.pkl")
    
    # Download AF files
    download_af_files(all)

if __name__ == '__main__':
    main()
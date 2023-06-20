import os
import pandas as pd
from tqdm import tqdm
from variables import CAZY_EXPANDED, AF_FILES

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
        os.system(f'curl {model_url} -o {AF_FILES}/{uniprot_id}.pdb')

    # Remove failed AF files (< 200 bytes)
    for file in tqdm(os.listdir(AF_FILES)):
        size = os.path.getsize(f'{AF_FILES}/{file}')
        if size < 300:
            os.remove(f'{AF_FILES}/{file}')

def main():
    '''File flow.'''

        # Input supreme frame
    df = pd.read_pickle(f'{CAZY_EXPANDED}/AA_supreme')
    
    # Download AF files
    download_af_files(df)

if __name__ == '__main__':
    main()
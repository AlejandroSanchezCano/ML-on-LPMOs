'''
Domain Parser for AlphaFold Models.

DPAM automatically recognizes globular domains from AlphaFold models based on 
predicted aligned errors (PAE), inter-residue distances in 3D structures, and
ECOD domains found by sequence (HHsuite) and structural (DALI) similarity 
searches. Parsing LPMOs domains is extremely useful to excise the enzymatic
core from the rest of the protein. The result of the domain parsing is stored 
in a file with extension .finalDPAM.domains. This information is extracted and 
stored in a dictionary with UniProt IDs as keys and domain starting and ending
points as values as a list of 2-length tuples.
'''

import os
import re
import pickle
from tqdm import tqdm
from ..config.config_parser import config

def main():
    '''Program flow.'''

    # Initialize output variable
    uniprot2domains = {}

    # Iterate over proteins to download
    input_proteins = os.listdir(config['AF_all'])
    for pdb_file in tqdm(input_proteins):
        
        uniprot = pdb_file.replace('.pdb', '')
        
        # Important location files
        path_pdb = f"{config['AF_all']}/{uniprot}.pdb"
        path_cif = f"{config['AF_cif']}/{uniprot}.cif"
        path_pae = f"{config['AF_pae']}/{uniprot}.json"

        # Download mmCIF
        url_cif = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.cif'
        os.system(f"curl {url_cif} -o {path_cif}")

        # Download PAE plot
        url_pae = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-predicted_aligned_error_v4.json'
        os.system(f"curl {url_pae} -o {path_pae}")
        
        # Move everthing to a DPAM result directory
        os.system(f'mkdir {config["DPAM_results"]}/{uniprot}')
        os.system(f'cp {path_pdb} {config["DPAM_results"]}/{uniprot}')
        os.system(f'cp {path_cif} {config["DPAM_results"]}/{uniprot}')
        os.system(f'cp {path_pae} {config["DPAM_results"]}/{uniprot}')

        # Change directory. For some reason the tool doesn't work otherwise
        os.chdir(f'{config["DPAM_results"]}/{uniprot}')

        # Run DPAM
        os.system(f'python ../../../../DPAM/DPAM.py {uniprot}.pdb {uniprot}.json {uniprot} . 16 ../../../../DPAM/datadir')

        # Parse DPAM output
        dpam_output = f'{config["DPAM_results"]}/{uniprot}/{uniprot}.finalDPAM.domains'
        with open(dpam_output, 'r') as handle:
            results = handle.read()
        domains = re.findall('[0-9]+-[0-9]+', results)
        
        # Store domains
        uniprot2domains[uniprot] = [list(map(int, domain.split('-')))\
                                    for domain in domains]
        with open(f'{config["data"]}/uniprot2domains.pkl', 'wb') as handle:
            pickle.dump(uniprot2domains, handle)

if __name__ == '__main__':
    main()

    
    
    
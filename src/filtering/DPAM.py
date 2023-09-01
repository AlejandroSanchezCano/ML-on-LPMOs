import os
import re
import pickle
from tqdm import tqdm
from ..config.config_parser import config

def download_pca_cif():

    # Initialize output variable
    uniprot2domains = {}

    # Iterate over proteins to download
    input_proteins = os.listdir(config['AF_all'])[2000:3000]
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
        dpam_output = f'{uniprot}.finalDPAM.domains'
        with open(dpam_output, 'r') as handle:
            results = handle.read()
        domains = re.findall('[0-9]+-[0-9]+', results)
        print(domains)

        # Store domains
        uniprot2domains[uniprot] = [map(int, domain.split('-'))\
                                    for domain in domains]
        with open(f'{config["data"]}/uniprot2domains.pkl', 'wb') as handle:
            pickle.dump(uniprot2domains, handle)

        

    return uniprot2domains


def main():
    '''Program flow.'''

    # Download mmCIF and PAE files
    download_pca_cif()

if __name__ == '__main__':
    main()

    
    
    
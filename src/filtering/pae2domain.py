
import os
import pickle
from tqdm import tqdm
from ..config.config_parser import config


def main():
    '''Program flow.'''
    # Import DPAM predictions
    with open(f'{config["data"]}/uniprot2domains.pkl', 'rb') as handle:
        uniprot2domains = pickle.load(handle)

    # Import uniprot2family
    with open(f'{config["data"]}/uniprot2family.pkl', 'rb') as handle:
        uniprot2family = pickle.load(handle)

    # Iterate over PAE files
    for uniprot in tqdm(uniprot2domains):

        if uniprot2family[uniprot] != 'AA14':
            continue

        # Run PAE to domains
        output_file = f'{config["DPAM_results"]}/{uniprot}/{uniprot}.alternative'
        os.system(f'python -m src.filtering.pae {config["AF_pae"]}/{uniprot}.json --output_file {output_file}')

        # Update output file
        new_domains = []
        with open(output_file, 'r') as file:
            for line in file.readlines():
                line = line.rstrip(',\n').split(',')
                new_domains.append([int(line[0]), int(line[-1])])

        # Update uniprot2domains if the PAE's prediction is better than DPAM
        if uniprot2domains[uniprot] == []:
            continue
        else:
            dpam = uniprot2domains[uniprot][0]

        if dpam[1] - dpam[0] > 250 or dpam[1] - dpam[0] < 100:
            uniprot2domains[uniprot] = new_domains
        
    # Store uniprot2domains
    with open(f'{config["data"]}/uniprot2domains.pkl', 'wb') as handle:
        pickle.dump(uniprot2domains, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    main()
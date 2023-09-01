'''
Crop core.
Use DPAM domain predictions to calculate the enzymatic core of the LPMO 
accessions. The length of the core is plot as an easy way of measuring 
the number of erroneous predictions. Those predictiosn with a core too long
or too short for an LPMO (arounf 165 aa) will be disregarded.
'''

import os
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from ..config.config_parser import config
from .parse_pdb import AlphaFoldStructure

def main():
    '''Program flow.'''

    # Import DPAM predictions
    with open(f'{config["data"]}/uniprot2domains.pkl', 'rb') as handle:
        uniprot2domains = pickle.load(handle)

    # Filter DPAM predictions
    filtered_uniprot2domains = {uniprot : uniprot2domains[uniprot]\
                                for uniprot in os.listdir(config['AF_his1'])}
    
    # Calculate enzymatic domain length
    lengths = []
    for domains in filtered_uniprot2domains.values():
        core_length = domains[0][1] - domains[0][0]
        lengths.append(core_length)
    
    # Plot enzymatic domain length
    ax = sns.histplot(lengths, kde = True)
    ax.set(xlabel = 'Core length')
    plt.savefig(f'{config["plots"]}/core_lengths.png', transparent=True)

    # Save core
    for uniprot, domains in filtered_uniprot2domains.items():
        core_length = domains[0][1] - domains[0][0]
        if core_length > 100 and core_length < 250:
            structure = AlphaFoldStructure(f'{config["AF_his1"]}/{uniprot}.pdb')
            structure.rewrite_range(structure.position[0], domains[0][1] - 1)

if __name__ == '__main__':
    main()
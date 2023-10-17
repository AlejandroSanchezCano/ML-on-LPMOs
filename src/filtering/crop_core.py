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
from tqdm import tqdm
import matplotlib.pyplot as plt
from ..config.config_parser import config
from .parse_pdb import AlphaFoldStructure

def main():
    '''Program flow.'''

    # Import DPAM predictions
    with open(f'{config["data"]}/uniprot2domains.pkl', 'rb') as handle:
        uniprot2domains = pickle.load(handle)

    # Manually examined accessions that DPAM gives no anwser at all of
    uniprot2domains['A0A448SQZ1'] = [(0, 1)]
    uniprot2domains['A0A127PCZ1'] = [(0, 1)]
    uniprot2domains['C5LGY5'] = [(0, 1)]
    uniprot2domains['W2HYN7'] = [(23, 176)]
    uniprot2domains['A0A447QUY5'] = [(0, 1)]
    uniprot2domains['A0A7T8HGE8'] = [(0, 1)]
    uniprot2domains['T0S300'] = [(0, 1)]
    uniprot2domains['W2P9R8'] = [(0, 1)]
    uniprot2domains['W2HZ06'] = [(23, 186)]
    uniprot2domains['A0A4P7KTX1'] = [(0, 1)]
    uniprot2domains['E5ADG7'] = [(0, 1)]
    uniprot2domains['D0NHK1'] = [(0, 1)]
    uniprot2domains['W2FS78'] = [(23, 181)]
    uniprot2domains['W2ZIM1'] = [(0, 1)]
    uniprot2domains['W2ZIM7'] = [(0, 1)]
    uniprot2domains['W2P5K3'] = [(0, 1)]
    uniprot2domains['R7Q626'] = [(0, 1)]
    uniprot2domains['A0A7U2R058'] = [(0, 1)]
    uniprot2domains['A0A5C1RII2'] = [(0, 1)]
    uniprot2domains['W2JPB3'] = [(0, 1)]
    uniprot2domains['G2QWF1'] = [(0, 1)]
    uniprot2domains['T0R0V2'] = [(0, 1)]
    uniprot2domains['W2YX34'] = [(0, 1)]
    uniprot2domains['A0A5K1JRS6'] = [(0, 1)]
    uniprot2domains['A0A7G5JGN7'] = [(0, 1)]

    # Filter DPAM predictions
    filtered_uniprot2domains = {}
    for uniprot in os.listdir(config['AF_his1']):
        uniprot = uniprot.replace('.pdb', '')
        filtered_uniprot2domains[uniprot] = uniprot2domains[uniprot]
    
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
    for uniprot, domains in tqdm(filtered_uniprot2domains.items()):
        core_length = domains[0][1] - domains[0][0]
        if core_length > 100 and core_length < 250:
            structure = AlphaFoldStructure(f'{config["AF_his1"]}/{uniprot}.pdb')
            structure.rewrite_range(
                f"{config['AF_core']}/{structure.id}.pdb",
                (structure.positions[0], domains[0][1] - 1)
                )
    
    ##### Let's try to save AA17 proteins
    # AlphaFold IDs -> corresponding families
    with open(f'{config["data"]}/uniprot2family.pkl', 'rb') as dic:
        dic = pickle.load(dic)

if __name__ == '__main__':
    main()
'''
Plot filtering.

Makes two stacked bar plots (unzoomed and zoomed) about the number of entries:
- GenBank
- UniProt
- AlphaFold
- AlphaFold filtered
- PDB 
'''

import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from ..config.config_parser import config

def main():
    '''Program flow.'''

    # Import data frame with all the LPMOs
    df = pd.read_pickle(f'{config["CAZy_expanded"]}/all.pkl')

    # Import uniprot to families dictionary
    with open(f'{config["data"]}/uniprot2family.pkl', 'rb') as handle:
        uniprot2family = pickle.load(handle)

    # LPMO families
    families = df['Family'].unique()

    # Initialize plotting object
    plot = {family : [] for family in families} 

    # Add total GenBank and UniProt IDs per family
    for family in families:
        sub_df = df[df['Family'] == family]
        plot[family].append(len(sub_df['GenBank'].dropna().drop_duplicates()))
        plot[family].append(len(sub_df['UniProt'].dropna().drop_duplicates()))

    # Add total AlphaFold IDs per family 
    af_proteins = [protein.replace('.pdb', '')\
                   for protein in os.listdir(config['AF_all'])] 
    filtered_uniprot2family = {uniprot:uniprot2family[uniprot]\
                               for uniprot in af_proteins}
    for family in families:
        plot[family].append(list(filtered_uniprot2family.values()).count(family))

    # Add filtered AlphaFold IDs per family
    af_proteins = [protein.replace('.pdb', '')\
                   for protein in os.listdir(config['AF_core'])] 
    filtered_uniprot2family = {uniprot:uniprot2family[uniprot]\
                               for uniprot in af_proteins}

    for family in families:
        plot[family].append(list(filtered_uniprot2family.values()).count(family))
    
    # Add total PDB IDs per family
    pdbs = [0, 42, 43, 2, 7, 2, 4, 1, 4]
    for fam, pdb in zip(plot, pdbs):
        plot[fam].append(pdb)
    
    # Plot unzoomed
    fig, ax = plt.subplots()
    width = 0.5
    legend = [
        'GenBank IDs', 
        'UniProt IDs', 
        'AlphaFold structures', 
        'Filtered AlphaFold structures', 
        'PDB structures'
        ]
    for label, y in zip(legend, zip(*list(plot.values()))):
        p = ax.bar(families, y, width, label=label)
    plt.title("Statistics about the number of entries")
    plt.xlabel('LPMO families')
    plt.ylabel('Number of entries')
    #plt.axhline(1000, color = 'black')
    ax.legend(loc="best")
    plt.savefig(f'{config["plots"]}/unzoomed_filtering.png', transparent=True)
    plt.show()

    # Plot zoomed
    fig, ax = plt.subplots()
    width = 0.5
    legend = [
        'GenBank IDs', 
        'UniProt IDs', 
        'AlphaFold structures', 
        'Filtered AlphaFold structures', 
        'PDB structures'
        ]
    for label, y in zip(legend, zip(*list(plot.values()))):
        p = ax.bar(families, y, width, label=label)
    plt.xlabel('LPMO families')
    plt.ylabel('Number of entries')
    plt.title("Statistics about the number of entries")
    ax.legend(loc="best")
    plt.ylim([0, 1000])
    plt.savefig(f'{config["plots"]}/zoomed_filtering.png', transparent=True)
    plt.show()

    print(plot)

if __name__ == '__main__':
    main()




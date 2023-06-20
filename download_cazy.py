''' Download CAZy info

This file uses the CAZy class to download or parse the information from CAZy

Functions
---------
main()
'''

import pandas as pd
from cazy import CAZy
from variables import CAZY_DATA

def main():
    '''File flow'''
    for family in [9, 10, 11, 13, 14, 15, 16, 17]: # Add AA0 manually
        # Instantiate LPMO family object
        cazyme = CAZy('AA', family)

        # Retrieve data
        print(f'Fetch structural data for family {cazyme.cazyme}')
        structures = cazyme.get_structures()
        print(f'Fetch characterized data for family {cazyme.cazyme}')
        characterized = cazyme.get_characterized()
        print(f'Fetch sequence data for family {cazyme.cazyme}')
        sequences = cazyme.get_sequences()
        print(f'Fetch full sequence data for family {cazyme.cazyme}')
        full_sequences = cazyme.get_full_sequences()
    
        # Store data
        structures.to_pickle(f'{CAZY_DATA}/CAZy_structures/{cazyme.cazyme}')
        characterized.to_pickle(f'{CAZY_DATA}/CAZy_characterized/{cazyme.cazyme}')
        sequences.to_pickle(f'{CAZY_DATA}/CAZy_downloadable_files/{cazyme.cazyme}')
        full_sequences.to_pickle(f'{CAZY_DATA}/CAZy_expanded/{cazyme.cazyme}')

if __name__ == '__main__':
    main()

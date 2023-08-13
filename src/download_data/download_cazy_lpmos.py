''' 
Download CAZy info

This file uses the CAZy class to download or parse the information from CAZy

Functions
---------
main()
'''

from cazy import CAZy
from ..config.config_parser import config

def main():
    '''Program flow'''
    for family in [0, 9, 10, 11, 13, 14, 15, 16, 17]:
        # Instantiate LPMO family object
        cazyme = CAZy('AA', family)

        # Retrieve data
        if family != 0: # AA0 has no structure or characterized HTML page
            print(f'Fetch structural data for family {cazyme.cazyme}')
            structures = cazyme.get_structures()
            print(f'Fetch characterized data for family {cazyme.cazyme}')
        characterized = cazyme.get_characterized()
        print(f'Fetch sequence data for family {cazyme.cazyme}')
        sequences = cazyme.get_sequences()
        print(f'Fetch full sequence data for family {cazyme.cazyme}')
        full_sequences = cazyme.get_full_sequences()
    
        # Store data
        structures.to_pickle(f"{config['CAZy_structures']}/{cazyme.cazyme}")
        characterized.to_pickle(f"{config['CAZy_characterized']}/{cazyme.cazyme}")
        sequences.to_pickle(f"{config['CAZy_downloadable_files']}/{cazyme.cazyme}")
        full_sequences.to_pickle(f"{config['CAZy_expanded']}/{cazyme.cazyme}")

if __name__ == '__main__':
    main()

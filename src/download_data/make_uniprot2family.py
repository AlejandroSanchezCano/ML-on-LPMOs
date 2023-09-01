'''
Make uniprot : family dictionary

Make a dictionary with key the UniProt IDs and values the LPMO families they
correspond to.

Functions
---------
main()
'''

import pickle
import pandas as pd
from ..config.config_parser import config

def main():
    '''Program flow.'''
    df = pd.read_pickle(f"{config['CAZy_expanded']}/all.pkl")
    dic = dict(zip(df['UniProt'], df['Family']))

    with open(f"{config['data']}/uniprot2family.pkl", 'wb') as file:
        pickle.dump(dic, file)

if __name__ == '__main__':
    main()

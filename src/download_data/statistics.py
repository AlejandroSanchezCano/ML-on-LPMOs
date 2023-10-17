from ..config.config_parser import config
import pandas as pd

'''
a, b = 0, 0
for family in ['AA0', 'AA9', 'AA10', 'AA11', 'AA13', 'AA14', 'AA15', 'AA16', 'AA17']:
    df = pd.read_pickle(f'{config["CAZy_expanded"]}/{family}')
    
    df = df.drop_duplicates('UniProt')
    
    s = len(df['GenBank'].dropna())
    d = len(no_dup['GenBank'].dropna())
    a += s
    b += d
    print(family, s, d)
print(a, b)
'''

import pickle
import os

with open(f'{config["data"]}/uniprot2family.pkl', 'rb') as handle:
    dic = pickle.load(handle)

s = 0
for uniprot in os.listdir(config["AF_core"]):
    uniprot = uniprot.replace('.pdb', '')
    if dic[uniprot] == 'AA10' or dic[uniprot] == 'AA9':
        s += 1
print(s)
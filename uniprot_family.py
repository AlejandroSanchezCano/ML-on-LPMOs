import pickle
import pandas as pd
from variables import CAZY_EXPANDED, CAZY_DATA

df = pd.read_pickle(f'{CAZY_EXPANDED}/AA_supreme')
dic = dict(zip(df['UniProt'], df['Family']))

with open(f'{CAZY_DATA}/uniprot_family.pkl', 'wb') as f:
    pickle.dump(dic, f)

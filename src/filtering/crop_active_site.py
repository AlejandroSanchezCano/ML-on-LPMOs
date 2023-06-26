'''
Crop active site by using AlphaFoldStructure's feature of selecting the 
neighbouring residues, those residues that fall in the sphere of given radius 
from a residue as the center point 
'''

import os
from tqdm import tqdm
from variables import AF_HIS1, AF_ACTIVE_SITE
from parse_structures import AlphaFoldStructures

for protein in tqdm(os.listdir(f'{AF_HIS1}')):
    structure = AlphaFoldStructure(f'{AF_HIS1}/{protein}')
    active_site = structure.neighbours(
        center = ('H', 0),
        radius = 20
        )[1]
    structure.rewrite_full(f'{AF_ACTIVE_SITE}/{protein}', active_site)
    

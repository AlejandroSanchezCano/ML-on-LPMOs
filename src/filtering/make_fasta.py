import os
import argparse
from tqdm import tqdm
from variables import AF_CORE, FASTA
from parse_structures import AlphaFoldStructure

# Get protein sequence in FASTA format
fastas = []
for protein in tqdm(os.listdir(AF_CORE)):
    fasta = AlphaFoldStructure(f'{AF_CORE}/{protein}').fasta()
    fastas.append(fasta)

# Make FASTA file
with open(f'{FASTA}/core.fasta', 'w') as file:
    for seq in fastas:
        file.write(seq)


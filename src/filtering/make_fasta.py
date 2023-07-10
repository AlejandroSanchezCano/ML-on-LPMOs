import os
import argparse
from tqdm import tqdm
from variables import AF_FILES, FASTA
from parse_structures import AlphaFoldStructure

# Get protein sequence in FASTA format
fastas = []
for protein in tqdm(os.listdir(AF_FILES)):
    fasta = AlphaFoldStructure(f'{AF_FILES}/{protein}').fasta()
    fastas.append(fasta)

# Make FASTA file
with open(f'{FASTA}/files.fasta', 'w') as file:
    for seq in fastas:
        file.write(seq)


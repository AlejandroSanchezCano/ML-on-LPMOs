import os
from tqdm import tqdm
from variables import AF_HIS1, FASTA
from parse_structures import AlphaFoldStructure

# Get protein sequence in FASTA format
fastas = []
for protein in tqdm(os.listdir(AF_HIS1)):
    fasta = AlphaFoldStructure(f'{AF_HIS1}/{protein}').fasta()
    fastas.append(fasta)

# Make FASTA file
with open(f'{FASTA}/His1.fasta', 'w') as file:
    for seq in fastas:
        file.write(seq)


import os
import argparse
from tqdm import tqdm
from .parse_pdb import AlphaFoldStructure

def handle_arguments() -> argparse.Namespace:
    '''
    Handles the arguments passed via command line.

    Return
    ------
    args : argparse.Namespace
        Arguments
    '''
    
    parser = argparse.ArgumentParser(
        prog = 'make_fasta',
        description = 'Makes FASTA file of AlphaFold structure .pdb files',
        epilog = 'See ya'
        )
    
    # Add arguments
    parser.add_argument(
        '-i', '--input_dir',
        type = str,
        help = 'Input directory containing protein .pdb files'
        )
    
    parser.add_argument(
        '-o', '--output_file',
        type = str,
        help = 'Output path of the FASTA file'
        )

    # Arguments from Namespace object
    return parser.parse_args()


def make_fasta(args : argparse.Namespace):
    '''
    Converts the .pdb structure files from the input directory in 
    AlphaFoldStructure objects. Then it stores their protein sequence in a 
    FASTA file with the UniProt ID as the sequence header.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments
    '''

    # Get protein sequence in FASTA format
    fastas = []
    for protein in tqdm(os.listdir(args.input_dir)):
        fasta = AlphaFoldStructure(f'{args.input_dir}/{protein}').fasta()
        fastas.append(fasta)

    # Make FASTA file
    with open(args.output_file, 'w') as file:
        for seq in fastas:
            file.write(seq)

def main():
    '''Program flow.'''

    # Arguments
    args = handle_arguments()

    # Make fasta
    make_fasta(args)
    

if __name__ == '__main__':
    main()


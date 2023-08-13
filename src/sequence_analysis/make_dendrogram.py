import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import Align

def handle_arguments() -> argparse.Namespace:
    '''
    Handles the arguments passed via command line.

    Return
    ------
    args : argparse.Namespace
        Arguments
    '''
    
    parser = argparse.ArgumentParser(
        prog = 'make_dendrogram',
        description = 'Compute pairwise alignment and build similarity matrix',
        epilog = 'See ya'
        )
    
    # Add arguments
    parser.add_argument(
        '-i', '--input_fasta',
        type = str,
        help = 'Input FASTA file path'
        )
    
    parser.add_argument(
        '-o', '--output_file',
        type = str,
        help = 'Output path of the sequence similarity matrix'
        )

    # Arguments from Namespace object
    return parser.parse_args()

def pairwise_alignment(args : argparse.Namespace):
    '''
    Reads the input FASTA file and computes global pairwise alignment among all
    sequences. The resulting sequence similarity matrix will be stored to be
    used as input for plotting a dendrogram.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments

    Returns
    -------
    None
    '''

    # List all sequences
    with open(args.input_fasta, 'r') as file:
        sequences = file.read().split('\n\n')[:-1]
        seqs = [fasta.split('\n')[1] for fasta in sequences]
        labels = [fasta[1:].split('\n')[0] for fasta in sequences]

    # Configure pairwise alignemnt
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = Align.substitution_matrices.load("PAM70")
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -0.5

    # Compute pairwise alignment -> O(n^2)
    n_seqs = len(seqs)
    matrix = np.zeros([n_seqs, n_seqs])
    for i in tqdm(range(n_seqs)):
        for j in range(i, n_seqs):
            score = aligner.align(seqs[i], seqs[j]).score
            matrix[i][j] = int(score)

    # Triangular matrix -> full matrix
    X = matrix + matrix.T - np.diag(np.diag(matrix))

    # Normalized matrix
    X_norm = X/np.diag(X)[:, None]

    # Save matrix by protein length
    matrix = pd.DataFrame(X_norm, columns = labels, index = labels)
    matrix.to_pickle(args.output_file)

def main():
    '''Progress flow.'''

    # Arguments
    args = handle_arguments()

    # Compute pairwise alignment values
    pairwise_alignment(args)

if __name__ == '__main__':
    main()
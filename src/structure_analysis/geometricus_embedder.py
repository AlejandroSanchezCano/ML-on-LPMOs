import os
import pickle
import argparse
import numpy as np
import geometricus as gm
from variables import CAZY_DATA, STRUCTURE_CLUSTERING

# Argument handling
parser = argparse.ArgumentParser(
                    prog = 'geometricus_embedder',
                    description = 'Runs Geometricus on a set of protein structures',
                    epilog = 'See ya'
                    )

parser.add_argument('--proteins', '-p', 
                    type = str, 
                    default = None,
                    help = 'Set of proteins')
parser.add_argument('--resolution', '-r', 
                    type = float, 
                    default = 1.,
                    help = 'Coarse-grainness of the shapemers')
parser.add_argument('--threads', '-t', 
                    type = int, 
                    default = 4,
                    help = 'Number of threads use in multithreading')
#structures_folder, resolution, n_threads = parser.parse_args()
args = parser.parse_args()
structures_folder = args.proteins
resolution = args.resolution
n_threads = args.threads

# Run Geometricus
invariants, errors = gm.get_invariants_for_structures(
    structures_folder, 
    n_threads = n_threads
    )
model = gm.ShapemerLearn.load()
shapemer_class = gm.Geometricus.from_invariants(
    invariants, 
    model = model, 
    protein_keys = [structure for structure in os.listdir(structures_folder)],
    resolution = resolution
    )
shapemer_count_matrix = shapemer_class.get_count_matrix()

# Normalizationfor protein length
shapemer_sum = np.sum(shapemer_count_matrix, axis = 1)
normalized_matrix = shapemer_count_matrix/shapemer_sum[:, None]

# Save normalized embedding
file_name = structures_folder.split('/')[-1]
np.save(f'{STRUCTURE_CLUSTERING}/{file_name}_resolution{resolution}.npy', normalized_matrix)
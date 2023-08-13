'''
Crop enzymatic core by using AlphaFoldStructure's pLDDT-based feature of 
identifying protein domains. For LPMOs, the enzymatic core is always the first
one. Because this identifying method is not perfect, the distribution of core
length is plotted to discard those cores with abnormally small or big size. 
'''




https://www.google.com/url?sa=i&url=https%3A%2F%2Fwww.biorxiv.org%2Fcontent%2F10.1101%2F2022.09.22.509116v1.full&psig=AOvVaw0ngEtzqHfsHfS8lJ9Wxrsq&ust=1690714033223000&source=images&cd=vfe&opi=89978449&ved=0CBIQjhxqFwoTCMD_9pbfs4ADFQAAAAAdAAAAABAE
https://github.com/tristanic/pae_to_domains


import os
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import matplotlib.pyplot as plt
from parse_structures import AlphaFoldStructure
from variables import AF_HIS1, PLOTS, AF_CORE, CAZY_EXPANDED

# Initialize necessary variables
structures, core_lengths = [], []
# Loop over filtered proteins
for protein in tqdm(os.listdir(AF_HIS1)):
    # Parse AF file
    structure = AlphaFoldStructure(f'{AF_HIS1}/{protein}')
    # Get core  (first domain)
    core = structure.domains(threshold = 70)[0]
    core_length = core[1] - core[0]
    core_lengths.append(core_length)
    # Store core
    if core_length > 120 and core_length < 280:
        structure.rewrite_range(
            output_file = f'{AF_CORE}/{protein}', 
            domains = [core])

# Plot core length distribution to choose cutoffs
ax = sns.histplot(core_lengths, kde = True)
ax.set(xlabel = 'Core length')
ax.axvline(120, color = 'red')
ax.axvline(280, color = 'red')
plt.savefig(f'{PLOTS}/core_lengths.png', transparent=True)

# Filter supreme data frame
supreme_df = pd.read_pickle(f'{CAZY_EXPANDED}/AA_supreme')
core_df = pd.DataFrame({'UniProt': [structure.id for structure in structures]})
filtered = pd.merge(
    left = supreme_df,
    right = core_df,
    how = 'inner',
    on = 'UniProt',
)
core_df.to_pickle(f'{CAZY_EXPANDED}/AA_core')
'''
Plot filtering

Makes two stacked bar plots (unzoomed and zoomed) about the number of entries:
- GenBank
- UniProt
- AlphaFold
- AlphaFold filtered
- PDB 
'''

import os
import pandas as pd
import matplotlib.pyplot as plt
from variables import CAZY_EXPANDED, AF_FILES, PLOTS

'''Program flow.'''

# Import data frames
supreme_df = pd.read_pickle(f'{CAZY_EXPANDED}/AA_supreme')
filtered_df = pd.read_pickle(f'{CAZY_EXPANDED}/AA_filtered')

# LPMO families
families = supreme_df['Family'].unique()

# Initialize plotting object
plot = {family : [] for family in families}

# Add total GenBank and UniProt IDs per family
for family in families:
    sub_dataframe = supreme_df[supreme_df['Family'] == family]
    plot[family].append(len(sub_dataframe['GenBank'].dropna().drop_duplicates()))
    plot[family].append(len(sub_dataframe['UniProt'].dropna().drop_duplicates()))

# Get total AlphaFold IDs per family via inner join
af_proteins = [protein.replace('.pdb', '') for protein in os.listdir(AF_FILES)]
af_proteins_df = pd.DataFrame({'UniProt': af_proteins})
total_af = pd.merge(
    left = supreme_df,
    right = af_proteins_df,
    how = 'inner',
    on = 'UniProt',
)

# Add total AlphaFold IDs per family 
for family in families:
    sub_dataframe = total_af[total_af['Family'] == family]
    plot[family].append(len(sub_dataframe['UniProt'].dropna().drop_duplicates()))

# Add filtered AlphaFold IDs per family
for family in families:
    sub_dataframe = filtered_df[filtered_df['Family'] == family]
    plot[family].append(len(sub_dataframe['UniProt'].dropna().drop_duplicates()))

# Add total PDB IDs per family
pdbs = [0, 42, 43, 2, 7, 2, 4, 1, 4]
for fam, pdb in zip(plot, pdbs):
    plot[fam].append(pdb)

# Plot unzoomed
fig, ax = plt.subplots()
width = 0.5
legend = [
    'GenBank IDs', 
    'UniProt IDs', 
    'AlphaFold structures', 
    'Filtered AlphaFold structures', 
    'PDB structures'
    ]
for label, y in zip(legend, zip(*list(plot.values()))):
    p = ax.bar(families, y, width, label=label)
plt.title("Statistics about the number of entries")
plt.xlabel('LPMO families')
plt.ylabel('Number of entries')
plt.axhline(1000, color = 'black')
ax.legend(loc="best")
plt.savefig(f'{PLOTS}/unzoomed_filtering.png', transparent=True)
plt.show()

# Plot zoomed
fig, ax = plt.subplots()
width = 0.5
legend = [
    'GenBank IDs', 
    'UniProt IDs', 
    'AlphaFold structures', 
    'Filtered AlphaFold structures', 
    'PDB structures'
    ]
for label, y in zip(legend, zip(*list(plot.values()))):
    p = ax.bar(families, y, width, label=label)
plt.xlabel('LPMO families')
plt.ylabel('Number of entries')
plt.title("Statistics about the number of entries")
ax.legend(loc="best")
plt.ylim([0, 1000])
plt.savefig(f'{PLOTS}/zoomed_filtering.png', transparent=True)
plt.show()

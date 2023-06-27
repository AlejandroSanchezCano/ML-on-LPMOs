import os
import umap
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from variables import CAZY_DATA, STRUCTURE_CLUSTERING

# Import shapemer count matrix
normalized_matrix = np.load(f'{STRUCTURE_CLUSTERING}/AF_His1_resolution4.0.npy')

df = pd.DataFrame(normalized_matrix)

# Add AF
ids = [structure.replace('.pdb', '') for structure in os.listdir(f'{CAZY_DATA}/AF_His1')]
df['ID'] = ids

# Add family
with open(f'{CAZY_DATA}/uniprot_family.pkl', 'rb') as dic:
    dic = pickle.load(dic)

df['Family'] = df['ID'].apply(lambda id : dic[id])
print(df)


# UMAP
#reducer = umap.UMAP(metric = "cosine", n_components = 2)
#space = reducer.fit_transform(normalized_matrix)

# PCA
pca = PCA()
space = pca.fit_transform(normalized_matrix)
space = pd.DataFrame(space)
space['ID'] = ids
space['Family'] = space['ID'].apply(lambda id : dic[id])
space = space.sort_values(by = ['Family'], key = lambda x: x.apply(lambda y: int(y[2:])))
print(space)
sns.scatterplot(data = space, x = 0, y = 1, hue= 'Family')
plt.savefig('test.png', transparent = True)

'''


# Scree plot
PC_values = np.arange(pca.n_components_) + 1
plt.plot(PC_values, pca.explained_variance_ratio_, 'ro-', linewidth=2)
plt.title('Scree Plot')
plt.xlabel('Principal Component')
plt.ylabel('Proportion of Variance Explained')
plt.show()
'''
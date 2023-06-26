import numpy as np
from variables import *
r1 = np.load(f'{STRUCTURE_CLUSTERING}/AF_His1_resolution1.0.npy')
r3 = np.load(f'{STRUCTURE_CLUSTERING}/AF_His1_resolution3.0.npy')

c = r1 == r3
print(c)
print(sum(c)[300:])
print(c.shape)
print(r1)
print(r3)

print(sum(r1 == 0))
import numpy as np
import pandas as pd
class SignalP():

    def __init__(self, path):
         self.df = None
         self._parse(SignalP._open(path))
        
    
    def _open(path):
        with open(path, 'r') as file:
            for line in file.readlines():
                if not line.startswith('#'):
                    yield line
    
    def _parse(self, lines):
        matrix = []
        for index, line in enumerate(lines):
            print(index, line)
            elements = line.split()
            if elements[1] == 'OTHER':
                matrix.append(elements[0:4] + [None]*2)
            else:
                selected_elements = [elements[i] for i in (0, 1, 2, 3, 6, 9)]
                matrix.append(selected_elements)
    
        self.df = pd.DataFrame(
            matrix, 
            columns = ['UniProt', 'Prediction', 'SP_prob', 
                       'Other_prob', 'Split', 'Split_prob']
            )
        



#s = SignalP('Data/SignalP/test.txt')

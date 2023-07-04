'''
Parse structures

This file holds hte AlphaFoldStructure class, which parses AlphaFold .pdb files
and retrieves relevant information
'''
import os
import re
import math
import numpy as np

#TODO: - Create an PDB class and a Structure class that is the father of the PDB and AF classes

class AlphaFoldStructure():
    '''
    Class that represents an AlphaFold structure.

    Attributes
    ----------
    path : str
        Absolute path of the .pdb file    
    id : str
        UniProt ID 
    full_seq : str
        Complete protein sequence
    seq : str
        Sequence held in the .pdb file
    pLDDT : np.array
        pLDDT confidence score
    coordinates : np.array
        Per atom coordinates
    coordinates_per_residue : np.array
        Per residue coordinates
    positions : np.array
        Residue index

    Methods
    -------
    __init__
    __str__
    __repr__
    __len__
    __parsing
    domains
    fasta
    neighbours
    rewrite_full
    rewrite_range
    
    Subclass
    --------
    AlphaFoldStructureError
    '''
    
    # Genetic code static variable
    _genetic_code = {
                    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 
                    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
                    }
    
    class AlphaFoldStructureError(Exception):
        '''Custom AlphaFoldStructure exception'''
        pass
    
    def __init__(self, path: str):
        '''
        AlphaFoldStructure constructor. A call is made to the __parsing()
        method to get parse the .pdb file upon instantiation.          

        Parameters
        ----------
        path : str
            Location of the .pdb file.

        '''
        self.path = os.path.abspath(path)
        self.id = path.split('/')[-1].replace('.pdb', '')
        self.full_seq = ''
        self.seq = ''
        self.pLDDT = []
        self.coordinates = []
        self.coordinates_per_residue = []
        self.positions = []
        self.__parsing()

    def __str__(self):
        return f'AlphaFold structure {self.id}'
    
    def __repr__(self):
        return f'AlphaFold structure {self.id}'
    
    def __len__(self):
        return len(self.seq)
    
    def __parsing(self):
        '''
        Parses .pdb file and store the relevant information in the instance 
        attributes.

        Returns
        -------
        None.

        '''
        # Initiate necessary variables
        start_atom_section = True
        coordinates_per_current = []
        
        # Open .pdb file
        with open(self.path, 'r') as file:
            lines = file.readlines()
            
        # Read lines
        for line in lines:
                
            # Get UniProt ID
            if line.startswith('TITLE'):
                self.uniprot = line.split()[-1][1:-1]
                
            # Get sequence
            if line.startswith('SEQRES'):
                peptide_3 = line.split()[4:]
                peptide_1 = ''.join(map(
                    lambda x: AlphaFoldStructure._genetic_code[x]
                    , peptide_3
                    ))
                self.full_seq += peptide_1
            
            # Get pLDDT score, position, coordinates and coordinates_per_residue
            if line.startswith('ATOM'):
                # Regex coordinates and pLDDT
                floats = re.findall(r'-?[0-9]+\.[0-9]+', line)
                *coordinates, _, pLDDT = list(map(float, floats))
                self.coordinates.append(coordinates)
                coordinates_per_current.append(coordinates)
                # Regex position
                current_position = int(re.search(r' A *([0-9]+)', line).group(1))
                # Regex aminoacid
                aminoacid = re.search(r' ([A-Z]{3}) ', line).group(1)
                # Per residue
                if start_atom_section:
                    # Update position flags
                    old_position = current_position
                    start_atom_section = False
                    # Update attributes that don't change within a residue 
                    self.positions.append(current_position)
                    self.pLDDT.append(pLDDT)
                    self.seq += AlphaFoldStructure._genetic_code[aminoacid]

                elif current_position != old_position:
                    # Update position flags
                    old_position = current_position
                    # Update attributes that don't change within a residue 
                    self.positions.append(current_position)
                    self.pLDDT.append(pLDDT)
                    self.seq += AlphaFoldStructure._genetic_code[aminoacid]
                    # Update that change within a residue 
                    self.coordinates_per_residue.append(np.mean(
                            coordinates_per_current[:-1],
                            axis = 0
                            ))
                    coordinates_per_current = [coordinates_per_current[-1]]
            
            # Add coordinates_per_current of final aminoacid
            if line.startswith('TER'):
                self.coordinates_per_residue.append(np.mean(
                        coordinates_per_current,
                        axis = 0
                        ))
                        
        # Convert to numpy arrays
        self.positions = np.array(self.positions)
        self.pLDDT = np.array(self.pLDDT)
        self.coordinates = np.matrix(self.coordinates)
        self.coordinates_per_residue = np.matrix(self.coordinates_per_residue)
    
    def domains(self, threshold: int = 70, consecutive: int = 3) -> list:
        '''
        Identify domains in the .pdb file based on pLDDT scores. When all 
        residues in a strecth of consecutive residues has a pLDDT score that 
        crosses a predetermined threshold, this is detected as the beginning or
        end of a domain. Domain recognition is based upon the experience that 
        domains in LPMOs (high pLDDT) are separated by stretches of low pLDDT
        values.
        While .pdb files indexes are 1-based, the indexes this method works 
        with are 0-based.

        Parameters
        ----------
        threshold : int, optional
            Threshold that marks the beginning or end of a domain. 
            The default is 70.
        consecutive : int, optional
            Consecutive number of aminoacids that must have a score above or 
            below the threshold to trigger the domain sensing. The default is 3.

        Returns
        -------
        list
            Domain regions -> [(3, 100), (120, 232)]

        '''
        # Initialize necessary variables
        cleaned_seq = []
        new_domain = True
        start_domain = 'not initialized'
        
        for index, pLDDT in enumerate(self.pLDDT[:-consecutive]):
            # Start of high-confidence domain 
            if new_domain and all(self.pLDDT[index : index + consecutive] >= threshold):
                start_domain = self.positions[index] - 1
                new_domain = False
            # End of high-confidence domain 
            elif not new_domain and all(self.pLDDT[index : index + consecutive] <= threshold):
                end_domain = self.positions[index] - 1
                cleaned_seq.append((start_domain, end_domain))
                new_domain = True

        # Account for failed structures
        if cleaned_seq == []:
            # Account for structures with full core
            if np.mean(self.pLDDT) > 95:
                return [(0, len(self) - 1)]
            # Account for structures with no CBM
            elif start_domain != 'not initialized':
                return [(start_domain, len(self) - 1)]
            else:
                return []
        else:
            return cleaned_seq
        
    def fasta(self) -> str:
        '''
        Returns sequence in fasta format

        Returns
        -------
        str
            Sequence in fasta format.

        '''
        return f'>{self.uniprot}\n{self.seq}\n\n'
    
    def neighbours(self,
                   center: tuple[str, int] = ('H', 0),
                   radius: int = 20
                   ) -> tuple[str,list[int]]:
        '''
        Find the residues that are spatially close to a residue. The default is 
        His1. This is done by calculating the Euclidean distance between the
        center residue and all the other residues in the protein.

        Parameters
        ----------
        center : tuple[str, int], optional
            Residue and its index that will act as the center of the sphere.
            The default is ('H', 0).
        radius : int, optional
            Radius of the sphere. The default is 20.

        Raises
        ------
        AlphaFoldStructureError
            When if the specified residue that acts as the center of the sphere 
            is not present in the sequence.

        Returns
        -------
        sequence : str
            Protein sequence with residues as uppercase if they neighbour the
            center residue.
        
        positions : list[int]
            Index of the residues that neighbor the center residue.
        '''
        # Raise error if the specified residue that acts as the center of
        # the sphere is not present in the sequence
        aminoacid, position = center
        if self.seq[position] != aminoacid:
            raise AlphaFoldStructure.AlphaFoldStructureError(
                f'The aminoacid {aminoacid} is not present in the position {position} in the sequence: {self.seq}'
                )
        # Get residues close to the center residue
        else:           
            sequence = ''
            positions = []
            center = self.coordinates_per_residue[0, :].tolist()[0]
            for name, coordinates, position in zip(
                    self.seq, 
                    self.coordinates_per_residue,
                    self.positions):
                residue = coordinates.tolist()[0] 
                distance = math.dist(center, residue)
                if distance <= radius:
                    sequence += name.upper()
                    positions.append(position)
                else:
                    sequence += name.lower()
                    
        return sequence, positions
    
    def rewrite_full(self, output_file: str, positions: list[int]) -> None:
        '''
        Rewrites the .pdb file with a subset of the residues

        Parameters
        ----------
        output_file : str
            Path to store new .pdb file version.
        positions : list[int]
            0-based indexes of the residues to keep.

        Returns
        -------
        None
        '''
        with open(self.path, 'r') as in_file, open(output_file, 'w') as out_file:
            for line in in_file:
                if line.startswith('ATOM'):
                    # Regex position
                    position = -1 + int(re.search(r' A *([0-9]+)', line).group(1))             
                    if position in positions:
                        out_file.write(line)
                else:
                    out_file.write(line)
                    
    def rewrite_range(self, 
                      output_file: str, 
                      domains: list[tuple[int, int]]
                      ) -> None:
        '''
        Rewrites the .pdb file with the subset of the residues specified by the
        domains range.

        Parameters
        ----------
        output_file : str
            Path to store new .pdb file version.
        domains : list[tuple[int, int]]
            List with the initial and final positions (as a tuple) of a range
            of residues (like a domain).

        Returns
        -------
        None.
        
        '''
        # Transform domains (list of tuples of ints) in positions (list of ints)
        positions = []
        for start, end in domains:
            positions += list(range(start, end))
        
        # Call regular rewrite
        self.rewrite_full(output_file, set(positions))
        
        
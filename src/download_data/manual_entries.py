'''
Add manual entries.

This file holds the ManualEntries class, which shares functionalities with 
the CAZy class with the intention to add new entries automatically provided
their UniProt ID.
TO DO: see what happens with failed IDs (UniProt ids without genbank id like 
G3XAP7)

Classes
-------
ManualEntries

'''

import re
import pandas as pd
from tqdm import tqdm
from bioservices import UniProt
from collections import defaultdict
from Bio import Entrez, SeqIO, GenBank


class ManualEntries():
    '''
    Class that represents a CAZy family or subfamily.

    Attributes
    ----------
    families : list[str]
        List of families where each manual entry belongs to .
    domains : list[str]
        List of life domains (i.e. Eukaryota, Bacteria, Virus) where each 
        manual entry belongs to.
    uniprot_ids : list[str]
        List of UniProt IDs to incorporate.
    df : pd.DataFrame 
        Info specified and retrieved from the manual entries in a pandas
        data frame format.


    Methods
    -------
    __init__
    __repr__
    __str__
    add_genbank
    add_species
    add_protein_fasta
    add_dna_fasta
    add_pdb
    add_manual_entries

    Class
    --------
    ManualEntriesError
    
    '''  

    class ManualEntriesError(Exception):
        '''Custom ManuealEntries exception'''
        pass

    def __init__(self, 
                 families : list[str], 
                 domains : list[str], 
                 uniprot_ids: list[str]
                 ):
        '''
        MabualEntries class constructor.
        The inputs are checked to be of the same length.

        Parameters
        ----------
        families : list[str]
            List of families where each manual entry belongs to .
        domains : list[str]
            List of life domains (i.e. Eukaryota, Bacteria, Virus) where each 
            manual entry belongs to.
        uniprot_ids : list[str]
            List of UniProt IDs to incorporate.
        '''       
        if len(families)* 3 != len(families) + len(domains) + len(uniprot_ids):
            raise ManualEntries.ManualEntriesError('Inputs are not of the same length')
        else:
            self.families = families
            self.domains = domains
            self.uniprot_ids = uniprot_ids
            self.df = pd.DataFrame({
                'Family': self.families,
                'Domain': self.domains, 
                'UniProt': self.uniprot_ids
                })
            
    def __str__(self):
        return f'{len(self.families)} new manual entries'
    
    def __repr__(self):
        return f'{len(self.families)} new manual entries'
            
    def add_genbank(self):
        '''
        Maps UniProt IDs to GenBank IDs and adds them to the data frame in
        construction.

        Returns
        -------
        None
        '''
        # Uniprot class instantiation
        u = UniProt()

        # Map UniProt - GenBank
        uniprot_results = u.mapping(
            fr = 'UniProtKB_AC-ID', 
            to = 'EMBL-GenBank-DDBJ_CDS', 
            query = self.uniprot_ids
            )
        uniprot2genbank = defaultdict(list)
        for result in uniprot_results['results']:
            uniprot2genbank[result['from']] = result['to']
        
        # Add to self.df
        self.df['GenBank'] = uniprot2genbank.values() 
        self.df = self.df.explode('GenBank')
        self.genbank_ids = self.df['GenBank'].to_list()
    
    def add_species(self):
        '''
        Fetches the GenBank file of each manual entry, parses it to find the
        species it belongs to and adds it to the data frame in construction.
        
        Returns
        -------
        None
        '''
        # NCBI requests require stating an email
        Entrez.email = ""

        # Fetch GenBank file
        self.species = []
        for genbank_id in tqdm(self.genbank_ids):
            handle = Entrez.efetch(
                db="protein", 
                id=genbank_id, 
                rettype="gb", 
                retmode="text"
                )
            entry = GenBank.read(handle)
            
            # Parse GenBank file
            for feature in entry.features:
                for qualifier in feature.qualifiers:
                    key = qualifier.key
                    if key == '/organism=':
                        strain = eval(qualifier.value)
                    elif key == '/strain=':
                        strain += ' ' + eval(qualifier.value)
            self.species.append(strain)
        
        # Add to self.df
        self.df['Species'] = self.species
        
    def add_protein_fasta(self):
        '''
        Add protein sequences to sequences data frame via GenBank ID - sequence
        mapping

        Returns
        -------
        None
        '''

        # NCBI requests require stating an email
        Entrez.email = ""

        # Input sequences
        sequences = self.genbank_ids
        
        # Protein sequences
        handle_protein = Entrez.efetch(
            db = "protein", 
            id = sequences, 
            rettype = "fasta", 
            retmode = "text"
            )
        SeqIO_protein = SeqIO.parse(handle_protein, 'fasta')
        
        # Make protein FASTA format
        counter = -1
        protein_sequences = []
        with tqdm(total = len(sequences)) as pbar: # The progress bar needs to be regulated manually
            for protein in SeqIO_protein:
                counter += 1
                # Repeat process with the next GenBank protein ID if there's no match
                while sequences[counter] != protein.id:
                    counter += 1
                    pbar.update(1)
                    protein_sequences.append(None)
                    # Security measure to avoid infinite loop
                    if counter == len(sequences):
                        raise Exception('A forever while loop has been stopped')
                # If there is a match, add FASTA
                else:
                    fasta_protein = f">{protein.description}\n{protein.seq}\n\n"
                    protein_sequences.append(fasta_protein)
                    pbar.update(1)
        
        # Add protein sequences to data frame
        self.df['Protein sequence'] = protein_sequences
        
        return self.df   
    
    def add_dna_fasta(self):
        '''
        Add DNA sequences to sequences data frame via GenBank ID - sequence
        mapping

        Returns
        -------
        None
        '''

        # NCBI requests require stating an email
        Entrez.email = ""

        # Input sequences
        sequences = self.genbank_ids

        # DNA sequences
        handle_dna = Entrez.efetch(
            db="protein", 
            id=sequences, 
            rettype='fasta_cds_na', 
            retmode="text"
            )
        SeqIO_dna = SeqIO.parse(handle_dna, 'fasta')

        def protein_id_from_dna_fasta(dna):
            '''
            Extract the protein ID from the fasta description. 

            Parameters
            ----------
            description : str
                Sequence description.

            Returns
            -------
            str
                protein ID.

            '''
            match = re.search(r'protein_id=([A-Z0-9_\.]+)\]', dna.description)
            return match.group(1)

        # Make DNA FASTA format
        counter = -1
        dna_sequences = []
        with tqdm(total = len(sequences)) as pbar: # The progress bar needs to be regulated manually
            for dna in SeqIO_dna:
                counter += 1
                # Repeat process with the next GenBank protein ID if there's no match
                while sequences[counter] != protein_id_from_dna_fasta(dna):
                    counter += 1
                    pbar.update(1)
                    dna_sequences.append(None)
                    # Security measure to avoid infinite loop
                    if counter == len(sequences):
                        raise Exception('A forever while loop has been stopped')
                # If there is a match, add FASTA
                else:
                    fasta_dna = f">{dna.description}\n{dna.seq}\n\n"
                    dna_sequences.append(fasta_dna)
                    pbar.update(1)

        # Add DNA sequences to data frame
        self.df['DNA sequence'] = dna_sequences
        
        return self.df

    def add_pdb(self):
        '''
        Adds UniProt IDs to sequences data frame via GenBank - UniProt mapping

        Returns
        -------
        None
        '''
        
        # Uniprot class instantiation
        u = UniProt()

        # Mapping UniProt - PDB
        uniprot_results = u.mapping(
            fr = 'UniProtKB_AC-ID', 
            to = 'PDB', 
            query = self.uniprot_ids
            )
        
        # Parse mapping results
        uniprot2pdb = defaultdict(list)
        for result in uniprot_results['results']:
            uniprot_id = result['from']
            pdb_id = result['to']
            uniprot2pdb[uniprot_id].append(pdb_id)
        
        # Add PDB IDs to data frame
        has_pdb = lambda uniprot: None if uniprot not in uniprot2pdb\
            else uniprot2pdb[uniprot]
        self.df['PDB'] = self.df['UniProt'].apply(has_pdb)

    def add_manual_entries(self, supreme_df):
        '''
        Parameters
        ----------
        supreme_df : pd.DataFrame
            Concatenated data frame of all the individual LPMO families data
            frames.
        Returns
        -------
        Data frame with the added information of the manual entries.
        '''
        return pd.concat([supreme_df, self.df], ignore_index = True)
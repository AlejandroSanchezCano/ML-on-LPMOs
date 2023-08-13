'''
InterPro results parser

This file holds the SignalP class, which uses the InterPro REST API to make
calls of the LPMO-associated entry IDs in order to select the domain range of
each entry ID per protein.

Class
-----
InterPro
'''

import math
import time
import pickle
import requests
from tqdm import tqdm
from typing import Tuple
from collections import defaultdict
from ..config.config_parser import config


class InterPro():
    '''
    Class that represents an entry ID in InterPro.
    This class is able to parse InterPro API responses exclusively to retrive
    the domain range of the entry ID in the entry proteins, but it could be 
    extended to the point of becoming an InterPro REST API parser. 

    Attributes
    ----------
    id : str
        Entry ID. Some entry ID examples can be IPR005103 (InterPro), cd21137
        (CDD) or G3DSA:2.70.50.70 (CATH-Gene3D). The ID format varies depending
        on the member database supported by InterPro that it belongs to.

    Methods
    -------
    __init__
    get_domains
    '''

    def __init__(self, name : str):
        '''
        SignalP constructor

        Parameters
        ----------
        name : str
            Entry ID
        '''
        self.id = name
    
    def get_domains(self) -> dict[str, Tuple[str, str]]:
        '''
        Uses the self.id attribute to make two requests to the InterPro API.
        The fist one has the objective to get the number of total entries, 
        which will be use in the second request as input for the progress bar.
        The json outputed by the REST API will be parsed and the domain info
        stored.
        Each call to the InterPro API is interspaced by 1 second of
        waiting time and will expect 200 protein results (batch_size).

        Parameters
        ----------
        None.

        Returns
        -------
        domains : dict[str, Tuple[str, str]]
            UniProt ID as keys and start and end of the domain/region as a 
            tuple.
        '''

        # Handle database
        if self.id.startswith('IPR'):
            db = 'interpro'
        elif self.id.startswith('G3DSA'):
            db = 'Cathgene3d'
        elif self.id.startswith('cd'):
            db = 'Cdd'

        # Get number of total entries
        url_core = 'https://www.ebi.ac.uk/interpro/api/protein/uniprot/entry'
        url_number_total_entries = f'{db}/{self.id}?page_size=1'
        url = url_core + url_number_total_entries
        request = requests.get(url)
        json = request.json()
        total = json['count']

        # Get domains
        domains = defaultdict(list)
        batch_size = 200
        url_domains = f'{db}/{self.id}?page_size={batch_size}'
        url = url_core + url_domains
        for _ in tqdm(range(math.ceil(total/batch_size))):
            request = requests.get(url)
            json = request.json()
            for accession in json['results']:
                uniprot = accession['metadata']['accession']
                
                # Error if accession['entries'] is a list of multiple elements
                entries = accession['entries']
                if len(entries) > 1:
                    raise Exception(f'UniProt accession {uniprot} has multiple entries: {entries}')
                
                # Take into account entries with multiple entry_protein_locations
                for entry_protein_location in entries[0]['entry_protein_locations']:
                
                    # Spaced domains result in multiple fragments
                    fragments = entry_protein_location['fragments']
                    # Get start and end of domain per fragment
                    domains[uniprot] += [(fragment['start'],  fragment['end'])\
                                          for fragment in fragments]
                    
            # Next batch
            url = json['next']
            time.sleep(1)

        return domains

def main():
    '''Program flow.'''

    # nothing for AA14 and AA17
    aa9 = InterPro('IPR005103').get_domains()
    aa10 = InterPro('IPR004302').get_domains()
    aa11 = InterPro('G3DSA:2.70.50.70').get_domains()
    aa13 = InterPro('cd21137').get_domains()

    # Unify dictionaries
    domains = aa13 | aa11 | aa10 | aa9
    
    # Store merged dictionary
    with open(f'{config["InterPro"]}/interpro.pkl', 'wb') as handle:
        pickle.dump(domains, handle)

if __name__ == '__main__':
    main()

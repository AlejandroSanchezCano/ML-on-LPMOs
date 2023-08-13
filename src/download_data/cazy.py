''' 
CAZy parser

This file holds the CAZy class, which serves as a perser and an API to connect
with CAZy database. Some additional features are also available such as 
retrieving DNA and protein FASTA sequences and UniProt ID from the GenBank ID, 
as well as mapping the UniProt IDs to the PDB.

Classes
-------
CAZy
'''

import re
import bs4
import requests
import pandas as pd
from tqdm import tqdm
from Bio import Entrez, SeqIO
from bioservices import UniProt
from collections import defaultdict
from typing import Optional, Union, Callable, Any, Generator

class CAZy():
    '''
    Class that represents a CAZy family or subfamily.

    Attributes
    ----------
    enzyme_class : str
        CAZyme class.
    family : Union[str, int]
        CAZyme family.
    subfamily : Union[str, int], optional
        CAZyme subfamily. The default is None.
    cazyme : str
        CAZyme string with the information of enzyme class, family and
        subfamily as CAZy works with -> e.g. GH5_2.

    Methods
    -------
    __init__
    __repr__
    __str__
    __get_soup_class
    __validate_enzyme_class
    __validate_sub_family
    __get_soup_structure
    __parse_table_row
    __parse_HTML
    get_structures
    __get_soup_characterized
    get_characterized
    __parse_txt
    get_sequences
    __add_protein_fasta
    __add_dna_fasta
    __add_uniprot
    __add_pdb
    get_full_sequences
    
    Class
    --------
    CAZyError
    
    '''
    
    # CAZy class dict
    __cazymes = {
        'GH': 'Glycoside-Hydrolases',
        'GT': 'GlycosylTransferases',
        'PL': 'Polysaccharide-Lyases',
        'CA': 'Carbohydrate-Esterases',
        'AA': 'Auxiliary-Activities',
        'CBM': 'Carbohydrate-Binding-Modules'   
    }
    
    class CAZyError(Exception):
        '''Custom CAZy exception'''
        pass

    def __init__(self, 
                 enzyme_class: str, 
                 family: Union[str, int], 
                 subfamily: Optional[Union[str, int]] = None
                 ):
        '''
        CAZy class constructor.
        A call to the validation methods is made to ensure validity of input
        attribute data.

        Parameters
        ----------
        enzyme_class : str
            CAZyme class.
        family : Union[str, int]
            CAZyme family.
        subfamily : Union[str, int], optional
            CAZyme subfamily. The default is None.

        '''       
        # Assign attributes
        self.enzyme_class = enzyme_class
        self.family = str(family)
        if subfamily:
            self.subfamily = str(subfamily)
            self.cazyme = f'{self.enzyme_class}{self.family}_{self.subfamily}'
        else:
            self.cazyme = f'{self.enzyme_class}{self.family}'
        
        # Validate input
        self.__validate_enzyme_class()
        self.__validate_sub_family()
        
    def __repr__(self):
        return f'CAZy object of CAZyme {self.cazyme}'
    
    def __str__(self):
        return f'CAZy object of CAZyme {self.cazyme}'
    
    def __get_soup_class(self) -> bs4.BeautifulSoup:
        '''
        Perform request of the CAZyme class main HTML page and return it as
        a beautiful soup object to be processed and parsed.    

        Raises
        ------
        CAZyError
            When the number of requests have exceeded what the web server's
            resources can handle and a HTTP 429 error is returned from the
            request.

        Returns
        -------
        soup : BeautifulSoup
            Souped html page.

        '''
        # Get souped HTML text
        enzyme_class = CAZy.__cazymes[self.enzyme_class]
        url = f"http://www.cazy.org/{enzyme_class}.html"
        soup_class = bs4.BeautifulSoup(
            requests.get(url).text, features= 'lxml'
            )

        # 'Too Many Requests' error
        if 'Too Many Requests' in soup_class.text:
            raise CAZy.CAZyError('Error from www.CAZy.org: Status 429, too many requests (try again soon)')
        
        return soup_class
    
    def __valid_input(self) -> tuple[list[str], dict[str:list[int]]]:
        '''
        Parses CAZyme class HTML page to get a list of valid input families and
        subfamilies. This is done by finding the 'a' (link) HTML elements, 
        filtering for the ones that are in class_family(_subfamily) format, get
        unique values, and sort them.

        Returns
        -------
        sorted_unique_families : list[str]
            Valid families -> ['1', '2', '3']
        
        subfamilies_available : list[str]
            Valid subfamilies -> {'GH12':['1', '2', '3']}
        
        '''
        # Get links in soup
        soup_class = self.__get_soup_class()
        links = soup_class.findAll('a')
        
        # Find families and subfamilies available in the links
        families_available = [] 
        subfamilies_available = defaultdict(list)
        pattern = rf'{self.enzyme_class}([0-9]*)_?([0-9]*)'
        for link in links:
            if 'href' in link.attrs:
                match = re.findall(pattern, link.attrs['href'])
                if match:
                    family, subfamily = match[0]
                    if subfamily:
                        class_family = f'{self.enzyme_class}{family}'
                        subfamilies_available[class_family].append(subfamily)
                    else:
                        families_available.append(family)
                        
        # Unique families
        unique_families = set(families_available)
        # Sort unique families
        sorted_unique_families = sorted(unique_families, key = int)
        
        return sorted_unique_families, subfamilies_available

    def __validate_enzyme_class(self) -> None:
        '''
        Validates input enzyme class: GH, GT, PL, CA, AA, CBM

        Raises
        ------
        CAZyError
            When invalid enzyme class.

        Returns
        -------
        None

        '''
        if self.enzyme_class not in CAZy.__cazymes:
            valid_enzyme_classes = ', '.join(CAZy.__cazymes.keys())
            raise CAZy.CAZyError(f'{self.enzyme_class} is an invalid enzyme class, choose from: {valid_enzyme_classes}')

    def __validate_sub_family(self) -> None:
        '''
        Validates input family and subfamily

        Raises
        ------
        CAZyError
            when invalid family or subfamily.

        Returns
        -------
        None
        
        '''
        # Valid input
        families, subfamilies = self.__valid_input()
        
        # Invalid family
        if self.family not in families:
            options = ", ".join(families)
            raise CAZy.CAZyError(f'{self.family} is an invalid family for {self.enzyme_class}, choose from: {options}')
            
        # Has subfamily attributed been specified?
        elif hasattr(self, 'subfamily'):
            class_family = self.enzyme_class + self.family
            # No subfamily
            if class_family not in subfamilies:
                raise CAZy.CAZyError(f'Family {class_family} has no subfamilies')
            # Invalid subfamily
            elif self.subfamily not in subfamilies[class_family]:
                subfamilies = subfamilies[class_family]
                options = ", ".join(subfamilies)
                raise CAZy.CAZyError(f'{self.subfamily} is an invalid subfamily for {class_family}, choose from: {options}')

###############################################################################
##########################    STRUCTURE CAZYMES    ############################
###############################################################################

    def __get_soup_structure(self) -> bs4.BeautifulSoup:
        '''
        Perform request of the CAZy family structure HTML page and return it as
        a beautiful soup object to be processed and parsed.  

        Raises
        ------
        CAZyError
            When the number of requests have exceeded what the web server's
            resources can handle and a HTTP 429 error is returned from the
            request.

        Returns
        -------
        soup_structure : BeautifulSoup
            Souped html page.

        '''
        
        # Get souped HTML text
        url = f"http://www.cazy.org/{self.cazyme}_structure.html"
        soup_structure = bs4.BeautifulSoup(
            requests.get(url).text, features= 'lxml'
            )

        # 'Too Many Requests' error
        if 'Too Many Requests' in soup_structure.text:
            raise CAZy.CAZyError('Error from www.CAZy.org: Status 429, too many requests (try again soon)')
        
        return soup_structure
    
    def __parse_table_row(tr: bs4.element.Tag) -> list[Union[str, list, tuple]]:
        '''
        Parses individual row (tr) objects. This function identifies if there 
        is a subtable in the cell (td) and calls itself recursively. If not, 
        retrieves the elements that the cell contains (as list or tuple if 
        multiple elements or as str if only 1)
        
        Parameters
        ----------
        tr : bs4.element.Tag
            Table row as a bs4 object.

        Returns
        -------
        list[Union[str, list, tuple]]
            Table row as a list.

        '''
        # Instantiate variable that will hold the row's cell contents
        row = []
        
        # Loop over all the table data cell (td) in the row (tr)
        tds = tr.findAll('td', recursive = False)
        for td in tds:
            # Detect table in the td
            sub_table = td.find('table')
            if sub_table:
                sub_result = []
                # Find trs in subtable
                sub_trs = sub_table.findAll('tr')
                # Find tds in trs
                for sub_tr in sub_trs:
                    # Apply same function
                    sub_result.append(CAZy.__parse_table_row(sub_tr))
                
                # Adjust how result look like depending on the number of sub_trs
                if len(sub_result) > 1:
                    row += list(zip(*sub_result))
                else:
                    row += sub_result[0]
                    
            # Parse normally if there's no table in the td        
            else:
                # 1. Extract text if several links/entries
                if len(td.findAll('a')) > 1:
                    text = [td_link.text for td_link in td.findAll('a')]
                # 2. Some GenBank IDs do not have a link
                elif re.match(genbank_pattern := r'[a-zA-Z]+_?[0-9]+\.[0-9]{1}', td.text):
                    text = re.findall(genbank_pattern, td.text)
                # 3. Remove leading/trailing characters
                else:
                    text = td.get_text(strip = True).rstrip(' #')
                
                # Append the td's text to the row
                row.append(text)
                
        return row
    
    def __parse_HTML(self, soup: bs4.BeautifulSoup) -> pd.DataFrame:
        '''
        Parses a souped HTML table and retrieves it in a pandas data frame
        format.

        Parameters
        ----------
        soup : bs4.BeautifulSoup
            Souped html page.

        Returns
        -------
        df : pd.DataFrame
            HTML table as data frame.

        '''
        # Get main HTML table
        HTML_table = soup.find(attrs = {"id":'pos_onglet', "class":'listing'})
        # HTML table as matrix
        matrix = []
        # Sometimes royaume tr is not the first tr
        royaume = None
        # Loop over the table rows
        trs = HTML_table.findAll('tr', recursive = False)
        for tr in trs:
            
            # Find kingdom
            if tr.get('class') == ['royaume']:
                royaume = tr.get_text(strip = True)
                continue
            
            if royaume:
                # Find column names (disregarding empty tds)
                if tr.get('id') == 'line_titre' and\
                    len(tr.findAll('td', recursive = False)) > 1:
                        titles = ['Domain'] + CAZy.__parse_table_row(tr)
                        # Remove trailing empty elements that sometimes appear
                        titles = [title for title in titles if title]
                        
                # Parse content rows
                else:
                    row = [royaume] + CAZy.__parse_table_row(tr)[:len(titles) - 1]
                    matrix.append(row)
                
        # Matrix -> data frame
        df = pd.DataFrame(matrix, columns = titles)     
        
        # The cells that come from parsing a subtable are retrieved as tuples
        # By identifying the columns that have tuples we can explode them.
        is_tuple = lambda cell: isinstance(cell, tuple)
        column_has_tuple = lambda column: column.apply(is_tuple)
        explode = [column for column in df if any(column_has_tuple(df[column]))]
        
        df = df.explode(explode) if explode != [] else df
        
        return df
    
    def get_structures(self) -> pd.DataFrame:
        '''
        Returns structural information of a CAZy (sub)family

        Returns
        -------
        df : pd.DataFrame
            CAZy (sub)family's structure information

        '''
        soup = self.__get_soup_structure()
        df = self.__parse_HTML(soup)
        return df
         
###############################################################################
########################    CHARACTERIZED CAZYMES    ##########################
###############################################################################

    def __get_soup_characterized(self) -> bs4.BeautifulSoup:
        '''
        Perform request of the CAZy characterized structure HTML page and return it as
        a beautiful soup object to be processed and parsed. 

        Raises
        ------
        CAZyError
            When the number of requests have exceeded what the web server's
            resources can handle and a HTTP 429 error is returned from the
            request.

        Returns
        -------
        soup_characterized : BeautifulSoup
            Souped html page.

        '''
        # Get souped HTML text
        url = f"http://www.cazy.org/{self.cazyme}_characterized.html"
        soup_characterized = bs4.BeautifulSoup(
            requests.get(url).text, features= 'lxml'
            )

        # 'Too Many Requests' error
        if 'Too Many Requests' in soup_characterized.text:
            raise CAZy.CAZyError('Error from www.CAZy.org: Status 429, too many requests (try again soon)')
        
        return soup_characterized
    
    def get_characterized(self) -> pd.DataFrame:
        '''
        Returns information of a CAZy (sub)family's entries that have been
        experimentally characterized

        Returns
        -------
        df : pd.DataFrame
            CAZy (sub)family's characterized information

        '''
        soup = self.__get_soup_characterized()
        df = self.__parse_HTML(soup)
        return df

###############################################################################
##########################    SEQUENCES CAZYMES    ############################
###############################################################################
    
    def __parse_txt(lines):
        names = ['Family', 'Domain', 'Species', 'GenBank']
        columns = zip(*[line.rstrip().split('\t') for line in lines])
        df = pd.DataFrame({name:column for name, column in zip(names, columns)})
        
        return df
        
    def get_sequences(self, local_path = None, download_path =  None) -> pd.DataFrame:
        'only yse download if local_path is nule'
        if local_path and download_path:
            raise CAZy.CAZyError("Only use download_path 'download_path' argument when the 'local_path' argument is not specified ")
        if local_path:
            with open(local_path, 'r') as file:
                lines = file.readlines()
                df = CAZy.__parse_txt(lines)
        else:
            url = f"http://www.cazy.org/IMG/cazy_data/{self.cazyme}.txt"
            text = requests.get(url).text
            lines = text.split('\n')[:-1]
            df = CAZy.__parse_txt(lines)
            if download_path:
                with open(download_path, 'w') as file:
                    file.write(text)
                    
        return df
            
###############################################################################
#######################   EXPAND SEQUENCES CAZYMES    #########################
###############################################################################

    def __add_protein_fasta(df: pd.DataFrame) -> pd.DataFrame:
        '''
        Add protein sequences to sequences data frame via GenBank ID - sequence
        mapping. Chunking is necessary to avoid weird situations like some
        IDs that should map but don't and then they generate bugs.


        Parameters
        ----------
        df : pd.DataFrame
            Sequences data frame.

        Raises
        ------
        CAZyError
            To stop an infinite loop.

        Returns
        -------
        df : pd.DataFrame
            Sequences data frame + protein sequences.

        '''
        # NCBI requests require stating an email
        Entrez.email = ""
        # Input sequences
        genbank = df['GenBank']
        # Initialize output variable
        protein_sequences = []
        # Chunk input
        chunked_genbank = [genbank[i:i + 100] for i in range(0, len(genbank), 100)]
        for chunk in tqdm(chunked_genbank):
            
            chunk = chunk.to_list()
            # Protein sequences
            handle_protein = Entrez.efetch(
                db = "protein", 
                id = chunk, 
                rettype = "fasta", 
                retmode = "text"
                )
            SeqIO_protein = SeqIO.parse(handle_protein, 'fasta')
        
            # Make protein FASTA format
            counter = -1

            for protein in SeqIO_protein:
                counter += 1
                # Repeat process with the next GenBank protein ID if there's no match
                while chunk[counter] != protein.id:
                    counter += 1
                    protein_sequences.append(None)
                    # Security measure to avoid infinite loop
                    if counter == len(chunk):
                        raise CAZy.CAZyError('A forever while loop has been stopped')
                # If there is a match, add FASTA
                else:
                    fasta_protein = f">{protein.description}\n{protein.seq}\n\n"
                    protein_sequences.append(fasta_protein)

        # Add protein sequences to data frame
        df['Protein sequence'] = protein_sequences
        return df

    def __add_dna_fasta(df: pd.DataFrame) -> pd.DataFrame:
        '''
        Add DNA sequences to sequences data frame via GenBank ID - sequence
        mapping. Chunking is necessary to avoid weird situations like some
        IDs that should map but don't and then they generate bugs.

        Parameters
        ----------
        df : pd.DataFrame
            Sequences data frame.

        Raises
        ------
        CAZyError
            To stop an infinite loop.

        Returns
        -------
        df : pd.DataFrame
            Sequences data frame + DNA sequences.

        '''
        
        def get_id_from_fasta(description: str) -> str:
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
            match = re.search(r'protein_id=([A-Z0-9_\.]+)\]', description)
            return match.group(1)
        
        # NCBI requests require stating an email
        Entrez.email = ""
        # Input sequences
        genbank = df['GenBank']
        # Initialize output variable
        dna_sequences = []
        # Chunk input
        chunked_genbank = [genbank[i:i + 100] for i in range(0, len(genbank), 100)]
        for chunk in tqdm(chunked_genbank):
            
            chunk = chunk.to_list()
            # Protein sequences
            handle_dna = Entrez.efetch(
                db = "protein", 
                id = chunk, 
                rettype = "fasta_cds_na", 
                retmode = "text"
                )
            SeqIO_dna = SeqIO.parse(handle_dna, 'fasta')
        
            # Make protein FASTA format
            counter = -1

            for dna in SeqIO_dna:
                counter += 1
                # Repeat process with the next GenBank protein ID if there's no match
                while chunk[counter] != get_id_from_fasta(dna.description):
                    counter += 1
                    dna_sequences.append(None)
                    # Security measure to avoid infinite loop
                    if counter == len(chunk):
                        raise CAZy.CAZyError('A forever while loop has been stopped')
                # If there is a match, add FASTA
                else:
                    fasta_dna = f">{dna.description}\n{dna.seq}\n\n"
                    dna_sequences.append(fasta_dna)

        # Add DNA sequences to data frame
        df['DNA sequence'] = dna_sequences
        return df
    
    def __add_uniprot(df: pd.DataFrame) -> pd.DataFrame:
        '''
        Adds UniProt IDs to sequences data frame via GenBank - UniProt mapping

        Parameters
        ----------
        df : pd.DataFrame
            Sequences data frame.

        Returns
        ------
        df : pd.DataFrame
            Sequences data frame + UniProt IDs.

        '''
        
        # Input sequences
        genbank = df['GenBank'].to_list()
        # Uniprot class instantiation
        u = UniProt()
        
        # Some small functions
        def genbank2uniprot(genbank: pd.Series) -> dict:
            '''
            Map GenBank IDs against UniProt to retrieve UniProt IDs

            Parameters
            ----------
            genbank : pd.Series
                GenBank IDs

            Returns
            -------
            dict
                Result of the mapping against UniProt.

            '''
            uniprot_results = u.mapping(
                fr = 'EMBL-GenBank-DDBJ_CDS', 
                to = 'UniProtKB', 
                query = genbank
                )
            return uniprot_results
        
        def divide_process_in_chunks(
                process: Callable[..., Any], 
                input_object: list, 
                chunk_size: int = 500
                ) -> Generator[list[str], None, None]:
            '''
            Divides a process/function into chunks to avoid overloading.

            Parameters
            ----------
            process : Callable[list, Any]
                Process or function that is to be chunked.
            input_object : list
                Process input to be sliced in different chunks.
            chunk_size : int, optional
                Size of the resulting chunks in which the input_object will be 
                divided. The default is 500.

            Yields
            ------
            Generator[list[str]]
                Result of the chunked process.

            '''
            chunked_input = [input_object[i:i + chunk_size] for i in range(0, len(input_object), chunk_size)]
            for chunk in chunked_input:
                yield process(chunk)

        # Get UniProt IDs from the nested dictionaries
        uniprot_ids = defaultdict(list)
        for uniprot_dict in divide_process_in_chunks(genbank2uniprot, genbank):
            # Manage resulted IDs
            for result in uniprot_dict['results']:
                genbank_id = result['from']
                uniprot_id = result['to']['primaryAccession']
                uniprot_ids[genbank_id].append(uniprot_id)

        # Add UniProt IDs to data frame
        has_uniprot = lambda genbank: None if genbank not in uniprot_ids\
            else uniprot_ids[genbank]
        df['UniProt'] = df['GenBank'].apply(has_uniprot)
        
        # Explode column if the length of the lists == 1 (1 to 1 mapping)
        if len(df['UniProt'].dropna()) == len(df['UniProt'].dropna().explode()):
            df = df.explode('UniProt')
        
        return df
    
    def __add_pdb(df: pd.DataFrame) -> pd.DataFrame:
        '''
        Adds PDB IDs to sequences data frame via UniProt ID - PDB ID mapping

        Parameters
        ----------
        df : pd.DataFrame
            Sequences data frame.

        Returns
        ------
        df : pd.DataFrame
            Sequences data frame + PDB IDs.


        '''
        # Uniprot class instantiation
        u = UniProt()
        # Mapping UniProt - PDB
        uniprot_results = u.mapping(
            fr = 'UniProtKB_AC-ID', 
            to = 'PDB', 
            query = df['UniProt'].dropna()
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
        df['PDB'] = df['UniProt'].apply(has_pdb)
    
        return df
    
    def get_full_sequences(self) -> pd.DataFrame:
        '''
        Adds protein sequences, DNA sequences, UniProt IDs and PDB IDs to the
        initial sequence data frame containing information about the domain, 
        family, strain and GenBank ID

        Returns
        -------
        df : pd.DataFrame
            Complete sequence data frame.

        '''
        df = self.get_sequences()
        print('Add protein sequences')
        df = CAZy.__add_protein_fasta(df)
        print('Add DNA sequences')
        df = CAZy.__add_dna_fasta(df)
        print('Add UniProt IDs')
        df = CAZy.__add_uniprot(df)
        print('Add PDB IDs')
        df = CAZy.__add_pdb(df)
        
        return df


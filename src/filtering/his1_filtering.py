'''
Filter and disregard AlphaFold .pdb structures that are not LPMOs. For that, 
the histidine brace will be used as feature for spotting LPMOs. Different
strategies can be employed:
- signalp = SignalP to identify signal peptides.
- domains = custom method in which the presence of His is checked in the 
  transition from signal peptide to protein core. This transition is found 
  leveraging the change of pLDDT values from low (signal peptide) to low 
  (protein core).
- neighbors = custom method in which, for each His in the beginning of the 
  protein, the presence of the other residues of the histidine brace (i.e. His 
  and Trp or Tyr) are checked.

Functions
---------
handle_arguments
starting_residue
manual_check_starts_h
manual_check_not_m_nor_h
his_around_cut
signalp
neighbors
domains
interpro
compare
store
'''

import os
import pickle
import argparse
import itertools
import pandas as pd
from tqdm import tqdm
from .signalP import SignalP
from typing import Tuple
from .parse_pdb import AlphaFoldStructure
from ..config.config_parser import config

def handle_arguments() -> argparse.Namespace:
    '''
    Handles the arguments passed via command line

    Return
    ------
    args : argparse.Namespace
        Arguments
    '''
    
    parser = argparse.ArgumentParser(
        prog = 'his1_filtering',
        description = 'Filters out non-LPMO structures and crops out signal peptide',
        epilog = 'See ya'
        )
    
    # Add arguments
    parser.add_argument(
        '-m', '--method',
        nargs = '+', 
        choices = ['signalp', 'neighbors', 'domains', 'interpro'],
        help = 'Performs any or a combination of different methods to identify signal peptides'
        )
    
    parser.add_argument(
        '-c', '--comparison',
        action = 'store_true',
        help = 'Compares the different methods to identify signal peptides'
        )

    # Arguments from Namespace object
    return parser.parse_args()

def starting_residue(
        structures: list[AlphaFoldStructure]
        ) -> Tuple[
                   list[AlphaFoldStructure], 
                   list[AlphaFoldStructure], 
                   list[AlphaFoldStructure]
                   ]:
    '''
    Separate proteins on whether the starting residue is Met, His or neither.
    Takes 0 secs
    
    Parameters
    ----------
    structures : AlphaFoldStructure
        List of AlphaFoldStructure objects

    Returns
    -------
    starts_m : list[AlphaFoldStructure]
        List of AlphaFoldStructure objects that start with Met
        
    starts_h : list[AlphaFoldStructure]
        List of AlphaFoldStructure objects that start with His
        
    not_m-nor_h : list[AlphaFoldStructure]
        List of AlphaFoldStructure objects that start with not Met nor His
    
    '''
    starts_m, starts_h, not_m_nor_h = [], [], []
    for structure in structures:
        # Get protein seq
        seq = structure.seq
        # Allocate based on starting residue
        if seq.startswith('M'):
            starts_m.append(structure)
        elif seq.startswith('H'):
            starts_h.append(structure)
        else:
            not_m_nor_h.append(structure)
    
    return starts_m, starts_h, not_m_nor_h

def manual_check_starts_h() -> list[AlphaFoldStructure]:
    '''
    The proteins starting with His are manually checked to have LPMO
    characteristics. The ones that pass the check are stored.

    Return 
    ------
    checked_starts_h : list[AlphaFoldStructure]
        Strcutures that pass the check
    '''

    # Manually inspect starts_h:
    starts_h_checked = [
    # 'A0A080ZKX3.pdb', # signal + ? cbm (before ATG)
    'A0A0S2UQK8.pdb', # no signal + 3 cbm 
    'A0A167KNI7.pdb', # no signal + 1 cbm
    'A0A167KNY8.pdb', # no signal + no cbm
    'A0A167KNZ4.pdb', # no signal + no cbm
    'A0A167KNZ9.pdb', # no signal + no cbm
    'A0A167KP22.pdb', # no signal + no cbm
    'A0A1C9ZMC3.pdb', # no signal + no cbm (terrible structure)
    'A0A1C9ZMC5.pdb', # no signal + 1 cbm
    'A0A1C9ZP88.pdb', # no signal + no cbm
    'A0A1C9ZUN9.pdb', # no signal + ? cbm
    'G3XAP7.pdb', # no signal + no cbm
    # 'M5DEQ0.pdb', # wtf is this?????
    # 'M5DEQ1.pdb', # wtf is this?????
    'R9WW03.pdb', # no signal + no cbm
    'V9PCY8.pdb', # no signal + ? cbm
    'V9PD63.pdb' # no signal + ? cbm
    ]

    # Save and store proteins that pass the check
    checked_starts_h = []
    for protein in starts_h_checked:
        structure = AlphaFoldStructure(f"{config['AF_all']}/{protein}")
        checked_starts_h.append(structure)
        indeces_structure = list(range(len(structure)))
        structure.rewrite_full(
            f"{config['AF_his1']}/{structure.id}.pdb", 
            indeces_structure
            )

    return checked_starts_h

def manual_check_not_m_nor_h() -> list[AlphaFoldStructure]:
    '''
    The proteins not starting with His nor Met are manually checked to have 
    LPMO characteristics

    Returns
    -------
    structures : list(AlphaFoldStructure)
        Structures that pass the check
    '''

    not_h_nor_m_checked = [
    # 'A0A059WU06.pdb', # partial protein (None) -> no signal + ? cbm (no homolog at all)
    'A0A081AW22.pdb', # signal + ? cbm
    'A0A0X9QXD6.pdb', # signal + 2 cbm
    'A0A109PR77.pdb', # signal + 2 cbm
    'A0A195BNI8.pdb', # signal + 1 cbm (uniprot says they are transmmebrane)
    'A0A288W730.pdb', # signal + no cbm
    # 'A0A385XGM1.pdb', # wtf is this????? -> after looking at identical proteins, CBM
    # 'A0A385XJ64.pdb', # wtf is this????? -> after looking at identical proteins, CBM
    # 'A0A3G2JNR6.pdb', # partial protein (None) -> no signal + ? cbm
    # 'A0A3G2JQB2.pdb', # partial protein (A0A395P178 193/218 identity) -> no signal + ? cbm
    # 'A0A3G2JRQ0.pdb', # partial protein (G9N0U1 225/226 identity) -> no signal + ? cbm
    # 'A0A5K1JUG4.pdb', # wtf is this????? ->  after looking at identical protein, == A0A5K1JTU3, which is already in df
    # 'A0A7R9BE57.pdb', # wtf is this????? -> after looking at identical proteins, CBM
    # 'A0A7T8H239.pdb', # I guess is a wrong folding -> after BLASTing, still wtf
    # 'A0A7T8K853.pdb', # wtf is this????? / I guess is a wrong folding -> after BLASTing, still wtf
    # 'A0A896INF4.pdb', # wtf is this????? -> after BLASTing == A0A2J7PJZ3, NOT in df
    # 'A7SSA3.pdb', # partial protein (XP_001625522.3)-> no signal + no cbm
    # 'B2AL94.pdb', # partial protein (A0A090CN22) -> no signal + no cbm
    # 'B2B346.pdb', # partial protein (A0A090D9H3) -> no signal + no cbm
    # 'C1KU36.pdb', # partial protein (A0A8S9A495, without structure) -> no signal + no cbm
    # 'C3KGQ8.pdb', # extended protein (M9MS70, alreade in df) -> sigal + no cbm
    # 'E9H5W4.pdb', # partial protein (A0A8J2WML7 187/206 identity but no tiene structure)-> signal + no cbm
    # 'F8QXM9.pdb', # partial protein (A0A423T9X9 97/126 identity) -> no signal + no cbm 
    # 'G2Q9N9.pdb', # partial protein (A0A4Q4Z315 122/152 identity) -> no signal + no cbm
    # 'G2QVT1.pdb', # partial protein -> signal + 1 cbm
    # 'G2RCI4.pdb', # partial protein -> signal + no cbm (but very similar to A0A3S4EX74 and Q2HET7, which are not in the df)
    # 'G4TYC7.pdb', # wtf is this????? -> after BLASTing still wtf
    # 'G5A5X7.pdb', # partial protein -> signal + ? cbm 
    # 'G5A5X9.pdb', # partial protein -> signal + no cbm 
    # 'G5AB86.pdb', # partial protein -> signal + no cbm 
    # 'G5AF67.pdb', # wtf is this????? -> after BLASTing prob CBM
    # 'G5AF69.pdb', # wtf is this????? -> after BLASTing prob CBM
    # 'H6S427.pdb', # partial protein -> signal + ? cbm 
    # 'M5DBI8.pdb', # partial protein (None) -> no signal + no cbm
    # 'M5DCP5.pdb', # partial protein (G2Q9T3, already in df 92/105 identity) -> no signal + no cbm 
    # 'M5DCP8.pdb', # partial protein (A0A2B7WVF0, not df 72/96 identity) -> no signal + no cbm -> this is the proof that there may LPMO not in CAZy, since M5DCP8 and A0A2B7WVF0 both share homoology with 7A8V, which is a characterized LPMO
    # 'M5DCQ2.pdb', # partial protein (C8V530, not in df 86/101 identity) -> no signal + no cbm
    # 'M5DD99.pdb', # partial protein (G2QHR1, in df 73/97 identity) -> no signal + no cbm
    # 'M5DDA9.pdb', # partial protein (Q2GSI4, not in df 81/98 identity) -> no signal + no cbm
    # 'M5DDB1.pdb', # wtf is this????? -> after BLASTing prob partial protein badly folded
    # 'V4BDB3.pdb', # partial protein -> signal + no cbm
    # 'V5N5H9.pdb', # partial protein (Q9RJY2, in df) -> no signal + no cbm
    # 'V9ER26.pdb', # prob not an LPMO, only a ig-like protein since it does not have an initial H and nor have it the 100% indentical proteins
    # 'V9ES95.pdb', # partial protein (A0A329RRC8 235/262 identity) -> no signal + ? CBM
    # 'V9FU93.pdb', # partial protein -> signal + ? CBM
    # 'V9FZ32.pdb', # partial protein (A0A081B0W1 and A0A0W8CVQ3, both not in df) -> no signal + ? CBM
    # 'V9PCZ9.pdb', # partial protein (== G0S5P8, not in df) -> no signal + 1 CBM
    # 'V9PD45.pdb', # partial protein (== G0SDM5, not in df) -> no signal + no CBM
    # 'W2GFG9.pdb', # partial protein -> signal + ? CBM
    # 'W2GFH9.pdb', # partial protein (A0A0W8D1Q2, not in df 253/257 identity) -> no signal + ? CBM
    # 'W2GFQ6.pdb', # partial protein (== W2MZB3, it's in df) -> no signal + ? CBM
    # 'W2GHI9.pdb', # partial protein (A0A0W8D1D3 not in df 204/206 identity) -> no signal + ? CBM
    # 'W2HL20.pdb', # partial protein (A0A081B0X0 and A0A0W8CVV8, both not in df) -> no signal + ? CBM
    # 'W2IDY3.pdb', # partial protein -> signal + ? CBM
    # 'W2ILZ6.pdb', # prob not an LPMO, only a ig-like protein since it does not have an initial H and nor have it the 100% indentical proteins
    # 'W2JTW2.pdb', # same result as W2HL20
    # 'W2LDN6.pdb', # partial protein (A0A0W8DIU2, not in df) -> no signal + no CBM
    # 'W2MYY2.pdb', # partial protein -> signal + ? CBM
    # 'W2N176.pdb', # same result as W2GFH9
    # 'W2VSR0.pdb', # partial protein -> signal + ? CBM
    # 'W2XP19.pdb', # partial protein -> signal + ? CBM
    # 'W2YX14.pdb', # partial protein -> signal + ? CBM
    # 'W4FXF0.pdb'  # partial protein -> signal + 1 CBM
    ]

    # Get structures that manually pass the check
    structures = []
    for protein in not_h_nor_m_checked:
        structure = AlphaFoldStructure(f"{config['AF_all']}/{protein}")
        structures.append(structure)
    
    return structures

def hist_around_cut(
    seq: str, 
    cut: int, 
    window_length: int
    ) -> int:
    '''
    The AlphaFoldStructure.domains() method is based purely on pLDDT scores,
    and it is often the case that there is a His but a 3 or less residues 
    apart from when it had been indicated by the method. It is also useful 
    to adjust the slightly off results of SignalP or verify the results from
    InterPro. This function checks if there is a His in the vicinity (window 
    length) of a given position (cut) of a sequence (seq) and returns its 
    index.

    Parameters
    ----------
    seq : str
        Protein sequence.
    cut : int
        Point of the sequence that represents the cut of a peptidase.
    window_length : int
        Length of the window around the cut point to look for His.

    Returns
    -------
    0 if His is not found, its index otherwise.

    '''
    # Adjust window in margins
    if cut - window_length < 0:
        window = seq[0 : cut + window_length]
    elif cut + window_length > len(seq):
        window = seq[cut - window_length : len(seq)]
    else:
        window = seq[cut - window_length : cut + window_length]

    # Explore vicinities for His1
    if seq[cut] == 'H': 
        return cut
    elif (h_found := window.find('H')) == -1:
        return None
    else:
        return cut + h_found - window_length
   
def signalp(structures : list[AlphaFoldStructure]) -> dict[str, int]:
    '''
    The resulting data frame from having run SignalP 6.0 on a list of protein 
    structure's sequences. The resulting gets parsed to know whether a sequence 
    has been predicted to have a signal peptide, and if so, where. Furthermore, 
    the positive results are validated. A UniProt entry will be ruled out if a 
    His is not identified in a window of 3 residues upstream and downstream of 
    the predicted cut point. This is done by the function "his_around_cut". In 
    summary, only when a sequence is detected by SignalP to have a signal 
    peptide and a His is located in the vicinities of the predicted cut, the 
    structure will pass the filter.

    Parameters
    ----------
    structures : list[AlphaFoldStructure]
        List of AlphaFold structures

    Returns
    ----------
    filter_result : dict[str, int]
        Dictionary containing the UniProt IDs paired with the sequence index of
        the identified His. When SignalP has predicted 'OTHER' or when 
        predicting 'SP' no His has been found around the cut point, the UniProt
        ID is paired with the value 0.
    '''

    # Parse SignalP results
    df = SignalP(f"{config['SignalP']}/prediction_results.txt").df

    # Instantiate return variable
    other = df[df['Prediction'] == 'OTHER'].index.to_list()
    filter_result = {uniprot:None for uniprot in other}

    # Iterate over non-other proteins
    cuts = df[df['Prediction'] != 'OTHER']['CUT']
    for uniprot, cut in tqdm(cuts.items(), total = len(cuts)):
        # Get sequence
        structure = AlphaFoldStructure(f"{config['AF_all']}/{uniprot}.pdb")
        seq = structure.seq

        # Search for His int the vicinities around cut
        filter_result[uniprot] = hist_around_cut(
            seq = seq, 
            cut = cut, 
            window_length = 3
            ) 
    
    return filter_result

def neighbors(structures : list[AlphaFoldStructure]) -> dict[str, int]:
    '''
    For each His found in the N-terminal region of the protein, the neighboring 
    residues are identified. If another His and Trp/Tyr are in the neighboring
    residues, the histidine brace (i.e. active site) has been identified.

    Parameters
    ----------
    structures : list[AlphaFoldStructure]
        List of AlphaFold structures
    
    Returns
    -------
    filter_result : dict[str, int]
        Dictionary containing the UniProt IDs paired with the sequence index of
        the identified His belonging to a histidine brace. When no His 
        belonging to a histidine brace is found, the UniProt ID is paired with 
        the value None.
    '''
    # Instantiate return variable
    filter_result = {}

    # Iterate over each residue of each structure
    for structure in tqdm(structures):
        for index, residue in enumerate(structure.seq[:70]):
            # Get neighboring residues
            if residue == 'H':
                neighbors = structure.neighbours(
                    center = ('H', index), 
                    radius = 11
                ).index.to_list()
                
                # Identify His in neghboring residues´
                neighbors = set(neighbors)
                if ('H' in neighbors and 'Y' in neighbors) or\
                   ('H' in neighbors and 'F' in neighbors):
                    filter_result[structure.id] = index
                    break
                elif structure.pLDDT[index] > 90:
                    filter_result[structure.id] = None
                    break
                else:
                    filter_result[structure.id] = None

    return filter_result

def domains(structures : list[AlphaFoldStructure]) -> dict[str, int]:
    '''
    Use AlphaFoldStructure.domains() method to identify the point where the
    signal peptide (normally of low quality) meets the core (normally of high)
    quality. Those structures who have His around that point pass the filter.
    
    Parameters
    ----------
    structures : list[AlphaFoldStructure]
        List of AlphaFoldStructure objects.

    Returns
    -------
    filter_result : dict[str, int]
        Dictionary containing the UniProt IDs paired with the sequence index of
        the identified His belonging to a histidine brace. When no His 
        belonging to a histidine brace is found, the UniProt ID is paired with 
        the value None.
    '''

    # Initiate output variables
    filter_result = {}
    
    for structure in tqdm(structures):
        # Get domains
        domain = structure.domains(threshold = 75)
        # Account for failed structures
        if domain == []:
            filter_result[structure.id] = None
        # Account for only core structures
        elif domain[0][0] == 0:
            filter_result[structure.id] = None
        # Account for structures with identified domains
        else:
            seq = structure.seq
            cut = domain[0][0]
            # Account for structures with His immediatelly after the cut
            if seq[cut] == 'H':
                filter_result[structure.id] = cut
            # Account for structures with His but around the cut
            elif (h_loc := hist_around_cut(seq, cut, 3)):
                filter_result[structure.id] = h_loc
            # Account for wtf structures
            else:
                filter_result[structure.id] = None
    
    return filter_result

def interpro(structures : list[AlphaFoldStructure]) -> dict[str, int]:
    '''
    InterPro entries include several proteins that share the same domain or 
    belong to the same familie thanks to a common domain-like region. Each
    protein has a description of where this domain/region is. This information
    has been previously retrieved by the InterPro class and this function finds
    whether the proteins are in the InterPro entries associated with LPMOs. If
    so, get the domain range and validates the presence of a His residue at the
    start.

    Parameters
    ----------
    structures : list[AlphaFoldStructure]
        List of AlphaFold structures
    
    Returns
    -------
    filter_result : dict[str, int]
        Dictionary containing the UniProt IDs paired with the sequence index of
        the identified His belonging to a histidine brace. When no His 
        belonging to a histidine brace is found, the UniProt ID is paired with 
        the value None.
    '''
    
    # Initiate output variables
    filter_result = {}

    # Open pickled dict
    with open(f'{config["InterPro"]}/interpro.pkl', 'rb') as handle:
        domains = pickle.load(handle)

    for structure in tqdm(structures):
        # Get the location of His1 predicted by InterPro
        try:
            start = domains[structure.id][0][0]
        # Except the UniProt ID is not in InterPro
        except IndexError:
            filter_result[structure.id] = None

        # Validate start  
        else:
            filter_result[structure.id] = hist_around_cut(
                seq = structure.seq, 
                cut = start - 1, 
                window_length = 3
                ) 

    return filter_result

def compare(results : dict[str, int]) -> pd.DataFrame:
    '''
    The UniProt entries where the signal peptide can be successfully (i.e. a
    His1 residue has been found) removed are stored in the results dictionary
    with the index of the His1 wheread those who cannot are stored with a None
    value.
    This function counts the number of not None per method and per combination
    of method and returns them in data frame format. The columns are all the 
    method combination and the two columns are the number of UniProt IDs that 
    pass and don't pass the filter. It's ordered from more to less counts.

    Paramaters
    ----------
    results : dict[str, int]
        Dictionary containing the UniProt IDs paired with the sequence index of
        the identified His belonging to a histidine brace. When no His 
        belonging to a histidine brace is found, the UniProt ID is paired with 
        the value None.

    Returns
    -------
    df : pd.DataFrame
        The columns are all the method combination and the two columns are
        'PASS' and 'NO PASS', so each value represents the number of UniProt 
        IDs that pass and don't pass the filter. It's ordered from more to less 
        counts.
    '''

    # Initialize output variable
    df = {}

    # Calculate number of entries and methods
    n_entries = len(list(results.values())[0])
    n_methods = len(results)

    # Get all possible method combination 
    for n in range(n_methods):
        for combination in itertools.combinations(results, n + 1): 
            
            # For each method in the combination, get the entries that pass the filter
            set_pass = set()
            for method in combination:
                result = results[method]
                yes_pass = set([uniprot for uniprot, his1 in result.items()\
                                if his1 is not None])
                set_pass = set_pass | yes_pass

            # Get total number of entries that pass and don't the filter
            name_combination = ' + '.join(combination)
            length_pass = len(set_pass)
            df[name_combination] = (length_pass, n_entries - length_pass)

    # Construct output data frame
    df = pd.DataFrame.from_dict(
        data = df, 
        orient = 'index', 
        columns = ('PASS', 'NO PASS')
        )
    df = df.sort_values(by = 'PASS', ascending = False)

    return df
            
def store(results : dict[str, int]):
    '''
    Stores the structure .pdb files of the UniProt IDs that pass the filter 
    specified by the combination of the filtering methods chosen.

    Paramaters
    ----------
    results : dict[str, int]
        Dictionary containing the UniProt IDs paired with the sequence index of
        the identified His belonging to a histidine brace. When no His 
        belonging to a histidine brace is found, the UniProt ID is paired with 
        the value None.

    Returns
    -------
    None
    '''

    # UniProt : his1_index
    uniprot_his1 = {}

    # Iterate over the methods' results
    for result in results.values():
        # Get the UniProt : His1 index pair of thos UniProt that pass the filter
        yes_pass = {uniprot : his1 for uniprot, his1 in result.index() if his1 is not None}
        # Update the output dictionary
        uniprot_his1 = uniprot_his1 | yes_pass

    # Store structures
    for structure, his1_index in uniprot_his1.index():
        structure.rewrite_range(
            f"{config['AF_his1']}/{structure.id}.pdb", 
            (his1_index, len(structure))
            )
        
def main():
    '''Program flow.'''
    
    # Arguments
    args = handle_arguments()

    # Proteins to evaluate    
    structures = [AlphaFoldStructure(f"{config['AF_all']}/{pdb_file}")\
                   for pdb_file in tqdm(os.listdir(config['AF_all']))]
    
    # Check starting residue
    starts_m, starts_h, not_h_nor_m = starting_residue(structures)
    manual_check_starts_h()
    valid = manual_check_not_m_nor_h()
    structures = starts_m + valid
    
    # Available methods
    methods = {
        'signalp' : signalp,
        'neighbors' : neighbors,
        'domains' : domains,
        'interpro' : interpro
    }
    
    # Run the passed methods
    method_results = {}
    for method in args.method:
        method_results[method] = methods[method](structures)
    
    # Perform comparison or store the results?
    if args.comparison:
        # Compare methods
        comparison = compare(method_results)
        print(comparison) 
    else:
        # Store results
        store(method_results)
     
if __name__ == '__main__':
    main()
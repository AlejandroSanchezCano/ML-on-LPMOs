'''
Filter and disregard AlphaFold .pdb structures that:
    - Have an overall bad quality structure.
    - Do not start with His after removing the signal peptide based on 
      N-terminal pLDDT, and therefore not considered as LPMOs.

Functions
---------
download_af_files
bad_quality
starting_residue
explore_starts_m
'''

import os
import argparse
import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import Union, Tuple
from Bio import BiopythonParserWarning
from parse_pdb import AlphaFoldStructure
from ..config.config_parser import config
warnings.simplefilter('ignore', BiopythonParserWarning)

def handle_arguments():
    '''
    Handles the arguments passed via command line

    Return
    ------
    args : argparse.Namespace
        Arguments
    '''

    # Parser
    parser = argparse.ArgumentParser(
        prog = 'filtering',
        description = 'Filters out non-LPMO structures or with low quality',
        epilog = 'See ya'
        )
    
    # Add arguments
    parser.add_argument(
        '-l', '--length', 
        type = int, 
        default = 20,
        help = 'Peptide length resulting from in-silico fragmentation'
        )
    
    parser.add_argument(
        '-c', '--cutoff',
        tyoe
    )

    # Get argument values
    args = parser.parse_args()
    return args

def bad_quality(
        peptide_length: int,
        cutoff: Union[int, float] 
        ) -> Tuple[list[str], list[str]]:    
    '''
    Impose a filtering based on overall poor quality structure. To do so, the
    protein will be diveded into peptides of length peptide_length. Then,
    the average pLDDT per peptide is calculated. The protein that do not have
    a peptide with an average pLDDT above a cutoff will not subject of future
    analysis.
    Takes 6 min for 4000 proteins.

    Parameters
    ----------
    peptide_length : int, optional
        Length of the peptides that the protein is divided into. The default is
        20
        
    cutoff : Union[int, float], optional
        pLDDT cutoff used to assess peptide quality. The default is 90.
    
    Returns
    -------
    kept : list[AlphaFoldStructure]
        List of AlphaFoldStructure objects that pass the quality filtering.
    
    removed : list[AlphaFoldStructure]
        List of AlphaFoldStructure objects that NOT pass the quality filtering.

    '''

    kept, removed = [], []
    for structure in tqdm(os.listdir(config['AF_all'])):
        # Get per-residue pLDDT
        af_structure = AlphaFoldStructure(f"{config['AF_all']}/{structure}")
        quality = af_structure.pLDDT
        # Divide protein in peptides -> matrix of pLDDT values
        peptides = [
            quality[position : position + peptide_length] \
            for position in range(0, len(quality), peptide_length)
            ]
        # Calculate mean pLDDT per peptide
        mean_pLDDT_per_peptide = np.array([peptide.mean()\
                                           for peptide in peptides])
        
        # Assess quality > cutoff
        if any(mean_pLDDT_per_peptide >= cutoff):
            kept.append(af_structure)
        else:
            removed.append(af_structure)
    
    return kept, removed

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
    for structure in tqdm(structures):
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
        structure = AlphaFoldStructure(f'{AF_FILES}/{protein}')
        checked_starts_h.append(structure)
        indeces_structure = list(range(len(structure)))
        structure.rewrite_full(
            f'{AF_HIS1}/{structure.id}.pdb', 
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
        structure = AlphaFoldStructure(f'{AF_FILES}/{protein}')
        structures.append(structure)
    
    return structures

def explore_starts_m(
        structures: list[AlphaFoldStructure]
        ) -> Tuple[list[AlphaFoldStructure], dict]:
    '''
    Use AlphaFoldStructure.domains() method to identify the point where the
    signal peptide (normally of low quality) meets the core (normally of high)
    quality. Those structures who have His around that point pass the filter
    and are stored in variables.ALPHAFOLD_STARTS_WITH_HIS file
    Takes 45 secs wirh 3866 structures.
    
    Parameters
    ----------
    structures : list[AlphaFoldStructure]
        List of AlphaFoldStructure objects.

    Returns
    -------
    success : [list[AlphaFoldStructure]
        List of AlphaFoldStructure objects that have His after removing the
        signal peptide
    
    fails : dict
    '''
    
    def hist_around_cut(
            seq: str, 
            cut: int, 
            window_length: int
            ) -> Union[bool, int]:
        '''
        The AlphaFoldStructure.domains() method is based purely on pLDDT scores,
        and it is often the case that there is a His but a 3 or less residues 
        apart from when it had been indicated by the method. This function 
        checks if there is a His in the vicinity (window length) of a given 
        position (cut) of a sequence (seq) and returns its index.

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
        False if His is not found, its index otherwise.

        '''
        window = seq[cut - window_length : cut + window_length]
        h_found = window.find('H')
        if h_found == -1:
            return False
        else:
            return cut + h_found - window_length
    
    # Add 1 protein from starts_h
    added = [AlphaFoldStructure(f'{AF_FILES}/A0A080ZKX3.pdb')]

    # Initiate output variables
    success = []
    fails = {
        'no_domains_identified' : [],
        'only_core': [],
        'no_idea': []
        }
    
    for structure in tqdm(structures + added):
        # Get domains
        domain = structure.domains(threshold = 75)
        # Account for failed structures
        if domain == []:
            fails['no_domains_identified'].append(structure)
        # Account for only core structures
        elif domain[0][0] == 0:
            fails['only_core'].append(structure)
        # Account for structures with identified domains
        else:
            seq = structure.seq
            cut = domain[0][0]
            # Account for structures with His immediatelly after the cut
            if seq[cut] == 'H':
                success.append(structure)
                structure.rewrite_range(
                    f'{AF_HIS1}/{structure.id}.pdb', 
                    [(cut, len(structure))]
                    )
            # Account for structures with His but around the cut
            elif (h_loc := hist_around_cut(seq, cut, 3)):
                success.append(structure)
                structure.rewrite_range(
                    f'{AF_HIS1}/{structure.id}.pdb', 
                    [(h_loc, len(structure))]
                    )
            # Account for wtf structures
            else:
                fails['no_idea'].append(structure)
    
    return success, fails

def custom_vs_signalp():
    pass

def filter_supreme_df(structures: list[AlphaFoldStructure]) -> None:

    # Import supreme data frame
    supreme_df = pd.read_pickle(f'{CAZY_EXPANDED}/AA_supreme')

    # Get AlphaFold IDs
    filtered_proteins = [structure.id for structure in structures]
    filtered_proteins = pd.DataFrame({'UniProt' : filtered_proteins})

    # Inner join
    filtered = pd.merge(
        left = supreme_df,
        right = filtered_proteins,
        how = 'inner',
        on = 'UniProt',
    )

    # Store filtered data frame
    filtered.to_pickle(f'{CAZY_EXPANDED}/AA_filtered')

        
def main():
    '''Program flow.'''

    # Argument parsing
    args = handle_arguments()

    # Filtering by overall poor structural quality
    print('Filtering bad quality proteins')
    kept, removed = bad_quality(peptide_length = 20, cutoff = 90)
    
    # Separate proteins on whether the starting residue is Met, His or neither
    print(f'Analyzing starting residue of {len(kept)} proteins')
    starts_m, starts_h, not_h_nor_m = starting_residue(kept) 
    
    # Manually check starts_h proteins
    print(f'Check {len(starts_h)} proteins that start with His')
    checked_starts_h = manual_check_starts_h()

    # Manually check not_m_nor_h proteins
    print(f'Check {len(not_h_nor_m)} proteins that do not start with His nor Met')
    checked_not_m_nor_h = manual_check_not_m_nor_h()

    # Remove signal peptide
    print(f'Remove signal peptide of {len(starts_m + checked_not_m_nor_h)} proteins')
    success, fails = explore_starts_m(starts_m + checked_not_m_nor_h)

    # Compare with SignalP
    print("Comparing custom method with SingalP's")
    both, none, custom, signalp = custom_vs_signalp(success, fails)

    # Make filtered data frame
    filtered_df = filter_supreme_df(success + checked_starts_h)

if __name__ == '__main__':
    main()
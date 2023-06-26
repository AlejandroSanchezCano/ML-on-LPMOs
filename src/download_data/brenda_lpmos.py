'''
BRENDA parser

Parses BRENDA information related to LPMOs, compare it to the existent LPMOs in
the information retrieved from CAZy and add the missing parts.

Functions
---------
LPMOs_from_BRENDA
LPMOs_from_CAZy
BRENDA_vs_CAZy
modify_df

'''

import requests
import pandas as pd
from tqdm import tqdm
from bs4 import BeautifulSoup
from bioservices import UniProt
from variables import CAZY_EXPANDED


def LPMOs_from_BRENDA() -> dict:
    '''Parse BRENDA HTML tables corresponding to LPMO enzymatic activities.
    
    Returns
    -------
    brenda dic : dict
        Dictionary with LPMO EC as keys and a list with protein fasta text as
        values

    '''
    
    u = UniProt()
    
    brenda_dic = {
        '1.14.99.53' : [],
        '1.14.99.54' : [],
        '1.14.99.55' : [],
        '1.14.99.56' : []
        }
    
    # Loop over LPMO-corresponding EC numbers
    for ec in tqdm(brenda_dic):
        link = f'https://www.brenda-enzymes.org/all_enzymes.php?ecno={ec}&table=Sequence#TAB'
        soup = BeautifulSoup(requests.get(link).text, features = 'lxml')
        table = soup.find("div", attrs={"id": "tab32"})
        rows = table.findAll('div', attrs = {'class' : ["row rgrey1", 'row rgrey2']})
        
        # Loop over the elements in each EC table
        for row in rows:
            for cell in row.find('div', attrs = {'class' :'cell'}):
                # Get UniProt ID from table
                uniprot_id = cell.get_text(strip = True).replace('pBLAST', '')
                # Get sequence from UniProt ID
                seq = u.get_fasta(uniprot_id)
                # Get UniProt ID from sequence
                get_uniprot_from_fasta = lambda fasta: fasta.split('|')[1]
                # Append
                brenda_dic[ec].append(get_uniprot_from_fasta(seq))
        
    return brenda_dic

def LPMOs_from_CAZY() -> pd.DataFrame:
    '''
    Merges all the sequence data frames into a single one containing info from
    all the families

    Returns
    -------
    df_supreme : pd.DataFrame
        Data frame containing info from all the families.
    '''
    
    # Dict with all the data frame objects
    dfs = {
        'AA0' : None,
        'AA9' : None,
        'AA10' : None,
        'AA11' : None,
        'AA13' : None,
        'AA14' : None, 
        'AA15' : None,
        'AA16' : None,
        'AA17' : None       
        }
    
    # Record individual data frames
    for df_name in dfs:
        df = pd.read_pickle(f'{CAZY_EXPANDED}/{df_name}')
        dfs[df_name] = df

    # Create supreme data frame
    df_supreme = pd.concat(list(dfs.values()), ignore_index=True)
    
    return df_supreme

def BRENDA_vs_CAZy(brenda_lpmos: dict, cazy_lpmos: pd.DataFrame) -> set:
    '''
    Return UniProt IDs present in BRENDA database but not in CAZy database.

    Parameters
    ----------
    brenda_lpmos : dict
        LPMOs in BRENDA database.
    cazy_lpmos : pd.DataFrame
        LPMOs in CAZy database.

    Returns
    -------
    set
        UniProt IDs present in BRENDA database but not in CAZy database..

    '''
    
    brenda_uniprot = [lpmo for ec in brenda_lpmos.values() for lpmo in ec]
    cazy_uniprot = cazy_lpmos['UniProt'].dropna()
    
    return set(brenda_uniprot) - set(cazy_uniprot)

def modify_df(df: pd.DataFrame) -> pd.DataFrame:
    '''
    Add entries from BRENDA missing in the supreme data frame and modify the
    entries where information (provided from BRENDA) is lacking

    Parameters
    ----------
    df : pd.DataFrame
        Data frame with info of all LPMOs.

    Returns
    -------
    df : pd.DataFrame
        Data frame with info of all LPMOs + BRENDA.

    '''
    
    # New entries
    # {'A0A8J8W0G0', # no structure
    #  'A0A8J8W6X8', # no structure
    #  'D0VWZ9', # PDB version of existing entry
    #  'G0R6T8', # 100% identical proteins already in df
    #  'G0RVK1', # not LMPMO -> BRENDA says it's LPMO but UniProt says it's GH
    #  'G3XAP7', # NEW but has not GenBank ID
    #  'Q5B027', # Already in db with another ID (A0A1U8QPP4)
    #  'Q5B1W7', # Already in db with another ID (A0A1U8QN05)
    #  'Q7RWN7', # Updated version of an existing ID
    #  'Q7SA19'} # Updated version of an existing ID 
    
    # 'D0VWZ9' (not in df) is 'G2RGE5' (in df) without signal peptide 
    loc = df[df['UniProt'] == 'G2RGE5'].index.item()
    df.at[loc, 'PDB'] = ['3EII', '3EJA']
    # 'G0R6T8' (not in df) == 'A0A4D6DD21', 'H2KXF9', 'O14405' (in df)
    loc = df[df['UniProt'] == 'A0A4D6DD21'].index.item()
    df.at[loc, 'PDB'] =['5O2W', '5O2X']
    loc = df[df['UniProt'] == 'H2KXF9'].index.item()
    df.at[loc, 'PDB'] =['5O2W', '5O2X']
    loc = df[df['UniProt'] == 'O14405'].index.item()
    df.at[loc, 'PDB'] =['5O2W', '5O2X']
    # It turns out that 'Q7RWN7' is the updated GenBank version of an existing ID
    loc = df[df['GenBank'] == 'EAA26873.1'].index.item()
    df.at[loc, 'GenBank'] = 'EAA26873.2'
    df.at[loc, 'Protein sequence'] = '>EAA26873.2 endoglucanase II [Neurospora crassa OR74A]\nMRSTLVTGLIAGLLSQQAAAHATFQALWVDGADYGSQCARVPPSNSPVTDVTSNAMRCNTGTSPVAKKCPVKAGSTVTVEMHQQANDRSCSSEAIGGAHYGPVLVYMSKVSDAASADGSSGWFKIFEDTWAKKPSSSSGDDDFWGVKDLNSCCGKMQVKIPSDIPAGDYLLRAEVIALHTAASAGGAQLYMTCYQISVTGGGSATPATVSFPGAYKSSDPGILVDIHSAMSTYVAPGPAVYSGGSSKKAGSGCVGCESTCKVGSGPTGTASAVPVASTSAAAGGGGGGGSGGCSVAKYQQCGGTGYTGCTSCASGSTCSAVSPPYYSQCV\n\n'
    df.at[loc, 'DNA sequence'] = '>lcl|CM002238.1_cds_EAA26873.2_1 [gene=gh61-5] [locus_tag=NCU08760] [protein=endoglucanase II] [protein_id=EAA26873.2] [location=join(5242656..5242904,5243111..5243800,5243863..5243916)] [gbkey=CDS]\nATGCGGTCCACTCTTGTCACCGGCCTCATCGCCGGCCTACTCTCCCAACAAGCCGCCGCCCACGCCACCTTCCAAGCCCTTTGGGTCGATGGTGCCGATTATGGCTCGCAATGCGCTCGCGTCCCTCCTTCCAACTCCCCCGTCACCGATGTGACTAGCAATGCCATGAGGTGTAACACGGGAACTTCGCCCGTTGCGAAGAAGTGCCCTGTCAAGGCGGGAAGTACGGTCACTGTTGAGATGCACCAGCAAGCAAATGACCGCTCCTGTTCCTCTGAAGCCATCGGTGGCGCTCACTACGGTCCCGTCCTCGTGTATATGTCCAAGGTCTCCGACGCCGCCTCCGCCGACGGTTCCTCTGGCTGGTTCAAGATCTTTGAGGACACCTGGGCCAAGAAGCCCTCCAGCTCCTCGGGCGACGATGATTTCTGGGGCGTCAAAGACCTCAACTCGTGCTGCGGCAAGATGCAGGTCAAGATCCCCTCGGACATCCCCGCGGGTGACTATCTCCTCCGTGCCGAGGTTATCGCGCTCCATACCGCCGCAAGCGCGGGAGGTGCCCAGTTGTACATGACCTGCTACCAGATCTCCGTTACCGGTGGTGGCTCCGCTACCCCGGCGACTGTCAGCTTTCCTGGTGCCTACAAGAGCTCCGACCCTGGTATCCTCGTTGACATCCACAGTGCCATGAGCACCTACGTCGCCCCCGGACCGGCTGTGTACTCGGGTGGAAGCTCCAAGAAGGCCGGAAGCGGCTGCGTGGGCTGCGAGTCTACTTGCAAGGTTGGCTCCGGCCCGACTGGAACTGCTTCTGCCGTCCCTGTTGCGAGCACGTCGGCGGCTGCTGGTGGTGGAGGCGGTGGTGGGAGCGGTGGCTGCAGCGTTGCAAAGTATCAGCAGTGTGGTGGAACCGGCTATACCGGGTGCACATCCTGCGCTTCCGGATCCACCTGCAGCGCTGTCTCACCTCCTTATTACTCCCAGTGTGTCTAA\n\n'
    df.at[loc, 'UniProt'] = 'Q7RWN7'
    # It turns out that 'Q7SA19' is the updated GenBank version of an existing ID
    loc = df[df['GenBank'] == 'EAA33178.1'].index.item()
    df.at[loc, 'GenBank'] = 'EAA33178.2'
    df.at[loc, 'Protein sequence'] = '>EAA33178.2 endoglucanase IV [Neurospora crassa OR74A]\nMKTFATLLASIGLVAAHGFVDNATIGGQFYQPYQDPYMGSPPDRISRKIPGNGPVEDVTSLAIQCNADSAPAKLHASAAAGSTVTLRWTIWPDSHVGPVITYMARCPDTGCQDWTPSASDKVWFKIKEGGREGTSNVWAATPLMTAPANYEYAIPSCLKPGYYLVRHEIIALHSAYSYPGAQFYPGCHQLQVTGSGTKTPSSGLVSFPGAYKSTDPGVTYDAYQAATYTIPGPAVFTC\n\n'
    df.at[loc, 'DNA sequence'] = '>lcl|CM002239.1_cds_EAA33178.2_1 [locus_tag=NCU07898] [protein=endoglucanase IV] [protein_id=EAA33178.2] [location=complement(join(1412641..1412684,1412787..1413039,1413178..1413477,1413564..1413590,1413729..1413821))] [gbkey=CDS]\nATGAAGACCTTTGCGACTCTTTTGGCTTCCATCGGCCTGGTGGCCGCTCACGGCTTTGTTGATAACGCCACTATTGGTGGTCAGTTTTATCAACCGTACCAGGACCCCTACATGGGCAGCCCCCCCGATCGAATCTCTCGTAAGATTCCCGGCAACGGCCCCGTCGAAGACGTCACTTCCCTCGCCATTCAGTGCAACGCCGACTCAGCCCCGGCCAAGCTTCATGCGTCCGCCGCCGCCGGATCGACTGTCACTTTGCGCTGGACCATTTGGCCCGACTCGCACGTGGGACCCGTCATCACCTACATGGCCCGCTGTCCCGACACGGGGTGCCAGGACTGGACCCCTAGCGCCAGTGATAAGGTGTGGTTCAAGATTAAGGAAGGTGGGAGGGAGGGAACGAGTAATGTTTGGGCTGCTACCCCCCTCATGACCGCCCCGGCCAACTACGAGTACGCCATCCCGTCCTGCCTCAAGCCCGGTTACTATCTGGTTAGGCACGAGATCATTGCGCTGCACAGCGCCTACTCTTATCCTGGTGCTCAGTTCTACCCGGGATGCCATCAGTTGCAGGTGACAGGTTCGGGAACCAAGACGCCCAGCTCGGGACTGGTCAGTTTCCCGGGCGCGTACAAGAGTACTGATCCGGGGGTTACTTATGATGCTTACCAGGCTGCCACTTATACCATCCCCGGTCCTGCTGTGTTTACTTGCTAA\n\n'
    df.at[loc, 'UniProt'] = 'Q7SA19'
    df.at[loc, 'PDB'] = '4EIS'
    # Add 'G3XAP7' row
    loc = len(df)
    df.loc[loc] = [
        'AA9', 
        'Eukaryota',
        'Thermoascus aurantiacus',
        None,
        '>tr|G3XAP7|G3XAP7_THEAU Endo-beta-1,4-glucanase D OS=Thermoascus aurantiacus OX=5087 PE=1 SV=1\nHGFVQNIVIDGKNYGGYLVNQYPYMSNPPEVIAWSTTATDLGFVDGTGYQTPDIICHRGAKPGALTAPVSPGGTVELQWTPWPDSHHGPVINYLAPCNGDCSTVDKTQLEFFKIAESGLINDDNPPGIWASDNLIAANNSWTVTIPTTIAPGNYVLRHEIIALHSAQNQDGAQNYPQCINLQVTGGGSDNPAGTLGTALYHDTDPGILINIYQKLSSYIIPGPPLYTG\n\n',
        None,
        'G3XAP7',
        ['2YET', '3ZUD', '7PU1', '7PZ3', '7PZ4', '7PZ5', '7PZ6', '7PZ7', '7PZ8', '7Q1K']
        ]
    
    return df

def main():
    '''Progress flow.'''
    
    # Parse LPMO BRENDA
    brenda_lpmos = LPMOs_from_BRENDA()
    
    # Import CAZy data frame
    cazy_lpmos = LPMOs_from_CAZY()
    
    # Compare BRENDA df vs CAZy df
    new_entries = BRENDA_vs_CAZy(brenda_lpmos, cazy_lpmos)

    # Modify df
    CAZy_plus_BRENDA = modify_df(cazy_lpmos)
    
    # Store
    CAZy_plus_BRENDA.to_pickle(f'{CAZY_EXPANDED}/AA_supreme')

if __name__ == '__main__':
    main()
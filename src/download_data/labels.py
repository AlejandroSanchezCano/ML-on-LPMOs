'''
Labels database.

CAZy and the review Frommhagen, M. et. al, 2019 were both used to accrue 
information about the LPMOs that have been characterized to be used as labels
in supervised ML experiments. Some descriptive statistics are also calculated 
and plotted to get a clearer idea of the frequency of characterized LPMOs and
features (regioselectivity and substrate specificity) per family.

Functions
---------
make_database
characterized_vs_all
add_new_entries
plot_regioselectivity
plot_substrate_specificity
plot_per_family
main
'''

import os
import pandas as pd
from tqdm import tqdm
from typing import Union
from .manual_entries import ManualEntries
from ..config.config_parser import config
from ..filtering.parse_pdb import AlphaFoldStructure


def make_database() -> pd.DataFrame:
    '''
    CAZy and the review Frommhagen, M. et. al, 2019 were both used to accrue 
    information about the LPMOs that have been characterized to be used as 
    labels in supervised ML experiments.

    Returns
    -------
    db : pd.DataFrame
        Characterized LPMOs
    '''

    # Manually populate the database
    db = [
        # AA9
        ('AA9', 'AfAA9B', 'Q4WP32', 'C1/C4', ('cellulose'), (21, 250), False, ('5x6a', '6h1z', '6ha5', '6haq')),
        ('AA9', 'AN3046', 'A0A1U8QMG7', 'C1', ('cellulose', 'xyloglucan'), (17, 238), False, None),
        ('AA9', 'AN1602.2', 'Q5BCX8', 'C4', ('cellulose'), (18, 231), True, None),
        ('AA9', 'CsLPMO9A', 'A0A6C0PI79', 'C1', ('cellulose'), (20, 245), True, None),
        ('AA9', 'CsLPMO9', 'M2RAI8', 'C1', ('cellulose'), (18, 235), False, ('7exk')),
        ('AA9', 'CsLPMO9B', 'M2QEW2', 'C1', ('cellulose'), (17, 243), True, None),
        ('AA9', 'CvAA9A', 'A0A223GEC9', 'C1/C4', ('cellulose', 'xyloglucan', 'glucomannan'), (22, 248), False, ('5nlt', '6ydc', '6ydd', '6yde', '6ydf')), 
        ('AA9', 'FgLPMO9A', 'I1REU9', 'C1/C4', ('cellulose', 'xyloglucan'), (21, 241), False, None),
        ('AA9', 'GtLPMO9A-2', 'A0A1C9ZMC5', 'C1/C4', ('cellulose', 'xyloglucan', 'glucomannan'), (0, 230), True, None),
        ('AA9', 'GtLPMO9B', 'F8T947', 'C1/C4', ('cellulose'), (20, 252), False, None),
        # ('AA9', 'GcLPMO9A', 'A0A0J9XKT0', 'C1/C4', ('cellulose', 'xyloglucan')), # no structure
        ('AA9', 'GcLPMO9B', 'A0A0J9XL55', 'C1/C4', ('cellulose', 'xyloglucan'), (18, 240), True, None), 
        ('AA9', 'HiLPMO9B', 'W4KMP1', 'C1', ('cellulose'), (19, 244), False, ('5nns')),
        ('AA9', 'HiLPMO9H', 'W4K498', 'C1', ('cellulose'), (17, 245), True, None),
        ('AA9', 'HiLPMO9I', 'W4K8M0', 'C4', ('cellulose', 'glucomannan'), (19, 243), True, None),
        ('AA9', 'Ls(AA9)A', 'A0A0S2GKZ1', 'C4', ('cellulose', 'xyloglucan', 'glucomannan', 'xylan'), (19, 254), False, ('5acf', '5acg', '5ach', '5aci', '5acj', '5n04', '5n05', '5nkw', '5nln', '5nlo', '5nlp', '5nlq', '5nlr', '5nls', '6ydg', '7nim', '7nin', '7pqr', '7ptz', '7pxi', '7pxj', '7pxk', '7pxl', '7pxm', '7pxn', '7pxr', '7pxs', '7pxt', '7pxu', '7pxv', '7pxw', '7pyd', '7pye', '7pyf', '7pyg', '7pyh', '7pyi', '7pyl', '7pym', '7pyn', '7pyo', '7pyp', '7pyq', '7pyu', '7pyw', '7pyx', '7pyy', '7pyz', '7py0', '8e1W')),
        ('AA9', 'MtLPMO9A', 'A0A0H4K9X4', 'C1/C4', ('cellulose', 'xyloglucan', 'xylan'), (17, 235), False, None),
        ('AA9', 'MtLPMO9B', 'A0A1C9CXI1', 'C1', ('cellulose'), (18, 267), True, None),
        ('AA9', 'MtLPMO9C', 'A0A1C9CXI0', 'C4', ('cellulose', 'xyloglucan'), (15, 237), False, None),
        ('AA9', 'MtLPMO9D', 'A0A218MJF1', 'C1', ('cellulose'), (17, 245), False, None),
        ('AA9', 'MtLPMO9E/J', 'G2Q7A5', 'C4', ('cellulose', 'xyloglucan'), (17, 246), False, None),
        ('AA9', 'MtLPMO9F', 'G2Q9F7', 'C4', ('cellulose', 'xyloglucan', 'glucomannan'), (15, 241), True, None),         
        ('AA9', 'MtLPMO9G', 'G2QK49', 'C4', ('cellulose'), (15, 235), False, None),
        ('AA9', 'MtLPMO9H', 'G2Q9T3', 'C1/C4', ('cellulose', 'xyloglucan'), (19, 265), True, None),
        ('AA9', 'MtLPMO9I', 'G2Q774', 'C1', ('cellulose'), (20, 242), False, None),
        ('AA9', 'MtLPMO9L', 'G2QJT0', 'C1', ('cellulose'), (17, 224), False, None),
        ('AA9', 'MYCTH_112089', 'G2QI82', 'C1', ('cellulose'), (17, 232), False, None),
        ('AA9', 'NCU00836', 'Q7SCJ5', 'C1', ('cellulose'), (18, 220), True, None),
        ('AA9', 'NcLPMO9D', 'Q1K8B6', 'C4', ('cellulose'), (15, 238), False, ('4eir', '7t5c', '7t5e')),
        ('AA9', 'NcLPMO9J', 'Q7SHD9', 'C1', ('cellulose'), (20, 265), True, None),
        ('AA9', 'NcLPMO9A', 'Q7S439', 'C4', ('cellulose'), (15, 237), True, ('5foh')),
        ('AA9', 'NcLPMO9C', 'Q7SHI8', 'C4', ('cellulose', 'xyloglucan', 'glucomannan'), (16, 243), True, ('4d7u', '4d7v')),
        ('AA9', 'NcLPMO9F', 'Q1K4Q1', 'C1', ('cellulose'), (17, 231), False, ('4qi8')),
        ('AA9', 'NCU07760', 'Q7S111', 'C1/C4', ('cellulose'), (18, 265), True, None),
        ('AA9', 'NcLPMO9M', 'Q7SA19', 'C1', ('cellulose'), (16, 238), False, ('4eis')),
        ('AA9', 'NcLPMO9E', 'Q7RWN7', 'C1', ('cellulose'), (20, 271), True, None),
        ('AA9', 'PsLPMO9A', 'A0A167KNZ4', 'C1/C4', ('cellulose'), (0, 213), False, None),
        ('AA9', 'PsLPMO9B', 'A0A167KNY8', 'C4', ('cellulose'), (0, 222), False, None),
        ('AA9', 'PcLPMO9D', 'H1AE14', 'C1', ('cellulose'), (18, 235), False, ('4b5q')),
        ('AA9', 'PaLPMO9A', 'B2B629', 'C1/C4', ('cellulose'), (20, 269), True, None),
        ('AA9', 'PaLPMO9B', 'B2AVF1', 'C1/C4', ('cellulose'), (19, 259), True, None),
        ('AA9', 'PaLPMO9D', 'B2ARG6', 'C1', ('cellulose'), (21, 269), False, None),
        ('AA9', 'PaLPMO9E', 'B2ATL7', 'C1', ('cellulose'), (19, 222), True, None),
        ('AA9', 'PaLPMO9F', 'B2B403', 'C1', ('cellulose'), (16, 253), False, None),
        ('AA9', 'PaLPMO9H', 'B2ADG1', 'C1/C4', ('cellulose', 'xyloglucan', 'glucomannan'), (15, 243), True, None),
        ('AA9', 'TaLPMO9A', 'G3XAP7', 'C1/C4', ('cellulose', 'xylan'), (0, 228), False, ('2yet', '3zud', '7pu1', '7pz3', '7pz4', '7pz5', '7pz6', '7pz7', '7pz8', '7q1k')),
        ('AA9', 'HjLPMO9A', 'O14405', 'C1/C4', ('cellulose'), (21, 268), True, None),
        ('AA9', 'TtLPMO9E', 'G2RGE5', 'C1', ('cellulose'), (18, 226), False, None),

        # AA10 
        ('AA10', 'AsLPMO10A', 'B6EQB6', 'C1', ('chitin'), (24, 205), True, None),
        ('AA10', 'AsLPMO10B', 'B6EQJ6', 'C1', ('chitin'), (25, 217), True, ('7okr')),
        ('AA10', 'BaAA10A', 'Q9F9Q5', 'C1', ('chitin'), (27, 206), False, None),
        ('AA10', 'BlAA10A', 'D0EW67', 'C1', ('chitin'), (54, 226), False, None),
        ('AA10', 'BtLPMO10A-FL', 'D0EW65', 'C1', ('chitin'), (40, 212), True, None),
        ('AA10', 'BtLPMO10A', 'A0A0C5K362', 'C1', ('chitin'), (34, 202), False, None),
        ('AA10', 'CfiLPMO10', 'F4H6A3', 'C1/C4', ('cellulose'), (33, 218), True, None),
        ('AA10', 'CflaLPMO10A', 'D5UGB1', 'C1', ('cellulose'), (33, 229), True, None),
        ('AA10', 'CflaLPMO10B', 'D5UGA8', 'C1/C4', ('cellulose'), (35, 221), True, None),
        ('AA10', 'CflaLPMO10C', 'D5UH31', 'C1/C4', ('cellulose'), (36, 219), True, None),
        ('AA10', 'CflaLPMO10D', 'D5UHY1', 'C1', ('chitin'), (42, 190), True, None),
        ('AA10', 'CjLPMO10A', 'B3PJ79', 'C1', ('chitin'), (36, 218), True, ('5fqj', '6z40', '6z41')),
        ('AA10', 'CjLPMO10B', 'B3PDT6', 'C1', ('cellulose'), (24, 230), True, None),
        ('AA10', 'EfAA10A', 'Q838S1', 'C1', ('chitin'), (28, 194), False, ('4a02', '4alc', '4ale', '4alq', '4alr', '4als', '4alt')),
        ('AA10', 'HcAA10-2', 'Q2SNS3', 'C1', ('cellulose'), (23, 229), True, None),
        ('AA10', 'JdLPMO10A', 'C7R4I0', 'C1', ('chitin'), (31, 173), True, ('5aa7', '5vg0', '5vg1')),
        ('AA10', 'KpLPMO10A', 'A0A4P7DN87', 'C1/C4', ('cellulose', 'chitin', 'xylan'), (38, 224), False, None), #C1 on chitin
        ('AA10', 'LmLPMO10', 'Q8Y4H4', 'C1', ('chitin'), (27, 193), True, None),
        ('AA10', 'Plu2352', 'Q7N4I5', 'C1', ('chitin'), (25, 199), False, ('6t5z')),
        ('AA10', 'SmAA10A', 'O83009', 'C1', ('chitin'), (27, 193), False, ('2bem', '2ben', '2lhs')),
        # ('AA10', '22kDa protein', 'Q59930', 'C1', ('chitin'), (27, 171), False, None), # AWFUL STRUCTURE
        ('AA10', 'SamLPMO10B', 'A3KIM2', 'C1', ('chitin'), (30, 172), False, None),
        ('AA10', 'SamLPMO10C', 'A3KKC4', 'C1', ('cellulose'), (34, 238), True, None),
        ('AA10', 'ScLPMO10A', 'Q9S296', 'C1', ('chitin'), (33, 213), True, None),
        ('AA10', 'ScLPMO10B', 'Q9RJC1', 'C1/C4', ('cellulose'), (42, 228), False, ('4oy6', '4oy8')),
        ('AA10', 'ScLPMO10C', 'Q9RJY2', 'C1', ('cellulose'), (34, 238), True, ('4oy7', '6f7e')),
        ('AA10', 'SgLPMO10A', 'B1VNK5', 'C1', ('cellulose', 'chitin'), (34, 238), True, None),
        ('AA10', 'SgLPMO10F', 'B1VN59', 'C1', ('chitin'), (30, 171), False, None),
        ('AA10', 'SliLPMO10E', 'A0A7U9DRA2', 'C1', ('chitin'), (29, 201), False, None),
        ('AA10', 'TtAA10A', 'C5BKQ9', 'C1/C4', ('cellulose'), (24, 228), True, ('6rw7')),
        ('AA10', 'TfLPMO10A', 'Q47QG3', 'C1/C4', ('cellulose'), (36, 222), False, ('5uiz')),
        ('AA10', 'TfAA10B', 'Q47PB9', 'C1', ('cellulose'), (31, 223), True, None),
        ('AA10', 'VcAA10B', 'A9Y370', 'C1', ('chitin'), (23, 202), True, None),
        # ('AA10', 'LpsAA10A', 'A0A8K1XN59', 'C1', ('cellulose')), # no structure
        ('AA10', 'lytic chitin monooxygenase', 'W5QLL4', 'C1', ('chitin'), (24, 206), False, ('6if7')),
        # ('AA10', 'fusolin (ACV034)', 'O70709', 'C1', ('chitin'), (16, 246), False, ('4yn1')), no AF structure

        # AA11
        ('AA11', 'AfAA11A', 'B0XZD3', 'C1', ('chitin'), (18, 219), False, ('7p3u')),
        ('AA11', 'A0AA11', 'Q2UA85', 'C1', ('chitin'), (19, 235), True, ('4mah', '4mai')),

        # AA13
        ('AA13', 'AnAA13', 'A0A1U8QN05', 'C1', ('chitin'), (18, 254), True, None),
        ('AA13', 'A0AA13', 'Q2U8Y3', 'C1', ('chitin'), (46, 279), False, ('4opb', '5lsv', '5t7j', '5t7k', '5t7n', '6tbq', '6tbr', '6tc4')),
        ('AA13', 'AtLPMO13A', 'Q0CGA6', 'C1', ('chitin'), (17, 250), True, None),
        ('AA13', 'NcAA13', 'Q6MWQ3', 'C1', ('chitin'), (18, 273), True, None),

        # AA14
        ('AA14', 'PcAA14A', 'A0A2I6QB00', 'C1', ('xylan'), (18, 299), False, None),
        ('AA14', 'PcAA14B', 'A0A2I6QAZ5', 'C1', ('xylan'), (21, 290), False, None),

        # AA15
        ('AA15', 'AaAA15', 'W4GRV7', 'C1', ('chitin'), (19, 218), True, None),
        ('AA15', 'TdAA15A', 'A0A2N8U5K6', 'C1', ('chitin', 'cellulose'), (20, 213), False, None),
        ('AA15', 'TdAA15B', 'A0A2N8U5I8', 'C1', ('chitin', 'cellulose'), (19, 212), False, None),

        # AA16
        ('AA16', 'AaAA16', 'A0A1L9X7U6', 'C1', ('cellulose'), (19, 192), False, None),

        # AA17
        ('AA17', 'PiAA17A', 'A0A833WPX8', 'C4', ('homogalacturonan'), (23, 197), False, None),
        ('AA17', 'PiAA17A', 'D0N2F6', 'C4', ('homogalacturonan'), (23, 195), False, None),
        ('AA17', 'PiAA17C', 'D0N2F7', 'C4', ('homogalacturonan'), (23, 194), False, ('6z5y'))
    ]


    # Make data frame
    columns = (
        'Family', 'Name', 'UniProt', 'Regioselectivity',\
        'Substrate specificity', 'Core', 'CBM', 'PDB'
        )

    df = pd.DataFrame(dict(zip(columns, zip(*db))))

    # Store data frame
    df.to_pickle(f"{config['supervised']}/labels.pkl")

    return df

def download(db : pd.DataFrame):
    '''
    Parameters
    ----------
    db : pd.DataFrame
        Characterized LPMOs
    
    Returns
    -------
    None
    '''

    for uniprot, core in tqdm(zip(db['UniProt'], db['Core']), total = len(db)):

        # Download AlphaFold file
        version = 'v4'
        model_url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_{version}.pdb'
        os.system(f"curl {model_url} -o {config['AF_labels']}/{uniprot}.pdb")

        # Crop to core
        structure = AlphaFoldStructure(f"{config['AF_labels']}/{uniprot}.pdb")
        structure.rewrite_range(
            f"{config['AF_labels_core']}/all/{uniprot}.pdb", 
            core
            )

def C1_vs_non_C1(regioselectivty : str) -> int:
    '''
    Converts the three-class regioselectivity problem into a binary problem by
    assigning 1s to C1-degrading LPMOs and 0s to C4 and C1/C4-degrading LPMOs

    Parameters
    ----------
    regioselectivity : str
        LPMO regioselectivity: either C1, C4 or C1/C4
    
    Returns
    -------
    1 if C1, 0 else
    '''

    if regioselectivty == 'C1':
        return 1
    else: 
        return 0

def C4_vs_non_C4(regioselectivty : str) -> int:
    '''
    Converts the three-class regioselectivity problem into a binary problem by
    assigning 1s to C4-degrading LPMOs and 0s to C1 and C1/C4-degrading LPMOs

    Parameters
    ----------
    regioselectivity : str
        LPMO regioselectivity: either C1, C4 or C1/C4
    
    Returns
    -------
    1 if C4, 0 else
    '''

    if regioselectivty == 'C4':
        return 1
    else: 
        return 0

def C1_vs_C4(regioselectivty : str) -> Union[int, None]:
    '''
    Converts the three-class regioselectivity problem into a binary problem by
    assigning 1s to C1-degrading LPMOs and 0s to C4-degrading LPMOs, 
    disregarding C1/C4-degrading enzymes.

    Parameters
    ----------
    regioselectivity : str
        LPMO regioselectivity: either C1, C4 or C1/C4
    
    Returns
    -------
    1 if C1, 0 if C4, None else.
    '''

    if regioselectivty == 'C1':
        return 1
    elif regioselectivty == 'C4':
        return 0
    else: 
        return None

def cellulose_vs_chitin(ss : str) -> Union[str, None]:
    '''
    Converts the multiclass regioselectivity problem into a binary problem by
    assigning 1s to cellulose-degrading LPMOs and 0s to chitin-degrading LPMOs,
    disregarding non-cellulose- and non-chitin degrading enzymes.

    Parameters
    ----------
    ss : str
        LPMO substrate specificity: it can be chitin, cellulose, xyloglycan, 
        homogelacturonan, etc.
    
    Returns
    -------
    1 if cellulose, 0 if chitin, None if both or none of them.
    '''
    if 'cellulose' in ss and 'chitin' in ss:
        return None
    elif 'cellulose' in ss:
        return 1
    elif 'chitin' in ss:
        return 0
    else: 
        return None

def prepare_labels_for_ML(db : pd.DataFrame) -> pd.DataFrame:
    '''
    Prepare the characterized LPMO database for the use of machine learning 
    algorithms. This is done by assigning 1s and 0s to the different classes of
    a problem, removing the entries that will not be used by a specific problem,
    and storing the results independently.

    Parameters
    ----------
    db : pd.DataFrame
        Characterized LPMOs
    
    Returns
    -------
    None.
    '''

    # Problem : binarizer function
    filters = {
        'C1' : ('Regioselectivity', C1_vs_non_C1),
        'C4' : ('Regioselectivity', C4_vs_non_C4),
        'C1C4' : ('Regioselectivity', C1_vs_C4),
        'SS' : ('Substrate specificity', cellulose_vs_chitin)
    }

    # Iterate over all the problems available
    for problem, (column, binarizer) in filters.items():

        # Apply binarizer function
        db['y'] = db[column].apply(lambda x : binarizer(x))

        # Remove rows that are not needed (None)
        db_dropnaed = db.dropna(subset = 'y')

        # Copy AF files to specific directory (e.g. all --> SS)
        for uniprot in db_dropnaed['UniProt']:
            os.system(f"cp {config['AF_labels_core']}/all/{uniprot}.pdb \
                  {config['AF_labels_core']}/{problem}")
            
        # Store labels database
        db_dropnaed.to_pickle(f"{config['supervised']}/{problem}.pkl")

def characterized_vs_all() -> set:
    '''
    Takes the UniProt IDs of the LPMO entries from CAZy and compares them with 
    those of the characterized entries. This ensures that all the characterized
    LPMOs are taken into account for subsequent analysis.

    Parameters
    ----------
    None

    Returns
    -------
    not_in_cazy : set
        Entries present in the characterized database but not in CAZy
    '''
    labels = pd.read_pickle(f"{config['supervised']}/labels.pkl")
    all = pd.read_pickle(f"{config['CAZy_expanded']}/all.pkl")
    unique_labels = set(labels['UniProt'])
    unique_all = set(all['UniProt'])
    not_in_cazy = unique_labels - unique_all

    return not_in_cazy

def add_new_entries():
    ''' 
    Use the ManualEntries class to update the data frame with all the
    retrieved information of LPMOs

    Returns
    -------
    None
    '''

    # Create ManualEntries object
    manual = ManualEntries(
        families = ['AA9', 'AA9', 'AA9', 'AA16'],
        domains = ['Eukaryota', 'Eukaryota', 'Eukaryota', 'Eukaryota'],
        uniprot_ids = ['W4K498', 'A0A0J9XL55', 'W4K8M0', 'A0A1L9X7U6']
    )

    # Add information
    manual.add_genbank()
    manual.add_species()
    manual.add_protein_fasta()
    manual.add_dna_fasta()
    manual.add_pdb()

    # Append to supreme df
    all = pd.read_pickle(f"{config['CAZy_expanded']}/all.pkl")
    all_v2 = manual.add_manual_entries(all)
    all_v2.to_pickle(f"{config['CAZy_expanded']}/all.pkl")

def plot_regioselectivity(db: pd.DataFrame):
    '''
    Plot regioselectivity counts. Because all chitin-degrading enzymes have
    C1 regioselectivity, a better overview of the regioselectivity counts 
    would be to remove these LPMOs.
    FIX: plot1 == plot2 but counts_chitin != counts_no_chitin, so something
         is not working with the pandas plot method.

    Parameters
    ----------
    db : pd.DataFrame
        Characterized LPMOs

    Returns
    -------
    None
    '''

    # With chitin
    counts_chitin = db['Regioselectivity'].value_counts()
    plot1 = counts_chitin.plot(
        kind = 'bar', 
        rot = 'horizontal', 
        title = 'Regioselectivity counts (with chitin)',
        xlabel = 'Regioselectivity',
        ylabel = 'Number of entries'
        )
    fig1 = plot1.get_figure()
    fig1.savefig(
        f"{config['plots']}/regio_counts_chitin.png", 
        transparent = True
        )
    
    # Without chitin
    db_no_chitin = db[db['Substrate specificity'] != ('chitin')]
    counts_no_chitin = db_no_chitin['Regioselectivity'].value_counts()
    plot2 = counts_no_chitin.plot(
        kind = 'bar', 
        rot = 'horizontal', 
        title = 'Regioselectivity counts (without chitin)',
        xlabel = 'Regioselectivity',
        ylabel = 'Number of entries'
        )
    fig2 = plot2.get_figure()
    fig2.savefig(
        f"{config['plots']}/regio_counts_no_chitin.png", 
        transparent = True
        )

def plot_substrate_specificity(db: pd.DataFrame):
    '''
    Plot substrate specificity counts. 

    Parameters
    ----------
    db : pd.DataFrame
        Characterized LPMOs
        
    Returns
    -------
    None
    '''

    counts = db['Substrate specificity'].explode().value_counts()
    plot = counts.plot(
        kind = 'bar', 
        fontsize = 7, 
        rot = 'horizontal',
        title = 'Substrate specificity counts',
        xlabel = 'Substrate specificity',
        ylabel = 'Number of entries'
        )
    fig = plot.get_figure()
    fig.savefig(
        f"{config['plots']}/substrate_counts.png", 
        transparent = True
        )

def plot_per_family(db: pd.DataFrame):
    '''
    Plot regioselectivity and substrate specificity per family.

    Parameters
    ----------
    db : pd.DataFrame
        Characterized LPMOs
        
    Returns
    -------
    None
    '''

    # Plot regioselectivity per family
    regio_counts = db[['Family', 'Regioselectivity']].value_counts()
    lpmo_sort = lambda index: [int(family[2:]) for family in index]
    regio_per_family = regio_counts.unstack().sort_index(key = lpmo_sort)
    plot1 = regio_per_family.plot(
        kind ='bar', 
        stacked = True,
        rot = 'horizontal',
        title = 'Regioselectivity per family',
        xlabel = 'Family',
        ylabel = 'Number of entries'
        )
    fig1 = plot1.get_figure()
    fig1.savefig(
        f"{config['plots']}/regio_per_family.png", 
        transparent = True
        )
    
    # Plot substrate specificity per family
    substrate_counts = db[
        ['Family', 
         'Substrate specificity']
         ].explode('Substrate specificity').value_counts()
    lpmo_sort = lambda index: [int(family[2:]) for family in index]
    substrate_per_family = substrate_counts.unstack().sort_index(key=lpmo_sort)
    plot1 = substrate_per_family.plot(
        kind ='bar', 
        stacked = True,
        rot = 'horizontal',
        title = 'Substrate specificity per family',
        xlabel = 'Family',
        ylabel = 'Number of entries'
        )
    fig1 = plot1.get_figure()
    fig1.savefig(
        f"{config['plots']}/substrate_per_family.png", 
        transparent = True
        )

def main():
    '''Program flow.'''

    # Labels database
    db = make_database()

    # Download labels
    download(db)

    # Prepare for machine learning
    prepare_labels_for_ML(db)

    # Is there any new LPMO not in CAZy or BRENDA?
    new_entries = characterized_vs_all()

    # Add new entries
    if new_entries:
        print('Adding new entries:', new_entries)
        add_new_entries()

    # Plot statistics
    plot_regioselectivity(db)
    plot_substrate_specificity(db)
    plot_per_family(db)   

if __name__ == '__main__':
    main()
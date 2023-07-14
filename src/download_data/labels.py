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

import pandas as pd
import matplotlib.pyplot as plt
from .manual_entries import ManualEntries
from ..config.config_parser import config


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
        ('AA9', 'AfAA9B', 'Q4WP32', 'C1/C4', ('cellulose')),
        ('AA9', 'AN3046', 'A0A1U8QMG7', 'C1', ('cellulose', 'xyloglucan')),
        ('AA9', 'AN1602.2', 'Q5BCX8', 'C4', ('cellulose')),
        ('AA9', 'CsLPMO9A', 'A0A6C0PI79', 'C1', ('cellulose')),
        ('AA9', 'CsLPMO9', 'M2RAI8', 'C1', ('cellulose')),
        ('AA9', 'CsLPMO9B', 'M2QEW2', 'C1', ('cellulose')),
        ('AA9', 'CvAA9A', 'A0A223GEC9', 'C1/C4', ('cellulose', 'xyloglucan', 'glucomannan')), 
        ('AA9', 'FgLPMO9A', 'I1REU9', 'C1/C4', ('cellulose', 'xyloglucan')),
        ('AA9', 'GtLPMO9A-2', 'A0A1C9ZMC5', 'C1/C4', ('cellulose', 'xyloglucan', 'glucomannan')),
        ('AA9', 'GtLPMO9B', 'F8T947', 'C1/C4', ('cellulose')),
        # ('AA9', 'GcLPMO9A', 'A0A0J9XKT0', 'C1/C4', ('cellulose', 'xyloglucan')), # no structure
        ('AA9', 'GcLPMO9B', 'A0A0J9XL55', 'C1/C4', ('cellulose', 'xyloglucan')), # add
        ('AA9', 'HiLPMO9B', 'W4KMP1', 'C1', ('cellulose')),
        ('AA9', 'HiLPMO9H', 'W4K498', 'C1', ('cellulose')), # add
        ('AA9', 'HiLPMO9I', 'W4K8M0', 'C4', ('cellulose', 'glucomannan')), # add
        ('AA9', 'Ls(AA9)A', 'A0A0S2GKZ1', 'C4', ('cellulose', 'xyloglucan', 'glucomannan', 'xylan')),
        ('AA9', 'MtLPMO9A', 'A0A0H4K9X4', 'C1/C4', ('cellulose', 'xyloglucan', 'xylan')),
        ('AA9', 'MtLPMO9B', 'A0A1C9CXI1', 'C1', ('cellulose')),
        ('AA9', 'MtLPMO9C', 'A0A1C9CXI0', 'C4', ('cellulose', 'xyloglucan')),
        ('AA9', 'MtLPMO9D', 'A0A218MJF1', 'C1', ('cellulose')),
        ('AA9', 'MtLPMO9E/J', 'G2Q7A5', 'C4', ('cellulose', 'xyloglucan')),
        ('AA9', 'MtLPMO9F', 'G2Q9F7', 'C4', ('cellulose', 'xyloglucan', 'glucomannan')),         
        ('AA9', 'MtLPMO9G', 'G2QK49', 'C4', ('cellulose')),
        ('AA9', 'MtLPMO9H', 'G2Q9T3', 'C1/C4', ('cellulose', 'xyloglucan')),
        ('AA9', 'MtLPMO9I', 'G2Q774', 'C1', ('cellulose')),
        ('AA9', 'MtLPMO9L', 'G2QJT0', 'C1', ('cellulose')),
        ('AA9', 'MYCTH_112089', 'G2QI82', 'C1', ('cellulose')),
        ('AA9', 'NCU00836', 'Q7SCJ5', 'C1', ('cellulose')),
        ('AA9', 'NcLPMO9D', 'Q1K8B6', 'C4', ('cellulose')),
        ('AA9', 'NcLPMO9J', 'Q7SHD9', 'C1', ('cellulose')),
        ('AA9', 'NcLPMO9A', 'Q7S439', 'C4', ('cellulose')),
        ('AA9', 'NcLPMO9C', 'Q7SHI8', 'C4', ('cellulose', 'xyloglucan', 'glucomannan')),
        ('AA9', 'NcLPMO9F', 'Q1K4Q1', 'C1', ('cellulose')),
        ('AA9', 'NCU07760', 'Q7S111', 'C1/C4', ('cellulose')),
        ('AA9', 'NcLPMO9M', 'Q7SA19', 'C1', ('cellulose')),
        ('AA9', 'NcLPMO9E', 'Q7RWN7', 'C1', ('cellulose')),
        ('AA9', 'PsLPMO9A', 'A0A167KNZ4', 'C1/C4', ('cellulose')),
        ('AA9', 'PsLPMO9B', 'A0A167KNY8', 'C4', ('cellulose')),
        ('AA9', 'PcLPMO9D', 'H1AE14', 'C1', ('cellulose')),
        ('AA9', 'PaLPMO9A', 'B2B629', 'C1/C4', ('cellulose')),
        ('AA9', 'PaLPMO9B', 'B2AVF1', 'C1/C4', ('cellulose')),
        ('AA9', 'PaLPMO9D', 'B2ARG6', 'C1', ('cellulose')),
        ('AA9', 'PaLPMO9E', 'B2ATL7', 'C1', ('cellulose')),
        ('AA9', 'PaLPMO9F', 'B2B403', 'C1', ('cellulose')),
        ('AA9', 'PaLPMO9H', 'B2ADG1', 'C1/C4', ('cellulose', 'xyloglucan', 'glucomannan')),
        ('AA9', 'TaLPMO9A', 'G3XAP7', 'C1/C4', ('cellulose', 'xylan')),
        ('AA9', 'HjLPMO9A', 'O14405', 'C1/C4', ('cellulose')),
        ('AA9', 'TtLPMO9E', 'G2RGE5', 'C1', ('cellulose')),

        # AA10 
        ('AA10', 'AsLPMO10A', 'B6EQB6', 'C1', ('chitin')),
        ('AA10', 'AsLPMO10B', 'B6EQJ6', 'C1', ('chitin')),
        ('AA10', 'BaAA10A', 'Q9F9Q5', 'C1', ('chitin')),
        ('AA10', 'BlAA10A', 'D0EW67', 'C1', ('chitin')),
        ('AA10', 'BtLPMO10A-FL', 'D0EW65', 'C1', ('chitin')),
        ('AA10', 'BtLPMO10A', 'A0A0C5K362', 'C1', ('chitin')),
        ('AA10', 'CfiLPMO10', 'F4H6A3', 'C1/C4', ('cellulose')),
        ('AA10', 'CflaLPMO10A', 'D5UGB1', 'C1', ('cellulose')),
        ('AA10', 'CflaLPMO10B', 'D5UGA8', 'C1/C4', ('cellulose')),
        ('AA10', 'CflaLPMO10C', 'D5UH31', 'C1/C4', ('cellulose')),
        ('AA10', 'CflaLPMO10D', 'D5UHY1', 'C1', ('chitin')),
        ('AA10', 'CjLPMO10A', 'B3PJ79', 'C1', ('chitin')),
        ('AA10', 'CjLPMO10B', 'B3PDT6', 'C1', ('cellulose')),
        ('AA10', 'EfAA10A', 'Q838S1', 'C1', ('chitin')),
        ('AA10', 'HcAA10-2', 'Q2SNS3', 'C1', ('cellulose')),
        ('AA10', 'JdLPMO10A', 'C7R4I0', 'C1', ('chitin')),
        ('AA10', 'KpLPMO10A', 'A0A4P7DN87', 'C1/C4', ('cellulose', 'chitin', 'xylan')), #C1 on chitin
        ('AA10', 'LmLPMO10', 'Q8Y4H4', 'C1', ('chitin')),
        ('AA10', 'plu2352', 'Q7N4I5', 'C1', ('chitin')),
        ('AA10', 'SmAA10A', 'O83009', 'C1', ('chitin')),
        ('AA10', '22kDa protein', 'Q59930', 'C1', ('chitin')),
        ('AA10', 'SamLPMO10B', 'A3KIM2', 'C1', ('chitin')),
        ('AA10', 'SamLPMO10C', 'A3KKC4', 'C1', ('cellulose')),
        ('AA10', 'ScLPMO10A', 'Q9S296', 'C1', ('chitin')),
        ('AA10', 'ScLPMO10B', 'Q9RJC1', 'C1/C4', ('cellulose')),
        ('AA10', 'ScLPMO10C', 'Q9RJY2', 'C1', ('cellulose')),
        ('AA10', 'SgLPMO10A', 'B1VNK5', 'C1', ('cellulose', 'chitin')),
        ('AA10', 'SgLPMO10F', 'B1VN59', 'C1', ('chitin')),
        ('AA10', 'SliLPMO10E', 'A0A7U9DRA2', 'C1', ('chitin')),
        ('AA10', 'TtAA10A', 'C5BKQ9', 'C1/C4', ('cellulose')),
        ('AA10', 'TfLPMO10A', 'Q47QG3', 'C1/C4', ('cellulose')),
        ('AA10', 'TfAA10B', 'Q47PB9', 'C1', ('cellulose')),
        ('AA10', 'VcAA10B', 'A9Y370', 'C1', ('chitin')),
        ('AA10', 'LpsAA10A', 'A0A8K1XN59', 'C1', ('cellulose')),
        ('AA10', 'lytic chitin monooxygenase', 'W5QLL4', 'C1', ('chitin')),
        ('AA10', 'fusolin (ACV034)', 'O70709', 'C1', ('chitin')),

        # AA11
        ('AA11', 'AfAA11A', 'B0XZD3', 'C1', ('chitin')),
        ('AA11', 'A0AA11', 'Q2UA85', 'C1', ('chitin')),

        # AA13
        ('AA13', 'AnAA13', 'A0A1U8QN05', 'C1', ('chitin')),
        ('AA13', 'A0AA13', 'Q2U8Y3', 'C1', ('chitin')),
        ('AA13', 'AtLPMO13A', 'Q0CGA6', 'C1', ('chitin')),
        ('AA13', 'NcAA13', 'Q6MWQ3', 'C1', ('chitin')),

        # AA14
        ('AA14', 'PcAA14A', 'A0A2I6QB00', 'C1', ('xylan')),
        ('AA14', 'PcAA14B', 'A0A2I6QAZ5', 'C1', ('xylan')),

        # AA15
        ('AA15', 'AaAA15', 'W4GRV7', 'C1', ('chitin')),
        ('AA15', 'TdAA15A', 'A0A2N8U5K6', 'C1', ('chitin', 'cellulose')),
        ('AA15', 'TdAA15B', 'A0A2N8U5I8', 'C1', ('chitin', 'cellulose')),

        # AA16
        ('AA16', 'AaAA16', 'A0A1L9X7U6', 'C1', ('cellulose')), #add

        # AA17
        ('AA17', 'PiAA17A', 'A0A833WPX8', 'C4', ('homogalacturonan')),
        ('AA17', 'PiAA17A', 'D0N2F6', 'C4', ('homogalacturonan')),
        ('AA17', 'PiAA17C', 'D0N2F7', 'C4', ('homogalacturonan')),
    ]

    # Make data frame
    columns = ['Family', 'Name', 'UniProt', 'Regioselectivity', 'Substrate specificity']
    df = pd.DataFrame(dict(zip(columns, zip(*db))))

    # Store data frame
    df.to_pickle(f"{config['supervised']}/labels.pkl")

    return df

def characterized_vs_all() -> set:
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
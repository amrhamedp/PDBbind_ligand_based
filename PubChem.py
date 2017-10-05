import os
import requests
import pandas as pd
from lxml import etree
import pubchempy


def get_AIDs(gid):
    """Use Gene ID number to get list of PubChem Assay ID"""

    req = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/geneid/%s/aids/TXT' % gid)
    aids = req.text.split()

    return aids


def download_aid_csv(aid, directory):
    """Download information from PubChem to directory"""

    f = os.path.join(directory, '%s.csv' % aid)
    if not os.path.isfile(f):
        os.system('wget -qO %s https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/CSV' % (f, aid))


def gene_id_to_data_frame(gene_id, aids, bioactivity_type):
    """Prepare data frame for target

    Parameters
    ----------
    gene_id: str
        Gene ID.

    Returns
    -------
    data: DataFrame
        Columns: gene_id, puchem_aid, pubchem_sid, bioactivity, value.

    """
    bioactivity = bioactivity_type + ' standard value'
    qualifier = bioactivity_type + ' qualifier'

    gene_dict = {'gene_id': [],
                 'pubchem_aid': [],
                 'pubchem_cid': [],
                 'pubchem_sid': [],
                 'smiles': [],
                 'bioactivity': [],
                 'qualifier': [],
                 'value': [],
                 'unit': []}

    for aid in aids:

        try:
            data = pd.read_csv('./aid_files/%s.csv' % aid, index_col=0)

            if 'PUBCHEM_SID' in data.columns and bioactivity in data.columns and\
               qualifier in data.columns and 'RESULT_UNIT' in data.index and '1' in data.index:
                for idx in list(data['1':].index):
                    gene_dict['gene_id'].append(gene_id)
                    gene_dict['pubchem_aid'].append(aid)
                    gene_dict['pubchem_cid'].append(data['PUBCHEM_CID'][idx])
                    gene_dict['pubchem_sid'].append(data['PUBCHEM_SID'][idx])
                    gene_dict['bioactivity'].append(bioactivity_type)
                    gene_dict['qualifier'].append(data[qualifier][idx])
                    gene_dict['value'].append(data[bioactivity][idx])

                    smiles = get_smiles(str(int(data['PUBCHEM_CID'][idx])))
                    print(smiles)
                    if smiles:
                        gene_dict['smiles'].append(smiles)
                    else:
                        gene_dict['smiles'].append('Unspecified')
                        print('blop')

                    if data[bioactivity]['RESULT_UNIT'] == 'MICROMOLAR':
                        gene_dict['unit'].append('uM')
                    elif data[bioactivity]['RESULT_UNIT'] == 'NANOMOLAR':
                        gene_dict['unit'].append('nM')
                    elif data[bioactivity]['RESULT_UNIT'] == 'PICOMOLAR':
                        gene_dict['unit'].append('pM')
                    elif data[bioactivity]['RESULT_UNIT'] == '%':
                        gene_dict['unit'].append('%')
                    else:
                        gene_dict['unit'].append('Unspecified')
        except:
            continue

    return pd.DataFrame(gene_dict)[['gene_id', 'pubchem_aid', 'pubchem_cid', 'pubchem_sid',
                                    'smiles', 'bioactivity', 'qualifier', 'value', 'unit']]



def create_data_frame(gene_id, aids, bioactivity_types, directory='', overwrite=False):
    """Create and save data frames to csv

    Parameters
    ----------
    gene_ids: list of stings
        List of Gene IDs.

    bioactivity_types: list of strings
        List of bioactivities.

    directory: str
        Directory where files will be saved.

    overwrite: bool
        Overwrites your files if True.

    """


    f = os.path.join(directory, 'pubchem_%s.csv' % gene_id)
    if overwrite or not os.path.isfile(f):

        data = pd.concat(gene_id_to_data_frame(gene_id, aids, bioactivity)
                         for bioactivity in bioactivity_types)
        data.to_csv(f, index=False)


def get_smiles(cid):
    """Get smiles from PubChem using Compound ID"""

    comp = pubchempy.Compound.from_cid(cid)
    return comp.canonical_smiles


def convert_unit(df, old_unit, new_unit, factor=1):
    """Replace old unit with new unit in bioactives dataframe
    Parameters
    ----------
    df : DataFrame
        DataFrame created by gene_id_to_data_frame function.
    old_unit : str
        Unit you want to convert.
    new_unit : str
        Unit to convert to.
    factor : float
        Factor to convert units.
    """

    idxs = df[df['unit'] == old_unit].index
    for idx in idxs:
        val = df['value'][idx]
        df.set_value(index=idx, col='value', value=val * factor)
        df.set_value(index=idx, col='unit', value=new_unit)

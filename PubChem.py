import os
import requests
import pandas as pd
from lxml import etree


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


def gene_id_to_data_frame(gene_id, aids):
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

    gene_dict = {'gene_id': [],
                 'pubchem_aid': [],
                 'pubchem_sid': [],
                 'bioactivity': [],
                 'value': []}

    for aid in aids:

        try:
            data = pd.read_csv('./aid_files/%s.csv' % aid, index_col=0)

            if 'PUBCHEM_SID' in data.columns and 'IC50' in data.columns and '1' in data.index:
                gene_dict['gene_id'].append(gene_id)
                gene_dict['pubchem_aid'].append(aid)
                gene_dict['pubchem_sid'].append(data['PUBCHEM_SID']['1'])
                gene_dict['bioactivity'].append('IC50')
                gene_dict['value'].append(data['IC50']['1'])
        except:
            continue

    return pd.DataFrame(gene_dict)

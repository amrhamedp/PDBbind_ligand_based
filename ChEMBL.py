import os
import sys
import requests
import pandas as pd
from lxml import etree
import chembl_webresource_client as chembl

targets = chembl.TargetResource()

def get_uniprot_id(pdb_id):
    """Use PDB ID to find Uniprot ID"""

    with requests.get('http://www.rcsb.org/pdb/rest/describeMol?structureId=%s' % pdb_id) as req:
        req_tree = etree.XML(req.text)
        uniprot_id = req_tree.find('.//accession').get('id')

    return uniprot_id


def get_chembl_id(uniprot_id):
    """Use Uniprot ID to find ChEMBL target ID"""

    try:
        chembl_id = targets.get(uniprot=uniprot_id)['chemblId']
        return chembl_id
    except:
        print('No chemblId for %s' % uniprot_id, file=sys.stderr)


def smiles_chembl(chembl_ids):
    """Use compound IDs from ChEMBL to find SMILES

    Parameters
    ----------
    chembl_ids: list of strings
        List of ChEMBL IDs.

    Returns
    -------
    smiles: list of strings
        List of smiles strings.
    """

    smiles = []
    compounds = chembl.CompoundResource()
    for chembl_id in chembl_ids:
        try:
            smiles.append(compounds.get(chembl_id)['smiles'])
        except:
            smiles.append('Unspecified')

    return smiles


def chembl_to_data_frame(uniprot_id):
    """Prepare data frame for target

    Parameters
    ----------
    uniprot_id: str
        Uniprot ID.

    Returns
    -------
    data: DataFrame
        Columns: uniprot_id, chembl_id, smiles, bioactivity_type, operator, value, units.

    """

    chembl_id = get_chembl_id(uniprot_id)
    if not chembl_id:
        return None

    activities = targets.bioactivities(chembl_id=chembl_id)
    chemblid = []
    act_type = []
    operator = []
    value = []
    units = []
    for act in activities:
        if act['value'] != 'Unspecified':
            chemblid.append(act['ingredient_cmpd_chemblid'])
            act_type.append(act['bioactivity_type'])
            operator.append(act['operator'])
            value.append(act['value'])
            units.append(act['units'])

    smiles = smiles_chembl(chemblid)
    data = pd.DataFrame({'chembl_id': chemblid, 'bioactivity_type': act_type,
                         'operator': operator, 'value': value, 'units': units,
                         'smiles': smiles, 'uniprot_id': uniprot_id})

    return data[['uniprot_id', 'chembl_id', 'smiles',
                 'bioactivity_type', 'operator', 'value', 'units']]


def create_data_frames(uniprot_ids, directory='', overwrite=False):
    """Create and save data frames to csv

    Parameters
    ----------
    uniprot_ids: list of stings
        List of strings: uniprot IDs.

    directory: str
        Directory where files will be saved.

    overwrite: bool
        Overwrites your files if True.

    """

    for uniprot_id in uniprot_ids:

        f = os.path.join(directory, 'uniprot', 'chembl_%s.csv' % uniprot_id)
        if overwrite or (not overwrite and not os.path.isfile(f)):

            data = chembl_to_data_frame(uniprot_id)
            if data is not None:
                data.to_csv(f, index=False)


def convert_unit(df, old_unit, new_unit, factor=1):
    """Replace old unit with new unit in bioactives dataframe

    Parameters
    ----------
    df : DataFrame
        DataFrame created by chembl_to_data_frame function.

    old_unit : str
        Unit you want to convert.

    new_unit : str
        Unit to convert to.

    factor : float
        Factor to convert units.

    """

    idxs = df[df['units'] == old_unit].index
    for idx in idxs:
        val = df['value'][idx]
        df.set_value(index=idx, col='value', value=val * factor)
        df.set_value(index=idx, col='units', value=new_unit)


def find_all_units(uniprot_ids, directory=''):
    """Find all units in all dataframes"""

    units = []
    for uniprot_id in uniprot_ids:

        f = os.path.join(directory, 'uniprot', 'chembl_%s.csv' % uniprot_id)
        if os.path.isfile(f):
            data = pd.read_csv(f)
            units += data['units'].tolist()

    return set(units)

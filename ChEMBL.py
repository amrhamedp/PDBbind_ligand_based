import os
import requests
import pandas as pd
from lxml import etree
import chembl_webresource_client as chembl

targets = chembl.TargetResource()

def get_uniprot_id(pdb_id):
    """Use pdb ID to find uniprot ID"""

    with requests.get('http://www.rcsb.org/pdb/rest/describeMol?structureId=%s' % pdb_id) as req:
        req_tree = etree.XML(req.text)
        uniprot_id = req_tree.find('.//accession').get('id')

    return uniprot_id


def get_chembl_id(uniprot_id):
    """Use uniprot ID to find chembl target ID"""

    try:
        chembl_id = targets.get(uniprot=uniprot_id)['chemblId']
        return chembl_id
    except:
        print('No chemblId for', uniprot_id)


def smiles_chembl(chembl_ids):
    """Use compound IDs from ChEMBL to find SMILES"""

    smiles = []
    compounds = chembl.CompoundResource()
    for chembl_id in chembl_ids:
        try:
            smiles.append(compounds.get(chembl_id)['smiles'])
        except:
            smiles.append('Unspecified')

    return smiles


def prepare_data_frame(uniprot_id):
    """Prepare data frame for target"""

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
    """Create and save data frames to csv"""

    for uniprot_id in uniprot_ids:

        f = os.path.join(directory, 'uniprot', 'chembl_%s.csv' % uniprot_id)
        if overwrite or (not overwrite and not os.path.isfile(f)):

            data = prepare_data_frame(uniprot_id)
            if data is not None:
                data.to_csv(f, index=False)


def replace(bs, old_unit, new_unit, factor=1):
    """Replace old unit with new unit in bioactives dataframe"""

    idxs = bs[bs['units'] == old_unit].index
    for idx in idxs:
        val = bs['value'][idx]
        bs.set_value(index=idx, col='value', value=float(val) * factor)
        bs.set_value(index=idx, col='units', value=new_unit)


def find_all_units(uniprot_ids, directory=''):
    """Find all units in all dataframes"""

    units = []
    for uniprot_id in uniprot_ids:

        f = os.path.join(directory, 'uniprot', 'chembl_%s.csv' % uniprot_id)
        if os.path.isfile(f):
            data = pd.read_csv(f)
            units += data['units'].tolist()

    return set(units)


def change_units(bs):

    replace(bs, old_unit='10\'4mol/L', new_unit='10\'4M', factor=1)
    replace(bs, old_unit='l mol-1', new_unit='/M', factor=1)
    # /nM
    replace(bs, old_unit='/M', new_unit='/nM', factor=1 / 10e9)
    replace(bs, old_unit='/uM', new_unit='/nM', factor=1 / 10e3)
    replace(bs, old_unit='10\'2/M', new_unit='/nM', factor=1 / 10e7)
    replace(bs, old_unit='10\'4/M', new_unit='/nM', factor=1 / 10e5)
    replace(bs, old_unit='10\'5/M', new_unit='/nM', factor=1 / 10e4)
    replace(bs, old_unit='10\'7/M', new_unit='/nM', factor=1 / 10e2)
    replace(bs, old_unit='10\'8/M', new_unit='/nM', factor=1 / 10e1)
    replace(bs, old_unit='10\'9/M', new_unit='/nM', factor=1)
    # nM
    replace(bs, old_unit='M', new_unit='nM', factor=10e9)
    replace(bs, old_unit='mM', new_unit='nM', factor=10e6)
    replace(bs, old_unit='uM', new_unit='nM', factor=10e3)
    replace(bs, old_unit='pM', new_unit='nM', factor=1 / 10e3)
    replace(bs, old_unit='10\'5mM', new_unit='nM', factor=10e11)
    replace(bs, old_unit='10\'4M', new_unit='nM', factor=10e13)
    replace(bs, old_unit='10\'10M', new_unit='nM', factor=10e19)

    idxs = bs[bs['units'] == '/nM'].index
    for idx in idxs:
        val = bs['value'][idx]
        bs.set_value(index=idx, col='value', value=1 / float(val))
        bs.set_value(index=idx, col='units', value='nM')

    return bs


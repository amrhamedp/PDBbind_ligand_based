import os
import requests
import pandas as pd
from lxml import etree
from chembl_webresource_client import *


def get_uniprot_id(pdbbind_id):
    """Use pdbbind id to find uniprot id"""

    req = requests.get('http://www.rcsb.org/pdb/rest/describeMol?structureId=%s' % pdbbind_id)
    req_tree = etree.XML(req.text)
    uniprot_id = req_tree.find('.//accession').get('id')
    req.close()

    return uniprot_id


def get_chembl_id(uniprot_id):
    """Use uniprot id to find chembl id"""

    targets = TargetResource()
    chembl_id = targets.get(uniprot=uniprot_id)['chemblId']

    return chembl_id


def bioactivities_chembl(chembl_id):
    """Find all bioactives of target from chembl"""

    targets = TargetResource()
    activities = targets.bioactivities(chembl_id)

    return activities


def smiles_chembl(chembl_ids):
    """Use ids from ChEMBL to find SMILES (simplified molecular-input line-entry system)"""
    smiles = []
    compounds = CompoundResource()
    for chembl_id in chembl_ids:
        try:
            smiles.append(compounds.get(chembl_id)['smiles'])
        except Exception as e:
            print(chembl_id, e)
            smiles.append('Unspecified')

    return smiles


def prepare_data_frame(chembl_id, smiles=False):
    """Prepare data frame with information about activities from ChEMBL"""

    activities = bioactivities_chembl(chembl_id=chembl_id)
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

    if smiles:
        smiles = smiles_chembl(chemblid)
        data = pd.DataFrame({'chembl id': chemblid, 'bioactivity type': act_type,
                             'operator': operator, 'value': value, 'units': units,
                             'smiles': smiles})
        return data[['chembl id', 'smiles', 'bioactivity type', 'operator', 'value', 'units']]
    else:
        data = pd.DataFrame({'chembl id': chemblid, 'bioactivity type': act_type,
                             'operator': operator, 'value': value, 'units': units})
        return data[['chembl id', 'bioactivity type', 'operator', 'value', 'units']]


def save_data_frame(pdbbind_ids, smiles=False, directory=''):
    """Save data frames to csv"""

    for pdbbind_id in pdbbind_ids:
        f = os.path.join(directory, 'chembl_%s.csv' % pdbbind_id)
        if not os.path.isfile(f):
            uniprot_id = get_uniprot_id(pdbbind_id=pdbbind_id)
            chembl_id = get_chembl_id(uniprot_id=uniprot_id)

            data = prepare_data_frame(chembl_id, smiles=smiles)
            data.to_csv(f)

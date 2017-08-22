import requests
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
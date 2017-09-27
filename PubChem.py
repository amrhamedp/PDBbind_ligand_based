import os
import requests
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


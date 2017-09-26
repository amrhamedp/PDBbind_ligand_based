import os
import requests
from lxml import etree


def get_gene_name(uniprot_id):
    """Use Uniprot ID to get gene name"""

    try:
        req_tree = etree.parse('http://www.uniprot.org/uniprot/%s.xml' % uniprot_id)
        gene_name = req_tree.find('.//{http://uniprot.org/uniprot}gene/{http://uniprot.org/uniprot}name').text

        return gene_name

    except:
        print('No gene name for', uniprot_id)
        return None


def get_AIDs(gene_name):
    """Use gene name to get list of PubChem Assay ID"""

    req = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/genesymbol/%s/aids/TXT' % gene_name)
    aids = req.text.split()

    return aids


def download_aid_csv(aid, directory):
    """Download information from PubChem to directory"""

    f = os.path.join(directory, '%s.csv' % aid)
    if not os.path.isfile(f):
        os.system('wget -qO %s https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/CSV' % (f, aid))


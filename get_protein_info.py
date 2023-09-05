# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 15:59:38 2023

@author: mast527
"""

import requests

from dash_bio.utils import protein_reader
import urllib.request as urlreq


from Bio import Entrez

def search_db(database='protein', search="Aplysia AND Pleurin", idtype='acc', **kwargs):
    '''
    input:
        database (str) --database to search.  default to protein  
        Here are the options:
            ['pubmed', 'protein', 'nuccore', 'ipg', 'nucleotide',
             'structure', 'genome', 'annotinfo', 'assembly',
             'bioproject', 'biosample', 'blastdbinfo', 'books',
             'cdd', 'clinvar', 'gap', 'gapplus', 'grasp', 'dbvar',
             'gene', 'gds', 'geoprofiles', 'homologene', 'medgen', 
             'mesh', 'nlmcatalog', 'omim', 'orgtrack', 'pmc', 'popset',
             'proteinclusters', 'pcassay', 'protfam', 'pccompound',
             'pcsubstance', 'seqannot', 'snp', 'sra', 'taxonomy',
             'biocollections', 'gtr'] 
            
        search (str) -- the term to search database 
    kwargs:
        retstart -- Sequential index of the first UID in the retrieved set 
                    to be shown in the XML output (default=0, corresponding
                    to the first record of the entire set). This parameter can
                    be used in conjunction with retmax to download an arbitrary
                    subset of UIDs retrieved from a search.
        retmax -- Total number of UIDs from the retrieved set to be shown in the XML output
        
        idtype --Specifies the type of identifier to return for sequence
                 databases (nuccore, popset, protein). By default, ESearch 
                 returns GI numbers in its output. id type to ‘acc’
        
        usehistory
        WebEnv
        query_key
        
        example: record = search_protein_db2(search='325296977', idtype='acc')
        '''
    Entrez.email = "david.mast@pnnl.gov"
    handle = Entrez.esearch(db=database, term=search, **kwargs)
    record=Entrez.read(handle)
    handle.close()
    print(record)
    return record



def get_record(id=''):
    '''kwargs:
        db, id, retype, retmode'''
    Entrez.email = "david.mast@pnnl.gov"
    try:
        #handle = Entrez.efetch(**kwargs)

        handle = Entrez.efetch(db="protein", id=id, s='8', rettype="gt", retmode="text",)
        read = handle.read()
        handle.close()
    except Exception as err:
        print()
        print(err)
        print('something went wrong!')
        
    else:
        print(read)
    
    
        return read
 



def get_uniprot_data(protein_id):
    base_url = "https://www.uniprot.org/uniprot/"
    response_format = ".json"  # Retrieve data in JSON format

    url = base_url + protein_id + response_format

    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception if there was an error
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error: {e}")
        print('Invalid protein ID. Please enter a valid uniprot id')
        raise requests.exceptions.RequestException

def get_uniprot_fasta(protein_id):
    
    url2 = f"https://www.uniprot.org/uniprot/{protein_id}.fasta"
    fasta_str = urlreq.urlopen(url2).read().decode('utf-8')
    sequence = protein_reader.read_fasta(datapath_or_datastring=fasta_str, is_datafile=False)[0]['sequence']

def get_protein_sequence(protein_id:str):
    protein_sequence = get_uniprot_data(protein_id)['sequence']['value']
    return protein_sequence
if __name__ == '__main__':
    x = get_protein_sequence('P10636')
    x1 = get_record('NP_001191504.1')
# # Example usage
# protein_id = "P12345"  # Replace with the UniProt ID of the protein you want to retrieve
# protein_data = get_uniprot_data(protein_id)

# if protein_data:
#     print(protein_data)
   
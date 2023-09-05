# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 10:22:20 2023

@author: mast527
"""

import gzip
from pyteomics import fasta
import os


def get_sequence_from_fasta(uniprot_accession:str, FASTA_gz_file:str):
    '''
    retrieves a protein sequence from a fasta file
    Input:
        uniprot_accession -- accesion for uniprot (ex:"A0PK00")
        FASTA_gz_file -- filepath to database (ex: human_proteome-2023.01.26-18.32.50.86.fasta.gz
    
    Output:
        tuple((sequence description, sequence))'''
    c = {"white": '\x1b[1;37m ',
         "red": '\x1b[1;31m',
         "yellow": '\x1b[1;33m',
         "blue": '\x1b[1;34m',
         "green": '\x1b[1;32m'}
    
    seq_descriptions = None
    with gzip.open(FASTA_gz_file, mode='rt') as gzfile:
        for description, sequence in fasta.FASTA(gzfile):
            if uniprot_accession in description:
                print(f"{c['blue']} {description} {c['yellow']} {uniprot_accession}{c['white']}.")
                seq_descriptions = tuple((description, {uniprot_accession:sequence}))
                break
            #seq_descriptions.append(description)
            else:
                pass
    
    if seq_descriptions is not None:
        print(f'Done! {c["yellow"]} {seq_descriptions[0]} {c["white"]} sequence obtained!')
        return seq_descriptions
    else: 
        print(f"{c['red']} Sequence aquisition FAILED:{c['white']} accession '{uniprot_accession}' not in FASTA")
        return None


def get_sequence_from_list(accession_list:list, FASTA_gz_file:str):
    '''obtain a list of sequence accessions of a fasta file
    
    Input:
        accession_list (list) -- a list of uniprot accessions as strings
                                    (ex: ["P10636", "A0PK00", "P10636-1"])
        FASTA_gz_file (str) -- filepath to the FASTA.gz file (ex: "./human_proteome-2023.01.26-18.32.50.86.fasta.gz" )
                
        
    
    '''
    
    f = get_sequence_from_fasta
    
    c = {"white": '\x1b[1;37m ',
         "red": '\x1b[1;31m',
         "yellow": '\x1b[1;33m',
         "blue": '\x1b[1;34m',
         "green": '\x1b[1;32m'}
    
    if os.path.isfile(FASTA_gz_file) and isinstance(accession_list,list):
            
        print("True, this is a file")
        return [f(peptide,FASTA_gz_file) for peptide in accession_list]
    
    else: 
        if os.path.isfile(FASTA_gz_file) is False:
            print(f"{c['red']} {FASTA_gz_file} not in current working directory or not a valid file c['white]")
            
    
        if isinstance(accession_list, list) is False:
            print("accession_list must be type: list")

        return None


    
# if __name__ == '__main__':
#     FASTA = r"human_proteome_UP000005640_9606.fasta.gz"
#     FASTA_dir = r'.\protein_databases'
#     FASTA_path = os.path.join(FASTA_dir, FASTA)
#     tau = get_sequence_from_fasta(uniprot_accession='P10636', FASTA_gz_file=FASTA_path)
    
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 08:09:05 2023

@author: mast527
"""

import get_protein_sequence_from_fasta as gpsff
import re
import os
import string_permutations



FASTA = r"human_proteome_UP000005640_9606.fasta.gz"
FASTA_dir = r'.\protein_databases'
FASTA_path = os.path.join(FASTA_dir, FASTA)
tau = gpsff.get_sequence_from_fasta(uniprot_accession='P10636', FASTA_gz_file=FASTA_path)
tau = tau[1]['P10636']





def find_pattern_in_prot(protein, s):
    
    l = string_permutations.get_randomized_strings(s)
    l1 = []
    r = []
    for i in l:
        x = re.compile(i)
        x1 = x.findall(protein)
        if len(x1) >=1:
            l1.append(i)
            r.append(x1)
        else:
            pass
    return dict(zip(l1,r))
    
def split_protein_sequence(protein_sequence):
    """
    Split a protein sequence into all possible sequence lengths.
    
    Args:
        protein_sequence (str): The input protein sequence.

    Returns:
        list: A list of all possible subsequences of varying lengths.
    """
    subsequences = []
    sequence_length = len(protein_sequence)

    for start in range(sequence_length):
        for end in range(start + 1, sequence_length + 1):
            subsequence = protein_sequence[start:end]
            subsequences.append(subsequence)

    return subsequences


protein_sequence = "MAAGGTRAFAGGPEAGGVGGR"
subsequences = split_protein_sequence(protein_sequence)

for subsequence in subsequences:
    print(subsequence)   
    
    
def max_peptide_length():
    pass




def get_centered_peptides(expression: str, peptide_list: list, all_possible=False, return_original_seqs=False):
    '''
    Generates aligned peptides based on a regex pattern
    Input: 
        expression (str) -- regular expression to align peptides
        peptide_list (list) -- a list of strings containing peptide sequences
    
    Output:
        a list containing aligned peptide sequences with string legnth of len(regex)
    '''
    if isinstance(expression, str) is False:
        raise TypeError("expression must be a string")

    elif isinstance(peptide_list, list) == False:
        raise TypeError("expression must be a list")

    pat = re.compile(expression)
    l_new = []
    l_new_seq = []
    l_pat = []

    for peptide in peptide_list:

        m = pat.findall(peptide)

        if len(m) > 0:
            if all_possible:
                for _ in m:
                    l_new.append(_)
                    l_new_seq.append(peptide)
                    l_pat.append(pat.pattern)
            else:
            
                m = m[0]
                l_new.append(m)
                l_pat.append(pat.pattern)

        else:
            pass
    
    
    if return_original_seqs:
        return l_new, l_new_seq
    
    return l_new
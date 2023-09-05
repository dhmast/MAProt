# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 09:16:34 2023

@author: mast527
"""

import re
from dataclasses import dataclass





@dataclass
class CenteredPeptideParams:
    '''
    regex_list : 
        - A list of regular expressions to find in protein
    min_seq_length : 
        -the minimum sequence length requested
    '''
    regex_list: list
    min_seq_length: int
    



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


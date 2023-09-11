# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 09:16:34 2023

@author: mast527
"""

import re
from dataclasses import dataclass
from collections import OrderedDict
import string_permutations
import random

def generate_random_peptide_list(num_peptides, min_length=6, max_length=15):
    
    '''
    
    Generate a list of random peptide sequences with specified characteristics.

    Args:
        num_peptides (int): The number of random peptides to generate.
        min_length (int, optional): The minimum length of a peptide sequence (default is 6).
        max_length (int, optional): The maximum length of a peptide sequence (default is 15).

    Returns:
        list: A list of randomly generated peptide sequences.

    A peptide sequence is composed of amino acids selected from a predefined set of
    amino acids. This function generates 'num_peptides' random peptide sequences,
    each with a length between 'min_length' and 'max_length' (inclusive).

    Example:
    >>> generate_random_peptide_list(10)
    ['PKMENL', 'YFIVLNT', 'ATWHTESM', 'RDYQLMST', 'VLEKFGHYR', 'MVESWI', 'VKSFIYD',
    'QMVFWIGY', 'RVTLMYCQI', 'SNMWVI']
    
    '''
    # Define the amino acids you want to include
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    peptide_list = []
    for _ in range(num_peptides):
        peptide_length = random.randint(min_length, max_length)
        peptide = ''.join(random.choice(amino_acids) for _ in range(peptide_length))
        peptide_list.append(peptide)
    return peptide_list


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

def get_centered_peptides(expression: str, peptide_list: list, all_possible=False, return_counts=False, return_original_seqs=False):
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
    
    
    if return_counts:
        all_possible = True
    
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
    
    if return_counts:
        # Create an ordered dictionary to store element counts
        element_counts = OrderedDict()
        
        # Use a list comprehension to count occurrences and populate the ordered dictionary
        [element_counts.update({element: element_counts.get(element, 0) + 1}) for element in l_new]
        
        # Convert the ordered dictionary to a list of tuples
        #ordered_counts_list = list(element_counts.items())
        
        sequences_list = list(element_counts.keys())
        counts_list = [_ for i, _ in element_counts.items()]
        return sequences_list, counts_list
        

    if return_original_seqs:
        return l_new, l_new_seq
    
    return l_new


def check_regex(regex:str):
    
    pass

# class CenteredPeptides:
    
#     def __init__(self, sequence_list:list, regex:str):
        
#         self.sequence_list = sequence_list
#         self.regex = regex
#         self.centered_peptides = []
        
#     def get_centered_peptide_combos(self):
#         self.ceterered_peptides = get_centered_peptides(self.regex, self.sequence_list)
        
    



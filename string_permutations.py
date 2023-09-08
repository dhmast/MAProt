# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 09:16:34 2023

@author: mast527
This module contains  functions for generating lists of string permutations. 
The get_randomized_strings() function takes in a single string, and then randomizes it until
all possible combinations are identified. 

It is much faster than using itertools to find permutations.
Example:
    In: string = 'ABC'
    In: x = get_randomized_strings(string)
    In: print(x)
    Out: ['BAC', 'CBA', 'ACB', 'BCA', 'CAB', 'ABC'] 
    
    
    
"""

import re
import random
import math
from collections import Counter
from itertools import permutations
import matplotlib.pyplot as plt
def num_perms_wo_repeats(regular_expression: str, other: bool = False) -> int:
    """
    Computes the number of unique permutations of a string with repeated characters.

    This function takes a string as input and computes the number of unique
    permutations that can be formed by rearranging its characters, even if
    some of the characters are repeated. If the `other` flag is set to True,
    the function also returns a tuple containing a list of the unique
    characters in the input string and their corresponding counts.

    Args:
        regular_expression (str): The input string to permute.
        other (bool): Optional flag indicating whether to return additional
                      information about the input string. Default is False.

    Returns:
        int or Tuple[int, List[str], List[int]]: The number of unique permutations
        of the input string. If `other` is True, also returns a tuple containing
        a list of the unique characters in the input string and their corresponding
        counts.
    """
    if not isinstance(regular_expression, str):
        raise TypeError("Input must be a string.")
    
    counts = Counter(regular_expression)
    factorial_counts = [math.factorial(count) for count in counts.values()]
    unique_combos = math.factorial(len(regular_expression)) // math.prod(factorial_counts)
    
    if other:
        chars = list(counts.keys())
        counts = list(counts.values())
        return unique_combos, chars, counts
    else:
        return unique_combos


def num_perms_w_repeats(s):
    
    '''
    computes the number of permutations with repeats.
    
    Parameters: 
        s : (str)
            - a string to 
    '''
    n = len(s)
    r = len(list(set(s)))
    if n < r:
        return 0  # There are no permutations if r > n
    else:
        return math.factorial(n) // math.factorial(n - r)

def get_string_permutations_itertools(s:str):
    '''
    itertools implementation
    s : (str)
        -string to get combos of. ex "....KR..."
    
    
    '''
    num_perms = num_perms_w_repeats(s)
    if num_perms >= 1000:
        print('num permutationss too high!!')
        return
    nums = list(s)
    perms = list(permutations(nums))
    perms = list(set([''.join(perm) for perm in perms]))
 
   
    return perms


def randomize_string(input_string):
    '''
    Parameters: 
        s : (str)
            - a string to shuffle
    '''    
    # Convert the string to a list of characters for shuffling
    char_list = list(input_string)

    # Shuffle the list of characters randomly
    random.shuffle(char_list)

    # Convert the shuffled list back to a string
    randomized_string = ''.join(char_list)
    
    return randomized_string

def randomize_string_by_number(input_string):
    char_list = list(input_string)
    n = len(char_list)

    # Create a shuffled index list
    shuffled_indices = list(range(n))
    random.shuffle(shuffled_indices)

    # Reorder characters based on shuffled indices
    randomized_string = ''.join(char_list[i] for i in shuffled_indices)

    return randomized_string



def randomize_string_no_repeats(input_string):
    '''slowest option'''
    char_list = list(input_string)
    n = len(char_list)
    
    while True:
        shuffled_chars = random.sample(char_list, n)
        randomized_string = ''.join(shuffled_chars)
        yield randomized_string



def get_randomized_strings(s, by=None, create_plot=False):
    '''
    

    Parameters
    ----------
    s : string
        regex pattern which will be randmized.
    by : str, optional
        which function will be used to randomized the string. The default is None.
        if by='number' than the randomize_string_by_number() function will be used
    create_plot : TYPE, optional
        DESCRIPTION. The default is False.

    Raises
    ------
    ValueError
        value error raised if the number of permutations is too large. 

    Returns
    -------
    hits : list
        a list of all possible permutations without repeated orders of the input string.
    
    Example
    -------
        In: string = 'ABC'
        In: x = get_randomized_strings(string)
        In: print(x)
        Out: ['BAC', 'CBA', 'ACB', 'BCA', 'CAB', 'ABC'] 
    '''
    total_perms = num_perms_wo_repeats(s)
    if total_perms > 400_000:
        print('WAY TOO MANY PERMUTATIONS! try a different string!')
        raise ValueError(f'{total_perms}...This number of permutations is too high!')
        
    print('expected number of permutations without repeats: ', total_perms)
    hits = set()
    num_hits = []
    num_iters = []
    h_to_i = []
    running = True
    i=0
    
    if by =='number':
        function = randomize_string_by_number
    
    else:
        function = randomize_string
    
    # randomize string by the
    while running:
        i+=1
        hit = function(s)
        hits.add(hit)
       
        num_hits.append(len(hits)/total_perms)
        num_iters.append(i)
        h_to_i.append(len(hits)/i)
        
        if len(hits) == total_perms:
            print(f'total number of permutations reached! at {i} iterations')
            running = False
            break
        elif len(hits)/i <= 0.05:
            print('computation time exceeded maximum')
            running = False
            break
            
    hits = list(hits)
    if create_plot:
        plt.scatter(num_iters, num_hits, label='number of hits/total permutations')
        plt.scatter(num_iters, h_to_i, label='number of hits / number iterations')
        plt.legend()
        plt.xlabel('# iterations')
        
    print('calculated_number of combos: ', len(hits))
    print('hits-to-iterations ratio : ', f'{len(hits)/i:0.2f}')
    return hits



def find_pattern_in_prot(protein, s):
    '''
    Parameters
    ----------
    protein : str
        protein sequence to search patterns
    s : str
        the regex string to generate permuations of, and search in sequence. ex: ...KR...

    Returns
    -------
    dict
        A dictionary of all regex permutations and the matched sequences in the protein.

    '''
    l = get_randomized_strings(s)
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


 

# def get_randomized_strings_no_repeats(s):
#     hits = set()
#     num_hits = []
#     num_iters = []
#     running = True
#     i=0
#     total_shuffles = 0
#     shuffler = randomize_string_no_repeats(s)

#     while running:
#         i += 1
#         shuffled = next(shuffler)
        
#         if shuffled not in hits:
#             hits.add(shuffled)
#             total_shuffles += 1
#         num_hits.append(len(hits))
#         num_iters.append(i)
#         if total_shuffles >= num_perms_wo_repeats(s):
#             running = False
#             break 
#     plt.scatter(num_iters, num_hits)
#     return hits

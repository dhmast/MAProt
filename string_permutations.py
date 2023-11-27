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
from itertools import permutations, product
import matplotlib.pyplot as plt

def sub_leu_ile_permutations(s: str):
    '''
    Generate all permutations of peptide sequences by swapping 'I' with 'L' and 'L' with 'I'.

    Parameters:
        s : str
            The input peptide sequence containing 'I' and 'L'.

    Returns:
        generator
            Yields strings containing all possible permutations of 'I' and 'L' in the peptide.

    Example:
        original_string = 'KTKIELDF'

        ## Generate all combinations
        >>> combinations = set(sub_leu_ile_permutations(original_string))

        >>> for combo in combinations:
        ...     print(combo)
        
        ## Output
        KTKIELDF
        KTKLELDF
        KTKIELDF
        KTKLELDF
    '''
    replacements = {'I': ['I', 'L'], 'L': ['I', 'L']}
    possible_combinations = [replacements.get(c, [c]) for c in s]
    for combo in product(*possible_combinations):
        yield ''.join(combo)


def peptide_substitution_permutations(s:str, replacements:dict):
    '''
    substitutes characters in a string based on a dict
    
    Parameters:
        s : str
            The input peptide sequence.
        
        replacements : dict
            The mapping for the replacements. 
            Ex 1.  {'I': ['I', 'L'], 'L': ['I', 'L']}
            Ex 2. {'A':['X'], 'F':['U']} #will replace all 'A' with 'X' and all 'F' with 'U'
        
    Returns:
        generator
            Yields strings containing all possible permutations of the mapping in the peptide.

    Example:
        >>> original_string = 'KTKIELDF'
        >>> replacement_dict = {'I': ['I', 'L'], 'L': ['I', 'L']}
        
        ## Generate all combinations
        >>> combinations = set(sub_leu_ile_permutations(original_string))

        >>> for combo in combinations:
        ...     print(combo)
        
        ## Output
        KTKIELDF
        KTKLELDF
        KTKIELDF
        KTKLELDF

    '''
    
    possible_combinations = [replacements.get(c, [c]) for c in s]
    for combo in product(*possible_combinations):
        yield ''.join(combo)

def num_perms_wo_repeats2(regular_expression: str, r, other: bool = False) -> int:

    if not isinstance(regular_expression, str):
            raise TypeError("Input must be a string.")
        
    counts = Counter(regular_expression)
    #factorial_counts = [math.factorial(count) for count in counts.values()]
    #unique_combos = math.factorial(len(regular_expression)) // math.prod(factorial_counts)
    unique_combos = math.factorial(len(regular_expression)) // math.factorial(len(regular_expression)- r)
    if other:
        chars = list(counts.keys())
        counts = list(counts.values())
        return unique_combos, chars, counts
    else:
        return unique_combos

    


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




def replace_regex_chars(input_string):
    '''
    replaces bracketed amino acids in a regex with a number, for shuffling
    Parameters:
        input_string : str
            the input regex string.
    
    returns:
        output_string : str
            The string with replaced characters
        substitution_mapping : dict
            the dictionary for reversing the substitution
    example:
        

        In: replace_regex_chars("....[AB]..[C]...")
        
        Out: ("....1..2...", {'0': '[AC]', '1': '[B]'}
    '''
    # remove start of line and end of line characters
    if input_string[len(input_string)-1] == '$':
        input_string = input_string.rstrip('$')
    if input_string[0] == '^':
        input_string = input_string.lstrip('^')
    
    
    replacement_regex = r'\[\w*\]'
    replacement = re.compile(replacement_regex)
    chars_to_replace = list(replacement.findall(input_string))
    
    # Define a list of unique replacement values
    replacement_values = list(range(len(chars_to_replace)))  # Add more values if needed
    replacement_values = [str(i) for i in replacement_values]
    substitution_mapping = dict(zip(replacement_values, chars_to_replace))
    
    # Use re.sub() with the callback function to perform the substitution
    output_string = re.sub(r'\[\w*\]',lambda x: replacement_values.pop(0) , input_string)
    return output_string, substitution_mapping


def reverse_replace_regex_chars(output_string, substitution_mapping):
    '''
    reverses the replace_regex_chars() function. Returns the original regex string
    '''
    # Use re.sub() with the callback function to perform the reverse substitution
    reversed_string = re.sub(r'\d',lambda x: substitution_mapping.get(x.group()), output_string)
    return reversed_string

def check_regex_input(input_string):
    invalid_chars = ',!@#%&*))/ -+'
    checks = dict(
        has_brackets = False,
        has_start = False,
        has_end = False,
    )
    
    bracket_regex = r'\[\w*\]'
    m = re.compile(bracket_regex)
    
    if m.search(input_string):
        checks.update(dict(has_brackets=True))
        bracketed = m.findall(input_string)
        print('number of substitutions to be made = ', len(m.findall(input_string)))
        
    if input_string[0] == '^':
        checks.update(dict(has_start=True))
    
    if input_string[len(input_string)-1] == '$':
        checks.update(dict(has_end=True))
    
    return checks
    


def convert_regex_input(input_string):
    
    
    #check for bracketed words
    bracket_regex = r'\[\w*\]'
    m = re.compile(bracket_regex)
    if m.search(input_string):
        
        print('Input string contains match')
        try:
            output_string, output_mapping = replace_regex_chars(input_string)
        except TypeError as err:
            print(err)
            raise TypeError('replace_regex_chars() function encountered an error. Inspect the input_string')
        return output_string, output_mapping
    else:
        return input_string, False
    


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
    checks = check_regex_input(s)
    s, output_mapping = convert_regex_input(s)
    print(s)
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
    if output_mapping:
        hits = [reverse_replace_regex_chars(_, output_mapping) for _ in hits]
        
    if checks.get('has_start'):
        hits = ['^'+hit for hit in hits]
    if checks.get('has_end'):
        hits = [hit+'$' for hit in hits]
         
    return hits



def get_randomized_strings_length_n(s, r=5, by=None, create_plot=False, stop_after=None):
    
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
    checks = check_regex_input(s)
    s, output_mapping = convert_regex_input(s)
    print(s)
    total_perms = num_perms_wo_repeats2(s, r)
    if stop_after is None:

        if total_perms > 400_000_000_000:
            print('WAY TOO MANY PERMUTATIONS! try a different string!')
            raise ValueError(f'{total_perms}...This number of permutations is too high!')
    else:
        if isinstance(stop_after, int):
            print(f'You stopped after finding {stop_after} permutations')
            total_perms = stop_after
        else:
            pass
            
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
        hits.add(hit[:r])
       
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
    if output_mapping:
        hits = [reverse_replace_regex_chars(_, output_mapping) for _ in hits]
        
    if checks.get('has_start'):
        hits = ['^'+hit for hit in hits]
    if checks.get('has_end'):
        hits = [hit+'$' for hit in hits]
         
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

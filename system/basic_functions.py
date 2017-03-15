# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 00:14:13 2017

@author: Demetris
"""

def join_list(gene_list):
    'gets a gene list and returns a string'
    return ''.join(gene_list)
    
def unpack_fasta(fasta_string):
    'gets a fasta string and returns its id and the dna sequence'
    data = list(fasta_string)
    
    #find where the genome is located in the string
    search = len('>Rosalind_')
    not_found = True
    while not_found:
        if data[search].isdigit():
            search +=1
        else:
            not_found = False
            
    extract_id = data[1:search]
    extract_dna = data[search:]
    
    return join_list(extract_id), join_list(extract_dna)  
    
def convert_to_bases(value):
    'gets a numpy array with numbers from zero to 3 and converts them to bases'
    
    if value == 0:
        return 'A'
    elif value == 1:
        return 'C'
    elif value == 2:
        return 'G'
    elif value == 3:
        return 'T'
    else:
        raise ValueError('the parameter numpy array should only contain 0,1,2,3')
    
def convert_base_to_number(base):
    'gets a base and converts it to numbers from 0 to 3'
    if base == 'A':
        return 0
    elif base == 'C':
        return 1
    elif base == 'G':
        return 2
    elif base == 'T':
        return 3
    else:
        raise ValueError('invalid base')
        
#print(convert_to_bases(3))
#print(convert_base_to_number('T'))
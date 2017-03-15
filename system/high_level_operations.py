# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 22:07:53 2017

@author: Demetris
"""
import numpy as np
from basic_functions import join_list, unpack_fasta, convert_to_bases, convert_base_to_number

def count_bases(dna_list):
    'gets a dna string and counts the number of occureneces of each base'
    dna_array = np.zeros(shape = 4)
    for base in dna_list:
        transcribed_base = convert_base_to_number(base)
        dna_array[transcribed_base]+=1
        
    return dna_array
    
dna = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'      
#print(count_bases(dna))


def transcription(t):
    ''' gets a dna string and returns the rna list'''
    
    data = list(t)
    
    i = 0
    for nucl in data:
        if nucl == 'T':
            data[i] = 'U'
        i +=1
    
    return data
    
a = 'GATGGAACTTGACTACGTAAATT'
#print(join_list(transcription(a)))

def get_reverse_complement(t):
    '''gets a dna string and returns the reverse of the other part of the 
       dna chain that correspnds to it
    '''
    data = list(t)
    data.reverse()
    result = []    
    
    for nucl in data:
        if nucl == 'A':
            result.append('T')
        elif nucl == 'C':
            result.append('G')     
        elif nucl == 'G':
            result.append('C')
        elif nucl == 'T':
            result.append('A')
        else:
            raise ValueError('invalid gene')
        
    return result

#print(unpack_fasta(c[0])) 
#print(unpack_fasta(f[0]))

b = 'AAAACCCGGT'
#print(join_list(get_reverse_complement(b)))

def identify_dna(string_list):
    'gets a list of strings that represent dna in FASTA format and returns the one with higher gc content'
    max_content = ('nothing',0)
    
    for dna_string in string_list:
        counter = 0.0
        extract = unpack_fasta(dna_string)
        extract_dna = list(extract[1])        
        
        #count bases
        bases_num = count_bases(extract_dna)
        counter = bases_num[1] + bases_num[2] + 0.0 #the 0.0 is to make counter a float
        
        gc_content = float(100 * counter /len(extract_dna))
        
        if gc_content > max_content[1]:
            max_content = (extract[0], gc_content)
            
    return max_content

c = ['>Rosalind_6404CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG',
     '>Rosalind_5959CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC',
     '>Rosalind_0808CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT']   
#print(identify_dna(c))

def get_pattern_location(s,t):
    'gets two dna strings and returns the locations of t as a substring of s'
    a = len(t)
    answer = []   
    
    #search through all possible combinations of s to see if they contain t
    for i in range(len(s)-a +1):
        if s[i:i+a] == t:
            answer.append(i+1)
            
    return answer
    
d = 'GATATATGCATATACTT'
e = 'ATAT'
#print(get_pattern_location(d,e))

def get_profile_matrix(string_list):
    'gets a list of FASTA strings and returns a numpy array of their profile matrix'
    extract = []
    length = len(unpack_fasta(string_list[0])[1])
    dna_array = np.zeros(shape = (4,length))
    
    #extract first the useful information out of the fasta format
    for dna_string in string_list:
        extract_id, extract_dna = unpack_fasta(dna_string)
        extract.append((tuple(extract_dna)))
    
    #count the number of occureneces of each base    
    for dna_tuple in extract:
        counter = 0
        for base in dna_tuple:
            transcribed_base = convert_base_to_number(base)
            dna_array[transcribed_base][counter]+=1
            counter +=1
    
    return dna_array
    
def get_consensus_string(profile_matrix):
    'gets a profile matrix numpy array and returns the consensus string'
    #find positions of maximum numbers in the profile_matrix 
    max_array = np.argmax(profile_matrix, axis = 0)
    
    result = []

    for x, value in np.ndenumerate(max_array):
        number = convert_to_bases(value)
        result.append(number)
        
    return join_list(result)

f = ['>Rosalind_1ATCCAGCT', '>Rosalind_2GGGCAACT','>Rosalind_3ATGGATCT','>Rosalind_4AAGCAACC',
     '>Rosalind_5TTGGAACT','>Rosalind_6ATGCCATT','>Rosalind_7ATGGCACT']
#print(get_profile_matrix(f))
#print(get_consensus_string(get_profile_matrix(f)))
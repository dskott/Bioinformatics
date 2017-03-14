# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 22:07:53 2017

@author: Demetris
"""
import numpy as np

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
    
#print(unpack_fasta(c[0])) 
#print(unpack_fasta(f[0]))

b = 'AAAACCCGGT'
#print(join_list(get_reverse_complement(b)))

def identify_dna(string_list):
    'gets a list of strings that represent dna in FASTA format'
    max_content = ('nothing',0)
    
    for dna_string in string_list:
        counter = 0.0
        extract = unpack_fasta(dna_string)
        extract_dna = list(extract[1])        
        
        for nucl in extract_dna:
            if nucl == 'G':
                counter +=1
            elif nucl == 'C':
                counter +=1
        
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

def count_bases(dna_list):
    'gets a list of lists that contain the bases of' 

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
         
f = ['>Rosalind_1ATCCAGCT', '>Rosalind_2GGGCAACT','>Rosalind_3ATGGATCT','>Rosalind_4AAGCAACC',
     '>Rosalind_5TTGGAACT','>Rosalind_6ATGCCATT','>Rosalind_7ATGGCACT']
#print(get_profile_matrix(f))
print(get_consensus_string(get_profile_matrix(f)))
#print(convert_to_bases(np.array([0, 3, 2, 1, 0, 0, 1, 3])))
#print(convert_base_to_number('T'))
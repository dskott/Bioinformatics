# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 22:07:53 2017

@author: Demetris
"""

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
    extract_id = data[1:14]
    extract_dna = data[14:]
    
    return join_list(extract_id), join_list(extract_dna)
    
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
#print(unpack_fasta(s)) 
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
print(get_pattern_location(d,e))
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
dna = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
data = list(dna)

result = [0]*4

for nucl in data:
    if nucl == 'A':
        result[0]+=1
    elif nucl == 'C':
        result[1]+=1     
    elif nucl == 'G':
        result[2]+=1
    elif nucl == 'T':
        result[3]+=1 
    else:
        raise ValueError('invalid gene')
            
print(result)

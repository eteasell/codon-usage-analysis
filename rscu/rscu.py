'''
The purpose of this file iis to calculate RSCU values for specific stretches of genes. This is not necessarily full genes
or just windows.
'''

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from lib.aminoacids import AA_TO_CODONS, CODON_TO_AA, STOP_CODONS
from collections import Counter

class RSCU():
    def __init__(self, sequence: str):
        self.sequence = sequence
        self.rscu_values = self.calculate_rscu()

    def calculate_rscu(self):
        '''
        This function calculates the RSCU of each codon within the given sequence.
        If the length of the given sequence is not a multiple of three, the function will return an empty dictionary.
        This function does NOT assume that the sequence starrts with ATG, since it can also calculate RSCU for windows of genes.
        
        RSCU_ij = (x_ij) / ((1/n_i) * sum over j = 1 to n_i of x_ij)
        '''
        sequence = self.sequence
        seq_len = len(sequence)
        
        if seq_len % 3 != 0:
            print("Invalid sequence: length is not a multiple of three")
            return {}
        
        rscu_vals = {codon: 0 for codon in CODON_TO_AA.keys()}
        
        codons = [sequence[i:i+3] for i in range(0, seq_len, 3)]
        
        # Checking that all codons in sequence are real codons
        if not all((codon in CODON_TO_AA.keys()) or (codon in STOP_CODONS) for codon in codons):
            print("Invalid sequence: contains non-standard codons")
            return {}
        
        # Dictionary containing the number of occurances for each codon in sequence
        counts = Counter(codons)
        
        # iterate through all possible amino acids (regardless of presence in sequence)
        for aa in AA_TO_CODONS.keys():
            # codons_i is the list of j codons for the given amino acid
            codons_i = AA_TO_CODONS[aa]
            
            # n_i is the number of codons for this amino acid i
            n_i = len(codons_i)
            
            # x_i is a numpy array to hold the counts x values for each codon j
            x_i = np.zeros(n_i)
            
            # initialize sum over n_i codons for this ith amino acid
            sum_i = 0
            
            # for each of the j codons, add count value to sum and keep count in x_i
            for j, codon in enumerate(codons_i):
                count_val = counts[codon]
                x_i[j] = count_val
                sum_i += count_val
                
            # normalize sum of x values by n_i  
            normalized_sum = sum_i / n_i
            
            # divide each x_ij by the normalized sum (if non-zero), otherwise all values are zero anyways
            if normalized_sum != 0:
                rscu = x_i / normalized_sum
            else:
                rscu = x_i
                
            # for each codon for this amino acid, add rscu to dict
            for j, codon in enumerate(codons_i):
                rscu_vals[codon] = rscu[j]
                    
        return rscu_vals
    
    
    def pretty_print(self):
        for codon in self.rscu_values.keys():
            print(f"{codon}: {str(self.rscu_values[codon])}")
    
  
    
def extract_window(record: SeqRecord, w: int, d: int, cs) -> str | None:
    '''
    record: the BioSeq sequence to be separated into window and full sequence
    w: window size (must be a multiple of 3)
    d: distance downstream from cleavage site
    cs: cleavage site
    '''
    assert w % 3 == 0
    assert isinstance(d, int)
    
    sequence = str(record.seq)
    if len(record.seq) % 3 != 0:
        print(f"Sequence {record.id} length is not a multiple of three.")
        return None

    if cs+d+w > len(sequence):
        print(f"Sequence {record.id} too short for window choice.")
        return None
    
    return sequence[cs+d:cs+d+w]
        

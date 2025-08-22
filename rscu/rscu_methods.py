'''
The purpose of this file iis to calculate RSCU values for specific stretches of genes. This is not necessarily full genes
or just windows.
'''

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from lib.aminoacids import AA_TO_CODONS_MULTI_CODON_FAMILIES, CODON_TO_AA, CODON_TO_AA_MULTI_CODON_FAMILIES
from collections import Counter

def calculate_rscu(codons: list[str]) -> dict[str, float]:
    '''
    codons: a list of codons in the sequence
    
    This function calculates the RSCU of each codon within the given sequence.
    ASSUMES THAT THE SEQUENCE VALIDITY HAS BEEN CHECKED BY ensure_valid_orf.
    
    RSCU_ij = (x_ij) / ((1/n_i) * sum over j = 1 to n_i of x_ij)
    '''
    assert isinstance(codons, list)
    
    rscu_vals = {codon: 0 for codon in CODON_TO_AA_MULTI_CODON_FAMILIES.keys()}
    
    # Dictionary containing the number of occurances for each codon in sequence
    counts = Counter(codons)
    
    # iterate through all possible amino acids (regardless of presence in sequence)
    for aa in AA_TO_CODONS_MULTI_CODON_FAMILIES.keys():
        # codons_i is the list of j codons for the given amino acid
        codons_i = AA_TO_CODONS_MULTI_CODON_FAMILIES[aa]
        
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

def rscu_to_string(rscu: dict[str, float]) -> str:
    '''
    rscu: a dictionary of codon to RSCU value
    
    This function returns a string representation of the RSCU values in a readable format.
    '''
    output = []
    for codon, value in rscu.items():
        output.append(f"{codon}: {value:.2f}")
    return "\n".join(output)


    
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

def extract_window_and_complement(record: SeqRecord, w: int, d: int, cs) -> str | None:
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
    
    return sequence[cs+d:cs+d+w], sequence[0:cs+d] + sequence[cs+d+w:]


def extract_SP_and_complement(record: SeqRecord, cs) -> str | None:
    '''
    record: the BioSeq sequence to be separated into window and full sequence
    cs: cleavage site
    '''
    
    sequence = str(record.seq)
    
    if len(record.seq) % 3 != 0:
        print(f"Sequence {record.id} length is not a multiple of three.")
        return None
    
    return sequence[0:cs], sequence[cs:]


def ensure_valid_orf(sequence: str | None) -> tuple[bool, list[str]]: 
    if sequence is None: return False, []
    
    seq_len = len(sequence)
    
    if seq_len % 3 != 0:
        print("Invalid sequence: length is not a multiple of three")
        return False, []
    
    codons = [sequence[i:i+3] for i in range(0, seq_len, 3)]
    
    # Checking that all codons in sequence are real codons
    if not all(codon in CODON_TO_AA.keys() for codon in codons):
        print("Invalid sequence: contains non-standard codons")
        return False, []
    
    return True, codons

    
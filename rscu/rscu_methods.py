'''
This file contains helper funnctions to calculate RSCU values for given sequences.
'''

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from lib.aminoacids import AA_TO_CODONS_MULTI_CODON_FAMILIES, CODON_TO_AA, CODON_TO_AA_MULTI_CODON_FAMILIES
from collections import Counter

def run_rscu_analysis(records: list[SeqRecord], 
                      arr_cs_n: np.ndarray, 
                      arr_npet_n: np.ndarray, 
                      w_n: int, 
                      UseCS: bool, 
                      UseSpAsWindow: bool):
    '''
    Run RSCU analysis on a list of SeqRecord objects.
    records: list of SeqRecord objects
    arr_cs: numpy array of cleavage site positions (in nucleotides, multiple of 3), must be same dim as records
    arr_npet: numpy array of npet lengths (in nucleotides, multiple of 3), must be same dim as records
    w_n: window size in nucleotides, multiple of 3
    UseCS: whether to use cleavage site for window extraction
    UseSpAsWindow: whether to use signal peptide as window instead of downstream ROI
    
    Returns: tuple of (window_rscu, full_rscu) where each is a dict of codon to RSCU value
    '''
    
    complement_orfs = []
    windows = []

    # for each record in fasta, extract ROI and complemetn and concatenate with the rest of the ORFs
    for i, record in enumerate(records): # record is SeqRecord object
        
        # find corresponding cs and npet values for this record
        cs_n = int(arr_cs_n[i])
        npet_n = int(arr_npet_n[i])
        
        # extract window and complement
        if UseSpAsWindow:
            extraction = extract_SP_and_complement(record, cs_n)
        else:
            # if UseCS is false, then we set cs_n to 0 to be able to extract windows from within the SP regions
            if not UseCS: cs_n = 0
            extraction = extract_window_and_complement(record, w_n, npet_n, cs_n)
        if extraction is not None:
            window_str, complement_str = extraction
        else:
            print(f"Invalid extraction for record {record.id}. Skipping.")
            continue
        
        valid_window, codons_window = ensure_valid_orf(window_str)
        valid_comp, codons_comp = ensure_valid_orf(complement_str)
        
        # only accept valid windows and full ORFS
        if window_str is None or not valid_window or not valid_comp:
            print(f"Invalid ORF for record {record.id}: complement valid: {valid_comp}, window valid: {valid_window}")
            continue
        else:
            # concatenate codons to full lists
            windows += codons_window
            complement_orfs += codons_comp

    window_rscu = calculate_rscu(windows)
    full_rscu = calculate_rscu(complement_orfs)
    
    return window_rscu, full_rscu
    

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
    record: the BioSeq sequence to be separated into window and complement sequence
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
    '''
    Check a sequence (string) of codons to ensure that it is a valid ORF.
    A valid ORF is defined as:
    - Length is a multiple of 3
    - All codons are real codons (in CODON_TO_AA)
    Returns a tuple of (is_valid, codon_list)
    ''' 
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

    
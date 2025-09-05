'''
Helper methods for working with genetic data
'''

from Bio.SeqRecord import SeqRecord
from lib.aminoacids import CODON_TO_AA


def extract_window_and_complement(record: SeqRecord, w: int, d: int, cs) -> str | None:
    '''
    record: the BioSeq sequence to be separated into window and complement sequence
    w: window size, must be a multiple of three (in nucleotides) 
    d: distance downstream from cleavage site (in nucleotides)
    cs: cleavage site (in nucleotides)
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

    
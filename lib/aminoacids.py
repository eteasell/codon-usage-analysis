'''
Dictionaries for amino acids info etc.
'''

AA_CHARGES = {
    # non-polar
    'G': 0, 
    'A': 0, 
    'V': 0,
    'L': 0,
    'I': 0,
    'M': 0,
    'P': 0,
    # aromatic
    'F': 0,
    'Y': 0,
    'W': 0,
    # polar, non-charged
    'S': 0,
    'T': 0,
    'C': 0,
    'N': 0,
    'Q': 0,
    # basic (+)
    'K': +1,
    'R': +1,
    'H': +1,
    # acidic (-)
    'D': -1,
    'E': -1,
    # extra values in fasta
    '*': 0
}

'''
Ala A
Arg R
Asn N
Asp D
Cys C
Gln Q
Glu E
Gly G
His H
Ile I
Leu L
Lys K
Met M
Phe F
Pro P
Ser S
Thr T
Trp W
Tyr Y
Val V
Stop
'''

# Excluding stop codons
AA_TO_CODONS = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAT', 'AAC'],
    'D': ['GAT', 'GAC'],
    'C': ['TGT', 'TGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'],
    'M': ['ATG'],
    'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG']
}

# Excluding stop codons
CODON_TO_AA = {
    'GCT': 'A', 
    'GCC': 'A', 
    'GCA': 'A', 
    'GCG': 'A',
    'CGT': 'R', 
    'CGC': 'R', 
    'CGA': 'R', 
    'CGG': 'R', 
    'AGA': 'R', 
    'AGG': 'R',
    'AAT': 'N', 
    'AAC': 'N',
    'GAT': 'D', 
    'GAC': 'D',
    'TGT': 'C', 
    'TGC': 'C',
    'CAA': 'Q', 
    'CAG': 'Q',
    'GAA': 'E', 
    'GAG': 'E',
    'GGT': 'G', 
    'GGC': 'G', 
    'GGA': 'G', 
    'GGG': 'G',
    'CAT': 'H', 
    'CAC': 'H',
    'ATT': 'I', 
    'ATC': 'I', 
    'ATA': 'I',
    'TTA': 'L', 
    'TTG': 'L', 
    'CTT': 'L', 
    'CTC': 'L', 
    'CTA': 'L', 
    'CTG': 'L',
    'AAA': 'K', 
    'AAG': 'K',
    'ATG': 'M',
    'TTT': 'F', 
    'TTC': 'F',
    'CCT': 'P', 
    'CCC': 'P', 
    'CCA': 'P', 
    'CCG': 'P',
    'TCT': 'S', 
    'TCC': 'S', 
    'TCA': 'S', 
    'TCG': 'S', 
    'AGT': 'S', 
    'AGC': 'S',
    'ACT': 'T', 
    'ACC': 'T', 
    'ACA': 'T', 
    'ACG': 'T',
    'TGG': 'W',
    'TAT': 'Y', 
    'TAC': 'Y',
    'GTT': 'V', 
    'GTC': 'V', 
    'GTA': 'V', 
    'GTG': 'V'
}

STOP_CODONS = ['TAA', 'TAG', 'TGA']

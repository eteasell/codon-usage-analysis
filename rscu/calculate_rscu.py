'''
The goal of this file is to run the RSCU calculations for the entire gene and within window for each gene in tsv_path.

Run this file by: python -m rscu.calculate_rscu (or else you will get the notorious python-importing-modules issues)
'''

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord  
from ast import literal_eval  
from rscu.rscu import RSCU, extract_window
import json
from rscu import PROJ_DIR

def load_sequence_data(tsv_path, fasta_path) -> tuple[pd.DataFrame, list]:
    df = pd.read_table(tsv_path, index_col='sys_name')
    records = list(SeqIO.parse(fasta_path, "fasta"))
    return df, records

def save_json(rscu_values, path):
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(rscu_values, f, ensure_ascii=False, indent=4)
    
   
# TODO: Set window length and downstream distance 
# Note that this code is running on nucleotide sequences, not amino acids, so make sure that the window and distances
# correspond to those chosen for amino acid analyses by mulitplying by 3
w = 21*3
d = 40*3

tsv_path = PROJ_DIR / 'data/filtered/sp_regions_filtered.tsv'
fasta_path = PROJ_DIR / 'data/filtered/dna_filtered.fasta'

df, records = load_sequence_data(tsv_path=tsv_path, fasta_path=fasta_path)

values = {}

# record is SeqRecord object
for record in records:
    #find corresponding row in tsv file
    row = df.loc[record.id]
    # parse cregion as tuple
    cregion = literal_eval(row['c-region'])
    # extract cleavage site as second entry in tuple
    # note that from this data file, the cleavage site corresponds to the amino acid position
    cs_aa = cregion[1]
    # cs_n is the position of the cleavage site within the nucleotide sequence
    cs_n = cs_aa * 3
    
    # extract window
    window_str = extract_window(record, w, d, cs_n)
    
    if window_str is None:
        continue
    
    # Calculate RSCU for both window and full gene
    full_rscu = RSCU(str(record.seq))
    
    if full_rscu.rscu_values == {}:
        print(f"errors for {record.id} full gene")
        
    window_rscu = RSCU(window_str)
    
    if window_rscu.rscu_values == {}:
        print(f"errors for {record.id} window")
    
    values[record.id] = {
        'full': full_rscu.rscu_values,
        'window': window_rscu.rscu_values
    }

path = PROJ_DIR / "data/rscu_values.json"
save_json(values, path)

################ Additionally, format data as a dataframe ##################
# Dataframe should be indexed on sys_name, with 218 rows, and 122 (2 x 61) columns for each codon (window and background)
codon_labels = [
    "GCT-w", "GCT-b", "GCC-w", "GCC-b", "GCA-w", "GCA-b", "GCG-w", "GCG-b",
    "CGT-w", "CGT-b", "CGC-w", "CGC-b", "CGA-w", "CGA-b", "CGG-w", "CGG-b",
    "AGA-w", "AGA-b", "AGG-w", "AGG-b", "AAT-w", "AAT-b", "AAC-w", "AAC-b",
    "GAT-w", "GAT-b", "GAC-w", "GAC-b", "TGT-w", "TGT-b", "TGC-w", "TGC-b",
    "CAA-w", "CAA-b", "CAG-w", "CAG-b", "GAA-w", "GAA-b", "GAG-w", "GAG-b",
    "GGT-w", "GGT-b", "GGC-w", "GGC-b", "GGA-w", "GGA-b", "GGG-w", "GGG-b",
    "CAT-w", "CAT-b", "CAC-w", "CAC-b", "ATT-w", "ATT-b", "ATC-w", "ATC-b",
    "ATA-w", "ATA-b", "TTA-w", "TTA-b", "TTG-w", "TTG-b", "CTT-w", "CTT-b",
    "CTC-w", "CTC-b", "CTA-w", "CTA-b", "CTG-w", "CTG-b", "AAA-w", "AAA-b",
    "AAG-w", "AAG-b", "ATG-w", "ATG-b", "TTT-w", "TTT-b", "TTC-w", "TTC-b",
    "CCT-w", "CCT-b", "CCC-w", "CCC-b", "CCA-w", "CCA-b", "CCG-w", "CCG-b",
    "TCT-w", "TCT-b", "TCC-w", "TCC-b", "TCA-w", "TCA-b", "TCG-w", "TCG-b",
    "AGT-w", "AGT-b", "AGC-w", "AGC-b", "ACT-w", "ACT-b", "ACC-w", "ACC-b",
    "ACA-w", "ACA-b", "ACG-w", "ACG-b", "TGG-w", "TGG-b", "TAT-w", "TAT-b",
    "TAC-w", "TAC-b", "GTT-w", "GTT-b", "GTC-w", "GTC-b", "GTA-w", "GTA-b",
    "GTG-w", "GTG-b"
]

mat = np.zeros((218,122))
gene_list = []

for i, gene in enumerate(values.keys()):
    gene_list.append(gene)
    
    gene_i = values[gene]
    
    window = dict(gene_i['window'])
    full = gene_i['full']
    
    for j, codon in enumerate(window.keys()):
        mat[i, 2*j] = window[codon]
        mat[i, 2*j + 1] = full[codon]
        

mat = np.round(mat,decimals=3)
df = pd.DataFrame(mat, index=gene_list, columns=codon_labels)

print(df)

df_path = PROJ_DIR / "data/rscu_values.tsv"
df.to_csv(df_path, sep="\t", index=True)
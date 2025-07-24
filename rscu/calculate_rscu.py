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

tsv_path = PROJ_DIR / 'data/sp_regions_slowmode.tsv'
fasta_path = PROJ_DIR / 'data/sp_dna_slowmode.fasta'

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
        'full': window_rscu.rscu_values,
        'window': full_rscu.rscu_values
    }

path = PROJ_DIR / "data/rscu_values.json"
save_json(values, path)
    

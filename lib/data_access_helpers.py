'''
Static functions for accessing data files.
'''

import pandas as pd
from Bio import SeqIO

def load_elongation_rates():
    path = '../data/daoduc_song_2018_rates.csv'
    return pd.read_csv(path, index_col='Codon')
    

def load_rscu_data(path):
    return pd.read_csv(path, sep='\t', index_col='codon')
    

def save_rscu(path, window_rscu, full_rscu):
    
    df = pd.DataFrame({
        "window": window_rscu,
        "complement": full_rscu
    })
    df.index.name = 'codon'
    df.to_csv(path, index=True, sep='\t')
    print(f"Saved RSCU data to {path}")
    
    
def load_fasta_records(path: str, genes: list[str] | None):
    records = list(SeqIO.parse(path, 'fasta'))
    
    if genes is not None:
        filtered_records = []
        for record in records:
            if record.id in genes:
                filtered_records.append(record)
        return filtered_records
    else:
        return records
                
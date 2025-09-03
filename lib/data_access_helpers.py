'''
Static functions for accessing data files.
'''

import pandas as pd

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
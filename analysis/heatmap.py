'''
The goal of this file to create a heatmap for the number of positively and negatively charged amino acids within a 
window of size n starting at distance d from the protein's SP cleavage site.
'''

import pandas as pd
from Bio import SeqIO
from ast import literal_eval
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from lib.aminoacids import AA_CHARGES
from analysis import PROJ_DIR

# window size
w = 21

# distance from cleavage site
d = 40

df_path = PROJ_DIR / 'data/filtered/sp_regions_filtered.tsv'
df = pd.read_table(df_path, index_col='sys_name')
print(df)

fasta_path = PROJ_DIR / 'data/filtered/proteins_filtered.fasta'
records = list(SeqIO.parse(fasta_path, "fasta"))
print(f"fasta records length: {len(records)}")

cregions = df['c-region'].apply(literal_eval)

rows = []
too_short = []

for record in records:
    cregion = cregions[record.id]
    
    #cleavage site
    cs = int(cregion[1])
    
    seq = str(record.seq)

    if cs+d+w > len(seq):
        print(f"Window too long for protein {record.id}")
        too_short.append(record.id)
        continue
    
    window = seq[cs+d:cs+d+w]
    sp_region = seq[0:cs]
    sp_downstream = seq[0:cs+d+w]
    downstream_only = seq[cs:]
    
    # TODO: choose from the options line 75 to 78, then don't forget to update the figure title
    seq_region = window
    charges = [AA_CHARGES[aa] for aa in seq_region]
    
    pos_count = charges.count(+1)
    neg_count = charges.count(-1)
    
    region = len(seq_region)
    
    pos_freq = pos_count / region
    neg_freq = neg_count / region
    
    row = {
        'sys_name': record.id,
        'cleavage_site': cs,
        'pos_freq': pos_freq,
        'neg_freq': neg_freq
    }
    
    rows.append(row)
    
print(f"Too Short Proteins: {len(too_short)}")
print(too_short)
    
df_heat = pd.DataFrame(rows)

print(f"Max freqs (pos and neg): {df_heat['pos_freq'].max()} and {df_heat['neg_freq'].max()}")

# Define bin width and count
bin_count = 10
freq_max = 0.5
bin_width = freq_max / bin_count  # = 0.05

# Discretize into 10 bins from 0.0 to 0.5
df_heat['neg_bin'] = np.floor(df_heat['neg_freq'] / bin_width).astype(int)
df_heat['pos_bin'] = np.floor(df_heat['pos_freq'] / bin_width).astype(int)

# Clip any values >= 0.5 into bin 9
df_heat['neg_bin'] = df_heat['neg_bin'].clip(0, bin_count - 1)
df_heat['pos_bin'] = df_heat['pos_bin'].clip(0, bin_count - 1)

print(df_heat)

# Create heatmap data
heatmap_data = (
    df_heat.groupby(['neg_bin', 'pos_bin'])
    .size()
    .unstack(fill_value=0)
    .reindex(index=range(bin_count), columns=range(bin_count), fill_value=0)
)

# Create bin labels for plotting
bin_edges = np.linspace(0, freq_max, bin_count + 1)
bin_labels = [f"{bin_edges[i]:.2f}–{bin_edges[i+1]:.2f}" for i in range(bin_count)]
heatmap_data.index = bin_labels  # neg_freq (y-axis)
heatmap_data.columns = bin_labels  # pos_freq (x-axis)

print(heatmap_data)

plt.figure(figsize=(8, 6))
sns.heatmap(heatmap_data, annot=True, fmt="d", cmap="YlGnBu").invert_yaxis()

# TODO: change title
plt.title(f"Amino acid charge frequencies within window (w={w}, d={d})")
plt.xlabel("Frequency of (+) charges")
plt.ylabel("Frequency of (-) charges")
plt.tight_layout()
plt.show()
'''
Running a paired t-test
'''

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel
from analysis import PROJ_DIR

# first load our RSCU data
rscu_path = PROJ_DIR / "data/rscu_values.tsv"
df = pd.read_table(rscu_path, index_col=0)

# Now we just have a numpy array to work with
arr = df.values

t_results = np.zeros((61), dtype=float)
p_results = np.zeros((61), dtype=float)

for i in range(0, 61):
    w_col = arr[:, 2*i]
    b_col = arr[:, 2*i + 1]
    
    t_stat, p_val = ttest_rel(b_col, w_col)
    
    t_results[i] = t_stat
    p_results[i] = p_val
    
# Save results by codon

codons = [
    "GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "CGA", "CGG", "AGA", "AGG",
    "AAT", "AAC", "GAT", "GAC", "TGT", "TGC", "CAA", "CAG", "GAA", "GAG",
    "GGT", "GGC", "GGA", "GGG", "CAT", "CAC", "ATT", "ATC", "ATA", "TTA",
    "TTG", "CTT", "CTC", "CTA", "CTG", "AAA", "AAG", "ATG", "TTT", "TTC",
    "CCT", "CCC", "CCA", "CCG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
    "ACT", "ACC", "ACA", "ACG", "TGG", "TAT", "TAC", "GTT", "GTC", "GTA",
    "GTG"
]

data = {
    't-statistic': t_results,
    'p-value': p_results
}

df = pd.DataFrame(data, index=codons)
out_path = PROJ_DIR / "data/stats/t_test.csv"
df.to_csv(out_path, index=True)
    




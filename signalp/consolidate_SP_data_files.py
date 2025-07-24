'''
The goal is to collect the results from running SignalP on the seven batches into one csv included standard name, systematic name, 
sgid, confidence, and location of SP region
'''

import pandas as pd
from gffpandas import gffpandas as gffpd
from signalp import PROJ_DIR

gff_files = [PROJ_DIR / f'data/signalp/out_trans_{i}/output.gff3' for i in range(7)]

dfs = []

for path in gff_files:
    gdf = gffpd.read_gff3(path)
    df = gdf._read_gff3_to_df()
    dfs.append(df)

combined_df = pd.concat(dfs, ignore_index=True)

seqs = combined_df['seq_id'].str.strip()

# Extract sys_name, std_name, and SGDID from seq_id
sys_name = seqs.str.split().str[0]
std_name = seqs.str.split().str[1]
sgdid = seqs.str.extract(r"(SGDID_[^\s,]+)")[0]

# Combine with selected fields
final_df = pd.concat([
    sys_name.rename("sys_name"),
    std_name.rename("std_name"),
    sgdid.rename("SGDID"),
    combined_df[['start', 'end', 'score']]
], axis=1)

# Save to a tab-separated file
final_df.to_csv(PROJ_DIR / "data/all_sp_loci.tsv", sep="\t", index=False)

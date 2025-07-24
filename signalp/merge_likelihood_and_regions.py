'''
The goal is to merge the separate rows per region in the sp_proteins_region_output.gff file into one row per protein in a tsv
'''

import pandas as pd
from gffpandas import gffpandas as gffpd
from signalp import PROJ_DIR

output_path = PROJ_DIR / 'data/sp_proteins_out/output.gff3'

gdf = gffpd.read_gff3(output_path)
df_meta = gdf._read_gff3_to_df()
print(df_meta)

seqs = df_meta['seq_id'].str.strip()

# Extract sys_name, std_name, and SGDID from seq_id
sys_name = seqs.str.split().str[0]
std_name = seqs.str.split().str[1]
sgdid = seqs.str.extract(r"(SGDID_[^\s,]+)")[0]

# Add parsed fields to the DataFrame
df_meta['sys_name'] = sys_name
df_meta['std_name'] = std_name
df_meta['SGDID'] = sgdid

# Keep only relevant columns
df_meta = df_meta[['sys_name', 'std_name', 'SGDID', 'score']]
df_meta = df_meta.rename(columns={'score':'prediction_likelihood'})

print(df_meta)

'''
Now we want to access the regions file and add in the n,h,c region boundaries
'''
region_path = PROJ_DIR / 'data/sp_proteins_out/region_output.gff3'

gdf = gffpd.read_gff3(region_path)
df_regions = gdf._read_gff3_to_df()

# Extract SGDID from seq_id
df_regions['SGDID'] = df_regions['seq_id'].str.extract(r"(SGDID_[^\s,]+)")

# Create (start, end) tuple
df_regions['region'] = list(zip(df_regions['start'], df_regions['end']))

# Keep only relevant columns
df_regions = df_regions[['SGDID', 'type', 'region']]

# Pivot so that each type becomes a column
df_regions = df_regions.pivot(index='SGDID', columns='type', values='region').reset_index()

# reorder columns
df_regions = df_regions[['SGDID', 'n-region', 'h-region', 'c-region']]
print(df_regions.columns)

df_merged = pd.merge(df_meta, df_regions, on="SGDID", how="left")

# reorder columns
df_merged = df_merged[['SGDID', 'sys_name', 'std_name', 'prediction_likelihood', 'n-region', 'h-region', 'c-region']]

print(df_merged)

# flatten columns and save
out_path = PROJ_DIR / 'data/sp_regions_slowmode.tsv'
df_merged.to_csv(out_path, sep="\t", index=False)
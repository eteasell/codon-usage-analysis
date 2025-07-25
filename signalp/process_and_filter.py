'''
This file takes the returned data from the signalp run (slowmode), filters out low confidences 
(if overall likelihood or cleavage-site position likelihood is less than 0.95), and consolidates into 
one tsv file with the following fields:

sys_name	std_name	prediction_likelihood	n-region	h-region	c-region    cs      cs_likelihood
'''

import pandas as pd
from gffpandas import gffpandas as gffpd
from signalp import PROJ_DIR

'''
First, read the slowmode output signalp file to get gene names and prediction likelihood
I am choosing to use this file over the prediction_results.txt file below because this file is more easily understood by 
Pandas with minimal manipulation, as opposed to using regex as below.
'''
signalp_out = PROJ_DIR / 'data/signalp/sp_proteins_out/output.gff3'

gdf = gffpd.read_gff3(signalp_out)
df_signalp_out = gdf._read_gff3_to_df()

seqs = df_signalp_out['seq_id'].str.strip()

# Extract sys_name, std_name, and SGDID from seq_id
sys_name = seqs.str.split().str[0]
std_name = seqs.str.split().str[1]
sgdid = seqs.str.extract(r"(SGDID_[^\s,]+)")[0]

# Add parsed fields to the DataFrame
df_signalp_out['sys_name'] = sys_name
df_signalp_out['std_name'] = std_name
df_signalp_out['SGDID'] = sgdid

# Keep only relevant columns
df_signalp_out = df_signalp_out[['sys_name', 'std_name', 'SGDID', 'score']]
df_signalp_out = df_signalp_out.rename(columns={'score':'prediction_likelihood'})

print(df_signalp_out)

# Now we want to access the regions file and add in the n,h,c region boundaries

region_path = PROJ_DIR / 'data/signalp/sp_proteins_out/region_output.gff3'

gdf = gffpd.read_gff3(region_path)
df_regions = gdf._read_gff3_to_df()

# Extract sys_name from seq_id
df_regions['sys_name'] = df_regions['seq_id'].apply(lambda x: x.split(" ")[0])

# Create (start, end) tuple
df_regions['region'] = list(zip(df_regions['start'], df_regions['end']))

# Keep only relevant columns
df_regions = df_regions[['sys_name', 'type', 'region']]

# Pivot so that each type becomes a column
df_regions = df_regions.pivot(index='sys_name', columns='type', values='region').reset_index()

# merge with data from the original output file
df_merged = pd.merge(df_signalp_out, df_regions, on="sys_name", how="left")

# reorder columns
df_merged = df_merged[['sys_name', 'std_name', 'SGDID', 'prediction_likelihood', 'n-region', 'h-region', 'c-region']]

print(df_merged)

'''
Below, I access the prediction_results.txt file which contains information about the cleavage site and likelihood.
The code below merges this data into the dataframe
'''

txt_path = PROJ_DIR / 'data/signalp/sp_proteins_out/prediction_results.txt'

df_tsv = df_merged.set_index('sys_name')

print(df_tsv)

df_txt = pd.read_csv(
    txt_path,
    sep="\t",
    comment="#",
    header=None,
    engine="python",
    names=["ID", "Prediction", "OTHER", "SP(Sec/SPI)", "CS Position"]
)

# manipulating the parsed txt file
df_txt["sys_name"] = df_txt["ID"].str.split().str[0]
df_txt["cs"] = df_txt["CS Position"].str.extract(r"CS pos: (\d+)-")[0].astype(float)
df_txt["cs_likelihood"] = df_txt["CS Position"].str.extract(r"Pr: ([0-9.]+)")[0].astype(float)

# keep only the nicely formatted columns
df_cs = df_txt[["sys_name", "cs", "cs_likelihood"]].set_index("sys_name")

# combine CS data with ongoing dataframe
df_combined = df_tsv.merge(df_cs, left_index=True, right_index=True, how="left")

# Filtering step for outliers: remove any genes with either overall likelihood or CS position likelihood less than 0.95
df_filtered = df_combined[df_combined['prediction_likelihood'] >= 0.95]
df_filtered_on_cs = df_filtered[df_filtered['cs_likelihood'] >= 0.95]

OUT_PATH = PROJ_DIR / "data/filtered/sp_regions_filtered.tsv"
df_filtered_on_cs.to_csv(OUT_PATH, sep="\t", index=True)
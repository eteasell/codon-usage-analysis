'''
The goal is to filter the original fasta files to include only the selected genes by SignalP ()
'''

from Bio import SeqIO
import pandas as pd
from signalp import PROJ_DIR

dna_path = PROJ_DIR / "data/orf_genomic_all.fasta"
protein_path = PROJ_DIR / "data/signalp/orf_trans_all.fasta"
tsv_path = PROJ_DIR / "data/filtered/sp_regions_filtered.tsv"

sys_name_values = pd.read_csv(tsv_path, sep="\t")["sys_name"].tolist()

print(sys_name_values)
print(len(sys_name_values))

# filter DNA file
filtered_dna = []

for record in SeqIO.parse(dna_path, "fasta"):
    sys_name = str(record.id.split(" ")[0])
    if sys_name in sys_name_values:
        filtered_dna.append(record)
        
print(len(filtered_dna))
out_path = PROJ_DIR / f"data/filtered/dna_filtered.fasta"
SeqIO.write(filtered_dna, out_path, "fasta")

# filter Proteins file
filtered_proteins = []

for record in SeqIO.parse(protein_path, "fasta"):
    sys_name = str(record.id.split(" ")[0])
    if sys_name in sys_name_values:
        filtered_proteins.append(record)
        
print(len(filtered_proteins))
out_path = PROJ_DIR / f"data/filtered/proteins_filtered.fasta"
SeqIO.write(filtered_proteins, out_path, "fasta")
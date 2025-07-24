'''
The goal is to filter the original all-protein fasta file to include only the 312 hits from the SignalP 6.0 fast mode run
'''

from Bio import SeqIO
import pandas as pd
from signalp import PROJ_DIR

fasta_path = PROJ_DIR / "data/orf_genomic_all.fasta"
tsv_path = PROJ_DIR / "data/sp_regions_slowmode.tsv"

sys_name_values = pd.read_csv(tsv_path, sep="\t")["sys_name"].tolist()

print(sys_name_values)
print(len(sys_name_values))

filtered = []

for record in SeqIO.parse(fasta_path, "fasta"):
    sys_name = str(record.id.split(" ")[0])
    if sys_name in sys_name_values:
        filtered.append(record)
        
print(len(filtered))
out_path = PROJ_DIR / f"data/sp_rna_slowmode.fasta"
SeqIO.write(filtered, out_path, "fasta")
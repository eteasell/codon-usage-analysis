from Bio import SeqIO
from pathlib import Path
from signalp import PROJ_DIR

# specific to the current yeast genome file
number_of_proteins = 6722

path = PROJ_DIR / "data/signalp/orf_trans_all.fasta"

bins = {i: [] for i in range(7)}

def num_to_rec(i):
    return (i // 1000)

j = 0
for i,record in enumerate(SeqIO.parse(path, "fasta")):
    bins[num_to_rec(i)].append(record)
    j = i
print(j)
    
for key in bins:
    print(len(bins[key]))
    out_path = PROJ_DIR / f"data/signalp/orf_trans_{key}.fasta"
    SeqIO.write(bins[key], out_path, "fasta")
    
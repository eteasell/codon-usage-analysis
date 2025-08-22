import sys
import os

# Add the directory containing your modules to Python path
sys.path.append('/Users/ellateasell/Research/CodonUsageBias/code')

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from rscu_methods import extract_window, ensure_valid_orf, calculate_rscu, rscu_to_string, extract_window_and_complement
from rscu import PROJ_DIR 
import pandas as pd
import numpy as np
from lib.aminoacids import AA_TO_CODONS_MULTI_CODON_FAMILIES, SLOW, FAST


def run_rscu_analysis(records: list[SeqRecord], df: pd.DataFrame, w_n, d_n, csUsed):
    if csUsed:
        print(f"Using cleavage site with window size {w_n} in nucelotides and distance {d_n} in nucelotides from CS")
    else:
        print(f"Not using cleavage site with window size {w_n} in nucleotides and distance {d_n} in nucleotides from start")
        

    full_orfs = []
    windows = []

    # for each record in fasta, extract ROI and full and concatenate with the rest of the ORFs
    for record in records: # record is SeqRecord object
        
        #find corresponding row in tsv file
        row = df.loc[record.id]
        
        # extract cleavage site
        # note that from this data file, the cleavage site corresponds to the amino acid position
        cs_aa = row['cs']
        # cs_n is the position of the cleavage site within the nucleotide sequence
        cs_n = int(cs_aa * 3)
        
        # extract window and complement
        if not csUsed: cs_n = 0
        extraction = extract_window_and_complement(record, w_n, d_n, cs_n)
        if extraction is not None:
            window_str, comeplement_str = extraction
        else:
            print(f"Invalid extraction for record {record.id}. Skipping.")
            continue
        valid_window, codons_window = ensure_valid_orf(window_str)
        valid_comp, codons_comp = ensure_valid_orf(comeplement_str)
        
        # only accept valid windows and full ORFS
        if window_str is None or not valid_window or not valid_comp:
            print(f"Invalid ORF for record {record.id}: complement valid: {valid_comp}, window valid: {valid_window}")
            continue
        else:
            windows += codons_window
            full_orfs += codons_comp

    window_rscu = calculate_rscu(windows)
    full_rscu = calculate_rscu(full_orfs)
    
    return window_rscu, full_rscu

def sum_of_abs_differences_fast_slow(rscu_w, rscu_b, codons):
    diffs = abs(rscu_w - rscu_b)
    slow_sum = sum(diffs[i] for i, codon in enumerate(codons) if codon in SLOW)
    fast_sum = sum(diffs[i] for i, codon in enumerate(codons) if codon in FAST)
    return slow_sum, fast_sum



##################################

# load SP ORFS (genomic)
all_fasta_path = PROJ_DIR / 'data/filtered/dna_filtered.fasta'
records = list(SeqIO.parse(all_fasta_path, 'fasta'))

# load SP data
tsv_path = PROJ_DIR / 'data/filtered/sp_regions_filtered.tsv'
df = pd.read_table(tsv_path, index_col='sys_name')

max_diff_slow = 0
best_w_aa_slow = 0
best_d_aa_slow = 0

max_diff_fast = 0
best_w_aa_fast = 0
best_d_aa_fast = 0

max_diff_overall = 0
best_w_aa_overall = 0
best_d_aa_overall = 0

for w_aa in range(10, 50, 5): # searching from window size (in amino acids) from 10, 15, 20, ...,50
    for d_aa in range(0, 100, 10): # searching in distance from CS (in amino acids) from 0, 10, 20, ..., 100
        
        w_n = w_aa * 3
        d_n = d_aa * 3
        csUsed = True
        
        window_rscu, full_rscu = run_rscu_analysis(records, df, w_n, d_n, csUsed)
        
        rscu_w = np.array(list(window_rscu.values()))
        rscu_b = np.array(list(full_rscu.values()))
        
        codons = list(window_rscu.keys())
        
        slow_sum, fast_sum = sum_of_abs_differences_fast_slow(rscu_w, rscu_b, codons)
        
        if slow_sum > max_diff_slow:
            max_diff_slow = slow_sum
            best_w_aa_slow = w_aa
            best_d_aa_slow = d_aa
            
        if fast_sum > max_diff_fast:
            max_diff_fast = fast_sum
            best_w_aa_fast = w_aa
            best_d_aa_fast = d_aa
            
        if slow_sum + fast_sum > max_diff_overall:
            max_diff_overall = slow_sum + fast_sum
            best_w_aa_overall = w_aa
            best_d_aa_overall = d_aa
            
            
   
print("")
print("------------------------------")
print("Best parameters found for SLOW CODONS:")         
print(f"Best window size (in amino acids): {best_w_aa_slow}, Best distance from CS (in amino acids): {best_d_aa_slow}, Max sum of abs differences: {max_diff_slow:.3f}")

print("------------------------------")
print("Best parameters found for FAST CODONS:")
print(f"Best window size (in amino acids): {best_w_aa_fast}, Best distance from CS (in amino acids): {best_d_aa_fast}, Max sum of abs differences: {max_diff_fast:.3f}")

print("------------------------------")
print("Best parameters found overall:")
print(f"Best window size (in amino acids): {best_w_aa_overall}, Best distance from CS (in amino acids): {best_d_aa_overall}, Max sum of abs differences: {max_diff_overall:.3f}")
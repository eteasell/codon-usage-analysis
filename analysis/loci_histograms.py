'''
Preliminary statistics for SP loci within genes
'''

import pandas as pd
from ast import literal_eval
import matplotlib.pyplot as plt
from analysis import PROJ_DIR

path = PROJ_DIR / 'data/sp_regions_slowmode.tsv'
df = pd.read_table(path)

print(df)

# Note to self: literal_eval is a useful method that turns strings of python code into actual python code
n_region_interval = df['n-region'].apply(literal_eval)
print(n_region_interval)

print(type(n_region_interval[0]))

n_region_starts = n_region_interval.apply(lambda x: x[0])
print(type(n_region_starts[0]))

# plt.hist(n_region_starts, bins=30, edgecolor='black')
# plt.xlabel('Start Position')
# plt.ylabel('Frequency')
# plt.title('Histogram of Region Start Positions')
# plt.show()

h_region_interval = df['h-region'].apply(literal_eval)
h_length = h_region_interval.apply(lambda x: x[1]-x[0])

c_region_interval = df['c-region'].apply(literal_eval)
c_length = c_region_interval.apply(lambda x: x[1]-x[0])

n_length = n_region_interval.apply(lambda x: x[1]-x[0])

total_lengths = h_length + c_length + n_length

# plt.hist(n_length, edgecolor='black')
# plt.xlabel('Length of n-region (amino acids)')
# plt.ylabel('Frequency')
# plt.title('Histogram of n-region lengths')
# plt.show()

# plt.hist(h_length, edgecolor='black')
# plt.xlabel('Length of h-region (amino acids)')
# plt.ylabel('Frequency')
# plt.title('Histogram of h-region lengths')
# plt.show()

plt.hist(c_length, edgecolor='black')
plt.xlabel('Length of c-region (amino acids)')
plt.ylabel('Frequency')
plt.title('Histogram of c-region lengths')
plt.show()

# plt.hist(total_lengths, edgecolor='black')
# plt.xlabel('Length of Signal Peptide (amino acids)')
# plt.ylabel('Frequency')
# plt.title('Histogram of Signal Peptide Lengths')
# plt.show()
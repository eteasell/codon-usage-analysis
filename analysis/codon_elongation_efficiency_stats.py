'''
Compare elongation efficiencies from Duc KD, Song YS (2018) to the Subramanian et al. 2022 database
'''
from analysis import PROJ_DIR

codon_db = PROJ_DIR / 'literature/Subramanian et al. 2022/nuclear_codon_statistics.tsv'
plos_db = PROJ_DIR / 'literature/Dao Duc & Song 2018/supp_data_updated.xlsx'

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df_codon = pd.read_table(codon_db, skiprows=8)
df_plos = pd.read_excel(plos_db, sheet_name=3)

df_merged = pd.merge(df_plos, df_codon, left_on="Codon", right_on="CODON", how="left")

df = df_merged[['Codon', 'Amino acid', 'median_elongation_rate', 'Preferred Codon']]

print(df)

# Get unique amino acids
amino_acids = df["Amino acid"].unique()

# Set up color palette
palette = {True: '#4C72B0', False: '#DD8452'}

# Create subplots: 5 rows × 4 columns = 20 panels
fig, axes = plt.subplots(nrows=5, ncols=4, figsize=(20, 15))
axes = axes.flatten()  # Make it easy to index

# Plot each amino acid separately
for i, aa in enumerate(sorted(amino_acids)):
    if i >= len(axes):
        break
    ax = axes[i]
    df_aa = df[df["Amino acid"] == aa].copy()

    # Sort codons by rate
    df_aa = df_aa.sort_values("median_elongation_rate", ascending=False)

    # Plot
    sns.barplot(
        data=df_aa,
        x="Codon",
        y="median_elongation_rate",
        hue="Preferred Codon",
        palette=palette,
        ax=ax
    )

    ax.set_title(aa)
    ax.set_ylabel("Elongation Rate")
    ax.set_xlabel("")
    ax.tick_params(axis='x', rotation=45)
    ax.legend_.remove()  # Remove individual legends to avoid clutter

# Hide unused subplot panels if < 20 amino acids
for j in range(i+1, len(axes)):
    axes[j].axis('off')

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, title="Preferred", loc='lower center', ncol=2)

fig.tight_layout(rect=[0, 0.05, 1, 1]) 
plt.show()
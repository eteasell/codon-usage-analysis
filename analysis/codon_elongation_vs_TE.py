'''
Compare elongation efficiencies from Duc KD, Song YS (2018) to the Pechmann et al. 2014 database
'''
from analysis import PROJ_DIR

pechmann_db = PROJ_DIR / 'literature/Pechmann et al. 2014/41594_2014_BFnsmb2919_MOESM20_ESM.xlsx'
plos_db = PROJ_DIR / 'literature/Dao Duc & Song 2018/supp_data.xlsx'
codon_db = PROJ_DIR / 'literature/Subramanian et al. 2022/nuclear_codon_statistics.tsv'

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df_pechmann = pd.read_excel(pechmann_db, skiprows=12)
df_plos = pd.read_excel(plos_db, sheet_name=3)
df_codon = pd.read_table(codon_db, skiprows=8)

# Below is DAO DUC vs PECHMANN
df_merged = pd.merge(df_plos, df_pechmann, left_on="Codon", right_on="Codon", how="left")

df = df_merged[['Codon', 'median_elongation_rate', 'S. cerevisiae']]
df = df.rename(columns={'S. cerevisiae': "translation_efficiency"})

print(df)

sns.scatterplot(data=df, x="median_elongation_rate", y="translation_efficiency")

# Below is DAO DUC vs SUBRAMANIAN
df_merged = pd.merge(df_plos, df_codon, left_on="Codon", right_on="CODON", how="left")
df = df_merged[['Codon', 'median_elongation_rate', 'RSCU']]

print(df)

sns.scatterplot(data=df, x="median_elongation_rate", y="RSCU")

# Below is PECHMANN vs SUBRAMANIAN
df_merged = pd.merge(df_pechmann, df_codon, left_on="Codon", right_on="CODON", how="left")
df = df_merged[['Codon', 'S. cerevisiae', 'RSCU']]
df = df.rename(columns={'S. cerevisiae': "translation_efficiency"})

print(df)

sns.scatterplot(data=df, x="translation_efficiency", y="RSCU")

plt.show()
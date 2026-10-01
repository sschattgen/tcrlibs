# this started with Supplementary File 2-P25NSTD02_Sequence_20250416
# VCV nucleotide sequences were set as sequence 
# now we'll trim the targets to match the read lengths.

#%%
import pandas as pd

#%%
df = pd.read_csv('/Users/sschattg/Downloads/HA/library_table.csv')

# %%
seed_seq = 'gaggacctgaacaaggtgtttcctccagaggtggccgtgttc'.upper()
seed_seq_len = len(seed_seq)
keep_length = 150 - seed_seq_len

# %%
# iterate through the Sequence column, find the seed_seq, and then pull the characters of keep_length preceeding it.
split_seqs = df['Sequence'].apply(lambda x: x.split(seed_seq))

# %%
target_seqs = list()
for i in range(len(split_seqs)):
    target_seqs.append(split_seqs[i][0][-keep_length:])

# %%
df2 = df.copy()
df2 = df2[['ID','Sequence']]
df2['Sequence'] = target_seqs
df2.to_csv('/Users/sschattg/Downloads/HA/library_table_clean.csv', index=False)
# %%

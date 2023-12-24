#!/usr/bin/env/python

import pandas as pd
import numpy as np
import os

def find_erroneous_and_missing_buscos(syn, merian_refs):
    number_seqs = []
    syn_results = pd.read_csv(syn, sep='\t', low_memory=False)
    syn_results = syn_results.rename(columns={"#marker": "busco"})
    n2_results = syn_results[['busco', 'n2_seq']].dropna()
    n2_results.head(100)
    n2_results = n2_results.merge(merian_refs, how="inner", on="busco")
    dominant_merians = n2_results.groupby(['n2_seq'])['merian'].agg(pd.Series.mode).to_frame()
    dominant_merians = dominant_merians.rename(columns={"merian": "dominant_merian"})
    n2_results = n2_results.merge(dominant_merians, how="inner", on="n2_seq")
    print(n2_results.loc[~(n2_results['merian'] == n2_results['dominant_merian'])])
    error_buscos = np.where(n2_results.merian != n2_results.dominant_merian)
    print(error_buscos)
    missing_buscos = set(merian_refs['busco']) -  set(n2_results['busco'])
    number_seqs.append(len(n2_results['n2_seq'].unique()))
    return(error_buscos, missing_buscos, number_seqs)
# %%
merian_refs = pd.read_csv('Final_Merians_based_on_syn_n2_from_r2_m5_full_table.tsv', sep='\t', 
header=None, usecols=[0,2], names=['busco', 'merian'])

# %%
file_dir = '/lustre/scratch123/tol/teams/blaxter/projects/lepidoptera_genomics/cw22/Leps_200/Analysis/syngraph/reviewers_response/subsampled_pickles/'
files = os.listdir(file_dir) #'.' = directory path
table_files = [k for k in files if '.table.tsv' in k]
print(len(table_files))
# %%
total_error, total_missing, total_number_seqs = [], [], []
for table in table_files:
    table_path = file_dir + table
    error, missing, number_seqs = find_erroneous_and_missing_buscos(table_path, merian_refs)
    total_missing.extend(list(missing))
    total_error.extend([error])
    total_number_seqs.append(number_seqs)
# %%
# count average number occurences per busco in totaL_missing
pd.Series(total_missing).value_counts().mean() # average 12% of time missing
# %%
# count number occurences per busco in total_error
pd.Series(total_error) # empty

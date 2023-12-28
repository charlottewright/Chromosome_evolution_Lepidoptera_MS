#!/usr/bin/python

#%%
import random
import pandas as pd
from itertools import product
import seaborn as sns
import statistics
import math
import matplotlib.pyplot as plt
from statistics import mean
from statistics import median
# %%
# Load the CSV data
df = pd.read_csv('single_copy_orthogroups_with_ALGs_assigned_Merians_accurately.tsv', sep='\t', header=None)
df = df.iloc[:, [3, 0, 1, 5]]  # Select columns 4, 1, 2, and 6
df.columns = ['BCS_LG', 'Aquee_prot', 'Mcinx_prot', 'Merian']

# BLGs (bilarerian LGs) transformations
df['BLG'] = df['BCS_LG'].str.replace('A1a', 'A1_mixed').str.replace('A1b', 'A1_mixed').str.replace('A1_mixed', 'A1axA1b')
df['BLG'] = df['BLG'].str.replace('Ea', 'E_mixed').str.replace('Eb', 'E_mixed').str.replace('E_mixed', 'EaxEb')
df['BLG'] = df['BLG'].str.replace('Qa', 'Q_mixed').str.replace('Qb', 'Q_mixed').str.replace('Q_mixed', 'QaxQb')
df['BLG'] = df['BLG'].str.replace('Qc', 'Q_mixed').str.replace('Qd', 'Q_mixed').str.replace('Q_mixed', 'QcxQd')

# Group and summarize by BLG and Merian
df = df.groupby(['BLG', 'Merian'])['BCS_LG'].count().reset_index()
df.columns = ['BLG', 'Merian', 'matches']

# Create a df with all possible combos of Merin & BLGs
all_combinations = pd.DataFrame(list(product(df['BLG'].unique(), df['Merian'].unique())), columns=['BLG', 'Merian'])
expanded_df = pd.merge(df, all_combinations, on=['BLG', 'Merian'], how='right')
expanded_df['matches'].fillna(0, inplace=True)

# Group by BLG and Merian to calculate totals
df = expanded_df.groupby('BLG')['matches'].sum().reset_index()
df.columns = ['BLG', 'total_BLG']
expanded_df = pd.merge(expanded_df, df, on='BLG', how='left')
df = expanded_df.groupby('Merian')['matches'].sum().reset_index()
df.columns = ['Merian', 'total_Merian']
df = pd.merge(expanded_df, df, on='Merian', how='left')

# Calculate proportions
df['prop_BLG'] = (df['matches'] / df['total_BLG']) * 100
df['prop_Merian'] = (df['matches'] / df['total_Merian']) * 100

BLGs = df.drop_duplicates(subset=["BLG"])
BLG_values = BLGs['BLG'].to_list()
BLG_weights = BLGs['total_BLG'].to_numpy()
merian2num_markers =  df.set_index('Merian').to_dict()['total_Merian']
# %%
def simulate_variance(df, number_simulations, weights):
  merian2sig = {}
  for m in list(df['Merian'].unique()): # for each merian
  #for m in ['M13']:
    random_var = []
    number_markers = int(merian2num_markers[m])
    if number_markers < 3:
      continue
    else:
      total_sims = number_simulations + 1
      for i in range(1,total_sims,1):
        if weights == "False":
          randomList = random.choices(BLG_values, k=number_markers) # no weights
        else:
          randomList = random.choices(BLG_values, BLG_weights, k=number_markers)
        randomCounts = pd.Series(randomList).value_counts().to_list()
        if len(randomCounts) != len(BLG_values):
          number_missing_zeros = len(BLG_values) -len(randomCounts)
          for i in range(1,(number_missing_zeros + 1), 1):
            randomCounts.append(0) # need to account for zeros i.e. BLGs that weren't sampled at all
        var = statistics.variance(randomCounts)
        random_var.append(var)
      print(m, ' done!')
      df_merian = df[df['Merian']== m]
      obs_var = statistics.variance(df_merian['matches'])
      sorted_random_var = sorted(random_var)
      percentile = 99.99
      value_at_percentile_cutoff = sorted_random_var[math.ceil(number_simulations * (percentile / 100))]
      if obs_var > value_at_percentile_cutoff:
        merian2sig[m] = "True"
      else:
          merian2sig[m] = "False"
      plt.figure()
      plot_title = str(m) + ' Number_markers=' + str(number_markers) + ' Obs var is less than ' +str(percentile)
      hist_plot = sns.histplot(data=random_var)
      hist_plot.axvline(x = obs_var, color="red")
      hist_plot.set(title=plot_title)
      plot_filename = str(m) + '_unequal_weights_' + weights + '.pdf'
      plt.savefig(plot_filename, format='pdf')
  return(merian2sig)

# %%
#equal_weights = simulate_variance(df, 'False')
unequal_weights = simulate_variance(df, 100000, 'True')
# %%
unequal_weights_df = pd.DataFrame(unequal_weights, index=[0]).transpose()
unequal_weights_df = unequal_weights_df.reset_index()
unequal_weights_df.columns = ['Merian', 'unequal_weight']
unequal_weights_df

#equal_weights_df = pd.DataFrame(equal_weights, index=[0]).transpose()
#equal_weights_df = equal_weights_df.reset_index()
#equal_weights_df.columns = ['Merian', 'equal_weight']
#equal_weights_df
#merged_df = pd.merge(equal_weights_df, unequal_weights_df,  on=['Merian'])
#merged_df['compare'] = "agree"
#merged_df.loc[(merged_df['equal_weight']!=merged_df['unequal_weight']), 'compare'] = 'disagree'
#merged_df['equal_weight'].value_counts()
# %%
unequal_weights_df['unequal_weight'].value_counts()
# %%
df, obs_var, random_var = simulate_variance(df, 'True') # True = unequal weights
print(mean(random_var))
print(median(random_var))
# %%
sorted_random_var = sorted(random_var)
percentile =99
value_at_percentile_cutoff = sorted_random_var[math.ceil(100000 * (percentile / 100))]
print(obs_var, value_at_percentile_cutoff)

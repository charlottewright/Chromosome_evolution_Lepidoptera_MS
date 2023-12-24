#%%
# #!/usr/bin/env python 
from os import spawnlp, supports_follow_symlinks
from sre_parse import parse_template
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import re
from collections import Counter
# %%

def parse_tables(mapped_rearrangements_file):
    all_rearrangenemnts = pd.read_csv(mapped_rearrangements_file, header=0, sep = "\t")
    all_fusions = all_rearrangenemnts[all_rearrangenemnts.Event == "fusion"] # remove fissions
    return(all_fusions)

def filter_fusions(all_fusions):
    filt_fusions = all_fusions[all_fusions.Node != "n6"] # LCA of Ditrysia
    filt_fusions = filt_fusions[filt_fusions.Node != "n1"] # LCA of trichoptera
    filt_fusions = filt_fusions[~filt_fusions.Node.str.contains("Limnephilus")] # trichoptera
    return(filt_fusions)

def reformat_ancient_fusion(filt_fusions): # this code replaces each instance of 'M17, M20' with 'M17.M20' such that the fusion is now made of n-1 Merian elements
    for index, row in filt_fusions.iterrows():
        Merians = row['Merians']
        Merians = re.sub(r'([^\s\w\.\,]|_)+', '', Merians)
        Merians = str(Merians).split(',')
        for Merian in Merians:
            if Merian == "M17":
                Merians.remove("M17")
        for i in range(len(Merians)):
            if Merians[i] == "M20":
                Merians[i] = "M17.M20"
                print('Updated merians:', Merians)
        Merians = re.sub(r'([^\s\w\.\,]|_)+', '', str(Merians))
        Merians = Merians.replace(' ','')
        filt_fusions.at[index,'Merians'] = Merians
    return(filt_fusions)

def make_merian_pairs_frequency_dict(filt_fusions):# Make a dict of the freq of each pair of Merians involved in a fusion
    freq_pairs_dict = {}
    for index, row in filt_fusions.iterrows():
        Merians = row['Merians']
        Merians = re.sub(r'([^\s\w\.\,]|_)+', '', Merians)
        Merians = str(Merians).split(',')
        if len(Merians) == 2: # only get fusions composed of two elements
            Merian_combos = [(Merians[i],Merians[j]) for i in range(len(Merians)) for j in range(i+1, len(Merians))]
            for i in Merian_combos: # add both orientations to freq_pairs_dict
                order_1 = sorted(i, reverse=True)
                order_2 = sorted(order_1)
                order_1, order_2 = tuple(order_1), tuple(order_2)
                if (order_1 in freq_pairs_dict): 
                    freq_pairs_dict[order_1] += 1
                else: 
                    freq_pairs_dict[order_1] = 1
                if (order_2 in freq_pairs_dict): 
                    freq_pairs_dict[order_2] += 1
                else: 
                    freq_pairs_dict[order_2] = 1
    return(freq_pairs_dict)

def make_individual_merian_frequency_dict(filt_fusions): # Make a dict containing the freq that each individual Merian is involved in a fusion
    freq_fusion = {}
    for index, row in filt_fusions.iterrows():
        Merians = row['Merians']
        Merians = re.sub(r'([^\s\w\.\,]|_)+', '', Merians)
        Merians = str(Merians)
        Merians = Merians.split(',')
        Merian_values = []
        if len(Merians) == 2: # only get fusions composed of 2 Merians 
            for Merian in Merians:
                # Need this line to remove whitespace
                Merian = Merian.strip()
                Merian_values.append(Merian)
            for k in Merians:
                if (k in freq_fusion): 
                    freq_fusion[k] += 1
                else: 
                    freq_fusion[k] = 1
    return(freq_fusion)

def convert_freq_pairs_dict_to_matrix(freq_pairs_dict):
    df = pd.DataFrame(freq_pairs_dict, index=[0]) # when passing scalar values for columns, you have to pass an index
    df = df.T # this makes rows the freqs and columns the Merian combos - so need to swap them around (transpose)
    df = df.reset_index() # move Merians from index to a column, this also splits each Merian into a seperate column
    df.columns = ['Merian_1', 'Merian_2', 'Freq'] # rename columns
    df = df.sort_values(by ='Merian_1' )
    df = df.sort_values(by ='Merian_2' )
    df_copy = df # group by Merian_1 and Merian_2, get the average
    df_copy = df_copy.groupby(['Merian_1', 'Merian_2']).mean()
    df_copy = df_copy.unstack(level=0) # unstack the indexes, and we’ll have our table
    df_copy = df_copy.fillna(0) # See lots of 'NaN' i.e. no fusion ever seen
    return(df_copy)

def plot_matrix(df_copy):
    fig, ax = plt.subplots(figsize=(11, 9))
    sb.heatmap(df_copy, cmap="Reds")  # plot heatmap
    ax.invert_yaxis()
    plt.show()
    return(plt)

def write_fusion_pairs_and_total_fusions_per_merian(freq_fusion, df_copy, prefix): # save output
    freq_fusion_df = pd.DataFrame(freq_fusion, index=[0])
    freq_fusion_df = freq_fusion_df.T
    freq_fusion_df = freq_fusion_df.reset_index()
    freq_fusion_df.columns = ['Merian', 'Freq']
    freq_fusion_df.to_csv(str(prefix + "_total_fusions_per_Merian.tsv"), header=True, index = False, sep = '\t')
    df_copy.to_csv(str(prefix + "_total_fusions_per_pair_of_Merians.tsv"), header=True, index = False, sep = '\t')
    return(freq_fusion_df)
#%%
# Define inputs
mapped_rearrangements_file = '/lustre/scratch123/tol/teams/blaxter/projects/lepidoptera_genomics/cw22/Leps_200/Analysis/LFSF/final_analysis/noComplex_281022/mapped_fusions_fissions_noComplex_181022.tsv'
prefix = '011122'

# Run functions
all_fusions = parse_tables(mapped_rearrangements_file,)
filt_fusions = filter_fusions(all_fusions)
filt_fusions = reformat_ancient_fusion(filt_fusions)
freq_pairs_dict = make_merian_pairs_frequency_dict(filt_fusions)
freq_fusion = make_individual_merian_frequency_dict(filt_fusions)
df_copy = convert_freq_pairs_dict_to_matrix(freq_pairs_dict)
plt = plot_matrix(df_copy)
freq_fusion_df = write_fusion_pairs_and_total_fusions_per_merian(freq_fusion, df_copy, prefix) # save output

# %%

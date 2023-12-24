
#!/usr/bin/env python3
#%%
import code
import random
import pandas as pd
from itertools import product
import seaborn as sns
import statistics
import math
import matplotlib.pyplot as plt
from ete3 import Tree
import sys
import os
#%%
# import tree
# with open(sys.argv[1], 'r') as treefile:
tree_file_path = 'r2_m5.newick.txt'
t = Tree(tree_file_path, format=1) # format=1 allows node labels to be read in
#print(t.get_ascii(show_internal=True))
#%%
all_species = list(t.get_leaf_names()) # get all tip names
species_to_remove = ["Micropterix_aruncella", "Limnephilus_lunatus", "Limnephilus_marmoratus","Limnephilus_rhombicus", "Glyphotaelius_pellucidus"]
ditrysia_species = [i for i in all_species if i not in species_to_remove]
print(len(ditrysia_species))
# count number species remaining
number_species =  len(ditrysia_species)
number_species_to_pick = int(number_species / 2) # make a whole number

if not os.path.exists('subsampled_lists'):
   os.makedirs('subsampled_lists')

if not os.path.exists('subsampled_trees'):
   os.makedirs('subsampled_trees')

#%%
n=1
for n in range(1,101,1):
    test_tree = t.copy()
    random_species = random.sample(ditrysia_species, number_species_to_pick) # pick a set of species at random
    for i in random_species:
        species = test_tree.search_nodes(name=i)[0]
        species.delete()
    ancestor = test_tree.get_common_ancestor("Limnephilus_lunatus", "Limnephilus_marmoratus","Limnephilus_rhombicus", "Glyphotaelius_pellucidus")
    test_tree.set_outgroup(ancestor)
    test_tree.write(format=0, outfile="subsampled_trees/subsampled_tree_" + str(n) + ".nw") # format=0 means no node labels - needed for syngraph
    all_species = list(test_tree.get_leaf_names())
    print(len(all_species))
    print(all_species[5])
    output_list = 'subsampled_lists/subsampled_list_' + str(n) + '.tsv'
    retained_species = set(all_species) - set(random_species) # want to keep those that have been retained in tree
    with open(output_list, 'w') as f:
        for i in retained_species:
            f.write(i +"\n")
        for i in species_to_remove:
            f.write(i + "\n")

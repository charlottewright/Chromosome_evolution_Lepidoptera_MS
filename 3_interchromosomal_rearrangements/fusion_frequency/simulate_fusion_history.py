#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: simulate_lepidopteran_fusions.py -s <INT> -f <INT> -t <STR> [-h]

  [Options]
    -s, --simulations <INT>                     Number of simulations
    -f, --fusions <INT>                         Number of fusions to uniformly sample across tree 
    -t, --tree <STR>                            Newick tree
    -h, --help                                  Show this message

"""

import sys
from docopt import docopt
import collections
import copy
import numpy
from scipy import stats
import random
import ete3
import math

random.seed(44)
numpy.random.seed(44)

def get_tree(tree_f):
	tree = ete3.Tree(str(tree_f))
	#print(tree.get_ascii())
	#print("")
	return tree

args = docopt(__doc__)
s_arg = int(args['--simulations'])
f_arg = int(args['--fusions'])
newick_path = args['--tree']

if s_arg < 10_000:
	sys.exit("[X] Increase simulations to at least 10_000")

tree = get_tree(newick_path) # read in newick tree

branches = [] # for every branch record whether it is external (1) or internal (0)
branch_lengths = [] # record branch lengths as this is the relative probability of a fusion landing there

trichopterans = {'Limnephilus_rhombicus', 'Limnephilus_marmoratus', 'Glyphotaelius_pellucidus', 
'Hydropsyche_tenuis', 'Limnephilus_lunatus'}

for node in tree.traverse(strategy="preorder"): # traverse the tree to collect the info mentioned above
	if not node.is_leaf() and not node.is_root():

		descendants = set()
		for descendant in node.iter_leaf_names():
			descendants.add(descendant)
		if descendants.intersection(trichopterans): # we ignore branches in the trichopteran part of the tree
			pass
		else:

			for child in node.get_children():

				branch_length = node.get_distance(child)
				branch_lengths.append(branch_length)

				if child.is_leaf():
					branches.append(1) # external branch
				else:
					branches.append(0) # internal branch

sim_results = []

for sim in range(0, s_arg): # for each sim do a random sampling and record how often fusions fall on external branches

	sample = random.choices(branches, weights=branch_lengths, k=f_arg)
	sim_results.append(sum(sample))

sim_results = sorted(sim_results)

# print some info about how many of the fusions we expect to find on external branches

for percentile in [0.01, 0.1, 1, 5, 10, 50, 90, 95, 99, 99.9, 99.99]:

	print("{}th percentile = {}".format(percentile, sim_results[math.ceil(s_arg * (percentile / 100))]))

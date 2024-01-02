# Section 3: Distribution of fusion and fission events across Lepidoptera

[**fusion_frequency**](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/tree/main/3_interchromosomal_rearrangements/fusion_frequency>)

Code and data to generate:
1. Boxplots of the distribution in proportional chromosome length within each Merian element. A TSV was plotted (prop_length_per_merian_per_species_141232.tsv).
2. Matrix of fusion events between pairs of Merian elements, where the shade of red indicates the total number of fusion events per Merian element. A TSV contained the frequency of fusions for each pair of Merian elements (freq_of_fusion_per_pair_of_merians_141232.tsv).
3. Barchart of the number of autosome-autosome and sex chromosome-autosome fusion events that each Merian element is involved in. A TSV contained the freuqency of each type of fusion events (number_autosome_autosome_vs_autosome_sex_fusions_per_Merian_141223.tsv).

These three TSVs were derived from two TSV files (011122_total_fusions_per_Merian.tsv and 011122_total_fusions_per_pair_of_Merians.tsv) which were generated using a python script (Make_heatmap_fusions.py). The three resulting TSV files were plotted using an R script (Fig3.R) to generate Figure 3a-c.

Also includes code required to compare the number of internal fusions to those on external nodes (simulate_fusion_history.py).

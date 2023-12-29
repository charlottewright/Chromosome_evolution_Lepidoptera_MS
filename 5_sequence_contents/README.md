# Section 5: Understanding biases in chromosomal fusions in Lepidoptera

[**sequence_patterns**](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/tree/main/5_sequence_contents/sequence_patterns>)

Code and data used to generate Figure 5a-c which contain scatterplots of repeat density, percent synteny and percent single copy orthologues per chromosome compared to proportional chromosome length for a set of species.
The data needed for these plots are the TSV files found in [Supplementary Table 10](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/blob/main/sup_tables/tableS10_chromosome_statistics.tsv>) and [Supplementary Table 8](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/blob/main/sup_tables/TableS8_species_statistics.tsv>) .
These TSV files were also used to generate scatterplots of GC, GC3 and coding density compared to proportional chromosome length for a set of species in Extended Data Figure 4. These TSV files are plotted using R scripts (Fig5A-Fig5C.R and functions_feature_patterns_plots.R and calculate_correlation_strengths_for_fig5.R).

[**repeat_distributions**](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/tree/main/5_sequence_contents/repeat_distributions>)

Code and data used to generate Extended Data Figure 5 which contains scatterplots of scaled repeat density for each major class of repeat. The TSV containing the repeat density per repeat class per Merian element scaled by the average repeat density of a given species (scaled_repeat_density_per_repeat_class_per_Merian_151223.tsv) was generated from a TSV (combined_repeat_density_per_family_per_chr.tsv) and plotted using an R script (plot_repeat_density_per_repeat_class_vs_chr_length_supplementary.R). 

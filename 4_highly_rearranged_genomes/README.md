# Section 4: Extensive rearrangemets in eight independent lineages

[**merian_paints**](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/tree/main/4_highly_rearranged_genomes/merian_paints>)

Code and data to generate the merian paints in Figure 4a, Figure 4b and Supplementary Figure 9.
Also requires the chromosome assignments to Merian elements from 
[Supplementary Table 10](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/blob/main/sup_tables/tableS10_chromosome_statistics.tsv>) and the assignmets of orthologues to Merian elements from [Supplementary Table 4](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/blob/main/sup_tables/TableS4_Merian_element_definitions.tsv>) and the 
[lepidopteran phylogeny](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/blob/main/1_genomes/phylogeny_210Leps_5Trichop.treefile>).
These data are plotted using the R scripts (Fig4.R, Pairwise_oxford_plots_supplementary.R, [functions_busco_painter.R](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/blob/main/2_merian_elements/functions_busco_painter.R>)
,[functions_tree_plots.R](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/blob/main/1_genomes/functions_tree_plots.R>). 

[**oxford_plots**](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/tree/main/4_highly_rearranged_genomes/oxford_plots>)

Code to generate the Oxford plots of gene order in *Lysandra coridon* compared to *Polyommatus icarus* in Supplementary Figure 9. This requires a TSV file containing the assignmet of orthologues to Merian elements from [Supplementary Table 4](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/blob/main/sup_tables/TableS4_Merian_element_definitions.tsv>) as well as two TSV files which contain the positions of each orthologue in each species (Lysandra_coridon.tsv and Polyommatus_icarus.tsv). These TSV files are plotted using an R scipt (Pairwise_oxford_plots_supplementary.R)


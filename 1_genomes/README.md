# Section 1: Over two hundred chromosomally-complete lepidopteran genomes

[**annotated_phylogeny**](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/tree/main/1_genomes/annotated_phylogeny>)

Contains code and files to generate the phylogeny coloured by the presence of interchromosomal rearrangements shown in Fig. 1a as well as the bar charts of genome size and chromosome number variation. The data needed for this plot are Supplementary Tables 1 and 6, which are plotted using R scripts (Fig1.R and functions_tree_plots.R).

[**lepidoptera_phylogeny**](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/tree/main/1_genomes/lepidoptera_phylogeny>)

Contains the newick file (phylogeny_210Leps_5Trichop.treefile) (derived from ASTRAL, with branch lengths from IQ-TREE) associated with the lepidoptera phylogeny shown in Figure 1A.

[**feature_distribution**](<https://github.com/charlottewright/Chromosome_evolution_Lepidoptera_MS/tree/main/1_genomes/feature_distribution>)

Contains the code to generate the distributions of GC, repeat density and coding density along chromosomes in 100 kb windows. BED files containing GC proportion, repeat density and coding density (combined_GC_per_chunk_100_per_chr.tsv.zip, combined_repeat_counts_per_100kb.tsv.zip, combined_GC_per_100kb.tsv.zip and combined_DivByThree.cds_density_counts_per_100kb.tsv.zip) are plotted using an R script (Plot_features_along_chr_supfig.R) which generates Exteded Data Figure 1a-d.


# Section 1: Over two hundred chromosomally-complete lepidopteran genomes

[**annotated_phylogeny**](<>)

Contains code and files to generate the phylogeny coloured by the presence of interchromosomal rearrangements shown in Fig. 1a as well as the bar charts of genome size and chromosome number variation. The data needed for this plot are Supplementary Tables 1 and 6, which are plotted using R scripts (Fig1.R and functions_tree_plots.R).

[**lepidoptera_phylogeny**](<>)

Contains the concatenated alignment (supermatrix.fa) and the resulting newick file (derived from ASTRAL, with branch lengths from IQ-TREE) associated with the lepidoptera phylogeny shown in Figure 1A.

[**feature_distribution**](<>)

Contains the code to generate the distributions of GC, repeat density and coding density along chromosomes in 100 kb windows. BED files containing GC proportion, repeat density and coding density are plotted using an R script (Plot_features_along_chr_supfig.R).


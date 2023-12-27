# Section 2: Thirty-two ancestral lepidopteran linkage groups

[**ancient_fusion**](<>)

Contains code and files to generate the annotated phylogeny and associated Oxford plots for Figure 2a.
The phylogeny can be found <>. A set of TSVs contain the locations of BUSCO genes in each species along with the Merian element assignment (busco_2_Merian_for_species_plotted_in_fig2A_151223.tsv and lep_odb10_busco_tables_all_spp.tar.gz). 
These files are plotted using R scripts (Fig2A.R and functions_oxford_plots.R).

[merian_paints](<>)

Contains code and data to generate the merian paints in Figure 2b. A set of TSVs contain the locations of BUSCO genes in each species along with the Merian element assignment (busco_2_Merian_for_species_plotted_in_fig2B_1223.tsv) derived from https://github.com/charlottewright/lep_busco_painter/. 
These files are plotted using R scripts (Fig2B.R and functions_busco_painter.R).

[compare_ALGs](<>)

Contains code and data to compare Merian elements to bilaterian linkage groups (BLGs). A TSV contains the assignment of each BUSCO gene to both Merian elements and bilaterian linkage groups (single_copy_orthogroups_with_Metazoa_ALGs_assigned_Merians_accurately.tsv).
The information from this TSV is summarised in two TSV files which are plotted as bar charts in Extended Data Figure 2a and 2c (Number_orthologues_per_Merian_also_assigned_to_BLGs_151223 and Number_orthologues_per_BLG_also_assigned_to_Merians_151223.tsv respectively). 
A heatmap is plotted using the TSV containing the proportion of each Merian element-defining orthologues assiged to each BLG in Extended Data Figure 2b
(Prop_Merian_othrologues_assigned_to_each_BLG_151223). The level of variation in assignment of orthologues to BLGs vs Merian elements was tested using a python script (simulate_BLG_vs_Merian_variation.py). These files are plotted using an R script (Compare_Merians_to_bilaterian_ALGs.R).

[fusion_frequency](<>)

Contains code and files to generate the scatter plot of proportional chromosome length compared to the number of fusion events shown in Extended Data Figure 3.
The number of fusions per Merian element were obtained using calculate_fusion_freq_per_merian.py and then summarised to generate the plot using an R script (make_source_data_for_chr_length_vs_num_fusions_per_Merian.R)
The resulting TSV contains the number of fusion events per Merian element and the average proportional length of each Merian element (Num_fusions_per_Merian_vs_average_prop_length.tsv).

### Figure 4
## 'Chromosome evolution in Lepidoptera'
## Wright et al., 2022
## cw22@sanger.ac.uk

library(phytools)
library(ggtree)
library(tidyverse)
library(scales)
library(ggpubr)
library(cowplot)

## set input path
setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/')

## load functions
source('functions_busco_painter.R')
source("functions_oxford_plots.R") 
source("functions_tree_plots.R")

## import data
assignments <- read.csv('../Sup_tables/Sup_table_all_chr.tsv', sep='\t')[,c(2,19)]
tree <- read.newick('../Data/busco2phylo/supermatrix_LG_G4_rooted_including_Hydropsyche.treefile')

## read in busco files
buscos_file_path <- '../Data/buscopaint/'
busco_dir <- "/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Data/AGORA/oxford_plot_analysis/"
Merian_assignments_ref <- read_busco('../Data/syngraph/Final_Merians_based_on_syn_n2_from_r2_m5_full_table.tsv')
shared_busco_tsv_suffix <- ".tsv"

## format assignments
colnames(assignments) <- c('query_chr', 'rearrangement_status')
Merian_assignments_ref <- select(Merian_assignments_ref, -c(start, end))
colnames(Merian_assignments_ref) <- c('busco_id', 'Merian')

## specify species to be plotted
target_species <- c("Aporia_crataegi","Pieris_brassicae","Pieris_napi","Pieris_rapae","Anthocharis_cardamines", 
             "Lysandra_coridon","Lysandra_bellargus","Polyommatus_icarus")

## make cladograms for subsets of species
p_spp <- c("Aporia_crataegi","Pieris_brassicae","Pieris_napi","Pieris_rapae","Anthocharis_cardamines")
l_spp <- c("Lysandra_coridon","Lysandra_bellargus","Polyommatus_icarus")

p_tree <- filter_tree(tree, p_spp) 
l_tree <- filter_tree(tree, l_spp) 

tip_label_size = 14
clade_label_fontsize = 13
info_size = 12
info_colour = "grey31"

pierini_tree <- ggtree(p_tree) +
  geom_text2(aes(subset=(node==1)), label = "italic('Anthocharis cardamines')", parse=TRUE, size=tip_label_size, hjust=0.4, vjust=1) +
  geom_text2(aes(subset=(node==1)), label = "'n = 30'", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=2.5, colour=info_colour) +
  geom_text2(aes(subset=(node==2)), label = "italic('Aporia crataegi')", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=1) +
  geom_text2(aes(subset=(node==2)), label = "'n = 25'", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=2.5, colour=info_colour) +
  geom_text2(aes(subset=(node==3)), label = "italic('Pieris brassicae')", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=1) +
  geom_text2(aes(subset=(node==3)), label = "'n = 15'", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=2.5, colour=info_colour) +
  geom_text2(aes(subset=(node==4)), label = "italic('Pieris napi')", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=1) +
  geom_text2(aes(subset=(node==4)), label = "'n = 25'", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=2.5,colour=info_colour) +
  geom_text2(aes(subset=(node==5)), label = "italic('Pieris rapae')", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=1) +
  geom_text2(aes(subset=(node==5)), label = "'n = 25'", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=2.5, colour=info_colour) + 
  geom_text2(aes(subset=(node==1)), label = "'1 fusion'", parse=TRUE, size=info_size, hjust=-0.1, vjust=-3, colour=info_colour) +
  geom_text2(aes(subset=(node==7)), label = "'18 fusions, 38 fissions'", parse=TRUE, size=info_size, hjust=-0.1, vjust=-0.1,colour=info_colour) +
  geom_text2(aes(subset=(node==2)), label = "' 39 fusions \n 13 fissions'", parse=TRUE, size=info_size, hjust=0, vjust=-2,colour=info_colour)  +
  geom_text2(aes(subset=(node==8)), label = "'35 fusions, 9 fissions'", parse=TRUE, size=info_size, hjust=-0.5, vjust=-0.1, colour=info_colour)  +
  geom_text2(aes(subset=(node==3)), label = "'10 fusions'", parse=TRUE, size=info_size, hjust=-0.1, vjust=-0.8,colour=info_colour) +
  geom_text2(aes(subset=(node==7)), label='Pierini', hjust = 0.5, vjust=2, color="black", size=clade_label_fontsize) +
  geom_point2(aes(subset=(node==7)), shape=21, size=10, fill="red")

lysandra_tree <- ggtree(l_tree) +
  geom_text2(aes(subset=(node==1)), label = "italic('Lysandra bellargus')", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=1)  +
  geom_text2(aes(subset=(node==1)), label = "'n = 45'", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=2.5,colour=info_colour) +
  geom_text2(aes(subset=(node==2)), label = "italic('Lysandra coridon')", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=1) +
  geom_text2(aes(subset=(node==2)), label = "'n = 90'", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=2.5,colour=info_colour) +
  geom_text2(aes(subset=(node==3)), label = "italic('Polyommatus icarus')", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=1) +
  geom_text2(aes(subset=(node==3)), label = "'n = 23'", parse=TRUE, size=tip_label_size, hjust=0.5, vjust=2.5,colour=info_colour) +
  geom_text2(aes(subset=(node==5)), label = "'15 fissions'", parse=TRUE, size=tip_label_size, hjust=0, vjust=-0.5, colour=info_colour) +
  geom_text2(aes(subset=(node==1)), label = "'6 fissions'", parse=TRUE, size=tip_label_size, hjust=1.2, vjust=-1, colour=info_colour) +
  geom_text2(aes(subset=(node==2)), label = "' 52 fissions\n 1 fusion'", parse=TRUE, size=tip_label_size, hjust=0, vjust=-1, colour=info_colour) +
  geom_text2(aes(subset=(node==3)), label = "'1 fusion'", parse=TRUE, size=tip_label_size, hjust=1.2, vjust=-0.5, colour=info_colour) +
  geom_point2(aes(subset=(node==5)), shape=21, size=10, fill="red")

upper_panel_p <- pierini_tree + layout_dendrogram() 
upper_panel_l <- lysandra_tree + layout_dendrogram() 

## make merian paint per species

## first edit assignments for L. coridon
# necessary as some chr have have fewer merian-defining loci than the detection threshold used (<= 17)
# manually edit to correct for this, using a lower threshold of 10 loci
# HG992076.1 (M27) - edit to split (15 buscos on HG992070.1)
# HG992123.1 - edit to split (14 buscos on HG992126.1)
# HG992070.1 - add in as M27 split ( 15 buscos map to M27)
# HG992106.1 - add in as M21 split (16 buscos map to M21)
# HG992126.1 - add in as M28 split (14 buscos map to M28)
# leaves two chr as genuinly intact - HG992057.1 (M31) & HG992077.1 (M29)

assignments <- assignments %>% mutate(rearrangement_status = ifelse(query_chr == "HG992076.1", "split", rearrangement_status))
assignments <- assignments %>% mutate(rearrangement_status = ifelse(query_chr == "HG992123.1", "split", rearrangement_status))
assignments <- assignments %>% 
  add_row(query_chr="HG992070.1", rearrangement_status = "split") %>%
  add_row(query_chr="HG992106.1", rearrangement_status = "split") %>%
  add_row(query_chr="HG992126.1", rearrangement_status = "split")

target_species_locations <- make_location_set(target_species, assignments)
xtext_size = 23

Aporia_plot <- get_busco_paint('Aporia_crataegi', target_species_locations) + theme(legend.position = "bottom") +rremove("xlab") + ggtitle(label="", subtitle= "") + theme(axis.text.x =element_text(size=xtext_size))
Pieris_bras_plot <- get_busco_paint('Pieris_brassicae', target_species_locations) + rremove("xlab") +  theme(legend.position = "none") + ggtitle(label="", subtitle= "") + theme(axis.text.x =element_text(size=xtext_size))
Pieris_napi_plot <- get_busco_paint('Pieris_napi', target_species_locations) + rremove("xlab") +  theme(legend.position = "none") + ggtitle(label="", subtitle= "") + theme(axis.text.x =element_text(size=xtext_size))
Pieris_rapae_plot <- get_busco_paint('Pieris_rapae', target_species_locations) + rremove("xlab") +  theme(legend.position = "none") + ggtitle(label="", subtitle= "") + theme(axis.text.x =element_text(size=xtext_size))
Anthocharis_plot <- get_busco_paint('Anthocharis_cardamines', target_species_locations) + rremove("xlab") +  theme(legend.position = "none") + ggtitle(label="", subtitle= "") + theme(axis.text.x =element_text(size=xtext_size))
Lysandra_cor_plot <- get_busco_paint('Lysandra_coridon', target_species_locations, 4) + rremove("xlab") +  theme(legend.position = "none") + ggtitle(label="", subtitle= "")+ theme(axis.text.x =element_text(size=xtext_size))
Lysandra_bel_plot <- get_busco_paint('Lysandra_bellargus', target_species_locations) + rremove("xlab") +  theme(legend.position = "none") + ggtitle(label="", subtitle= "")+ theme(axis.text.x =element_text(size=xtext_size))
Polyommatus_plot <- get_busco_paint('Polyommatus_icarus', target_species_locations) + rremove("xlab") +  theme(legend.position = "none") + ggtitle(label="", subtitle= "")+ theme(axis.text.x =element_text(size=xtext_size))

# testing
lower_panel_p <- Anthocharis_plot + Aporia_plot + Pieris_rapae_plot + Pieris_bras_plot + Pieris_napi_plot +  
  plot_layout(guides = "collect", nrow=1)
Lysandra_cor_plot + Lysandra_bel_plot +  Polyommatus_plot +  plot_layout(guides = "collect")
# end of testing
lower_panel_p <- plot_grid(Anthocharis_plot, Aporia_plot, Pieris_rapae_plot,Pieris_bras_plot, Pieris_napi_plot, nrow=1)
lower_panel_l <- plot_grid(Polyommatus_plot, Lysandra_bel_plot, Lysandra_cor_plot,nrow=1)

lower_panel_p <- annotate_figure(lower_panel_p, bottom = text_grob("Chr position (Mb)",size=30)) +  theme(plot.margin = unit(c(1, 0, 0, 0), "cm"))
lower_panel_l <- annotate_figure(lower_panel_l, bottom = text_grob("Chr position (Mb)",size=30)) +  theme(plot.margin = unit(c(1, 0, 0, 0), "cm"))

## combine cladograms and merian paints into one plot 
Pierini_panel <- (upper_panel_p / lower_panel_p) + plot_layout(heights = c(1, 4))
Lysandra_panel <- (upper_panel_l / lower_panel_l) + plot_layout(heights = c(1, 4))

## save figures
fig_height = 20
fig_width = 20

ggsave(plot=Pierini_panel, '../Figures/Main_text/Fig4A.pdf', device='pdf', width = fig_width, height = fig_width, dpi = 300, units = "in")
ggsave(plot=Lysandra_panel, '../Figures/Main_text/Fig4B.pdf', device='pdf', width = fig_width, height = fig_height, dpi = 300, units = "in")

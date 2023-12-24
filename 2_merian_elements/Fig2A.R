## Figure 2A
## 'Chromosome evolution in Lepidoptera'
## Wright et al., 2022
## cw22@sanger.ac.uk

library(tidyverse) 
library(patchwork)
library(ggplot2)
library(dplyr)
library(phytools)
library(ggtree)

## set input paths
setwd("/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/")

## load functions
source("/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera_MS/scripts/functions_oxford_plots.R") 

## import data
busco_dir <- "/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Data/BUSCOs/All/"
Merian_assignments_ref <- read_busco('../Data/syngraph/Final_Merians_based_on_syn_n2_from_r2_m5_full_table.tsv')
tree <- read.newick('../Data/busco2phylo/supermatrix_LG_G4_rooted_including_Hydropsyche.treefile')

## filter tree to just keep one representative of Trichoptera (L. marmomatus)
tips_to_drop <- c('Tinea_semifulvella', 'Limnephilus_rhombicus', 'Limnephilus_lunatus', 'Hydropsyche_tenuis')
for (i in 1:length(tips_to_drop)){
  tree <- drop.tip(tree, tips_to_drop[i])
}

## draw tree
p <- ggtree(tree) 
collapsed_p <-scaleClade(p, 216, scale = 0.003) %>% collapse(node=216, 'max', fill="white", colour="black")
tip_label_size = 10
clade_label_fontsize = 10
collapsed_p <- collapsed_p + 
  geom_text2(aes(subset=(node==1)), label = "italic('Glyphotaelius pellucidus')", parse=TRUE, size=tip_label_size, hjust=-.02) +
  geom_text2(aes(subset=(node==2)), label = "italic('Limnephilus marmoratus')", parse=TRUE, size=tip_label_size, hjust=-.02) +
  geom_text2(aes(subset=(node==3)), label = "italic('Micropterix aruncella')", parse=TRUE, size=tip_label_size,, hjust=-.02) +
  geom_text2(aes(subset=(node==4)), label = "italic('Tinea trinotella')", parse=TRUE, size=tip_label_size, hjust=-.02)

collapsed_p <- collapsed_p + xlim(-0.22,1.1) # this alters the spacing around plot

lep_colour <- '#3a0ca3'
ditrysian_colour <- '#4361ee'
summary_cladogram <- collapsed_p + 
  geom_text2(aes(subset=(node==216)), label = "bold('APODITRYSIA')", parse=TRUE, size=clade_label_fontsize, nudge_x=0.52, nudge_y=0.4) +
  geom_text2(aes(subset=(node==216)), label = "italic('E.g. Diarsia rubi')", parse=TRUE, size=tip_label_size, nudge_x=0.52, nudge_y=0.2) +
  geom_text2(aes(subset=(node==213)), label='TRICHOPTERA', nudge_y = -0.13, nudge_x=-0.2, color="black", fontface=2, size=clade_label_fontsize) +
  geom_text2(aes(subset=(node==214)), label='n=32',  nudge_x=0.1, color=lep_colour, fontface=2, size=10) +
  geom_text2(aes(subset=(node==214)), label='LEPIDOPTERA', nudge_y = 0.13, nudge_x=-0.2, color=lep_colour, fontface=2, size=clade_label_fontsize) +
  geom_text2(aes(subset=(node==215)), label='DITRYSIA', nudge_y = -0.2, nudge_x=-0.15, color=ditrysian_colour, fontface=2, size=clade_label_fontsize) +
  geom_text2(aes(subset=(node==215)), label='n=31',  nudge_x=0.1, color=ditrysian_colour, fontface=2, size=10) +
  geom_text2(aes(subset=(node==215)), label='M17-M20 fusion', nudge_y = 0.13, nudge_x=-0.2, color="black", fontface=1, size=clade_label_fontsize) +
  geom_point2(aes(subset=(node==214)), shape=21, size=10, fill=lep_colour) +
  geom_point2(aes(subset=(node==215)), shape=21, size=10, fill=ditrysian_colour) +
  geom_text2(aes(subset=(node==4)), label = 'n = 30', size=tip_label_size, nudge_y=-0.18, nudge_x=0.14, colour="Dimgrey") +
  geom_text2(aes(subset=(node==3)), label = 'n = 31', size=tip_label_size, nudge_y=-0.18, nudge_x=0.25, colour="Dimgrey") +
  geom_text2(aes(subset=(node==2)), label = 'n = 30', size=tip_label_size, nudge_y=-0.18, nudge_x=0.28, colour="Dimgrey") +
  geom_text2(aes(subset=(node==1)), label = 'n = 30', size=tip_label_size, nudge_y=-0.18, nudge_x=0.28, colour="Dimgrey") +
  geom_text2(aes(subset=(node==216)), label = 'n = 31', size=tip_label_size, nudge_y=0.05, nudge_x=0.52, colour="Dimgrey") +
  geom_text2(aes(subset=(node==3)), label = 'M11-MZ fusion', fontface=1, size=clade_label_fontsize, nudge_y = 0.13, nudge_x=-0.2, colour="black") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  ggtitle("A") + theme(title = element_text(face="bold", size=18))

## make oxford plots
shared_busco_tsv_suffix <- ".tsv"
Merian_assignments_ref <- select(Merian_assignments_ref, -c(start, end))
colnames(Merian_assignments_ref) <- c('busco_id', 'Merian')

buscofiles <- list.files(busco_dir, full.names = T,
                         pattern = shared_busco_tsv_suffix)
names(buscofiles) <- make.names(sub(paste0(".+//(.+)", 
                                           shared_busco_tsv_suffix), "\\1", 
                                    buscofiles))
buscos <- map_df(buscofiles, read_busco,
                 .id = "assembly")  
buscos <- merge(buscos, Merian_assignments_ref, by='busco_id')

buscos <- buscos %>%
  mutate(short_species_id = sub(busco_dir, "\\1", assembly),  
         multispecies_sequence = paste0(short_species_id, "_", Sequence))

filt_buscos <- filter_buscos(buscos, multispecies_sequence)
ancient_fusion_buscos <- filt_buscos[(filt_buscos$Merian == 'M17') | (filt_buscos$Merian == 'M20'),]

ancient_fusion_buscos <- ancient_fusion_buscos %>% group_by(multispecies_sequence, Merian) %>% mutate(count = n()) %>% ungroup()
ancient_fusion_buscos <- ancient_fusion_buscos %>% filter(count >5)

subset_spp <- c('Melitaea_cinxia', 'Diarsia_rubi', 'Tinea_trinotella', 'Micropterix_aruncella', 'Limnephilus_rhombicus', 'Glyphotaelius_pellucidus')
subset_spp_buscos <- ancient_fusion_buscos[ancient_fusion_buscos$short_species_id %in% subset_spp, ]

## draw oxford plots
custom_colors <- c('#f72585', '#7209b7') %>% rep(16)

a <- multispecies_oxford_plot(ancient_fusion_buscos, "Melitaea_cinxia", "Diarsia_rubi", custom_colors) + xlab(NULL) + ylab(NULL) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
b <- multispecies_oxford_plot(ancient_fusion_buscos, "Melitaea_cinxia", "Tinea_trinotella", custom_colors) + xlab(NULL) + ylab(NULL) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
c <- multispecies_oxford_plot(ancient_fusion_buscos, "Melitaea_cinxia", "Micropterix_aruncella", custom_colors) + xlab(NULL)+ ylab(NULL) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
d <- multispecies_oxford_plot(ancient_fusion_buscos, "Melitaea_cinxia", "Limnephilus_rhombicus", custom_colors)  + xlab(NULL)+ ylab(NULL) +  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
e <- multispecies_oxford_plot(ancient_fusion_buscos, "Melitaea_cinxia", "Glyphotaelius_pellucidus", custom_colors)  + xlab(NULL) + ylab(NULL) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
                                                                                                                                                      
e <- e + theme(axis.title.x=element_text(size=25)) +
  xlab("Melitaea cinxia") # add x-axis title

# plot cladogram and oxford plots together
oxford_plots <- a + b + c + d + e + plot_layout(ncol = 1)
Fig2A <- summary_cladogram + oxford_plots + plot_layout(widths = c(3,1))

## save Fig2A
ggsave(plot=Fig2A, '../Figures/Fig2A', device='pdf', width = 15, height = 15, dpi = 300, units = "in")


# save busco table as an output tsv to add to github repo
write.table(buscos, file = "../../Chromosome_evolution_Lepidoptera_MS/data/busco_2_Merian_for_all_217_species_151223.tsv", row.names=FALSE, sep="\t", quote = FALSE)
# save subset_spp_buscos as an output tsv to add to source_data
write.table(subset_spp_buscos, file = "../../Chromosome_evolution_Lepidoptera_MS/data/busco_2_Merian_for_species_plotted_in_fig2A_151223.tsv", row.names=FALSE, sep="\t", quote = FALSE)


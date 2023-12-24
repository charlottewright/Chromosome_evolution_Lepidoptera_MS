## Figure 5A-5C
## 'Chromosome evolution in Lepidoptera'
## Wright et al., 2022
## cw22@sanger.ac.uk

library(patchwork)
library(ggplot2)f
library(dplyr)
library(tidyr)

## set input paths
setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/')

## load functions to assist with plotting
source('functions_feature_patterns_plots.R')

## import data
correlations_table <- read.csv('../Sup_tables/Sup_table_correlations_strength.tsv', sep='\t')
features <- read.csv('../Sup_tables/Sup_table_all_chr.tsv', sep='\t')

## format data
features$status[features$assigned_merian == "M17+M20"] <- "ancestral_Ditrysia"
ancestral_chromosomes <- features %>% filter(status %in% c("ancestral","ancestral_Ditrysia")) # only keep ancestral merians
ancestral_autosomes <- ancestral_chromosomes %>% filter(assigned_merian!= "MZ")
chr_to_keep <- ancestral_autosomes %>% group_by(species) %>%  dplyr::mutate(num_ancestral = n()) %>% ungroup() %>% filter(num_ancestral>=10) 
# note that the resulting df causes issues with ggplot2, instead use it to filter the original df 
chr_to_keep <- chr_to_keep$chr
ancestral_features_filt_noMZ <- ancestral_autosomes[ancestral_autosomes$chr %in% chr_to_keep,]

# combine df with correlation values
significant_species_repeat <- (correlations_table %>% filter(rep_pvalue < 0.05))$species
significant_species_cds <- (correlations_table %>% filter(cds_pvalue < 0.05))$species
significant_species_synteny <- (correlations_table %>% filter(synteny_pvalue < 0.05))$species
significant_species_orthologs <- (correlations_table %>% filter(ortho_pvalue < 0.05))$species
significant_species_gc3 <- (correlations_table %>% filter(gc3_pvalue < 0.05))$species
significant_species_gc <- (correlations_table %>% filter(gc_pvalue < 0.05))$species

ancestral_features_filt_noMZ$repeat_sig <- "FALSE"
ancestral_features_filt_noMZ$repeat_sig[ancestral_features_filt_noMZ$species %in% significant_species_repeat] <- "TRUE"

ancestral_features_filt_noMZ$cds_sig <- "FALSE"
ancestral_features_filt_noMZ$cds_sig[ancestral_features_filt_noMZ$species %in% significant_species_cds] <- "TRUE"

ancestral_features_filt_noMZ$synteny_sig <- "FALSE"
ancestral_features_filt_noMZ$synteny_sig[ancestral_features_filt_noMZ$species %in% significant_species_synteny] <- "TRUE"

ancestral_features_filt_noMZ$orthologs_sig <- "FALSE"
ancestral_features_filt_noMZ$orthologs_sig[ancestral_features_filt_noMZ$species %in% significant_species_orthologs] <- "TRUE"

ancestral_features_filt_noMZ$gc3_sig <- "FALSE"
ancestral_features_filt_noMZ$gc3_sig[ancestral_features_filt_noMZ$species %in% significant_species_gc3] <- "TRUE"

ancestral_features_filt_noMZ$gc_sig <- "FALSE"
ancestral_features_filt_noMZ$gc_sig[ancestral_features_filt_noMZ$species %in% significant_species_gc] <- "TRUE"

sig_levels = c("TRUE","FALSE")

## convert to factor for consistent colour matching
ancestral_features_filt_noMZ$repeat_sig = factor(ancestral_features_filt_noMZ$repeat_sig, levels=sig_levels)
ancestral_features_filt_noMZ$orthologs_sig = factor(ancestral_features_filt_noMZ$orthologs_sig, levels=sig_levels)
ancestral_features_filt_noMZ$cds_sig = factor(ancestral_features_filt_noMZ$cds_sig, levels=sig_levels)
ancestral_features_filt_noMZ$orthologs_sig = factor(ancestral_features_filt_noMZ$orthologs_sig, levels=sig_levels)
ancestral_features_filt_noMZ$synteny_sig = factor(ancestral_features_filt_noMZ$synteny_sig, levels=sig_levels)
ancestral_features_filt_noMZ$gc3_sig = factor(ancestral_features_filt_noMZ$gc3_sig, levels=sig_levels)
ancestral_features_filt_noMZ$gc_sig = factor(ancestral_features_filt_noMZ$gc_sig, levels=sig_levels)

## make plots
Bombycoidea_plot <- plot_superfamily("Bombycoidea", ancestral_features_filt_noMZ, "repeat_density","repeat_sig") +
  ggtitle("Repeat density") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = "0.5"))
Geometroidea_plot <- plot_superfamily("Geometroidea", ancestral_features_filt_noMZ,"repeat_density","repeat_sig")
Pyraloidea_plot <- plot_superfamily("Pyraloidea", ancestral_features_filt_noMZ,"repeat_density","repeat_sig")
Noctuoidea_plot <- plot_superfamily("Noctuoidea", ancestral_features_filt_noMZ, "repeat_density","repeat_sig")
Gelechioidea_plot <- plot_superfamily("Gelechioidea", ancestral_features_filt_noMZ,"repeat_density","repeat_sig")
Papilionoidea_plot <- plot_superfamily("Papilionoidea", ancestral_features_filt_noMZ,"repeat_density","repeat_sig")
Sesioidea_plot <- plot_superfamily("Sesioidea", ancestral_features_filt_noMZ,"repeat_density","repeat_sig")
Tortricoidea_plot_with_axes <- plot_superfamily_with_axes("Tortricoidea", ancestral_features_filt_noMZ,"repeat_density","repeat_sig")

Bombycoidea_cds_plot <- plot_superfamily("Bombycoidea", ancestral_features_filt_noMZ,"cds_density", "cds_sig") +
  ggtitle("Coding density") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = "0.5"))
Pyraloidea_cds_plot <- plot_superfamily("Pyraloidea", ancestral_features_filt_noMZ,"cds_density", "cds_sig")
Noctuoidea_cds_plot <- plot_superfamily("Noctuoidea", ancestral_features_filt_noMZ,"cds_density", "cds_sig")
Gelechioidea_cds_plot <- plot_superfamily("Gelechioidea", ancestral_features_filt_noMZ,"cds_density", "cds_sig")
Geometroidea_cds_plot <- plot_superfamily("Geometroidea", ancestral_features_filt_noMZ,"cds_density", "cds_sig")
Papilionoidea_cds_plot <- plot_superfamily("Papilionoidea", ancestral_features_filt_noMZ,"cds_density", "cds_sig")
Sesioidea_cds_plot <- plot_superfamily("Sesioidea", ancestral_features_filt_noMZ,"cds_density", "cds_sig")
Tortricoidea_cds_plot_with_axies <- plot_superfamily_with_axes("Tortricoidea", ancestral_features_filt_noMZ,"cds_density", "cds_sig")

Bombycoidea_synteny_plot <- plot_superfamily("Bombycoidea", ancestral_features_filt_noMZ,"synteny", "synteny_sig") +
  ggtitle("Synteny") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = "0.5"))
Pyraloidea_synteny_plot <- plot_superfamily("Pyraloidea", ancestral_features_filt_noMZ, "synteny","synteny_sig")
Noctuoidea_synteny_plot <- plot_superfamily("Noctuoidea", ancestral_features_filt_noMZ, "synteny","synteny_sig")
Gelechioidea_synteny_plot <- plot_superfamily("Gelechioidea", ancestral_features_filt_noMZ, "synteny","synteny_sig")
Geometroidea_synteny_plot <- plot_superfamily("Geometroidea", ancestral_features_filt_noMZ, "synteny","synteny_sig")
Papilionoidea_synteny_plot <- plot_superfamily("Papilionoidea", ancestral_features_filt_noMZ, "synteny","synteny_sig")
Sesioidea_synteny_plot <- plot_superfamily("Sesioidea", ancestral_features_filt_noMZ, "synteny","synteny_sig")
Tortricoidea_synteny_plot_with_axies <- plot_superfamily_with_axes("Tortricoidea", ancestral_features_filt_noMZ, "synteny","synteny_sig")

Bombycoidea_orthologs_plot <- plot_superfamily("Bombycoidea", ancestral_features_filt_noMZ, "prop_single_copy", "orthologs_sig") +
  ggtitle("Single copy orthologs") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = "0.5"))
Pyraloidea_orthologs_plot  <- plot_superfamily("Pyraloidea", ancestral_features_filt_noMZ, "prop_single_copy", "orthologs_sig")
Noctuoidea_orthologs_plot  <- plot_superfamily("Noctuoidea", ancestral_features_filt_noMZ, "prop_single_copy", "orthologs_sig")
Gelechioidea_orthologs_plot  <- plot_superfamily("Gelechioidea", ancestral_features_filt_noMZ, "prop_single_copy", "orthologs_sig")
Geometroidea_orthologs_plot  <- plot_superfamily("Geometroidea", ancestral_features_filt_noMZ, "prop_single_copy", "orthologs_sig")
Papilionoidea_orthologs_plot <- plot_superfamily("Papilionoidea", ancestral_features_filt_noMZ, "prop_single_copy", "orthologs_sig")
Sesioidea_orthologs_plot  <- plot_superfamily("Sesioidea", ancestral_features_filt_noMZ, "prop_single_copy", "orthologs_sig")
Tortricoidea_orthologs_plot_with_axies <- plot_superfamily_with_axes("Tortricoidea", ancestral_features_filt_noMZ, "prop_single_copy", "orthologs_sig")

Bombycoidea_gc3_plot <- plot_superfamily("Bombycoidea", ancestral_features_filt_noMZ, "gc3_per", "gc3_sig") +
  ggtitle("GC3") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = "0.5"))
Pyraloidea_gc3_plot  <- plot_superfamily("Pyraloidea", ancestral_features_filt_noMZ, "gc3_per", "gc3_sig")
Noctuoidea_gc3_plot  <- plot_superfamily("Noctuoidea", ancestral_features_filt_noMZ, "gc3_per", "gc3_sig")
Gelechioidea_gc3_plot  <- plot_superfamily("Gelechioidea", ancestral_features_filt_noMZ, "gc3_per", "gc3_sig")
Geometroidea_gc3_plot  <- plot_superfamily("Geometroidea", ancestral_features_filt_noMZ, "gc3_per", "gc3_sig")
Papilionoidea_gc3_plot <- plot_superfamily("Papilionoidea", ancestral_features_filt_noMZ, "gc3_per", "gc3_sig")
Sesioidea_gc3_plot  <- plot_superfamily("Sesioidea", ancestral_features_filt_noMZ, "gc3_per", "gc3_sig")
Tortricoidea_gc3_plot_with_axies <- plot_superfamily_with_axes("Tortricoidea", ancestral_features_filt_noMZ, "gc3_per", "gc3_sig")

Bombycoidea_gc_plot <- plot_superfamily("Bombycoidea", ancestral_features_filt_noMZ, "gc", "gc_sig") +
  ggtitle("GC") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = "0.5"))
Pyraloidea_gc_plot <- plot_superfamily("Pyraloidea", ancestral_features_filt_noMZ, "gc", "gc_sig")
Noctuoidea_gc_plot <- plot_superfamily("Noctuoidea", ancestral_features_filt_noMZ, "gc", "gc_sig")
Gelechioidea_gc_plot <- plot_superfamily("Gelechioidea", ancestral_features_filt_noMZ, "gc", "gc_sig")
Geometroidea_gc_plot <- plot_superfamily("Geometroidea", ancestral_features_filt_noMZ, "gc", "gc_sig")
Papilionoidea_gc_plot <- plot_superfamily("Papilionoidea", ancestral_features_filt_noMZ, "gc", "gc_sig")
Sesioidea_gc_plot  <- plot_superfamily("Sesioidea", ancestral_features_filt_noMZ, "gc", "gc_sig")
Tortricoidea_gc_plot_with_axies <- plot_superfamily_with_axes("Tortricoidea", ancestral_features_filt_noMZ, "gc", "gc_sig")
Papilionoidea_gc_plot <- Papilionoidea_gc_plot + ylab("Feature %") + theme(axis.title = element_text(size=15))

fig5_pt1 <- wrap_plots(Bombycoidea_plot,Bombycoidea_synteny_plot, Bombycoidea_orthologs_plot,ncol=3) +
  plot_annotation(tag_levels = 'A')                        
                           
fig5_pt2 <- wrap_plots(Geometroidea_plot,Geometroidea_synteny_plot, Geometroidea_orthologs_plot,
                       Noctuoidea_plot,Noctuoidea_synteny_plot, Noctuoidea_orthologs_plot,
                       Pyraloidea_plot,Pyraloidea_synteny_plot,Pyraloidea_orthologs_plot,
                       Papilionoidea_plot,Papilionoidea_synteny_plot, Papilionoidea_orthologs_plot,
                       Gelechioidea_plot,Gelechioidea_synteny_plot,Gelechioidea_orthologs_plot,
                       Sesioidea_plot,Sesioidea_synteny_plot,Sesioidea_orthologs_plot,
                       Tortricoidea_plot_with_axes,Tortricoidea_synteny_plot_with_axies,Tortricoidea_orthologs_plot_with_axies, ncol=3)
  
fig5A-C <- ggarrange(fig5_pt1,fig5_pt2, ncol=1, heights=c(1,5))                 

ggsave(plot=fig5A, path = "../Figures/", 
       filename = "Fig5A-Fig5C.pdf", width=8, height=10)

## make corresponding supplementary figure - GC, GC3 and coding density
SupFig_pt1 <- wrap_plots(Bombycoidea_gc_plot, Bombycoidea_gc3_plot,Bombycoidea_cds_plot,ncol=3) +
  plot_annotation(tag_levels = 'A')                         

SupFig_pt2 <- wrap_plots(Geometroidea_gc_plot,Geometroidea_gc3_plot,Geometroidea_cds_plot,
                         Noctuoidea_gc_plot,Noctuoidea_gc3_plot,Noctuoidea_cds_plot,
                         Pyraloidea_gc_plot,Pyraloidea_gc3_plot,Pyraloidea_cds_plot,
                         Papilionoidea_gc_plot,Papilionoidea_gc3_plot,Papilionoidea_cds_plot,
                         Gelechioidea_gc_plot,Gelechioidea_gc3_plot,Gelechioidea_cds_plot,
                         Sesioidea_gc_plot,Sesioidea_gc3_plot,Sesioidea_cds_plot,
                         Tortricoidea_gc_plot_with_axies,Tortricoidea_gc3_plot_with_axies,Tortricoidea_cds_plot_with_axies, ncol=3)

SupFig_plot <- ggarrange(SupFig_pt1,SupFig_pt2, ncol=1, heights=c(1,5))                 

ggsave(plot=SupFig_plot, path = "../Figures/Supplementary/", 
       filename = "SF_GC_and_CDS_analysis_spearmans_smooth_with_dots.png", width=8, height=10)

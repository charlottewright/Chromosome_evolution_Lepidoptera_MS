library(dplyr)
library(ggplot2)
library(patchwork)

# Aim: Plot repeat density across fused chr and compare to non-fused chr for sup fig
setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/')

# functions
plot_fused_vs_ancestral <- function(df, merian_info_per_segment, spp1, spp2, merian1, merian2){
  merian1_df <- merian_info_per_segment %>% filter(species==spp1) %>% filter(merian_segment == merian1)
  merian1_fusion_coord <- merian1_df$prop_stop # get pos of first merian ending 
  merian2_df <- merian_info_per_segment %>% filter(species==spp1) %>% filter(merian_segment == merian2)
  merian2_fusion_coord <- merian2_df$prop_stop # get pos of first merian ending
  fusion_coord <- min(merian1_fusion_coord, merian2_fusion_coord) # use whichever is smaller i.e. which merian is first along chr
  
  
  spp1_df <- df[df$species == spp1,] 
  spp2_df <- df[df$species == spp2,] 
  spp1_df <- spp1_df[spp1_df$assigned_merian == paste(merian1, ',', merian2, sep = ''),]
  spp2_df <- spp2_df[spp2_df$assigned_merian %in% c(merian1, merian2),]
  spp1_df$plot_colour <- "green" # sanity check - shouldn't see any green
  chr_lengths <- spp2_df %>% group_by(chr) %>% summarise(Value = max(prop_midpos))
  ratio_plots <- chr_lengths[1,2] / chr_lengths[2,2]
  
  if (merian1_fusion_coord < merian2_fusion_coord) { # if merian1 is before merian2 along the chr
    spp1_df$plot_colour[spp1_df$prop_midpos < merian2_fusion_coord] <- merian2 
    spp1_df$plot_colour[spp1_df$prop_midpos < merian1_fusion_coord] <- merian1
  } else { # else if merian2 is before merian1
    spp1_df$plot_colour[spp1_df$prop_midpos < merian1_fusion_coord] <- merian1
    spp1_df$plot_colour[spp1_df$prop_midpos < merian2_fusion_coord] <- merian2 
  }
  
  spp2_df$plot_colour <- spp2_df$assigned_merian
  
  merian1_max <- max(spp1_df$prop_midpos, na.rm=TRUE)
  merian2_max <- max(spp2_df$prop_midpos, na.rm = TRUE)
  xlim_max <- max(merian1_max, merian2_max)
  
  
  spp1_plot <- ggplot(data=spp1_df, aes(x=prop_midpos, y=repeat_density)) + 
    geom_point(size=0.4, alpha=0.2, aes(colour=plot_colour)) + facet_wrap(~chr,  scales = "free_x") +
    geom_smooth(method ="loess", span=0.3, colour="black") +ylim(15, 80) + theme_bw() + my_labels +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    scale_color_manual(values=c("blue","red"), name="Merian element") + 
    theme(axis.title = element_text(size=8)) + # draw line at fusion coord
    geom_vline(xintercept = fusion_coord, linetype="dotted",  
               color = "blue", size=1.5)  + ggtitle(spp1) + theme(plot.title = element_text(size=9, face = "italic"))
  
  spp2_df$plot_colour <- factor(spp2_df$plot_colour, levels=c("blue", merian1, merian2))
  spp2_chr1_df <- spp2_df[spp2_df$assigned_merian == merian1,]
  spp2_chr2_df <- spp2_df[spp2_df$assigned_merian == merian2,]
  spp2_chr1_max <- max(spp2_chr1_df$prop_midpos)
  spp2_chr2_max <- max(spp2_chr2_df$prop_midpos)
  
  spp2_chr1_plot <- ggplot(data=spp2_chr1_df, aes(x=prop_midpos, y=repeat_density)) + 
    geom_point(size=0.4, alpha=0.2, aes(colour=plot_colour)) +
    geom_smooth(method ="loess", span=0.6, colour="black") + theme_bw() + ylim(15, 80) + my_labels +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +     
    scale_color_manual(values=c("blue", "red"), breaks=c(merian1, merian2), name="Merian element") + 
    theme(axis.title = element_text(size=8),
          plot.title = element_text(size=9, face = "italic"),
          legend.position = "none") + ylab("") 
  
  spp2_chr2_plot <- ggplot(data=spp2_chr2_df, aes(x=prop_midpos, y=repeat_density)) + 
    geom_point(size=0.4, alpha=0.2, aes(colour=plot_colour)) +
    geom_smooth(method ="loess", span=0.6, colour="black") + theme_bw() + ylim(15, 80) + my_labels +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) + 
    scale_color_manual(values=c("blue", "red"), breaks =c(merian1, merian2),name="Merian element") + 
    theme(axis.title = element_text(size=8)) +
    theme(plot.title = element_text(size=9, face = "italic"),
          legend.position = "none") + ylab("") 
  
  
  if (spp2_chr1_max > spp2_chr2_max) { # if merian1 is before merian2 along the chr
    spp2_chr2_plot <- spp2_chr2_plot + theme(axis.text.y= element_blank(),
                                             axis.ticks.y = element_blank())
    unfused_plot <- spp2_chr1_plot + ggtitle(spp2) +spp2_chr2_plot + plot_layout(ncol = 2,  widths=c(ratio_plots, 1)) 
  } else { # else if merian2 is before merian1
    spp2_chr1_plot <- spp2_chr1_plot + theme(axis.text.y= element_blank(),
                                             axis.ticks.y = element_blank())
    unfused_plot <- spp2_chr2_plot + ggtitle(spp2) +spp2_chr1_plot + plot_layout(ncol = 2,  widths=c(ratio_plots,1)) 
    
  }
  plot <- spp1_plot / unfused_plot +  plot_layout(ncol = 1,  heights=c(1,1)) 
  return(plot)
}

plot_fused_only <- function(df, merian_info_per_segment, spp1, merian1, merian2){
  merian1_df <- merian_info_per_segment %>% filter(species==spp1) %>% filter(merian_segment == merian1)
  merian1_fusion_coord <- merian1_df$prop_stop # get pos of first merian ending 
  merian2_df <- merian_info_per_segment %>% filter(species==spp1) %>% filter(merian_segment == merian2)
  merian2_fusion_coord <- merian2_df$prop_stop # get pos of first merian ending
  fusion_coord <- min(merian1_fusion_coord, merian2_fusion_coord) # use whichever is smaller i.e. which merian is first along chr
  
  
  spp1_df <- df[df$species == spp1,] 
  spp1_df <- spp1_df[spp1_df$assigned_merian == paste(merian1, ',', merian2, sep = ''),]
  spp1_df$plot_colour <- "green" # sanity check - shouldn't see any green
  
  if (merian1_fusion_coord < merian2_fusion_coord) { # if merian1 is before merian2 along the chr
    spp1_df$plot_colour[spp1_df$prop_midpos < merian2_fusion_coord] <- merian2 
    spp1_df$plot_colour[spp1_df$prop_midpos < merian1_fusion_coord] <- merian1
  } else { # else if merian2 is before merian1
    spp1_df$plot_colour[spp1_df$prop_midpos < merian1_fusion_coord] <- merian1
    spp1_df$plot_colour[spp1_df$prop_midpos < merian2_fusion_coord] <- merian2 
  }
  
  merian1_max <- max(spp1_df$prop_midpos, na.rm=TRUE)
  xlim_max <- max(merian1_max)
  
  spp1_plot <- ggplot(data=spp1_df, aes(x=prop_midpos, y=repeat_density)) + 
    geom_point(size=1, alpha=0.5, aes(colour=plot_colour)) +
    geom_smooth(method ="loess", span=0.6, colour="black", se=F) +ylim(15, 80) + theme_bw() + my_labels +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    scale_color_manual(values=c("blue","red"), name="Merian element") + 
    theme(axis.title = element_text(size=8)) + # draw line at fusion coord
    geom_vline(xintercept = fusion_coord, linetype="dotted",  
               color = "blue", size=2) + xlim(0,xlim_max) + 
    theme(plot.title = element_text(size=9, face = "italic")) +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  return(spp1_plot)
}
all_windows <- read.csv('../Data/combined/per_100kb/combined_repeat_counts_per_100kb.tsv', sep='\t', header=FALSE)

all_chromosomes <- read.csv('../Sup_tables/Sup_table_all_assigned_chr_no_complex.tsv', sep='\t')
all_chromosomes <- all_chromosomes %>% select(assigned_merian, species, chr,genome_size)
merian_info_per_segment <- read.csv('../Data/lep_fusion_split2/final_analysis/Fusion_segment_merian_info_adjusted.tsv', sep='\t', header=FALSE)
colnames(merian_info_per_segment) <- c('species', 'chr', 'merian_segment','start', 'stop')

# Clean data
all_windows <- all_windows[,c(1,2,3,4,5)]
colnames(all_windows) <- c('species', 'chr', 'start', 'stop', 'repeat_density')
all_windows$midpos <- (all_windows$stop + all_windows$start)/2
all_windows$midpos <- all_windows$midpos/1000000
all_windows$repeat_density <- all_windows$repeat_density*100
all_windows <- all_windows[,c(1,2,5,6)]


# Make one df
all_windows <- merge(all_windows, all_chromosomes, by=c("species", "chr"))
all_windows$genome_size <- all_windows$genome_size / 1000000 # convert to Mb
all_windows$prop_midpos <- (all_windows$midpos / all_windows$genome_size)*100 # convert midpos to prop length
merian_info_per_segment <- all_chromosomes %>% 
  select(species, genome_size) %>% 
  unique() %>%
  merge(merian_info_per_segment)


merian_info_per_segment$prop_stop <- (merian_info_per_segment$stop / merian_info_per_segment$genome_size)*100
my_labels <- list(ylab("Repeat density (%)"),
                  xlab("Proportional chr length (%)"))



# Recent fusions
a <- plot_fused_vs_ancestral(all_windows, merian_info_per_segment, "Agrochola_circellaris", "Agrochola_macilenta", "M30", "M5")
b <- plot_fused_vs_ancestral(all_windows, merian_info_per_segment, "Dendrolimus_kikuchii", "Dendrolimus_punctatus", "M31", "MZ")
c <- plot_fused_vs_ancestral(all_windows, merian_info_per_segment, "Agonopterix_subpropinquella", "Agonopterix_arenella", "M19", "M21")

recent <-  plot_spacer() + a + b +  c + plot_layout(ncol = 1, heights=c(0.1,1,1,1))
recent


ggsave(plot=recent, "../Figures/Supplementary/SupFig_Repeat_density_in_recent_fusions.pdf", width = 5, height = 8, units = "in", limitsize=FALSE)

# Example of an old fusion
old <- plot_fused_only(all_windows, merian_info_per_segment, "Aphantopus_hyperantus","M14", "M29")
ggsave(plot=old, "../Figures/Supplementary/SupFig_Repeat_density_in_old_fusion.pdf", width = 6, height = 6, units = "in", limitsize=FALSE)
ggsave(plot=old, "../Figures/Supplementary/SupFig_Repeat_density_in_old_fusion.png", width = 6, height = 6, units = "in", limitsize=FALSE)




# Ancient fusions
e <- plot_fused_vs_ancestral(all_windows, merian_info_per_segment, "Tinea_trinotella", "Micropterix_aruncella", "M17", "M20")
f <- plot_fused_vs_ancestral(all_windows, merian_info_per_segment, "Polyommatus_icarus", "Colias_croceus", "M12", "M29")
g <- plot_fused_vs_ancestral(all_windows, merian_info_per_segment, "Apotomis_betuletana", "Acleris_emargana", "M13", "M5")
h <- plot_fused_vs_ancestral(all_windows, merian_info_per_segment, "Deilephila_porcellus", "Bombyx_mori", "M24", "M25")


# output plotted tsv tables to save as source data
get_plotted_df_fused_only <- function(df, merian_info_per_segment, spp1, merian1, merian2){
  merian1_df <- merian_info_per_segment %>% filter(species==spp1) %>% filter(merian_segment == merian1)
  merian1_fusion_coord <- merian1_df$prop_stop # get pos of first merian ending 
  merian2_df <- merian_info_per_segment %>% filter(species==spp1) %>% filter(merian_segment == merian2)
  merian2_fusion_coord <- merian2_df$prop_stop # get pos of first merian ending
  fusion_coord <- min(merian1_fusion_coord, merian2_fusion_coord) # use whichever is smaller i.e. which merian is first along chr
  
  
  spp1_df <- df[df$species == spp1,] 
  spp1_df <- spp1_df[spp1_df$assigned_merian == paste(merian1, ',', merian2, sep = ''),]
  spp1_df$plot_colour <- "green" # sanity check - shouldn't see any green
  
  if (merian1_fusion_coord < merian2_fusion_coord) { # if merian1 is before merian2 along the chr
    spp1_df$plot_colour[spp1_df$prop_midpos < merian2_fusion_coord] <- merian2 
    spp1_df$plot_colour[spp1_df$prop_midpos < merian1_fusion_coord] <- merian1
  } else { # else if merian2 is before merian1
    spp1_df$plot_colour[spp1_df$prop_midpos < merian1_fusion_coord] <- merian1
    spp1_df$plot_colour[spp1_df$prop_midpos < merian2_fusion_coord] <- merian2 
  }
  
  merian1_max <- max(spp1_df$prop_midpos, na.rm=TRUE)
  xlim_max <- max(merian1_max)

  return(spp1_df)
}

make_fused_vs_ancestral_df <- function(df, merian_info_per_segment, spp1, spp2, merian1, merian2){
  merian1_df <- merian_info_per_segment %>% filter(species==spp1) %>% filter(merian_segment == merian1)
  merian1_fusion_coord <- merian1_df$prop_stop # get pos of first merian ending 
  merian2_df <- merian_info_per_segment %>% filter(species==spp1) %>% filter(merian_segment == merian2)
  merian2_fusion_coord <- merian2_df$prop_stop # get pos of first merian ending
  fusion_coord <- min(merian1_fusion_coord, merian2_fusion_coord) # use whichever is smaller i.e. which merian is first along chr
  
  spp1_df <- df[df$species == spp1,] 
  spp2_df <- df[df$species == spp2,] 
  spp1_df <- spp1_df[spp1_df$assigned_merian == paste(merian1, ',', merian2, sep = ''),]
  spp2_df <- spp2_df[spp2_df$assigned_merian %in% c(merian1, merian2),]
  spp1_df$plot_colour <- "green" # sanity check - shouldn't see any green
  chr_lengths <- spp2_df %>% group_by(chr) %>% summarise(Value = max(prop_midpos))
  ratio_plots <- chr_lengths[1,2] / chr_lengths[2,2]
  
  if (merian1_fusion_coord < merian2_fusion_coord) { # if merian1 is before merian2 along the chr
    spp1_df$plot_colour[spp1_df$prop_midpos < merian2_fusion_coord] <- merian2 
    spp1_df$plot_colour[spp1_df$prop_midpos < merian1_fusion_coord] <- merian1
  } else { # else if merian2 is before merian1
    spp1_df$plot_colour[spp1_df$prop_midpos < merian1_fusion_coord] <- merian1
    spp1_df$plot_colour[spp1_df$prop_midpos < merian2_fusion_coord] <- merian2 
  }
  
  spp2_df$plot_colour <- spp2_df$assigned_merian
  spp2_df$plot_colour <- factor(spp2_df$plot_colour, levels=c("blue", merian1, merian2))
  spp2_chr1_df <- spp2_df[spp2_df$assigned_merian == merian1,]
  spp2_chr2_df <- spp2_df[spp2_df$assigned_merian == merian2,]
  output_df <- rbind(spp1_df, spp2_chr1_df, spp2_chr2_df)
  return(output_df)
}
# prop_midpos, y=repeat_density,plot_colour)
old_df <- get_plotted_df_fused_only(all_windows, merian_info_per_segment, "Aphantopus_hyperantus","M14", "M29")
nrow(old_df)
write.table(old_df, file = "../../Chromosome_evolution_Lepidoptera_MS/data/repeat_density_across_fusion_chr_in_Aphantopus_hyperantus_151223.tsv", row.names=FALSE, sep="\t", quote = FALSE)


supfig_A_df <- make_fused_vs_ancestral_df(all_windows, merian_info_per_segment, "Agrochola_circellaris", "Agrochola_macilenta", "M30", "M5")
supfig_B_df <- make_fused_vs_ancestral_df(all_windows, merian_info_per_segment, "Dendrolimus_kikuchii", "Dendrolimus_punctatus", "M31", "MZ")
supfig_C_df <- make_fused_vs_ancestral_df(all_windows, merian_info_per_segment, "Agonopterix_subpropinquella", "Agonopterix_arenella", "M19", "M21")

write.table(supfig_A_df, file = "../../Chromosome_evolution_Lepidoptera_MS/data/repeat_density_across_fusion_chr_in_Agrochola_circellaris_151223.tsv", row.names=FALSE, sep="\t", quote = FALSE)
write.table(supfig_B_df, file = "../../Chromosome_evolution_Lepidoptera_MS/data/repeat_density_across_fusion_chr_in_Dendrolimus_kikuchii_151223.tsv", row.names=FALSE, sep="\t", quote = FALSE)
write.table(supfig_C_df, file = "../../Chromosome_evolution_Lepidoptera_MS/data/repeat_density_across_fusion_chr_in_Agonopterix_subpropinquella_151223.tsv", row.names=FALSE, sep="\t", quote = FALSE)


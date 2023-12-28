## Figure 2B
## 'Chromosome evolution in Lepidoptera'
## Wright et al., 2022
## cw22@sanger.ac.uk

library(tidyverse)
library(ggpubr) 

## set input paths
setwd("/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera_MS/")

## load functions
source('2_merian_elements/functions_busco_painter.R')

## import data
assignments <- read.csv('sup_tables/tableS10_chromosome_statistics.tsv', sep='\t')[,c(2,19)] # 6233
# NB: currently misses H. sara - need to fix
colnames(assignments) <- c('query_chr', 'rearrangement_status')
buscos_file_path <- '../Chromosome_evolution_Lepidoptera/Data/buscopaint/'

## specify set of species to plot
spp_list <- c('Micropterix_aruncella', 'Melitaea_cinxia', 'Bombyx_mori', 'Spilosoma_lubricipeda', 'Pieris_napi',
              'Diarsia_rubi', 'Tinea_trinotella')

## custom copy of function with tweaked text sizes for making this figure
busco_paint_coloured_chr <- function(spp_df, num_col, title, karyotype){
  merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
  rearrangement_order = c('ancestral', 'fusion') # leaving out 'split' as only want to distinguish between ancestral vs rearranged
  col_palette <- hue_pal()(32) 
  merian_colours <- c(
    "MZ" = col_palette[1], "M1" = col_palette[2], "M2" = col_palette[3], 
    "M3"= col_palette[4], "M4" = col_palette[5], "M5" = col_palette[6], 
    "M6" = col_palette[7], "M7" = col_palette[8], "M8" =  col_palette[9], 
    "M9" = col_palette[10], "M10" = col_palette[11], "M11" = col_palette[12], 
    "M12" = col_palette[13], "M13" = col_palette[14], "M14" = col_palette[15], 
    "M15" = col_palette[16],"M16" = col_palette[17], "M17" = col_palette[18], 
    "M18" = col_palette[19], "M19" = col_palette[20], "M20" = col_palette[21], 
    "M21" = col_palette[22], "M22" = col_palette[23], "M23" = col_palette[24], 
    "M24" = col_palette[25],"M25" = col_palette[26], "M26" = col_palette[27], 
    "M27" = col_palette[28],"M28" =  col_palette[29], "M29" = col_palette[30], 
    "M30"= col_palette[31], "M31" = col_palette[32], "self"= "grey")
  rearrangement_colours <- c("ancestral" = "black", "fusion"= "red", "split" = "red")
  
  sub_title <- paste("n =", karyotype) 
  the_plot <- ggplot(data = spp_df) +
    scale_colour_manual(values=merian_colours, aesthetics=c("fill"), breaks=merian_order) +
    scale_colour_manual(values=rearrangement_colours, aesthetics=c("colour"), breaks=rearrangement_order, labels=c('Ancestral', 'Rearranged')) + # labels are custom text for legend elements. Breaks are the elements you want in the legend.
    geom_rect(aes(xmin=start, xmax=length, ymax=0, ymin =12, colour=rearrangement_status), fill="white") + #colour="black" if don't want to apply red/black to boxes
    geom_rect(aes(xmin=position-2e4, xmax=position+2e4, ymax=0, ymin =12, fill=status_f)) + # was fill=status_f
    facet_wrap(query_chr_f ~., ncol=num_col, strip.position="right") + guides(scale="none") + # was query_chr - query_chr_f allows order by status (i.e. put fusion/split chr at top)
    xlab("Position (Mb)") +
    scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
    scale_y_continuous(breaks=NULL) + labs(fill = "Merian element") +
    theme(strip.text.y = element_blank(), 
          strip.background = element_blank()) +
    theme(strip.text.x = element_text(margin = margin(0,0,0,0, "cm")), 
          panel.background = element_rect(fill = "white", colour = "white"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.line.x = element_line(color="black", size = 0.5)) +
    theme(legend.position="right") + 
    ggtitle(label=title, subtitle= sub_title)  + 
    theme(plot.title = element_text(hjust = 0.5, face="italic", size=25),
          plot.subtitle = element_text(hjust = 0.5, size=25)) + 
    # strip.text.y = element_text(angle = 0)) +
    guides(fill=guide_legend("Merian Element"), color = "none")  +
    guides(colour=guide_legend("Rearrangement status"), color="none") +
    theme(axis.text.x = element_text(size=20))
  return(the_plot)
}

location_set <- data_frame()
for (i in spp_list){ # read in busco data
  print(i)
  args1 <- paste(buscos_file_path, 'all_spp/', i, '_complete_location.tsv', sep='')
  args2 <- paste(buscos_file_path, 'fai/', i, '.fasta.fai', sep='')
  args3 <- i
  locations_filt <- prepare_data(args1, args2, args3)
  locations_filt$spp <- args3
  location_set <- rbind(location_set, locations_filt)
}
 
location_set <- merge(location_set, assignments, by='query_chr') # merge with fusion/fission status

# make a plot per species
Microperix_paint <- get_busco_paint('Micropterix_aruncella', location_set) + rremove("xlab") +  theme(legend.position = "none")
Melitaea_paint <- get_busco_paint('Melitaea_cinxia', location_set) + rremove("xlab") +  theme(legend.position = "none")
Bombyx_paint <- get_busco_paint('Bombyx_mori', location_set) + rremove("xlab") +  theme(legend.position = "none")
Spilosoma_paint <- get_busco_paint('Spilosoma_lubricipeda', location_set) + rremove("xlab") +  theme(legend.position = "none")
Diarsia_paint <- get_busco_paint('Diarsia_rubi', location_set) + rremove("xlab") +  theme(legend.position = "none")
Tinea_paint <- get_busco_paint('Tinea_trinotella', location_set) + rremove("xlab") +  theme(legend.position = "none")
Pieris_paint <- get_busco_paint('Pieris_napi', location_set) + rremove("xlab") +  theme(legend.position = "none")

legend_paint <- get_busco_paint('Pieris_napi', location_set) + rremove("xlab") +
  theme(legend.key.size = unit(1.1, 'cm'),
        legend.text = element_text(size=18), 
        legend.title = element_text(size=20)) 

legend <- get_legend(legend_paint)
  
# combine all species into one plot
paints_combined <- plot_grid(Microperix_paint,Tinea_paint, Diarsia_paint, Melitaea_paint, Pieris_paint,legend, nrow=1)
paints_combined <- annotate_figure(paints_combined, bottom = text_grob("Chr position (Mb)",size=25)) +  theme(plot.margin = unit(c(1, 0, 0, 0), "cm"))
paints_combined <- paints_combined + draw_figure_label(label="B", position = c("top.left"), size=40, fontface = 2)

## save figure
ggsave('../Figures/Fig2B.pdf', device='pdf', width = 16, height = 16, dpi = 300, units = "in")


# save subset_spp_buscos as an output tsv to add to source_data
subset_spp <- c('Micropterix_aruncella', 'Melitaea_cinxia', 'Bombyx_mori', 'Spilosoma_lubricipeda', 'Diarsia_rubi', 'Tinea_trinotella', 'Pieris_napi')
subset_spp_buscos <- location_set[location_set$short_species_id %in% subset_spp, ]
write.table(subset_spp_buscos, file = "2_merian_elements/busco_2_Merian_for_species_plotted_in_fig2B_1223.tsv", row.names=FALSE, sep="\t", quote = FALSE)

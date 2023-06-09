### Funcions for making busco paints in R ####

# dependenices:
library(tidyverse)
library(scales)
# How to run if want to plot without chr colour (all black)
#args1 <- '../Data/buscopaint/all_spp/Bombyx_mori_complete_location.tsv'
#args2 <- '../Data/buscopaint/fai/Bombyx_mori.fasta.fai'
#args3 <- 'Bombyx_mori'

#locations_filt <- prepare_data(args1, args2, args3)
#num_contigs <- as.character(length(unique(locations_filt$query_chr)))
#subset_merians <- set_colour_mapping(locations_filt)
#b <- busco_paint(locations_filt, 1, args3, num_contigs, 'none')

# With chr coloured (ancestral vs fusion vs fission)
#locations_filt <- prepare_data(args1, args2, args3)
#num_contigs <- as.character(length(unique(locations_filt$query_chr)))
#b <- busco_paint_coloured_chr(locations_filt, 1, args3, num_contigs, 'none')

make_location_set <- function(spp_list, assignments){ # wrapper for prepare_data
  location_set <- tibble()
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
  return(location_set)
}


prepare_data <- function(args1, args2, args3){
  locations <- read_tsv(args1)
  contig_lengths <- read_tsv(args2, col_names=FALSE)
  colnames(contig_lengths) <- c('Seq', 'length', 'offset', 'linebases', 'linewidth')
  plot_title <- args3
  # Format location data 
  locations <- locations %>% filter(!grepl(':', query_chr))
  locations <- merge(locations, contig_lengths, by.x="query_chr", by.y="Seq")
  locations$Length <- locations$length *1000000
  locations$start <- 0
  locations_filt <- locations  %>% group_by(query_chr) %>%   # Filter df to only keep query_chr with >2 BUSCO to remove shrapnel
    mutate(n_busco = n()) %>% # Make a new column reporting number BUSCOs per query_chr
    ungroup() %>%
    filter(n_busco > 3)
  total_contigs <- length(unique(locations$query_chr))# Num query_chr before filtering
  #print(total_contigs) # 104 - number of contigs before filtering by number of BUSCOs
  filt_contigs <- length(unique(locations_filt$query_chr)) #   Num query_chr after filtering
  # print(filt_contigs) # 31 - Number of contigs after filtering by num_BUSCOs>3 
  num_removed_contigs <- length(unique(locations$query_chr)) - length(unique(locations_filt$query_chr))   # Num query_ chr removed by filtering
  # print(num_removed_contigs) # Number of contigs removed by filtering
  merian_order <- c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
  locations_filt$status_f =factor(locations_filt$status, levels=merian_order)

  return(locations_filt)
}

set_colour_mapping <- function(locations_filt){ # Set mapping of Merian2colour
  merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
  colour_palette <- append(hue_pal()(32), 'grey')
  status_merians <- unique(locations_filt$status)
  subset_merians <- subset(colour_palette, merian_order %in% status_merians)
  return(subset_merians)
}

# Adapted plot - plots each chr as a box of correct length - all chr are black 
busco_paint <- function(spp_df, num_col, title, karyotype,legend_pos='right'){
  sub_title <- paste("n contigs =", karyotype)
  the_plot <- ggplot(data = spp_df) +
    scale_colour_manual(values=subset_merians, aesthetics=c("colour", "fill")) +
    geom_rect(aes(xmin=start, xmax=length, ymax=0, ymin =12, colour="black"), colour="black", fill="white") + #colour="black" if don't want to apply red/black to boxes
    geom_rect(aes(xmin=position-2e4, xmax=position+2e4, ymax=0, ymin =12, fill=status_f)) + 
    facet_wrap(query_chr ~., ncol=num_col, strip.position="right") + guides(scale="none") + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) +
    xlab("Position (Mb)") +
    scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
    scale_y_continuous(breaks=NULL) + labs(fill = "Merian element") +
    theme(strip.text.y = element_blank()) + # uncomment if want to remove facet i.e. contig labels
    #   theme(strip.text.y = element_text(angle = 0),
    #        strip.background = element_blank()) +
    theme(strip.text.x = element_text(margin = margin(0,0,0,0, "cm")), 
          panel.background = element_rect(fill = "white", colour = "white"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.line.x = element_line(color="black", size = 0.5)) +
    theme(legend.position=legend_pos) + ggtitle(label=title, subtitle= sub_title)  + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    theme(plot.title=element_text(face="italic")) +
    guides(fill=guide_legend("Merian Element"), color = "none")
  return(the_plot)
}


busco_paint_coloured_chr <- function(spp_df, num_col, title, karyotype){
  merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
  rearrangement_order = c('ancestral', 'fusion') # leaving out 'split' as only want to distinguish between ancestral vs rearranged
  col_palette <- hue_pal()(32) # was hue_pal
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

  sub_title <- paste("n =", karyotype) # need to fix 'mycolors' = change to rearrangement_colours?
  the_plot <- ggplot(data = spp_df) +
    scale_colour_manual(values=merian_colours, aesthetics=c("fill"), breaks=merian_order) +
    scale_colour_manual(values=rearrangement_colours, aesthetics=c("colour"), breaks=rearrangement_order, labels=c('Ancestral', 'Rearranged')) + # labels are custom text for legend elements. Breaks are the elements you want in the legend.
    geom_rect(aes(xmin=start, xmax=length, ymax=0, ymin =12, colour=rearrangement_status), fill="white") + #colour="black" if don't want to apply red/black to boxes
    geom_rect(aes(xmin=position-2e4, xmax=position+2e4, ymax=0, ymin =12, fill=status_f)) + # was fill=status_f
    facet_wrap(query_chr_f ~., ncol=num_col, strip.position="right") + guides(scale="none") + # was query_chr - query_chr_f allows order by status (i.e. put fusion/split chr at top)
    xlab("Position (Mb)") +
    scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
    scale_y_continuous(breaks=NULL) + labs(fill = "Merian element") +
   theme(strip.text.y = element_blank(), # comment if want to remove facet i.e. contig labels
        #  strip.text.y = element_text(angle = 0),
          strip.background = element_blank()) +
    theme(strip.text.x = element_text(margin = margin(0,0,0,0, "cm")), 
          panel.background = element_rect(fill = "white", colour = "white"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.line.x = element_line(color="black", size = 0.5)) +
    theme(legend.position="right") + 
    ggtitle(label=title, subtitle= sub_title)  + 
    theme(plot.title = element_text(hjust = 0.5, face="italic", size=20),
          plot.subtitle = element_text(hjust = 0.5, size=15)) + 
         # strip.text.y = element_text(angle = 0)) +
    guides(fill=guide_legend("Merian Element"), color = "none")  +
    guides(colour=guide_legend("Rearrangement status"), color="none")
  return(the_plot)
}

get_busco_paint <- function(query_spp, location_set, num_cols=1){
  filt_buscos <- location_set %>% filter(spp == query_spp)
  num_contigs <- as.character(length(unique(filt_buscos$query_chr)))
  fusion_buscos <- filt_buscos[filt_buscos$rearrangement_status == "fusion",] # find fused chr
  split_buscos <- filt_buscos[filt_buscos$rearrangement_status == "split",] # find fission chr
  ancestral_buscos <- filt_buscos[filt_buscos$rearrangement_status == "ancestral",] # find ancestral chr
  fusion_chr <- unique(fusion_buscos$query_chr)
  fission_chr <- unique(split_buscos$query_chr)
  ancestral_chr <- unique(ancestral_buscos$query_chr)
  chr_levels <- c(fusion_chr, fission_chr, ancestral_chr) # combine all chr but ordered by status
  filt_buscos$query_chr_f =factor(filt_buscos$query_chr, levels=chr_levels) # set chr order as order for plotting
  # the_plot <- busco_paint(filt_buscos, 1, args3, num_contigs, 'none') # uncomment if want without coloured chr (i.e. all black)
  the_plot <- busco_paint_coloured_chr(filt_buscos, num_cols, query_spp, num_contigs) # potentially change colour palette
  return(the_plot)
}



library(tidyverse)
library(scales)

# functions from https://github.com/charlottewright/lep_busco_painter/blob/main/plot_buscopainter.R
prepare_data_with_index <- function(args1, args2){
  locations <- read_tsv(args1, col_types=cols())
  contig_lengths <- read_tsv(args2, col_names=FALSE, col_types = cols())
  colnames(contig_lengths) <- c('Seq', 'length')
  locations <- locations %>% filter(!grepl(':', query_chr)) # format location data 
  locations <- locations %>% filter(!grepl(':', assigned_chr)) # format location data 
  locations <- merge(locations, contig_lengths, by.x="query_chr", by.y="Seq")
  locations$start <- 0
  return(locations)
}

filter_buscos <- function(locations, minimum){ # minimum of buscos to be present
  locations_filt <- locations  %>% 
    group_by(query_chr) %>%   # filter df to only keep query_chr with >=3 buscos to remove shrapnel
    mutate(n_busco = n()) %>% # make a new column reporting number buscos per query_chr
    ungroup() %>%
    filter(n_busco >= minimum)
  return(locations_filt)
}


paint_merians_all <- function(spp_df, num_col, title){
  colour_palette <- append(hue_pal()(27), 'grey')
  merian_order <- c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'unassigned')
  spp_df$assigned_chr_f =factor(spp_df$assigned_chr, levels=merian_order)
  chr_levels <- subset(spp_df, select = c(query_chr, length)) %>% unique() %>% arrange(desc(length))
  chr_levels <- chr_levels$query_chr
  spp_df$query_chr_f =factor(spp_df$query_chr, levels=chr_levels) # set chr order as order for plotting
  the_plot <- ggplot(data = spp_df) +
    scale_colour_manual(values=colour_palette, aesthetics=c("colour", "fill")) +
    geom_rect(aes(xmin=start, xmax=length, ymax=0, ymin =12), colour="black", fill="white") + 
    geom_rect(aes(xmin=position, xmax=position+0.5, ymax=0, ymin =12, fill=assigned_chr_f)) +
    facet_wrap(query_chr_f ~., ncol=num_col, strip.position="left") + guides(scale="none") + 
    theme(strip.text.y.left = element_text(angle=0)) +
    xlab("Relative position") +
    scale_y_continuous(breaks=NULL) + 
    ggtitle(label=title)  +
    guides(fill=guide_legend("Merian element"), color = "none") +
    busco_paint_theme
  return(the_plot)
}

busco_paint_theme <- theme(legend.position="right",
                           strip.text.x = element_text(margin = margin(0,0,0,0, "cm")), 
                           panel.background = element_rect(fill = "white", colour = "white"), 
                           panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(),
                           axis.line.x = element_line(color="black", size = 0.5),
                           axis.text.x = element_text(size=10),
                           axis.title.x = element_text(size=10),
                           strip.text.y = element_text(angle=0),
                           strip.background = element_blank(),
                           plot.title = element_text(hjust = 0.5, face="italic", size=20))

setwd("/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/")

args1 <- '../Data/AGORA/anc212_complete_locations.tsv'
args2 <- '../Data/AGORA/anc212_CAR_lengths.tsv'
args3 <- 'anc212'
minimum <- 10
locations <- prepare_data_with_index(args1, args2)
locations$spp <- args3
locations_filt <- filter_buscos(locations, minimum)
SupFig <- paint_merians_all(locations_filt, 3, '')

ggsave(plot=SupFig, '../Figures/SupFig_merian_paint_agora_anc212.pdf', device='pdf', width = 10, height = 8, dpi = 300, units = "in")

# add asterix to moved BUCO

locations_syngraph<-locations[(locations$status =="unassigned"),] 
locations_all <- filter_buscos(locations_syngraph, 0)
locations_2 <- filter_buscos(locations_syngraph, 2)
locations_5 <- filter_buscos(locations_syngraph, 5)
locations_10 <- filter_buscos(locations_syngraph, 10)

nrow(locations_all) # 4112
nrow(locations_2) # 3093
nrow(locations_5) # 1957
nrow(locations_10) # 1033

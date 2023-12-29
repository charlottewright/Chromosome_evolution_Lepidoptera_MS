


library(wesanderson)

make_pairwise_dataframe <- function(df1, df2){
  df1$midpos_x <- df1$midpos
  df2$midpos_y <- df2$midpos
  pairwise_df = merge(df1, df2, by="busco_id")
  return(pairwise_df)
}

# correct function
make_oxford_plot <- function(df){
  the_plot <- ggplot() +  # scales=free, space=free is what you want to get boxes scaled by length :) 
    facet_grid(Sequence.x_f~Sequence.y_f, space = "free", scales="free")  + 
    geom_point(data = df, aes(x = midpos_y, y = midpos_x, color=Merian), size=0.1) + # was Sequence.x_f 
    theme_bw() + 
    theme(legend.position = "none", axis.ticks = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_rect(color = "#999999", size=0.3)) + 
    theme(panel.spacing = unit(0, "lines")) +
    theme(strip.background = element_blank()) +
    theme(strip.text.y = element_text(angle = 0,size=6)) +
    theme(strip.text.x = element_text(angle = 90,size=6)) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
    geom_point(data=df, aes(x=sequence_y_min, y=sequence_x_min), color="white") + # needed to scale each box
    geom_point(data=df, aes(x=sequence_y_max, y=sequence_x_max), color="white")     # needed to scale each box
  return(the_plot)
}



scale_lengths <- function(buscos){
  buscos$start <- buscos$start / 1000000
  buscos$end <- buscos$end / 1000000
  buscos$midpos <- (buscos$start + buscos$end)/2
  return(buscos)
}

make_pairwise_dataframe <- function(df1, df2){
  df1$midpos_x <- df1$midpos
  df2$midpos_y <- df2$midpos
  pairwise_df = merge(df1, df2, by="busco_id")
  return(pairwise_df)
}

filter_buscos_x <- function(buscos, threshold){
  filtered <- group_by(buscos, Sequence.x) %>%
    mutate(nGenes = n(),
           mxGpos = max(midpos_x)) %>%
    ungroup() %>%
    filter(nGenes >threshold)
  return(filtered)
}

filter_buscos_y <- function(buscos, threshold){
  filtered <- group_by(buscos, Sequence.y) %>%
    mutate(nGenes = n(),
           mxGpos_y = max(midpos_y)) %>%
    ungroup() %>%
    filter(nGenes >threshold)
  return(filtered)
}


setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/')
Merian_assignments_ref <- read_busco('../Data/syngraph/Final_Merians_based_on_syn_n2_from_r2_m5_full_table.tsv')
Merian_assignments_ref <- select(Merian_assignments_ref, -c(start, end))
colnames(Merian_assignments_ref) <- c('busco_id', 'Merian')

Lysandra_coridon <- read_busco('../Data/BUSCOs/All/Lysandra_coridon.tsv')
Polyommatus_icarus <- read_busco('..//Data/BUSCOs/All/Polyommatus_icarus.tsv')

Lysandra_coridon <- scale_lengths(Lysandra_coridon)
Polyommatus_icarus <- scale_lengths(Polyommatus_icarus)

ilLysCorr_vs_ilPolIcar <- make_pairwise_dataframe(Lysandra_coridon, Polyommatus_icarus)
ilLysCorr_vs_ilPolIcar <- filter_buscos_x(ilLysCorr_vs_ilPolIcar,2) %>% filter_buscos_y(1)
ilLysCorr_vs_ilPolIcar <- merge(ilLysCorr_vs_ilPolIcar, Merian_assignments_ref, by="busco_id")
merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
ilLysCorr_vs_ilPolIcar$Merian <- factor(ilLysCorr_vs_ilPolIcar$Merian, levels = merian_order)

# Make columns to define limits of each facet
ilLysCorr_vs_ilPolIcar <- ilLysCorr_vs_ilPolIcar %>% group_by(Sequence.x) %>% mutate(sequence_x_max = max(end.x)) %>% ungroup()
ilLysCorr_vs_ilPolIcar <- ilLysCorr_vs_ilPolIcar %>% group_by(Sequence.y) %>% mutate(sequence_y_max = max(end.y)) %>% ungroup()
ilLysCorr_vs_ilPolIcar$sequence_y_min <- 0
ilLysCorr_vs_ilPolIcar$sequence_x_min <- 0


# Default order of chr is by size 
sequence_x_order <- ilLysCorr_vs_ilPolIcar %>% arrange(mxGpos) %>% select(Sequence.x) %>% unique() %>% pull()
sequence_y_order <- ilLysCorr_vs_ilPolIcar %>% arrange(mxGpos_y) %>% select(Sequence.y) %>% unique() %>% pull()
sequence_x_order <- rev(sequence_x_order)

# not sure about OW569334.1
sequence_y_order <- c("OW569337.1", "OW569323.1", "OW569335.1", "OW569338.1", "OW569331.1", "OW569340.1", "OW569332.1", 
                      "OW569333.1", "OW569336.1", "OW569326.1", "OW569328.1",
                      "OW569329.1","OW569342.1", "OW569321.1", "OW569327.1",
                      "OW569339.1", "OW569324.1", "OW569325.1", "OW569322.1",
                      "OW569330.1","OW569341.1","OW569334.1", "OW569320.1")

ilLysCorr_vs_ilPolIcar$Sequence.x_f <- factor(ilLysCorr_vs_ilPolIcar$Sequence.x, levels = sequence_x_order)
ilLysCorr_vs_ilPolIcar$Sequence.y_f <- factor(ilLysCorr_vs_ilPolIcar$Sequence.y, levels = sequence_y_order)

ilLysCorr_vs_ilPolIcar_plot <- make_oxford_plot(ilLysCorr_vs_ilPolIcar) + 
  xlab("Polyommatus icarus") +   ylab("Lysandra coridon") +
  theme(axis.title = element_text(face="italic"),
        legend.position = "bottom") + 
  guides(colour=guide_legend(override.aes = list(size=4), title="Merian element", nrow=3,byrow=TRUE)) 

ggsave(plot=ilLysCorr_vs_ilPolIcar_plot, '../Figures/Supplementary/SupFig_Oxford_plot_lysandra_coridon_vs_polyommatus_icarus.png', width=10, height=10)



### Funcions for making Oxford Plots ####
# Function for loading buscos
read_busco <- function(buscoFile){
  busco <- read_tsv(buscoFile,
                    col_names = c("busco_id", "Status", "Sequence",
                                  "start", "end", "strand", "Score", "Length",
                                  "OrthoDB_url", "Description"),
                    col_types = c("ccciicdicc"),
                    comment = "#") %>%
    filter(Status == "Complete") %>%
    select(busco_id, Sequence, start, end)
  
  return(busco)
}

# Function for creating multispecies coordinates index from busco_tibble
# busco_tibble needs to have columns "multispecies_sequence" and "end"
buscos_to_multi_sp_index <- function(busco_tb){
  multi_sp_index <- group_by(busco_tb, multispecies_sequence) %>%
    summarise(size = max(end), .groups = "drop") %>%
    arrange(size) %>%
    #    arrange(multispecies_sequence) %>%
    mutate(multispecies_sequence = factor(multispecies_sequence),
           multispecies_index = cumsum(size) - size) %>%
    select(multispecies_sequence, multispecies_index, size)
  
  return(multi_sp_index)
}

# Function for multispecies Oxford plot
multispecies_oxford_plot <- function(busco_tbl, sp_x, sp_y, colour_palette=hue_pal()(32)){
  
  # if(is.null(colour_palette))
  #  colour_palette <- hue_pal()(32) # specify colours - default is hue_pal i.e. ggplot default colours
  
  merian_colors <- colour_palette
  # manually specifying busco_tbl:
  #busco_tbl <- filt_buscos
  #sp_x <- x_species
  #sp_y <- y_species
  
  x_busco <- filter(busco_tbl, assembly %in% sp_x)
  y_busco <- filter(busco_tbl, assembly %in% sp_y) # comment out if want to plot many spp on y-axis
#  y_busco <- busco_tbl # uncomment if want to plot many spp on y axis (i.e. dont filter table)
  # Make buscos table of just x_spp and y_spp
  buscos_to_cmpr <- left_join(x_busco, y_busco,
                              by = "busco_id")
  
  # Make df of 'multispecies_sequence' vs 'end' for spp_x
  multi_sp_index_x <- select(buscos_to_cmpr, 
                             multispecies_sequence = multispecies_sequence.x,
                             end = end.x) 
  # Order sequences by size (max end pos) then convert to additive length (i.e. index pos)
  multi_sp_index_x <- buscos_to_multi_sp_index(multi_sp_index_x)
  # Do same for spp_y
  multi_sp_index_y <- select(buscos_to_cmpr, 
                             multispecies_sequence = multispecies_sequence.y,
                             end = end.y) %>% buscos_to_multi_sp_index()
  
  
  # Merge with Merian info??
  merian_order = c('self', 'MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30')
  buscos_to_cmpr <- left_join(buscos_to_cmpr, Merian_assignments_ref, by="busco_id")
  buscos_to_cmpr$Merian_f = factor(buscos_to_cmpr$Merian, levels=merian_order)
  
  # Extrapolate USCOs to multispecies coordinates
  multispecies_buscos <- left_join(buscos_to_cmpr, multi_sp_index_x,
                                   by = c("multispecies_sequence.x" = "multispecies_sequence")) %>%
    left_join(multi_sp_index_y,
              by = c("multispecies_sequence.y" = "multispecies_sequence")) %>%
    mutate(multispecies_start.x = start.x + multispecies_index.x,
           multispecies_start.y = start.y + multispecies_index.y) %>%
    select(busco_id, multispecies_start.x, multispecies_start.y,
           multispecies_sequence.x, multispecies_sequence.y, Merian_f)
  
  # Main plot
  p <- ggplot(multispecies_buscos, aes(x = multispecies_start.x,
                                       y = multispecies_start.y, colour=Merian_f)) +  scale_colour_manual(values=merian_colors) +
    geom_point(alpha = 1, size=0.7) +
   scale_x_continuous(breaks = multi_sp_index_x$multispecies_index + (multi_sp_index_x$size)/2,
        #              labels = multi_sp_index_x$multispecies_sequence, # uncomment if want labels of chr on y-axis
                      expand = c(0, 0)) +
   scale_y_continuous(breaks = multi_sp_index_y$multispecies_index + (multi_sp_index_y$size)/20,
                  #   labels = multi_sp_index_y$multispecies_sequence, # uncomment if want labels of chr on y-axis
                    expand = c(0, 0),
                    guide=guide_axis(angle=90))  +
    geom_vline(xintercept=multi_sp_index_x$multispecies_index, size = 1, colour="black") + #was size=0.2, colour="grey"
    geom_hline(yintercept=multi_sp_index_y$multispecies_index, size = 1, colour="black") + # was size=0.2, colour="grey"
    theme_bw() +
   theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
   #    axis.text.y = element_text(angle=90), 
        axis.text = element_blank(), # uncomment if want labels of chr on y-axis
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.title = element_text(face = "italic")) +
    labs(x=sp_x, y=sp_y) 
  return(p)
}

filter_buscos <- function(buscos, multispecies_sequence){
  filtered <- group_by(buscos, multispecies_sequence) %>%
    mutate(nGenes = n(),
           mxGpos = max(start)) %>%
    ungroup() %>%
    filter(nGenes >2)
  return(filtered)
}

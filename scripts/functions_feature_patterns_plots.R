
# Functions associated with plotting feature patterns along chr
# Used in 'Fig5_feature_patterns.R'
plot_superfamily <- function(superfamily_query, dataframe, feature, corr_feature){
  axis_text_size = 12
  group.colors <- c("TRUE" = "#00A08A", "FALSE" = "#F98400") # needed because if a df only has 'FALSE", it would default to the first colour rather than second
  dataframe$feature <- dataframe[, feature] 
  dataframe$corr_feature <- dataframe[, corr_feature]
  superfamily_query <- dataframe[dataframe$superfamily == superfamily_query, ]
  the_plot <- ggplot(data=superfamily_query, aes(x=prop_length, y=feature,group=species)) +
    #geom_point(aes(alpha=0.2,colour=corr_feature),size=0.2) + 
    #facet_grid(~species) +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.text = element_text(size=axis_text_size)) +
    scale_color_manual(values = group.colors) +
    geom_smooth(span=1.2, se=FALSE, aes(color=corr_feature), size=0.4) +
    #  geom_smooth(method = "lm", se=FALSE, size=0.4, aes(color=corr_feature)) + # "loess
    labs(x="",y="") +
    theme(plot.margin = margin(0.001, 0.001, 0.001, 0.001, "cm")) +
    scale_x_continuous(limits = c(1,5.7)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) 
  return(the_plot)
}

plot_superfamily_with_MZ <- function(superfamily_query, dataframe, feature, corr_feature,ancestral_MZ){
  axis_text_size = 12
  Z_df <- ancestral_MZ[ancestral_MZ$superfamily == superfamily_query,]
  group.colors <- c("TRUE" = "#00A08A", "FALSE" = "#F98400") # needed because if a df only has 'FALSE", it would default to the first colour rather than second
  dataframe$feature <- dataframe[, feature] 
  Z_df$feature <- Z_df[, feature] 
  dataframe$corr_feature <- dataframe[, corr_feature]
  superfamily_query <- dataframe[dataframe$superfamily == superfamily_query, ]
  the_plot <- ggplot(data=superfamily_query, aes(x=prop_length, y=feature,group=species)) +
    #geom_point(aes(alpha=0.2,colour=corr_feature),size=0.2) + 
    #facet_grid(~species) +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.text = element_text(size=axis_text_size)) +
    scale_color_manual(values = group.colors) +
    geom_smooth(span=1.2, se=FALSE, aes(color=corr_feature), size=0.4) +
    #  geom_smooth(method = "lm", se=FALSE, size=0.4, aes(color=corr_feature)) + # "loess
    labs(x="",y="") +
    theme(plot.margin = margin(0.001, 0.001, 0.001, 0.001, "cm")) +
    scale_x_continuous(limits = c(1,5.7)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
  geom_point(data=Z_df, aes(x=prop_length, y=feature), color="#F21A00",size=0.5) # uncomment if want red dots for Z-chr
  
  return(the_plot)
}

plot_superfamily_with_axes <- function(superfamily_query, dataframe, feature, corr_feature){
  axis_text_size = 12
  group.colors <- c("TRUE" = "#00A08A", "FALSE" = "#F98400")
  dataframe$feature <- dataframe[, feature] 
  dataframe$corr_feature <- dataframe[, corr_feature]
  superfamily_query <- dataframe[dataframe$superfamily == superfamily_query, ]
  the_plot <- ggplot(data=superfamily_query, aes(x=prop_length, y=feature,group=species)) +
    # geom_point(aes(alpha=0.2,colour=corr_feature),size=1) + 
    #facet_grid(~species) +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.text = element_text(size=axis_text_size)) +
    scale_color_manual(values = group.colors) +
    geom_smooth(span=1.2, se=FALSE, aes(color=corr_feature), size=0.4) +
    #    geom_smooth(method = "lm", se=FALSE, size=0.4, aes(color=corr_feature)) + # loess
    labs(x="",y="") +
    theme(plot.margin = margin(0.001, 0.001, 0.001, 0.001, "cm")) +
    scale_x_continuous(limits = c(1,5.7)) 
  return(the_plot)
}
plot_superfamily_with_MZ_and_axes <- function(superfamily_query, dataframe, feature, corr_feature,ancestral_MZ){
  axis_text_size = 12
  Z_df <- ancestral_MZ[ancestral_MZ$superfamily == superfamily_query,]
  group.colors <- c("TRUE" = "#00A08A", "FALSE" = "#F98400")
  dataframe$feature <- dataframe[, feature] 
  Z_df$feature <- Z_df[, feature] 
  dataframe$corr_feature <- dataframe[, corr_feature]
  superfamily_query <- dataframe[dataframe$superfamily == superfamily_query, ]
  the_plot <- ggplot(data=superfamily_query, aes(x=prop_length, y=feature,group=species)) +
    # geom_point(aes(alpha=0.2,colour=corr_feature),size=1) + 
    #facet_grid(~species) +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.text = element_text(size=axis_text_size)) +
    scale_color_manual(values = group.colors) +
    geom_smooth(span=1.2, se=FALSE, aes(color=corr_feature), size=0.4) +
    #    geom_smooth(method = "lm", se=FALSE, size=0.4, aes(color=corr_feature)) + # loess
    labs(x="",y="") +
    theme(plot.margin = margin(0.001, 0.001, 0.001, 0.001, "cm")) +
    scale_x_continuous(limits = c(1,5.7)) 
  # +
  # geom_point(data=Z_df, aes(x=prop_length, y=feature), color="#F21A00",size=0.5) # uncomment if want red dots for Z-chr
  return(the_plot)
}

get_average_feature_density_per_species <- function(df, feature_column){
  df$feature_column <- df[, feature_column]
  df$feature_mult <- df$feature_column * df$length
  total_length_per_species <- df %>%
    dplyr:::group_by(species) %>% filter(!is.na(col1)) %>%
    dplyr:::summarise(total_length = sum(length))
  
  sum_feature_per_spp <- df %>%
    dplyr:::group_by(species) %>% filter(!is.na(feature_mult)) %>%
    dplyr:::summarise(total_feature_mult = sum(feature_mult)) 
  
  feature_scaled <- merge(sum_feature_per_spp, total_length_per_species, by="species")
  feature_scaled$feature_scale_factor <- feature_scaled$total_feature_mult / feature_scaled$total_length
  feature_scaled <- feature_scaled[,c(1,4)]
  colnames(feature_scaled) <- c('species', paste(feature_column, "_scale_factor", sep=""))
  df <-  merge(df, feature_scaled, by="species", all.x=TRUE)
  return(df)
}

make_panel_across_superfamily <- function(df, superfamily){
  group.colors <- c("Autosome" = "#B4B4B4", "neoZ" = "#046C9A", "Z"="#F21A00")
  subset_df <- df[df$superfamily == superfamily,]
  subset_df_Z <- filter(subset_df, grepl("MZ",merians))
  subset_df$Chromosome <- "Autosome"
  subset_df$Chromosome[(subset_df$merians %in% subset_df_Z$merians)]<- "neoZ"
  subset_df$Chromosome[subset_df$merians == "MZ"] <- "Z"
  the_plot <-  ggplot(subset_df, aes(x=prop_length, y=repeat_density)) + geom_point(aes(colour=Chromosome)) + 
    facet_wrap(~species, scales="free",ncol=6) + theme_bw() +  
    theme(strip.text.x = element_text(size = 7)) +
    scale_color_manual(values=group.colors) +
    guides(color = guide_legend(override.aes = list(size = 5))) + xlab("Prop chr length (%)") + ylab("Repeat density (%)")
  return(the_plot)
}

make_half_panel <- function(subset_df){
  group.colors <- c("autosome" = "#B4B4B4", "neoZ" = "#046C9A", "Z"="#F21A00")
  subset_df_Z <- filter(subset_df, grepl("MZ",merians))
  subset_df$Chromosome <- "autosome"
  subset_df$Chromosome[(subset_df$merians %in% subset_df_Z$merians)]<- "neoZ"
  subset_df$Chromosome[subset_df$merians == "MZ"] <- "Z"
  the_plot <-  ggplot(subset_df, aes(x=prop_length, y=repeat_density)) + geom_point(aes(colour=Chromosome)) + 
    facet_wrap(~species, scales="free",ncol=6) + theme_bw() +  
    theme(strip.text.x = element_text(size = 7)) +
    scale_color_manual(values=group.colors) +
    guides(color = guide_legend(override.aes = list(size = 5))) + xlab("Prop chr length (%)") + ylab("Repeat density (%)")
  return(the_plot)
}
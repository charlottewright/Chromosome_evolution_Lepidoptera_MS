
library(dplyr)
library(ggplot2)

## plotting function
plot_two_df_features <- function(df1_autosomes, df1_MZ, df2_autosomes, df2_MZ){
  p <- ggplot(data=df1_autosomes, aes(x=length_mb, y=intergenic_size_mb, color=chr_type)) + geom_point() +
    theme_classic() + geom_smooth(method = "loess", se=FALSE) +
    geom_point(data = df1_MZ, aes(x=length_mb, y=intergenic_size_mb, color=chr_type)) + 
    labs(y="Mean intergenic distance (Mb)", x="Chr length (Mb)") + labs(color="") +
    scale_color_manual(values=c("black", "red")) 
  
  q <- ggplot(data=df1_autosomes, aes(x=length_mb, y=intron_lengths, color=chr_type)) + geom_point() +
    theme_classic() + geom_smooth(method = "loess", se=FALSE) +
    geom_point(data = df1_MZ, aes(x=length_mb, y=intron_lengths, color=chr_type)) + 
    labs(y="Mean intron length (bp)", x="Chr length (Mb)") + labs(color="") +
    scale_color_manual(values=c("black", "red"))  
  
  r <- ggplot(data=df1_autosomes, aes(x=length_mb, y=intron_numbers, color=chr_type)) + geom_point() +  
    theme_classic() + geom_smooth(method = "loess", se=FALSE) +
    geom_point(data = df1_MZ, aes(x=length_mb, y=intron_numbers, color=chr_type)) + 
    labs(y="Mean number of introns per gene", x="Chr length (Mb)") + labs(color="") +
    scale_color_manual(values=c("black", "red"))
  
  s <- ggplot(data=df2_autosomes, aes(x=length_mb, y=intergenic_size_mb, color=chr_type)) + geom_point() +
    theme_classic() + geom_smooth(method = "loess", se=FALSE) +
    geom_point(data = df2_MZ, aes(x=length_mb, y=intergenic_size_mb, color=chr_type)) + 
    labs(y="Mean intergenic distance (Mb)", x="Chr length (Mb)") + labs(color="") +
    scale_color_manual(values=c("black", "red")) 
  
  t <- ggplot(data=df2_autosomes, aes(x=length_mb, y=intron_lengths, color=chr_type)) + geom_point() +
    theme_classic() + geom_smooth(method = "loess", se=FALSE) +
    geom_point(data = df2_MZ, aes(x=length_mb, y=intron_lengths, color=chr_type)) + 
    labs(y="Mean intron length (bp)", x="Chr length (Mb)") + labs(color="") +
    scale_color_manual(values=c("black", "red"))  
  
  u <- ggplot(data=df2_autosomes, aes(x=length_mb, y=intron_numbers, color=chr_type)) + geom_point() +  
    theme_classic() + geom_smooth(method = "loess", se=FALSE) +
    geom_point(data = df2_MZ, aes(x=length_mb, y=intron_numbers, color=chr_type)) + 
    labs(y="Mean number of introns per gene", x="Chr length (Mb)") + labs(color="") +
    scale_color_manual(values=c("black", "red"))
  
  plots_combined <- p + q + r +s + t + u + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  return(plots_combined)
}

make_df <- function(df_intergenic, df_intron_lengths, df_intron_numbers, assignments){
  colnames(df_intergenic) <- c('chr', 'intergenic_size')
  colnames(df_intron_lengths) <- c('chr', 'intron_lengths')
  colnames(df_intron_numbers) <- c('chr', 'intron_numbers')
  
  df <- merge(df_intergenic, assignments, by="chr") %>% merge(df_intron_lengths, by="chr") %>% merge(df_intron_numbers, by="chr")
  
  df$intergenic_size_mb <- df$intergenic_size / 100000
  df$length_mb <- df$length / 100000
  
  df$chr_type <- "Autsosomal Merian element"
  df$chr_type[df$assigned_merian == "MZ"] <- "Merian element Z"
  return(df)
}

setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Data/genes/')

M_intergenic <- read.csv('Melitaea_cinxia_average_intergenic_distance_per_chr.tsv', sep='\t', header=FALSE)
M_intron_lengths <- read.csv('Melitaea_cinxia_average_intron_length_per_chr.tsv', sep='\t', header=FALSE)
M_intron_numbers <- read.csv('Melitaea_cinxia_average_number_introns_per_chr.tsv', sep='\t', header=FALSE)
assignments <- read.csv('../../../Chromosome_evolution_Lepidoptera_MS/sup_tables/tableS10_chromosome_statistics.tsv', sep='\t')[,c(3,4,6,20)]

D_intergenic <- read.csv('Diarsia_rubi_average_intergenic_distance_per_chr.tsv', sep='\t', header=FALSE)
D_intron_lengths <- read.csv('Diarsia_rubi_average_intron_length_per_chr.tsv', sep='\t', header=FALSE)
D_intron_numbers <- read.csv('Diarsia_rubi_average_number_introns_per_chr.tsv', sep='\t', header=FALSE)

D_df <- make_df(D_intergenic, D_intron_lengths, D_intron_numbers, assignments)
D_autosomes = D_df[D_df$assigned_merian != "MZ",]
D_MZ = D_df[D_df$assigned_merian == "MZ",]

M_df <- make_df(M_intergenic, M_intron_lengths, M_intron_numbers, assignments)
M_autosomes = M_df[M_df$assigned_merian != "MZ",]
M_MZ = M_df[M_df$assigned_merian == "MZ",]

D_plot <- plot_features(D_autosomes, D_MZ)
M_plot <- plot_features(M_autosomes, M_MZ)
both_plot <- plot_two_df_features(D_autosomes, D_MZ, M_autosomes, M_MZ)

ggsave(plot=both_plot, '../../Figures/Response_to_reviewers/Prop_length_vs_intergic_and_intron_stats_D_rubi_and_M_cinxia.pdf', device = 'pdf', width = 9, height = 7, units = "in", limitsize=FALSE) #  save figure
ggsave(plot=both_plot, '../../Figures/Response_to_reviewers/Prop_length_vs_intergic_and_intron_stats_D_rubi_and_M_cinxia.png', device = 'png', width = 9, height = 7, units = "in", limitsize=FALSE) #  save figure

## calculate correlation strengths
corr_intergenic <- cor.test(x=D_autosomes$length, y=D_autosomes$intergenic_size_mb, method = 'spearman')
corr_intron_size <- cor.test(x=D_autosomes$length, y=D_autosomes$intron_lengths, method = 'spearman')
corr_intron_number <- cor.test(x=D_autosomes$length, y=D_autosomes$intron_numbers, method = 'spearman')

corr_intergenic
corr_intron_size
corr_intron_number

  
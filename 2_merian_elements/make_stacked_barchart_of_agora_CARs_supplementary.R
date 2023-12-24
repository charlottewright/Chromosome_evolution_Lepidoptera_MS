library(dplyr)
library(ggplot2)

setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/')

df <- read.csv('../Data/AGORA/comparison_merians_to_agora_anc212.tsv', sep='\t', header=TRUE)

head(df)

merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')

df$merian <- factor(df$merian,levels = merian_order)

SupFig_stacked <- ggplot(df, aes(fill=CAR, y=count, x=merian)) + 
  geom_bar(position="stack", stat="identity") + 
  geom_col(color = "#0099f9", fill = "#ffffff") +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + labs(y='Number of BUSCOs', x='Merian element')

ggsave(plot=SupFig_stacked, '../Figures/SupFig_stacked_barchart_agora_anc212.pdf', device='pdf', width = 10, height = 8, dpi = 300, units = "in")
ggsave(plot=SupFig_stacked, '../Figures/SupFig_stacked_barchart_agora_anc212.png', device='png', width = 10, height = 8, dpi = 300, units = "in")


df_longest_CARs <- df %>% filter(count >= 15)



library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(rstatix)
setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/')

# read in orthofinder output
df <- read.csv('Data/ALG_comparisons/single_copy_orthogroups_with_ALGs_assigned_Merians_accurately.tsv', sep='\t', header=FALSE)[,c(4,1,2,6)]
colnames(df) <- c('BCS_LG', 'Aquee_prot', 'Mcinx_prot', 'Merian')

df$BLG <- df$BCS_LG
df$BLG <- gsub("A1a", "A1_mixed", df$BLG)
df$BLG <- sub("A1b", "A1_mixed", df$BLG)
df$BLG <- sub("A1_mixed", "A1axA1b", df$BLG)

df$BLG <- gsub("Ea", "E_mixed", df$BLG)
df$BLG <- gsub("Eb", "E_mixed", df$BLG)
df$BLG <- gsub("E_mixed", "EaxEb", df$BLG)

df$BLG <- gsub("Qa", "Q_mixed", df$BLG)
df$BLG <- gsub("Qb", "Q_mixed", df$BLG)
df$BLG <- gsub("Q_mixed", "QaxQb", df$BLG)

df$BLG <- gsub("Qc", "Q_mixed", df$BLG)
df$BLG <- gsub("Qd", "Q_mixed", df$BLG)
df$BLG <- gsub("Q_mixed", "QcxQd", df$BLG)

df <- df %>% group_by(BLG, Merian) %>% summarise('matches' = n())

# create a df with all possible combos of Merin & BCS
all_combinations <- expand.grid(BLG = unique(df$BLG), 
                                Merian = unique(df$Merian))


expanded_df <- merge(df, all_combinations, by = c("BLG", "Merian"), all = TRUE)
expanded_df[is.na(expanded_df$matches), "matches"] <- 0
df <- expanded_df

df <- df %>% group_by(BLG) %>% mutate('total_BLG' = sum(matches)) %>% ungroup()
df <- df %>% group_by(Merian) %>% mutate('total_Merian' = sum(matches)) %>% ungroup()

df$prop_BLG <- (df$matches / df$total_BLG) *100
df$prop_Merian <- (df$matches / df$total_Merian) *100

df_filt <- df %>% filter(prop_BLG > 5)

merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31')
df$Merian <- factor(df$Merian, levels=merian_order)

BLG_order <- c('A1axA1b', 'EaxEb',	'D',	'F',	'C1',	'G','H',	'K',	'N',	'L',	'M',	'B1',	'I',	'B2',	'P',	'J2',
                'O1',	'J1',	'O2',	'B3',	'QaxQb','A2',	'QcxQd','C2',	'R')
df$BLG <- factor(df$BLG, levels=BLG_order)

plot_prop_BLG <- ggplot(df, aes(x=Merian, y=BLG, fill=prop_BLG)) + geom_tile() + scale_fill_gradient(low="white", high="blue") + ggtitle("Proportion of BLG") +
  theme(axis.text.x = element_text(angle=90))
plot_prop_Merian <- ggplot(df, aes(x=Merian, y=BLG, fill=prop_Merian)) + geom_tile() + scale_fill_gradient(low="white", high="blue") + ggtitle('Proportion of Merian') +
  theme(axis.text.x = element_text(angle=90))

heatmap_plots <- plot_prop_BLG + plot_prop_Merian


df_stats <- df %>% select(c(BLG, Merian, matches))
data_wide <- spread(df_stats, Merian, matches)
data_wide <- data.frame(data_wide)
rownames(data_wide) <- data_wide$BLG
data_wide <- data_wide[,2:31]

Merian_df <- df %>% select('Merian', 'total_Merian') %>% distinct()
BLG_df <- df %>% select('BLG', 'total_BLG') %>% distinct()
Merian_df$stat_status <- "Significant"
false_merians <- c('M13', 'M18', 'M29', 'M31', 'M30','M20', 'M12', 'M16', 'M28', 'MZ', 'M15') # output from python simulation script
Merian_df <- Merian_df %>%  mutate(stat_status = ifelse(Merian %in% false_merians, "Non-significant", stat_status))
Merian_df <- Merian_df %>% mutate(stat_status = ifelse(total_Merian < 15, "Excluded", stat_status))
Merian_df$stat_status <- as.factor(Merian_df$stat_status, levels=c("Significant", "Non-significant", "Excluded"))

merian_boxplot <- ggplot(data=Merian_df, aes(x=Merian, y=total_Merian), fill="dimgrey") + 
  theme_classic()  +
  geom_bar(stat="identity", position="stack") +
  labs(y="Number orthologs per Merian", fill="Estimate of variance") +
  theme(axis.text.x= element_blank(),
        axis.title.x = element_blank()) 

var_legend <- cowplot::get_legend(merian_boxplot) 

BLG_boxplot <- ggplot(data=BLG_df, aes(x=BLG, y=total_BLG)) + 
  theme_classic()  +
  geom_bar(stat="identity", position="stack", fill="dimgrey") +
  labs(y="Number orthologues per BLG") +
  theme(axis.text.y= element_blank(),
        axis.title.y = element_blank()) + coord_flip()

plot_prop_Merian <- ggplot(df, aes(x=Merian, y=BLG, fill=prop_Merian)) + 
  geom_tile() + scale_fill_gradient2(low="white", high="#254482") +
  theme(axis.text.x = element_text(angle=90),
        legend.position = "left") +
  labs(fill = "Proportion of Merian element")

sup_fig <- plot_spacer() + merian_boxplot+theme(legend.position = "none")  + var_legend +  plot_spacer() + plot_prop_Merian + BLG_boxplot +
  plot_layout(widths = c(0.5, 2, 1), heights = c(1, 2))

ggsave(plot=sup_fig, 'Figures/Supplementary/Merians_vs_BLgs.pdf', device = 'pdf', width = 10, height = 6.6 , units = "in", limitsize=FALSE) #  save figure

# test for correlation overall between BLGs and Merians
# need to have simulate flag as not enough memory to run in full
test <- fisher.test(data_wide, simulate.p.value=TRUE)
test$p.value # 0.000499 - reject null hypothesis - there is a sig relationship between the two variables!
chisq.test(data_wide)$expected
# see may obs < 5
# if try:
chisq.test(data_wide)
# In chisq.test(data_wide) : Chi-squared approximation may be incorrect
# means theres obs <5 so should use Fisher'
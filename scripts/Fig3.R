### Figure 3
## 'Chromosome evolution in Lepidoptera'
## Wright et al., 2022
## cw22@sanger.ac.uk

library(reshape2) 
library(ggrepel) 
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(forcats)
library(patchwork)

## set input path
setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/')

## import data
freq_merians <- read.csv('../Results/Syngraph_vs_LFSF/final_analysis/011122_total_fusions_per_Merian.tsv', sep='\t', header=TRUE)
freq_merian_pairs <- read.csv('../Results/Syngraph_vs_LFSF/final_analysis/011122_total_fusions_per_pair_of_Merians.tsv', sep='\t', header=TRUE, skip=1)
chr_lengths <- read.csv('../Sup_tables/Sup_table_all_chr.tsv', sep='\t')[,c(1,2,3,18,21)]
assignments <- read.csv('../Sup_tables/Sup_table_all_chr.tsv', sep='\t')[,c(1,2,19,20)]
colnames(assignments) <- c('species', 'chr', 'status', 'assigned_ref_chr')
complex_spp <- c('Apeira_syringaria', 'Leptidea_sinapis', 'Brenthis_ino', 'Lysandra_bellargus', 'Lysandra_coridon', 'Melinaea_marsaeus_rileyi', 'Melinaea_menophilus', 'Operophtera_brumata', 'Philereme_vetulata', 'Pieris_brassicae', 'Pieris_napi', 'Pieris_rapae', 'Aporia_crataegi', 'Tinea_semifulvella')
assignments <- assignments[!assignments$species %in% complex_spp,] # these analyses do not include complex species
genome_sizes <- read.csv('../Sup_tables/Table_S1_Assembly_Information_040523.tsv', sep='\t', header=TRUE)[,c(1,9)]


## format data
names(freq_merian_pairs)[9] <- "M17+M20"
assignments[assignments  == "M17,M20"] <- "M17+M20"
colnames(genome_sizes) <- c('species', 'genome_size')
header <- colnames(freq_merian_pairs)
row.names(freq_merian_pairs) <- header

# filter data
trichoptera_species <- c('Glyphotaelius_pellucidus', 'Limnephilus_lunatus', 'Limnephilus_marmoratus', 'Limnephilus_rhombicus')
chr_lengths <- chr_lengths[!chr_lengths$species %in% trichoptera_species,]
LFSF <- LFSF[!LFSF$Node %in% c(trichoptera_species, 'n1'),]

################################################################################
## get counts of data - this isn't used for anything :-)
LFSF <- read.csv('../Results/Syngraph_vs_LFSF/final_analysis/mapped_fusions_fissions_noComplex_181022.tsv', sep='\t', header=TRUE)
LFSF_subsets <- read.csv('../Results/Syngraph_vs_LFSF/final_analysis/noComplex_181022_subset_fusions.tsv', sep='\t', header=TRUE)
LFSF_subsets <- LFSF_subsets[c(1:2)]
LFSF_subsets$remove <- "TRUE"
LFSF <- merge(LFSF, LFSF_subsets, by=c("Merians", "Node"), all.x=TRUE)  # this makes >1 entry for M17.M20 as >once in LFSF_subsets
# check what this dataframe is meant to have??
# it seems to have the pre-subset fusions, not post-subset i.e. M17,M20 not M17,M20,M21??

#LFSF[is.na(LFSF)] = "FALSE"
#LFSF_filt <- LFSF[LFSF$remove == "FALSE",]
#nrow(LFSF_filt) # 168
#Z_fusions <- filter(LFSF_filt, grepl('MZ', Merians))

#Fusion_counts <- LFSF %>% group_by(Merians) %>% count() 
#Filt_fusion_counts <- LFSF_filt %>% group_by(Merians) %>% count() 
#Z_counts <- Z_fusions %>% group_by(Merians) %>% count()
#sum(Z_counts$n) # 33 independent Z-fusions
#sum(Fusion_counts$n) # 199 total rearrangements
#sum(Filt_fusion_counts$n) # 168 completely independent fusions

################################################################################

# convert to matrix, then to dataframe
freq_merian_pairs <- as.matrix(freq_merian_pairs)
freq_merian_pairs <- melt(freq_merian_pairs)
colnames(freq_merian_pairs) <- c("Merian_1", "Merian_2", "Freq")
merian_order <- c("M29", "M31", "M25", "M30", "M27", "M28", "M26", "M24", "M19","M13", 
                  "M14", "M11", "M23", "M10", "M15", "M22", "M21", "M4", "M6", "M16", 
                  "M18", "M7", "M12", "M5", "M3", "M9", "M8", "M1", "M17+M20", "M2","MZ")
reverse_merian_order <- rev(merian_order)

freq_merian_pairs$Merian_1_f <- factor(freq_merian_pairs$Merian_1, levels = merian_order)
freq_merian_pairs$Merian_2_f <- factor(freq_merian_pairs$Merian_2, levels = merian_order)

## get prop_chr_length for ancestral chr only
prop_chr_length <- merge(chr_lengths, assignments, by=c('species', 'chr'), all.x = TRUE)
prop_chr_length_M17_M20 <- prop_chr_length[prop_chr_length$assigned_ref_chr == 'M17,M20',] %>% filter_all(any_vars(!is.na(.)))
prop_chr_length_M17_M20[prop_chr_length_M17_M20 == "fusion"] <- "ancestral"
prop_chr_length <- prop_chr_length[prop_chr_length$status == 'ancestral',] %>% filter_all(any_vars(!is.na(.))) # filter out fusions/fissions 
prop_chr_length <- rbind(prop_chr_length, prop_chr_length_M17_M20)
prop_chr_length <- prop_chr_length[!prop_chr_length$assigned_ref_chr %in% c('M17', 'M20'),] # removes unfused M17/M20 i.e. in M. aruncella
prop_chr_length$Merian_f = factor(prop_chr_length$assigned_ref_chr, levels=merian_order)


### Check Z-autosome vs autosome-autosome fusions ###
# Aim: check for differences between Z-autosome vs autosome-autosome fusions
df <- freq_merian_pairs[,c(1,2,3)]
Z_autosome_fusions <- df[df$Merian_1 == "MZ",] # can just take from Merian_1 as matrix was refleted i.e. M1,M2 same as M2,M1 count
autosome_autosome_fusions <- df[(df$Merian_1 != "MZ") & (df$Merian_2 != "MZ"),] 

autosome_autosome_fusions <- autosome_autosome_fusions %>% group_by(Merian_1) %>% 
  summarise(Freq = sum(Freq))
autosome_autosome_fusions$Merian_1 <- as.character(autosome_autosome_fusions$Merian_1)
Z_autosome_fusions$Merian_2 <- as.character(Z_autosome_fusions$Merian_2)
Z_autosome_fusions$Merian_f = factor(Z_autosome_fusions$Merian_2, levels=merian_order)
autosome_autosome_fusions$Merian_f = factor(autosome_autosome_fusions$Merian_1, levels=merian_order)
Z_autosome_fusions$Type <- "Z-autosome"
autosome_autosome_fusions$Type <- "Autosome-autosome"
Z_autosome_fusions <- Z_autosome_fusions[, 2:5]
colnames(Z_autosome_fusions) <- c("Merian_1", "Freq" ,"Merian_f", "Type")
Z_autosome_fusions$Freq <- as.numeric(Z_autosome_fusions$Freq)
Z_autosome_fusions[nrow(Z_autosome_fusions) + 1,] <- c("MZ", sum(Z_autosome_fusions$Freq), "MZ", "Z-autosome")
# Combine Z & autosome-autosome fusions
Z_vs_autosome_fusions <- rbind(Z_autosome_fusions, autosome_autosome_fusions)  
Z_vs_autosome_fusions$Freq<- as.numeric(Z_vs_autosome_fusions$Freq)

# set plotting parameters
x_text_size = 13
x_title_size = 25

# plot as a stacked barchart
Z_vs_autosome_barchart <- ggplot(data=Z_vs_autosome_fusions, aes(x=Merian_f, y=Freq, fill=Type)) + 
  theme_classic()  +
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=c("lightgrey", "Dimgrey")) + 
  labs(y="Number of fusion events", x="Merian element") +
  coord_flip() +   scale_x_discrete(position = "top") +
  theme(axis.text = element_text(size = x_text_size),
        legend.position = "left",
        axis.title=element_text(size=x_title_size))

# plot heatmap
heatmap_plot <- ggplot(freq_merian_pairs, aes(Merian_1_f, Merian_2_f, fill=Freq)) + 
  geom_tile() +   scale_fill_gradient2(low="navy", mid="white", high="red",guide="legend", 
                                       midpoint=0, limits=range(freq_merian_pairs$Freq)) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=x_text_size),
        axis.text.x = element_text(angle=90, hjust=0.95),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  guides(fill=guide_legend(title="Frequency"))

# make boxplots of length distribution
length_boxplot <- ggplot(prop_chr_length, aes(x=Merian_f, y=prop_length, colour=Merian_f))  + 
  geom_boxplot(outlier.shape = NA) + 
  theme_classic() +coord_flip()  +
  labs(y="Proportional chr length (%)", x=("Merian element")) +
  theme(legend.position = "none", 
        axis.title = element_text(size = x_title_size),
        axis.text = element_text(size=x_text_size)) 
length_boxplot

# get number of observations per merian
num_obs_df <- prop_chr_length %>% group_by(assigned_ref_chr) %>% summarise(num_obs = n()) %>% ungroup()
num_obs_df <- num_obs_df %>% arrange(factor(assigned_ref_chr, levels = rev(merian_order)))

length_boxplot <- ggplot(prop_chr_length, aes(x=Merian_f, y=prop_length, colour=Merian_f))  + 
  geom_boxplot(outlier.shape = NA) + 
  theme_classic() +coord_flip()  +
  labs(y="Proportional chr length (%)", x=("Merian element")) +
  theme(legend.position = "none", 
        axis.title = element_text(size = x_title_size),
        axis.text = element_text(size=x_text_size)) +
  scale_x_discrete(breaks=rev(merian_order), label=function(x) paste0(x, ' (n=', num_obs_df$num_obs, ')')) # label categories with Merian and number of obs


# combine plots in Fig3 and save
Fig3A <- length_boxplot + labs(subtitle = "A") +  theme(legend.position = 'none')
Fig3B_and_Fig3C <- heatmap_plot + labs(subtitle = "B") + Z_vs_autosome_barchart + labs(subtitle = "C") +
  plot_layout(widths = c(1.4, 1), guides = "collect") &  theme(legend.position = 'top')

Fig3 <- Fig3A + Fig3B_and_Fig3C +  plot_layout(ncol=2, widths=c(1,2.4))  & theme(plot.subtitle = element_text(size=15, face=2))
ggsave(plot=Fig3, '../Figures/Main_text/NEE_submission2/Fig3_061223.pdf', device = 'pdf', width = 20, height = 10, units = "in", limitsize=FALSE) #  save figure

# output plotted tsv tables to save as source data
fig3a_plotted_source_data <- merge(prop_chr_length, num_obs_df, by='assigned_ref_chr')
write.table(fig3a_plotted_source_data, file = "../../Chromosome_evolution_Lepidoptera_MS/data/prop_length_per_merian_per_species_141232.tsv", row.names=FALSE, sep="\t", quote = FALSE)
write.table(freq_merian_pairs, file = "../../Chromosome_evolution_Lepidoptera_MS/data/freq_of_fusion_per_pair_of_merians_141232.tsv", row.names=FALSE, sep="\t", quote = FALSE)
write.table(Z_vs_autosome_fusions, file = "../../Chromosome_evolution_Lepidoptera_MS/data/number_autosome_autosome_vs_autosome_sex_fusions_per_Merian_141223.tsv", row.names=FALSE, sep="\t", quote = FALSE)

# Test for correlation between avg_prop_length and frequency of fusions

## calculate the average prop_chr_length per Merian
mean_prop_length <- prop_chr_length %>% group_by(assigned_ref_chr) %>% summarise_at(vars(prop_length), list(name = mean))
colnames(mean_prop_length) <- c('Merian', 'avg_prop_length')

corr <- cor.test(x=freq_merians$Freq, y=freq_merians$avg_prop_length, method = 'spearman')

# plot supplementary figures - frequency of fusions per merian vs avg_prop_length
freq_merians <- merge(freq_merians, mean_prop_length, by='Merian') # only considering pairwise fusions

corr_plot <- ggplot(data=freq_merians, aes(x=avg_prop_length, y=Freq, label=Merian)) + geom_point(size=4) +  
  theme_classic() + geom_text_repel(size=6,min.segment.length = 0, seed = 42, box.padding = 0.5) +
  ylab('Number of fusion events') + xlab('Proportional chromosome length (%)') +
  theme(axis.text=element_text(size=16, color='black'), axis.title = element_text(size=20))

ggsave(plot=corr_plot, '../Figures/Supplementary/Avg_Merian_length_vs_freq_fusions.pdf', device = 'pdf', width = 12, height = 10, units = "in", limitsize=FALSE) #  save figure
### Figure 5
## 'Chromosome evolution in Lepidoptera'
## Wright et al., 2022
## cw22@sanger.ac.uk

library(dplyr)
library(ggplot2)
library(patchwork)

## set input path
setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/')

## import data
repeat_per_segment <- read.csv('../Data/lep_fusion_split2/final_analysis/combined_repeat_counts_per_chr_and_Merian.tsv', sep='\t', header=FALSE)
merian_info_per_segment <- read.csv('../Data/lep_fusion_split2/final_analysis/Fusion_segment_merian_info_adjusted.tsv', sep='\t', header=FALSE)
all_chromosomes <- read.csv('../Sup_tables/Sup_table_all_chr.tsv', sep='\t')  # read in all_chromosomes (including unassigned + complex species)

## format data
all_chromosomes$genome_repeat_density <- all_chromosomes$genome_repeat_density * 100
all_chromosomes$scaled_repeat_density <- all_chromosomes$scaled_repeat_density/100
colnames(repeat_per_segment) <- c('species', 'chr', 'start', 'stop', 'repeat_density')
colnames(merian_info_per_segment) <- c('species', 'chr', 'merian_segment','start', 'stop')
repeat_per_segment <- merge(repeat_per_segment, merian_info_per_segment, by=c("chr", "start", "stop","species"), all.x=TRUE)
repeat_per_segment$length <- repeat_per_segment$stop - repeat_per_segment$start
all_chromosomes$start <- 0
all_chromosomes$stop <- all_chromosomes$length

## identify fusion chromosomes
fusion_chr_df_complete<- all_chromosomes[all_chromosomes$status == "fusion",]
fusion_chr_df <- fusion_chr_df_complete[,c('species','chr','assigned_merian','status','start','stop')] # just take the columns needed for now

## make a lookup table that maps each species to its genome size
species2genome_size <- all_chromosomes[,c('species', 'genome_size')]
species2genome_size <- unique(species2genome_size)

## make a lookup table that maps each species to its average repeat density
species2genome_repeat_density <- all_chromosomes[,c('species', 'genome_repeat_density')]
species2genome_repeat_density <- unique(species2genome_repeat_density)

## make a lookup table that maps each Merian element to its average scaled repeat density 
# using averages_per_feature generated from 'Fig5_stats_on_feature_trends.R'
averages_per_feature <- read.csv('../Sup_tables/SupTable_average_feature_density_per_Merian.tsv', sep='\t', header =TRUE)
Merian2Average_scaled_repeat_density <- averages_per_feature[,c('assigned_merian', 'avg_scaled_repeat_density')]
colnames(Merian2Average_scaled_repeat_density) <- c("merian_segment", "avg_scaled_repeat_density")

## make a lookup table that maps each Merian element to its average proportional length
Merian2Average_length <- averages_per_feature[,c('assigned_merian', 'average_prop_length_per_Merian')]
colnames(Merian2Average_length) <- c("merian_segment", "avg_prop_length_per_Merian")

## make a lookup table that maps each chromosome to its proportional length
chr2prop_length <- fusion_chr_df_complete[,c("chr", "prop_length")]
colnames(chr2prop_length) <- c("chr", "complete_chr_prop_length")

# remove non-fused chromosomes from segmented_df (i.e. get rid of ancestral/split)
fusion_chr <- fusion_chr_df$chr
fusion_chr <- fusion_chr[!is.na(fusion_chr)] # remove NAs
repeat_per_segment <- repeat_per_segment[repeat_per_segment$chr %in% fusion_chr,]

# combine repeat information with Merian information
repeat_per_segment <- merge(repeat_per_segment, fusion_chr_df, all.x=TRUE, by=c('chr','species','start','stop'))
repeat_per_segment["status"][is.na(repeat_per_segment["status"])] <- "segment" # if part of chr, call a "segment"

# add in genome size and average repeat density info
repeat_per_segment <- merge(repeat_per_segment, species2genome_size, by="species", all.x=TRUE)
repeat_per_segment <- merge(repeat_per_segment, species2genome_repeat_density, by="species", all.x=TRUE)

# calculate prop_length & scaled_repeat_density
repeat_per_segment$prop_length <- repeat_per_segment$length / repeat_per_segment$genome_size
repeat_per_segment$scaled_repeat_density <- repeat_per_segment$repeat_density / repeat_per_segment$genome_repeat_density
repeat_per_segment$repeat_density <- repeat_per_segment$repeat_density*100

small_merians <- c("M29", "M31", "M25", "M30") # these are the four smallest merians on average

# make a dataframe for each small Merian, then combine together. Make a list of all such chromosomes
M29_df <- repeat_per_segment %>% filter(grepl('M29', assigned_merian))
M31_df <- repeat_per_segment %>% filter(grepl('M31', assigned_merian))
M25_df <- repeat_per_segment %>% filter(grepl('M25', assigned_merian))
M30_df <- repeat_per_segment %>% filter(grepl('M30', assigned_merian))

small_merian_fusions <- rbind(M29_df, M31_df,M25_df, M30_df)
small_merians_fusions_chr <- small_merian_fusions$chr

# find all Z-autosome fusions
Z_chr_df <- repeat_per_segment %>% filter(grepl('MZ', assigned_merian))
Z_chr <- unique(Z_chr_df$chr)

# make final dataframe for plotting - all autosome-autosome fusions, containing the 4 smallest Merians and two-part fusions only
small_merian_fusions <- repeat_per_segment[repeat_per_segment$chr %in% small_merians_fusions_chr,]
small_merian_fusions <- small_merian_fusions[!small_merian_fusions$chr %in% Z_chr,]

small_merian_fusions <- small_merian_fusions %>% group_by(chr) %>% mutate(n_row = n()) # count number entries per chr
small_merian_fusions <- small_merian_fusions[small_merian_fusions$n_row == '3',] # filter to just keep simple two-part fusions (i.e. original chr + two fusion entries)
complex_spp <- c('Apeira_syringaria', 'Leptidea_sinapis', 'Brenthis_ino', 'Lysandra_bellargus', 'Lysandra_coridon', 'Melinaea_marsaeus_rileyi', 'Melinaea_menophilus', 'Operophtera_brumata', 'Philereme_vetulata', 'Pieris_brassicae', 'Pieris_napi', 'Pieris_rapae', 'Aporia_crataegi', 'Tinea_semifulvella')
small_merian_fusions <- small_merian_fusions[!small_merian_fusions$species %in% complex_spp,] # for these analyses - preferable to not have complex spp

# for each fusion, calculate the 'large' and 'small' Merian involved ("relative_size")
segments_df <- small_merian_fusions[small_merian_fusions$status == "segment",]
complete_df <- small_merian_fusions[small_merian_fusions$status == "fusion",]
segments_df$relative_size <- "Large/Small"
complete_df$relative_size <- "Complete"
segments_df <- segments_df %>% group_by(chr) %>% mutate(max_segment_length = max(length))
segments_df$relative_size[segments_df$length == segments_df$max_segment_length] <- "Large"
segments_df$relative_size[segments_df$length != segments_df$max_segment_length] <- "Small"
small_merian_fusions <- rbind(segments_df, complete_df)

# combine raw data (observed) with Merian averages (expected)
small_merian_fusions <- merge(small_merian_fusions, Merian2Average_scaled_repeat_density, all.x=TRUE, by="merian_segment")
small_merian_fusions$scaled_repeat_density <- small_merian_fusions$scaled_repeat_density * 100
small_merian_fusions$prop_length = small_merian_fusions$prop_length*100
small_merian_fusions$repeat_delta <- (small_merian_fusions$avg_scaled_repeat_density - small_merian_fusions$scaled_repeat_density)*100

#avg_scaled_repeat_density = 0.9809362
# scaled_repeat_density = 1.0164027

small_merian_fusions <- merge(small_merian_fusions, Merian2Average_length, all.x=TRUE, by="merian_segment")
small_merian_fusions <- merge(small_merian_fusions, chr2prop_length, all.x=TRUE, by ="chr")
small_merian_fusions$length_delta <- small_merian_fusions$complete_chr_prop_length - small_merian_fusions$avg_prop_length_per_Merian

# conduct a two-tailed paired t-test to determine whether there's a difference in mean for small vs large segments
small_segment_deltas <- small_merian_fusions[small_merian_fusions$relative_size == "Small",]
large_segment_deltas <- small_merian_fusions[small_merian_fusions$relative_size == "Large",]
small_segment_deltas <- small_segment_deltas$repeat_delta
large_segment_deltas <- large_segment_deltas$repeat_delta
t.test(small_segment_deltas, large_segment_deltas, alternative = "two.sided", var.equal = FALSE)

## make plots
axis_text_size = 14

# example of a recent fusion: Agrochola_circellaris
target_spp_df <- all_chromosomes[all_chromosomes$species == "Agrochola_circellaris",]
target_spp_df <- target_spp_df[target_spp_df$prop_length > 0.02,]
target_spp_df <- target_spp_df[!grepl("MZ", target_spp_df$assigned_merian),] # remove MZ (regardless of fused or not)
target_spp_fusions <- small_merian_fusions[small_merian_fusions$species == "Agrochola_circellaris",]

# make copies of dataframe for source data
Fig5D_source_data_1 <- target_spp_df
Fig5D_source_data_2 <- target_spp_fusions

target_spp_fusions$tempvar <- "Recent fusion" # do this to add title as a strip later

Fig5D <- ggplot() + 
  geom_point(data=target_spp_df, color='black', fill="grey", shape=21, size=5, aes(x=prop_length, y=repeat_density)) +
  geom_point(data=target_spp_fusions, color='black', shape=21, size=5, aes(x=prop_length, y=repeat_density, fill=factor(relative_size))) +
  theme_bw() +
  scale_fill_manual(name='Relative size of chr',
                    values=c("black", "#B40F20","#46ACC8")) + 
  labs(x="Proportional chr length (%)", y="Repeat density (%)") +
  theme(axis.title = element_text(size=axis_text_size), 
        axis.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.text = element_text(size=axis_text_size))

Fig5D <- Fig5D + facet_grid(. ~ tempvar) +
  theme(strip.background = element_rect(fill="lightgrey"),
        strip.text = element_text(size=15, colour="black"))

# example of an old fusion: Aphantopus_hyperantus
target_spp_df <- all_chromosomes[all_chromosomes$species == "Aphantopus_hyperantus",]
target_spp_df <- target_spp_df[target_spp_df$prop_length > 0.05,]  # remove unplaced scaffs
target_spp_df <- target_spp_df[!grepl("MZ", target_spp_df$assigned_merian),] # remove MZ (regardless of fused or not)
target_spp_fusions <- small_merian_fusions[small_merian_fusions$species == "Aphantopus_hyperantus",]

# make copies of dataframe for source data
Fig5E_source_data_1 <- target_spp_df
Fig5E_source_data_2 <- target_spp_fusions

target_spp_fusions$tempvar <- "Old fusion" # do this to add title as a strip later

Fig5E <- ggplot() + 
  geom_point(data=target_spp_df, color='black', fill="grey", shape=21, size=5, aes(x=prop_length, y=repeat_density)) +
  geom_point(data=target_spp_fusions, color='black', shape=21, size=5, aes(x=prop_length, y=repeat_density, fill=factor(relative_size))) +
  theme_bw() +
  scale_fill_manual(name='Relative size of chr',
                    values=c("black", "#B40F20","#46ACC8")) + 
  labs(x="Proportional chr length (%)", y="Repeat density (%)") +
  theme(axis.title = element_text(size=axis_text_size), 
        axis.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.text = element_text(size=axis_text_size))

Fig5E <- Fig5E + facet_grid(. ~ tempvar) + # add title
  theme(strip.background = element_rect(fill="lightgrey"),
        strip.text = element_text(size=axis_text_size, colour="black"))

## make a boxplot of all large vs small Merian comparisons 
segments_only <- small_merian_fusions[small_merian_fusions$relative_size != "Complete",]
segments_only$tempvar <- "All fusions"
segments_only$relative_size <- factor(segments_only$relative_size , levels=c("Small", "Large"))

Fig5F <- ggplot(data=segments_only, aes(x=relative_size, y=repeat_delta, fill=relative_size)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#46ACC8", "#B40F20")) + 
  geom_jitter(color="black",alpha=0.5, size=2) +
  labs(x="Proportional length (%)", y=expression(Delta("Expected for Merian - observed repeat density"))) +
  theme_bw() + xlab("") +
  theme(axis.title = element_text(size=axis_text_size), 
        axis.text = element_text(size=axis_text_size),
        legend.position="none")   +
  geom_boxplot(outlier.shape = NA) # just plot outliers

Fig5F <- Fig5F + facet_grid(. ~ tempvar) + # add title
  theme(strip.background = element_rect(fill="lightgrey"),
        strip.text = element_text(size=axis_text_size, colour="black")) 

## Pierini
Pierini_spp <- c("Anthocharis_cardamines", "Pieris_rapae", "Pieris_napi", "Pieris_brassicae", "Aporia_crataegi")
Pierini <- all_chromosomes[all_chromosomes$species %in% Pierini_spp,]
Pierini <- Pierini[!grepl("CAJ", Pierini$chr),] # remove unplaced scaffolds (all have a 'CAJ' prefix)- one chr is in range of unplaced scaffs, thus use this instead
Pierini_noMZ <- Pierini[!grepl("MZ", Pierini$assigned_merian),] # remove MZ (regardless of fused or not)
Pierini_noMZ <- Pierini_noMZ %>% filter(prop_length > 0.2) # remove mito

## Lysandra
Lysandra_spp <- c("Lysandra_bellargus","Lysandra_coridon", "Polyommatus_icarus")
Lysandra <- all_chromosomes[all_chromosomes$species %in% Lysandra_spp,]
Lysandra <- Lysandra[!grepl("CAJ", Lysandra$chr),] # remove unplaced scaffolds
Lysandra_noMZ <- Lysandra[!grepl("MZ", Lysandra$assigned_merian),] # remove MZ (regardless of fused or not)
Lysandra_noMZ <- Lysandra_noMZ %>% filter(prop_length > 0.2) # remove mito

Pierini_noMZ$tempvar <- "Pierini"
Lysandra_noMZ$tempvar <- "Lysandra"

Fig5G <- ggplot(Pierini_noMZ, aes(x=prop_length, y=repeat_density)) + 
  geom_point(color='black', shape=21, size=4, aes(fill=factor(species))) + theme_bw() +
  scale_fill_manual(name='Pierini species',
                     values=c('Anthocharis_cardamines'='darkgrey', 'Aporia_crataegi'='green', 
                              'Pieris_rapae'='orange', 'Pieris_napi'='purple', 'Pieris_brassicae'='blue')) +
  labs(x="Proportional chr length (%)", y="Repeat density (%)") +
  theme(legend.position = "right",
        axis.title = element_text(size=axis_text_size), 
        axis.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.text = element_text(size=axis_text_size, face='italic'))

Fig5G <- Fig5G + facet_grid(. ~ tempvar) +
  theme(strip.background = element_rect(fill="lightgrey"),
        strip.text = element_text(size=axis_text_size, colour="black"))

Fig5H <- ggplot(Lysandra_noMZ, aes(x=prop_length, y=repeat_density)) + 
  geom_point(color='black', shape=21, size=4, aes(fill=factor(species))) + theme_bw() +
  scale_fill_manual(name='Lysandra species', values=c('Polyommatus_icarus'='darkgrey', 'Lysandra_coridon'='pink', 'Lysandra_bellargus'='lightblue')) +
  labs(x="Proportional chr length (%)", y="Repeat density (%)") +
  theme(legend.position = "right",
        axis.title = element_text(size=axis_text_size), 
        axis.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.text = element_text(size=axis_text_size, face='italic'))

Fig5H <- Fig5H + facet_grid(. ~ tempvar) +
  theme(strip.background = element_rect(fill="lightgrey"),
        strip.text = element_text(size=axis_text_size, colour="black"))

## make figure
design <- "
123
456
"  
Fig5D_to_Fig5H <- Fig5D + Fig5E + Fig5F + Fig5G + Fig5H + guide_area() + plot_layout(design=design, guides = "collect")

ggsave(plot=Fig5D_to_Fig5H, '../Figures/Fig5D-Fig5H.pdf', device = 'pdf', width = 13, height = 10, units = "in", limitsize=FALSE)


# output plotted tsv tables to save as source data
# Fig 5D - Agrochola_circellaris
cols_1_keep <- c('chr',	'assigned_merian','prop_length',	'repeat_density', 'status')
Fig5D_source_data_1 <- Fig5D_source_data_1 %>% select(all_of(cols_1_keep))
Fig5D_source_data_1$relative_size <- 'NA'
cols_2_keep <- c('chr',	'assigned_merian','prop_length',	'repeat_density', 'status', 'relative_size')
Fig5D_source_data_2 <- Fig5D_source_data_2 %>% select(all_of(cols_2_keep))
Fig5D_source_data <- rbind(Fig5D_source_data_1, Fig5D_source_data_2)
write.table(Fig5D_source_data, file = "../Chromosome_evolution_Lepidoptera_MS/data/Fig5d_recent_fusion_repeat_density_141232.tsv", row.names=FALSE, sep="\t", quote = FALSE)

# Fig 5E - Aphantopus_hyperantus
Fig5E_source_data_1 <- Fig5E_source_data_1 %>% select(all_of(cols_1_keep))
Fig5E_source_data_1$relative_size <- 'NA'
cols_2_keep <- c('chr',	'assigned_merian','prop_length',	'repeat_density', 'status', 'relative_size')
Fig5E_source_data_2 <- Fig5E_source_data_2 %>% select(all_of(cols_2_keep))
FigE_source_data <- rbind(Fig5E_source_data_1, Fig5E_source_data_2)
write.table(Fig5E_source_data, file = "../Chromosome_evolution_Lepidoptera_MS/data/Fig5e_recent_fusion_repeat_density_141232.tsv", row.names=FALSE, sep="\t", quote = FALSE)

# Fig 5F - observed_vs_expected_repeat_density_per_segment_141232
# data=segments_only, aes(x=relative_size, y=repeat_delta, fill=relative_size
Fig5F_source_data <- segments_only
cols_keep <- c('chr',	'merian_segment',	'species',	'relative_size',	'repeat_delta')
Fig5F_source_data <- Fig5F_source_data %>% select(all_of(cols_keep))
write.table(Fig5F_source_data, file = "../Chromosome_evolution_Lepidoptera_MS/data/observed_vs_expected_repeat_density_per_segment_141232", row.names=FALSE, sep="\t", quote = FALSE)
head(segments_only)
# Fig 5G (Pierini)
write.table(Pierini_noMZ, file = "../Chromosome_evolution_Lepidoptera_MS/data/Repeat_density_vs_prop_length_Pierini_141232.tsv", row.names=FALSE, sep="\t", quote = FALSE)

# Fig 5H (Lysandra)
write.table(Lysandra_noMZ, file = "../Chromosome_evolution_Lepidoptera_MS/data/Repeat_density_vs_prop_length_Lysandra_141232.tsv", row.names=FALSE, sep="\t", quote = FALSE)

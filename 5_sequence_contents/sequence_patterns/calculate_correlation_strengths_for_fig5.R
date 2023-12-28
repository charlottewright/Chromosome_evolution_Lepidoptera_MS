## Calculating strength of feature correlations
## 'Chromosome evolution in Lepidoptera'
## Wright et al., 2022
## cw22@sanger.ac.uk

library(dplyr)
library(tidyr)

## set input paths
setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/')

## load functions to assist with plotting
source('functions_feature_patterns_plots.R')

## import data
# pre-Hsara was:
#features <- read.csv('../Data/combined/summarised_features_per_chr_270223.tsv', sep='\t')
features <- read.csv('../Data/combined/summarised_features_per_chr_100523.tsv', sep='\t')
metadata <- read.csv('../Sup_tables/Table_S1_Assembly_Information_040523.tsv', sep='\t', header=TRUE)[,c(1,7,9)]
assignment_info <- read.csv('../Data/lep_fusion_split2/final_analysis/chr_assignments_all_spp_final_analysis.tsv', sep='\t', col.names = c('species', 'chr', 'status', 'merians'), header=FALSE)

assignment_info$chr[assignment_info$chr == "SUPER_10"] <- "OU860789.1"
assignment_info$chr[assignment_info$chr == "SUPER_11"] <- "OU860790.2"
assignment_info$chr[assignment_info$chr == "SUPER_12"] <- "OU860791.1"
assignment_info$chr[assignment_info$chr == "SUPER_13"] <- "OU860792.2"
assignment_info$chr[assignment_info$chr == "SUPER_14"] <- "OU860793.1"
assignment_info$chr[assignment_info$chr == "SUPER_15"] <- "OU860794.2"
assignment_info$chr[assignment_info$chr == "SUPER_16"] <- "OU860795.2"
assignment_info$chr[assignment_info$chr == "SUPER_17"] <- "OU860796.2"
assignment_info$chr[assignment_info$chr == "SUPER_18"] <- "OU860797.2"
assignment_info$chr[assignment_info$chr == "SUPER_19"] <- "OU860798.1"
assignment_info$chr[assignment_info$chr == "SUPER_20"] <- "OU860799.1"
assignment_info$chr[assignment_info$chr == "SUPER_Z"] <- "OU860800.1"
assignment_info$chr[assignment_info$chr == "SUPER_1"] <- "OU860780.1"
assignment_info$chr[assignment_info$chr == "SUPER_2"] <- "OU860781.1"
assignment_info$chr[assignment_info$chr == "SUPER_3"] <- "OU860782.2"
assignment_info$chr[assignment_info$chr == "SUPER_4"] <- "OU860783.2"
assignment_info$chr[assignment_info$chr == "SUPER_5"] <- "OU860784.1"
assignment_info$chr[assignment_info$chr == "SUPER_6"] <- "OU860785.2"
assignment_info$chr[assignment_info$chr == "SUPER_7"] <- "OU860786.1"
assignment_info$chr[assignment_info$chr == "SUPER_8"] <- "OU860787.2"
assignment_info$chr[assignment_info$chr == "SUPER_9"] <- "OU860788.1"

## merge features, metadata and assignment information
colnames(metadata) <- c('species', 'superfamily', 'genome_size')
features <- merge(features, metadata, by="species", all.x=TRUE) %>%
  merge(assignment_info, by=c("chr","species"),all.x=TRUE)

## format data
features$assigned_merian[features$assigned_merian == "M17,M20"] <- "M17+M20"
features$prop_length <- (features$length / features$genome_size)*100
features$length_mb <- (features$length)/1000000
features$cds_density <- features$cds_density*100
features$repeat_density <- features$repeat_density*100
features$genome_repeat_density <- features$genome_repeat_density*100

features$prop_single_copy <- (features$num_single_copy_genes / (features$num_multi_copy_genes + features$num_other_genes + features$num_single_copy_genes))*100
features$total_orthologs <- features$num_multi_copy_genes + features$num_other_genes + features$num_single_copy_genes
features <- transform(features, gc3_per = as.numeric(gc3_per))
features$scaled_gc3 <- features$gc3_per / features$genome_gc3_per

## before filtering any chromosome out, calculate average feature density per genome
features <- get_average_feature_density_per_species(features, "cds_density")
features <- get_average_feature_density_per_species(features, "prop_single_copy")
features <- get_average_feature_density_per_species(features, "gc")

features$scaled_cds_density <- features$cds_density / features$cds_density_scale_factor
features$scaled_ortholog_prop <- features$prop_single_copy / features$prop_single_copy_scale_factor
features$scaled_gc <- features$gc / features$gc_scale_factor
features$scaled_repeat_density <- features$repeat_density / features$genome_repeat_density

## filter to just keep chromosomes that have been assigned to Merians with confidence, save table
assigned_chr_features <- features[!is.na(features$assigned_merian),] 
write.table(assigned_chr_features, file='../Sup_tables/Sup_table_all_chr.tsv', quote=FALSE, sep='\t', row.names = FALSE)

complex_spp <- c('Apeira_syringaria', 'Leptidea_sinapis', 'Brenthis_ino', 'Lysandra_bellargus', 'Lysandra_coridon', 'Melinaea_marsaeus_rileyi', 'Melinaea_menophilus', 'Operophtera_brumata', 'Philereme_vetulata', 'Pieris_brassicae', 'Pieris_napi', 'Pieris_rapae', 'Aporia_crataegi', 'Tinea_semifulvella')
assigned_chr_features_noComplex <- assigned_chr_features[!assigned_chr_features$species %in% complex_spp,] # for these analyses - preferable to not have complex spp

write.table(assigned_chr_features_noComplex, file='../Sup_tables/Sup_table_all_assigned_chr_no_complex.tsv', quote=FALSE, sep='\t', row.names = FALSE)

## find superfamilies with >=5 species
filtered_superfamilies <- metadata %>% group_by(superfamily) %>% dplyr:::mutate(nSpecies = n()) %>%
  ungroup() %>% filter(nSpecies >= 5)
superfamilies_keep <- unique(filtered_superfamilies$superfamily) # find superfamilies with >=5 species

features$status[features$assigned_merian == "M17+M20"] <- "ancestral_Ditrysia"
ancestral_chromosomes <- features %>% filter(status %in% c("ancestral","ancestral_Ditrysia")) # only keep ancestral merians
ancestral_features_noMZ <- ancestral_chromosomes %>% filter(assigned_merian!= "MZ")
chr_to_keep <- ancestral_features_noMZ %>% group_by(species) %>%  dplyr::mutate(num_ancestral = n()) %>% ungroup() %>% filter(num_ancestral>=10) 
# note that the resulting df causes issues with ggplot2, instead use it to filter the original df 
chr_to_keep <- chr_to_keep$chr
ancestral_features_filt_noMZ <- ancestral_features_noMZ[ancestral_features_noMZ$chr %in% chr_to_keep,]
ancestral_MZ <- ancestral_chromosomes[ancestral_chromosomes$assigned_merian == "MZ",]
# note filtering for 10 obs makes no difference to number of species, but good to check

# calculate strength of correlation
species_list <- unique(ancestral_features_filt_noMZ$species)
cds_rows <- ancestral_features_filt_noMZ[!is.na(ancestral_features_filt_noMZ$cds_density),]
species_with_cds_info <- unique(cds_rows$species)

repeat_rows <- ancestral_features_filt_noMZ[!is.na(ancestral_features_filt_noMZ$repeat_density),]
species_with_repeats_info <- unique(repeat_rows$species)

repeat_cor_per_species <- data.frame(matrix(nrow = 0, ncol = length(1:3))) 
synteny_cor_per_species <- data.frame(matrix(nrow = 0, ncol = length(1:3)))
coding_cor_per_species <- data.frame(matrix(nrow = 0, ncol = length(1:3)))
orthologs_cor_per_species <- data.frame(matrix(nrow = 0, ncol = length(1:3)))
gc3_cor_per_species <- data.frame(matrix(nrow = 0, ncol = length(1:3)))
gc_cor_per_species <- data.frame(matrix(nrow = 0, ncol = length(1:3)))

colnames(repeat_cor_per_species) <- c('species', 'rep_corr', 'rep_pvalue')
colnames(synteny_cor_per_species) <- c('species', 'synteny_corr', 'synteny_pvalue')
colnames(coding_cor_per_species) <- c('species', 'cds_corr', 'cds_pvalue')
colnames(orthologs_cor_per_species) <- c('species', 'ortho_corr', 'ortho_pvalue')
colnames(gc3_cor_per_species) <- c('species', 'gc3_corr', 'gc3_pvalue')
colnames(gc_cor_per_species) <- c('species', 'gc_corr', 'gc_pvalue')

ancestral_features_filt_noMZ$repeat_density <- as.numeric(ancestral_features_filt_noMZ$repeat_density)

for (i in species_with_repeats_info){
  subset_data <- ancestral_features_filt_noMZ[ancestral_features_filt_noMZ$species == i,]
  cor_strength <- cor(subset_data$prop_length, subset_data$repeat_density, method = "spearman", use = "complete.obs")
  p_value <- cor.test(subset_data$prop_length, subset_data$repeat_density,method="spearman")$p.value
  repeat_cor_per_species[nrow(repeat_cor_per_species) + 1,] <- c(i, cor_strength, p_value)
}

for (i in species_list){
  subset_data <- ancestral_features_filt_noMZ[ancestral_features_filt_noMZ$species == i,]
  cor_strength <- cor(subset_data$prop_length, subset_data$synteny, method = "spearman", use = "complete.obs")
  p_value <- cor.test(subset_data$prop_length, subset_data$synteny, method="spearman")$p.value
  synteny_cor_per_species[nrow(synteny_cor_per_species) + 1,] <- c(i, cor_strength, p_value)
  cor_strength <- cor(subset_data$prop_length, subset_data$gc, method = "spearman", use = "complete.obs")
  p_value <- cor.test(subset_data$prop_length, subset_data$gc, method="spearman")$p.value
  gc_cor_per_species[nrow(gc_cor_per_species) + 1,] <- c(i, cor_strength, p_value)
}

for (i in species_with_cds_info){
  subset_data <- ancestral_features_filt_noMZ[ancestral_features_filt_noMZ$species == i,]
  cor_strength <- cor(subset_data$prop_length, subset_data$cds_density, method = "spearman", use = "complete.obs")
  p_value <- cor.test(subset_data$prop_length, subset_data$cds_density, method="spearman")$p.value
  coding_cor_per_species[nrow(coding_cor_per_species) + 1,] <- c(i, cor_strength, p_value)
  cor_strength <- cor(subset_data$prop_length, subset_data$prop_single_copy, method = "spearman", use = "complete.obs")
  p_value <- cor.test(subset_data$prop_length, subset_data$prop_single_copy, method="spearman")$p.value
  orthologs_cor_per_species[nrow(orthologs_cor_per_species) + 1,] <- c(i, cor_strength, p_value)
  cor_strength <- cor(subset_data$prop_length, subset_data$gc3_per, method = "spearman", use = "complete.obs")
  p_value <- cor.test(subset_data$prop_length, subset_data$gc3_per, method="spearman")$p.value
  gc3_cor_per_species[nrow(gc3_cor_per_species) + 1,] <- c(i, cor_strength, p_value)
}

repeat_cor_per_species <- transform(repeat_cor_per_species, rep_pvalue = as.numeric(rep_pvalue))
synteny_cor_per_species <- transform(synteny_cor_per_species, synteny_pvalue = as.numeric(synteny_pvalue))
coding_cor_per_species <- transform(coding_cor_per_species, cds_pvalue = as.numeric(cds_pvalue))
orthologs_cor_per_species <- transform(orthologs_cor_per_species, ortho_pvalue = as.numeric(ortho_pvalue))
gc3_cor_per_species <- transform(gc3_cor_per_species, gc3_pvalue = as.numeric(gc3_pvalue))
gc_cor_per_species <- transform(gc_cor_per_species, gc_pvalue = as.numeric(gc_pvalue))

nrow(repeat_cor_per_species %>% filter(rep_pvalue < 0.05)) # (Spearmans: 180/193 - 93%)
nrow(synteny_cor_per_species %>% filter(synteny_pvalue < 0.05)) # (Spearmans: 132/193 - 68%) 
nrow(coding_cor_per_species %>% filter(cds_pvalue < 0.05)) # (Spearmans: 1/184  0.5% -ve correlation; 33/184 - 18%) 
nrow(orthologs_cor_per_species %>% filter(ortho_pvalue < 0.05)) # (Spearmans: 174/184 - 95%)
nrow(gc3_cor_per_species %>% filter(gc3_pvalue < 0.05)) # (Spearman's: 93/184 - 50% total; 48% -ve correlation as 4 are +ve)
nrow(gc_cor_per_species %>% filter(gc_pvalue < 0.05)) # (Spearman's: 163/193 - 84%; all are -ve correlation)

correlations_table <- merge(gc_cor_per_species, repeat_cor_per_species, all.x=TRUE)
correlations_table <- merge(correlations_table, synteny_cor_per_species, all.x=TRUE)
correlations_table <- merge(correlations_table, coding_cor_per_species, all.x=TRUE)
correlations_table <- merge(correlations_table, orthologs_cor_per_species, all.x=TRUE)
correlations_table <- merge(correlations_table, gc3_cor_per_species, all.x=TRUE)


write.table(correlations_table, file='../Sup_tables/Sup_table_correlations_strength.tsv', quote=FALSE, sep='\t', row.names = FALSE)

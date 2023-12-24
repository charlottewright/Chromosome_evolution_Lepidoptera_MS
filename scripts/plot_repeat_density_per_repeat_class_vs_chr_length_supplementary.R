## Extended Data Fig. 5
## 'Chromosome evolution in Lepidoptera'
## Wright et al., 2022
## cw22@sanger.ac.uk


# This script makes the Extended Data Figure 5 and outputs the accompanying source data
# Relatioship between the mean density of each repeat class vs proportional length, per Meria element across the dataset

library(dplyr)
library(tidyverse)

repeats <- read.csv('sup_tables/combined_repeat_density_per_family_per_chr.tsv', sep='\t', header=FALSE)[,c(1,2,3,9)]
colnames(repeats) <- c('species', 'repeat_type', 'chr', 'repeat_prop')
repeats <- merge(repeats, ancestral_chromosomes, by=c('chr', 'species'))

for (i in(unique(repeats$repeat_type))) {
  
  subsetted_df <- repeats[repeats$repeat_type == i,]
  subsetted_df$repeat_prop_scaled <- subsetted_df$repeat_prop / subsetted_df$genome_repeat_density
  subsetted_df_stats <- get_average_and_sdev_per_Merian(subsetted_df, "repeat_prop_scaled")
  
  avg <- subsetted_df_stats[,1]
  sdev <-  subsetted_df_stats[,2]
  
  
pdf(paste0('../Chromosome_evolution_Lepidoptera/Figures/Supplementary/', i,'_density_per_Merian_vs_prop_length.pdf'), width=8, height=6)

layout <- plot(lengths_vector, avg,
               ylim=range(c(avg-sdev, avg+sdev)),
               pch=19, cex=2, xlab="Average proportional length of Merian element (%)", ylab="Average repeat density (%)",
               col = ifelse(1:length(lengths_vector) == 31, "red", "black"),
               main=i
)
  
#  draw errorbars using the 'arrows' function by specifying "arrowheads"
arrows(lengths_vector, avg-sdev, lengths_vector, avg+sdev, length=0.05, angle=90, code=3,
       ifelse(1:length(lengths_vector) == 31, "red", "black"))
#legend("topright", legend=c("Autosomal Merian element", "Merian element Z"), 
legend("topright", legend=c("Autosomal Merian element", "Merian element Z"),
       col=c("black", "red"),  lty=1:1, cex=1)
dev.off()
}

# output plotted tsv tables to save as source data
# Extended Data 5
source_data_E5 <- repeats
source_data_E5$repeat_prop_scaled <- source_data_E5$repeat_prop / source_data_E5$genome_repeat_density
cols_keep <- c('chr', 'species', 'assigned_merian', 'prop_length','repeat_type', 'repeat_density', 'genome_repeat_density','scaled_repeat_density')
source_data_E5 <- source_data_E5 %>% select(all_of(cols_keep))
head(source_data_E5)
write.table(source_data_E5, file = "../Chromosome_evolution_Lepidoptera_MS/data/scaled_repeat_density_per_repeat_class_per_Merian_151223.tsv", row.names=FALSE, sep="\t", quote = FALSE)

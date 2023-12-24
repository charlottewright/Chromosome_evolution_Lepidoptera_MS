
# This script was the make a source data file for a sup fifure that shows prop length of each Merian vs number of fusion events it was involved in.

setwd("/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts")

merians <- read.csv('../Chromosome_evolution_Lepidoptera_MS/sup_tables/TableS9_Merian_statistics.tsv', sep='\t', header=TRUE)
merians <- merians[,c(1:2)]

fusions <- read.csv('../Chromosome_evolution_Lepidoptera_MS/data/011122_total_fusions_per_Merian.tsv', sep='\t', header=TRUE)
colnames(fusions) <- c('Merian_element', 'Freq_fusions')

fusions <- merge(fusions, merians, by="Merian_element")
# output plotted tsv tables to save as source data
# Extended Data 3
write.table(fusions, file = "../Chromosome_evolution_Lepidoptera_MS/data/Num_fusions_per_Merian_vs_freq_fusions.tsv", row.names=FALSE, sep="\t", quote = FALSE)

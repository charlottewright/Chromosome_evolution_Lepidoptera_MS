
setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Data/genes/')
df <- read.csv('Melitaea_cinxia.gene_and_repeat_stats.tsv', sep='\t', header=TRUE)


library(ggplot2)
library(patchwork)

df$length <- df$length / 1000000
df$gene_length <- df$gene_length / 1000000
df$repeat_length <- df$repeat_length / 1000000
df$intron_length <- df$intron_length / 1000000
df$non_genic_length <- df$non_genic_length / 1000000

df_Z <- df[df$chr == "HG992209.1",]
df_autosomes <- df[df$chr != "HG992209.1",]

gene_count <- ggplot(data=df_autosomes, aes(x=length, y=gene_count)) + geom_point() + theme_classic() + labs(x="Chr length (Mb)", y="Gene count") +
  geom_smooth(method = "loess", se=FALSE, col="black") +
  geom_point(data=df_Z, aes(x=length, y=gene_count), color="red")
repeat_count <- ggplot(data=df_autosomes, aes(x=length, y=repeat_count)) + geom_point() + theme_classic() + labs(x="Chr length (Mb)", y="Repeat count") +
  geom_smooth(method = "loess", se=FALSE, col="black") +
  geom_point(data=df_Z, aes(x=length, y=repeat_count), color="red")
gene_span <- ggplot(data=df_autosomes, aes(x=length, y=gene_length)) + geom_point() + theme_classic() + labs(x="Chr length (Mb)", y="Gene span") +
  geom_smooth(method = "loess", se=FALSE, col="black") +
  geom_point(data=df_Z, aes(x=length, y=gene_length), color="red")
repeat_span <- ggplot(data=df_autosomes, aes(x=length, y=repeat_length)) + geom_point() + theme_classic() + labs(x="Chr length (Mb)", y="Repeat span") +
  geom_smooth(method = "loess", se=FALSE, col="black") +
  geom_point(data=df_Z, aes(x=length, y=repeat_length), color="red")
intron_span <- ggplot(data=df_autosomes, aes(x=length, y=intron_length)) + geom_point() + theme_classic() + labs(x="Chr length (Mb)", y="Intron span") +
  geom_smooth(method = "loess", se=FALSE, col="black") +
  geom_point(data=df_Z, aes(x=length, y=intron_length), color="red")
non_genic_span <- ggplot(data=df_autosomes, aes(x=length, y=non_genic_length)) + geom_point() + theme_classic() + labs(x="Chr length (Mb)", y="Non-genic span") +
  geom_smooth(method = "loess", se=FALSE, col="black") +
  geom_point(data=df_Z, aes(x=length, y=non_genic_length), color="red")

feature_plots <- gene_count + repeat_count + gene_span + repeat_span + intron_span + non_genic_span
feature_plots <- feature_plots + plot_annotation(tag_levels = "A")
ggsave(plot=feature_plots, "../../Figures/Response_to_reviewers/Melitaea_cinxia_gene_and_repeat_counts.png", width = 9, height = 7)  

# check for correlation strengths
corr_gene_count <- cor.test(x=df_autosomes$length, y=df_autosomes$gene_count, method = 'spearman')
corr_repeat_count <- cor.test(x=df_autosomes$length, y=df_autosomes$repeat_count, method = 'spearman')
corr_gene_span <- cor.test(x=df_autosomes$length, y=df_autosomes$gene_length, method = 'spearman')
corr_repeat_span <- cor.test(x=df_autosomes$length, y=df_autosomes$repeat_length, method = 'spearman')
corr_intron_span <- cor.test(x=df_autosomes$length, y=df_autosomes$intron_length, method = 'spearman')
corr_non_genic_span <- cor.test(x=df_autosomes$length, y=df_autosomes$non_genic_length, method = 'spearman')


corr_gene_count
corr_repeat_count
corr_gene_span
corr_repeat_span
corr_intron_span
corr_non_genic_span

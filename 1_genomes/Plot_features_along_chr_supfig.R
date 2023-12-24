

# Aim: SupFig - plot gc, repeat and coding density in 100kb windows, plus GC in 100 chunks per chr
setwd('/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/')


genes <- read.csv('../Data/combined/per_100kb/combined_DivByThree.cds_density_counts_per_100kb.tsv', sep='\t', col.names = c('species', 'chr', 'start', 'stop', 'gene_density'))
repeats <- read.csv('../Data/combined/per_100kb/combined_repeat_counts_per_100kb.tsv', sep='\t', col.names = c('species', 'chr', 'start', 'stop', 'repeat_density'))
gc <- read.csv('../Data/combined/per_100kb/combined_GC_per_100kb.tsv', sep='\t', col.names = c('species', 'chr', 'start', 'stop', 'gc_per'))
assignments <- read.csv('../Data/lep_fusion_split2/final_analysis/chr_assignments_all_spp_final_analysis.tsv', sep='\t')
gc_chunks <- read.csv('../Data/combined/combined_GC_per_chunk_100_per_chr.tsv', sep='\t', col.names=c('species', 'chr', 'start', 'stop', 'gc_per'))

# tidy data
colnames(assignments) <- c('species',	'query_chr', 'rearrangement_status',	'Merian')

# filter dataframes to just keep reps from superfamily cladogram ('Fig5_superfamily_cladogram.R)
# Species listed here in order of appearance in superfamily_cladogram in Fig5
superfam_reps <- c("Deilephila_porcellus", "Agriopis_aurantiaria", "Abrostola_tripartita",
                   "Acentria_ephemerella", "Anthocharis_cardamines", 
                   "Blastobasis_adustella", "Bembecia_ichneumoniformis", "Apotomis_betuletana")



# Acleris_emargana and Agonopterix_arenella don't have MZ
# Acleris_emargana - Tortricoidea, try Apotomis_betuletana instead (all tortricids have a MZ fusion)
# Agonopterix_arenella - Gelechioidea, try Blastobasis_adustella instead


gc <- gc %>% filter(species %in% superfam_reps) %>% 
  mutate(midpos = (start+stop)/2)  %>%
  select(c('species', 'chr', 'midpos', 'gc_per'))

repeats <- repeats %>% filter(species %in% superfam_reps) %>% 
  mutate(midpos = (start+stop)/2)  %>%
  select(c('species', 'chr', 'midpos', 'repeat_density'))

genes <- genes %>% filter(species %in% superfam_reps) %>% 
  mutate(midpos = (start+stop)/2)  %>%
  select(c('species', 'chr', 'midpos', 'gene_density'))

gc_chunks$start <- as.numeric(gc_chunks$start)
gc_chunks <- gc_chunks %>% filter(species %in% superfam_reps) %>% 
  mutate(midpos = (start+stop)/2)  %>%
  select(c('species', 'chr', 'midpos', 'gc_per'))

# Add chunk number to gc_chunks for plotting
gc_chunks <- gc_chunks %>% group_by(species, chr) %>% mutate(chunk_n = row_number())

# combine dataframes into one for plotting
feature_df <- merge(gc, repeats, by=c('species', 'chr', 'midpos')) %>% merge(genes, by =c('species', 'chr', 'midpos'))

# find which chr in M1 in each species
M1 <- assignments %>% filter(Merian == "M1" & species %in% superfam_reps)
M1_chr <- M1$query_chr
all_ancestral_chr <- assignments %>% filter(rearrangement_status == "ancestral" &  species %in% superfam_reps)
all_ancestral_chr <- all_ancestral_chr$query_chr
MZ_chr <- assignments %>% filter(Merian == "MZ" & species %in% superfam_reps)
MZ_chr <- MZ_chr$query_chr

# filter combined dataframe to just keep M1
feature_df <- feature_df %>% filter(chr %in% M1_chr)
feature_df$midpos <- feature_df$midpos / 1000000 # convert to Mb
gc_chunks <- gc_chunks %>% filter(chr %in% all_ancestral_chr)
feature_df <- feature_df %>% filter(gc_per < 0.42)
                    

# Set species as factor using order given in superfam_reps
feature_df$species <- factor(feature_df$species, levels = superfam_reps)
gc_chunks$species <- factor(gc_chunks$species, levels = superfam_reps)
gc_chunks_MZ <- gc_chunks %>% filter(chr %in% MZ_chr) # filter to just keep MZ


gc_plot <- ggplot(feature_df, aes(x=midpos, y=gc_per)) + geom_point(alpha=0.1, colour="#78B7C5") + facet_wrap(~species, ncol=1, scales="free_y") +
  geom_smooth(method='loess', span=0.8, colour='#3B9AB2')  + theme_pubr() + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(face="bold")) + 
  labs(x="Chr position (Mb)", y="GC %") +
  ggtitle('A')

repeat_plot <- ggplot(feature_df, aes(x=midpos, y=repeat_density)) + geom_point(alpha=0.1, colour="#78B7C5") + facet_wrap(~species, ncol=1, scales="free_y") +
  geom_smooth(method='loess', span=0.8, colour="#3B9AB2")  + theme_pubr() + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(face="bold")) + 
  labs(x="Chr position (Mb)", y="Repeat density %") +
  ggtitle('B')

coding_plot <- ggplot(feature_df, aes(x=midpos, y=gene_density)) + geom_point(alpha=0.1, colour="#78B7C5") + facet_wrap(~species, ncol=1, scales="free_y") +
  geom_smooth(method='loess', span=0.8, colour="#3B9AB2")  + theme_pubr() + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(face="bold")) +
  labs(x="Chr position (Mb)", y="Coding density %") +
  ggtitle('C')

gc_chunk_plot <- ggplot() + geom_point(data = gc_chunks, aes(x=chunk_n, y=gc_per), alpha=0.1, colour="#78B7C5") + facet_wrap(~species, ncol=1, scales="free_y") +
  geom_smooth(method='loess', span=0.8, colour="#3B9AB2")  + theme_pubr() + 
  theme(strip.background = element_blank(),
       strip.text.x = element_blank(),
       plot.title = element_text(face="bold")) +
  labs(x="Relative chr position", y="GC  %") +
  ggtitle('D')   + geom_point(data = gc_chunks_MZ, aes(x=chunk_n, y=gc_per), alpha=0.5, colour="#F21A00")

feature_plot <- gc_plot + repeat_plot + coding_plot + gc_chunk_plot + plot_layout(ncol=4)

ggsave(plot=feature_plot, '../Figures/Supplementary/SupFig_feature_patterns_along_chr.pdf', device='pdf', width=10, height=10)

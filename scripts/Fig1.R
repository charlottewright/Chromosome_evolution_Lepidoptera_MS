## Figure 1
## 'Chromosome evolution in Lepidoptera'
## Wright et al., 2022
## cw22@sanger.ac.uk

library(phytools)
library(treeio)
library(dplyr)
library(ggtree)
library(ggplot2) 
library(aplot)

## set input path
setwd("/Users/cw22/Documents/R_work/Chromosome_evolution_Lepidoptera/Scripts/")

## import data
tree <- read.newick('../Data/busco2phylo/supermatrix_LG_G4_rooted_including_Hydropsyche.treefile')
LFSF <- read.csv('../Results/Syngraph_vs_LFSF/final_analysis/mapped_fusions_fissions_noComplex_181022.tsv', sep='\t', header=TRUE)
metadata <- read.csv('../Sup_tables/Table_S1_Assembly_Information_040523.tsv', sep='\t', header=TRUE)[,c(1,6,7,8,9,12)]

## read in functions associated with manipulating and plotting trees
source("functions_tree_plots.R")

## format metadata
metadata$Genome_size <- as.numeric(as.character(metadata$Genome_size))/ 1000000 # convert to numeric to have a continuous variable on axes
metadata$Species <- as.character(metadata$Species)
metadata$Chr_number <- as.numeric(metadata$Chr_number)

# remove trichopteran species and species that failed QC from metadata and tree 
tips_to_drop <- c("Limnephilus_lunatus", "Limnephilus_marmoratus", "Limnephilus_rhombicus", "Glyphotaelius_pellucidus", "Hydropsyche_tenuis","Cnaphalocrocis_medinalis", "Zerene_cesonia") # Remove "Cnaphalocrocis_medinalis" from dataframe - not in tree! (lacks M30)
metadata <- metadata[!metadata$Species %in% tips_to_drop, ]

for (i in 1:length(tips_to_drop)){
  tree <- drop.tip(tree, tips_to_drop[i])
}

# remove rearrangements that involve trichopterans
LFSF <- LFSF[!grepl("Limnephilus", LFSF$Node),]
LFSF <- LFSF[!grepl("Limnephilus", LFSF$Tips),] 
LFSF <- LFSF[!grepl("Glyphotaelius", LFSF$Node),]
LFSF <- LFSF[!grepl("Glyphotaelius", LFSF$Tips),] 

# check all species are in both the metadata and in the tree
a = metadata$Species
b = unique(a)
c = tree$tip.label
d = setdiff(c, b)
d # should be empty

## set levels for famililes and superfamilies
superfamily_levels <- c("Noctuoidea", "Geometroidea", "Bombycoidea", "Drepanoidea", 
                        "Pyraloidea", "Papilionoidea", "Gelechioidea", "Carposinoidea", 
                        "Pterophoroidea", "Sesioidea", "Zygaenoidea", "Cossoidea", "Tortricoidea", 
                        "Yponomeutoidea", "Tineoidea", "Micropterigoidea")

family_levels <- c("Noctuidae", "Erebidae", "Nolidae", "Notodontidae", "Geometridae", 
                   "Sphingidae","Bombycidae", "Lasiocampidae", "Drepanidae",
                   "Crambidae", "Pyralidae", "Nymphalidae", "Lycaenidae", "Pieridae",
                   "Hesperiidae", "Papilionidae", "Blastobasidae", "Depressariidae",
                   "Pterophoridae", "Carposinidae", "Sesiidae", "Zygaenidae", "Cossidae",
                   "Tortricidae", "Ypsolophidae", "Plutellidae", "Tineidae", "Micropterigidae")

metadata$Superfamily <- factor(metadata$Superfamily, levels = superfamily_levels)
metadata$Family <- factor(metadata$Family, levels = family_levels)

## specify aesthetics for plotting
colour_palette <- c('Dimgrey', 'darkgrey', '#f6ae2d', '#f26419', '#86bbd8') 

col1 <- colour_palette[1]  # colour1 for genome size / chr number plots
col2 <- colour_palette[2] # colour2 for genome size / chr number plots
colors_family <- c(col1, col2) %>% rep(14)

default_brancolour <- 'grey25'
rearranged_colour <- colour_palette[3]
complex_colour <- colour_palette[4]
ancient_fusion_colour <- colour_palette[5]


## prepare lists of nodes to plot
list_nodes <- LFSF$Node # get all nodes with a rearrangement 
list_nodes <- unique(list_nodes) # remove dups
list_spp <- Filter(function(x) !any(grepl("^n", x)), list_nodes) # get list of tips with a rearrangement
list_internal_nodes <- Filter(function(x) any(grepl("^n", x)), list_nodes) # get list of internal nodes with a rearrangement
list_internal_nodes <- Filter(function(x) !any(grepl("n6$", x)), list_internal_nodes) # remove n6

## get internal nodes with rearrangements
internal_nodes_wrt_tree <- list()
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Aricia_agestis', 'Aricia_artaxerxes'))) #n227 = 298
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Stenoptilia_bipunctidactyla', 'Marasmarcha_lunaedactyla'))) #n78 = 256
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Hedya_salicella', 'Apotomis_turbidana', 'Apotomis_betuletana', 'Notocelia_uddmanniana', 'Epinotia_nisella', 'Pammene_fasciana','Grapholita_molesta', 'Leguminivora_glycinivorella', 'Cydia_splendana'))) #n22 = 232
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Ypsolopha_sequella', 'Ypsolopha_scabrella'))) # "n20"    = 227
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Lycaena_phlaeas', 'Glaucopsyche_alexis', 'Celastrina_argiolus', 'Plebejus_argus', 'Polyommatus_icarus', 'Cyaniris_semiargus', 'Aricia_agestis', 'Aricia_artaxerxes'))) #n106 = 293 
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Eilema_sororcula', 'Eilema_depressum'))) #n305 = 378
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Agonopterix_subpropinquella', 'Agonopterix_arenella', 'Carcina_quercana', 'Blastobasis_lacticolella', 'Blastobasis_adustella'))) # n45 = 250 
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Manduca_sexta', 'Mimas_tiliae', 'Laothoe_populi', 'Hemaris_fuciformis', 'Deilephila_porcellus', 'Hyles_vespertilio', 'Hyles_euphorbiae'))) # n156 = 353
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Thera_britannica', 'Chloroclysta_siterata'))) # n303 = 346
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Thymelicus_sylvestris', 'Ochlodes_sylvanus', 'Hesperia_comma'))) # n110 = 304
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Papilio_machaon', 'Iphiclides_podalirius'))) # n64 = 307
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Pararge_aegeria', 'Erebia_aethiops' ,'Erebia_ligea', 'Hipparchia_semele', 'Melanargia_galathea', 'Lasiommata_megera'))) # n143 = 270
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Dendrolimus_punctatus', 'Dendrolimus_kikuchii'))) # n138 = 359
internal_nodes_wrt_tree <- append(internal_nodes_wrt_tree, getMRCA(tree, c('Acleris_emargana', 'Acleris_sparsana','Pandemis_cinnamomeana', 'Apotomis_betuletana', 'Apotomis_turbidana', 'Hedya_salicella', 'Cydia_splendana','Leguminivora_glycinivorella','Grapholita_molesta','Pammene_fasciana', 'Epinotia_nisella','Notocelia_uddmanniana'))) # n17 = 229

## Find species that possess internal node rearrangements
internal_nodes <- dplyr::filter(LFSF, grepl("^n",Node))
list_spp_with_internal_nodes <- internal_nodes$Tips

## find species that possess species-specific rearrangements
tips_with_fusions <- data.frame()
for (i in 1:length(list_spp)) {
    node_number <- (which_node(list_spp[[i]]))
    new_row <- data.frame(species = list_spp[i], node = node_number)
    tips_with_fusions <- rbind(tips_with_fusions, new_row)
}

# create a dataframe of the default color for branches
branchColor <- data.frame(node = 1:Nnode2(tree), colour = default_brancolour)
tips_with_fusions$color = rearranged_colour # Now change colour to 'orange' for TIPS with tip-specific rearrangements
branchColor <- merge(branchColor, tips_with_fusions, by="node", all.x=TRUE)
branchColor <- branchColor[,c(1,4)]
branchColor[["color"]][is.na(branchColor[["color"]])] <- default_brancolour

rownames(branchColor) <- branchColor$node # needed to make index same as node numbers
internal_nodes_2_colour <- c('288', '284', '285', '286', '287', '288', '289', '290', '291', '242', '241', '243', '344', '345', '346', '347', '348', '295', '261', '262', '263', '264', '265', '266', '220', '221', '222', '223', '224', '225', '226', '227', '228', '229')
tip_list = internal_nodes[1, "Tips"]
tip_list <- gsub(" ", "", tip_list)
tip_list <- as.list(strsplit(tip_list, ",")[[1]])

# make a list of all tips of an internal node with a fusion
for (row in 2:nrow(internal_nodes)) {
  tip_spp <- internal_nodes[row, "Tips"]
  tip_spp <- gsub(" ", "", tip_spp)
  tip_spp <- as.list(strsplit(tip_spp, ",")[[1]])
  tip_list <- append(tip_list, tip_spp)
}

tip_list <- unique(tip_list) # remove dups

for (i in 1:length(tip_list)){
  matching_node <- which_node(tip_list[i]) # first look up the node
  branchColor[matching_node, 2] <- rearranged_colour
  
}

for (i in 1:length(internal_nodes_2_colour)){
  branchColor[internal_nodes_2_colour[i], 2] <- rearranged_colour
}

## colour the "complex" species
complex_spp <- c('Apeira_syringaria', 'Leptidea_sinapis', 'Brenthis_ino', 'Lysandra_bellargus', 'Lysandra_coridon', 'Melinaea_marsaeus_rileyi', 'Melinaea_menophilus', 'Operophtera_brumata', 'Philereme_vetulata', 'Pieris_brassicae', 'Pieris_napi', 'Pieris_rapae', 'Aporia_crataegi', 'Tinea_semifulvella')

complex_tips <- data.frame()
for (i in 1:length(complex_spp)){
  node_num <- which(tree$tip.label==complex_spp[i])
  branchColor[node_num, 2] <- complex_colour
}

getMRCA(tree, c('Pieris_brassicae', 'Pieris_napi', 'Pieris_rapae', 'Aporia_crataegi')) # 254
complex_internal_nodes <- c('255', '256') # internal nodes to node 254 i.e. Pierids

for (i in 1:length(complex_internal_nodes)){
  branchColor[complex_internal_nodes[i], 2] <- complex_colour
}

### Find nodes where shared complex rearrangements are present ###
getMRCA(tree, c('Pieris_brassicae', 'Pieris_napi', 'Pieris_rapae', 'Aporia_crataegi')) # 254
getMRCA(tree, c('Lysandra_bellargus', 'Lysandra_coridon')) # 290 
getMRCA(tree, c('Melinaea_marsaeus_rileyi', 'Melinaea_menophilus'))  # 282
  
### Add family/superfamily labels ###
family <- metadata[,c(1,4)]
superfamily <- metadata[,c(1,3)]
family <- family %>% group_by(Family) %>% mutate(families_grouped = paste0(Species, collapse = ",")) # get list of tips belonging to each family
family <- family[,2:3] %>% distinct() # remove duplicate rows. Just have [family, spp]
superfamily <- superfamily %>% group_by(Superfamily) %>% mutate(superfamilies_grouped = paste0(Species, collapse = ",")) # get list of tips belonging to each family
superfamily <- superfamily[,2:3] %>% distinct() # remove duplicate rows. Just have [family, spp]

family_labels <- data.frame()
superfamily_labels <- data.frame()

# First do for family
for (row in 1:nrow(family)) {
  family_name <- family[row, "Family"]
  species <- family[row, "families_grouped"]
  species <- as.character(species)
  species <- as.list(strsplit(species, ",")[[1]])
  num_species <- length(species)
  if (num_species > 1) {
    node_number <-  getMRCA(tree, unlist(species))
    new_row <- data.frame(family = family_name, node = node_number)
    family_labels <- rbind(family_labels, new_row)
  }
    else{
      node_number <- which_node(species)
      new_row <- data.frame(family = family_name, node = node_number)
      family_labels <- rbind(family_labels, new_row)
}}

# Then for superfamily
for (row in 1:nrow(superfamily)) {
  superfamily_name <- superfamily[row, "Superfamily"]
  species <- superfamily[row, "superfamilies_grouped"]
  species <- as.character(species)
  species <- as.list(strsplit(species, ",")[[1]])
  num_species <- length(species)
  if (num_species > 1) {
    node_number <-  getMRCA(tree, unlist(species))
    new_row <- data.frame(superfamily = superfamily_name, node = node_number)
    superfamily_labels <- rbind(superfamily_labels, new_row)
  }
  else{
    node_number <- which_node(species)
    new_row <- data.frame(superfamily = superfamily_name, node = node_number)
    superfamily_labels <- rbind(superfamily_labels, new_row)
  }}


# plot tree and annotate with family labels
# colour branches according to rearrangement status and add add circles to denote nodes with shared rearrangements
annotated_family_tree <-ggtree(tree, size = 0.5) %<+%
  branchColor +
  aes(colour = I(color)) +
  geom_tippoint(size = 0.1, color = default_brancolour) +
  #  geom_tiplab(size = 2, color = default_brancolour) +
  geom_treescale(x = 0, y = 15, width=0.1) +
  #geom_treescale(x = 0, y = -5, width = 0.1, fontsize = 3, linesize = 0.5)  +  # was y=200
  geom_cladelabel(node=373, label="Noctuidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=360, label="Erebidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=351, label="Notodontidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=343, label="Sphingidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=341, label="Bombycidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=314, label="Geometridae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=310, label="Drepanidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=306, label="Pyralidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=300, label="Crambidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=297, label="Papilionidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=292, label="Hesperiidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=283, label="Lycaenidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=258, label="Nymphalidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=251, label="Pieridae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=245, label="Pterophoridae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=243, label="Blastobasidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=240, label="Depressariidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=232, label="Sesiidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=219, label="Tortricidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=215, label="Ypsolophidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=213, label="Tineidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=148, label="Nolidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=138, label="Lasiocampidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=33, label="Carposinidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=27, label="Zygaenidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=26, label="Cossidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=4, label="Plutellidae", color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=1, label="Micropterigidae", color="black", align=TRUE, fontsize=2) +
  geom_point2(aes(subset=(node==288)), shape=21, size=2, fill=rearranged_colour, colour=default_brancolour) + 
  geom_point2(aes(subset=(node==246)), shape=21, size=2, fill=rearranged_colour)  + 
  geom_point2(aes(subset=(node==222)), shape=21, size=2, fill=rearranged_colour, colour=default_brancolour) +  #specify black outline as else is orange as is a branch that's orange
  geom_point2(aes(subset=(node==217)), shape=21, size=2, fill=rearranged_colour) + 
  geom_point2(aes(subset=(node==283)), shape=21, size=2, fill=rearranged_colour) + 
  geom_point2(aes(subset=(node==368)), shape=21, size=2, fill=rearranged_colour) + 
  geom_point2(aes(subset=(node==240)), shape=21, size=2, fill=rearranged_colour) +  
  geom_point2(aes(subset=(node==343)), shape=21, size=2, fill=rearranged_colour) + 
  geom_point2(aes(subset=(node==336)), shape=21, size=2, fill=rearranged_colour) +  
  geom_point2(aes(subset=(node==294)), shape=21, size=2, fill=rearranged_colour) + 
  geom_point2(aes(subset=(node==297)), shape=21, size=2, fill=rearranged_colour) + 
  geom_point2(aes(subset=(node==260)), shape=21, size=2, fill=rearranged_colour) + 
  geom_point2(aes(subset=(node==349)), shape=21, size=2, fill=rearranged_colour) + 
  geom_point2(aes(subset=(node==219)), shape=21, size=2, fill=rearranged_colour) + 
  geom_point2(aes(subset=(node==254)), shape=21, size=2, fill=complex_colour, colour=default_brancolour) + # nodes with shared complex rearrangements
  geom_point2(aes(subset=(node==290)), shape=21, size=2, fill=complex_colour, colour=default_brancolour) + 
  geom_point2(aes(subset=(node==282)), shape=21, size=2, fill=complex_colour, colour=default_brancolour) 

# plot tree and annotate with superfamily labels
# colour branches according to rearrangement status and add add circles to denote nodes with shared rearrangements
annotated_superfamily_tree <-ggtree(tree, size = 0.5) %<+%
  branchColor +
  aes(colour = I(color)) +
  geom_tippoint(size = 0.1, color = default_brancolour) +
  #  geom_tiplab(size = 2, color = default_brancolour) +
  geom_treescale(x = 0, y = 15, width=0.1) +  
  geom_cladelabel(node=350, label=toupper("Noctuidae"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=299, label=toupper("Pyraloidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=219, label=toupper("Tortricoidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=240, label=toupper("Gelechioidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=314, label=toupper("Geometroidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=248, label=toupper("Papilionoidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=232, label=toupper("Sesioidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=341, label=toupper("Bombycoidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=33, label=toupper("Carposinoidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=245, label=toupper("Pterophoroidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=310, label=toupper("Drepanoidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=1, label=toupper("Micropterigoidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=215, label=toupper("Yponomeutoidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=213, label=toupper("Tineoidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=26, label=toupper("Cossoidea"), color="black", align=TRUE, fontsize=2) +
  geom_cladelabel(node=27, label=toupper("Zygaenoidea"), color="black", align=TRUE, fontsize=2) +
  geom_point2(aes(subset=(node==212)), shape=21, size=2, fill=ancient_fusion_colour) + # ancient fusion label
  geom_point2(aes(subset=(node==288)), shape=21, size=2, fill=rearranged_colour, colour=default_brancolour) + 
  geom_point2(aes(subset=(node==246)), shape=21, size=2, fill=rearranged_colour) + 
  geom_point2(aes(subset=(node==222)), shape=21, size=2, fill=rearranged_colour, colour=default_brancolour) +  #specify black outline as else is orange as is a branch that's orange
  geom_point2(aes(subset=(node==217)), shape=21, size=2, fill=rearranged_colour) +
  geom_point2(aes(subset=(node==283)), shape=21, size=2, fill=rearranged_colour) + 
  geom_point2(aes(subset=(node==368)), shape=21, size=2, fill=rearranged_colour) +
  geom_point2(aes(subset=(node==240)), shape=21, size=2, fill=rearranged_colour) +
  geom_point2(aes(subset=(node==343)), shape=21, size=2, fill=rearranged_colour) +
  geom_point2(aes(subset=(node==336)), shape=21, size=2, fill=rearranged_colour) +
  geom_point2(aes(subset=(node==294)), shape=21, size=2, fill=rearranged_colour) +
  geom_point2(aes(subset=(node==297)), shape=21, size=2, fill=rearranged_colour) +
  geom_point2(aes(subset=(node==260)), shape=21, size=2, fill=rearranged_colour) +
  geom_point2(aes(subset=(node==349)), shape=21, size=2, fill=rearranged_colour) +
  geom_point2(aes(subset=(node==219)), shape=21, size=2, fill=rearranged_colour) +
  geom_point2(aes(subset=(node==254)), shape=21, size=2, fill=complex_colour, colour=default_brancolour) + # nodes with shared complex rearrangements
  geom_point2(aes(subset=(node==290)), shape=21, size=2, fill=complex_colour, colour=default_brancolour) + 
  geom_point2(aes(subset=(node==282)), shape=21, size=2, fill=complex_colour, colour=default_brancolour)

## plot chr_numbers
chr_number_plot <- ggplot(metadata, aes(Species, Chr_number)) + geom_col(aes(fill=Family), colour="white") +
  coord_flip() + theme_tree2() + theme(legend.position='none') + ylab("Haploid number (n)") +
  scale_fill_manual(values=colors_family) + ylim(0, 90)

## plot genome_sizes
genome_size_plot <- ggplot(metadata, aes(Species, Genome_size)) + geom_col(aes(fill=Family),colour="white") + 
  coord_flip() + theme_tree2() + theme(legend.position='none') + ylab("Genome size (Mb)") +
  scale_fill_manual(values=colors_family) +ylim(0, 2300) 

## plot tree with family labels alongside chr_numbers barchart                                                                                                                                       
Fig1A_and_Fig1B <- chr_number_plot %>% insert_left(annotated_family_tree, width=5) 
## plot tree with superfamily labels alongside genome_size barchart                                                                                                                                       
Fig1C <- genome_size_plot %>% insert_left(annotated_superfamily_tree, width=5)

## make histogram of chr_number vs no.species 
chr_number_histogram <- ggplot(metadata, aes(x=Chr_number)) + geom_histogram(binwidth = 1, fill=col2, boundary=0) + 
  theme_classic() + theme(axis.title = element_text(size = 20)) + scale_x_continuous(breaks = c(0, seq(0, 90, 5))) +
  xlab("Haploid chromosome number") + ylab("Number of species") 

## make histogram of genome_size vs no.species 
genome_size_histogram <- ggplot(metadata, aes(x=Genome_size))  + geom_histogram(binwidth = 50, fill=col2, boundary=0) + 
  theme_classic() + theme(axis.title = element_text(size = 20)) + scale_x_continuous(breaks = c(0, seq(0, 2500,250)), expand=c(0,0)) +
  xlab("Genome size (Mb)") + ylab("Number of species")


## save plots
ggsave(plot=Fig1A_and_Fig1B, "../Figures/Fig1A_and_Fig1B.pdf", width = 15, height = 9, units = "in", limitsize=FALSE)
ggsave(plot=Fig1C, "../Figures/Fig1C.pdf", width = 15, height = 9, units = "in", limitsize=FALSE)
ggsave(plot=genome_size_histogram, '../Figures/Supplementary/Histogram_genome_size.png', width = 9, height = 5, units = "in", limitsize=FALSE)
ggsave(plot=chr_number_histogram, '../Figures/Supplementary/Histogram_chromosome_number.png', width = 9, height = 5, units = "in", limitsize=FALSE)

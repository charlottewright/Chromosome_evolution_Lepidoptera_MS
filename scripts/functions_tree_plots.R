# Functions associated with manipulating and plotting trees

which_node <- function(tip_spp){ # lookup node number in tree for a given species
  node_number <- which(tree$tip.label==tip_spp)
  return(node_number)
}

filter_tree <- function(full_tree, tips_to_keep){ # filter tree to only keep a set of tips
  all_tips <- full_tree$tip.label
  tips_to_drop <- setdiff(all_tips, tips_to_keep)
  tree <- full_tree
  for (i in 1:length(tips_to_drop)){
    tree <- drop.tip(tree, tips_to_drop[i])
  } 
  return(tree)
}

# Michael Matschiner, 2015-07-09

# Load libraries.
library(ape)

# Get the command line arguments.
args <- commandArgs(trailingOnly = TRUE)
input_tree_file_name <- args[1]
output_tree_file_name <- args[2]
list_file_name <- args[3]

# Read the tree.
tree <- read.nexus(input_tree_file_name)

# Read the list of taxa to keep.
keep <- read.table(list_file_name)

# Remove all tips not in the list.
tree <- drop.tip(tree, tree$tip.label[-match(keep[,1], tree$tip.label)])

# Write the pruned tree to file.
write.tree(tree, file=output_tree_file_name)

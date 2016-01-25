# Michael Matschiner, 2015-07-07

# Load required libraries.
library(ape)

# Get the command lines arguments.
args <- commandArgs(trailingOnly = TRUE)
input_tree_file_name <- args[1]
tip_name <- args[2]
output_tree_file_name <- args[3]

# Read all trees.
tree <- read.nexus(input_tree_file_name)

# Prune the specified tip from all trees.
tree <- drop.tip(tree,tip=tip_name)

# Write the pruned tree to file.
write.tree(tree,output_tree_file_name)
# Michael Matschiner, 2015-07-07

# Load required libraries.
library(phytools)

# Get the command lines arguments.
args <- commandArgs(trailingOnly = TRUE)
input_tree_file_name <- args[1]
tip_name <- args[2]
output_tree_file_name <- args[3]

# Read all trees.
trees <- read.nexus(input_tree_file_name)

# Prune the specified tip from all trees.
trees <- lapply(trees,drop.tip,tip=tip_name)
class(trees) <- "multiPhylo"

# Write the pruned tree to file.
write.tree(trees,output_tree_file_name)
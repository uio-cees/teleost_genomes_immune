# Michael Matschiner, 2015-07-29

# Load libraries.
library(ape)

# Get the command line arguments.
args <- commandArgs(trailingOnly = TRUE)
input_tree_file_name <- args[1]
output_tree_file_name <- args[2]

# Read the tree.
tree <- read.nexus(input_tree_file_name)

# Remove all ENSEMBL species.
tree <- drop.tip(tree, "Gasacu")
tree <- drop.tip(tree, "Takrub")
tree <- drop.tip(tree, "Tetnig")
tree <- drop.tip(tree, "Orenil")
tree <- drop.tip(tree, "Orylat")
tree <- drop.tip(tree, "Poefor")
tree <- drop.tip(tree, "Xipmac")

# Write the pruned tree to file.
write.nexus(tree, file=output_tree_file_name, translate=F)

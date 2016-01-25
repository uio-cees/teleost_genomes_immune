# Michael Matschiner, 2015-07-10

# Load libraries.
library(ape)

# Get the command line arguments.
args <- commandArgs(trailingOnly = TRUE)
beast2_bins_input_tree_file_name <- args[1]
beast2_bins_output_tree_file_name <- args[2]
raxml_bins_input_tree_file_name <- args[3]
raxml_bins_output_tree_file_name <- args[4]
beast2_single_input_tree_file_name <- args[5]
beast2_single_output_tree_file_name <- args[6]
raxml_single_input_tree_file_name <- args[7]
raxml_single_output_tree_file_name <- args[8]

# Root and ladderize all trees.
tree <- read.tree(beast2_bins_input_tree_file_name)
tree <- root(tree, outgroup=c("Danrer"), node=NULL, resolve.root=FALSE)
tree <- ladderize(tree, right=FALSE)
write.tree(tree, beast2_bins_output_tree_file_name)
tree <- read.tree(raxml_bins_input_tree_file_name)
tree <- root(tree, outgroup=c("Danrer"), node=NULL, resolve.root=FALSE)
tree <- ladderize(tree, right=FALSE)
write.tree(tree, raxml_bins_output_tree_file_name)
tree <- read.nexus(beast2_single_input_tree_file_name)
tree <- root(tree, outgroup=c("Danrer"), node=NULL, resolve.root=FALSE)
tree <- ladderize(tree, right=FALSE)
write.tree(tree, beast2_single_output_tree_file_name)
tree <- read.tree(raxml_single_input_tree_file_name)
tree <- root(tree, outgroup=c("Danrer"), node=NULL, resolve.root=FALSE)
tree <- ladderize(tree, right=FALSE)
write.tree(tree, raxml_single_output_tree_file_name)

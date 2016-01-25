# Michael Matschiner, 2015-07-09

# Load libraries.
library(phytools)
require(methods)

# Get the command line arguments.
args <- commandArgs(trailingOnly = TRUE)
tree_file_name <- args[1]
table_file_name <- args[2]
log_copy_number_phenogram_file_name <- args[3]

# Read the tree.
tree <- read.tree(tree_file_name)

# Read the copy number table.
table <- read.table(table_file_name)
log_copy_numbers <- table$V2
names(log_copy_numbers) <- table$V1

# Produce the phenogram for log copy numbers.
pdf(log_copy_number_phenogram_file_name, height=7, width=7)
phenogram(tree, log_copy_numbers, spread.labels=TRUE, spread.cost=c(1,0), link=20, offset=1, fsize=0.5)
dev.off()

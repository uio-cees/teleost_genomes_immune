# Michael Matschiner, 2015-07-17

# Load the phytools library.
library(phytools)
require(methods)

# Get the command lines arguments.
args <- commandArgs(trailingOnly = TRUE)
tree_file_name <- args[1]
table_file_name <- args[2]
plot_file_name <- args[3]

# Read the phylogenetic tree.
phy <- read.tree(tree_file_name)

# Read the table with trait information.
table <- read.table(table_file_name,header=T,row.names = "Species")

# Open the output pdf.
pdf(plot_file_name, width=7, height=7)

# Produce a phylomorphospace plot.
phylomorphospace(phy, table, xlab="log(MHC I copy number)", ylab="Net diversification")

# Turn off device.
dev.off()
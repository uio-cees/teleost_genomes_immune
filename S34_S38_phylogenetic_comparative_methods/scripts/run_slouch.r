# Michael Matschiner, 2015-07-17.

# Load the slouch source code.
source("source_slouch.r")

# Get the command line arguments.
args <- commandArgs(trailingOnly = TRUE)
table_file_name <- args[1]
plot_file_name <- args[2]

# Read the slouch input table.
table <- read.table(table_file_name)
attach(table)

# Define the range of values for h (the phylogenetic half-life) and vy (the stationary variance).
h <- seq(11,200,length.out=190)
vy <- seq(0.11,2.0,length.out=190)

# Open a pdf file for the slouch plot.
pdf(plot_file_name, width=7, height=7)

# Run slouch's function model.fit.
model.fit(ancestor, time, h, vy, trait, trait.me, fixed.fact=regime, support=1)

# Cleanup.
detach(table)
dev.off()
# Michael Matschiner, 2015-09-22

# Load libraries.
library(diversitree)

# Set file names.
tree_file_name <- "unresolved.tre"
state_file_name <- "resolved.txt"
unresolved_table_file_name <- "unresolved.txt"
resultsfile_name <- "results.txt"

# Read the tree.
tree <- read.tree(tree_file_name)

# Read the table with states of resolved clades.
resolved_table <- read.table(state_file_name,header=T)
resolved <- resolved_table$state
names(resolved) <- resolved_table$tip.label

# Read the table with information for unresolved clades.
unresolved_table <- read.table(unresolved_table_file_name,header=T)
unresolved <-as.data.frame(unresolved_table)

# Create a likelihood function for this data.
lik <- make.bisse(tree, resolved, unresolved=unresolved)
lik <- set.defaults(lik, root = ROOT.GIVEN, root.p = c(1, 0))

# Add constraints to the model.
lik_more_constrained <- constrain(lik, lambda0 ~ lambda1, mu0 ~ mu1)
lik_less_constrained <- constrain(lik, mu0 ~ mu1)

# Select a good starting point for parameter optimization for the more constrained model.
start_more_constrained <- c(0.08, 0.01, 0.01, 0.001)
names(start_more_constrained) <- c("lambda1", "mu1", "q01", "q10")

# Find the maximum likelihood parameter combination for the more constrained model.
mle_more_constrained <- find.mle(lik_more_constrained, start_more_constrained, method="subplex")
l1 <- sprintf(paste(mle_more_constrained$par.full, sep=" "))
l2 <- sprintf(paste(mle_more_constrained$lnLik, sep=" "))

# Use the maximum likelihood parameter combination for the more constrained model as the starting point for the less constrained model.
start_less_constrained <- c(mle_more_constrained$par[1],mle_more_constrained$par[1],mle_more_constrained$par[2],mle_more_constrained$par[3],mle_more_constrained$par[4])
mle_less_constrained <- find.mle(lik_less_constrained, start_less_constrained, method="subplex")
l3 <- sprintf(paste(mle_less_constrained$par.full, sep=" "))
l4 <- sprintf(paste(mle_less_constrained$lnLik, sep=" "))

# Write results to file.
fileConn<-file(resultsfile_name)
writeLines(c(l1,l2,l3,l4), fileConn)
close(fileConn)

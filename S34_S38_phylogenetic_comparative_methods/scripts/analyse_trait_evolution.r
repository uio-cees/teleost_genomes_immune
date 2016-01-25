# Michael Matschiner, 2015-07-15

# Load required libraries.
library(pGLS)
library(geiger)
library(caper)
library(picante)
library(surface)
library(igraph)
library(phytools)
library(bayou)
library(OUwie)

# Get the command line arguments.
args <- commandArgs(trailingOnly = TRUE)
tree_file_name <- args[1]
data_file_name <- args[2]
plot_file_name <- args[3]

# Open the output pdf.
pdf(plot_file_name, width=7, height=7)

# Read the pruned phylogenetic tree.
tree <- read.tree(tree_file_name)
# plot.phylo(tree)
# nodelabels()

# Read the table of natural-log-transformed copy number estimates.
data <- read.table(data_file_name, row.names=1)
traits <- data[[1]]
names(traits) <- rownames(data)
stderrs <- data[[2]]
names(stderrs) <- rownames(data)

# # Use the phylosignal function of the picante package to test for phylogenetic signal, measured with the K statistic.
# K ~ 1.0 suggests Brownian motion, K < 1 indicates less resemblance among relatives than expected under Brownian motion,
# and K > 1 indicates more resemblance among relatives than expected under Brownian motion (Blomberg et al., 2003, Evolution, 57:717-745).
cat("\n----------------------------------------picante----------------------------------------\n\n")
phylosignal(traits, tree)

# Use the phylosig function of the phytools package to test for phylogenetic signal, measured with Pagel's lambda.
# Pagel's lambda is 0 <= lambda <= 1, with lambda ~ 0 indicating no phylogenetic signal, and lambda ~ 1 indicating as much
# phylogenetic signal as expected under Brownian motion (Pagel, 1999, Nature, 401:877-884).
cat("\n----------------------------------------phytools----------------------------------------\n\n")
phylosig(tree, traits, method="lambda", test=TRUE)

# Compare the fit of standard single-regime models, using geiger.
cat("\n----------------------------------------geiger----------------------------------------\n\n")
bm.model <- fitContinuous(tree, traits, model="BM")
white.model <- fitContinuous(tree, traits, model="white")
ou.model <- fitContinuous(tree, traits, model="OU")
eb.model <- fitContinuous(tree, traits, model="EB")
results <- c(bm.model$opt$aicc,white.model$opt$aicc,ou.model$opt$aicc,eb.model$opt$aicc)
names(results) <- c("bm","white_noise","ou","eb")
results

# Use the surface package to identify OU regimes in the phylogeny.
cat("\n----------------------------------------surface----------------------------------------\n\n")
surface_tree <- nameNodes(tree)
surface_data <- data
surface_data$V3 <- NULL
olist <- convertTreeData(surface_tree, surface_data)
otree <- olist[[1]]
.pardefault <- par()
plot(otree, node.names=TRUE, cex=0.5)
par(.pardefault)
odata <- olist[[2]]
fwd <- surfaceForward(otree, odata, aic_threshold = -1, exclude = 0, verbose = FALSE, plotaic=FALSE, sample_shifts=TRUE)
k1 <- length(fwd)
bwd <- surfaceBackward(otree, odata, starting_model = fwd[[k1]], aic_threshold = 0, only_best = TRUE, verbose = FALSE, plotaic = FALSE)
k2 <- length(bwd)
surfaceSummary(fwd, bwd)
surfaceTreePlot(surface_tree, cex=0.5, bwd[[k2]], labelshifts = TRUE)
surfaceAICPlot(fwd, bwd)

# Start surface with a OU model that has two regime shifts, one within Percomorphaceae, and one within Gadiformes.
cat("\n----------------------------------------surface: ou2----------------------------------------\n\n")
ou2_fwd <- startingModel(otree, odata, shifts=c("58"="b","31"="c"))
k3 <- length(ou2_fwd)
ou2_bwd <- surfaceBackward(otree, odata, starting_model = ou2_fwd[[k3]], aic_threshold = 0, only_best = TRUE, verbose = FALSE, plotaic = FALSE)
k4 <- length(ou2_bwd)
surfaceSummary(ou2_fwd, ou2_bwd)
surfaceTreePlot(surface_tree, cex=0.5, ou2_bwd[[k4]], labelshifts = TRUE)
surfaceAICPlot(ou2_fwd, ou2_bwd)

# Start surface with a OU model with multiple potential shifts.
cat("\n----------------------------------------surface: ou3----------------------------------------\n\n")
ou3_fwd <- startingModel(otree, odata, shifts=c("58"="b","31"="c","5"="d"))
k5 <- length(ou3_fwd)
ou3_bwd <- surfaceBackward(otree, odata, starting_model = ou3_fwd[[k5]], aic_threshold = 0, only_best = TRUE, verbose = FALSE, plotaic = FALSE)
k6 <- length(ou3_bwd)
surfaceSummary(ou3_fwd, ou3_bwd)
surfaceTreePlot(surface_tree, cex=0.5, ou3_bwd[[k6]], labelshifts = TRUE)
surfaceAICPlot(ou3_fwd, ou3_bwd)

# Start surface with a OU model with multiple potential shifts.
cat("\n----------------------------------------surface: ou4----------------------------------------\n\n")
ou4_fwd <- startingModel(otree, odata, shifts=c("58"="b","31"="c","17"="d","5"="e"))
k7 <- length(ou4_fwd)
ou4_bwd <- surfaceBackward(otree, odata, starting_model = ou4_fwd[[k7]], aic_threshold = 0, only_best = TRUE, verbose = FALSE, plotaic = FALSE)
k8 <- length(ou4_bwd)
surfaceSummary(ou4_fwd, ou4_bwd)
surfaceTreePlot(surface_tree, cex=0.5, ou4_bwd[[k8]], labelshifts = TRUE)
surfaceAICPlot(ou4_fwd, ou4_bwd)

# Use bayou to identify OU regimes in the phylogeny.
cat("\n----------------------------------------bayou----------------------------------------\n\n")

# Prepare a prior using bayou's make.prior function. This also produces a plot of all prior probabilities.
prior <- make.prior(tree, dists=list(dalpha="dlnorm", dsig2="dlnorm", dsb="dsb", dk="cdpois",
	dtheta="dnorm"), param=list(dalpha=list(meanlog=-5, sdlog=2),
	dsig2=list(meanlog=-1, sdlog=2), dk=list(lambda=6, kmax=200),
	dsb=list(bmax=Inf,prob=tree$edge.length), dtheta=list(mean=mean(traits), sd=2)))

# Use six plots per figure.
par(mfrow=c(2,3))

# Run a first bayou mcmc chain.
fit1 <- bayou.mcmc(tree, traits, SE=stderrs, model="OU", prior, ngen=2000000, new.dir=TRUE, plot.freq=400000, ticker.freq=200000)

# Report mcmc run statistics.
fit1

# Summarize the first bayou mcmc chain.
chain1 <- load.bayou(fit1, save.Rdata=FALSE, cleanup=TRUE)
chain1 <- set.burnin(chain1, 0.3)
out <- summary(chain1)

# Reset to use one plots per figure.
par(mfrow=c(1,1))

# Produce a tree figure for the first bayou mcmc chain.
plotSimmap.mcmc(tree, chain1, burnin=0.3, circle=TRUE, fsize=0.4)

# Produce a phenogram for the first bayou mcmc chain.
phenogram.density(tree, traits, chain=chain1, burnin=0.3, pp.cutoff=0.1)

# Run a second bayou mcmc chain.
fit2 <- bayou.mcmc(tree, traits,SE=stderrs, model="OU", prior, ngen=2000000, new.dir=TRUE, plot.freq=NULL, ticker.freq=200000)

# Summarize the second bayou mcmc chain.
chain2 <- load.bayou(fit2, save.Rdata=FALSE, cleanup=FALSE)
chain2 <- set.burnin(chain2, 0.3)

# Assess convergence by comparison of the two chains, using Gelman's R-statistic for the three parameters lnL, alpha, and sig2.
RlnL <- gelman.R("lnL", chain1=chain1, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain1, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain1, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
L1 <- Lposterior(chain1,tree, burnin=0.3)
L2 <- Lposterior(chain2,tree, burnin=0.3)
plot(L1$pp,L2$pp, xlim=c(0,1), ylim=c(0,1), xlab="Chain 1", ylab="Chain 2")
curve(1*x, add=TRUE, lty=2)

# Turn off device.
dev.off()

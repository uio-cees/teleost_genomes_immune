# Michael Matschiner, 2015-07-06.

# Load required libraries.
library(BAMMtools)
library(coda)
require(methods)

# Get the command lines arguments.
args <- commandArgs(trailingOnly = TRUE)
tree_file_name <- args[1]
mcmc_out_file_name <- args[2]
event_data_file_name <- args[3]
prior_probs_file_name <- args[4]
rates_file_name <- args[5]
post_probs_file_bame <- args[6]
netdiv_best_plot_file_name <- args[7]
lambda_best_plot_file_name <- args[8]
mu_best_plot_file_name <- args[9]
set_plot_of_4_file_name <- args[10]
set_plot_of_9_file_name <- args[11]
bayes_factors_plot_file_name <- args[12]
cohort_plot_file_name <- args[13]
net_div_tree_file_name <- args[14]
branch_file_name <- args[15]
branch2_file_name <- args[16]

# Check convergence.
convergence <- FALSE
mcmcout <- read.csv(mcmc_out_file_name, header=T)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
if (effectiveSize(postburn$logLik) > 200 && effectiveSize(postburn$N_shifts) > 200) {
	convergence <- TRUE
}

# Load the tree the event data, and the results from simulations using the prior only.
if (convergence) {
	tree <- read.tree(tree_file_name)
	edata <- getEventData(tree, eventdata = event_data_file_name, burnin=0.1)
	priordata <- read.table(prior_probs_file_name,header=T,sep=",")
	priorshifts <- getBranchShiftPriors(tree, priordata)

	# Calculate mean and 95% HPDs of clade net diversification rates.
	background_rates <- getCladeRates(edata, node=80, nodetype="exclude")
	#background_rates <- getCladeRates(edata, node=87, nodetype="exclude")
	background_netdiv <- background_rates$lambda - background_rates$mu
	background_netdiv_mean <- mean(background_netdiv)
	background_netdiv_hpd <- HPDinterval(as.mcmc(background_netdiv))
	background_netdiv_hpd_lower <- background_netdiv_hpd[1]
	background_netdiv_hpd_upper <- background_netdiv_hpd[2]
	background_lambda <- background_rates$lambda
	background_lambda_mean <- mean(background_lambda)
	background_lambda_hpd <- HPDinterval(as.mcmc(background_lambda))
	background_lambda_hpd_lower <- background_lambda_hpd[1]
	background_lambda_hpd_upper <- background_lambda_hpd[2]
	background_mu <- background_rates$mu
	background_mu_mean <- mean(background_mu)
	background_mu_hpd <- HPDinterval(as.mcmc(background_mu))
	background_mu_hpd_lower <- background_mu_hpd[1]
	background_mu_hpd_upper <- background_mu_hpd[2]
	gadiform_rates <- getCladeRates(edata, node=88)
	#gadiform_rates <- getCladeRates(edata, node=95)
	gadiform_netdiv <- gadiform_rates$lambda - gadiform_rates$mu
	gadiform_netdiv_mean <- mean(gadiform_netdiv)
	gadiform_netdiv_hpd <- HPDinterval(as.mcmc(gadiform_netdiv))
	gadiform_netdiv_hpd_lower <- gadiform_netdiv_hpd[1]
	gadiform_netdiv_hpd_upper <- gadiform_netdiv_hpd[2]
	gadiform_lambda <- gadiform_rates$lambda
	gadiform_lambda_mean <- mean(gadiform_lambda)
	gadiform_lambda_hpd <- HPDinterval(as.mcmc(gadiform_lambda))
	gadiform_lambda_hpd_lower <- gadiform_lambda_hpd[1]
	gadiform_lambda_hpd_upper <- gadiform_lambda_hpd[2]
	gadiform_mu <- gadiform_rates$mu
	gadiform_mu_mean <- mean(gadiform_mu)
	gadiform_mu_hpd <- HPDinterval(as.mcmc(gadiform_mu))
	gadiform_mu_hpd_lower <- gadiform_mu_hpd[1]
	gadiform_mu_hpd_upper <- gadiform_mu_hpd[2]
	percomorphaceae_rates <- getCladeRates(edata, node=122)
	#percomorphaceae_rates <- getCladeRates(edata, node=129)
	percomorphaceae_netdiv <- percomorphaceae_rates$lambda - percomorphaceae_rates$mu
	percomorphaceae_netdiv_mean <- mean(percomorphaceae_netdiv)
	percomorphaceae_netdiv_hpd <- HPDinterval(as.mcmc(percomorphaceae_netdiv))
	percomorphaceae_netdiv_hpd_lower <- percomorphaceae_netdiv_hpd[1]
	percomorphaceae_netdiv_hpd_upper <- percomorphaceae_netdiv_hpd[2]
	percomorphaceae_lambda <- percomorphaceae_rates$lambda
	percomorphaceae_lambda_mean <- mean(percomorphaceae_lambda)
	percomorphaceae_lambda_hpd <- HPDinterval(as.mcmc(percomorphaceae_lambda))
	percomorphaceae_lambda_hpd_lower <- percomorphaceae_lambda_hpd[1]
	percomorphaceae_lambda_hpd_upper <- percomorphaceae_lambda_hpd[2]
	percomorphaceae_mu <- percomorphaceae_rates$mu
	percomorphaceae_mu_mean <- mean(percomorphaceae_mu)
	percomorphaceae_mu_hpd <- HPDinterval(as.mcmc(percomorphaceae_mu))
	percomorphaceae_mu_hpd_lower <- percomorphaceae_mu_hpd[1]
	percomorphaceae_mu_hpd_upper <- percomorphaceae_mu_hpd[2]
	l0 <- sprintf("")
	l1 <- sprintf("Net diversification")
	l2 <- sprintf("Background: %.3f (%.3f-%.3f)", background_netdiv_mean, background_netdiv_hpd_lower, background_netdiv_hpd_upper)
	l3 <- sprintf("Gadiformes: %.3f (%.3f-%.3f)", gadiform_netdiv_mean, gadiform_netdiv_hpd_lower, gadiform_netdiv_hpd_upper)
	l4 <- sprintf("Percomorphaceae: %.3f (%.3f-%.3f)", percomorphaceae_netdiv_mean, percomorphaceae_netdiv_hpd_lower, percomorphaceae_netdiv_hpd_upper)
	l5 <- sprintf("Speciation")
	l6 <- sprintf("Background: %.3f (%.3f-%.3f)", background_lambda_mean, background_lambda_hpd_lower, background_lambda_hpd_upper)
	l7 <- sprintf("Gadiformes: %.3f (%.3f-%.3f)", gadiform_lambda_mean, gadiform_lambda_hpd_lower, gadiform_lambda_hpd_upper)
	l8 <- sprintf("Percomorphaceae: %.3f (%.3f-%.3f)", percomorphaceae_lambda_mean, percomorphaceae_lambda_hpd_lower, percomorphaceae_lambda_hpd_upper)
	l9 <- sprintf("Extinction")
	l10 <- sprintf("Background: %.3f (%.3f-%.3f)", background_mu_mean, background_mu_hpd_lower, background_mu_hpd_upper)
	l11 <- sprintf("Gadiformes: %.3f (%.3f-%.3f)", gadiform_mu_mean, gadiform_mu_hpd_lower, gadiform_mu_hpd_upper)
	l12 <- sprintf("Percomorphaceae: %.3f (%.3f-%.3f)", percomorphaceae_mu_mean, percomorphaceae_mu_hpd_lower, percomorphaceae_mu_hpd_upper)
	fileConn<-file(rates_file_name)
	writeLines(c(l1,l2,l3,l4,l0,l5,l6,l7,l8,l0,l9,l10,l11,l12), fileConn)
	close(fileConn)

	# Calculate the posterior probabilities of numbers of rate shifts and write these to file.
	post_probs <- table(postburn$N_shifts) / nrow(postburn)
	write.table(post_probs,post_probs_file_bame,quote=F,col.names=c("N_shifts","prop"), row.names=F)

	# Produce plots for the single best shift configuration.
	best <- getBestShiftConfiguration(edata, prior=priorshifts)
	pdf(netdiv_best_plot_file_name, height=7, width=7)
	par(pin=c(7,7))
	best_plot <- plot.bammdata(best, spex="netdiv", lwd=2, pal=c("#268bd2","#6c71c4","#d33682","#dc322f","#cb4b16","#eac115"),labels=T,cex=0.5)
	addBAMMlegend(best_plot, corners = c(10, 12, 30, 40), nTicks=2, side=4, las=1, cex.axis=0.5)
	addBAMMshifts(best)
	dev.off()
	pdf(lambda_best_plot_file_name, height=7, width=7)
	par(pin=c(7,7))
	best_plot <- plot.bammdata(best, spex="s", lwd=2, pal=c("#268bd2","#6c71c4","#d33682","#dc322f","#cb4b16","#eac115"),labels=T,cex=0.5)
	addBAMMlegend(best_plot, corners = c(10, 12, 30, 40), nTicks=2, side=4, las=1, cex.axis=0.5)
	addBAMMshifts(best)
	dev.off()
	pdf(mu_best_plot_file_name, height=7, width=7)
	par(pin=c(7,7))
	best_plot <- plot.bammdata(best, spex="e", lwd=2, pal=c("#268bd2","#6c71c4","#d33682","#dc322f","#cb4b16","#eac115"),labels=T,cex=0.5)
	addBAMMlegend(best_plot, corners = c(10, 12, 30, 40), nTicks=2, side=4, las=1, cex.axis=0.5)
	addBAMMshifts(best)
	dev.off()

	# Produce plots for the sets of credible shift combinations.
	cset <- credibleShiftSet(edata, prior=priorshifts, BFcriterion=3)
	pdf(set_plot_of_4_file_name, height=7, width=10)
	set_plot <- plot.credibleshiftset(cset, spex="netdiv", plotmax=4, lwd=1.5, pal=c("#268bd2","#6c71c4","#d33682","#dc322f","#cb4b16","#eac115"), legend=T)
	dev.off()
	pdf(set_plot_of_9_file_name, height=7, width=10)
	set_plot <- plot.credibleshiftset(cset, spex="netdiv", plotmax=9, lwd=1.5, pal=c("#268bd2","#6c71c4","#d33682","#dc322f","#cb4b16","#eac115"), legend=T)
	dev.off()

	# Produce a plot with branch lengths according to Bayes Factors.
	bftree <- bayesFactorBranches(edata, priorshifts)
	bftree2 <- ladderize(bftree)
	bftree3 <- rotate(bftree2, 110)
	bftree4 <- bftree3#rotate(bftree3, 80)
	pdf(bayes_factors_plot_file_name, height=7, width=7)
	par(pin=c(7,7))
	plot.phylo(bftree4, cex=0.5, show.node.label=TRUE)
	# nodelabels()
	add.scale.bar(length=10,cex=0.5)
	dev.off()

	# Produce a plot for the macroevolutionary cohort analysis.
	cmat <- getCohortMatrix(edata)
	pdf(cohort_plot_file_name, height=7, width=7)
	cohorts(cmat, edata)
	dev.off()

	# Produce a plot with branch lengths according to mean rate estimates per branch.
	net_div_tree_object <- getMeanBranchLengthTree(edata, rate="ndr")
	net_div_tree <- net_div_tree_object$phy
	pdf(net_div_tree_file_name, height=7, width=7)
	par(pin=c(7,7))
	plot(net_div_tree, cex=0.5)
	dev.off()

	# Write a file to connect branch begin and end with branch rates.
	l0 <- sprintf("")
	l1 <- sprintf("Branch begin")
	l2 <- sprintf(paste(edata$begin, sep=" "))
	l3 <- sprintf("Branch end")
	l4 <- sprintf(paste(edata$end, sep=" "))
	l5 <- sprintf("Branch rate")
	l6 <- sprintf(paste(net_div_tree$edge.length, sep=" "))
	l7 <- sprintf("Branch tip label")
	l8 <- sprintf(paste(net_div_tree$tip.label, sep=" "))
	fileConn<-file(branch_file_name)
	writeLines(c(l1,l2,l0,l3,l4,l0,l5,l6,l0,l7,l8), fileConn)
	close(fileConn)

	# Write a file to list speciation and extinction rates per species.
	speciation_tree_object <- getMeanBranchLengthTree(edata, rate="speciation")
	speciation_tree <- speciation_tree_object$phy
	extinction_tree_object <- getMeanBranchLengthTree(edata, rate="extinction")
	extinction_tree <- extinction_tree_object$phy
	l0 <- sprintf("")
	l1 <- sprintf("Branch begin")
	l2 <- sprintf(paste(edata$begin, sep=" "))
	l3 <- sprintf("Branch end")
	l4 <- sprintf(paste(edata$end, sep=" "))
	l5 <- sprintf("Branch speciation")
	l6 <- sprintf(paste(speciation_tree$edge.length, sep=" "))
	l7 <- sprintf("Branch extinction")
	l8 <- sprintf(paste(extinction_tree$edge.length, sep=" "))
	l9 <- sprintf("Branch tip label")
	l10 <- sprintf(paste(speciation_tree$tip.label, sep=" "))
	fileConn<-file(branch2_file_name)
	writeLines(c(l1,l2,l0,l3,l4,l0,l5,l6,l0,l7,l8,l0,l9,l10), fileConn)
	close(fileConn)
}

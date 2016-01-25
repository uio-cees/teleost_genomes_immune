# Michael Matschiner, 2015-10-06

# Get the command line arguments.
args <- commandArgs(trailingOnly = TRUE)
input_table_file_name <- paste(args[1],"/summary.txt",sep="")
log_likelihood_difference_plot_file_name <- paste(args[1],"/log_likelihood_difference.pdf",sep="")
net_diversification_plot_file_name <- paste(args[1],"/net_diversification.pdf",sep="")

table <- read.table(input_table_file_name,header=T)
rownames(table) <- table$threshold

pdf(log_likelihood_difference_plot_file_name, height=7, width=7)
plot(log(table$threshold),table$median_log_likelihood_difference,type="l",xlim=c(2.302585093,4.094344562),ylim=c(0,70),xlab="MHC I copy number threshold",ylab="log(likelihood) difference",main="Trait-dependent speciation")
lines(log(table$threshold),table$lower_log_likelihood_difference,type="l",col="grey")
lines(log(table$threshold),table$upper_log_likelihood_difference,type="l",col="grey")
dev.off()

pdf(net_diversification_plot_file_name, height=7, width=7)
plot(log(table$threshold),table$mean_unconstrained_netdiv_low_copy_number,type="l",xlim=c(2.302585093,4.094344562),ylim=c(0.05,0.08),xlab="MHC I copy number threshold",ylab="Net diversification",main="Trait-dependent speciation")
lines(log(table$threshold),table$mean_unconstrained_netdiv_low_copy_number-table$stdev_unconstrained_netdiv_low_copy_number,type="l",col="grey")
lines(log(table$threshold),table$mean_unconstrained_netdiv_low_copy_number+table$stdev_unconstrained_netdiv_low_copy_number,type="l",col="grey")
lines(log(table$threshold),table$mean_unconstrained_netdiv_high_copy_number,type="l")
lines(log(table$threshold),table$mean_unconstrained_netdiv_high_copy_number-table$stdev_unconstrained_netdiv_high_copy_number,type="l",col="grey")
lines(log(table$threshold),table$mean_unconstrained_netdiv_high_copy_number+table$stdev_unconstrained_netdiv_high_copy_number,type="l",col="grey")
dev.off()

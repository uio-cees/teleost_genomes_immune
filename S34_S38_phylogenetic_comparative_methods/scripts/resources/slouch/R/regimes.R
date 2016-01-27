`regimes` <-
function (topology, times, regime.specs, term) {
		N <- length(term);
		reg <- set.of.regimes(topology,regime.specs);
		R <- length(reg);
		beta <- vector(R*N, mode="list");
		for (i in 1:N) {
			for (k in 1:R) {
				p <- pedigree(topology, term[i]);
				n <- length(p);
				beta[[i + N*(k-1)]] <- as.integer(regime.specs[p[1:(n-1)]] == reg[k]);
			}
		}    
		return(beta);
	}


`epochs` <-
function (topology, times, term) {
		N <- length(term);
		e <- vector(length=N,mode="list");
		for (k in 1:N) {
			p <- pedigree(topology,term[k]);
			e[[k]] <- times[p];	
		}
		return(e);
	}


`set.of.regimes` <-
function (topology, regime.specs) {
		n <- length(regime.specs);
		id <- seq(1,n)[topology > 0];       # find all non-root nodes
		reg <- sort(unique(regime.specs[id]));
		return(reg);		  
	}


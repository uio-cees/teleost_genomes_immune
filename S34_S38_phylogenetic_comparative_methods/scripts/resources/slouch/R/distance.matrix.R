`distance.matrix` <-
function (topology, times) {
  term <- terminal.twigs(topology);
  N <- length(term);
  dm <- matrix(data=0,nrow=N,ncol=N);
  dm[1,1] <- 0;
  for (i in 2:N) {
    pedi <- pedigree(topology,term[i]);
    for (j in 1:(i-1)) {
      pedj <- pedigree(topology,term[j]);
      for (k in 1:length(pedi)) {
        if (any(pedj == pedi[k])) break;
      }
      dm[j,i] <- dm[i,j] <- (times[term[i]]-times[pedi[k]]) + (times[term[j]]-times[pedi[k]]); 
    }
    dm[i,i] <- 0;
  }
  return(dm);
}


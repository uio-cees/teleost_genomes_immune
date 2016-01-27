`tsja` <-
function (topology, times) {
  term <- terminal.twigs(topology);
  N <- length(term);
  t.ja <- matrix(data=0,nrow=N,ncol=N);
  t.ja[1,1]=0
  for (i in 2:N) {
    pedi <- pedigree(topology,term[i]);
    for (j in 1:(i-1)) {
      pedj <- pedigree(topology,term[j]);
      for (k in 1:length(pedi)) {
        if (any(pedj == pedi[k])) break;
      }
      t.ja[j,i] <- t.ja[i,j] <- (times[term[j]]-times[pedi[k]]); 
    }
    t.ja[i,i] <- 0;
  }
  return(t.ja);
}


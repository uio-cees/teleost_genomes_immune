`branch.times` <-
function (topology, times) {
  term <- terminal.twigs(topology);
  N <- length(term);
  bt <- matrix(data=0,nrow=N,ncol=N);
  bt[1,1] <- times[term[1]];
  for (i in 2:N) {
    pedi <- pedigree(topology,term[i]);
    for (j in 1:(i-1)) {
      pedj <- pedigree(topology,term[j]);
      for (k in 1:length(pedi)) {
        if (any(pedj == pedi[k])) break;
      }
      bt[j,i] <- bt[i,j] <- times[pedi[k]];
    }
    bt[i,i] <- times[term[i]];
  }
  return(bt);
}


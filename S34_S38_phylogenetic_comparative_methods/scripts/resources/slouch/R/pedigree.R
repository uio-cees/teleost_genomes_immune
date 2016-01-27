`pedigree` <-
function (topology, k) {
  p <- k;
  k <- topology[k];
  while (k != 0) {
    p <- c(p, k);
    k <- topology[k];
  }
  return(p);
}


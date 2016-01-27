`arrange.tree` <-
function (root, topology) {
  k <- which(topology==root);
  n <- length(k);
  reltree <- rep(0,length(topology));
  reltree[root] <- 0.5;
  p <- list(NULL);
  if (n > 0) {
    m <- rep(0,n);
    for (j in 1:n) {
      p[[j]] <- arrange.tree(k[j],topology);
      m[j] <- length(which(p[[j]] != 0));
    }
    cm <- c(0,cumsum(m));
    for (j in 1:n) {
      reltree <- reltree + (cm[j]/sum(m))*(p[[j]] != 0) + (m[j]/sum(m))*p[[j]];
    }
  }
  return(reltree);
}


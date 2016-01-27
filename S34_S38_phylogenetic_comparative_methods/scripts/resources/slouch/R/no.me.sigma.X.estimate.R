no.me.sigma.X.estimate <- function (predictor, topology, times) {
  pt <- parse.tree(topology,times);
  n <- pt$N;
  v <- pt$bt;
  w <- matrix(data=1,nrow=pt$N,ncol=1);
  dat <- predictor[!is.na(predictor)];
  beta<-solve(t(w)%*%solve(v)%*%w)%*%(t(w)%*%solve(v)%*%dat)
  e<-dat-beta
  theta <- beta
  sigma <- sqrt((e %*% solve(v,e))/(n-1));
  dim(sigma) <- 1;
  return(list(as.numeric(theta), as.numeric(sigma)));
}

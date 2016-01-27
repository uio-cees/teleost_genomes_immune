sigma.X.estimate <- function (predictor,me.predictor, topology, times) {
  pt <- parse.tree(topology,times);
  n <- pt$N;
  v <- pt$bt;
  w <- matrix(data=1,nrow=pt$N,ncol=1);
  me<-diag(me.predictor[!is.na(me.predictor)])
  dat <- predictor[!is.na(predictor)];
  beta<-solve(t(w)%*%solve(v)%*%w)%*%(t(w)%*%solve(v)%*%dat)
  e<-dat-beta
  sigma<-as.numeric((t(e)%*%solve(v)%*%e)/(n-1))
  repeat{
   beta<-solve(t(w)%*%solve(v + me/sigma)%*%w)%*%(t(w)%*%solve(v + me/sigma)%*%dat)
   e<-dat-beta
   theta <- beta
   sigma1<-(t(e)%*%solve(v +me/sigma)%*%e)/(n-1)
   if(abs(as.numeric(sigma1)-sigma) <= 0.0000001*sigma) break
   sigma <- as.numeric(sigma1)
  }
  return(list(as.numeric(theta), as.numeric(sigma)));
}

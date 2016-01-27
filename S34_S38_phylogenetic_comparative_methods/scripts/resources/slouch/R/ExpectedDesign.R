`ExpectedDesign` <-
function(topology, times, alpha, random.cov, me.random.cov)
{
 #### Need to work in code for multiple predictors as well as for when we have a fixed factor i.e. generalize
  T <- times[terminal.twigs(topology)];	
  n<-length(T) 
  pred<-data.frame(random.cov)
  n.pred<-length(pred[1,]);
  pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred);
  me.pred<-matrix(data=me.random.cov[!is.na(me.random.cov)], ncol=n.pred);	
	
  D<-matrix(0,nrow=n,ncol=n.pred+1) # generalize for multiple predictors
  currDcol<-1
	s.T<-matrix(data=0, ncol=n.pred)  # PREDICTOR THETA
	
	for(i in 1:n.pred)
	{
		s.T[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[1]);
	}
		
  
  D[,1]<-1
	for(i in 2:(n.pred+1))
	{	
     D[,i]<-(1-(1-exp(-alpha*T))/(alpha*T))*s.T[,(i-1)]
# change for sigma.X.estimate(pred, me.pred, topology, times)[[1]]	
	}	
  return(D)	
}


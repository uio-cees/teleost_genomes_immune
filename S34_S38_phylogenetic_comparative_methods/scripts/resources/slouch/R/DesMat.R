`DesMat` <-
function(topology, times, random.cov, me.random.cov=NULL, alpha)
{
T <- times[terminal.twigs(topology)];
	pred<-data.frame(random.cov)
n.pred<-length(pred[1,]);
pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred);
	
n<-length(T)	 
  D<-matrix(0,nrow=n,ncol=n.pred+1)
  D[,1]<-1
	D[,2:(n.pred+1)]<-(1-(1-exp(-alpha*T))/(alpha*T))*pred  
  return(D)			
}


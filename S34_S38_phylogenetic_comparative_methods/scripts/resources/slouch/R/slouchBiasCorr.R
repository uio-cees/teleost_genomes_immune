`slouchbiascorr` <-
function(beta, topology, times, half_life_value, vy, response, me.response=NULL, random.cov, me.random.cov=NULL)
{
Y <- response[!is.na(response)]; 
N <- length(Y);	
	alpha<-log(2)/half_life_value;	
pred<-data.frame(random.cov)
n.pred<-length(pred[1,]);	

SlouchD<-DesMat(topology, times, random.cov, me.random.cov,alpha)	
		
DesPhyCov<-function(topology, times, DesignMatrix, random.cov, me.random.cov, alpha)
  {
		pred<-data.frame(random.cov)
		n.pred<-length(pred[1,]);
		pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred);	 
	  me.pred<-matrix(data=me.random.cov[!is.na(me.random.cov)], ncol=n.pred);
		kX<-n.pred	
	  
	  s.X<-matrix(data=0, ncol=n.pred)  # PREDICTOR SIGMA
	
	  for(i in 1:n.pred)
	  {
		  s.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[2]);
	  }
	 
		bt<-parse.tree(topology, times)$bt
		T <- times[terminal.twigs(topology)];
#s.X<-as.numeric(sigma.X.estimate(random.cov, me.random.cov, topology, times)[[2]])
		Vd<-matrix(0,ncol(DesignMatrix)*nrow(DesignMatrix),ncol(DesignMatrix)*nrow(DesignMatrix))
	  BMCovMat<-matrix(0, nrow=nrow(bt), nco=ncol(bt))
	  for(i in 1:n.pred)
	  {
		  BMCovMat<-BMCovMat+(bt*sqrt(s.X[,i]))  ### CHECK IF SIMPLE ADDITION OF VARIANCES IS OK BUG ???
	  }
	   
		Zeta<-(((1-(1-exp(-alpha*T))/(alpha*T)))*bt)^2
	  
	  
		Vd[((ncol(DesignMatrix)-kX)*nrow(DesignMatrix)+1):(ncol(DesignMatrix)*nrow(DesignMatrix)),((ncol(DesignMatrix)-kX)*nrow(DesignMatrix)+1):(ncol(DesignMatrix)*nrow(DesignMatrix))]<-Zeta*BMCovMat	
		return(Vd)	
	}
	
 
SlouchVd<-DesPhyCov(topology, times, SlouchD, random.cov, me.random.cov, alpha=alpha)
	 

DesMeCov<-function(topology, times, alpha, DesignMatrix, me.random.cov)  # equivalent to calculateDesignCovariance.dim1x1.simpReg
 {
  me.pred<-data.frame(me.random.cov)
  n.pred<-length(me.pred[1,]);	 
  kX<-n.pred	
  bt<-parse.tree(topology, times)$bt
  T <- times[terminal.twigs(topology)];
  Vu<-matrix(0,ncol(DesignMatrix)*nrow(DesignMatrix),ncol(DesignMatrix)*nrow(DesignMatrix))
  
  Zeta<-(((1-(1-exp(-alpha*T))/(alpha*T)))*bt)^2
	 MeCovMat<-matrix(0, nrow=nrow(bt), ncol=ncol(bt)) 	
	 for(i in 1:n.pred)
	 {
		 MeCovMat<-MeCovMat+Zeta*diag(me.pred[,i][!is.na(me.pred[,i])])  ### CHECK IF SIMPLE ADDITION OF VARIANCES IS OK???
	 }	 
	 
  Vu[((ncol(DesignMatrix)-kX)*nrow(DesignMatrix)+1):(ncol(DesignMatrix)*nrow(DesignMatrix)),((ncol(DesignMatrix)-kX)*nrow(DesignMatrix)+1):(ncol(DesignMatrix)*nrow(DesignMatrix))]<-MeCovMat
  return(Vu)
 } 
 
SlouchVu<-DesMeCov(topology, times, alpha, SlouchD, me.random.cov)  
SlouchED<-ExpectedDesign(topology, times, alpha, random.cov, me.random.cov) 
V<-ResidV(beta, topology, times, alpha, vy, response, me.response, random.cov, me.random.cov, mecov.random.cov=NULL, ultrametric=TRUE)
	
### OUTPUT ###
	

SlouchGLSestims<-pseudoinverse(t(SlouchD)%*%solve(V)%*%SlouchD)%*%t(SlouchD)%*%solve(V)%*%response[!is.na(response)]
SlouchGLSestimsBiasCorr<-BiasCor(SlouchGLSestims,SlouchD,SlouchED,V,SlouchVu,SlouchVd)
bK<- SlouchGLSestims / SlouchGLSestimsBiasCorr 
SlouchGLSestimsBiasVar<-pseudoinverse(diag(bK^2)*(t(SlouchD)%*%solve(V)%*%SlouchD))	### NEED TO CHECK THIS??
output<-matrix(data=0, nrow=(n.pred+1), ncol=4, dimnames=list(c("Intercept", if(n.pred==1) deparse(substitute(random.cov)) else colnames(random.cov)), c("B(biased)", "Bias factor", "B(unbiased)", "Std. Error")))

	output[,1]<-round(SlouchGLSestims,4)
	output[,2]<- round(bK,4)
	output[,3]<-round(SlouchGLSestimsBiasCorr,4)
	output[,4]<-round(sqrt(diag(SlouchGLSestimsBiasVar)),4)
	
	print(output)	
}


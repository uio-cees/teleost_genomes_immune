`DesPhyCov` <-
function(topology, times, DesignMatrix, random.cov, me.random.cov)
 {
	 pred<-data.frame(random.cov)
	 n.pred<-length(pred[1,]);
	 pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred);	 
 kX<-n.pred	
 bt<-parse.tree(topology, times)$bt
  T <- times[terminal.twigs(topology)];
  s.X<-as.numeric(sigma.X.estimate(random.cov, me.random.cov, topology, times)[[2]])
  Vd<-matrix(0,ncol(DesignMatrix)*nrow(DesignMatrix),ncol(DesignMatrix)*nrow(DesignMatrix))
  BMCovMat<-bt*sqrt(s.X)
  Zeta<-(((1-(1-exp(-alpha*T))/(alpha*T)))*bt)^2
  Vd[((ncol(DesignMatrix)-kX)*nrow(DesignMatrix)+1):(ncol(DesignMatrix)*nrow(DesignMatrix)),((ncol(DesignMatrix)-kX)*nrow(DesignMatrix)+1):(ncol(DesignMatrix)*nrow(DesignMatrix))]<-Zeta*BMCovMat	
  return(Vd)	
 }


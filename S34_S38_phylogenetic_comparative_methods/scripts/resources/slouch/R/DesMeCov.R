`DesMeCov` <-
function(topology, times, alpha, DesignMatrix, me.pred)  # equivalent to calculateDesignCovariance.dim1x1.simpReg
 {
 kX<-1	
  bt<-parse.tree(topology, times)$bt
  T <- times[terminal.twigs(topology)];
  Vu<-matrix(0,ncol(DesignMatrix)*nrow(DesignMatrix),ncol(DesignMatrix)*nrow(DesignMatrix))
  
  Zeta<-(((1-(1-exp(-alpha*T))/(alpha*T)))*bt)^2
  Vu[((ncol(DesignMatrix)-kX)*nrow(DesignMatrix)+1):(ncol(DesignMatrix)*nrow(DesignMatrix)),((ncol(DesignMatrix)-kX)*nrow(DesignMatrix)+1):(ncol(DesignMatrix)*nrow(DesignMatrix))]<-if(is.matrix(me.pred)) Zeta*me.pred else Zeta*diag(me.pred)
  return(Vu)
 }


`slouchtree.plot` <-
function (topology, times, names = NULL, regimes = NULL, cex = NULL, lwd=NULL, reg.col=NULL) {
  if(is.null(cex)) cex<-1;	
  if(is.null(lwd)) lwd<-1;	
    rx <- range(times);
  rxd <- 0.1*diff(rx);

  if (is.null(regimes))
    regimes <- factor(rep(1,length(topology)));

  levs <- levels(as.factor(regimes));
  palette <- rainbow(length(levs));

  for (r in 1:length(levs)) {
    y <- tree.layout(topology);
    x <- times;
    f <- which(topology > 0 & regimes == levs[r]);
    pp <- topology[f];
    X <- array(data=c(x[f], x[pp], rep(NA,length(f))),dim=c(length(f),3));
    Y <- array(data=c(y[f], y[pp], rep(NA,length(f))),dim=c(length(f),3));
    oz <- array(data=1,dim=c(2,1));
    X <- kronecker(t(X),oz);
    Y <- kronecker(t(Y),oz);
    X <- X[2:length(X)];
    Y <- Y[1:(length(Y)-1)];
    if(!is.null(regimes)) 
    {if(is.null(reg.col))	
    C <- rep(palette[r],length(X)) 
    }   
    {if(!is.null(reg.col)) 	
    C <- rep(reg.col[r],length(X)) 
    }    	   
     if (r > 1) par(new=TRUE);
    par(yaxt='n')
    par(bty="n")
    par(font="2")
    
    plot(X,Y,type='l',col=C,lwd=lwd,xlab='time',ylab='',xlim = rx + c(-rxd,rxd),ylim=c(0,1));
    if (!is.null(names))
      text(X[seq(1,length(X),6)],Y[seq(1,length(Y),6)],names[f],pos=4, cex=cex);
  }
  par(yaxt="s") #reset graphic parameter to default
  par(bty="o")
  par(font="1")
}


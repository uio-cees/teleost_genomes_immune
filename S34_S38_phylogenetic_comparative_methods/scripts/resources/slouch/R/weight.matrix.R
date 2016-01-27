`weight.matrix` <-
function (alpha, topology, times, N, regime.specs, fixed.cov, intercept) 
{
		
    N <- N
    reg <- set.of.regimes(topology, regime.specs)
    R <- length(reg)
    T <- times[terminal.twigs(topology)]
    ep <- epochs(topology, times, terminal.twigs(topology))
    beta <- regimes(topology, times, regime.specs, terminal.twigs(topology))
    W <- matrix(data = 0, nrow = N, ncol = R + 1, dimnames = list(c(), 
        c("Ya", as.character(set.of.regimes(topology, regime.specs)))))
    W[, 1] <- exp(-alpha * T)
    for (i in 1:N) {
        delta <- diff(exp(alpha * (ep[[i]] - T[i])))
        for (k in 1:R) {
            W[i, k + 1] <- -sum(delta * beta[[i + N * (k - 1)]])
        }
    }
    if (is.null(intercept)) 
        W <- W
    else {
        if (intercept == "root") {
            root.reg <- as.character(regime.specs[times == 0])
            nonroot.reg <- as.character(reg[reg != root.reg])
            int <- as.matrix(W[, 1] + W[, root.reg])
            colnames(int) = root.reg
            if (max(int[, 1]) <= 0.01) 
                W <- W
            else {
                W2 <- cbind(int, W[, nonroot.reg])
                W <- W2
            }
        }
        else W[, 1] <- intercept
    }
    if (max(W[, 1]) <= 0.01)  
        W <- W[, -1];
        
        if(!is.null(fixed.cov))
{	
 fixed.pred<-data.frame(fixed.cov);
  n.fixed.pred<-length(fixed.pred[1,]);
  fixed.pred<-matrix(data=fixed.pred[!is.na(fixed.pred)], ncol=n.fixed.pred);  
  W<-cbind(W, fixed.pred)
 }	
       
        return(W)
}


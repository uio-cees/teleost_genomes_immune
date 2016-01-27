`fast.svd` <-
function (m, tol) 
{
    n <- dim(m)[1]
    p <- dim(m)[2]
    EDGE.RATIO <- 2
    if (n > EDGE.RATIO * p) {
        return(psmall.svd(m, tol))
    }
    else if (EDGE.RATIO * n < p) {
        return(nsmall.svd(m, tol))
    }
    else {
        return(positive.svd(m, tol))
    }
}


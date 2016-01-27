`psmall.svd` <-
function (m, tol) 
{
    B <- crossprod(m)
    s <- svd(B, nu = 0)
    if (missing(tol)) 
	tol <- dim(B)[1] * max(s$d) * .Machine$double.eps
    Positive <- s$d > tol
    d <- sqrt(s$d[Positive])
    v <- s$v[, Positive, drop = FALSE]
    u <- m %*% v %*% diag(1/d, nrow = length(d))
    return(list(d = d, u = u, v = v))
}


`nsmall.svd` <-
function (m, tol) 
{
    B <- m %*% t(m)
    s <- svd(B, nv = 0)
    if (missing(tol)) 
	tol <- dim(B)[1] * max(s$d) * .Machine$double.eps
    Positive <- s$d > tol
    d <- sqrt(s$d[Positive])
    u <- s$u[, Positive, drop = FALSE]
    v <- crossprod(m, u) %*% diag(1/d, nrow = length(d))
    return(list(d = d, u = u, v = v))
}


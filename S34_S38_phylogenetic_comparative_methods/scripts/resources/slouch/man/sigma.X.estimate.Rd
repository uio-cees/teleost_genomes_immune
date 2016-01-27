\name{sigma.X.estimate}
\alias{sigma.X.estimate}
\title{Estimate Brownian Motion evolution model parameters}
\description{
  \code{sigma.X.estimate} can be used to estimate \eqn{\theta} (BM mean or estimated state at root of phylogeny) and \eqn{\sigma} squared (BM variance) for a univariate Brownian motion of evolution for a trait on a given phylogeny. The \code{\link{model.fit}} function uses it to estimate the Brownian Motion model parameters for the randomly evolving predictor variable /s.   
}
\usage{
sigma.X.estimate(predictor, me.predictor, topology, times)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predictor}{A predictor variable vector for the OUBM model) in the slouch tree data frame }
  \item{me.predictor}{A vector of the observation variance for the predictor variable}
  \item{topology}{A column in the \pkg{slouch} tree data frame named \code{ancestor} that specifies which node a current species or internal node is connected to. The ancestral node for the node at the root of the tree gets a numeric value of 0.  }
  \item{times}{A column in the slouch tree data frame named "time" specifying the cumulative time of the current node in relation to the root node, which is at time 0.}
}
\details{
  \code{sigma.X.estimate} returns restricted maximum likelihood estimates (REML), so may differ slightly from the e.g. \code{brown} function estimates in \pkg{ouch}, which uses maximum likelihood. The first number returned is theta and the second sigma squared.
}

\author{Jason Pienaar \email{jasonpienaar@gmail.com} }
\seealso{\code{\link{model.fit}}}
\examples{
data(Lio)
attach(Lio)
## e.g. estimate brownian motion mean and variance parameters for body temperature (Tb)
## no observation error variance (requires a place holder vector) 
me<-rep(0, times=30)
sigma.X.estimate(Tb, me, ancestor, time)
## include observation error variance for body temp (mv.tb)
sigma.X.estimate(Tb,mv.tb, ancestor, time)
detach(Lio)
}
\keyword{models}



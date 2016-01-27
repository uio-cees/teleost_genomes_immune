\name{make.tree}
\alias{make.tree}

\title{Create a bifurcating tree of n tip species.}
\description{
 A rudimentary function for creating symmetric,  bifurcating trees of n tips (allowed values are given below) in slouch tree format. It is intended for use with simulations to examine the influence of tree shape on estimation precision, bias and accuracy. Species names are randomly picked as different colors. The tree height is 1 and internal nodes are given as proportions of this. Note that \pkg{ape} has more sophisticated tree generating functions so an alternative is to create such trees in \pkg{ape} then convert them to \pkg{slouch} trees via \code{ape2ouch} and \code{ouch2slouch}. 
}
\usage{
make.tree(n.tips, shape = NULL, stretch = NULL, ouch = FALSE)
}
\arguments{
  \item{n.tips}{the number of tip species to create a tree for (must be larger than 3}
  \item{shape}{ the default \code{NULL} value will return a fully balanced tree, if \code{shape} is set to anything else, \code{make.tree} will return a fully unbalanced tree}
  \item{stretch}{A numeric value  to which the time vector can be raised. For example, for a value of 2, all times are squared, and the tip branches become longer relative to branches closer to the root. A value of 0.5 would give the opposite effect.}
  \item{ouch}{The default \code{FALSE} value will return a slouch tree data frame, if \code{ouch} is set to true, \code{make.tree} will return an \code{ouchtree} object for use with the \pkg{ouch} package}
}
\author{Jason Pienaar \email{jasonpienaar@gmail.com} }
\seealso{
\code{\link{slouchtree.plot}}, \code{\link{slouch2ouch}}, \code{\link{ouch2slouch}}  
}
\examples{
## make a slouch tree with 8 tips
eg.tree<-make.tree(8)
slouchtree.plot(eg.tree$ancestor, eg.tree$time)
## make non-symmetric tree with 8 tips
eg.tree<-make.tree(8, shape=1)
slouchtree.plot(eg.tree$ancestor, eg.tree$time)
## stretch branches at tips 
eg.tree<-make.tree(8, stretch=2)
slouchtree.plot(eg.tree$ancestor, eg.tree$time)
## stretch branches closer to root
eg.tree<-make.tree(8, stretch=0.25)
slouchtree.plot(eg.tree$ancestor, eg.tree$time)
## make an ouchtree object (requires package ouch to work)
\dontrun{eg.ouch.tree<-make.tree(8, ouch=TRUE)
class(eg.ouch.tree)}
}
\keyword{datagen}

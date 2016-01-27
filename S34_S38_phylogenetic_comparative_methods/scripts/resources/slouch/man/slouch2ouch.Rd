\name{slouch2ouch}
\alias{slouch2ouch}
\title{Convert a slouch tree data frame to an ouch tree object}
\description{
  \code{slouch2ouch} translates \pkg{slouch}'s data frame representation of a phylogenetic tree into \pkg{ouch}'s  \code{ouchtree} object. The \pkg{ouch} package includes a function \code{ape2ouch} that converts \pkg{ape} tree objects of class \code{phylo} into \pkg{ouch} tree objects, and the former includes a number of useful functions for tree pruning and plotting.   
}
\usage{
slouch2ouch(nodes, ancestors, times, labels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodes}{a column in the slouch tree data frame representing the node numbers}
  \item{ancestors}{a column in the slouch tree data frame representing the ancestral nodes for a given node}
  \item{times}{a column in the slouch data frame representing the times at which nodes occur}
  \item{labels}{the \code{species} names column in the slouch data frame }
}
\value{
  returns an \code{ouchtree} object for use with the \pkg{ouch} package}
\author{Jason Pienaar \email{jasonpienaar@gmail.com}}
\seealso{\code{\link{ouch2slouch}}, \code{\link{make.tree}}, \code{\link{slouchtree.plot}}}
\examples{
## convert a slouch tree data.frame into an ouchtree object
## requires ouch package so not run through R CMD check
\dontrun{eg.ouch.tree<-make.tree(8, ouch=TRUE)
class(eg.ouch.tree)
slouch.tree<-ouch2slouch(eg.ouch.tree)
class(slouch.tree)
ouch.tree<-slouch2ouch(slouch.tree$nodes, slouch.tree$ancestor, slouch.tree$time, slouch.tree$species)
class(ouch.tree)}
}
\keyword{manip}


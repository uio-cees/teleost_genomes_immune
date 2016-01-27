\name{ouch2slouch}
\alias{ouch2slouch}
\title{Convert an ouch tree object to a slouch tree data frame}
\description{
 \code{ouch2slouch} translates \pkg{ouch}'s \code{ouchtree} object representation of a phylogenetic tree into \pkg{slouch}'s data frame representation. The \pkg{ouch} package includes a function \code{ape2ouch} that converts \pkg{ape} tree objects of class \code{phylo} into \code{ouchtree} objects,  and the former includes a number of useful functions for tree pruning and plotting. The principle use of this function is to allow users to obtain a \pkg{slouch} tree from an original nexus or newick representation, possibly pruned via \pkg{ape} and \pkg{ouch}.   

}
\usage{
ouch2slouch(tree)
}
\arguments{
  \item{tree}{An object of class \code{ouchtree} created in the \pkg{ouch} package }
}

\value{
 Returns a \pkg{slouch} data frame representation of the tree }
\author{Jason Pienaar \email{jasonpienaar@gmail.com}}
\seealso{\code{\link{slouch2ouch}},\code{\link{make.tree}},\code{\link{slouchtree.plot}}}
\examples{
## convert an ouchtree object to a slouch tree data frame
## make an ouchtree object (requires ouch package so not run through R CMD check)
\dontrun{eg.ouch.tree<-make.tree(8, ouch=TRUE)
class(eg.ouch.tree)
slouch.tree<-ouch2slouch(eg.ouch.tree)
class(slouch.tree)
slouchtree.plot(slouch.tree$ancestor, slouch.tree$time)}
}

\keyword{manip}


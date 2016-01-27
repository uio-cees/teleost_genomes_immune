% --- Source file: slouchtree.plot.R ---
\name{slouchtree.plot}
\alias{slouchtree.plot}
\title{Plot a trees and niches encoded in slouch format as a phylogram}
\description{
Allows for plotting the tree and visualizing reconstructed ancestral states on the internal branches. 
}
\usage{
slouchtree.plot(topology, times, names = NULL, regimes = NULL, cex = NULL, lwd = NULL, reg.col = NULL)
}
\arguments{
  \item{topology}{A column in the \pkg{slouch} tree data.frame named ancestor that specifies which node a current species or internal node is connected to. The ancestral node for the node at the root of the tree gets a numeric value of 0.}
  \item{times}{A column in the slouch tree data frame specifying the cumulative time of the current node in relation to the root node, which occurs at time 0.}
  \item{names}{An optional vector of species names, internal nodes are encoded as NA}
  \item{regimes}{An optional vector of codes (numeric or character) that specify different categorical predictor variables (regimes). Usually these would also include a reconstruction of the states on the phylogeny so that the internal nodes are also coded according to which regime the branch leading back from the current node belongs to (using e.g. the \code{fitch} function).  
}
  \item{cex}{As in the R base package, \code{cex} controls text size.}
  \item{lwd}{As in the R base package, \code{lwd} controls line width.}
  \item{reg.col}{An optional vector of color names (as characters) that allows user-specified coloring of the different regimes.}
}
\details{
  \code{plot.tree} was borrowed and modified from the \pkg{ouch} package written by Marguerite Butler and Aaron King (2004), and is therefore similar the plotting functions in that package. 
}
\references{ 
Butler, M. A. & King. A. A. (2004) Phylogenetic comparative analysis: a modeling approach for adaptive evolution. \emph{The american Naturalist} \bold{164}: 683--695.  
}
\author{Jason Pienaar \email{jasonpienaar@gmail.com}}
\seealso{ \code{\link{fitch}},\code{\link{make.tree}},\code{\link{ouch2slouch}},\code{\link{slouch2ouch}}}
\examples{
data(Lio)
attach(Lio)
## plot tree topology
slouchtree.plot(ancestor, time)
## add species names, adjust text size and line width
slouchtree.plot(ancestor, time, species, cex=0.8, lwd=2)
## plot niche reconstructions
pars.climate<-fitch(Lio, Lio$climate, root="Tropical")
slouchtree.plot(ancestor, time, species, cex=0.8, lwd=2, regimes=pars.climate)
pars.del.temp<-fitch(Lio, Lio$temp, root="Warm", deltran=TRUE)
slouchtree.plot(ancestor, time, species, cex=0.8, lwd=2, regimes=pars.del.temp)
## choose own regime colors
slouchtree.plot(ancestor, time, species, cex=0.8, lwd=2, regimes=pars.del.temp, reg.col=c("blue", "red"))
detach(Lio)
}
\keyword{aplot}

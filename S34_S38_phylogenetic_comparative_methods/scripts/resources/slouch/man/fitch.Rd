% --- Source file: fitch.R ---
\name{fitch}
\alias{fitch}
\title{Unordered Fitch parsimony reconstruction of discrete character states}
\description{
Implements an unordered Fitch parsimony reconstruction of discrete niche     variables for use in the OU models where optima are modeled on discrete, categorical niche encodings. Allows for delayed and accelerated transformations to deal with ambiguities. 
}
\usage{
fitch(tree.data, niche, deltran = FALSE, acctran = FALSE, root = NULL)
}
\arguments{
  \item{tree.data}{The attached \pkg{slouch} data.frame (i.e. use \code{attach}), 	including all tree and taxon data}
  \item{niche}{The specific niche variable in the \pkg{slouch} data.frame to be 	reconstructed, entered as data.frame\code{$}niche }
  \item{deltran}{Implements a delayed transformation algorithm in order to 	deal with ambiguous nodes}
  \item{acctran}{Implements an accelerated transformation algorithm to deal 	with ambiguous nodes}
  \item{root}{An optional argument allowing the user to define a character 	state for the root (useful if the root node is ambiguously 		reconstructed)}
}
\details{
The \code{fitch} function is meant to be interactive, where the user acts on the advice given in the returned messages whilst attempting to reconstruct ancestral states. If the root node is ambiguous after an initial reconstruction (a message will be printed to the screen if this is the case), this needs to be set by the user using the \code{root = "state"} argument in the function call. Any remaining ambiguous nodes can then be dealt with by specifying \code{deltran} or \code{acctran ="TRUE"} in the function call 
}
\value{
The \code{fitch} function returns a vector of reconstructed character states, which can be used directly in the \code{\link{slouchtree.plot}} or \code{\link{model.fit}} functions of \pkg{slouch} if assigned to a variable name
}
\references{ 
Fitch, M.W. (1971) Defining the course of Evolution: Minimum change for a specific tree topology. \emph{Systematic Zoology} \bold{20}:406--416.

Swofford, D. L. and W.P. Maddison (1987) Reconstructing ancestral character states under Wagner parsimony. \emph{Mathematical Biosciences} \bold{87}: 199--229. 
}
\author{Jason Pienaar \email{jasonpienaar@gmail.com}}
\seealso{
\code{\link{slouchtree.plot}}, \code{\link{model.fit}}
}
\examples{
data(Lio)
attach(Lio)
## reconstruct climate niches
pars.climate<-fitch(Lio, Lio$climate)
## ambiguity at root, set to e.g. root="Tropical" as suggested by function
pars.climate<-fitch(Lio, Lio$climate, root="Tropical")
pars.climate
## reconstruct temperate niches
pars.temp<-fitch(Lio, Lio$temp)
## ambiguity at root and several internal nodes, set to e.g. root="Warm"
pars.temp<-fitch(Lio, Lio$temp, root="Warm")
## still 6 internal ambiguities, attempt deltran reconstruction
pars.del.temp<-fitch(Lio, Lio$temp, root="Warm", deltran=TRUE)
## no more ambiguities, note that actran reconstruction does not work on this example
detach(Lio)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

\name{slouch-package}
\alias{slouch-package}
\alias{slouch}
\docType{package}
\title{
Stochastic Linear Ornstein Uhlenbeck models for Comparative Hypotheses
}
\description{
the \pkg{slouch} package provides functions for fitting macro evolutionary 	models of trait evolution on a phylogeny as well as useful tree 	plotting, parsimony reconstruction, iterated GLS regression 		functions for incorporating observation error of continuous 		variables, and bias correction functions for when predictor 		variable error biases regression slopes. Macro evolutionary trait 	optima can be modeled as a single global optimum for all species, 	as a function of different fixed categorically encoded niches that 	can change over time for groups of species (similar to the 		\pkg{ouch} package, as a multiple or single predictor regression on 	fixed, continuously encoded measures of other traits or the 		external environment,  as a multiple or single predictor regression 	on randomly evolving continuous measures of the external 		environment, or as any combination of the preceding models. 		Observation error for any continuously measured variable, as well 	as more complex interactions between observation errors, such as 	covariances between response and predictor observation errors can 	be included in the model fitting functions. A would-be user should 	also consult the \pkg{ouch} page at: 
	http://tsuga.biology.lsa.umich.edu/ouch/.}
\details{
\tabular{ll}{
Package: \tab slouch\cr
Type: \tab Package\cr
Version: \tab 1.2-0\cr
Date: \tab 2011-12-28\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
The starting point for \pkg{slouch} analysis is a \R \code{data.frame} containing information for encoding a phylogenetic tree, vectors of trait values for fitting models to, vectors of categorical or continuous predictor variables to model optima on and optionally vectors of observation errors for all continuous variables. The \pkg{slouch} tree encoding can be obtained from newick / nexus files via the \pkg{ape} (\code{read.nexus}) and \pkg{ouch} (\code{ape2ouch}) packages using \code{\link{ouch2slouch}}, or can be encoded manually in a \code{.txt} or \code{.csv} file and read into \R using the standard \code{\link{read.table}} function of the base package. Columns for the other variables can then either be read in along with the manually encoded file, or added to the \code{\link{data.frame}} using e.g. the \R base package \code{\link{cbind}} command. Fixed, categorical niches can be mapped onto a phylogeny using an unordered \code{\link{fitch}} parsimony algorithm and viewed using \code{\link{slouchtree.plot}}. Macroevolutionary trait optima can be estimated by modeling them on any combination of fixed and randomly evolving predictor variables using the \code{\link{model.fit}} function. Biases due to predictor measurement error can be corrected for using the \code{\link{slouchbiascorr}} function, and non-phylogentically controlled regressions that account for observation errors can be implemented in \code{\link{np.regression}}.   
  
}
\author{
Jason Pienaar, Thomas F Hansen, Steven H Orzack and Kzrysztof Bartoszek 

Maintainer: Jason Pienaar <jasonpienaar@gmail.com>

}
\references{
Butler, M. A. & King. A. A. (2004) Phylogenetic comparative analysis: a modeling approach for adaptive evolution. \emph{The american Naturalist} \bold{164}: 683--695.  

Hansen, T. F. (1997) Stabilizing selection and the comparative analysis of adaptation. \emph{Evolution} \bold{51}: 1341--1351.

Hansen, T. F., Pienaar, J & Orzack, S, H. (2008) A comparative method for studying adaptation to a randomly evolving environment. \emph{Evolution} \bold{62}:1965--1967.

Labra, A., Pienaar, J. & Hansen, T. F.(2009) Evolution of thermophysiology in Liolaemus lizards: adaptation, phylogenetic inertia and niche tracking. \emph{The American Naturalist} \bold{174}: 204--220. 

Pienaar et al (in prep) An overview of comparative methods for testing adaptation to external environments.  

Hansen, T. F. and Bartoszeck, K. (in review) Interpreting the evolutionary regression: the interplay between observational and biological errors in phylogenetic comparative studies.
}
\keyword{ package }
\seealso{

\code{\link[<ape>:<ape>-package]{<ape>}}, \code{\link[<ouch>:<ouch>-package]{<ouch>}}
}

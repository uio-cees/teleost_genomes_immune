\name{model.fit}
\alias{model.fit}
\title{Main model fitting function for Slouch}
\description{
  Fits OU models with optima modeled on: no predictors (i.e. an intercept only regression, useful for estimating phylogenetic signal), fixed categorical predictors mapped onto the tree, fixed categorical and fixed continuous predictors (i.e. fixed predictor ANCOVA), fixed categorical and randomly evolving continuous predictors (i.e. mixed model ANCOVA), fixed continuous predictors (i.e. a fixed predictor regression), randomly evolving continuous predictors (i.e. a random predictor regression), any combination of fixed and randomly evolving continuous predictors (i.e. a mixed model regression). In all cases observation error for continuous predictors and continuous response variables can be incorporated, as well as observation error covariance  between the response and predictor variables. 
}
\usage{
model.fit(topology, times, half_life_values, vy_values, response, me.response = NULL, fixed.fact = NULL, fixed.cov = NULL,
me.fixed.cov = NULL, mecov.fixed.cov = NULL, random.cov = NULL, me.random.cov = NULL, mecov.random.cov = NULL, intercept = "root", ultrametric = TRUE, support = NULL, convergence = NULL)
}
\arguments{
  \item{topology}{A column in the slouch tree data.frame named 	\code{"ancestor"} that specifies which node a current species or internal node is connected to. The ancestral node for the node at the root of the tree gets a numeric value of 0.}
  \item{times}{A column in the slouch tree data frame named \code{"time"} specifying the cumulative time of the current node in relation to the root node, which occurs at time 0.}
  \item{half_life_values}{A user specified vector of candidate half-life values that are evaluated jointly with a user-specified vector of candidate Vy (stationary variance of the response trait) values by likelihood.}
  \item{vy_values}{A user specified vector of candidate Vy values that are evaluated jointly with a user-specified vector of candidate half-life values by likelihood.}
  \item{response}{The vector of trait values in the slouch data frame that are tested for adaptation to optima modeled on various niche variables.}
  \item{me.response}{An optional vector of observation error variances for the response variable.}
  \item{fixed.fact}{A categorical vector encoding fixed categorical niches to model the optima on. To gain the most out of the analysis, it is recommended that some form of ancestral state reconstruction of these variables on the phylogeny is attempted (i.e. using the \code{\link{fitch}} parsimony function, where the \code{fitch} object can be parsed directly to \code{model.fit}). Internal nodes will therefore also have niche state values. If a reconstruction is not attempted, all internal nodes can be encoded as a single ancestral niche.}
  \item{fixed.cov}{A continuous, numeric vector (or matrix for multiple predictors) of the same length as the number of tip species.  Encodes a numerical state for the optima to be be modeled on as a regression but assumes that the numeric state is fixed for a given lineage (i.e. has not evolved randomly on the phylogeny). To include multiple predictors use \code{cbind(covariate1, covariate2 .. etc)} in the function call.}
  \item{me.fixed.cov}{An optional vector (or matrix) of measurement error variance / s for the fixed covariate / s. To include multiple measurement variances,  use \code{cbind(m.e.1, m.e.2 .. etc)} in the function call. For multiple regression,  the number of entries must match the number of predictors.}
  \item{mecov.fixed.cov}{An optional vector (or matrix) of measurement error covariance / s for the fixed covariate / s with the response variable.To include multiple measurement covariances, use \code{cbind(m.e.cov1, m.e.cov2 .. etc)} in the function call. For multiple regression,  the number of entries must match the number of predictors.}
  \item{random.cov}{A continuous numeric vector (or matrix for multiple predictors)of the same length as the number of tip species.  Encodes a numerical state for each tip species for the optima to be be modeled on as a regression. Assumes that the numeric state evolves randomly on the phylogeny where each such predictor is modeled as a Brownian Motion process). To include multiple predictors use \code{cbind(covariate1, covariate2 .. etc)} in the function call.}
  \item{me.random.cov}{An optional vector (or matrix) of observation error variance / s for the randomly evolving covariate / s. To include multiple observation variances, use \code{cbind(m.e.1, m.e.2 .. etc)} in the function call. For multiple regression, the number of entries must match the number of predictors.}
  \item{mecov.random.cov}{An optional vector (or matrix) of observation error covariance / s for the random covariate / s with the response variable. To include multiple entries, use \code{cbind(m.e.cov1,m.e.cov2, .. etc)} in the function call. For multiple regression, the number of entries must match the number of predictors.}
  \item{intercept}{An optional control variable for specyfing how to model the intercept (ya) for the fixed, categorical niches (if included). It is rarely possible obtain a meaningful estimate of this parameter, as its coefficent often tends to 0. We recommend either combining this coefficient with the coefficent for the niche that mapped to the root of the tree the (default, where \code{intercept} is set to = \code{"root"}), or using some integer value (e.g. \code{intercept = 1}). To attempt to obtain an estimate, use \code{intercept = NULL}.}
  \item{ultrametric}{By default \code{model.fit} assumes an ultrametric tree for the mean structure of the model. By setting \code{ultrametric=FALSE}, the more complicated design matrix for the mean structure is used (see Hansen et al (2008), where two additional parameters need to be estimated). In either case the function will always use a variance matrix corresponding to the tree structure.}
  \item{support}{The support region plotted in the likelihood surface for the joint half-life, Vy vectors, by default includes everything within 2 log-likelihood units. This can be changed to any number of support units (n) by specifying \code{support = n}}
  \item{convergence}{The convergence criteria for the iterated GLS routine is by default 1e-06 but can be set to n by specifying \code{convergence = n} (note that the default value could be problematic if the variables themselves are in very small units).}
}
\details{
  The papers introducing the OU models for fixed, categorical predictors (Hansen 1997, Butler & King 2004) introduce and explain the theory behind these models. The paper introducing the OUBM model and user manual for SLOUCH v1.0 (Hansen et al, 2008) as well as an online appendix for Labra et al. (2011) explain the theory behind the models where the optima themselves evolve randomly and are modeled on continuous, randomly evolving predictors. Pienaar et al (submitted), provides an overview of the theory behind the OU models as well as descriptions of the extensions where fixed, categorical, fixed continuous  and randomly evolving continuous predictor variables can be combined to model optima. }
\value{
  \item{MODEL SUMMARY}{A list stating: the type of model being fit; the data set the model is being fit to; the response variable name; the predictor variable names (listed as NULL if not present in the model); the sample size (number of tip species); and the maximum tree height.}
  \item{GRID SEARCH PARAMETER SUPPORT}{A recurring list of the current half-life and Vy values being evaluated, the likelihood support for that half-life / Vy combination, and in the case that predictors are specified, the optimal mean or regression parameters (intercept terms (1 for ultrametric trees, 3 for non-ultrametric trees) and slope / s). This list serves two functions, the first being to identify the half-life and Vy values that fall within the 2-unit support region (along with the parameter estimates) and the second to provide a loop counter to determine the grid search progress. In addition this list can be used to determine for which combinations of half-life / Vy values, the likelihood cannot be determined and so aids in refining these vectors during the grid 			search.}
  \item{MODEL PARAMETERS}{A list of the best estimates for the parameters describing evolution on the phylogeny: the rate of adaptation towards changing optima (alpha); the phylogenetic half life (ln(2)/alpha) which is the time taken for the trait to evolve halfway towards it's optimum in expectation (a measure of the strength of phylogenetic inertia where its units are in terms of the branch lengths of the tree); the phylogenetic correction factor (another measure of phylogenetic inertia), used to correct the observed evolutionary regression to obtain the regression expected if adaptation were instantaneous); and the stationary variance (vy), the response variable secondary variance once the effects of the optima have been accounted for.}
  \item{PRIMARY OPTIMA}{A list of the best estimates for primary optima, where the parameters are in in terms of the model specified: for a model with no predictors, a single global optimum is returned; for a model with fixed, categorical niches estimated mean primary optima for those niches are returned; if fixed or randomly evolving covariates are included in the preceding model, these are returned as optimal regression slopes; if only fixed covariates are modeled, an optimal regression intercept and slope / s are returned; If only randomly evolving covariates are modeled, two sets of regressions are returned, an evolutionary and an optimal regression. The evolutionary regression is the observed, phlogenetically corrected regression slope (should be similar to estimates returned by phylogenetic GLS routines of other programs). The optimal regression is the slope corrected by the phylogenetic correction factor, and gives the expected slope if adaptation towards the randomly evolving optima were instantaneous. An optimal regression intercept (bo) is also returned based on the following assumption: the composite intercept K , which is compromised of the trait value at the root of the tree (ya), the predictor state at the root (x0) and the regression intercept (bo) can be decomposed by assuming that ya = xa (i.e. that the trait was at its optimum at the root of the tree). Note that if a model with fixed and randomly evolving covariates is specified, the fixed regression slope / s will appear in both the evolutionary and optimal regression slope sections, but will have the same value.}
  \item{MODEL FIT}{A list of model fit criteria for model selection and evaluation including: model support (Likelihood), Akaike information criteria (AIC), small sample size corrected AIC (AICc), the stricter Swartz information criteria (SIC), the phylogenetically corrected r square value for the percent variance explained and the total (SST) and error (SSE) sums of squares. A 3-dimensional perspective plot where the vectors of half-life and Vy values are given on the X and Y axes and the likelihood for each combination (scaled to 0 for the best estimate) on the Z axis is also returned in the plotting window.}
 
}
\references{ 
Butler, M. A. & King. A. A. (2004) Phylogenetic comparative analysis: a modeling approach for adaptive evolution. \emph{The american Naturalist} \bold{164}: 683--695.  

Hansen, T. F. (1997) Stabilizing selection and the comparative analysis of adaptation. \emph{Evolution} \bold{51}: 1341--1351.

Hansen, T. F., Pienaar, J & Orzack, S, H. (2008) A comparative method for studying adaptation to a randomly evolving environment. \emph{Evolution} \bold{62}:1965--1967.

Labra, A., Pienaar, J. & Hansen, T. F.(2009) Evolution of thermophysiology in Liolaemus lizards: adaptation, phylogenetic inertia and niche tracking. \emph{The American Naturalist} \bold{174}: 204--220. 

Pienaar et al (submitted) An overview of comparative methods for testing adaptation to external environments.  
 }
\author{Jason Pienaar \email{jasonpienaar@gmail.com}}
\note{ The model.fit function will not always work for a given combination of half-life and Vy values. For some combinations it is not possible to calculate the likelihood or the regression parameters due to numerical issues (e.g. for half-life values that are much larger than the actual values). In these cases, the trick is to refine these vectors (e.g. decrease or shift the range) until parameters can be estimated. Once a ballpark range of estimates and their support has been obtained by rough grids, the user should refine the increments between candidate half-life and Vy values until satisfied that round-off error is small enough.
}
\seealso{\code{\link{slouchbiascorr}}, \code{\link{slouchtree.plot}},\code{\link{sigma.X.estimate}},\code{\link{fitch}}, \code{\link{np.regression}}}
\examples{
data(Lio)
attach(Lio)
## NO OPTIMA (EVALUATING PHYLOGENTIC SIGNAL)

## Tsel modeled on grand mean (no predictors, no observation error) 
## mean structure assumed for ultrametric tree, variance structure for non-ultrametric
## half life is a measure of phylogenetic signal
h<-10:30*0.1
vy<-1:20
model.fit(ancestor, time, h, vy, Tsel)
## suggest strong phylogentis signal in Tsel (BM evolution)
## DO NOT YET conclude phylogentic inertia / BM evolution
## could be due to phylogentifc signal in the predictors defining Tsel's optima 

## test for phylogenetic signal in Tb
h<-0:6*0.1
vy<-1:15
model.fit(ancestor, time, h, vy, Tb)
## rough grid suggest half life at 0.1 (10 % tree height)
## more refined grid (h<-1:60*0.01) suggests 0.07 for half-life

## could also set non ultrametric design to estimate Ya and Bo
model.fit(ancestor, time, h, vy, Tb, ultrametric=FALSE)
## could also include observation error in this and preceding models for Tb
model.fit(ancestor, time, h, vy, Tb, me.response=mv.tb, ultrametric=FALSE)
## note how it quite seriously effects the best half-life and vy estimates

## MODELING OPTIMA 
## note these models take a long time to run depending on the size of the grid search (they are therefore not run through R CMD check)

## model Tsel optima on reconstructed climate niche
\dontrun{
h<-1:10*0.1
vy<-1:20*0.1
model.fit(ancestor, time, h, vy, Tsel, me.response=mv.tsel, fixed.fact=pars.climate)
## AICc = 93.8, r squared = 44 %, half-life = 0.4
## incorporate fixed covariate (e.g. assume Tsel scales mechanically with Tb) with observation error
model.fit(ancestor, time, h, vy, Tsel, me.response=mv.tsel, fixed.fact=pars.climate, fixed.cov=Tb, me.fixed.cov=mv.tb)
# AICc = 94.5, r squared = 49 %, half-life = 0.3
## assume Tb itself evolves randomly
model.fit(ancestor, time, h, vy, Tsel, me.response=mv.tsel, fixed.fact=pars.climate, random.cov=Tb, me.random.cov=mv.tb)
## AICc = 87, r squared = 64 %, half-life = 0.2
## model optima on multiple randomly evolving predictors (assume ultrametric mean structure)
model.fit(ancestor, time, h, vy, Tsel, me.response=mv.tsel, random.cov=cbind(Tb, To, Ts), me.random.cov=cbind(mv.tb, mv.to, mv.ts))
## AICc = 80, r squared = 70%, half-life = 0.1
## the best model is probably some combination of climate and the continuous predictors
## clearly the phylogenetic signal is not in Tsel, but in one of the predictors for its optima. 
## these models can also be run using ultrametric = false, and with any number and combination of predictors
}
detach(Lio)
}


\keyword{models}
\keyword{regression}

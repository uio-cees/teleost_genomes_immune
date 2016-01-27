\name{np.regression}
\alias{np.regression}
\title{Non-phylogentic iterated GLS Regression model fitting function for Slouch}
\description{
  Estimates regression model parameters when both response and predictor (or multiple predictors) have observation error variances and / or covariances. Does NOT take any phylogenetic information into account.
}
\usage{
np.regression(response, me.response = NULL, random.cov, me.random.cov = NULL, mecov.random.cov = NULL, convergence = NULL)
}
\arguments{
  \item{response}{The vector of trait values in the slouch data frame that are tested for adaptation to optima modeled on various niche variables.}
  \item{me.response}{An optional vector of observation error variances for the response variable.}
  \item{random.cov}{A continuous, numeric vector (or matrix for multiple predictors) of the same length as the number of tip species.  Encodes a numerical state for the response variable to be be modeled on as a regression. To include multiple predictors use 	cbind(covariate1, covariate2, .. etc) in the function call.}
  \item{me.random.cov}{An optional vector (or matrix) of observation error variance / s for the randomly evolving covariate / s. To include multiple entries, use \code{cbind(m.e.1, m.e.2, .. etc)} in the function call. For multiple regression, the number of entries must match the number of predictors.}
  \item{mecov.random.cov}{An optional vector (or matrix) of observation error covariance / s for the random covariate / s with the response variable.To include multiple entries, use \code{cbind(cov1, cov2 .. etc)} in the function call. For multiple regression, the number of entries must match the number of predictors.}
  \item{convergence}{The convergence criteria for the iterated GLS routine is by default 1e-06 but can be set to n by specifying \code{convergence = n} (note that the default value could be problematic if the variables themselves are in very small units).}
}
\details{
  \code{np.regression} works in the same way as most \pkg{slouch} functions, taking arguments directly from the \pkg{slouch} data frame. 
}
\value{
  \item{REGRESSION PARAMETERS}{A list of the best estimates for the non phylogenetically corrected regression intercept and slope / s along with the standard errors of these estimates.}
  \item{MODEL FIT}{A list of model fit criteria for model selection and evaluation including: model support (Likelihood), Akaike information criteria (AIC), small sample size corrected AIC (AICc), the stricter Swartz information criteria (SIC), the r square value for the percent variance explained and the total (SST) and error (SSE) sums of squares.}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jason Pienaar \email{jasonpienaar@gmail.com}}
\seealso{\code{\link{model.fit}}, \code{\link{slouchbiascorr}}}
\examples{
data(Lio)
attach(Lio)
## single predictor y~x regression, no obervation error
## example Tsel regressed on Tb (will give same estimates as base lm function)
np.regression(Tsel, random.cov=Tb)
## include observation error for response varible
np.regression(Tsel, me.response=mv.tsel, random.cov=Tb)
## include observation error for predictor variable
np.regression(Tsel, me.response=mv.tsel, random.cov=Tb)
np.regression(Tsel, me.response=mv.tsel, random.cov=Tb, me.random.cov=mv.tb)
## note could also include observation error covariance between predictor and response observation errors. not shown in this example.
## multiple regression, no observation error
## Tsel regressed on body temperature and substrate temperature
np.regression(Tsel, random.cov=cbind(Tb, Ts))
## multiple regression with observation errors for all variable
np.regression(Tsel, me.response=mv.tsel, random.cov=cbind(Tb, Ts), me.random.cov=cbind(mv.tb, mv.ts))
detach(Lio)
}
\keyword{models}
\keyword{regression}


% --- Source file: slouchbiascorr.R ---
\name{slouchbiascorr}
\alias{slouchbiascorr}

\title{Calculate phylogenetic regression bias caused by observation error in continuous predictors}
\description{
  \code{Slouchbiascorr} calculates the degree of bias in the regression slopes returned by \code{model.fit} when the predictor variables are associated with observation error variances.  
}
\usage{
slouchbiascorr(beta, topology, times, half_life_value, vy, response, me.response = NULL, random.cov, me.random.cov = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{ a user-created matrix of the optimal regression estimates returned by \code{model.fit}. Use \code{matrix(cbind(K, slope1, slope2, ..etc))} to create this matrix. (NB \code{matrix} must be used to get the correct dimensions)}
  \item{topology}{A column in the slouch tree data.frame named 	\code{"ancestor"} that specifies which node a current species or internal node is connected to. The ancestral node for the node at the root of the tree gets a numeric value of 0.}
  \item{times}{A column in the slouch tree data frame named \code{"time"} specifying the cumulative time of the current node in relation to the root node, which occurs at time 0.}
  \item{half_life_value}{The best half-life estimate returned by \code{model.fit} }
  \item{vy}{The best Vy estimate returned by \code{model.fit}}
  \item{response}{The vector of trait values in the slouch data frame that are tested for adaptation to optima modeled on various niche variables.}
  \item{me.response}{An optional vector of observation error variances for the response variable.}
  \item{random.cov}{The same random covariate / s used in code{model.fit}}
  \item{me.random.cov}{ the same observation error / s for the random covariate / s used in \code{model.fit} }
}
\details{
  Hansen and Bartoszeck (2012) provide the theoretical underpinnings for this function
}
\value{
  \code{slouchbiasscorr} returns a matrix where the first row column contains the uncorrected regression estimates, the second column contains the bias correction factors, the third column the bias corrected regression estimates and the last column the standard errors for these bias corrected estimates}
\references{
Hansen, T. F. and Bartoszeck, K. (in review) Interpreting the evolutionary regression: the interplay between observational and biological errors in phylogenetic comparative studies.
}
\author{Jason Pienaar \email{jasonpienaar@gmail.com} and Kzrysztof Bartoszeck \email{krzbar@chalmers.se}}
\note{The theory for bias corrections when observation errors are correlated with each other has not been developed yet, hence this option does not yet exist in \code{slouchbiascorr} 
}
\seealso{ \code{\link{model.fit}}}
\examples{
data(Lio)
attach(Lio)
## fit a model for Tsel where optima are modeled on Ts and Tb
## Ts and Tb include observation error which could bias the regression
h<-0.3
vy<-5
model.fit(ancestor, time, h, vy, Tsel, me.response=mv.tsel, random.cov=cbind(Ts,Tb), me.random.cov=cbind(mv.ts, mv.tb))
## optimal regression estiamtes are K=25.90157, Ts = 0.23788, Tb = 0.27057
## create beta matrix
beta<-matrix(cbind(25.90157,0.23788, 0.27057))
## use slouchbiascorr to determine and correct for bias
slouchbiascorr(beta, ancestor, time, 0.3, 5, Tsel, mv.tsel, cbind(Ts,Tb),  cbind(mv.ts, mv.tb))
detach(Lio)
}


\keyword{models}
\keyword{regression}


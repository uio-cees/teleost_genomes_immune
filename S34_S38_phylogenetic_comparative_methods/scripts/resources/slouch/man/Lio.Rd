\name{Lio}
\alias{Lio}
\docType{data}
\title{Liolaemus Lizards thermoregulatory data}
\description{
  This is part the Liolaemus lizard data set used in Labra et al. (2009). It was used to model evolving environmentally determined optima for the body temperature these organisms choose in a thermal gradient (Tsel). Tsel is assumed to be the optimal body temperature for a given species in their past and present environments. The \code{data.frame} is given in a \pkg{slouch} \code{data.frame} format, ready to use by \pkg{slouch}. It includes columns that describe the composite phylogenetic tree (\code{nodes}, \code{species}, \code{ancestor}, \code{time}. The rest of the columns are the thermoregulatory data itself, including various continuous measures (and their observation error variances) and categorical encodings of the lizards and their environments}
\usage{data(Lio)}
\format{
  A data frame with 59 observations (30 species and 29 internal nodes) on the following 16 variables.
  \describe{
    \item{\code{nodes}}{A numeric vector numbering the nodes (root node = 1)}
    \item{\code{species}}{Species names for extant species, internal nodes are encoded as missing values}
    \item{\code{ancestor}}{A numeric vector encoding ancestral nodes that \code{nodes} are connected to}
    \item{\code{time}}{Cumulative branch lengths in terms of nucleotide substitutions scaled by maximum tree height}
    \item{\code{Tsel}}{Selected body temperatures in degrees centigrade}
    \item{\code{mv.tsel}}{Observation error variances for \code{Tsel}}
    \item{\code{Tb}}{Mean body temperatures during heating and cooling experiments (in degrees centigrade)}
    \item{\code{mv.tb}}{Observation error variances for \code{Tb}}
    \item{\code{svl}}{Snout vent length in millimeters}
    \item{\code{mv.svl}}{Observation error variances for \code{svl}}
    \item{\code{Ts}}{Mean substrate temperatures}
    \item{\code{mv.ts}}{Observation error variances for \code{Ts}}
    \item{\code{To}}{Mean operative temperature}
    \item{\code{mv.to}}{Observation error variances for \code{To}}
    \item{\code{climate}}{a niche encoding factor with levels \code{Continental} \code{Mediterranean} \code{Oceanic} \code{Tropical} \code{anc}. Here \code{anc} means that a single niche is encoded for all ancestral lineages, these can be estimated by parsimony using \code{\link{fitch}}, or some other ancestral reconstruction program}
    \item{\code{temp}}{a factor with levels \code{Cold} \code{Warm} \code{anc}}
  }
}
\details{
  Refer to Labra et al. (2009) for details on how the measurements were taken and for other variables used in the analyses.
}
\source{
  Labra, A., Pienaar, J. & Hansen, T. F.(2009) Evolution of thermophysiology in Liolaemus lizards: adaptation, phylogenetic inertia and niche tracking. \emph{The American Naturalist} \bold{174}: 204--220.
}
\references{
  Lobo, F. (2001) A phylogenetic analysis of lizards of the Liolaemus chiliensis group (Iguania: Tropiduridae). \emph{Herpetological Journal} \bold{11}: 137--150. 
Schulte II, J. A., J. R. Macey, R. E. Espinoza, and A. Larson. (2000) Phylogenetic relationships in the iguanid lizard genus Liolaemus: multiple origins of  viviparous reproduction and evidence for recurring Andean vicariance and dispersal. \emph{Biological Journal Linnean Society} \bold{69}: 75--102.  
}
\examples{
data(Lio)
## attach the data set for easier use with slouch
}
\keyword{datasets}

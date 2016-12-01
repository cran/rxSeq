\name{nLogLik}
\alias{nLogLik}
\title{Negative log likelihood for coefficients provided in results of the fit using process function}
\description{
  Calculates negative log(likelihood) of an X chromosome joint TReC and ASE counts model at a given set of parameters
}
\usage{
nLogLik(res, rc, genei, hessian=FALSE)}
\arguments{
  \item{res}{result object from process function}
  \item{rc}{Read count data object created by readCounts function}
  \item{genei}{get results for i'th gene}
  \item{hessian}{a logical option whether to calculate a Hessian matrix, the default values is set to FALSE.}
}

\value{
  output - list(nll=-log.likelihood,hessian=hessian matrix)
}
\author{ Vasyl Zhabotynsky vasyl@unc.edu }
\seealso{\code{\link{process}}, \code{\link{rcA}}, \code{\link{readCounts}}.}
\keyword{ methods }

\examples{\dontrun{
# get negative-log likelihood at the given point:
nLogLik(res=trecase.X.out, rc=rcX, genei=1, hessian=TRUE)
}}

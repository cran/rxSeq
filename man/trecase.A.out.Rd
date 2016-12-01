\name{trecase.A.out}
\alias{trecase.A.out}
\title{Example of results produced by optimizing step using process function on autosomal genes. Structured as a list.}
\description{
  A list containing test results as well as parameter estimates for joint model evaluated by process function for autosomal genes.
}
\value{
  a list of following matrices (if there is only one sex, only the relevant tests and matrices are outputed) :
  \item{pvals}{matrix of p-values from description for each gene corresponding row}
  \item{coef.full}{matrix of full model fit coefficients, -log(likelihood at these coefficients),phi, theta (2 overdispersion parameters used)}  
  \item{coef.add}{matrix of additive restricted fit coefficients, -log(likelihood at these coefficients),phi, theta}
  \item{coef.poo}{matrix of parent of origin restricted fit coefficients, -log(likelihood at these coefficients),phi, theta}
  \item{coef.dom}{matrix of dominance restricted fit coefficients, -log(likelihood at these coefficients),phi, theta}
  \item{coef.same}{matrix of TReC=ASE additive restricted fit coefficients, -log(likelihood at these coefficients),phi, theta}
  \item{coef.ase.add}{matrix of ASE additive restricted fit coefficients, -log(likelihood at these coefficients),phi, theta}
  \item{coef.sex}{matrix of sex restricted fit coefficients, -log(likelihood at these coefficients),phi, theta}
  \item{coef.sex.add}{matrix of sex specific additive restricted fit coefficients, -log(likelihood at these coefficients),phi, theta}
  \item{coef.dev.dom}{matrix of dominance deviation for male restricted fit coefficients, -log(likelihood at these coefficients),phi, theta}
  \item{errorlist}{a list of errors}
  \item{hess.lst}{a list of heassian matrices, if parameter \bold{hessian} is set to TRUE}  
}
\author{ Vasyl Zhabotynsky vasyl@unc.edu }
\seealso{\code{\link{process}}, \code{\link{rcA}}, \code{\link{data.A}}.}
\keyword{ methods }

\examples{
names(trecase.A.out)
}

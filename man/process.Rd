\name{process}
\alias{process}
\title{Optimization wrapper, maximizing either the joint model of total (TReC) and allele specific (ASE) counts or just TReC}
\description{
  Performs optimization of one of four combinations: joint TReC and ASE or just TReC for autosome or X chromosome and tests with lrt test several hypotheses:
  additive, parent of origin, dominance, consistency of TreC and ASE additive effect, ASE only additive effect, 
  sex, sex specific additive, dominance deviation for males.
}
\usage{
process(rc, hessian=FALSE)}
\arguments{
  \item{rc}{an object of class readCounts.}
  \item{hessian}{a flag whether Hessian matrix for these genes should be calculated, by default set to FALSE}
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
\seealso{\code{\link{get.tausB}},\code{\link{nLogLik}}, \code{\link{data.X}}, \code{\link{data.A}}, \code{\link{rcA}}, \code{\link{readCounts}}.}
\keyword{ methods }

\examples{\dontrun{
# fitting X chromosome data example, for a full model, i.e. assuming we have allele specific reads:
trecase.A.out = process(rc=rcA)
names(trecase.A.out)
trecase.A.out$pval

#alternatively for X chromosome:
trecase.X.out = process(rc=rcX)
names(trecase.X.out)
trecase.X.out$pval
}}

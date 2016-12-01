\name{get.tausB}
\alias{get.tausB}
\title{Produce \emph{Xce} estimates for mice with allele specific reads}
\description{
  \emph{Xce} estimation for mice with allele specific reads.
}
\usage{
get.tausB(n, n0B, geneids, min.cnt=50, exclude.prop=.05, Xist.ID="ENSMUSG00000086503")}
\arguments{
  \item{n}{vector of allele specific counts for each mouse}  
  \item{n0B}{vector of allele specific counts for allele B}
  \item{geneids}{gene IDs}
  \item{min.cnt}{minimum number of allele specific counts}
  \item{exclude.prop}{minimum proportion of allele specific counts for each allele}
  \item{Xist.ID}{and ID of \emph{Xist}, to exclude it from estimating \emph{Xce}, since \emph{Xce} would 1-tausB}
}

\value{
  output - matrix of 4 rows: 
  \item{med.tauB}{taus estimated via median}
  \item{ave.tauB}{taus estimated via percent of allele B counts}
  \item{all.genes}{number of genes that had passed minimum count}
  \item{used.genes}{number of genes that had required percent of each allele}
  each column represent respective mouse.
}
\author{ Vasyl Zhabotynsky vasyl@unc.edu }
\seealso{\code{\link{process}},\code{\link{nLogLik}}, \code{\link{data.X}}, \code{\link{rcX}}.}
\keyword{ utilities }


\examples{
# Estimating XCE effect for each mouse for X chromosome
get.tausB(n=data.X$n, n0B=data.X$n0B, geneids=data.X$geneids)
}

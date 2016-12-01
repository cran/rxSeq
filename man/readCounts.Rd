\name{readCounts}
\alias{readCounts}
\title{A list object that should be used as input to optimization process function.}
\description{
  It should contain at least total read counts (TReC) and classification of crosses 1 to 8. To fit the full model should also have appropriate allele specific counts n and n0B. Also is used along with results of optimization as input to nLogLik function if one needs to calculate Hessian matrix.
}
\value{  
  \item{index}{vector defining the cross of the mouse, female - AB=1, BA=2, AA=3, BB=4, and male - AB=5, BA=6, AA=7, BB=8.
  If mice are of only one sex, AB=1, BA=2, AA=3, BB=4.}
  \item{y}{matrix of TReC counts. Note, the expected input assumes that inbred mice will be in the last columns of the table, after the last F1 mouse.}
  \item{n}{matrix of ASE counts for corresponding F1 mouse (classes 1,2,5,6) for corresponding genes.}  
  \item{n0B}{matrix of ASE counts belonging for allele B, for correponding genes and mice as in \bold{n}.}  
  \item{kappas}{A parameter, specifying as overall TReC for the mouse, on log scale}
  \item{tausB}{\emph{Xce} effect: expression of allele B relative to the overall allele specific count for each mouse. Set to NULL in autosomes.}
  \item{gene.switch}{For which genes \emph{Xce} effect should be switched. Null for autosomes.}
  \item{geneids}{ids of genes, if not provided, rownames of the matrix \bold{y} will be used}
  \item{chrom}{this field would be set to be "X" since this is dataset for chromosome X}
  \item{model}{set to be "full", can be modified to "short" to run a TReC oply model}
  \item{geneids}{Ensembl gene ids}  
  \item{tech.ctrl}{a list of overdispersion boundaries and log(2)}
}
\author{ Vasyl Zhabotynsky vasyl@unc.edu }
\keyword{ utilities }
\seealso{\code{\link{process}},\code{\link{nLogLik}}, \code{\link{simRX}}.}

\examples{
# see total read counts (TReC) for first 2 X chromosome genes of a data example:
rcA = readCounts(index=data.A$index, y=data.A$y[1:2,], n=data.A$n[1:2,], n0B=data.A$n0B[1:2,], 
                 kappas=data.A$kappas, geneids=data.A$geneids[1:2])
}
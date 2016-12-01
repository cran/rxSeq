\name{data.X}
\alias{data.X}
\title{Sample data example for X chromosome genes}
\description{
  This data set provides with example of experimental data for a subset of autosomal genes.
  The full model requires a combination of total read counts (\bold{y}) - all the reads belonging for a gene, 
  and finding out which of these reads we can specifically attribute to allele A or allele B - allele specific counts (\bold{n}),
  separately the reads attributed specifically to allele B (\bold{n0B}).
  Also, \bold{tausB} - is the \emph{Xce} effect for each F1 mouse, which specifies the proportion of allele specific reads belonging to allele B.
  Also, it includes the other data pieces to fit the model: \bold{kappas} - total number of counts for each mouse, on log scale, 
  index - specifying which cross each mouse belongs to, and geneids - Ensembl ids of genes.
  They, as well as the datasets simulated with \bold{simRX} can be fitted using \bold{proc.trecase.X} or \bold{proc.trec.X}.  
}
\value{  
  \item{index}{vector defining the cross of the mouse, female - AB=1, BA=2, AA=3, BB=4, and male - AB=5, BA=6, AA=7, BB=8.
  If mice are of only one sex, AB=1, BA=2, AA=3, BB=4.}
  \item{y}{matrix of TReC counts. Note, the expected input assumes that inbred mice will be in the last columns of the table, after the last F1 mouse.}
  \item{n}{matrix of ASE counts for corresponding F1 mouse (classes 1,2,5,6) for corresponding genes.}  
  \item{n0B}{matrix of ASE counts belonging for allele B, for correponding genes and mice as in \bold{n}.}  
  \item{kappas}{A parameter, specifying as overall TReC for the mouse, on log scale}
  \item{tausB}{\emph{Xce} effect: expression of allele B relative to the overall allele specific count for each mouse. Use some allele specific counts to estimate the effect.}
  \item{geneids}{ids of genes, if not provided, rownames of the matrix \bold{y} will be used}

  \item{index}{index defining which mouse belongs to which cross}
  \item{y}{modified total read counts}
  \item{n}{modified allele specific counts}
  \item{n0B}{modified allele specific coutns, belonging to allele B}
  \item{kappas}{offset, defining a library size for each mouse}
  \item{tausB}{\emph{Xce} effect for each mouse, for a given cross}
  \item{geneids}{Ensembl gene ids}  
}
\author{ Vasyl Zhabotynsky vasyl@unc.edu }
\keyword{ utilities }
\seealso{\code{\link{process}},\code{\link{readCounts}}.}

\examples{
# see total read counts (TReC) for first 2 X chromosome genes of a data example:
data.X$y[1:2,]
}
\name{rcX}
\alias{rcX}
\title{Reformatted data for chromosome X set to be used as input to process function}
\description{
  This is an object of type readCounts provides with example of experimental data for a subset of X chromosome genes.
  The full model requires a combination of total read counts (\bold{y}) - all the reads belonging for a gene, 
  and finding out which of these reads we can specifically attribute to allele A or allele B - allele specific counts (\bold{n}),
  separately the reads attributed specifically to allele B (\bold{n0B}).
  Also, \bold{tausB} - is the \emph{Xce} effect for each F1 mouse, which specifies the proportion of allele specific reads belonging to allele B.
  Also, it includes the other data pieces to fit the model: \bold{kappas} - total number of counts for each mouse, on log scale, 
  index - specifying which cross each mouse belongs to, and geneids - Ensembl ids of genes.
  They, as well as the datasets simulated with \bold{simRX} can be fitted using \bold{process} with appropriate options chrom="X" and field model to be either "full" or "short".  
}
\value{  

                genes.switch=genes.switch,
                geneids=geneids,
                
  \item{index}{vector defining the cross of the mouse, female - AB=1, BA=2, AA=3, BB=4, and male - AB=5, BA=6, AA=7, BB=8.
  If mice are of only one sex, AB=1, BA=2, AA=3, BB=4.}
  \item{y}{matrix of TReC counts. Note, the expected input assumes that inbred mice will be in the last columns of the table, after the last F1 mouse.}
  \item{n}{matrix of ASE counts for corresponding F1 mouse (classes 1,2,5,6) for corresponding genes.}  
  \item{n0B}{matrix of ASE counts belonging for allele B, for correponding genes and mice as in \bold{n}.}  
  \item{kappas}{A parameter, specifying as overall TReC for the mouse, on log scale}
  \item{tausB}{\emph{Xce} effect: expression of allele B relative to the overall allele specific count for each mouse. Use some allele specific counts to estimate the effect.}
  \item{gene.switch}{For which genes \emph{Xce} effect should be switched. Xist gene set to be switched in this set.}
  \item{geneids}{ids of genes, if not provided, rownames of the matrix \bold{y} will be used}
  \item{chrom}{this field would be set to be "X" since this is dataset for chromosome X}
  \item{model}{set to be "full", can be modified to "short" to run a TReC oply model}
  \item{geneids}{Ensembl gene ids}  
  \item{tech.ctrl}{a list of overdispersion boundaries and log(2)}
}
\author{ Vasyl Zhabotynsky vasyl@unc.edu }
\keyword{ utilities }
\seealso{\code{\link{process}},\code{\link{readCounts}},\code{\link{simRX}}.}

\examples{
# see total read counts (TReC) for first 2 X chromosome genes of a data example:
rcX
}
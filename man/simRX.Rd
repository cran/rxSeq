\name{simRX}
\alias{simRX}
\title{Produce simulated counts}
\description{
  This function is producing simulated counts for the joint model with Negative-Binomial distribution for TReC and Beta-Binomial for ASE counts. The simulated dataset should be reformatted to readCounts format to be used for optimization.
}
\usage{
simRX(b0f, b0m, b1f, b1m, beta_sex, beta_dom, beta_k=1, phi=1, theta=1, n=6, 
      mean.base.cnt=50, range.base.cnt=60, perc.ase=.35, n.simu=1E4, 
      is.X=FALSE, tauB=NULL, seed=NULL)}
\arguments{
  \item{b0f}{a female additive strain effect}  
  \item{b0m}{a male additive strain effect}  
  \item{b1f}{a female parent of origin effect}  
  \item{b1m}{a male parent of origin effect}  
  \item{beta_sex}{a sex effect}
  \item{beta_dom}{a dominance effect}
  \item{beta_k}{an effect associated with the library size kappas}
  \item{phi}{a Negative-Binomial overdispersion, default value is 1}
  \item{theta}{a Beta-Binomial overdispersion, default value is 1}
  \item{n}{a vector defining number of mice in each cross, default value is 6}
  \item{mean.base.cnt}{a target expected number of counts for the base group (with no effects), default value is 50}
  \item{range.base.cnt}{a range in which the expected number of counts  for the base group will vary, default value is 60}
  \item{perc.ase}{a percent reads that are allele-specific, default value is 35\%}
  \item{n.simu}{a number of simulations, default value is 1E4}
  \item{is.X}{a flag if the value to be simulated is X for chromosome (otherwise autosome), default value is FALSE}
  \item{tauB}{a value describing allelic imbalance - \emph{Xce} effect for the cross, default value is NULL, in which case 50\% will be simulated}
  \item{seed}{a random seed to be set, no set by default.}
}

\value{
  output - 3 matrices with one row - one gene, one column - one mouse: 
  \item{index}{vector defining the cross of the mouse, female - AB=1, BA=2, AA=3, BB=4, and male - AB=5, BA=6, AA=7, BB=8.
  If mice are of only one sex, AB=1, BA=2, AA=3, BB=4.}
  \item{y}{A matrix of total read counts}
  \item{n}{A matrix of allele specific counts}
  \item{n0B}{A matrix of allele specific counts associated with allele B}
  \item{kappas}{Offset parameter, given as overall TReC for the mouse.}
  \item{tausB}{In case of the simulating X chromosome the provided \emph{Xce} effect is returned: expression of allele B relative to the overall allele specific count for each mouse. }
}
\author{ Vasyl Zhabotynsky vasyl@unc.edu }
\seealso{\code{\link{process}},\code{\link{readCounts}}.}
\keyword{ utilities }

\examples{
# simulating autosomal data:
dat.A = simRX(b0f=.5, b0m=.6, b1f=.3, b1m=.4, beta_sex=.1, beta_dom=.1, n.simu=1E1)
names(dat.A)
# simulating autosomal data:
dat.X = simRX(b0f=.5, b0m=.6, b1f=.3, b1m=.4, beta_sex=.1, beta_dom=.1, n.simu=1E1, 
              is.X=TRUE, tauB=.3)
names(dat.X)
}

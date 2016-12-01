`ll.jRCI.sex.b1.A` = 
function(par, yi, ind.lst, X, sex, ni, ni0, xs, iphi, theta)
{
  logL = ll.tRCI.A(par0=c(par[1:5], par[-(1:4)]), yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, iphi=iphi)
  logL = logL + ll.aRC(par0=c(par[3], par[5]), ni[sex == 1], ni0[sex == 1], xs[sex == 1], theta)
  logL = logL + ll.aRC(par0=c(par[4], par[5]), ni[sex == -1], ni0[sex == -1], xs[sex == -1], theta)
  
  return(logL)
}

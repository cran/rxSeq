`ll.jRCI.b1.A2` = 
function(par, yi, ind.lst, X, sex=NULL, ni, ni0, xs, iphi, theta)
{
  seq = 1:4
  logL = ll.tRCI.A(par0=c(par[seq], 0, 0, par[-seq]), yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, iphi=iphi)
  logL = logL + ll.aRC(par0=c(par[3], 0), ni[sex == 1], ni0[sex == 1], xs[sex == 1], theta)
  logL = logL + ll.aRC(par0=c(par[4], 0), ni[sex == -1], ni0[sex == -1], xs[sex == -1], theta)

  return(logL)
}
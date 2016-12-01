`ll.jRCI.b0.A2` = 
function(par, yi, ind.lst, X, sex=NULL, ni, ni0, xs, iphi, theta)
{
  logL = ll.tRCI.A(par0=c(0, 0, 0, 0, par), yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, iphi=iphi)
  logL = logL + ll.aRC(par0=c(0, par[1]), ni[sex == 1], ni0[sex == 1], xs[sex == 1], theta)
  logL = logL + ll.aRC(par0=c(0, par[2]), ni[sex == -1], ni0[sex == -1], xs[sex == -1], theta)

  return(logL)
}
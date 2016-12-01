`ll.jRCI.b0.A1` = 
function(par, yi, ind.lst, X, ni, ni0, xs, iphi, theta)
{
  logL = ll.tRCI.A(par0=c(0, 0, par), yi=yi, ind.lst=ind.lst, X=X, twosex=FALSE, iphi=iphi)
  logL = logL + ll.aRC(par0=c(0, par[1]), ni, ni0, xs, theta)
  return(logL)
}
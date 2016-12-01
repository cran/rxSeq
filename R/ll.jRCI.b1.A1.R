`ll.jRCI.b1.A1` = 
function(par, yi, ind.lst, X, ni, ni0, xs, iphi, theta)
{
  seq = 1:2
  logL = ll.tRCI.A(par0=c(par[seq], 0, par[-seq]), yi=yi, ind.lst=ind.lst, X=X, twosex=FALSE, iphi=iphi)
  logL = logL + ll.aRC(par0=c(par[2], 0), ni, ni0, xs, theta)

  return(logL)
}
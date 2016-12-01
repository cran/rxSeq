`ll.jRCI.sex.A` = 
function(par, yi, ind.lst, X, ni, ni0, xs, iphi, theta)
{
  logL = ll.tRCI.A(par0=c(rep(par[1:3], each=2), par[-(1:3)]), yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, iphi=iphi)
  logL = logL + ll.aRC(par0=par[2:3], ni=ni, ni0=ni0, xs=xs, theta=theta)
  
  return(logL)
}

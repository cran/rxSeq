`ll.jRCI.sex.X` = 
function(par, yi, ind.lst, X, ni, ni0, xs, iphi, theta, l.tau.r, l.tau.a)
{
  logL = ll.tRCI.X(par0=c(par[1], par), yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, iphi=iphi, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
  logL = logL + ll.aRC.X(par0=par[2:3], ni=ni, ni0=ni0, xs=xs, theta=theta, l.tau.r=l.tau.r)
  
  return(logL)
}

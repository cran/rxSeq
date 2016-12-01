`ll.jRCI.trec.ase.X` = 
function(par, yi, ind.lst, X, twosex, ni, ni0, xs, iphi, theta, l.tau.r, l.tau.a)
{
  if(twosex){
    ind = c(1, 3)
    logL = ll.tRCI.X(par0=c(par[1:2], par[-2]), yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
  }else{
    ind = c(1, 2)    
    logL = ll.tRCI.X(par0=c(par[1], par), yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
  }
  logL = logL + ll.aRC.X(par0=par[ind], ni=ni, ni0=ni0, xs=xs, theta=theta, l.tau.r=l.tau.r)

  return(logL)
}
                     
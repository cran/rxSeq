`ll.jRCI.b1.X` = 
function(par, yi, ind.lst, X, twosex, ni, ni0, xs, iphi, theta, l.tau.r, l.tau.a)
{
  if(twosex){
    seq = 1:3
    ind = 3
  }else{
    seq = 1:2
    ind = 2
  }
  logL = ll.tRCI.X(par0=c(par[seq], 0, par[-seq]), yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
  logL = logL + ll.aRC.X(par0=c(par[ind], 0), ni=ni, ni0=ni0, xs=xs, theta=theta, l.tau.r=l.tau.r)

  return(logL)
}

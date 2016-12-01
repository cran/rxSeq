`ll.jRCI.b0.X` = 
function(par, yi, ind.lst, X, twosex, ni, ni0, xs, iphi, theta, l.tau.r, l.tau.a)
{
  if(twosex){
    zeroes = c(0, 0, 0)
  }else{
    zeroes = c(0, 0)
  }
  logL = ll.tRCI.X(par0=c(zeroes, par), yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
  logL = logL + ll.aRC.X(par0=c(0, par[1]), ni=ni, ni0=ni0, xs=xs, theta=theta, l.tau.r=l.tau.r)

  return(logL)
}

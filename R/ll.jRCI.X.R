`ll.jRCI.X` = 
function(par, yi, ind.lst, X, twosex, ni, ni0, xs, l.tau.r, l.tau.a, iphi, theta){
  logL = ll.tRCI.X(par0=par, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
  if(twosex){
    ase = c(3, 4)
  }else{
    ase = c(2, 3)
  }
  logL = logL + ll.aRC.X(par0=par[ase], ni=ni, ni0=ni0, xs=xs, theta=theta, l.tau.r=l.tau.r)

  return(logL)
}

`ll.jRCI.b0.A` = 
function(par, yi, ind.lst, X, twosex, sex=NULL, ni, ni0, xs, iphi, theta)
{
  if(twosex){
    zeroes = c(0, 0, 0, 0)
    logL = ll.aRC(par0=c(0, par[1]), ni[sex == 1], ni0[sex == 1], xs[sex == 1], theta)
    logL = logL + ll.aRC(par0=c(0, par[2]), ni[sex == -1], ni0[sex == -1], xs[sex == -1], theta)
  }else{
    zeroes = c(0, 0)
    logL = ll.aRC(par0=c(0, par[1]), ni, ni0, xs, theta)
  }
  logL = logL + ll.tRCI.A(par0=c(zeroes, par), yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi)

  return(logL)
}
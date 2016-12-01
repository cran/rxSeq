`ll.jRCI.b1.A` = 
function(par, yi, ind.lst, X, twosex, sex=NULL, ni, ni0, xs, iphi, theta)
{
  if(twosex){
    zeroes = c(0, 0)
    seq = 1:4
    logL = ll.aRC(par0=c(par[3], 0), ni[sex == 1], ni0[sex == 1], xs[sex == 1], theta)
    logL = logL + ll.aRC(par0=c(par[4], 0), ni[sex == -1], ni0[sex == -1], xs[sex == -1], theta)
  }else{
    zeroes = 0
    seq = 1:2
    logL = ll.aRC(par0=c(par[2], 0), ni, ni0, xs, theta)
  }
  logL = logL + ll.tRCI.A(par0=c(par[seq], zeroes, par[-seq]), yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi)

  return(logL)
}
`ll.aRC` = 
function(par0, ni, ni0, xs, theta){
  mu = par0[1] + par0[2] * xs
  lpi  = mu - log1p(exp(mu))
  logL =  -logBB(ni, ni0, exp(lpi), theta)
  
  return(logL)
}

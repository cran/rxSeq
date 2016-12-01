# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# negative log likelihood of ASE data from reciprical cross
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
`ll.aRC.X` = 
function(par0, ni, ni0, xs, theta, l.tau.r){#tausB
  mu = par0[1] + par0[2] * xs + l.tau.r
  lpi  = mu - log1p(exp(mu))
  logL =  -logBB(ni, ni0, exp(lpi), theta)
  
  return(logL)
}

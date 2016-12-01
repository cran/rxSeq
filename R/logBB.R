`logBB` = 
function(ni, ni0, pi, theta){
  alpha = pi/theta
  beta  = (1 - pi)/theta
  
  logL = sum(
          lchoose(ni, ni0) + lbeta(ni0 + alpha, ni - ni0 + beta) - lbeta(alpha, beta)
          )
  return(logL)
}
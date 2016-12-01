`loglikNB0` = 
function(subi, iphi, lmu, y, index){
  whch = (index == subi)
  yp = y[whch]
  lmup = lmu[whch]
  logL = sum(
  (lgamma(yp + iphi) -  lgamma(iphi) - lgamma(yp + 1) +  yp * lmup)  - 
  (iphi + yp) * log(iphi + exp(lmup)) + iphi * log(iphi)
  )
  
  return(logL)
}
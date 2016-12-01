`loglikNB2` = 
function(phi, lmu, y){
  iphi = 1/phi

  logL = sum(
  (lgamma(y + iphi) -  lgamma(iphi) - lgamma(y + 1) +  y * lmu)  - 
  (iphi + y) * log(iphi + exp(lmu)) + iphi * log(iphi)
  )
  
  return(logL)
}
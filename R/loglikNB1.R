`loglikNB1` = 
function(iphi, lmu, y){
  logL = sum(
  (lgamma(y + iphi) - lgamma(iphi) - lgamma(y + 1) +  y * lmu) - 
  (iphi + y) * log(iphi + exp(lmu)) + iphi * log(iphi)
  )
  
  return(logL)
}
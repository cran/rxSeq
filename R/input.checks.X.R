`input.checks.X` = 
function(coef, phi, index, y, kappas, theta=FALSE, n=NULL, n0B=FALSE, trecase=FALSE){
  classes = (unique(index))
  f1 = c(1, 2)
  if(!all(classes %in% 1:8)){
    stop("index must be one of the predefined crosses: for AB, BA, AA, BB, for female or one sex model respectively 1, 2, 3, 4, for male mice in two sex model 5, 6, 7, 8")
  }
  nclasses = length(classes)
  if((nclasses != 4)  &  (nclasses != 8)){
    stop("model requires all 4 crosses for optimization, index must have each class of 1:4 or 1:8")
  }
  ny = length(y)
  if(ny != length(index)){
    stop("index must be provided for each mouse, length(y) != length(index)")
  }
  if(ny != length(kappas)){
    stop("kappas must be provided for each mouse, length(y) != length(kappas)")
  }
  if(trecase){
    nn = length(n)
    nn0B = length(n0B)
    if(nn != nn0B){
      stop("allele specific counts mismatch: counts of allele B (n0B) don't match aggregate allele specific counts (n)")      
    }     
    f1.classes = unique(index[1:nn])
    if(!all(f1.classes %in% f1)){
      stop("data should be presented so that F1 mice with appropriate allele specific counts are listed prior to the mice without allele specific counts, such as male F1s or inbred mice")
    } 
    nya = sum(index %in% f1)
    if(nya != nn){
      stop("index doesn't correspond to the vector of allele specific reads (n), each element of n should have a corresponding class 1 or 2")
    }
    if(nclasses == 8  &  length(coef) != 9){
      stop("coef should be of length 9: add_F_trc, add_M_trc, add_F_ase, poo_F, beta_0, beta_kappa, beta_sex, beta_dom, beta_dev.dom")
    }
    if(nclasses == 4  &  length(coef) != 6){
      stop("coef should be of length 6: add_trc, add_ase, poo, beta_0, beta_kappa, beta_dom")
    }
  }else{
    if(nclasses == 8  &  length(coef) != 8){
      stop("coef should be of length 8: add_F_trc, add_M_trc, beta_0, beta_kappa, beta_sex, beta_dom, beta_dev.dom, poo")
    }
    if(nclasses == 4  &  length(coef) != 5){
      stop("coef should be of length 5: add_trc, beta_0, beta_kappa, beta_dom, beta_poo")
    }  
  }
}

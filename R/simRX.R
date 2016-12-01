`simRX` = 
function(b0f, b0m, b1f, b1m, beta_sex, beta_dom, beta_k=1, phi=1, theta=1, n=6, mean.base.cnt=50, range.base.cnt=60, perc.ase=.35, n.simu=1E4, is.X=FALSE, tauB=NULL, seed=NULL){  
  if(!is.null(seed))set.seed(seed)
  if(any(n<1))stop("all crosses must have at least one mouse")
  ns = rep(n, length=8)
  index = c(sapply(1:8, function(i){rep(i, ns[i])}))
  ind.lst = lapply(1:8, function(i){which(index %in% i)})
  #produce kappas
  min.cnt = mean.base.cnt - range.base.cnt/2
  max.cnt = mean.base.cnt + range.base.cnt/2
  kappas = log(runif(length(index), min=min.cnt, max=max.cnt))
  betas = c(beta_k, beta_sex, beta_dom)

  sex = rep(1, length(index))
  sex[index %in% 5:8] = -1
  dom = rep(0, length(index))
  dom[index %in% c(1:2, 5:6)] = 1
  X = cbind(kappas, sex, dom)
  
  if(is.X){
    F1s = index %in% c(1:2)
    if(is.null(tauB)){
      tausB = rep(.5, sum(F1s))
    }else{
      tausB = rep(tauB, sum(F1s))
    }
    l.tau.a = log(1 - tausB)
    l.tau.b = log(tausB)
    l.tau.r = l.tau.b - l.tau.a  
    sex.ase = sex[F1s]
    x = index[F1s]
    x[x == 2] = -1
    
    end = log(2) - log1p(exp(b1f))
    etas = rep(0, length(index))                                                                          #AxA, fem - base
    etas[ind.lst[[4]]] = b0f                                                                              #BxB, fem
    etas[ind.lst[[1]]] =       log1p(exp(l.tau.r[ind.lst[[1]]] + b0f + b1f)) + l.tau.a[ind.lst[[1]]] + end#AxB, fem
    etas[ind.lst[[2]]] = b1f + log1p(exp(l.tau.r[ind.lst[[2]]] + b0f - b1f)) + l.tau.a[ind.lst[[2]]] + end#BxA, fem
    etas[ind.lst[[5]]] = b1f + end                                                                        #AxB, mal
    etas[ind.lst[[6]]] = b0m + b1f + end                                                                  #BxA, mal
    etas[ind.lst[[7]]] = b1f + end                                                                        #AxA, mal
    etas[ind.lst[[8]]] = b0m + b1f + end                                                                  #BxB, mal
  }else{  
    F1s = index %in% c(1:2, 5:6)
    x = index[F1s]
    x[x == 5] = 1
    x[x == 2] = -1
    x[x == 6] = -1
    sex.ase = sex[F1s]
  
    etas = rep(0, length(index)) #base will be eta3 = 0 - AxA, fem
    ends =  - log1p(exp( - b1f))
    etas[ind.lst[[4]]] =        b0f #BxB, fem
    etas[ind.lst[[1]]] = -b1f + log1p(exp(b0f + b1f))   + ends#AxB, fem
    etas[ind.lst[[2]]] =        log1p(exp(b0f - b1f))   + ends#BxA, fem
    etas[ind.lst[[5]]] = -b1m + log1p(exp(b0m + b1m))   + ends#AxB, mal
    etas[ind.lst[[6]]] =        log1p(exp(b0m - b1m))   + ends#BxA, mal
    etas[ind.lst[[7]]] =        log1p(exp(-b1m))        + ends#AxA, mal
    etas[ind.lst[[8]]] =  b0m + log1p(exp(-b1m))        + ends#BxB, mal
  }

    mus = exp(c(X %*% betas) + etas)
    mu.ase = c(exp(b0f + b1f * x[sex.ase == 1]), exp(b0m + b1m * x[sex.ase == -1]))
    pi1    = mu.ase/(1 + mu.ase)
    aa1    = pi1/theta
    bb1    = 1/theta - aa1
      
    y = matrix(0, nrow=n.simu, ncol=length(index))
    n0 = n = matrix(0, nrow=n.simu, ncol=length(sex.ase))
         
    #ind.F1 = which(index %in% F1s)
    ind.F1 = which(F1s)
    for(k in 1:n.simu){
      if(k%%500 == 0)message("generating data for ", k, "'th simulation")
      
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      # simulate data: total read couont
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
      y[k, ] = rnegbin(length(mus), mus, 1/phi)
  
      n[k, ] = rbinom(length(ind.F1), y[k, ind.F1], perc.ase)
      
      for(i in 1:ncol(n)){
        if(n[k, i] == 0){
          n0[k, i] = 0
        }else{
          n0[k, i] = rbetabinom.ab(1, n[k, i], aa1[i], bb1[i])
        }
      }
    }
   rownames(y) = rownames(n) = rownames(n0) = paste("GID", 1:n.simu, sep="")
   crosses = c("AB_fem", "BA_fem", "AA_fem", "BB_fem", "AB_mal", "BA_mal", "AA_mal", "BB_mal")
   mouse.id = paste(crosses[index], "_ID", 1:length(index), sep="")
   colnames(y) = mouse.id
   colnames(n) = colnames(n0) = mouse.id[ind.F1]

   if(is.X)   
     out = list(index=index, y=y, n=n, n0B=n0, kappas=kappas,tausB=tausB, geneids=rownames(y))
   else
     out = list(index=index, y=y, n=n, n0B=n0, kappas=kappas, geneids=rownames(y))
   return(out)
}

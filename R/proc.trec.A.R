`proc.trec.A` = 
function(rc, hessian=FALSE){
  index = rc$index 
  y = rc$y
  kappas = rc$kappas
  geneids = rc$geneids

  if(is.null(geneids)){geneids = 1:nrow(y)}
  if(nrow(y) != length(geneids)){
    stop("geneids length is not equal to the number of rows of y")
  }
  if(ncol(y) != length(index)){
    stop("number of rows of y should match length of index")
  }  

  twosex = max(index)>4
  if(twosex){
    if(!all((1:8) %in% index)){
      stop("package expects to have both sexes of all 4 groups: AB, BA, AA, BB")      
    }
    ind.lst = lapply(1:8, function(i){which(index %in% i)})
  }else{
    if(!all((1:4) %in% index)){
      stop("package expects to have all 4 groups: AB, BA, AA, BB")
    }
    ind.lst = lapply(1:4, function(i){which(index %in% i)})
  }
  if(is.null(kappas)){
    kappas = colSums(y)
  }else{
    if(ncol(y) != length(kappas)){
      stop("number of columns of y should match length of kappas")
    }
  }
  #maxtrial = 10;
  eps.sh = 1E-2; maxiter.sh = 20
  eps = 1E-3; maxiter = 100
  iphi = 2; theta = .5

  ngenes = length(geneids)
  nsamples = length(index)
  intercept = rep(1, nsamples)
  dom = as.numeric(index %in% c(1, 2, 5, 6))
  poo = rep(0, length(index))
  poo[index %in% c(1, 5)] = 1
  poo[index %in% c(2, 6)] = -1

  if(twosex){
    sex.tot = intercept;sex.tot[index>4] = -1
    X = cbind(intercept, kappas, sex.tot, dom, dom * sex.tot, poo, poo * sex.tot)
    X.sex.dom = X[, -5]
    X.dom = X.sex.dom[, -4]
    X.sex.poo = X[, -7]
    X.poo = X.sex.poo[, -6]  
    X.sex = X[, c(-3, -5, -7)]  
  
    full.model = 2 + ncol(X)
    df=c(2, 2, 2, 4, 1, 1, 1)
    pvals = matrix(NA, nrow=ngenes, ncol=length(df))
    colnames(pvals) = c("pval_add", "pval_poo", "pval_dom", "pval_sex", "pval_sex.add", "pval_sex.poo", "pval_sex.dom")
    coef.full = matrix(NA, nrow=ngenes, ncol=full.model + 2)
    colnames(coef.full) = c("add_F_trc", "add_M_trc", "interc", "kappa", "sex", "dom", "sex:dom", "poo_F", "poo_M", "ll", "phi")
  }else{
    X = cbind(intercept, kappas, dom, poo)
    X.dom = X[, -3]
    X.poo = X[, -4]

    sex.ase = NULL
  
    full.model = 1 + ncol(X)
    df=c(1, 1, 1)
    pvals = matrix(NA, nrow=ngenes, ncol=length(df))
    colnames(pvals) = c("pval_add", "pval_poo", "pval_dom")
    coef.full = matrix(NA, nrow=ngenes, ncol=full.model + 2)
    colnames(coef.full) = c("add_trc", "interc", "kappa", "dom", "poo", "ll", "phi")
  }
  ll.ind = full.model + 1; phi.ind = ll.ind + 1
  rownames(pvals) = geneids
  rownames(coef.full) = geneids 

  coef.add =  coef.poo = coef.dom = coef.full
  coef.sex = coef.sex.add = coef.sex.poo = coef.sex.dom = coef.add
  message("processing ", ngenes, " genes")
  errorlist = NULL
  if(twosex){
    for(j in 1:ngenes){
      if(j%%100 == 0){
        message(j, "th gene")
      }
      #fit the smallest first 
      #4. sex effect: b0F' - b0M'=b0F - b0M=b1F - b1M=beta2=beta4=0
      inits.sex = c(0, -5, 1, 0, 0)
      res.sex = trec.sex.A(yi=y[j, ], ind.lst=ind.lst, X=X.sex, 
               start=inits.sex, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
  
      if(is.null(res.sex)){
        inits.poo = c(0, 0, -5, 1, 0, 0, 0)
        inits.add = c(-5, 1, 0, 0, 0, 0, 0)
        inits.dom = c(0, 0, -5, 1, 0, 0, 0)
      }else{
        inits.add = res.sex[3:9]
        inits.poo = res.sex[1:7]
        inits.dom = res.sex[c(1:5, 8:9)]
      }
      #2. POE
      res.poo = trec.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X.poo, twosex=twosex, 
               start=inits.poo, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      #1. Strain
      res.add = trec.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
              start=inits.add, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      #3. Dominance
      res.dom = trec.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, 
              start=inits.dom, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      #5. sex specific addain: b0F' - b0M'=b0F - b0M=0
      if(is.null(res.add)){inits.sex.add = c(0, -5, 1, 0, 0, 0, 0, 0)}else{inits.sex.add = res.add[2:9]}
      res.sex.add = trec.sex.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, 
               start=inits.sex.add, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      #6. sex specific poo: b1F=b1M
      if(is.null(res.poo)){inits.sex.poo = c(0, 0, -5, 1, 0, 0, 0, 0)}else{inits.sex.poo = res.poo[1:8]}
      res.sex.poo = trec.sex.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X.sex.poo, 
               start=inits.sex.poo, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      #7. sex specific Dominance
      if(is.null(res.dom)){inits.sex.dom = c(0, 0, -5, 1, 0, 0, 0, 0)}else{inits.sex.dom = res.dom[c(1:6, 8:9)]}
      res.sex.dom = trec.sex.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.sex.dom, 
             start=inits.sex.dom, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      
      #find better start for full    
      lls = list(res.add, res.poo, res.dom, res.sex, res.sex.add, res.sex.poo, res.sex.dom)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind, trecase=FALSE)                         
      #fit full
      res.full = trec.full.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             start=new.inits$inits, iphi=new.inits$iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      if(is.null(res.full)){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next          
      }
      iphi_i = res.full[phi.ind]
      #refit add
      inits.add = res.full[3:9]
      inits.poo = res.full[1:7]
      inits.dom = res.full[c(1:5, 8:9)]
      inits.sex = c(mean(res.full[1:2]), res.full[3:4], mean(res.full[6:7]), mean(res.full[8:9]))
      inits.sex.add = c(mean(res.full[1:2]), res.full[3:9])
      inits.sex.poo = c(res.full[1:7], mean(res.full[8:9]))
      inits.sex.dom = c(res.full[1:5], mean(res.full[6:7]), res.full[8:9])
      res.add1 = trec.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             start=inits.add, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit poo
      res.poo1 = trec.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X.poo, twosex=twosex, 
             start=inits.poo, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit dom
      res.dom1 = trec.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, 
             start=inits.dom, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit sex
      res.sex1 = trec.sex.A(yi=y[j, ], ind.lst=ind.lst, X=X.sex, 
             start=inits.sex, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #5. refit sex specific addain
      res.sex.add1 = trec.sex.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, 
             start=inits.sex.add, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #6. sex specific poo: b1F=b1M
      res.sex.poo1 = trec.sex.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X.sex.poo, 
             start=inits.sex.poo, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #7. sex specific Dominance
      res.sex.dom1 = trec.sex.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.sex.dom, 
             start=inits.sex.dom, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)

      res.add = get.better(lst=list(res.add, res.add1), ind=ll.ind, trecase=FALSE)
      res.poo = get.better(lst=list(res.poo, res.poo1), ind=ll.ind, trecase=FALSE)
      res.dom = get.better(lst=list(res.dom, res.dom1), ind=ll.ind, trecase=FALSE)
      res.sex = get.better(lst=list(res.sex, res.sex1), ind=ll.ind, trecase=FALSE)
      res.sex.add = get.better(lst=list(res.sex.add, res.sex.add1), ind=ll.ind, trecase=FALSE)
      res.sex.poo = get.better(lst=list(res.sex.poo, res.sex.poo1), ind=ll.ind, trecase=FALSE)
      res.sex.dom = get.better(lst=list(res.sex.dom, res.sex.dom1), ind=ll.ind, trecase=FALSE)
  
      #find better start for full
      lls = list(res.full, res.add, res.poo, res.dom, res.sex, res.sex.add, res.sex.poo, res.sex.dom)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind, trecase=FALSE)
      #fit full
      res.full1 = trec.full.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             start=new.inits$inits, iphi=new.inits$iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      res.full = get.better(lst=list(res.full, res.full1), ind=ll.ind, trecase=FALSE)
      coef.full[j, ] = res.full
      tests = test(full=res.full, short=lls[-1], df=df, ind=ll.ind, twosex=twosex, genei=j)
      pvals[j, ] = tests$pvals
      errorlist = c(errorlist, tests$error)
      if(!is.null(res.add)){coef.add[j, ] = res.add}
      if(!is.null(res.poo)){coef.poo[j, ] = res.poo}
      if(!is.null(res.dom)){coef.dom[j, ] = res.dom}
      if(!is.null(res.sex)){coef.sex[j, ] = res.sex}
      if(!is.null(res.sex.add)){coef.sex.add[j, ] = res.sex.add}
      if(!is.null(res.sex.poo)){coef.sex.poo[j, ] = res.sex.poo}
      if(!is.null(res.sex.dom)){coef.sex.dom[j, ] = res.sex.dom}
    }
    coef.full[, phi.ind] = 1/coef.full[, phi.ind]
    coef.add[, phi.ind] = 1/coef.add[, phi.ind]
    coef.poo[, phi.ind] = 1/coef.poo[, phi.ind]    
    coef.dom[, phi.ind] = 1/coef.dom[, phi.ind]    
    coef.sex[, phi.ind] = 1/coef.sex[, phi.ind]
    coef.sex.add[, phi.ind] = 1/coef.sex.add[, phi.ind]
    coef.sex.poo[, phi.ind] = 1/coef.sex.poo[, phi.ind]
    coef.sex.dom[, phi.ind] = 1/coef.sex.dom[, phi.ind]
    
    res = list(pvals=pvals, coef.full=coef.full, coef.add=coef.add, coef.poo=coef.poo, 
    coef.dom=coef.dom, coef.sex=coef.sex, coef.sex.add=coef.sex.add, coef.sex.poo=coef.sex.poo, 
    coef.sex.dom=coef.sex.dom, errorlist=errorlist)
  }else{
    for(j in 1:ngenes){
    #for(j in 1:2){
      if(j%%100 == 0){
        message(j, "th gene")
      }
      #fit the smallest first 
      #1. Strain
      inits.add = c(-5, 1, 0, 0)
      #yi=y[j, ];sex=sex.ase; ni=nis[j, ];ni0=ni0B[j, ];start=inits.add
      res.add = trec.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
              start=inits.add, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)     
      if(is.null(res.add)){
        inits.poo = c(0, -5, 1, 0)
        inits.dom = c(0, -5, 1, 0)
      }else{
        inits.poo = res.add[1:4]
        inits.dom = res.add[c(1:3, 5)]
      }
      #2. POE
      res.poo = trec.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X.poo, twosex=twosex, 
               start=inits.poo, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      #3. Dominance
      res.dom = trec.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, 
              start=inits.dom, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)        
      #find better start for full    
      lls = list(res.add, res.poo, res.dom)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind, trecase=FALSE)                         
      #fit full
      res.full = trec.full.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             start=new.inits$inits, iphi=new.inits$iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      if(is.null(res.full)){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next          
      }
      iphi_i = res.full[phi.ind]
      theta_i = res.full[ll.ind + 2]
      #refit add
      inits.add = res.full[2:5]
      inits.poo = res.full[1:4]
      inits.dom = res.full[c(1:3, 5)]
      res.add1 = trec.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             start=inits.add, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit poo
      res.poo1 = trec.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X.poo, twosex=twosex, 
             start=inits.poo, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit dom
      res.dom1 = trec.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, 
             start=inits.dom, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)

      res.add = get.better(lst=list(res.add, res.add1), ind=ll.ind, trecase=FALSE)
      res.poo = get.better(lst=list(res.poo, res.poo1), ind=ll.ind, trecase=FALSE)
      res.dom = get.better(lst=list(res.dom, res.dom1), ind=ll.ind, trecase=FALSE)
  
      #find better start for full
      lls = list(res.full, res.add, res.poo, res.dom)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind, trecase=FALSE)
      #fit full
      res.full1 = trec.full.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             start=new.inits$inits, iphi=new.inits$iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      res.full = get.better(lst=list(res.full, res.full1), ind=ll.ind, trecase=FALSE)
      coef.full[j, ] = res.full
      tests = test(full=res.full, short=lls[-1], df=df, ind=ll.ind, twosex=twosex, genei=j)
      pvals[j, ] = tests$pvals
      errorlist = c(errorlist, tests$error)
      if(!is.null(res.add)){coef.add[j, ] = res.add}
      if(!is.null(res.poo)){coef.poo[j, ] = res.poo}
      if(!is.null(res.dom)){coef.dom[j, ] = res.dom}
    }
    coef.full[, phi.ind] = 1/coef.full[, phi.ind]
    coef.add[, phi.ind] = 1/coef.add[, phi.ind]
    coef.poo[, phi.ind] = 1/coef.poo[, phi.ind]    
    coef.dom[, phi.ind] = 1/coef.dom[, phi.ind]    
    res = list(pvals=pvals, coef.full=coef.full, coef.add=coef.add, coef.poo=coef.poo, coef.dom=coef.dom, errorlist=errorlist)  
  }
  if(hessian){
    if(ngenes>100)warning("This option is recommended only for the extra analysis on a subset of genes.")    
    hess.lst = as.list(rep(NA, ngenes))
    for(j in 1:ngenes){
      tag = tryCatch({    
        nll = nLogLik.trec.A(coef=coef.full[j, 1:(ll.ind-1)], phi=coef.full[j, phi.ind], 
          rc, j, hessian=TRUE)
        hess.lst[[j]] = nll$hess
        0
      }, error=function(e) {
        warning("a problem with hessian")
        1
      })
      names(hess.lst) = geneids
    }
    res$hess.lst = hess.lst
  }  
  return(res)
}


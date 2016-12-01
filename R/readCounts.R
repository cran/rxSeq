'readCounts' = function(index, y, n=NULL, n0B=NULL, chrom="auto", kappas=NULL, tausB=NULL, 
                         genes.switch="ENSMUSG00000086503", geneids=rownames(y))
{
  tech.ctrl = list(iphi_l=0.01, iphi_u=1e5, theta_l=1e-5, theta_u=1e2, l2=log(2), maxtrial=10)
#  globalVariables("tech.ctrl")
  if(is.null(kappas)){
    kps = colSums(y)
  }else{
    kps = kappas
  }
  model = "short"
  if(!is.null(n)&!is.null(n0B)){
    model = "full"
  }
  if(ncol(y) != length(kps)){
    stop("number of columns of y should match length of kappas")
  }
        rc = list(
                index = index,
                y = y,
                n = n,
                n0B = n0B,
                kappas = kps,
                tausB = tausB,
                genes.switch=genes.switch,
                geneids=geneids,
                chrom = chrom,
                model = model,
                tech.ctrl = tech.ctrl
       )

        ## Set the name for the class
        class(rc) = append(class(rc),"ReadCounts")
        return(rc)
}
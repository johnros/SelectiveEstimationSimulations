#### Utility functions FCR of CQC after BH simulation ####

## Packages ####
library(selectiveCI)
library(lattice)
library(grid)
library(ggplot2)
# library(tractor.base)
# library(AnalyzeFMRI)



## Convert 1D variance to 1D FWHM
var2FWHM<- function(var){
  FWHM<- sqrt(8 * log(2) * var) 
  return(FWHM)
}
## Testing:
# Should return unit matrix:
# diag(var2FWHM(1)^2, 3) / (8 * log(2))



## Test if a CI covers a parameter value
coverage.test<- function(x, l, u){
  check<- NA
  if(is.numeric(x) && is.numeric(u) && is.numeric(l)){
    stopifnot(l<=u)
    check<-  x >= l && x<= u
  }
  return(check)
}
## Testing:
# coverage.test(1,2,5)
# coverage.test(x=4,l=2,u=5)
# coverage.test(x=-1,l=-2,u=5)
# coverage.test(x=3,l=-2,u=5)
# coverage.test(x=-5,l=-4,u=-7)
# coverage.test(x=-5,l=-7,u=-5.8)





## Error catching wrapper to CQC 
tryCQC_CI<- function(...){
  result<- list(lower=NA, upper=NA)
  try(result<- QuasiConventionalCI(...))
  return(result)
}
## Testing:
# tryCQC_CI(x=4, sigsq=1, lambda=2, cutoff=2, alpha=0.05)
# tryCQC_CI(x=1, sigsq=1, lambda=2, cutoff=2, alpha=0.05)





## Make 3D signal field
makeSignal3D<- function(dim.brain, n.sample, pi1, rho, ...){
  n.tests<- prod(dim.brain)
  
  scaled.signal<- sqrt(n.sample-3) * atanh(rho) 
  result<- array(rbinom(n.tests, 1, pi1) * scaled.signal, dim=dim.brain)
  
  return(result)
}  
## Testing:
# .test<- makeSignal3D(dim.brain=c(5,5,5), n.sample=16L, pi1=0.5, rho=0.5, msk=array(!FALSE, dim=c(5,5,5)))
# dim(.test)
# table(is.na(.test))
# levelplot(.test)






## Generate simulaton configuration
load(file='observed_sigma_voxels_Forma.RData')
generateConfig3D<-function(
  dim.brain, rho, voxdim=c(1,1,1), lambda=2, alpha=0.1, FDR=0.1, ksize=17, observed.sigma=observed.sigma.voxels, 
  sigma.scale, n.sample, pi1, cl.size=NA
){
  
  config<- list(    
    dim.brain=dim.brain, 
    rho=rho, 
    lambda=lambda, 
    alpha=alpha, 
    voxdim=voxdim, 
    ksize=ksize, 
    sigma=observed.sigma * sigma.scale,
    FDR=FDR,
    n.sample=n.sample, 
    pi1=pi1,
    cl.size=cl.size)
  
  return(config)
}
## Testing:
# generateConfig3D(dim.brain=c(5,5,5), rho=0.5, sigma.scale=1, n.sample=16L, pi1=0.2,cl.size=10)




## Generate field and compute FCR: 
ComputeFCRs3D<- function(
  dim.brain, voxdim, sigma, ksize, n.sample, pi1, rho, FDR, alpha, lambda, ... ){
  ### Sketch:
  # Generate signal
  # Generate noise
  # Compute B-H mask
  # Compute CQC CIs
  # Compute coverage
  
  result<- list(V=NA, R=NA, FCR=NA)
  
  signal.field<- makeSignal3D(dim.brain=dim.brain, n.sample=n.sample, pi1=pi1, rho=rho)    
  
  msk<- array(TRUE, dim=dim.brain)
  
  noise.field <- Sim.3D.GRF(d = dim.brain, 
                            voxdim = voxdim, 
                            sigma = sigma, 
                            ksize = ksize, 
                            mask = msk, 
                            type = "field")$mat
  
  observations<- signal.field + noise.field 
  
  pvals<- 2 * pnorm(abs(observations), lower.tail=FALSE)
  selection.ind<- p.adjust(pvals, method='BH') <= FDR
  n.selected<- sum(selection.ind)
  n.tested<- prod(dim.brain)
  pvalue.cutoff<- FDR / n.tested * n.selected
  selection.cutoff<- qnorm(1 - pvalue.cutoff / 2)
  
  
  if(n.selected>0) try({
    .cis<- sapply(observations[selection.ind], function(x) {
      tryCQC_CI(x=x, sigsq=1, lambda=lambda, cutoff=selection.cutoff, alpha=alpha)
    })
    
    coverage.ind<- mapply(coverage.test, 
                          x=signal.field[selection.ind],
                          l=.cis['lower',],
                          u=.cis['upper',])
    
    if(sum(!is.na(coverage.ind))<1) return(result)
    
    R<- length(coverage.ind)
    V<- sum(!na.omit(coverage.ind))
    observed.FCR<- ifelse(R>0, V/R, 0)
    result<- list(V=V, R=R, FCR=observed.FCR)
  })
  else{
    result<- list(V=0, R=0, FCR=0)
  }
  
  return(result)
}
## Testing:
# ComputeFCRs3D(
#   dim.brain=c(10,10,10), voxdim=c(1,1,1), sigma=diag(1,3), ksize=17, 
#   n.sample=16L, pi1=0.5, rho=0.5, FDR=0.1, alpha=0.1, lambda=2)
# do.call(ComputeFCRs3D, generateConfig3D(dim.brain=c(5,5,5), rho=0.5, sigma.scale=1, n.sample=16L, pi1=0.6))










## Replivate simualtion and return average FCR:
replicateFCRs3D<- function(replications, config, FUN) {
  
  FCRs<- matrix(nrow=replications, ncol=3, 
                dimnames=list(NULL, c('V','R','FCR')))
  
  for(i in 1:replications){
    FCRs[i,]<- unlist(do.call(FUN, config))
  }
  
  mean.FCR<- mean(FCRs[,'FCR'], na.rm=TRUE)
  sd.FCR<- sd(FCRs[,'FCR'], na.rm=TRUE)  
  length.FCR<- as.integer(sum(!is.na(FCRs[,'FCR'])))
  
  result<- c(
    mean.FCR= mean.FCR, 
    sd.FCR= sd.FCR, 
    length.FCR=length.FCR  )
  
  return(result)
}
## Testing:
#  .config<- generateConfig3D(dim.brain=c(10,10,10), rho=0.5, sigma.scale=1, n.sample=16L, pi1=0.6, cl.size=10)
#  replicateFCRs3D(replications=4L, config=.config, FUN=ComputeFCRs3D)







## Wrap replication of FCR to fit apply structure
wrapComputeFCRs3D<- function(x){
  config<- generateConfig3D(
    dim.brain=c(10,10,10), 
    rho=x[['rho']], 
    sigma.scale=x[['sigma.scale']], 
    n.sample=x[['n.sample']], 
    pi1=x[['pi1']])
  
  replicateFCRs3D(replications=x[['replications']], config=config, FUN=ComputeFCRs3D)
}



test.configurations<- expand.grid(
  rho= c(0, 0.2, 0.4),
  replications=2L, 
  sigma.scale= c(0.2, 0.7),
  pi1=c(0, 0.1, 0.5),
  n.sample=c(4, 16, 32),
  cl.size=c(1,2))
## Testing:
# wrapComputeFCRs3D(test.configurations[20,])
# wrapComputeFCRs3D_threshed(test.configurations[20,])
# apply(test.configurations, 1, wrapComputeFCRs3D)










## Packages ####
library(tractor.base)
library(selectiveCI)
library(grid)
library(ggplot2)
library(AnalyzeFMRI)






#### Return an index for each contigous excursion set
giveExcursionIndexes <- function(mask) {
  dims <- dim(mask)
  assignemt.matrix<- array(0L, dim = dims)
  assignment<- 1
  
  for(i in 1:dims[1]){
    for (j in 1:dims[2]){
      for(k in 1:dims[3]){
        # try(rm(i,j,k)); i<- 20; j<- 20; k<- 26
        # Check if already assigned:  			
        # If active, look for assignment in neighbouring voxels:
        if(!is.na(mask[i,j,k]) && mask[i,j,k]!=0) {
          # Assign new group number 
          assignemt.matrix[i,j,k]<- assignment 
          assignment<- assignment + 1 
          update.assignemnt<- FALSE
          
          # Check if neighbouring group already exists:										
          for( ii in (i-1):(i+1) ){
            for ( jj in (j-1):(j+1) ){
              for( kk in (k-1):(k+1) ){
                # try(rm(ii,jj,kk)); ii<- i-1; jj<- j-1; kk<- k-1
                if(ii==i && jj==j && kk==k) next() 
                if(ii>dims[1] || jj>dims[2] || kk>dims[3]) next()																
                if(ii<1 || jj<1 || kk<1) next()
                if(is.na(mask[ii,jj,kk])) next()								
                
                if(assignemt.matrix[ii,jj,kk]!=0L) {
                  assignemt.matrix[i,j,k]<- assignemt.matrix[ii,jj,kk]
                  update.assignemnt<- TRUE
                }
              }
            }
          }
          if(update.assignemnt) assignment<- assignment - 1					
        }									
      }				
    }
  }
  return(assignemt.matrix)
}
## Testing:
# A<- array(0, dim=c(4,4,4))
# A[c(1,2,4),c(1,3,4),c(1,3,4)]<- 1
# .test<- giveExcursionIndexes(A)
# levelplot(A)
# levelplot(.test)
# A<- array(1, dim=c(4,4,4))
# A<- array(TRUE, dim=c(4,4,4))
# .test<- giveExcursionIndexes(A)
# table(.test)
# A<- array(FALSE, dim=c(5,5,5))
# A[c(1,2,4),c(1,4),c(1,3,4)]<- TRUE
# A[c(2),c(2),c(2)]<- TRUE
# A
# .test<- giveExcursionIndexes(A)
# .test[,,]
# table(.test)
# levelplot(A)
# levelplot(.test)




## Apply cluster size threshing to given mask
getSizeThreshMask<- function(selection.mask, size, brain.mask){
  effective.mask<- selection.mask & brain.mask
  full.data.excursions<- giveExcursionIndexes(effective.mask)
  indexes<- table(full.data.excursions)
  if (length(indexes)==1L && names(indexes)=="1") large.excursion.indexes<- table(full.data.excursions)>size
  else large.excursion.indexes<- table(full.data.excursions)[-1]>=size
  full.data.size.mask<- apply(full.data.excursions, c(1,2,3), function(x) x %in% names(large.excursion.indexes)[large.excursion.indexes])
  return(full.data.size.mask)
}
## Testing:
# TODO: fix size masking
# debug(getSizeThreshMask)
# getSizeThreshMask(array(FALSE, dim=c(5,5,5)), 2, array(TRUE, dim=c(5,5,5)))
# getSizeThreshMask(array(TRUE, dim=c(5,5,5)), 1, array(TRUE, dim=c(5,5,5)))
# A<- array(FALSE, dim=c(5,5,5))
# A[c(1,2,4),c(1,4),c(1,3,4)]<- TRUE
# A[c(2),c(2),c(2)]<- TRUE
# levelplot(A)
# levelplot(getSizeThreshMask(A, 0, array(TRUE, dim=c(5,5,5))))
# levelplot(getSizeThreshMask(A, 1, array(TRUE, dim=c(5,5,5))))
# levelplot(getSizeThreshMask(A, 3, array(TRUE, dim=c(5,5,5))))
# levelplot(getSizeThreshMask(A, 4, array(TRUE, dim=c(5,5,5))))
# levelplot(getSizeThreshMask(A, 5, array(TRUE, dim=c(5,5,5))))










getPearsonBonferonniCutoff<- function(n, alpha, tests){
  # Note: alpha is expected to be *uncorrected*
  # n<- 100
  # alpha <- 0.05
  # tests<- 100
  pvale.cutoff<- alpha/tests
  zscore.cutoff<- qnorm(pvale.cutoff / 2, lower.tail=FALSE)
  pearson.cutoff<- zscore2pearson(z=zscore.cutoff, n=n)
  return(pearson.cutoff)
}
## Testing:
# getPearsonBonferonniCutoff(n=100, alpha=0.05, tests=100)
# getPearsonBonferonniCutoff(n=100, alpha=0.005, tests=100)
# getPearsonBonferonniCutoff(n=100, alpha=0.05, tests=1000)
# getPearsonBonferonniCutoff(n=10, alpha=0.05, tests=1000)


## Convert 1D variance to 1D FWHM
var2FWHM<- function(var){
  FWHM<- sqrt(8 * log(2) * var) 
  return(FWHM)
}
## Testing:
# Should return unit matrix:
# diag(var2FWHM(1)^2, 3) / (8 * log(2))



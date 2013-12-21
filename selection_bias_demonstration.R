#### Simulate distribution of correlation after selection ####



getRejectedDistribution <- function (resels, 
                                     n.subjects, 
                                     rho, 
                                     error.rate, 
                                     error.type='BH') {
  rho.sd<- 1 / sqrt(n.subjects-3)
  rho.z<- atanh(rho)
  correlations<- tanh(rnorm(resels, mean=rho.z, sd=rho.sd ))
  pvals<- 2 * pnorm(abs(atanh(correlations)), sd=rho.sd, lower.tail=FALSE)
  
  adjusted.pvals<- p.adjust(pvals, method=error.type)
  rejection.indicator<- adjusted.pvals < error.rate
  result<- correlations[rejection.indicator]
  return(result)  
}
## Testing:
# .testing<- getRejectedDistribution(resels=111, n.subjects=16, rho=0.2, error.rate=0.5)
# mean(.testing)


# Minimal number of resels was computed from Tom et al's data
# Maximal number of resels equal number of voxels
configurations<- expand.grid(
  resels=round(c(5e3, 1e4, 3e5),-1),              
  n.subjects=floor(seq(12, 100, length.out=6)),
  rho=seq(0,1, length.out=10),
  error.rate=c(0.05, 0.2), 
  error.type='bonf', stringsAsFactors=FALSE)



## Run simulation
selected.correlations<- list()
for(i in 1: nrow(configurations)){
  selected.correlations[[i]]<- do.call(getRejectedDistribution, as.list(configurations[i,]))
}



## Compute means:
selection.bias<- cbind(configurations, 
                       mean=sapply(selected.correlations, mean, na.rm=TRUE),
                       median=sapply(selected.correlations, median, na.rm=TRUE),
                       sd=sapply(selected.correlations, sd, na.rm=TRUE), 
                       n=sapply(selected.correlations, function(x) sum(!is.na(x))))
selection.bias$arm<- with(selection.bias, sd/sqrt(n)*2)


# To report SEs in paper:
max(selection.bias$arm, na.rm=TRUE)


## Vizualize results:
library(ggplot2)
error.rates<- unique(selection.bias$error.rate)
.subset<- subset(selection.bias, error.rate==error.rates[1])

base.plot<- ggplot(data=.subset) + xlab("True Correlation") + ylab("Mean of Selected Correlations")+
  geom_point(aes(y=mean , x=rho)) + 
  facet_grid(resels~ n.subjects , labeller = label_both) +
  geom_abline(intercept = 0, slope = 1, lty=1, col='grey')
base.plot







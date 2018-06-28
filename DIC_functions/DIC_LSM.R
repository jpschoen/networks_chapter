# clear workspace
rm(list=ls())
set.seed(19)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#libraries
library(latentnet)

# Load the florentine marriage network
data(florentine)

# Posterior parameter estimation:
p.flo <- ergmm(flomarriage ~ edges + kstar(2),
               control=ergmm.control(burnin=10000,
                                     sample.size= 10000,
                                     interval=5))
# MCMC diagnostics and posterior summaries:
mcmc.diagnostics(p.flo)
summary(p.flo)


#calculate DIC from ergmm object
ergmm_DIC = function(ergmm_object){
  L = ergmm_object$mcmc.mle$lpY
  llSum = 0
  n = length(ergmm_object$sample)
  for (s in 1:n){
    llSum = llSum + ergmm_object$sample$lpY[s]
  }
  P = 2 * (L - (llSum/n))
  DIC = -2*(L-P)
  DIC
}

ergmm_DIC(p.flo)



#bergm results for comparison 
#AIC: 112.1    
#BIC: 117.7   
#DIC -39.08868
#BERGM Betas
#edges  -1.67294  
#kstar2  0.01265

# discussion about problems with DIC
#http://rsos.royalsocietypublishing.org/content/5/3/171519
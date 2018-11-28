# functions for calculating DIC


#calculate DIC from LSM object
LSM_DIC = function(ergmm_object){
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

# calculate DIC for bergm
# function to calculate likelihood at each posterior bergm sample
ll.bergm <- function(bergm.posterior.sample,bergm.offset.formula,MCMC.samplesize=2000,nsteps=50){
  # bergm.posterior is a matrix with posterior draws in rows and parameters in columns
  # bergm.offset.formula is the formula on which bergm.posterior was inferred, but with
  ## each term enclosed in offset()
  # MCMC.samplesize sets the sample size used to estimate the likelihood
  # nsteps sets the number of bridges used to estimate the likelihood
  require(ergm)
  thetas <- bergm.posterior.sample
  lls <- numeric(nrow(thetas))
  for(i in 1:length(lls)){
    estML <- ergm(bergm.offset.formula,offset.coef=thetas[i,])
    lli <- logLik(estML,add=T,control=control.logLik.ergm(MCMC.samplesize=MCMC.samplesize,nsteps=nsteps))
    lls[i] <- as.numeric(lli$mle.lik)
  }
  thetabar <- apply(thetas,2,mean)
  estML <- ergm(bergm.offset.formula,offset.coef=thetabar)
  llbar <- logLik(estML,add=T,control=control.logLik.ergm(MCMC.samplesize=MCMC.samplesize,nsteps=nsteps))
  list(lls=lls,llbar=llbar)
}

# function to calculate DIC for bergm
dic.bergm <- function(ll.bergm.result){
  Dtbar <- -2*as.numeric(ll.bergm.result$llbar$mle.lik)
  Dbar <- -2*mean(ll.bergm.result$lls)
  2*Dbar-Dtbar
}

# clear workspace
rm(list=ls())
set.seed(19)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#libraries
library(Bergm)
library(ergm)

# Load the florentine marriage network
data(florentine)
# Posterior parameter estimation:
p.flo <- bergm(flomarriage ~ edges + kstar(2),
               burn.in = 50,
               aux.iters = 500,
               main.iters = 500,
               gamma = 1)
# MCMC diagnostics and posterior summaries:
bergm.output(p.flo)

# estimate the same model using ML
estML0 <- ergm(flomarriage ~ edges + kstar(2))

# evaluate the log-likelihood at the first set of parameter
# values sampled from bergm
estML <- ergm(flomarriage ~ offset(edges) + offset(kstar(2)),offset.coef=p.flo$Theta[1,])
system.time(ll <- logLik(estML,add=T,control=control.logLik.ergm(MCMC.samplesize=100000,nsteps=50)))
ll$mle.lik

# calculate DIC for bergm
# function to calculate likelihood at each posterior bergm sample
ll.bergm <- function(bergm.posterior.sample,bergm.offset.formula,MCMC.samplesize=100000,nsteps=50){
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

# simple example with values set too low
lls <- ll.bergm(p.flo$Theta,flomarriage ~ offset(edges) + offset(kstar(2)),MCMC.samplesize=100,nsteps=10)

# calculate the dic
dic.bergm(lls)

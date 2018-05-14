# clear workspace
rm(list=ls())
set.seed(19)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#libraries

# load data
load("ICPSR_21600_W1_1/DS0001/21600-0001-Data.rda") #questionaire
load("ICPSR_21600_W1_2/DS0002/21600-0002-Data.rda") #context
load("ICPSR_21600_W1_3/DS0003/21600-0003-Data.rda") #network variables
load("ICPSR_21600_W1_4/DS0004/21600-0004-Data.rda") #weights

# 4,431 should be final sample size

# subset and create network

# Wang et al code begins


# R script for ERGM estimation
model1 <- ergm(net~edges,constraints=~bd(attribs=sexattr,maxout=maxout))
model2 <- ergm(net~edges+mutual,constraints=~bd(attribs=sexattr,maxout=maxout))
model3 <- ergm(net~edges+mutual +absdiffcat('grade')+nodemix('female',base=1), constraints=~bd(attribs=sexattr,maxout=maxout))
model4 <- ergm(net~edges+mutual +absdiffcat('grade')+nodemix('female',base=1) +nodematch(class)+nodematch(clubs)+nodematch(sports), constraints=~bd(attribs=sexattr,maxout=maxout))
model5 <- ergm(net~edges+mutual+absdiffcat('grade')+nodemix('female',base=1) +nodematch(class)+nodematch(clubs)+nodematch(sports) +gwesp(δopt,fixed=T), constraints=~bd(attribs=sexattr,maxout=maxout))
# Application of convergence criterions (Model 5 as an example)
sf <- model5$sample-model5$sample.obs
t.k <- abs(apply(sf,2,mean))/apply(model5$sample,2,sd)
tconv.max <- sqrt(t(apply(sf,2,mean) %*% solve(as.matrix(cov(sf))) %*% apply(sf,2,mean)))
while (max(t.k)>0.1 | tconv.max>0.25) {
  par <- coef(model5)
  model5 <- ergm(net~edges+mutual+absdiffcat('grade')+nodemix('female',base=1)
                 +nodematch(class)+nodematch(clubs)+nodematch(sports)
                 +gwesp(δopt,fixed=T),
                 constraints=~bd(attribs=sexattr,maxout=maxout),
                 control=control.ergm(init=par,MCMLE.maxit=20000))
  sf <- model5$sample-model5$sample.obs
  t.k <- abs(apply(sf,2,mean)/apply(model5$sample,2,sd))
  tconv.max <- sqrt(t(apply(sf,2,mean) %*% solve(as.matrix(cov(sf))) %*% apply(sf,2,mean)))
}
#R script for ERGM simulation (Model 5 as an example)
net.fit<-model5
net.sim5<-simulate(net.fit,
                   constraints=~observed+bd(attribs=sexattr,minout=minout,maxout=maxout),
                   nsim=250)
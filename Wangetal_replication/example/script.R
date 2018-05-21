# clear workspace
rm(list=ls())
set.seed(19)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#We will need the following libraries:
library(network)
library(ergm)
library(gtools)

#This undirected adjacency matrix includes 'all' relationships between characters in Grey's Anatomy. 
#Read in the adj matrix without the first column
grey <- read.csv("Grey_sex.csv",header=T)[,-1]

# create some missingness
rc <- seq(3,44,3)
rc <- permutations(n=length(rc),r=2,v=rc,repeats.allowed=F)
for(i in 1:nrow(rc)){
  grey[rc[i,1],rc[i,2]] <-NA
}


#Read in node attributes, which include Gender, Race, Birth Year, Professional Position, Season he/she first appears, and astrological sign.
attributes <- read.csv("Grey_Attributes.csv",header=T)


#Create a network object using the sociomatrix and its attributes
ga.net<-network(grey, vertex.attrnames=colnames(attributes),
                directed=F, hyper=F, loops=F, multiple=F, bipartite=F)
ga.net[rbind(c(1,4),c(3,5))] <- NA

#Set vertex attributes
set.vertex.attribute(ga.net,names(attributes),attributes)


#Model 5
model5<-ergm(ga.net~edges+nodematch("sex")+degree(1)+nodematch("race")+
                absdiff("birthyear"),control=control.ergm(
                  MCMC.burnin=50000, MCMC.interval=5000))
summary(model5)

#model5 <- ergm(net~edges+mutual+absdiffcat('grade')+nodemix('female',base=1) +nodematch(class)+nodematch(clubs)+nodematch(sports) +gwesp(δopt,fixed=T), constraints=~bd(attribs=sexattr,maxout=maxout))

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

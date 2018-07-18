#Task: compare the DIC for LSM, SBM, and ERGM models


# clear workspace
rm(list=ls())
set.seed(19)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#libraries
library(latentnet)
library(igraph)
library(CIDnetworks)
library(Bergm)
library(ergm)
library(igraphdata)
library(intergraph)

#load data
data(florentine)
data(karate)
zach <- asNetwork(karate)
data('sampson')
rm(karate)
#load functions
source("DIC_functions.R")

# Set Model Parameters
burnin = 500
sample_size = 2000
interval = 5

#create dataframe to store values
results <- as.data.frame(matrix(nrow=3, ncol=3))
row.names(results) <- c("LSM", "SBM", "ERGM")
colnames(results) <- c("florentine", "zach", "samplike")

#pull out adjacency matrix for SBMs
flo_adj <- as.matrix.network(flomarriage,matrix.type="adjacency")
zach_adj <- as.matrix.network(zach,matrix.type="adjacency")
samplike_adj <- as.matrix.network(samplike,matrix.type="adjacency")

######### Florentine data ############

#LSM
flo_LSM <- ergmm(flomarriage ~ edges + kstar(2),
               control=ergmm.control(burnin=burnin,
                                     sample.size= sample_size,
                                     interval=interval))
results$florentine[1] <-LSM_DIC(flo_LSM)

#SBM
flo_SBM<- CID.Gibbs(flo_adj, components = c(SBM(2)), draws = sample_size, burnin = burnin,
                    thin = interval)
results$florentine[2] <-flo_SBM$DIC[1]

#ERGM
flo_ERGM <- bergm(flomarriage ~ edges + kstar(2),
               burn.in = burnin,
               aux.iters = sample_size,
               main.iters = sample_size,
               gamma = 1)
flo_lls <- ll.bergm(flo_ERGM$Theta,flomarriage ~ offset(edges) + offset(kstar(2)),
                MCMC.samplesize=100,nsteps=10)
results$florentine[3] <-dic.bergm(flo_lls)



######### zach data  ############

#LSM
zach_LSM <- ergmm(zach ~ edges + kstar(2),
                 control=ergmm.control(burnin=burnin,
                                       sample.size= sample_size,
                                       interval=interval))
results$zach[1] <-LSM_DIC(zach_LSM)

#SBM
zach_SBM<- CID.Gibbs(flo_adj, components = c(SBM(2)), draws = sample_size, burnin = burnin,
                    thin = interval)
results$zach[2] <-zach_SBM$DIC[1]

#ERGM
zach_ERGM <- bergm(zach ~ edges + kstar(2),
                  burn.in = burnin,
                  aux.iters = sample_size,
                  main.iters = sample_size,
                  gamma = 1)
zach_lls <- ll.bergm(zach_ERGM$Theta,zach ~ offset(edges) + offset(kstar(2)),
                MCMC.samplesize=100,nsteps=10)
results$zach[3] <-dic.bergm(zach_lls)
######### Sampson data ############

#LSM
samplike_LSM <- ergmm(samplike~edges + mutual,
                  control=ergmm.control(burnin=burnin,
                                        sample.size= sample_size,
                                        interval=interval))
results$samplike[1] <-LSM_DIC(samplike_LSM)


#SBM
samplike_SBM<- CID.Gibbs(samplike_adj, components = c(SBM(2)), draws = sample_size, burnin = burnin,
                     thin = interval)
results$samplike[2] <-samplike_SBM$DIC[1]


#ERGM
samplike_ERGM <- bergm(samplike~edges+mutual,
                   burn.in = burnin,
                   aux.iters = sample_size,
                   main.iters = sample_size,
                   gamma = 1)
samplike_lls <- ll.bergm(samplike_ERGM$Theta,samplike ~ offset(edges) + offset(mutual),
                     MCMC.samplesize=100,nsteps=10)
results$samplike[3] <-dic.bergm(samplike_lls)


save(results, file="comparison.Rdata")
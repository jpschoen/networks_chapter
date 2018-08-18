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
sim_n = 100

#create dataframe to store values
results <- as.data.frame(matrix(nrow=3, ncol=3))
row.names(results) <- c("LSM", "SBM", "ERGM")
colnames(results) <- c("florentine", "zach", "samplike")

#pull out adjacency matrix for SBMs
flo_adj <- as.matrix.network(flomarriage,matrix.type="adjacency")
zach_adj <- as.matrix.network(zach,matrix.type="adjacency")
samplike_adj <- as.matrix.network(samplike,matrix.type="adjacency")



#1.   Let SBM(D) be the SBM estimated on D, ERGM(D) be the ERGM estimated on D, 
#         and LSM(D) be the LSM estimated on D. D = {flo,zach,sampson}


start_time <- Sys.time()
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
                MCMC.samplesize=200,nsteps=10)
results$florentine[3] <-dic.bergm(flo_lls)

save(list=c("flo_ERGM", "flo_LSM", "flo_SBM", "flo_adj","flomarriage"), file="flo_results.Rdata")

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
end_time <- Sys.time()
end_time-start_time # 1.95 hours

save(results, file="comparison.Rdata")






#Comapre DIC

#2a.  Simulate M (e.g., 500) networks from SBM(D). For each simulated network, estimate LSM, SBM, and ERGM. See how frequently the DIC favors SBM.

#Simulate flo
CID_sim <- CID(n.nodes=nrow(flo_adj),
            intercept=summary(flo_SBM)$intercept[3],
            components=list(SBM(2, block.matrix=summary(flo_SBM)$SBMcid$block.matrix, 
                      membership=summary(flo_SBM)$SBMcid$membership)),
            #class.outcome="binary", #also gaussian or ordinal
            is.directed=TRUE,
            generate=TRUE)
#create network objects for model
net_sim <- asNetwork(as.data.frame(cbind(CID_sim$edge.list,CID_sim$outcome)),directed=FALSE)  # create network objects for ergmm and ergm
adj_sim <- as.matrix.network(net_sim,matrix.type="adjacency")

#estimate models for flo simulation
simtest_SBM <- CID.Gibbs(adj_sim, components = c(SBM(2)), draws = sample_size, burnin = burnin,
thin = interval)
simtest_LSM <- ergmm(net_sim ~ edges + kstar(2),
                    control=ergmm.control(burnin=burnin,
                                          sample.size= sample_size,
                                          interval=interval))
simtest_ERGM <- bergm(net_sim ~ edges + kstar(2),
                  burn.in = burnin,
                  aux.iters = sample_size,
                  main.iters = sample_size,
                  gamma = 1)
simtest_lls <- ll.bergm(simtest_ERGM,flomarriage ~ offset(edges) + offset(kstar(2)),
                    MCMC.samplesize=100,nsteps=10)

#record info

#Simulate  zach
CID_sim <- CID(n.nodes=nrow(zach_adj),
             intercept=summary(zach_SBM)$intercept[3],
             components=list(SBM(2, block.matrix=summary(zach_SBM)$SBMcid$block.matrix, 
                                 membership=summary(flo_SBM)$SBMcid$membership)),
             #class.outcome="binary", #also gaussian or ordinal
             is.directed=TRUE,
             generate=TRUE)
#create network objects for model
net_sim <- asNetwork(as.data.frame(cbind(CID_sim$edge.list,CID_sim$outcome)),directed=FALSE)  # create network objects for ergmm and ergm
adj_sim <- as.matrix.network(net_sim,matrix.type="adjacency")

#estimate models for flo simulation
simtest_SBM <- CID.Gibbs(adj_sim, components = c(SBM(2)), draws = sample_size, burnin = burnin,
                         thin = interval)
simtest_LSM <- ergmm(net_sim ~ edges + kstar(2),
                     control=ergmm.control(burnin=burnin,
                                           sample.size= sample_size,
                                           interval=interval))
simtest_ERGM <- bergm(net_sim ~ edges + kstar(2),
                      burn.in = burnin,
                      aux.iters = sample_size,
                      main.iters = sample_size,
                      gamma = 1)
simtest_lls <- ll.bergm(simtest_ERGM,flomarriage ~ offset(edges) + offset(kstar(2)),
                        MCMC.samplesize=100,nsteps=10)


# Simulate sampson
CID_sim <- CID(n.nodes=nrow(samplike_adj),
             intercept=summary(flo_SBM)$intercept[3],
             components=list(SBM(2, block.matrix=summary(flo_SBM)$SBMcid$block.matrix, 
                                 membership=summary(flo_SBM)$SBMcid$membership)),
             #class.outcome="binary", #also gaussian or ordinal
             is.directed=TRUE,
             generate=TRUE)
#create network objects for model
net_sim <- asNetwork(as.data.frame(cbind(CID_sim$edge.list,CID_sim$outcome)),directed=FALSE)  # create network objects for ergmm and ergm
adj_sim <- as.matrix.network(net_sim,matrix.type="adjacency")

#estimate models for flo simulation
simtest_SBM <- CID.Gibbs(adj_sim, components = c(SBM(2)), draws = sample_size, burnin = burnin,
                         thin = interval)
simtest_LSM <- ergmm(net_sim ~ edges + kstar(2),
                     control=ergmm.control(burnin=burnin,
                                           sample.size= sample_size,
                                           interval=interval))
simtest_ERGM <- bergm(net_sim ~ edges + kstar(2),
                      burn.in = burnin,
                      aux.iters = sample_size,
                      main.iters = sample_size,
                      gamma = 1)
simtest_lls <- ll.bergm(simtest_ERGM,flomarriage ~ offset(edges) + offset(kstar(2)),
                        MCMC.samplesize=100,nsteps=10)




#2b.  Simulate M (e.g., 500) networks from LSM(D). For each simulated network, 
#         estimate LSM, SBM, and ERGM. See how frequently the DIC favors LSM.




#2c.  Simulate M (e.g., 500) networks from ERGM(D). For each simulated network, 
#         estimate LSM, SBM, and ERGM. See how frequently the DIC favors ERGM.



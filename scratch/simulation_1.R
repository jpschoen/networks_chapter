# clear workspace
rm(list=ls())
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#libraries
library(latentnet)
library(igraph)
library(CIDnetworks)
library(Bergm)
library(ergm)
library(igraphdata)
library(intergraph)

# set parallelization for functions
library(foreach)
library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores-2)
registerDoParallel(cl)


#load data and initial models
load("flo_results.Rdata")
#load functions
source("DIC_functions.R")

#number for seed and results
a = 1
#set dataframe

# Set Model Parameters
burnin = 500
sample_size = 20
interval = 5
bergm_MCMC_n = 10
sim_n = 100
seed = a

#set seed
set.seed(seed)
# result names
result_names <- c("SBM_flo_SBM", "SBM_flo_LSM", "SBM_flo_ERGM",
                       "LSM_flo_SBM", "LSM_flo_LSM", "LSM_flo_ERGM",
                       "SBM_flo_SBM", "SBM_flo_LSM", "SBM_flo_ERGM")
dfr <-as.data.frame(matrix(nrow=1, ncol=9))
for(n in seq(10,20,100)){
s <- n-9
f <- n  
# run parallel loop for simulations
df <-foreach(i=s:f, .combine=rbind, .packages= c("CIDnetworks","latentnet", "Bergm", "ergm", "intergraph")) %dopar% {
#Comapre DIC
#create dataframe to store values
results_row <- as.data.frame(matrix(nrow=1, ncol=9))
#2a.  Simulate M (e.g., 500) networks from SBM(D). For each simulated network, estimate LSM, SBM, and ERGM. See how frequently the DIC favors SBM.
#Simulate flo
sim <- CID(n.nodes=nrow(flo_adj),
               intercept=summary(flo_SBM)$intercept[3],
               components=list(SBM(2, block.matrix=summary(flo_SBM)$SBMcid$block.matrix, 
                                   membership=summary(flo_SBM)$SBMcid$membership)),
               class.outcome="binary", #binary, gaussian, or ordinal
               is.directed=F,
               generate=TRUE)
#create network objects for model
el <- cbind(sim$edge.list,sim$outcome) # edgelist (a=origin, b=destination, c=weight)
mat<-matrix(0, nrow(flo_adj), nrow(flo_adj))
for(i in 1:NROW(el)) mat[ el[i,1], el[i,2] ] <- el[i,3]  # SEE UPDATE
adj_sim <- mat
net_sim <- as.network(x = adj_sim , # the network object
                      directed = F, # specify whether the network is directed
                      loops = FALSE, # do we allow self ties (should not allow them)
                      matrix.type = "adjacency" # the type of input
)  # create network objects for ergmm and ergm


#estimate models for flo simulation
simtest_SBM <- CID.Gibbs(adj_sim, components = c(SBM(2)), draws = sample_size, burnin = burnin,
                         thin = interval)
results_row[1] <- simtest_SBM$DIC[1]
simtest_LSM <- ergmm(net_sim ~ edges + kstar(2),
                     control=ergmm.control(burnin=burnin,
                                           sample.size= sample_size,
                                           interval=interval))
results_row[2] <-LSM_DIC(simtest_LSM)
simtest_ERGM <- bergm(net_sim ~ edges + kstar(2),
                      burn.in = burnin,
                      aux.iters = sample_size,
                      main.iters = sample_size,
                      gamma = 1)
simtest_lls <- ll.bergm(simtest_ERGM$Theta,net_sim ~ offset(edges) + offset(kstar(2)),
                        MCMC.samplesize=bergm_MCMC_n,nsteps=10)
results_row[3] <-dic.bergm(simtest_lls)



#2b.  Simulate networks from LSM(D). For each simulated network, 
#     estimate LSM, SBM, and ERGM. See how frequently the DIC favors LSM.
#Simulate flo
sim <- simulate(flo_LSM) 
#create network objects for model
#net_sim <- asNetwork(as.data.frame(cbind(sim$edge.list,sim$outcome)),directed=FALSE)  # create network objects for ergmm and ergm
adj_sim <- as.matrix.network(sim,matrix.type="adjacency")

#estimate models for flo simulation
simtest_SBM <- CID.Gibbs(adj_sim, components = c(SBM(2)), draws = sample_size, burnin = burnin,
                         thin = interval)
results_row[4] <- simtest_SBM$DIC[1]
simtest_LSM <- ergmm(sim ~ edges + kstar(2),
                     control=ergmm.control(burnin=burnin,
                                           sample.size= sample_size,
                                           interval=interval))
results_row[5] <-LSM_DIC(simtest_LSM)
simtest_ERGM <- bergm(sim ~ edges + kstar(2),
                      burn.in = burnin,
                      aux.iters = sample_size,
                      main.iters = sample_size,
                      gamma = 1)
simtest_lls <- ll.bergm(simtest_ERGM$Theta,sim ~ offset(edges) + offset(kstar(2)),
                        MCMC.samplesize=bergm_MCMC_n,nsteps=10)
results_row[6] <-dic.bergm(simtest_lls)



#2c.  Simulate networks from ERGM(D). For each simulated network, 
#     estimate LSM, SBM, and ERGM. See how frequently the DIC favors ERGM.
#Simulate flo
ergm_obj <- ergm(flomarriage ~ edges + kstar(2),coef=bergm.output(flo_ERGM)$statistics[,1]) #create ergm object to simulate from
sim <- simulate(ergm_obj)
#create network objects for model
#net_sim <- asNetwork(as.data.frame(cbind(sim$edge.list,sim$outcome)),directed=FALSE)  # create network objects for ergmm and ergm
adj_sim <- as.matrix.network(sim,matrix.type="adjacency")

#estimate models for flo simulation
simtest_SBM <- CID.Gibbs(adj_sim, components = c(SBM(2)), draws = sample_size, burnin = burnin,
                         thin = interval)
results_row[7] <- simtest_SBM$DIC[1]
simtest_LSM <- ergmm(sim ~ edges + kstar(2),
                     control=ergmm.control(burnin=burnin,
                                           sample.size= sample_size,
                                           interval=interval))
results_row[8] <-LSM_DIC(simtest_LSM)
simtest_ERGM <- bergm(sim ~ edges + kstar(2),
                      burn.in = burnin,
                      aux.iters = sample_size,
                      main.iters = sample_size,
                      gamma = 1)
simtest_lls <- ll.bergm(simtest_ERGM$Theta,sim ~ offset(edges) + offset(kstar(2)),
                        MCMC.samplesize=bergm_MCMC_n,nsteps=10)
results_row[9] <-dic.bergm(simtest_lls)
#return results
results_row
}
dfr <- rbind(dfr,df)
save(dfr, file=paste0("sim",a,"_results.Rdata"))
}

#stop cluster
stopCluster(cl)

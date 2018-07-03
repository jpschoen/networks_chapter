# clear workspace
rm(list=ls())
set.seed(19)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#libraries
library(latentnet)
library(igraph)
library(CIDnetworks)

# Load the florentine marriage network
data(florentine)

marriage_m <- as.matrix.network(flomarriage,matrix.type="adjacency")
business_m <- as.matrix.network(flobusiness,matrix.type="adjacency")


model1 <- CID.Gibbs(marriage_m, draws = 10000, burnin = 10000,
              thin = 5)

plot(model1)
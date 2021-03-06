# clear workspace
#rm(list=ls())
#set.seed(19)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#We will need the following libraries:
library(network)
library(ergm)
library(gtools)
library(latentnet)
source('functions.R')

#This undirected adjacency matrix includes 'all' relationships between characters in Grey's Anatomy. 
#Read in the adj matrix without the first column
grey <- read.csv("example/Grey_sex.csv",header=T)[,-1]

#create data with NAs
df_NA <- insert_NAs(data=grey, type = "partial", partial_level = .1, complete_level = .1)

#Read in node attributes, which include Gender, Race, Birth Year, Professional Position, Season he/she first appears, and astrological sign.
attributes <- read.csv("example/Grey_Attributes.csv",header=T)

#Create a network object using the sociomatrix and its attributes
ga.net<-network(df_NA, vertex.attrnames=colnames(attributes),
                directed=F, hyper=F, loops=F, multiple=F, bipartite=F)

#Set vertex attributes
set.vertex.attribute(ga.net,names(attributes),attributes)

#define model
model1 <-as.formula(ga.net~edges+nodematch("sex")+degree(1)+nodematch("race")+
                           absdiff("birthyear"))

#Estimate
net_sims <- ergm_predictions(model1)

#check accuracy
score_list <- return_accuracy(net_sims, grey)
mean(score_list)




# latent space model ####

#estimate
latent_model <- ergmm(model1,control=ergmm.control(burnin=10000,sample.size= 2000,interval=5))
#simulate
latent.sims <- simulate(latent_model, nsim=2)
#check accuracy
score_list_latent <- return_accuracy_l(latent.sims[[2]], grey, df_NA)
mean(score_list_latent)





# clear workspace
#rm(list=ls())
#set.seed(19)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#We will need the following libraries:
library(network)
library(ergm)
library(gtools)
source('functions.R')

#This undirected adjacency matrix includes 'all' relationships between characters in Grey's Anatomy. 
#Read in the adj matrix without the first column
grey <- read.csv("Grey_sex.csv",header=T)[,-1]

#Read in node attributes, which include Gender, Race, Birth Year, Professional Position, Season he/she first appears, and astrological sign.
attributes <- read.csv("Grey_Attributes.csv",header=T)

#Create a network object using the sociomatrix and its attributes
ga.net<-network(grey, vertex.attrnames=colnames(attributes),
                directed=F, hyper=F, loops=F, multiple=F, bipartite=F)

#Set vertex attributes
set.vertex.attribute(ga.net,names(attributes),attributes)

#define model
model1 <-as.formula(ga.net~edges+nodematch("sex")+degree(1)+nodematch("race")+
                           absdiff("birthyear"))

#create data with NAs
df_NAs <- insert_NAs()

#create data with NAs
df_NAs <- insert_NAs()

#create data with NAs
df_NAs <- insert_NAs()


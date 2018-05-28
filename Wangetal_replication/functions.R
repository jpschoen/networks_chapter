#1) We select some specified number, k, of observed edges (1s) and nulls (0s) from y-obs and replace 
#   them with NAs to produce a modified y-obsH. 

#complete: has covarate data, but missing all DV edges for some actors
#partial: some edges missing for some actors


# feed in adjaceny matrix
insert_NAs <-function(data, type = c("partial", "complete", "both"), partial_level = .1, complete_level = .1){
  if(type == "partial"){
    data_NA <- data
    pp <- seq(1,nrow(data),1)
    pp <- permutations(n=length(pp),r=2,v=pp,repeats.allowed=F)
    pp_sample <- pp[sample(nrow(pp), round(partial_level*nrow(pp))), ]
    for(i in 1:nrow(pp_sample)){
      data_NA[pp_sample[i,1],pp_sample[i,2]] <-NA
    }
    return(data_NA)
  }else if (type == "complete"){
    data_NA <- data
    cp_sample <- sample(nrow(data), round(partial_level*nrow(data)))
    for(i in cp_sample){
      data_NA[i,] <-NA
    }
    return(data_NA)
  }else if (type == "both"){
    data_NA <- data
    cp_sample <- sample(nrow(data_NA), round(partial_level*nrow(data_NA)))
    for(i in cp_sample){
      data_NA[i,] <-NA
    }    
    pp <- seq(1,nrow(data),1)
    pp <- permutations(n=length(pp),r=2,v=pp,repeats.allowed=F)
    pp_sample <- pp[sample(nrow(pp), round(partial_level*nrow(pp))), ]
    for(i in 1:nrow(pp_sample)){
      data_NA[pp_sample[i,1],pp_sample[i,2]] <-NA
    }
    return(data_NA)  
  }
}


# 2) Using the selected ERGM family, we estimate θH|yobsH. 
# 3) We simulate D multiple draws from YmisH|θH,yobsH. The number of draws 
#    (D) can be determined by the researcher 


ergm_predictions <- function(model){
  print("Running initial ERGM")
  model_a<-ergm(model,control=control.ergm(
                    MCMC.burnin=50000, MCMC.interval=5000))
  model_b<-ergm(model,control=control.ergm(
                   MCMC.burnin=50000, MCMC.interval=5000,init=coef(model_a),MCMLE.last.boost=1))
  print("Application of convergence criterions")
  sf <- model_b$sample-model_b$sample.obs
  t.k <- abs(apply(sf,2,mean))/apply(model_b$sample,2,sd)
  tconv.max <- sqrt(t(apply(sf,2,mean) %*% solve(as.matrix(cov(sf))) %*% apply(sf,2,mean)))
  while (max(t.k)>0.1 | tconv.max>0.25) {
    par <- coef(model_b)
    model_b <- ergm(model,control=control.ergm(
                       MCMC.burnin=50000, MCMC.interval=5000,init=coef(model_a),MCMLE.last.boost=1))
    sf <- model_b$sample-model_b$sample.obs
    t.k <- abs(apply(sf,2,mean)/apply(model_b$sample,2,sd))
    tconv.max <- sqrt(t(apply(sf,2,mean) %*% solve(as.matrix(cov(sf))) %*% apply(sf,2,mean)))
  }
  print("250 Network Simulations")
  net.fit<-model_b
  net.sim<-simulate(net.fit,
                     constraints=~observed,
                     nsim=250)
  return(net.sim)
}


#4) Using the simulated draws, we estimate the predictive accuracy of the imputation procedure via 
#   the mean fraction of held-out 1s and 0s (jointly or respectively) that are correctly imputed.

return_accuracy <-function(simulations, observed){
  score <-vector("numeric",length(simulations))
  for(i in 1:length(simulations)){
   test = as.matrix.network(simulations[[i]])-observed
   score[i] = sum(abs(test))/(nrow(observed)^2)
  }
  return(score)
}




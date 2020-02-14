# Generate 1000 Datasets for five different population parameter settings and each for nine different pool size
# in total 9*5*1000 datasets

library(stochprofML)

#setwd("SimulationStudies/1_GenerateDatasets/")
pool_sizes<-c(1,2,5,10,15,20,50, "mix", "mix2")
nn <- 1000


# in which parameter setting are you in?
for(whichSet in c("Set1","Set2","Set3","Set4","Set5")){
  #setwd(whichSet)
  numb_datasets<-9000
  set.seed(123)
  sample_seed <- sample(1:1000000, numb_datasets, replace=FALSE)
  counter_seed <- 1
  
  # Parameters of Set 1
  k = 50
  TY=2
  p = c(0.2, 0.8)
  mu.list <- list()
  m <- 1
  mu <- c(2, 0)
  mu.list[[m]] <- mu
  sigma <- 0.2
  
  # Change in Parameters for Set 2,3,4 and 5
  if(whichSet=="Set2"){
    p = c(0.1, 0.9)
  }else if(whichSet=="Set3"){
    p = c(0.4, 0.6)
  }else if(whichSet=="Set4"){
    mu <- c(2, 1)
  }else if(whichSet=="Set5"){
    sigma <- 0.5
  }
  
  # Which poolsize n do you want?
  for(j in pool_sizes){
    n<-j
    if(n=="mix"){
      geg_n <- c(rep(c(1,2,5,10), 12), 1, 2)
    }else if(n=="mix2"){
      geg_n <- c(rep(c(10,15,20,50), 12), 10, 15)
    }else{
      n = as.numeric(n)
    }
    
    
  #Generate 1000 dataset for each pool size
  generated_data <- list() 
  if(j %in% pool_sizes[8:9]){
    n = geg_n  
  }   
  
      for (i in 1:nn){
        set.seed(sample_seed[counter_seed])
  
        dataset <- matrix(0, ncol = m, nrow = k)
  
        stochprofML ::: set.model.functions("LN-LN")
  
        dataset[, m] <- r.sum.of.mixtures.LNLN(k,n,p,mu.list[[m]],rep(sigma,TY))
  
        name_data <- paste("Dataset",i, sep = "_")
        generated_data[[i]]<- dataset
        names(generated_data)[i]<-name_data
  
        counter_seed= counter_seed+1
      }
      #Save the list
      Var_name <- paste(whichSet, j, sep = "_")
      Var_rda <- paste(Var_name, "rda", sep = ".")
  
      save(generated_data, file = Var_rda)
  
  }
  setwd("..")
}
  
  


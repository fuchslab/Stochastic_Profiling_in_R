# Fit LNLN to all datasets of all datasets

library(stochprofML)
path<- getwd()

pool_sizes<-c(1,2,5,10,15,20,50, "mix", "mix2")
nn <- 1000

# in which parameter setting are you in?
for(whichSet in c("Set1","Set2","Set3","Set4","Set5")){
  
  set.seed(456)
  sample_seed <- sample(1:1000000, 9*nn, replace=FALSE)
  counter_seed <- 1
  
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
  name_set <- paste(whichSet,n, sep = "_")
  name_set1 <- paste(name_set, "rda", sep = ".")
  load(paste(path, name_set1, sep="/"))
  
  # Name for the saving the results
  Var_name_z <- paste(whichSet, whichSet, sep = "/")
  Var_name <- paste(Var_name_z, j, sep = "_")
  Var_name2 <- paste(Var_name, "fitting", sep = "_")
  Var_rda <- paste(Var_name2, "rda", sep = ".")
  
  # list for all 1000 fits
  result_fit_list<-list()
  
  # 1000 datasets for each pool-size setting
  if(j %in% pool_sizes[8:9]){
    n = geg_n  
  }
  for (i in 1:nn){
    
    #2 Population, 1 Gene should be fitted
    TY = 2
    m <- 1
    genes <- paste("Gene",1:m)
    
    # use specific dataset
    dataset<- generated_data[[i]]
    
    # set fit specific seed
    set.seed(sample_seed[counter_seed])
    
    sink( tempfile() )
    
    # Fitting of the Parameter
    result_fit_list[[i]] <- stochprof.loop(model="LN-LN",dataset=dataset,n=n,TY=TY,genenames=genes,fix.mu=F,loops=10,
                                           until.convergence=F,print.output=F,show.plots=F,plot.title="Entire Cluster",use.constraints=F)
    sink()
    
    # save results of Dataset 1 till 1000 in one list
    name_data <- paste("Dataset",i, sep = "_")
    names(result_fit_list)[i]<-name_data
    
    #set counter one up (in total total from 1 to 9000)
    counter_seed=counter_seed+1
    
    #save result for each setting
    save(result_fit_list, file = Var_rda)
    
  }
}
  


}

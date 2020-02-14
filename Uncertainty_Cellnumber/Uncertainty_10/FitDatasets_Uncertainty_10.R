# Fit LNLN to all datasets of all datasets

library(stochprofML)
path <- getwd()
# setwd('SimulationsStudien')

nn <- 1000

set.seed(987)
sample_seed <- sample(1:1e+06, 10 * nn, replace = FALSE)
counter_seed <- 1

# true n, load original dataset and original fits
n <- 10
load(file = paste(path, "2_FitDatasets/Set1/Set1_10_fitting.rda", sep = "/"))
load(file = paste(path, "1_GenerateDatasets/Set1/Set1_10.rda", sep = "/"))

whichSet <- "Set1"


for (j in c(5, 6, 7, 8, 9, 11, 12, 13, 14, 15)) {
  # Name for the saving the results
  Var_name <- paste(whichSet, "10", sep = "_")
  Var_name2 <- paste(Var_name, "fitting", j, sep = "_")
  Var_rda <- paste(Var_name2, "rda", sep = ".")
  
  result_fit_uncertain_list <- list()
  
  # 1000 datasets for each pool-size setting
  
  for (i in 1:nn) {
    
    # 2 Population, 1 Gene should be fitted
    TY <- 2
    m <- 1
    genes <- paste("Gene", 1:m)
    
    # use specific dataset
    dataset <- generated_data[[i]]
    
    # set fit specific seed
    set.seed(sample_seed[counter_seed])
    
    sink(tempfile())
    
    # Fitting of the Parameter
    result_fit_uncertain_list[[i]] <- stochprof.loop(model = "LN-LN", dataset = dataset, 
      n = j, TY = TY, genenames = genes, fix.mu = FALSE, loops = 10, until.convergence = FALSE, 
      print.output = FALSE, show.plots = FALSE, plot.title = "Entire Cluster", 
      use.constraints = FALSE)
    sink()
    
    # save results of Dataset 1 till 1000 in one list
    name_data <- paste("Dataset", i, sep = "_")
    names(result_fit_uncertain_list)[i] <- name_data
    
    # set counter one up (in total total from 1 to 9000)
    counter_seed <- counter_seed + 1
    
    # save result for each setting
    save(result_fit_uncertain_list, file = Var_rda)
    
  }
}




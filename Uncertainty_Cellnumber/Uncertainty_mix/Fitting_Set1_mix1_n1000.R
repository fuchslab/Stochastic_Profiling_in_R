# Use Dataset 1 from Set 1 underlying small mix vector Simulate 1000 similar n
# vectors and estimate parameters

library(stochprofML)

geg_n <- c(rep(c(1, 2, 5, 10), 12), 1, 2)

# Set parameters
i <- 1
TY <- 2
n <- geg_n
m <- 1
n_rep <- 1000
genes <- paste("Gene", 1:m)

# load dataset
Var_name <- paste("Set1", "mix", sep = "_")
Var_rda <- paste(Var_name, "rda", sep = ".")
load(file = Var_rda)

# all seeds that are used
set.seed(256)
sample_seed <- sample(1:1e+06, 9 * n_rep, replace = FALSE)
counter_seed <- 7001


result_fitting <- list()

for (j in 1:n_rep) {
  
  set.seed(j)
  n_sim <- rpois(length(n), n)
  n_sim[n_sim == 0] <- 1
  dataset <- generated_data[[1]]
  
  set.seed(sample_seed[counter_seed])
  sink(tempfile())
  result_fitting[[j]] <- stochprof.loop(model = "LN-LN", dataset = dataset, n = n_sim, 
    TY = TY, genenames = genes, fix.mu = F, loops = 10, until.convergence = F, 
    print.output = F, show.plots = F, plot.title = "Entire Cluster", use.constraints = F)
  sink()
  
  name_data <- paste("n_vector", j, sep = "_")
  names(result_fitting)[j] <- name_data
  
  save(result_fitting, file = "Fitting_mix1_n1000_complete.rda")
  
}


library(stochprofML)

set.seed(10)
k <- 1000
n <- 10
TY <- 2
p_2 <- c(0.6, 0.4)
p_3 <- c(0.5, 0.1, 0.4)
mu_1 <- 0.47
mu_2 <- -0.87
mu_3 <- 0.2
sigma_1 <- 0.03
sigma_2 <- 0.4
sigma_3 <- 0.7
lambda = 5


#### EXPLN 3 Pop
set.seed(35)
gene_EXPLN_3 <- r.sum.of.mixtures.EXPLN(k = k, n, p_3, c(mu_1,mu_2), rep(sigma_1, 2), lambda)

# Infer parameters
library(microbenchmark)

micro_stochprof_EXPLN_3 <- function(){
  stochprof.loop(model = "EXP-LN", 
                 dataset = matrix(gene_EXPLN_3, ncol = 1), n = n, TY = 3, 
                 genenames = "SimGene", fix.mu = FALSE, loops = 5,
                 until.convergence = TRUE, print.output = FALSE, show.plots = FALSE,
                 plot.title = "Simulated Gene",	use.constraints = FALSE)
}


set.seed(56786)
res_EXPLN_3 <- microbenchmark(micro_stochprof_EXPLN_3(), times = 5)

save(gene_EXPLN_2, res_EXPLN_2, file="EXPLN_2_5times.rda")

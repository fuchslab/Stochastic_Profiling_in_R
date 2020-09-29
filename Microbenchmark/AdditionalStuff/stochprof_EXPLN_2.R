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


#### EXPLN 2 Pop
set.seed(30)
gene_EXPLN_2 <- r.sum.of.mixtures.EXPLN(k = k, n, p_2, mu_1, sigma_1, lambda)


# Infer parameters
library(microbenchmark)

micro_stochprof_EXPLN_2 <- function(){
  stochprof.loop(model = "EXP-LN", 
                 dataset = matrix(gene_EXPLN_2, ncol = 1), n = n, TY = 2, 
                 genenames = "SimGene", fix.mu = FALSE, loops = 5,
                 until.convergence = TRUE, print.output = FALSE, show.plots = FALSE,
                 plot.title = "Simulated Gene",	use.constraints = FALSE)
}


set.seed(2345)
res_EXPLN_2 <- microbenchmark(micro_stochprof_EXPLN_2(), times = 1000)


save(gene_EXPLN_2, res_EXPLN_2, file="EXPLN_2.rda")
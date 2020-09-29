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


#### LNLN 3 Pop
set.seed(15)
gene_LNLN_3 <- r.sum.of.mixtures.LNLN(k = k, n, p_3, c(mu_1, mu_2, mu_3), rep(sigma_1, 3))

# Infer parameters
library(microbenchmark)

micro_stochprof_LNLN_3 <- function(){
  stochprof.loop(model = "LN-LN", 
                 dataset = matrix(gene_LNLN_3, ncol = 1), n = n, TY = 3, 
                 genenames = "SimGene", fix.mu = FALSE, loops = 5,
                 until.convergence = TRUE, print.output = FALSE, show.plots = FALSE,
                 plot.title = "Simulated Gene",	use.constraints = FALSE)
}

set.seed(456)
res_LNLN_3 <- microbenchmark(micro_stochprof_LNLN_3(), times = 5)

save(gene_LNLN_3, res_LNLN_3, file="LNLN_3_5times.rda")
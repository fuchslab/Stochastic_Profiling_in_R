library(stochprofML)

# Fix Distribution parameters to simulate data
k <- 1000
n <- 10
TY <- 3
p_3 <- c(0.1, 0.4, 0.5)
mu_1 <- 1.5
mu_2 <- -0.4
mu_3 <- -2.5
sigma <- 0.2
trials <- 1000


gene_LNLN_3 <- list()
res_LNLN_3 <- list()

# Generate 1000 datasets containing 3 Populations folowing the LNLN model. 
# All 1000 datasets have the same population parameters
# Estimate for each dataset once the parameters with stochprofML
for(i in 1 : trials){
set.seed(735 + i)
gene_LNLN_3[[i]] <- r.sum.of.mixtures.LNLN(k = k, n, p_3, c(mu_1, mu_2, mu_3), rep(sigma, 3))

set.seed(8974 + i)
res_LNLN_3[[i]] <-  stochprof.loop(model = "LN-LN", 
                 dataset = matrix(gene_LNLN_3[[i]], ncol = 1), n = n, TY = TY, 
                 genenames = "SimGene", fix.mu = FALSE, loops = 5,
                 until.convergence = TRUE, print.output = FALSE, show.plots = FALSE,
                 plot.title = "Simulated Gene",	use.constraints = FALSE)
save(gene_LNLN_3, res_LNLN_3, file="stochprofML_LNLN_3.rda")
}




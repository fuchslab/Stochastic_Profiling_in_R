# Generate Data
library(stochprofML)

set.seed(10)
k <- 1000
n <- 10
TY <- 2
p <- c(0.62, 0.38)
mu <- c(0.47, -0.87)
sigma <- 0.03
gene_LNLN <- r.sum.of.mixtures.LNLN(k = k, n, p, mu, rep(sigma, TY))



pdf(paste0("Simulated_Gene.pdf"), width = 7, height = 4)
x <- seq(from = min(gene_LNLN), to = max(gene_LNLN), length = 500)
stochprofML:::set.model.functions("LN-LN")
y <- d.sum.of.mixtures(x, n, p, mu, rep(sigma, TY), logdens = FALSE)
hist(gene_LNLN, main = paste("Simulated Gene"), breaks = 50,
    xlab = "Sum of mixtures of lognormals", ylab = "Density",
    freq = FALSE, col = "lightgrey")
lines(x, y, col = "blue", lwd = 2)
legend("topright", legend = "data generating pdf", col = "blue", lwd = 2, bty = "n")
dev.off()

# Infer parameters
library(stochprofML)

pdf(paste0("Simulated_Gene_4.pdf"))
par(mfrow(c(2,2)))
set.seed(20)
result <- stochprof.loop(model = "LN-LN", 
    dataset = matrix(gene_LNLN, ncol = 1), n = n, TY = TY, 
    genenames = "SimGene", fix.mu = FALSE, loops = 10,
    until.convergence = FALSE, print.output = FALSE, show.plots = TRUE,
    plot.title = "Simulated Gene",	use.constraints = FALSE)
dev.off()
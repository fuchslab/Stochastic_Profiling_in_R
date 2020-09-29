# Plot Figure 3

# load stochprofML results
load("stochprofML_LNLN_3.rda")

# write results in parameter matrices with confidence intervals
p_1 <- matrix(rep(0, 3 * 1000), nrow = 3)
p_2 <- matrix(rep(0, 3 * 1000), nrow = 3)
mu_1 <- matrix(rep(0, 3 * 1000), nrow = 3)
mu_2 <- matrix(rep(0, 3 * 1000), nrow = 3)
mu_3 <- matrix(rep(0, 3 * 1000), nrow = 3)
sigma <- matrix(rep(0, 3 * 1000), nrow = 3)

for (j in 1:1000) {
  p_1[2, j] <- res_LNLN_3[[j]]$mle[1]
  try(p_1[1, j] <- res_LNLN_3[[j]]$ci["p_1", 1], silent = TRUE)
  try(p_1[3, j] <- res_LNLN_3[[j]]$ci["p_1", 2], silent = TRUE)
  
  p_2[2, j] <- res_LNLN_3[[j]]$mle[2]
  try(p_2[1, j] <- res_LNLN_3[[j]]$ci["p_2", 1], silent = TRUE)
  try(p_2[3, j] <- res_LNLN_3[[j]]$ci["p_2", 2], silent = TRUE)
  
  mu_1[2, j] <- res_LNLN_3[[j]]$mle[3]
  try(mu_1[1, j] <- res_LNLN_3[[j]]$ci["mu_1_gene_SimGene", 1], silent = TRUE)
  try(mu_1[3, j] <- res_LNLN_3[[j]]$ci["mu_1_gene_SimGene", 2], silent = TRUE)
  
  mu_2[2, j] <- res_LNLN_3[[j]]$mle[4]
  try(mu_2[1, j] <- res_LNLN_3[[j]]$ci["mu_2_gene_SimGene", 1], silent = TRUE)
  try(mu_2[3, j] <- res_LNLN_3[[j]]$ci["mu_2_gene_SimGene", 2], silent = TRUE)
  
  mu_3[2, j] <- res_LNLN_3[[j]]$mle[5]
  try(mu_3[1, j] <- res_LNLN_3[[j]]$ci["mu_3_gene_SimGene", 1], silent = TRUE)
  try(mu_3[3, j] <- res_LNLN_3[[j]]$ci["mu_3_gene_SimGene", 2], silent = TRUE)
  
  sigma[2, j] <- res_LNLN_3[[j]]$mle[6]
  try(sigma[1, j] <- res_LNLN_3[[j]]$ci["sigma", 1], silent = TRUE)
  try(sigma[3, j] <- res_LNLN_3[[j]]$ci["sigma", 2], silent = TRUE)
}

#plot vioplots

pdf(paste0("Simstudy3Pop_LNLN.pdf"), width = 10, height = 4)
par(mfrow = c(1, 6), oma = c(2, 2, 2, 0), mar = c(3, 3, 3, 2), mgp = c(4, 1, 0), 
    bty = "n")

vioplot::vioplot(list(p_1[2,]), names = c("", ""), col = c("lightgrey"), ylab = "", yaxt = "n", axes = F, ylim = c(0, 1))
lines(c(0.5, 5.5), c(0.1, 0.1), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(p[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(p_2[2,]), names = c("", ""), col = c("lightgrey"), ylab = "", yaxt = "n", axes = F, ylim = c(0, 1))
lines(c(0.5, 5.5), c(0.4, 0.4), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(p[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(sigma[2,]), names = c(""), col = c("lightgrey"), ylab = "", yaxt = "n", axes = F, ylim = c(0, 5))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = c(0, 2, 4))
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)


vioplot::vioplot(list(mu_1[2,]), names = c(""), col = c("lightgrey"), ylab = "", yaxt = "n", axes = F, ylim = c(0, 10))
lines(c(0.5, 5.5), c(1.5, 1.5), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_2[2,]), names = c(""), col = c("lightgrey"), ylab = "", yaxt = "n", axes = F, ylim = c(-4, 4))
lines(c(0.5, 5.5), c(-0.4, -0.4), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_3[2,]), names = c(""), col = c("lightgrey"), ylab = "", yaxt = "n", axes = F)
lines(c(0.5, 5.5), c(-2.5, -2.5), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[3]), side = 2, line = 3, las = 1, cex = 1.5)
legend(x = 0, y = 10, legend = "true value", lw = 2, lty = 2, col = "#FF7F00", bty = "n", 
       xpd = NA, cex = 1.5)

title(paste0("Parameter estimates for 3 populations on 1000 datasets"), outer = TRUE, 
      cex.main = 2)

dev.off()

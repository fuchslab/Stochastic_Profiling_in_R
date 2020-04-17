# Use all 10-cell data sets from Set 1 and fit them with other pool sizes
n_uncertain <- c(5, 6, 7, 8, 9, 11, 12, 13, 14, 15)

# Set1 and parameters
whichSet <- "Set1"

p <- matrix(rep(0, 11 * 1000), nrow = 11)
mu_1 <- matrix(rep(0, 11 * 1000), nrow = 11)
mu_2 <- matrix(rep(0, 11 * 1000), nrow = 11)
sigma <- matrix(rep(0, 11 * 1000), nrow = 11)


for (n in n_uncertain) {
  Var_name <- paste(whichSet, "10", sep = "_")
  Var_name2 <- paste(Var_name, "fitting", n, sep = "_")
  Var_rda <- paste(Var_name2, "rda", sep = ".")
  load(file = Var_rda)
  
  i <- n - 4
  for (j in 1:1000) {
    p[i, j] <- result_fit_uncertain_list[[j]]$mle[1]
    mu_1[i, j] <- result_fit_uncertain_list[[j]]$mle[2]
    mu_2[i, j] <- result_fit_uncertain_list[[j]]$mle[3]
    sigma[i, j] <- result_fit_uncertain_list[[j]]$mle[4]
  }
}

# p-plots setwd('..')
load("SimulationStudies/2_FitDatasets/Set1/Set1_10_fitting.rda")

for (j in 1:1000) {
  p[6, j] <- result_fit_list[[j]]$mle[1]
  mu_1[6, j] <- result_fit_list[[j]]$mle[2]
  mu_2[6, j] <- result_fit_list[[j]]$mle[3]
  sigma[6, j] <- result_fit_list[[j]]$mle[4]
}


# Preparation done use colors FF7F00 #true -> orange color


# Create Plots True values
p_true <- c(0.2, 0.8)
mu_true <- c(2, 0)
sigma_true <- 0.2

# find 2.5%, 97,5% and median of the fitted parameters of all 1000 datasets for one cell number


p_all_3<-apply(p,1,quantile,probs=c(0.025,0.5,0.975))
mu1_all_3<-apply(mu_1,1,quantile,probs=c(0.025,0.5,0.975))
mu2_all_3<-apply(mu_2,1,quantile,probs=c(0.025,0.5,0.975))
sigma_all_3<-apply(sigma,1,quantile,probs=c(0.025,0.5,0.975))
# plot figure

pdf(paste0("Uncertainty_10.pdf"), width = 10, height = 6)

par(mfrow = c(2, 2), oma = c(2, 1, 2, 1), mar = c(4, 5, 3, 0), mgp = c(3, 1, 0), 
  bty = "n")

EnvStats::errorBar(x = 5:15, y = p_all_3[ 2,], lower = p_all_3[1,], upper = p_all_3[ 
  3,], incr = FALSE, gap = FALSE, bar.ends.size = 0.2, pch = 19, xlab = "", ylab = "", 
  ylim = c(0, 1), yaxt = "n", xaxt = "n", cex.axis = 1.3, cex.lab = 1.7)
axis(side = 2, at = seq(0, 1, by = 0.2), seq(0, 1, by = 0.2), cex.axis = 1.3, tick = TRUE, 
  las = 1)
axis(side = 1, at = 5:15, 5:15, cex.axis = 1.3, tick = TRUE)
abline(h = p_true[1], lw = 2, lty = 2, col = "#FF7F00")
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)
mtext("Pool size", side = 1, line = 3, las = 1, cex = 1.5)

EnvStats::errorBar(x = 5:15, y = mu1_all_3[ 2,], lower = mu1_all_3[ 1,], upper = mu1_all_3[ 
  3,], incr = FALSE, gap = FALSE, bar.ends.size = 0.2, pch = 19, xlab = "", ylab = "", 
  yaxt = "n", xaxt = "n", ylim = c(1, 3), cex.axis = 1.3, cex.lab = 1.7)
axis(side = 2, at = seq(1, 3, by = 1), seq(1, 3, by = 1), cex.axis = 1.3, tick = TRUE, 
  las = 1)
axis(side = 1, at = 5:15, 5:15, cex.axis = 1.3, tick = TRUE)
abline(h = mu_true[1], lw = 2, lty = 2, col = "#FF7F00")
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)
mtext("Pool size", side = 1, line = 3, las = 1, cex = 1.5)
legend(x = 10, y = 3.5, legend = "true value", lw = 2, lty = 2, col = "#FF7F00", 
  bty = "n", xpd = NA, cex = 1.3)

EnvStats::errorBar(x = 5:15, y = mu2_all_3[ 2,], lower = mu2_all_3[1,], upper = mu2_all_3[ 
  3,], incr = FALSE, gap = FALSE, bar.ends.size = 0.2, pch = 19, xlab = "", ylab = "", 
  yaxt = "n", xaxt = "n", ylim = c(-2.5,1.5), cex.axis = 1.3, cex.lab = 1.7)
axis(side = 2, at = seq(-2, 1, by = 1), seq(-2, 1, by = 1), cex.axis = 1.3, 
  tick = TRUE, las = 1)
axis(side = 1, at = 5:15, 5:15, cex.axis = 1.3, tick = TRUE)
abline(h = mu_true[2], lw = 2, lty = 2, col = "#FF7F00")
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)
mtext("Pool size", side = 1, line = 3, las = 1, cex = 1.5)

EnvStats::errorBar(x = 5:15, y = sigma_all_3[ 2,], lower = sigma_all_3[1,], upper = sigma_all_3[
  3,], incr = FALSE, gap = FALSE, bar.ends.size = 0.2, pch = 19, xlab = "", ylab = "", 
  yaxt = "n", xaxt = "n", ylim = c(0, 0.4), cex.axis = 1.3, cex.lab = 1.7)
axis(side = 1, at = 5:15, 5:15, cex.axis = 1.3, tick = TRUE)
axis(side = 2, at = seq(0, 0.4, by = 0.2), seq(0, 0.4, by = 0.2), cex.axis = 1.3, tick = TRUE, 
  las = 1)
abline(h = sigma_true[1], lw = 2, lty = 2, col = "#FF7F00")
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)
mtext("Pool size", side = 1, line = 3, las = 1, cex = 1.5)

dev.off()



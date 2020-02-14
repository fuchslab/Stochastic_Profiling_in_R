# Simstudy Uncertainty in n-vektor Plot results of the simulation studies for all
# Parameter sets and pool compositions

# setwd('UncertaintyCellnumber\Uncertainty_mix')


# Load fitted data for each Set and create parameter vectors Set1
load(file = "2_FitDatasets/Set1/Complete_fitting.rda")

p <- matrix(rep(0, 9 * 1000), nrow = 9)
mu_1 <- matrix(rep(0, 9 * 1000), nrow = 9)
mu_2 <- matrix(rep(0, 9 * 1000), nrow = 9)
sigma <- matrix(rep(0, 9 * 1000), nrow = 9)

p_CI_lower <- matrix(rep(0, 9 * 1000), nrow = 9)
mu_1_CI_lower <- matrix(rep(0, 9 * 1000), nrow = 9)
mu_2_CI_lower <- matrix(rep(0, 9 * 1000), nrow = 9)
sigma_CI_lower <- matrix(rep(0, 9 * 1000), nrow = 9)


p_CI_upper <- matrix(rep(0, 9 * 1000), nrow = 9)
mu_1_CI_upper <- matrix(rep(0, 9 * 1000), nrow = 9)
mu_2_CI_upper <- matrix(rep(0, 9 * 1000), nrow = 9)
sigma_CI_upper <- matrix(rep(0, 9 * 1000), nrow = 9)


for (i in 1:9) {
  
  for (j in 1:1000) {
    
    p[i, j] <- complete_fitting[[i]][[j]]$mle[1]
    try(p_CI_lower[i, j] <- complete_fitting[[i]][[j]]$ci["p_1", 1], silent = TRUE)
    try(p_CI_upper[i, j] <- complete_fitting[[i]][[j]]$ci["p_1", 2], silent = TRUE)
    
    mu_1[i, j] <- complete_fitting[[i]][[j]]$mle[2]
    try(mu_1_CI_lower[i, j] <- complete_fitting[[i]][[j]]$ci["mu_1_gene_Gene 1", 
      1], silent = TRUE)
    try(mu_1_CI_upper[i, j] <- complete_fitting[[i]][[j]]$ci["mu_1_gene_Gene 1", 
      2], silent = TRUE)
    
    mu_2[i, j] <- complete_fitting[[i]][[j]]$mle[3]
    try(mu_2_CI_lower[i, j] <- complete_fitting[[i]][[j]]$ci["mu_2_gene_Gene 1", 
      1], silent = TRUE)
    try(mu_2_CI_upper[i, j] <- complete_fitting[[i]][[j]]$ci["mu_2_gene_Gene 1", 
      2], silent = TRUE)
    
    sigma[i, j] <- complete_fitting[[i]][[j]]$mle[4]
    try(sigma_CI_lower[i, j] <- complete_fitting[[i]][[j]]$ci["sigma", 1], silent = TRUE)
    try(sigma_CI_upper[i, j] <- complete_fitting[[i]][[j]]$ci["sigma", 2], silent = TRUE)
  }
}

# p-plots
p_set1 <- data.frame(p_lower = c(t(p_CI_lower)), p = c(t(p)), p_upper = c(t(p_CI_upper)))
mu_1_set1 <- data.frame(mu_1_lower = c(t(mu_1_CI_lower)), mu_1 = c(t(mu_1)), mu_1_upper = c(t(mu_1_CI_upper)))
mu_2_set1 <- data.frame(mu_2_lower = c(t(mu_2_CI_lower)), mu_2 = c(t(mu_2)), mu_2_upper = c(t(mu_2_CI_upper)))
sigma_set1 <- data.frame(sigma_lower = c(t(sigma_CI_lower)), sigma = c(t(sigma)), 
  sigma_upper = c(t(sigma_CI_upper)))


# load Fit of mix1 of similar vectors

load("Fitting_mix1_n1000_complete.rda")

p_mix1 <- data.frame(p_lower = rep(0, 1000), p = rep(0, 1000), p_upper = rep(0, 1000))
mu_1_mix1 <- data.frame(mu_1_lower = rep(0, 1000), mu_1 = rep(0, 1000), mu_1_upper = rep(0, 
  1000))
mu_2_mix1 <- data.frame(mu_2_lower = rep(0, 1000), mu_2 = rep(0, 1000), mu_2_upper = rep(0, 
  1000))
sigma_mix1 <- data.frame(sigma_lower = rep(0, 1000), sigma = rep(0, 1000), sigma_upper = rep(0, 
  1000))

for (i in 1:1000) {
  p_mix1$p[i] <- result_fitting[[i]]$mle[1]
  try(p_mix1$p_lower[i] <- result_fitting[[i]]$ci[1, 1], silent = TRUE)
  try(p_mix1$p_upper[i] <- result_fitting[[i]]$ci[1, 2], silent = TRUE)
  mu_1_mix1$mu_1[i] <- result_fitting[[i]]$mle[2]
  try(mu_1_mix1$mu_1_lower[i] <- result_fitting[[i]]$ci[2, 1], silent = TRUE)
  try(mu_1_mix1$mu_1_upper[i] <- result_fitting[[i]]$ci[2, 2], silent = TRUE)
  mu_2_mix1$mu_2[i] <- result_fitting[[i]]$mle[3]
  try(mu_2_mix1$mu_2_lower[i] <- result_fitting[[i]]$ci[3, 1], silent = TRUE)
  try(mu_2_mix1$mu_2_upper[i] <- result_fitting[[i]]$ci[3, 2], silent = TRUE)
  sigma_mix1$sigma[i] <- result_fitting[[i]]$mle[4]
  try(sigma_mix1$sigma_lower[i] <- result_fitting[[i]]$ci[4, 1], silent = TRUE)
  try(sigma_mix1$sigma_upper[i] <- result_fitting[[i]]$ci[4, 2], silent = TRUE)
}

# load Fit of mix2 of similar vectors

load("Fitting_mix2_n1000_complete.rda")
p_mix2 <- data.frame(p_lower = rep(0, 1000), p = rep(0, 1000), p_upper = rep(0, 1000))
mu_1_mix2 <- data.frame(mu_1_lower = rep(0, 1000), mu_1 = rep(0, 1000), mu_1_upper = rep(0, 
  1000))
mu_2_mix2 <- data.frame(mu_2_lower = rep(0, 1000), mu_2 = rep(0, 1000), mu_2_upper = rep(0, 
  1000))
sigma_mix2 <- data.frame(sigma_lower = rep(0, 1000), sigma = rep(0, 1000), sigma_upper = rep(0, 
  1000))

for (i in 1:1000) {
  p_mix2$p[i] <- result_fitting[[i]]$mle[1]
  try(p_mix2$p_lower[i] <- result_fitting[[i]]$ci[1, 1], silent = TRUE)
  try(p_mix2$p_upper[i] <- result_fitting[[i]]$ci[1, 2], silent = TRUE)
  mu_1_mix2$mu_1[i] <- result_fitting[[i]]$mle[2]
  try(mu_1_mix2$mu_1_lower[i] <- result_fitting[[i]]$ci[2, 1], silent = TRUE)
  try(mu_1_mix2$mu_1_upper[i] <- result_fitting[[i]]$ci[2, 2], silent = TRUE)
  mu_2_mix2$mu_2[i] <- result_fitting[[i]]$mle[3]
  try(mu_2_mix2$mu_2_lower[i] <- result_fitting[[i]]$ci[3, 1], silent = TRUE)
  try(mu_2_mix2$mu_2_upper[i] <- result_fitting[[i]]$ci[3, 2], silent = TRUE)
  sigma_mix2$sigma[i] <- result_fitting[[i]]$mle[4]
  try(sigma_mix2$sigma_lower[i] <- result_fitting[[i]]$ci[4, 1], silent = TRUE)
  try(sigma_mix2$sigma_upper[i] <- result_fitting[[i]]$ci[4, 2], silent = TRUE)
}


########################################################################################## 




pdf(paste0("Simstudy3_mix1_2_result_vio.pdf"), width = 12, height = 5)

par(mfrow = c(1, 4), oma = c(2, 2, 2, 0), mar = c(3, 3, 3, 2), mgp = c(4, 1, 0), 
  bty = "n")

vioplot::vioplot(list(p_mix1$p, p_mix2$p), names = c(""), main = "", yaxt = "n", 
  ylab = "", col = "lightgrey")
axis(side = 2, at = seq(0, 1, by = 0.2), seq(0, 1, by = 0.2), cex.axis = 1.3, las = 1)
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)
lines(c(0.6, 1.4), c(0.2, 0.2), col = "#FF7F00", lw = 2, lty = 2, )
lines(c(0.6, 1.4), c(p_set1$p[7001], p_set1$p[7001]), col = "#1F78B4", lw = 2, lty = 2, 
  )
lines(c(1.6, 2.4), c(0.2, 0.2), col = "#FF7F00", lw = 2, lty = 2, )
lines(c(1.6, 2.4), c(p_set1$p[8001][1], p_set1$p[8001][1]), col = "#1F78B4", lw = 2, 
  lty = 2, )
axis(side = 1, at = c(1, 2), c("Pool:
1-2-5-10
cells", "Pool:
10-15-20-50
cells"), 
  cex.axis = 1.3, tick = FALSE, outer = TRUE)

vioplot::vioplot(list(mu_1_mix1$mu_1, mu_1_mix2$mu_1), names = c(""), main = "", 
  yaxt = "n", ylab = "", col = "lightgrey")
axis(side = 2, at = seq(0, 4, by = 1), seq(0, 4, by = 1), cex.axis = 1.3, las = 1)
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)
lines(c(0.6, 1.4), c(2, 2), col = "#FF7F00", lw = 2, lty = 2, )
lines(c(0.6, 1.4), c(mu_1_set1$mu_1[7001], mu_1_set1$mu_1[7001]), col = "#1F78B4", 
  lw = 2, lty = 2, )
lines(c(1.6, 2.4), c(2, 2), col = "#FF7F00", lw = 2, lty = 2, )
lines(c(1.6, 2.4), c(mu_1_set1$mu_1[8001], mu_1_set1$mu_1[8001]), col = "#1F78B4", 
  lw = 2, lty = 2, )
axis(side = 1, at = c(1, 2), c("Pool:
1-2-5-10
cells", "Pool:
10-15-20-50
cells"), 
  cex.axis = 1.3, tick = FALSE, outer = TRUE)

vioplot::vioplot(list(mu_2_mix1$mu_2, mu_2_mix2$mu_2), names = c(""), main = "", 
  yaxt = "n", ylab = "", col = "lightgrey", ylim = c(-40, 10))
axis(side = 2, at = seq(-40, 10, by = 10), seq(-40, 10, by = 10), cex.axis = 1.3, 
  las = 1)
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)
lines(c(0.6, 1.4), c(0, 0), col = "#FF7F00", lw = 2, lty = 2, )
lines(c(0.6, 1.4), c(mu_2_set1$mu_2[7001], mu_2_set1$mu_2[7001]), col = "#1F78B4", 
  lw = 2, lty = 2, )
lines(c(1.6, 2.4), c(0, 0), col = "#FF7F00", lw = 2, lty = 2, )
lines(c(1.6, 2.4), c(mu_2_set1$mu_2[8001], mu_2_set1$mu_2[8001]), col = "#1F78B4", 
  lw = 2, lty = 2, )

axis(side = 1, at = c(1, 2), c("Pool:
1-2-5-10
cells", "Pool:
10-15-20-50
cells"), 
  cex.axis = 1.3, tick = FALSE, outer = TRUE)

vioplot::vioplot(list(sigma_mix1$sigma, sigma_mix2$sigma), names = c(""), main = "", 
  yaxt = "n", ylab = "", col = "lightgrey", ylim = c(0, 1.4))
axis(side = 2, at = seq(0, 1.4, by = 0.2), seq(0, 1.4, by = 0.2), cex.axis = 1.3, 
  las = 1)
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)
lines(c(0.6, 1.4), c(0.2, 0.2), col = "#FF7F00", lw = 2, lty = 2, )
lines(c(0.6, 1.4), c(sigma_set1$sigma[7001], sigma_set1$sigma[7001]), col = "#1F78B4", 
  lw = 2, lty = 2, )
lines(c(1.6, 2.4), c(0.2, 0.2), col = "#FF7F00", lw = 2, lty = 2, )
lines(c(1.6, 2.4), c(sigma_set1$sigma[8001], sigma_set1$sigma[8001]), col = "#1F78B4", 
  lw = 2, lty = 2, )

axis(side = 1, at = c(1, 2), c("Pool:
1-2-5-10
cells", "Pool:
10-15-20-50
cells"), 
  cex.axis = 1.3, tick = FALSE, outer = TRUE)



legend(x = 0.5, y = 1.7, legend = c("true value", "estimate with true n-vetor"), 
  lty = c(2, 2), col = c("#FF7F00", "#1F78B4"), bty = "n", xpd = NA, cex = 1.3)


title(paste0("Parameter estimates for 1000 datasets"), outer = TRUE, cex.main = 2)


dev.off()


# Plot results of the simulation studies for all Parameter sets and pool
# compositions

# setwd('SimulationStudies/')

# Load fitted data for each Set and create parameter vectors Set1

n <- c(1, 2, 5, 10, 15, 20, 50, "mix of 1-2-5-10", "mix of 10-15-20-50")
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

# generate dataframes
p_set1 <- data.frame(p_lower = c(t(p_CI_lower)), p = c(t(p)), p_upper = c(t(p_CI_upper)))
mu_1_set1 <- data.frame(mu_1_lower = c(t(mu_1_CI_lower)), mu_1 = c(t(mu_1)), mu_1_upper = c(t(mu_1_CI_upper)))
mu_2_set1 <- data.frame(mu_2_lower = c(t(mu_2_CI_lower)), mu_2 = c(t(mu_2)), mu_2_upper = c(t(mu_2_CI_upper)))
sigma_set1 <- data.frame(sigma_lower = c(t(sigma_CI_lower)), sigma = c(t(sigma)), 
  sigma_upper = c(t(sigma_CI_upper)))



############################################## Set2

n <- c(1, 2, 5, 10, 15, 20, 50, "mix of 1-2-5-10", "mix of 10-15-20-50")
load(file = "2_FitDatasets/Set2/Complete_fitting.rda")


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
p_set2 <- data.frame(p_lower = c(t(p_CI_lower)), p = c(t(p)), p_upper = c(t(p_CI_upper)))
mu_1_set2 <- data.frame(mu_1_lower = c(t(mu_1_CI_lower)), mu_1 = c(t(mu_1)), mu_1_upper = c(t(mu_1_CI_upper)))
mu_2_set2 <- data.frame(mu_2_lower = c(t(mu_2_CI_lower)), mu_2 = c(t(mu_2)), mu_2_upper = c(t(mu_2_CI_upper)))
sigma_set2 <- data.frame(sigma_lower = c(t(sigma_CI_lower)), sigma = c(t(sigma)), 
  sigma_upper = c(t(sigma_CI_upper)))

############################################## Set3

n <- c(1, 2, 5, 10, 15, 20, 50, "mix of 1-2-5-10", "mix of 10-15-20-50")
load(file = "2_FitDatasets/Set3/Complete_fitting.rda")


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
p_set3 <- data.frame(p_lower = c(t(p_CI_lower)), p = c(t(p)), p_upper = c(t(p_CI_upper)))
mu_1_set3 <- data.frame(mu_1_lower = c(t(mu_1_CI_lower)), mu_1 = c(t(mu_1)), mu_1_upper = c(t(mu_1_CI_upper)))
mu_2_set3 <- data.frame(mu_2_lower = c(t(mu_2_CI_lower)), mu_2 = c(t(mu_2)), mu_2_upper = c(t(mu_2_CI_upper)))
sigma_set3 <- data.frame(sigma_lower = c(t(sigma_CI_lower)), sigma = c(t(sigma)), 
  sigma_upper = c(t(sigma_CI_upper)))

############################################## Set4

n <- c(1, 2, 5, 10, 15, 20, 50, "mix of 1-2-5-10", "mix of 10-15-20-50")
load(file = "2_FitDatasets/Set4/Complete_fitting.rda")



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

# 
p_set4 <- data.frame(p_lower = c(t(p_CI_lower)), p = c(t(p)), p_upper = c(t(p_CI_upper)))
mu_1_set4 <- data.frame(mu_1_lower = c(t(mu_1_CI_lower)), mu_1 = c(t(mu_1)), mu_1_upper = c(t(mu_1_CI_upper)))
mu_2_set4 <- data.frame(mu_2_lower = c(t(mu_2_CI_lower)), mu_2 = c(t(mu_2)), mu_2_upper = c(t(mu_2_CI_upper)))
sigma_set4 <- data.frame(sigma_lower = c(t(sigma_CI_lower)), sigma = c(t(sigma)), 
  sigma_upper = c(t(sigma_CI_upper)))

############################################## Set5

n <- c(1, 2, 5, 10, 15, 20, 50, "mix of 1-2-5-10", "mix of 10-15-20-50")
load(file = "2_FitDatasets/Set5/Complete_fitting.rda")


for (i in 1:9) {
  
  for (j in 1:1000) {
    # mle_complete[j]<-complete_fitting[[i]][[j]]$mle
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
p_set5 <- data.frame(p_lower = c(t(p_CI_lower)), p = c(t(p)), p_upper = c(t(p_CI_upper)))
mu_1_set5 <- data.frame(mu_1_lower = c(t(mu_1_CI_lower)), mu_1 = c(t(mu_1)), mu_1_upper = c(t(mu_1_CI_upper)))
mu_2_set5 <- data.frame(mu_2_lower = c(t(mu_2_CI_lower)), mu_2 = c(t(mu_2)), mu_2_upper = c(t(mu_2_CI_upper)))
sigma_set5 <- data.frame(sigma_lower = c(t(sigma_CI_lower)), sigma = c(t(sigma)), 
  sigma_upper = c(t(sigma_CI_upper)))
# Preparation done

########################################################################################################### Create Plots for different viewpoints/analysis etc

# use three colors in plots FF7F00 #true -> orange color 1F78B4 #estimated ->
# blue color 89BBBF # reference -> turquois color

# Part 1 # Impact of increasing Cellnumbers (for all sets) set 1


pdf(paste0("Simstudy1_set1_result_vio_sep.pdf"), width = 12, height = 7)
par(mfrow = c(4, 2), oma = c(2, 2, 2, 0), mar = c(3, 3, 3, 2), mgp = c(4, 1, 0), 
  bty = "n")

vioplot::vioplot(list(p_set1$p[1:1000], p_set1$p[1001:2000], p_set1$p[2001:3000], 
  p_set1$p[3001:4000], p_set1$p[7001:8000]), names = c("", "", "", "", ""), col = c(rep("lightgrey", 
  3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", axes = F, ylim = c(0, 0.6))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = c(0, 0.2, 0.4, 0.6))
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(p_set1$p[3001:4000], p_set1$p[4001:5000], p_set1$p[5001:6000], 
  p_set1$p[6001:7000], p_set1$p[8001:9000]), names = c("", "", "", "", ""), col = c("#89BBBF", 
  rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F)
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
legend(x = 4, y = 1.8, legend = "true value", lw = 2, lty = 2, col = "#FF7F00", bty = "n", 
  xpd = NA, cex = 1.5)
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_1_set1$mu_1[1:1000], mu_1_set1$mu_1[1001:2000], mu_1_set1$mu_1[2001:3000], 
  mu_1_set1$mu_1[3001:4000], mu_1_set1$mu_1[7001:8000]), names = c("", "", "", 
  "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", 
  axes = F)
lines(c(0.5, 5.5), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_1_set1$mu_1[3001:4000], mu_1_set1$mu_1[4001:5000], mu_1_set1$mu_1[5001:6000], 
  mu_1_set1$mu_1[6001:7000], mu_1_set1$mu_1[8001:9000]), names = c("", "", "", 
  "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F)
lines(c(0.5, 5.5), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_2_set1$mu_2[1:1000], mu_2_set1$mu_2[1001:2000], mu_2_set1$mu_2[2001:3000], 
  mu_2_set1$mu_2[3001:4000], mu_2_set1$mu_2[7001:8000]), names = c("", "", "", 
  "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", 
  axes = F, h = 0.3, ylim = c(-2, 1))
lines(c(0.5, 5.5), c(0, 0), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = c(-2, -1, 0, 1))
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_2_set1$mu_2[3001:4000], mu_2_set1$mu_2[4001:5000], mu_2_set1$mu_2[5001:6000], 
  mu_2_set1$mu_2[6001:7000], mu_2_set1$mu_2[8001:9000]), names = c("", "", "", 
  "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F, 
  h = 2, ylim = c(-50, 10))
lines(c(0.5, 5.5), c(0, 0), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(sigma_set1$sigma[1:1000], sigma_set1$sigma[1001:2000], sigma_set1$sigma[2001:3000], 
  sigma_set1$sigma[3001:4000], sigma_set1$sigma[7001:8000]), names = c("", "", 
  "", "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", 
  yaxt = "n", axes = F, ylim = c(0, 0.8))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 1, at = 1:5, c("single
cells", "2-cell
pools", "5-cell
pools", "10-cell
pools", 
  "1-2-5-
10-cell 
pools"), cex.axis = 1.5, tick = FALSE, outer = TRUE)
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(sigma_set1$sigma[3001:4000], sigma_set1$sigma[4001:5000], sigma_set1$sigma[5001:6000], 
  sigma_set1$sigma[6001:7000], sigma_set1$sigma[8001:9000]), names = c("", "", 
  "", "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", 
  axes = F, ylim = c(0, 1))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 1, at = 1:5, c("10-cell
pools", "15-cell
pools", "20-cell
pools", 
  "50-cell
pools", "10-15-20-
50-cell
pools"), cex.axis = 1.5, tick = FALSE, 
  outer = TRUE)
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)

title(paste0("Parameter estimates for different pool decompositions"), outer = TRUE, 
  cex.main = 2)
dev.off()

######################################################################################################## set 2


pdf(paste0("Simstudy1_set2_result_vio_sep.pdf"), width = 12, height = 7)
par(mfrow = c(4, 2), oma = c(2, 2, 2, 0), mar = c(3, 3, 3, 2), mgp = c(4, 1, 0), 
  bty = "n")

vioplot::vioplot(list(p_set2$p[1:1000], p_set2$p[1001:2000], p_set2$p[2001:3000], 
  p_set2$p[3001:4000], p_set2$p[7001:8000]), names = c("", "", "", "", ""), col = c(rep("lightgrey", 
  3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", axes = F, ylim = c(0, 0.6))
lines(c(0.5, 5.5), c(0.1, 0.1), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = c(0, 0.2, 0.4, 0.6))
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(p_set2$p[3001:4000], p_set2$p[4001:5000], p_set2$p[5001:6000], 
  p_set2$p[6001:7000], p_set2$p[8001:9000]), names = c("", "", "", "", ""), col = c("#89BBBF", 
  rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F)
lines(c(0.5, 5.5), c(0.1, 0.1), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
legend(x = 4, y = 1.6, legend = "true value", lw = 2, lty = 2, col = "#FF7F00", bty = "n", 
  xpd = NA, cex = 1.5)
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_1_set2$mu_1[1:1000], mu_1_set2$mu_1[1001:2000], mu_1_set2$mu_1[2001:3000], 
  mu_1_set2$mu_1[3001:4000], mu_1_set2$mu_1[7001:8000]), names = c("", "", "", 
  "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", 
  axes = F, ylim = c(0, 2.5))
lines(c(0.5, 5.5), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_1_set2$mu_1[3001:4000], mu_1_set2$mu_1[4001:5000], mu_1_set2$mu_1[5001:6000], 
  mu_1_set2$mu_1[6001:7000], mu_1_set2$mu_1[8001:9000]), names = c("", "", "", 
  "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F, 
  ylim = c(0, 3.5))
lines(c(0.5, 5.5), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = seq(0, 3, 1))
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_2_set2$mu_2[1:1000], mu_2_set2$mu_2[1001:2000], mu_2_set2$mu_2[2001:3000], 
  mu_2_set2$mu_2[3001:4000], mu_2_set2$mu_2[7001:8000]), names = c("", "", "", 
  "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", 
  axes = F, h = 0.5, ylim = c(-1, 0.2))
lines(c(0.5, 5.5), c(0, 0), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = seq(-1, 0.2, 0.2))
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_2_set2$mu_2[3001:4000], mu_2_set2$mu_2[4001:5000], mu_2_set2$mu_2[5001:6000], 
  mu_2_set2$mu_2[6001:7000], mu_2_set2$mu_2[8001:9000]), names = c("", "", "", 
  "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F, 
  h = 2, ylim = c(-40, 10))
lines(c(0.5, 5.5), c(0, 0), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(sigma_set2$sigma[1:1000], sigma_set2$sigma[1001:2000], sigma_set2$sigma[2001:3000], 
  sigma_set2$sigma[3001:4000], sigma_set2$sigma[7001:8000]), names = c("", "", 
  "", "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", 
  yaxt = "n", axes = F, ylim = c(0, 0.4))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 1, at = 1:5, c("single
cells", "2-cell
pools", "5-cell
pools", "10-cell
pools", 
  "1-2-5-
10-cell 
pools"), cex.axis = 1.5, tick = FALSE, outer = TRUE)
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(sigma_set2$sigma[3001:4000], sigma_set2$sigma[4001:5000], sigma_set2$sigma[5001:6000], 
  sigma_set2$sigma[6001:7000], sigma_set2$sigma[8001:9000]), names = c("", "", 
  "", "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", 
  axes = F, ylim = c(0, 1))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 1, at = 1:5, c("10-cell
pools", "15-cell
pools", "20-cell
pools", 
  "50-cell
pools", "10-15-20-
50-cell
pools"), cex.axis = 1.5, tick = FALSE, 
  outer = TRUE)
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)

title(paste0("Parameter estimates for different pool decompositions"), outer = TRUE, 
  cex.main = 2)
dev.off()

######################################################################################################## set 3


pdf(paste0("Simstudy1_set3_result_vio_sep.pdf"), width = 12, height = 7)
par(mfrow = c(4, 2), oma = c(2, 2, 2, 0), mar = c(3, 3, 3, 2), mgp = c(4, 1, 0), 
  bty = "n")

vioplot::vioplot(list(p_set3$p[1:1000], p_set3$p[1001:2000], p_set3$p[2001:3000], 
  p_set3$p[3001:4000], p_set3$p[7001:8000]), names = c("", "", "", "", ""), col = c(rep("lightgrey", 
  3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", axes = F, ylim = c(0, 1))
lines(c(0.5, 5.5), c(0.4, 0.4), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1))
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(p_set3$p[3001:4000], p_set3$p[4001:5000], p_set3$p[5001:6000], 
  p_set3$p[6001:7000], p_set3$p[8001:9000]), names = c("", "", "", "", ""), col = c("#89BBBF", 
  rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F)
lines(c(0.5, 5.5), c(0.4, 0.4), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
legend(x = 4, y = 1.8, legend = "true value", lw = 2, lty = 2, col = "#FF7F00", bty = "n", 
  xpd = NA, cex = 1.5)
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_1_set3$mu_1[1:1000], mu_1_set3$mu_1[1001:2000], mu_1_set3$mu_1[2001:3000], 
  mu_1_set3$mu_1[3001:4000], mu_1_set3$mu_1[7001:8000]), names = c("", "", "", 
  "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", 
  axes = F, ylim = c(1, 3))
lines(c(0.5, 5.5), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_1_set3$mu_1[3001:4000], mu_1_set3$mu_1[4001:5000], mu_1_set3$mu_1[5001:6000], 
  mu_1_set3$mu_1[6001:7000], mu_1_set3$mu_1[8001:9000]), names = c("", "", "", 
  "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F)
lines(c(0.5, 5.5), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_2_set3$mu_2[1:1000], mu_2_set3$mu_2[1001:2000], mu_2_set3$mu_2[2001:3000], 
  mu_2_set3$mu_2[3001:4000], mu_2_set3$mu_2[7001:8000]), names = c("", "", "", 
  "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", 
  axes = F, h = 0.5, ylim = c(-3, 2))
lines(c(0.5, 5.5), c(0, 0), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = seq(-3, 2, 1))
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_2_set3$mu_2[3001:4000], mu_2_set3$mu_2[4001:5000], mu_2_set3$mu_2[5001:6000], 
  mu_2_set3$mu_2[6001:7000], mu_2_set3$mu_2[8001:9000]), names = c("", "", "", 
  "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F, 
  h = 2, ylim = c(-50, 10))
lines(c(0.5, 5.5), c(0, 0), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(sigma_set3$sigma[1:1000], sigma_set3$sigma[1001:2000], sigma_set3$sigma[2001:3000], 
  sigma_set3$sigma[3001:4000], sigma_set3$sigma[7001:8000]), names = c("", "", 
  "", "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", 
  yaxt = "n", axes = F, ylim = c(0, 0.8))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 1, at = 1:5, c("single
cells", "2-cell
pools", "5-cell
pools", "10-cell
pools", 
  "1-2-5-
10-cell 
pools"), cex.axis = 1.5, tick = FALSE, outer = TRUE)
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(sigma_set3$sigma[3001:4000], sigma_set3$sigma[4001:5000], sigma_set3$sigma[5001:6000], 
  sigma_set3$sigma[6001:7000], sigma_set3$sigma[8001:9000]), names = c("", "", 
  "", "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", 
  axes = F, ylim = c(0, 1))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 1, at = 1:5, c("10-cell
pools", "15-cell
pools", "20-cell
pools", 
  "50-cell
pools", "10-15-20-
50-cell
pools"), cex.axis = 1.5, tick = FALSE, 
  outer = TRUE)
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)

title(paste0("Parameter estimates for different pool decompositions"), outer = TRUE, 
  cex.main = 2)
dev.off()


######################################################################################################## set 4




pdf(paste0("Simstudy1_set4_result_vio_sep.pdf"), width = 12, height = 7)
par(mfrow = c(4, 2), oma = c(2, 2, 2, 0), mar = c(3, 3, 3, 2), mgp = c(4, 1, 0), 
  bty = "n")

vioplot::vioplot(list(p_set4$p[1:1000], p_set4$p[1001:2000], p_set4$p[2001:3000], 
  p_set4$p[3001:4000], p_set4$p[7001:8000]), names = c("", "", "", "", ""), col = c(rep("lightgrey", 
  3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", axes = F, ylim = c(0, 1))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = seq(0, 1, 0.2))
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(p_set4$p[3001:4000], p_set4$p[4001:5000], p_set4$p[5001:6000], 
  p_set4$p[6001:7000], p_set4$p[8001:9000]), names = c("", "", "", "", ""), col = c("#89BBBF", 
  rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F)
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
legend(x = 4, y = 1.8, legend = "true value", lw = 2, lty = 2, col = "#FF7F00", bty = "n", 
  xpd = NA, cex = 1.5)
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_1_set4$mu_1[1:1000], mu_1_set4$mu_1[1001:2000], mu_1_set4$mu_1[2001:3000], 
  mu_1_set4$mu_1[3001:4000], mu_1_set4$mu_1[7001:8000]), names = c("", "", "", 
  "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", 
  axes = F, ylim = c(1, 3.5))
lines(c(0.5, 5.5), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_1_set4$mu_1[3001:4000], mu_1_set4$mu_1[4001:5000], mu_1_set4$mu_1[5001:6000], 
  mu_1_set4$mu_1[6001:7000], mu_1_set4$mu_1[8001:9000]), names = c("", "", "", 
  "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F, 
  ylim = c(1, 4))
lines(c(0.5, 5.5), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_2_set4$mu_2[1:1000], mu_2_set4$mu_2[1001:2000], mu_2_set4$mu_2[2001:3000], 
  mu_2_set4$mu_2[3001:4000], mu_2_set4$mu_2[7001:8000]), names = c("", "", "", 
  "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", 
  axes = F, h = 0.3, ylim = c(-1, 2))
lines(c(0.5, 5.5), c(1, 1), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = c(-1, 0, 1, 2))
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_2_set4$mu_2[3001:4000], mu_2_set4$mu_2[4001:5000], mu_2_set4$mu_2[5001:6000], 
  mu_2_set4$mu_2[6001:7000], mu_2_set4$mu_2[8001:9000]), names = c("", "", "", 
  "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F, 
  h = 2, ylim = c(-50, 10))
lines(c(0.5, 5.5), c(1, 1), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(sigma_set4$sigma[1:1000], sigma_set4$sigma[1001:2000], sigma_set4$sigma[2001:3000], 
  sigma_set4$sigma[3001:4000], sigma_set4$sigma[7001:8000]), names = c("", "", 
  "", "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", 
  yaxt = "n", axes = F, ylim = c(0, 0.6))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 1, at = 1:5, c("single
cells", "2-cell
pools", "5-cell
pools", "10-cell
pools", 
  "1-2-5-
10-cell 
pools"), cex.axis = 1.5, tick = FALSE, outer = TRUE)
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(sigma_set4$sigma[3001:4000], sigma_set4$sigma[4001:5000], sigma_set4$sigma[5001:6000], 
  sigma_set4$sigma[6001:7000], sigma_set4$sigma[8001:9000]), names = c("", "", 
  "", "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", 
  axes = F, ylim = c(0, 0.6))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 1, at = 1:5, c("10-cell
pools", "15-cell
pools", "20-cell
pools", 
  "50-cell
pools", "10-15-20-
50-cell
pools"), cex.axis = 1.5, tick = FALSE, 
  outer = TRUE)
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)

title(paste0("Parameter estimates for different pool decompositions"), outer = TRUE, 
  cex.main = 2)
dev.off()

######################################################################################################## set 5



pdf(paste0("Simstudy1_set5_result_vio_sep.pdf"), width = 12, height = 7)
par(mfrow = c(4, 2), oma = c(2, 2, 2, 0), mar = c(3, 3, 3, 2), mgp = c(4, 1, 0), 
  bty = "n")

vioplot::vioplot(list(p_set5$p[1:1000], p_set5$p[1001:2000], p_set5$p[2001:3000], 
  p_set5$p[3001:4000], p_set5$p[7001:8000]), names = c("", "", "", "", ""), col = c(rep("lightgrey", 
  3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", axes = F, ylim = c(0, 1))
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = seq(0, 1, 0.2))
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(p_set5$p[3001:4000], p_set5$p[4001:5000], p_set5$p[5001:6000], 
  p_set5$p[6001:7000], p_set5$p[8001:9000]), names = c("", "", "", "", ""), col = c("#89BBBF", 
  rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F)
lines(c(0.5, 5.5), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
legend(x = 4, y = 1.8, legend = "true value", lw = 2, lty = 2, col = "#FF7F00", bty = "n", 
  xpd = NA, cex = 1.5)
mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_1_set5$mu_1[1:1000], mu_1_set5$mu_1[1001:2000], mu_1_set5$mu_1[2001:3000], 
  mu_1_set5$mu_1[3001:4000], mu_1_set5$mu_1[7001:8000]), names = c("", "", "", 
  "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", 
  axes = F, ylim = c(0.5, 3.5))
lines(c(0.5, 5.5), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_1_set5$mu_1[3001:4000], mu_1_set5$mu_1[4001:5000], mu_1_set5$mu_1[5001:6000], 
  mu_1_set5$mu_1[6001:7000], mu_1_set5$mu_1[8001:9000]), names = c("", "", "", 
  "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F, 
  ylim = c(0, 4))
lines(c(0.5, 5.5), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_2_set5$mu_2[1:1000], mu_2_set5$mu_2[1001:2000], mu_2_set5$mu_2[2001:3000], 
  mu_2_set5$mu_2[3001:4000], mu_2_set5$mu_2[7001:8000]), names = c("", "", "", 
  "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", yaxt = "n", 
  axes = F, h = 0.3, ylim = c(-2, 1))
lines(c(0.5, 5.5), c(0, 0), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1, at = c(-2, -1, 0, 1))
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(mu_2_set5$mu_2[3001:4000], mu_2_set5$mu_2[4001:5000], mu_2_set5$mu_2[5001:6000], 
  mu_2_set5$mu_2[6001:7000], mu_2_set5$mu_2[8001:9000]), names = c("", "", "", 
  "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", axes = F, 
  h = 2, ylim = c(-50, 10))
lines(c(0.5, 5.5), c(0, 0), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(sigma_set5$sigma[1:1000], sigma_set5$sigma[1001:2000], sigma_set5$sigma[2001:3000], 
  sigma_set5$sigma[3001:4000], sigma_set5$sigma[7001:8000]), names = c("", "", 
  "", "", ""), col = c(rep("lightgrey", 3), "#89BBBF", "lightgrey"), ylab = "", 
  yaxt = "n", axes = F, ylim = c(0, 1.5))
lines(c(0.5, 5.5), c(0.5, 0.5), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 1, at = 1:5, c("single
cells", "2-cell
pools", "5-cell
pools", "10-cell
pools", 
  "1-2-5-
10-cell 
pools"), cex.axis = 1.5, tick = FALSE, outer = TRUE)
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)

vioplot::vioplot(list(sigma_set5$sigma[3001:4000], sigma_set5$sigma[4001:5000], sigma_set5$sigma[5001:6000], 
  sigma_set5$sigma[6001:7000], sigma_set5$sigma[8001:9000]), names = c("", "", 
  "", "", ""), col = c("#89BBBF", rep("lightgrey", 4)), ylab = "", yaxt = "n", 
  axes = F, ylim = c(0, 1.5))
lines(c(0.5, 5.5), c(0.5, 0.5), lw = 2, lty = 2, col = "#FF7F00")
axis(side = 1, at = 1:5, c("10-cell
pools", "15-cell
pools", "20-cell
pools", 
  "50-cell
pools", "10-15-20-
50-cell
pools"), cex.axis = 1.5, tick = FALSE, 
  outer = TRUE)
axis(side = 2, cex.axis = 1.3, las = 1)
mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)

title(paste0("Parameter estimates for different pool decompositions"), outer = TRUE, 
  cex.main = 2)
dev.off()


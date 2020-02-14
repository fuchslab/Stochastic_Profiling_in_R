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

# p-plots
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

# use two colors in plots FF7F00 #true -> orange color 1F78B4 #estimated -> blue
# color 89BBBF # reference -> turquois color

# Simstudy: Changing Parameters

# plot violinplots for each cellnumber
cellnumber_all <- c("single cells", "2-cell pools", "5-cell pools", "10-cell pools", 
  "15-cell pools", "20-cell pools", "50-cell pools", "mixed pools of single, 2-, 5- and 10-cells", 
  "mixed pools of 10-, 15-, 20-, and 50-cells")
n_short <- c(1, 2, 5, 10, 15, 20, 50, "1-2-5-10", "10-15-20-50")
bw_set <- NULL
ylim_intervall <- NULL

for (i in 1:length(n)) {
  index <- (((i - 1) * 1000 + 1):(i * 1000))
  
  cellnumber <- cellnumber_all[i]
  plotname <- paste("Simstudy2", bquote(.(n_short[i])), "vioplot.pdf", sep = "_")
  
  pdf(plotname, width = 12, height = 7)
  
  par(mfrow = c(3, 4), oma = c(2, 2, 2, 0), mar = c(3, 3, 4, 2), mgp = c(2, 1, 
    0))
  # mtext('Parameter estimates for 10-cell pools', side = 3, line = 0, outer =
  # TRUE)
  if (i > 4) {
    bw_set <- 1
    ylim_intervall <- c(-5, max(mu_2_set2$mu_2[index], mu_2_set1$mu_2[index], 
      mu_2_set3$mu_2[index]))
  }
  
  ## varying p p
  vioplot::vioplot(list(p_set2$p[index], p_set1$p[index], p_set3$p[index]), names = c("", 
    "", ""), col = c("lightgrey", "#89BBBF", "lightgrey"), at = c(1, 2, 4), axes = F, 
    ylab = "", xlab = "", pch = 20, yaxt = "n")
  lines(c(0.6, 4.4), c(0.06, 0.44), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 2, 4), c("p = 0.1", "p = 0.2", "p = 0.4"), cex.axis = 1.1, 
    tick = FALSE)
  mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)
  mtext(expression(paste("Varying ", p)), side = 3, line = 1, las = 1, cex = 1.5)
  
  
  ## mu_1
  vioplot::vioplot(list(mu_1_set2$mu_1[index], mu_1_set1$mu_1[index], mu_1_set3$mu_1[index]), 
    names = c("", "", ""), col = c("lightgrey", "#89BBBF", "lightgrey"), at = c(1, 
      2, 4), axes = F, ylab = "", xlab = "", pch = 20, main = "", yaxt = "n")
  lines(c(0.6, 4.4), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 2, 4), c("p = 0.1", "p = 0.2", "p = 0.4"), cex.axis = 1.1, 
    tick = FALSE)
  mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)
  
  ## mu_2
  vioplot::vioplot(list(mu_2_set2$mu_2[index], mu_2_set1$mu_2[index], mu_2_set3$mu_2[index]), 
    names = c("", "", ""), col = c("lightgrey", "#89BBBF", "lightgrey"), at = c(1, 
      2, 4), axes = F, ylab = "", xlab = "", pch = 20, main = "", yaxt = "n", 
    h = bw_set, ylim = ylim_intervall)
  lines(c(0.6, 4.4), c(0, 0), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 2, 4), c("p = 0.1", "p = 0.2", "p = 0.4"), cex.axis = 1.1, 
    tick = FALSE)
  mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)
  
  
  ## sigmna
  vioplot::vioplot(list(sigma_set2$sigma[index], sigma_set1$sigma[index], sigma_set3$sigma[index]), 
    names = c("", "", ""), col = c("lightgrey", "#89BBBF", "lightgrey"), at = c(1, 
      2, 4), axes = F, ylab = "", xlab = "", pch = 20, main = "", yaxt = "n")
  lines(c(0.6, 4.4), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 2, 4), c("p = 0.1", "p = 0.2", "p = 0.4"), cex.axis = 1.1, 
    tick = FALSE)
  legend(x = 2, y = max(sigma_set2$sigma[index], sigma_set1$sigma[index], sigma_set3$sigma[index]) + 
    0.4 * max(sigma_set2$sigma[index], sigma_set1$sigma[index], sigma_set3$sigma[index]), 
    legend = c("true value", "p = 0.2", expression(paste(mu[1], " =2")), expression(paste(mu[2], 
      " = 0")), expression(paste(sigma, " = 0.2"))), lty = c(2, NA, NA, NA, 
      NA), lw = c(2, NA, NA, NA, 
                   NA), col = c("#FF7F00", NA, NA, NA, NA), bty = "n", xpd = NA, cex = 1.2)
  mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)
  
  ##### varying mu_2
  
  ## p
  vioplot::vioplot(list(p_set1$p[index], p_set4$p[index]), names = c("", ""), col = c("#89BBBF", 
    "lightgrey"), at = c(1, 3), axes = F, ylab = "", xlab = "", pch = 20, yaxt = "n")
  lines(c(0.6, 3.4), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 3), c(expression(paste(mu[2], " = 0")), expression(paste(mu[2], 
    " = 1"))), cex.axis = 1.1, tick = FALSE)
  mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)
  mtext(expression(paste("Varying ", mu[2])), side = 3, line = 1, las = 1, cex = 1.5)
  
  ## mu_1
  vioplot::vioplot(list(mu_1_set1$mu_1[index], mu_1_set4$mu_1[index]), names = c("", 
    ""), col = c("#89BBBF", "lightgrey"), at = c(1, 3), axes = F, ylab = "", 
    xlab = "", pch = 20, main = "", yaxt = "n")
  lines(c(0.6, 3.4), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 3), c(expression(paste(mu[2], " = 0")), expression(paste(mu[2], 
    " = 1"))), cex.axis = 1.1, tick = FALSE)
  mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)
  
  ## mu_2
  vioplot::vioplot(list(mu_2_set1$mu_2[index], mu_2_set4$mu_2[index]), names = c("", 
    ""), col = c("#89BBBF", "lightgrey"), at = c(1, 3), axes = F, ylab = "", 
    xlab = "", pch = 20, main = "", yaxt = "n", h = bw_set, ylim = ylim_intervall)
  lines(c(-1, 5), c(-1, 2), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 3), c(expression(paste(mu[2], " = 0")), expression(paste(mu[2], 
    " = 1"))), cex.axis = 1.1, tick = FALSE)
  mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)
  
  
  ## sigma
  vioplot::vioplot(list(sigma_set1$sigma[index], sigma_set4$sigma[index]), names = c("", 
    ""), col = c("#89BBBF", "lightgrey"), at = c(1, 3), axes = F, ylab = "", 
    xlab = "", pch = 20, main = "", yaxt = "n")
  lines(c(0.6, 3.4), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 3), c(expression(paste(mu[2], " = 0")), expression(paste(mu[2], 
    " = 1"))), cex.axis = 1.1, tick = FALSE)
  mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)
  
  
  ##### varying sigma p
  vioplot::vioplot(list(p_set1$p[index], p_set5$p[index]), names = c("", ""), col = c("#89BBBF", 
    "lightgrey"), at = c(1, 3), axes = F, ylab = "", xlab = "", pch = 20, yaxt = "n")
  lines(c(0.6, 3.4), c(0.2, 0.2), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 3), c(expression(paste(sigma, " = 0.2")), expression(paste(sigma, 
    " = 0.5"))), cex.axis = 1.1, tick = FALSE)
  mtext(expression(p), side = 2, line = 3, las = 1, cex = 1.5)
  mtext(expression(paste("Varying ", sigma)), side = 3, line = 1, las = 1, cex = 1.5)
  
  
  ## mu_1
  vioplot::vioplot(list(mu_1_set1$mu_1[index], mu_1_set5$mu_1[index]), names = c("", 
    ""), col = c("#89BBBF", "lightgrey"), at = c(1, 3), axes = F, ylab = "", 
    xlab = "", pch = 20, main = "", yaxt = "n")
  lines(c(0.6, 3.4), c(2, 2), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 3), c(expression(paste(sigma, " = 0.2")), expression(paste(sigma, 
    " = 0.5"))), cex.axis = 1.1, tick = FALSE)
  mtext(expression(mu[1]), side = 2, line = 3, las = 1, cex = 1.5)
  
  ## mu_2
  vioplot::vioplot(list(mu_2_set1$mu_2[index], mu_2_set5$mu_2[index]), names = c("", 
    ""), col = c("#89BBBF", "lightgrey"), at = c(1, 3), axes = F, ylab = "", 
    xlab = "", pch = 20, main = "", yaxt = "n", h = bw_set, ylim = ylim_intervall)
  lines(c(0.6, 3.4), c(0, 0), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 3), c(expression(paste(sigma, " = 0.2")), expression(paste(sigma, 
    " = 0.5"))), cex.axis = 1.1, tick = FALSE)
  mtext(expression(mu[2]), side = 2, line = 3, las = 1, cex = 1.5)
  
  
  ## sigma
  vioplot::vioplot(list(sigma_set1$sigma[index], sigma_set5$sigma[index]), names = c("", 
    ""), col = c("#89BBBF", "lightgrey"), at = c(1, 3), axes = F, ylab = "", 
    xlab = "", pch = 20, main = "", yaxt = "n")
  lines(c(-1, 3.66666), c(-0.1, 0.6), lw = 2, lty = 2, col = "#FF7F00")
  axis(side = 2, cex.axis = 1.3, las = 1)
  axis(side = 1, at = c(1, 3), c(expression(paste(sigma, " = 0.2")), expression(paste(sigma, 
    " = 0.5"))), cex.axis = 1.1, tick = FALSE)
  mtext(expression(sigma), side = 2, line = 3, las = 1, cex = 1.5)
  
  
  title(paste0("Parameter estimates for ", cellnumber), outer = TRUE, cex.main = 2)
  
  dev.off()
}




# general assumptions:
# Parameters of Population A are based on 50 10-cell samples, 
# In there, Population A got an estimated fraction of p_A = 12%.
# Thus these parameters are based on ca. 60 cells.

# Parameters of Population B are based on 100 10-cell samples, 
# In there, Population B got an estimated fraction of p_B = 20%.
# Thus these parameters are based on ca. 100 cells.


library(MASS)
source("OVL_LN_LN.R")

# 60 und 200
n1 <- 60
n2 <- 200


# 60 und 60
n3 <- 60
n4 <- 60

pdf("Overlap_D.pdf",  width = 12, height = 3)


par(mfrow=c(1,2))

repeat_I <- 1000



# Case D_60_200

mu_A_D <- 2.1
mu_B_D <- 2.03
sigma_A_D <- 0.19
sigma_B_D <- 0.2

OVL_D <- OVL_LN_LN(mu_A_D, mu_B_D, sigma_A_D, sigma_B_D)

# generate new Datasets
A_D_1 <- list()
B_D_1 <- list()

set.seed(940)
for(i in 1:repeat_I){
  A_D_1[[i]] <- rlnorm(n1, meanlog = mean(mu_A_D, mu_B_D), sdlog = mean(sigma_A_D,sigma_B_D))
  B_D_1[[i]] <- rlnorm(n2, meanlog = mean(mu_A_D, mu_B_D), sdlog = mean(sigma_A_D,sigma_B_D))
}

# fit Parameter for each dataset
Para_A_D_1 <- list()
Para_B_D_1 <- list()
set.seed(943)
for(i in 1:repeat_I){
  Para_A_D_1[[i]] <-fitdistr(A_D_1[[i]],"lognormal")$estimate
  Para_B_D_1[[i]] <-fitdistr(B_D_1[[i]],"lognormal")$estimate
}

# Calculate Overlap for each pair
OVL_D_1_all <- c()
for(i in 1:repeat_I){
  OVL_D_1_all[i]<-OVL_LN_LN(Para_A_D_1[[i]][1],Para_B_D_1[[i]][1], Para_A_D_1[[i]][2], Para_B_D_1[[i]][2])
}



Sort_OVL_D_1 <- sort(OVL_D_1_all)
OVL_D_1_95 <- Sort_OVL_D_1[round(length(Sort_OVL_D_1)-length(Sort_OVL_D_1)*0.95)]

hist_D_1 <- hist(OVL_D_1_all,breaks = 100, plot = FALSE)
hist(OVL_D_1_all,breaks = 100, col = c(rep("grey",sum(hist_D_1$breaks<OVL_D_1_95)), rep("black",sum(hist_D_1$breaks>=OVL_D_1_95))), border = c(rep("grey",sum(hist_D_1$breaks<OVL_D_1_95)), rep("black",sum(hist_D_1$breaks>=OVL_D_1_95))),
     xlab = "Overlap", main ="60 and 200 single-cells")

lines(x=rep(OVL_D,2), y=c(0,1000), col = "#89BBBF", lwd= 3)




# Case D_60_60

# generate new Datasets
A_D_2 <- list()
B_D_2 <- list()

set.seed(940)
for(i in 1:repeat_I){
  A_D_2[[i]] <- rlnorm(n3, meanlog = mean(mu_A_D, mu_B_D), sdlog = mean(sigma_A_D,sigma_B_D))
  B_D_2[[i]] <- rlnorm(n4, meanlog = mean(mu_A_D, mu_B_D), sdlog = mean(sigma_A_D,sigma_B_D))
}

# fit Parameter for each dataset
Para_A_D_2 <- list()
Para_B_D_2 <- list()
set.seed(943)
for(i in 1:repeat_I){
  Para_A_D_2[[i]] <-fitdistr(A_D_2[[i]],"lognormal")$estimate
  Para_B_D_2[[i]] <-fitdistr(B_D_2[[i]],"lognormal")$estimate
}

# Calculate Overlap for each pair
OVL_D_2_all <- c()
for(i in 1:repeat_I){
  OVL_D_2_all[i]<-OVL_LN_LN(Para_A_D_2[[i]][1],Para_B_D_2[[i]][1], Para_A_D_2[[i]][2], Para_B_D_2[[i]][2])
}



Sort_OVL_D_2 <- sort(OVL_D_2_all)
OVL_D_2_95 <- Sort_OVL_D_2[round(length(Sort_OVL_D_2)-length(Sort_OVL_D_2)*0.95)]

hist_D_2 <- hist(OVL_D_2_all,breaks = 100, plot = FALSE)
hist(OVL_D_2_all,breaks = 100, col = c(rep("grey",sum(hist_D_2$breaks<OVL_D_2_95)), rep("black",sum(hist_D_2$breaks>=OVL_D_2_95))), border = c(rep("grey",sum(hist_D_2$breaks<OVL_D_2_95)), rep("black",sum(hist_D_2$breaks>=OVL_D_2_95))),
      xlab = "Overlap", main ="60 and 60 single-cells")

lines(x=rep(OVL_D,2), y=c(0,1000), col = "#89BBBF", lwd= 3)


dev.off()

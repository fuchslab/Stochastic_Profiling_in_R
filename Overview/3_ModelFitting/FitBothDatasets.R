# Apply stochprofML on both datasets to estimate the distribution parameters

library(stochprofML)

###########################################################################################
# Dataset A

load( "A_1000.rda")
set.seed(123)
result <- stochprof.loop(model="EXP-LN",dataset=A,n=10,TY=3,genenames=NULL,fix.mu=F,
                         loops=10,until.convergence=F,print.output=F,show.plots=F,
                         plot.title="Entire Cluster",use.constraints=F)

save(result, file="Result_A_1000.rda")


###########################################################################################
# Dataset B
set.seed(10)
n_vektor <- c(5,15,8,sample(2:20,996,replace=TRUE),10)

set.seed(123)
result <- stochprof.loop(model="EXP-LN",dataset=B,n=n_vektor,TY=3,genenames=NULL,fix.mu=F,
                         loops=10,until.convergence=F,print.output=F,show.plots=F,
                         plot.title="Entire Cluster",use.constraints=F)

save(result, file="Result_B_1000.rda")


########################################################################################################################
# n-cell density with Histogram of Data A
# Plot in a Histogram
dens_A <- stochprofML::mix.d.sum.of.mixtures.EXPLN((1:1000)/10, n.vector = 10, p.vector =  p.vector_A, mu.vector =  mu.vector_A, sigma.vector = sigma_A, lambda = lambda_A)

svg(paste0( "histA_fitdens.svg"),  width = 7, height = 8)
hist(A, breaks=seq(-0.5, 100, by = 1), freq = FALSE, col = "black", main = "Sample A", border = "black")
lines(density(A), col=rgb(187,224,227,maxColorValue = 255),lwd=3)
lines((1:1000)/10, dens_A,col= rgb(222,205,135,maxColorValue = 255),lwd=3)
legend("topright", legend =c("data generating pdf", "pdf of fitted model"), lwd=c(3,3), col = c(rgb(187,224,227,maxColorValue = 255), rgb(222,205,135,maxColorValue = 255)))
dev.off()

## n-cell density with Histogram of Data B
# Plot in a Histogram
dens_B <- stochprofML::mix.d.sum.of.mixtures.EXPLN((1:1500)/10, n.vector = n_vektor, p.vector =  p.vector_B, mu.vector =  mu.vector_B, sigma.vector = sigma_B,lambda = lambda_B)

svg(paste0( "histB_fitdens.svg"),  width = 7, height = 8)
hist(B, breaks=seq(-0.5, 150, by = 1), freq = FALSE, col = "black", main = "Sample B", border ="black")
lines(density(B), col=rgb(187,224,227,maxColorValue = 255),lwd=3)
lines((1:1500)/10, dens_B,col=rgb(222,205,135,maxColorValue = 255),lwd=3)
legend("topright", legend =c("data generating pdf", "pdf of fitted model"), lwd=c(3,3), col = c(rgb(187,224,227,maxColorValue = 255), rgb(222,205,135,maxColorValue = 255)))
dev.off()



save(dens_A, dens_B, A, B, file="Dens_fitsA_B.rda")


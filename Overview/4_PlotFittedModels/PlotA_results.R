# PLot the fitted models
###########################################################
# Datset A

load(file="Result_A_1000.rda")

result_A <- result
p.vector_A = c(result_A$mle["p_1"],result_A$mle["p_2"],1-result_A$mle["p_1"]-result_A$mle["p_2"])
mu.vector_A = c(result_A$mle["mu_1_gene_1"], result_A$mle["mu_2_gene_1"])
sigma_A  = rep(result_A$mle["sigma"], 2)
lambda_A= result_A$mle["lambda_1"]

#### 1 Pop Densities A single cell
Pop1_A <- dlnorm((1:10000)/100, meanlog =  mu.vector_A[1], sdlog = sigma_A[1])
Pop2_A <- dlnorm((1:10000)/100, meanlog =  mu.vector_A[2], sdlog = sigma_A[2])
Pop3_A <- dexp((1:10000)/100,rate = lambda_A)

svg(paste0( "Pop1A_dens.svg"),  width = 7, height = 8)
plot(((1:10000)/100)[1:3000],Pop1_A[1:3000], lwd = 3, type = "l", col =rgb(red = 137, green = 184, blue = 167, max = 255))
polygon(((1:10000)/100)[1:3000],Pop1_A[1:3000],col=rgb(red = 187, green = 224, blue = 227, max = 255), border = NA)
lines(((1:10000)/100)[1:3000],Pop1_A[1:3000], lwd = 3, col =rgb(red = 137, green = 184, blue = 167, max = 255))
dev.off()

svg(paste0( "Pop2A_dens.svg"),  width = 7, height = 8)
plot(((1:10000)/100)[1:1000],Pop2_A[1:1000], lwd = 3, type = "l", col =rgb(red = 45, green = 45, blue = 138, max = 255))
polygon(((1:10000)/100)[1:1000],Pop2_A[1:1000],col=rgb(red = 156, green = 156, blue = 223, max = 255), border = NA)
lines(((1:10000)/100)[1:1000],Pop2_A[1:1000], lwd = 3,  col =rgb(red = 45, green = 45, blue = 138, max = 255))
dev.off()

svg(paste0( "Pop3A_dens.svg"),  width = 7, height = 8)
plot(((1:10000)/100)[1:200],Pop3_A[1:200],  lwd = 3, type = "l", col =rgb(red = 27, green = 27, blue = 27, max = 255))
polygon(c(0,((1:10000)/100)[1:200]),c(0,Pop3_A[1:200]), col=rgb(red = 128, green = 128, blue = 128, max = 255), border = NA)
lines(((1:10000)/100)[1:200],Pop3_A[1:200],  lwd = 3, col =rgb(red = 27, green = 27, blue = 27, max = 255))
dev.off()

# Mix Density A single cell

Pop123_A <- stochprofML::d.sum.of.mixtures.EXPLN((1:2000)/100, n = 1, p.vector=p.vector_A, mu.vector =  mu.vector_A, sigma  = sigma_A, lambda= lambda_A, logden = F)

svg(paste0( "Pop123mixA_dens.svg"),  width = 7, height = 8)
plot(((1:2000)/100),Pop123_A,  lwd = 3, type = "l", col = "black")
polygon(c(0,((1:60)/100), 60/100),c(0,Pop123_A[1:60], 0), col = rgb(red = 128, green = 128, blue = 128, max = 255), border = NA)
polygon(c(60/100,((60:420)/100), 420/100),c(0,Pop123_A[60:420], 0), col = rgb(red = 156, green = 156, blue = 223, max = 255), border = NA)
polygon(c(60/100,(60:120)/100, 120/100),c(0,Pop123_A[60:120],0), col =rgb(red = 128, green = 128, blue = 128, max = 255), angle=45, dens = 50,border = NA, lwd = 3)
polygon(c(420/100,((420:2000)/100), 2000/100),c(0,Pop123_A[420:2000], 0), col = rgb(red = 187, green = 224, blue = 227, max = 255), border = NA)
polygon(c(420/100,(420:550)/100, 550/100),c(0,Pop123_A[420:550],0), col =rgb(red = 156, green = 156, blue = 223, max = 255), angle=45, dens = 50,border = NA, lwd = 3)
lines(((1:2000)/100),Pop123_A,  lwd = 3, col = "black")
dev.off()

# save densities A
save(Pop1_A, Pop2_A, Pop3_A, Pop123_A, file = "PopA_dens.rda")


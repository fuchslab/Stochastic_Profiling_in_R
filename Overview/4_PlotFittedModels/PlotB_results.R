# PLot the fitted models
###########################################################
# Datset B

load(file="Result_B_1000.rda")

result_B <- result

p.vector_B = c(result_B$mle["p_1"],result_B$mle["p_2"],1-result_B$mle["p_1"]-result_B$mle["p_2"])

mu.vector_B = c(result_B$mle["mu_1_gene_1"], result_B$mle["mu_2_gene_1"])
sigma_B  = rep(result_B$mle["sigma"], 2)
lambda_B= result_B$mle["lambda_1"]

# 1 Pop Densities A single cell
Pop1_B <- dlnorm((1:10000)/100, meanlog =  mu.vector_B[1], sdlog = sigma_B[1])
Pop2_B <- dlnorm((1:10000)/100, meanlog =  mu.vector_B[2], sdlog = sigma_B[2])
Pop3_B <- dexp((1:10000)/100,rate = lambda_B)

svg(paste0( "Pop1B_dens.svg"),  width = 7, height = 8)
plot(((1:10000)/100)[1:3000],Pop1_B[1:3000], lwd = 3, type = "l", col =rgb(red = 137, green = 184, blue = 167, max = 255))
polygon(((1:10000)/100)[1:3000],Pop1_B[1:3000],col=rgb(red = 187, green = 224, blue = 227, max = 255), border = NA)
lines(((1:10000)/100)[1:3000],Pop1_B[1:3000], lwd = 3, col =rgb(red = 137, green = 184, blue = 167, max = 255))
dev.off()

svg(paste0( "Pop2B_dens.svg"),  width = 7, height = 8)
plot(((1:10000)/100)[1:1000],Pop2_B[1:1000], lwd = 3, type = "l", col =rgb(red = 45, green = 45, blue = 138, max = 255))
polygon(((1:10000)/100)[1:1000],Pop2_B[1:1000],col=rgb(red = 156, green = 156, blue = 223, max = 255), border = NA)
lines(((1:10000)/100)[1:1000],Pop2_B[1:1000], lwd = 3,  col =rgb(red = 45, green = 45, blue = 138, max = 255))
dev.off()

svg(paste0( "Pop3B_dens.svg"),  width = 7, height = 8)
plot(((1:10000)/100)[1:200],Pop3_B[1:200],  lwd = 3, type = "l", col =rgb(red = 27, green = 27, blue = 27, max = 255))
polygon(c(0,((1:10000)/100)[1:200]),c(0,Pop3_B[1:200]), col=rgb(red = 128, green = 128, blue = 128, max = 255), border = NA)
lines(((1:10000)/100)[1:200],Pop3_B[1:200],  lwd = 3, col =rgb(red = 27, green = 27, blue = 27, max = 255))
dev.off()

# Mix Density B single cell

Pop123_B <- stochprofML::d.sum.of.mixtures.EXPLN((1:2000)/100, n = 1, p.vector=p.vector_B, mu.vector =  mu.vector_B, sigma  = sigma_B, lambda= lambda_B, logden = F)

svg(paste0( "Pop123mixB_dens.svg"),  width = 7, height = 8)
plot(((1:2000)/100),Pop123_B,  lwd = 3, type = "l", col = "black")

polygon(c(0,((1:60)/100), 60/100),c(0,Pop123_B[1:60], 0), col = rgb(red = 128, green = 128, blue = 128, max = 255), border = NA)

polygon(c(60/100,((60:420)/100), 420/100),c(0,Pop123_B[60:420], 0), col = rgb(red = 156, green = 156, blue = 223, max = 255), border = NA)
polygon(c(60/100,(60:120)/100, 120/100),c(0,Pop123_B[60:120],0), col =rgb(red = 128, green = 128, blue = 128, max = 255), angle=45, dens = 50,border = NA, lwd = 3)


polygon(c(420/100,((420:2000)/100), 2000/100),c(0,Pop123_B[420:2000], 0), col = rgb(red = 187, green = 224, blue = 227, max = 255), border = NA)
polygon(c(420/100,(420:550)/100, 550/100),c(0,Pop123_B[420:550],0), col =rgb(red = 156, green = 156, blue = 223, max = 255), angle=45, dens = 50,border = NA, lwd = 3)

lines(((1:2000)/100),Pop123_B,  lwd = 3, col = "black")
dev.off()

### save densities B
save(Pop1_B, Pop2_B, Pop3_B, Pop123_B, file = "PopB_dens.rda")




# Plot true distributions, model assumpiton in stochprofML

p.vector = c(0.3,0.6,0.1)

mu.vector = c(2, 0.8)
sigma  = c(0.4, 0.4)
lambda= 3


Pop1 <- dlnorm((1:10000)/100, meanlog =  mu.vector[1], sdlog = sigma[1])
Pop2 <- dlnorm((1:10000)/100, meanlog =  mu.vector[2], sdlog = sigma[2])
Pop3 <- dexp((1:10000)/100,rate = lambda)

svg(paste0( "Pop1_dens.svg"),  width = 7, height = 8)
plot(((1:10000)/100)[1:3000],Pop1[1:3000], lwd = 3, type = "l", col =rgb(red = 137, green = 184, blue = 167, max = 255))
polygon(((1:10000)/100)[1:3000],Pop1[1:3000],col=rgb(red = 187, green = 224, blue = 227, max = 255), border = NA)
lines(((1:10000)/100)[1:3000],Pop1[1:3000], lwd = 3, col =rgb(red = 137, green = 184, blue = 167, max = 255))
dev.off()

svg(paste0( "Pop2_dens.svg"),  width = 7, height = 8)
plot(((1:10000)/100)[1:1000],Pop2[1:1000], lwd = 3, type = "l", col =rgb(red = 45, green = 45, blue = 138, max = 255))
polygon(((1:10000)/100)[1:1000],Pop2[1:1000],col=rgb(red = 156, green = 156, blue = 223, max = 255), border = NA)
lines(((1:10000)/100)[1:1000],Pop2[1:1000], lwd = 3,  col =rgb(red = 45, green = 45, blue = 138, max = 255))
dev.off()

svg(paste0( "Pop3_dens.svg"),  width = 7, height = 8)
plot(((1:10000)/100)[1:200],Pop3[1:200],  lwd = 3, type = "l", col =rgb(red = 27, green = 27, blue = 27, max = 255))
polygon(c(0,((1:10000)/100)[1:200]),c(0,Pop3[1:200]), col=rgb(red = 128, green = 128, blue = 128, max = 255), border = NA)
lines(((1:10000)/100)[1:200],Pop3[1:200],  lwd = 3, col =rgb(red = 27, green = 27, blue = 27, max = 255))
dev.off()

### Mix

Pop123 <- stochprofML::d.sum.of.mixtures.EXPLN((1:2000)/100, n = 1, p.vector=p.vector, mu.vector =  mu.vector, sigma  = sigma, lambda= lambda, logden = F)

svg(paste0( "Pop123mix_dens.svg"),  width = 7, height = 8)
plot(((1:2000)/100),Pop123,  lwd = 3, type = "l", col = "black")
polygon(c(0,((1:60)/100), 60/100),c(0,Pop123[1:60], 0), col = rgb(red = 128, green = 128, blue = 128, max = 255), border = NA)
polygon(c(60/100,((60:420)/100), 420/100),c(0,Pop123[60:420], 0), col = rgb(red = 156, green = 156, blue = 223, max = 255), border = NA)
polygon(c(60/100,(60:120)/100, 120/100),c(0,Pop123[60:120],0), col =rgb(red = 128, green = 128, blue = 128, max = 255), angle=45, dens = 50,border = NA, lwd = 3)
polygon(c(420/100,((420:2000)/100), 2000/100),c(0,Pop123[420:2000], 0), col = rgb(red = 187, green = 224, blue = 227, max = 255), border = NA)
polygon(c(420/100,(420:550)/100, 550/100),c(0,Pop123[420:550],0), col =rgb(red = 156, green = 156, blue = 223, max = 255), angle=45, dens = 50,border = NA, lwd = 3)
lines(((1:2000)/100),Pop123,  lwd = 3, col = "black")
dev.off()


save(Pop1, Pop2, Pop3, Pop123, file = "Orig_Pop_dens.rda")

# First Try: Different distributions


mu_1 <- 2
mu_2 <- 2.7
sigma_1 <- 0.25
sigma_2 <- 0.2
x <- seq(0.1, 40, by = 0.1)

Distribution1 <- dlnorm(x, meanlog = mu_1, sdlog = sigma_1)
Distribution2 <- dlnorm(x, meanlog = mu_2, sdlog = sigma_2)
Overlap_1_2 <- OVL_LN_LN(mu_1, mu_2, sigma_1, sigma_2)

pdf("Example_Overlap.pdf", width = 12, height = 5)
par(mfrow = c(2, 2))
plot(x, Distribution1, type = "l", col = "#FF7F00", lwd = 3, main = paste("Overlap area:", round(bquote(.(Overlap_1_2))*100, 0), "%"), xlab = "Measurements", ylab = "Density", bty = "n")
lines(x, Distribution2, col = "#1F78B4", lwd = 3)
polygon(x, c(pmin(Distribution1, Distribution2)), col = "#FF7F00", border = NA)
polygon(x, c(pmin(Distribution1, Distribution2)), col = "#1F78B4", angle = 45, dens = 50, border = NA, lwd = 1)
lines(x, Distribution1, col = "#FF7F00", lwd = 3)
lines(x, Distribution2, col = "#1F78B4", lwd = 3)
legend("topright", c( expression(paste("LN(", mu[1], " = 2.00, ", sigma[1], " = 0.25)")), expression(paste("LN(", mu[2], " = 2.70, ", sigma[2], " = 0.20)")) ), col = c("#FF7F00", "#1F78B4"), lwd = c(3, 3), bty = "n")
fig_label("A", cex = 2) 
#function fig_label() at the end of this document

mu_3 <- 1.75
mu_4 <- 2.09
sigma_3 <- 0.19
sigma_4 <- 0.21
x <- seq(0.1, 40, by = 0.1)
Distribution3 <- dlnorm(x, meanlog = mu_3, sdlog = sigma_3)
Distribution4 <- dlnorm(x, meanlog = mu_4, sdlog = sigma_4)
Overlap_3_4 <- OVL_LN_LN(mu_3, mu_4, sigma_3, sigma_4)
plot(x, Distribution3, type = "l", col = "#FF7F00", lwd = 3, main = paste("Overlap area:", round(bquote(.(Overlap_3_4))*100, 0), "%"), xlab = "Measurements", ylab = "Density", bty = "n")
lines(x, Distribution4, col = "#1F78B4", lwd = 3)
polygon(x, c(pmin(Distribution3, Distribution4)), col = "#FF7F00", border = NA)
polygon(x, c(pmin(Distribution3, Distribution4)), col = "#1F78B4", angle = 45, dens = 50, border = NA, lwd = 1)
lines(x, Distribution3, col = "#FF7F00", lwd = 3)
lines(x, Distribution4, col = "#1F78B4", lwd = 3)
legend("topright", c( expression(paste("LN(", mu[1], " = 1.75, ", sigma[1], " = 0.19)")), expression(paste("LN(", mu[2], " = 2.09, ", sigma[2], " = 0.21)")) ), col = c("#FF7F00", "#1F78B4"), lwd = c(3, 3), bty = "n")
fig_label("B", cex = 2) 

mu_5 <- 1.89
mu_6 <- 2.07
sigma_5 <- 0.21
sigma_6 <- 0.25
x <- seq(0.1, 40, by = 0.1)
Distribution5 <- dlnorm(x, meanlog = mu_5, sdlog = sigma_5)
Distribution6 <- dlnorm(x, meanlog = mu_6, sdlog = sigma_6)
Overlap_5_6 <- OVL_LN_LN(mu_5, mu_6, sigma_5, sigma_6)
plot(x, Distribution5, type = "l", col = "#FF7F00", lwd = 3, main = paste("Overlap area:", round(bquote(.(Overlap_5_6))*100, 0), "%"), xlab = "Measurements", ylab = "Density", bty = "n")
lines(x, Distribution6, col = "#1F78B4", lwd = 3)
polygon(x, c(pmin(Distribution5, Distribution6)), col = "#FF7F00", border = NA)
polygon(x, c(pmin(Distribution5, Distribution6)), col = "#1F78B4", angle = 45, dens = 50, border = NA, lwd = 1)
lines(x, Distribution5, col = "#FF7F00", lwd = 3)
lines(x, Distribution6, col = "#1F78B4", lwd = 3)
legend("topright", c( expression(paste("LN(", mu[1], " = 1.89, ", sigma[1], " = 0.21)")), expression(paste("LN(", mu[2], " = 2.07, ", sigma[2], " = 0.25)")) ), col = c("#FF7F00", "#1F78B4"), lwd = c(3, 3), bty = "n")
fig_label("C", cex = 2) 

mu_7 <- 2.1
mu_8 <- 2.03
sigma_7 <- 0.19
sigma_8 <- 0.2
x <- seq(0.1, 40, by = 0.1)
Distribution7 <- dlnorm(x, meanlog = mu_7, sdlog = sigma_7)
Distribution8 <- dlnorm(x, meanlog = mu_8, sdlog = sigma_8)
Overlap_7_8 <- OVL_LN_LN(mu_7, mu_8, sigma_7, sigma_8)


plot(x, Distribution7, type = "l", col = "#FF7F00", lwd = 3, main = paste("Overlap area:", round(bquote(.(Overlap_7_8))*100, 0), "%"), xlab = "Measurements", ylab = "Density", bty = "n")
lines(x, Distribution8, col = "#1F78B4", lwd = 3)
polygon(x, c(pmin(Distribution7, Distribution8)), col = "#FF7F00", border = NA)
polygon(x, c(pmin(Distribution7, Distribution8)), col = "#1F78B4", angle = 45, dens = 50, border = NA, lwd = 1)
lines(x, Distribution7, col = "#FF7F00", lwd = 3)
lines(x, Distribution8, col = "#1F78B4", lwd = 3)
legend("topright", c( expression(paste("LN(", mu[1], " = 2.10, ", sigma[1], " = 0.19)")), expression(paste("LN(", mu[2], " = 2.03, ", sigma[2], " = 0.20)")) ), col = c("#FF7F00", "#1F78B4"), lwd = c(3, 3), bty = "n")
fig_label("D", cex = 2) 

dev.off()



###############################################################################################################################################################################
# function found at: https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}



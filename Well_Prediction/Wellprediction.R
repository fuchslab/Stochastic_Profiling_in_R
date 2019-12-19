##################################################
# load packages
library("stochprofML")
library(ggplot2)
library(gridExtra)
library(dplyr)
library(knitr)
library(RColorBrewer)
library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

####################################
# Functions
hpi <- function(probs_distr, j.vec.all , hpi_p){
  pos <- which.max(probs_distr)   
  MLE <- j.vec.all[pos,]
  MLE_v <- probs_distr[pos]
  hpi_v <- MLE_v
  pos_l <- pos
  pos_u <- pos
  
  while(hpi_v < hpi_p){
    if(pos_l > 1 & pos_u < length(probs_distr)){
      pos_l_temp <- pos_l-1
      pos_u_temp <- pos_u + 1
      prob_l_temp <- probs_distr[pos_l_temp]
      prob_u_temp <- probs_distr[pos_u_temp]
      decision <- which.max(c(prob_l_temp, prob_u_temp))
      if(decision == 1) pos_l <- pos_l_temp else pos_u <- pos_u_temp
      if(decision == 1) hpi_v <-hpi_v + prob_l_temp else hpi_v <-hpi_v + prob_u_temp
    } else if(pos_l > 1 & pos_u == length(probs_distr)){
      pos_l <- pos_l-1
      hpi_v <-hpi_v + probs_distr[pos_l]
    } else if( pos_u < length(probs_distr)){
      pos_u <- pos_u+1
      hpi_v <-hpi_v + probs_distr[pos_u]
    }
    
  }
  
  return(list("HPI"= hpi_v, "MLE"=MLE, "HPI_l"= j.vec.all[pos_l,], "HPI_u"= j.vec.all[pos_u,]))
  
}


###################################################################################################
p.combo.sum.of.mix <- function(y,n,p.vector.est=FALSE,mu.vector.est,sigma.vector.est){
  p.vector.est.is <- p.vector.est
  p.vector.est.is <- is.numeric(p.vector.est.is)
  j.combis <- stochprofML:::comb.summands(n,length(mu.vector.est))
  
  Dens_j <- matrix(ncol=dim(j.combis)[1], nrow=length(y))
  
  #p.vector not needed
  if(p.vector.est.is==FALSE){
    for(i in 1:nrow(j.combis)) {
      Dens_j[,i] <-  stochprofML:::d.sum.of.types.LNLN(y,j.combis[i,],mu.vector.est,sigma.vector.est,logdens=F)
    }
    Prob_j <- Dens_j/rowSums(Dens_j)
  }else{
    for (i in 1:nrow(j.combis))  {
      weight <- dmultinom(x=j.combis[i,],prob=p.vector.est,log=F)
      mixture.density <- stochprofML:::d.sum.of.types.LNLN(y,j.combis[i,],mu.vector.est,sigma.vector.est,logdens=F)
      Dens_j[,i] <- weight * mixture.density
    }
    Prob_j <- Dens_j/rowSums(Dens_j)
  }
  
  return(list("Prob_j"=Prob_j,"Dens_j"= Dens_j))
}

#########################################################################


p.vector <- c(0.2,0.8)
mu.vector <- c(2,0)
sigma.vector <- c(0.2,0.2)
k <- 100 #number of measurement samples
n <- 5 #the number of cells taken from each measurement sample
gene_nr <- 3
set_seeds <- c(1234,345)


#########################################################################
##############################################################################################
# start general part
library(stochprofML)
model="LN-LN"

###############################################################################################

set_a_seed <- set_seeds[1]
set.seed(set_a_seed)
N.matrix <- t(rmultinom(n = k, size = n, prob = p.vector))
###############################################################################################
## Generate Datasets (3 genes)
###############################################################################################

Dataset <- matrix(ncol=k, nrow = gene_nr)
for(i in 1:gene_nr){
  set.seed(set_a_seed+100+i)
  Dataset[i,] <- r.sum.of.mixtures.LNLN(k, n, p.vector, mu.vector, sigma.vector, N.matrix)
}

#Estimate Parameters
estimated.parameter <- list()
TY <- length(p.vector)
set_a_seed <- set_seeds[2]
for(i in 1:gene_nr){
  set.seed(set_a_seed+100+i)
  estimated.parameter[[i]] <- stochprof.loop(model=model,dataset=t(Dataset[i,,drop =FALSE]),n=n,TY=TY,print.output=FALSE,show.plots=FALSE)
}
###############################################################################################
# Calculate all the densities for each gene fit
###############################################################################################

j.vector.all<-stochprofML:::comb.summands(n=n, k=TY)
dens.all <- list()
for(i in 1:gene_nr){
  x= seq(round(min(Dataset[i,])),round(max(Dataset[i,])),(round(max(Dataset[i,]))-round(min(Dataset[i,])))/99)
  dens.all[[i]] <- matrix(nrow=2*dim(j.vector.all)[1]+3,ncol = length(x))
  dens.all[[i]][1,] <- x
  dens.all[[i]][2,] <- stochprofML:::d.sum.of.mixtures.LNLN(x,n,c(estimated.parameter[[i]]$mle[1],(1-estimated.parameter[[i]]$mle[1])),c(estimated.parameter[[i]]$mle[2],estimated.parameter[[i]]$mle[3]),rep(estimated.parameter[[i]]$mle[4],TY),logdens = F)
  for (j in 1:nrow(j.vector.all)){
    dens.all[[i]][j+2,]<- stochprofML:::d.sum.of.types.LNLN(x, j.vector=j.vector.all[j,], mu.vector=c(estimated.parameter[[i]]$mle[2],estimated.parameter[[i]]$mle[3]), sigma.vector=rep(estimated.parameter[[i]]$mle[4],TY), logdens = FALSE)
  }
  dens.all[[i]][nrow(j.vector.all)+3,] <- stochprofML:::d.sum.of.mixtures.LNLN(x,n,p.vector,mu.vector,sigma.vector,logdens = F)
  for (j in 1:nrow(j.vector.all)){
    dens.all[[i]][j+nrow(j.vector.all)+3,]<- stochprofML:::d.sum.of.types.LNLN(x, j.vector=j.vector.all[j,], mu.vector=mu.vector, sigma.vector=sigma.vector, logdens = FALSE)
  }
}

###############################################################################################
# Calculate all conditional densities for all possible well compositions and the Porbability,
# either with or without using the predicted p
###############################################################################################
cond.dens.prob.all_nop <- list()
cond.dens.prob.all_wp <- list()
cond.dens.prob.all.orig_nop <- list()
cond.dens.prob.all.orig_wp <- list()
cond.prob.allGenes_nop <- matrix(0,ncol = nrow(j.vector.all), nrow= length(Dataset[i,]))
cond.prob.allGenes_wp <- matrix(0,ncol = nrow(j.vector.all), nrow= length(Dataset[i,]))
cond.prob.allGenes.orig_nop <- matrix(0,ncol = nrow(j.vector.all), nrow= length(Dataset[i,]))
cond.prob.allGenes.orig_wp <- matrix(0,ncol = nrow(j.vector.all), nrow= length(Dataset[i,]))
for(i in 1:gene_nr){
  for (j in 1:nrow(j.vector.all)){
    cond.dens.prob.all_nop[[i]]<-p.combo.sum.of.mix(y=t(Dataset[i,,drop =FALSE]),n,p.vector.est=FALSE,mu.vector.est=c(estimated.parameter[[i]]$mle[2],estimated.parameter[[i]]$mle[3]),sigma.vector.est=rep(estimated.parameter[[i]]$mle[4],TY))
    cond.dens.prob.all_wp[[i]]<-p.combo.sum.of.mix(y=t(Dataset[i,,drop =FALSE]),n,p.vector.est=c(estimated.parameter[[i]]$mle[1],(1-estimated.parameter[[i]]$mle[1])),mu.vector.est=c(estimated.parameter[[i]]$mle[2],estimated.parameter[[i]]$mle[3]),sigma.vector.est=rep(estimated.parameter[[i]]$mle[4],TY))
    cond.dens.prob.all.orig_nop[[i]]<-p.combo.sum.of.mix(y=t(Dataset[i,,drop =FALSE]),n,p.vector.est=FALSE,mu.vector.est=mu.vector,sigma.vector.est=sigma.vector)
    cond.dens.prob.all.orig_wp[[i]]<-p.combo.sum.of.mix(y=t(Dataset[i,,drop =FALSE]),n,p.vector.est=p.vector,mu.vector.est=mu.vector,sigma.vector.est=sigma.vector)
  }
  cond.prob.allGenes_nop <- cond.prob.allGenes_nop + cond.dens.prob.all_nop[[i]]$Prob_j
  cond.prob.allGenes_wp <- cond.prob.allGenes_wp + cond.dens.prob.all_wp[[i]]$Prob_j
  cond.prob.allGenes.orig_nop <- cond.prob.allGenes.orig_nop + cond.dens.prob.all.orig_nop[[i]]$Prob_j
  cond.prob.allGenes.orig_wp <- cond.prob.allGenes.orig_wp + cond.dens.prob.all.orig_wp[[i]]$Prob_j
}
cond.prob.allGenes_nop <- cond.prob.allGenes_nop/gene_nr
cond.prob.allGenes_wp <- cond.prob.allGenes_wp/gene_nr
cond.prob.allGenes.orig_nop <- cond.prob.allGenes.orig_nop/gene_nr
cond.prob.allGenes.orig_wp <- cond.prob.allGenes.orig_wp/gene_nr

###############################################################################################
# Predict Cell Combos
###############################################################################################
# HPI done here only works properly for 2Pop

Predic_nop <- list()
Predic_wp <- list()
MLE_nop <- list()
MLE_wp <- list()
HPI_v_nop <- list()
HPI_v_wp <- list()
HPI_l_nop <- list()
HPI_l_wp <- list()
HPI_u_nop <- list()
HPI_u_wp <- list()
for(i in 1:gene_nr){
  Predic_nop[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  Predic_wp[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  MLE_nop[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI_l_nop[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI_u_nop[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI_v_nop[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = 1)
  MLE_wp[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI_l_wp[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI_u_wp[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI_v_wp[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = 1)
  for (j in 1:ncol(Dataset[i,,drop =FALSE])){
    Predic_nop[[i]][j,] <- colSums(j.vector.all * cond.dens.prob.all_nop[[i]]$Prob_j[j,])
    Predic_wp[[i]][j,] <- colSums(j.vector.all * cond.dens.prob.all_wp[[i]]$Prob_j[j,])
    mle_hpi95_nop <- hpi(probs_distr = cond.dens.prob.all_nop[[i]]$Prob_j[j,] , j.vec.all = j.vector.all, hpi_p = 0.95)
    MLE_nop[[i]][j,] <- mle_hpi95_nop$MLE
    HPI_l_nop[[i]][j,] <- mle_hpi95_nop$HPI_l
    HPI_u_nop[[i]][j,] <- mle_hpi95_nop$HPI_u
    HPI_v_nop[[i]][j,] <-  mle_hpi95_nop$HPI
    mle_hpi95_wp <- hpi(probs_distr = cond.dens.prob.all_wp[[i]]$Prob_j[j,] , j.vec.all = j.vector.all, hpi_p = 0.95)
    MLE_wp[[i]][j,] <- j.vector.all[which.max(cond.dens.prob.all_wp[[i]]$Prob_j[j,]), ]
    HPI_l_wp[[i]][j,] <- mle_hpi95_wp$HPI_l
    HPI_u_wp[[i]][j,] <- mle_hpi95_wp$HPI_u
    HPI_v_wp[[i]][j,] <-  mle_hpi95_wp$HPI
    
  }
}
# same for original parameters
Predic.orig_nop <- list()
Predic.orig_wp <- list()
MLE.orig_nop <- list()
MLE.orig_wp <- list()
HPI.orig_v_nop <- list()
HPI.orig_v_wp <- list()
HPI.orig_l_nop <- list()
HPI.orig_l_wp <- list()
HPI.orig_u_nop <- list()
HPI.orig_u_wp <- list()
for(i in 1:gene_nr){
  Predic.orig_nop[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  Predic.orig_wp[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  MLE.orig_nop[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI.orig_l_nop[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI.orig_u_nop[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI.orig_v_nop[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = 1)
  MLE.orig_wp[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI.orig_l_wp[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI.orig_u_wp[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI.orig_v_wp[[i]] <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = 1)
  for (j in 1:ncol(Dataset[i,,drop =FALSE])){
    Predic.orig_nop[[i]][j,] <- colSums(j.vector.all * cond.dens.prob.all.orig_nop[[i]]$Prob_j[j,])
    Predic.orig_wp[[i]][j,] <- colSums(j.vector.all * cond.dens.prob.all.orig_wp[[i]]$Prob_j[j,])
    mle_hpi95.orig_nop <- hpi(probs_distr = cond.dens.prob.all.orig_nop[[i]]$Prob_j[j,] , j.vec.all = j.vector.all, hpi_p = 0.95)
    MLE.orig_nop[[i]][j,] <- mle_hpi95.orig_nop$MLE
    HPI.orig_l_nop[[i]][j,] <- mle_hpi95.orig_nop$HPI_l
    HPI.orig_u_nop[[i]][j,] <- mle_hpi95.orig_nop$HPI_u
    HPI.orig_v_nop[[i]][j,] <-  mle_hpi95.orig_nop$HPI
    mle_hpi95.orig_wp <- hpi(probs_distr = cond.dens.prob.all.orig_wp[[i]]$Prob_j[j,] , j.vec.all = j.vector.all, hpi_p = 0.95)
    MLE.orig_wp[[i]][j,] <- j.vector.all[which.max(cond.dens.prob.all.orig_wp[[i]]$Prob_j[j,]), ]
    HPI.orig_l_wp[[i]][j,] <- mle_hpi95.orig_wp$HPI_l
    HPI.orig_u_wp[[i]][j,] <- mle_hpi95.orig_wp$HPI_u
    HPI.orig_v_wp[[i]][j,] <-  mle_hpi95.orig_wp$HPI
    
  }
}

#same for all genes combined
Predic.allGenes_nop <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
Predic.allGenes_wp <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
MLE.allGenes_nop <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_l.allGenes_nop <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_u.allGenes_nop <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_v.allGenes_nop <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = 1)
MLE.allGenes_wp <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_l.allGenes_wp <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_u.allGenes_wp <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_v.allGenes_wp <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = 1)
for (j in 1:ncol(Dataset[1,,drop =FALSE])){
  Predic.allGenes_nop[j,] <- colSums(j.vector.all * cond.prob.allGenes_nop[j,])
  Predic.allGenes_wp[j,] <- colSums(j.vector.all * cond.prob.allGenes_wp[j,])
  mle_hpi95.allGenes_nop <- hpi(probs_distr = cond.prob.allGenes_nop[j,] , j.vec.all = j.vector.all, hpi_p = 0.95)
  MLE.allGenes_nop[j,] <- mle_hpi95.allGenes_nop$MLE
  HPI_l.allGenes_nop[j,] <- mle_hpi95.allGenes_nop$HPI_l
  HPI_u.allGenes_nop[j,] <- mle_hpi95.allGenes_nop$HPI_u
  HPI_v.allGenes_nop[j,] <-  mle_hpi95.allGenes_nop$HPI
  mle_hpi95.allGenes_wp <- hpi(probs_distr = cond.prob.allGenes_wp[j,] , j.vec.all = j.vector.all, hpi_p = 0.95)
  MLE.allGenes_wp[j,] <- j.vector.all[which.max(cond.prob.allGenes_wp[j,]), ]
  HPI_l.allGenes_wp[j,] <- mle_hpi95.allGenes_wp$HPI_l
  HPI_u.allGenes_wp[j,] <- mle_hpi95.allGenes_wp$HPI_u
  HPI_v.allGenes_wp[j,] <-  mle_hpi95.allGenes_wp$HPI
  
}

#same for all genes combined of the original parameters
Predic.allGenes.orig_nop <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
Predic.allGenes.orig_wp <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
MLE.allGenes.orig_nop <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_l.allGenes.orig_nop <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_u.allGenes.orig_nop <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_v.allGenes.orig_nop <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = 1)
MLE.allGenes.orig_wp <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_l.allGenes.orig_wp <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_u.allGenes.orig_wp <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_v.allGenes.orig_wp <- matrix(nrow = ncol(Dataset[i,,drop =FALSE]), ncol = 1)
for (j in 1:ncol(Dataset[1,,drop =FALSE])){
  Predic.allGenes.orig_nop[j,] <- colSums(j.vector.all * cond.prob.allGenes.orig_nop[j,])
  Predic.allGenes.orig_wp[j,] <- colSums(j.vector.all * cond.prob.allGenes.orig_wp[j,])
  mle_hpi95.allGenes.orig_nop <- hpi(probs_distr = cond.prob.allGenes.orig_nop[j,] , j.vec.all = j.vector.all, hpi_p = 0.95)
  MLE.allGenes.orig_nop[j,] <- mle_hpi95.allGenes.orig_nop$MLE
  HPI_l.allGenes.orig_nop[j,] <- mle_hpi95.allGenes.orig_nop$HPI_l
  HPI_u.allGenes.orig_nop[j,] <- mle_hpi95.allGenes.orig_nop$HPI_u
  HPI_v.allGenes.orig_nop[j,] <-  mle_hpi95.allGenes.orig_nop$HPI
  mle_hpi95.allGenes.orig_wp <- hpi(probs_distr = cond.prob.allGenes.orig_wp[j,] , j.vec.all = j.vector.all, hpi_p = 0.95)
  MLE.allGenes.orig_wp[j,] <- j.vector.all[which.max(cond.prob.allGenes.orig_wp[j,]), ]
  HPI_l.allGenes.orig_wp[j,] <- mle_hpi95.allGenes.orig_wp$HPI_l
  HPI_u.allGenes.orig_wp[j,] <- mle_hpi95.allGenes.orig_wp$HPI_u
  HPI_v.allGenes.orig_wp[j,] <-  mle_hpi95.allGenes.orig_wp$HPI
  
}


# Check how often the true combo gets predicted
# Check how often the true combo gets predicted
True.predictions_nop <- c()
True.predictions_wp <- c()
True.predictions.orig_nop <- c()
True.predictions.orig_wp <- c()
True.predictions.allGenes_nop <- c()
True.predictions.allGenes_wp <- c()
True.predictions.allGenes.orig_nop <- c()
True.predictions.allGenes.orig_wp <- c()

for(i in 1:gene_nr){
  True.predictions_nop[i] <- sum(round(Predic_nop[[i]],0)== N.matrix)/ncol(N.matrix)
  True.predictions_wp[i] <- sum(round(Predic_wp[[i]],0)== N.matrix)/ncol(N.matrix)
  True.predictions.orig_nop[i] <- sum(round(Predic.orig_nop[[i]],0)== N.matrix)/ncol(N.matrix)
  True.predictions.orig_wp[i] <- sum(round(Predic.orig_wp[[i]],0)== N.matrix)/ncol(N.matrix)  
}
True.predictions.allGenes_nop <- sum(round(Predic.allGenes_nop,0)== N.matrix)/ncol(N.matrix)
True.predictions.allGenes_wp <- sum(round(Predic.allGenes_wp,0)== N.matrix)/ncol(N.matrix)  
True.predictions.allGenes.orig_nop <- sum(round(Predic.allGenes.orig_nop,0)== N.matrix)/ncol(N.matrix)
True.predictions.allGenes.orig_wp <- sum(round(Predic.allGenes.orig_wp,0)== N.matrix)/ncol(N.matrix)  


# Check how often the true combo gets predicted via MLE and how often it is contained in the interval
True.predictions.MLE_nop <- c()
True.predictions.MLE_wp <- c()
True.predictions.MLE_HPI_nop <- c()
True.predictions.MLE_HPI_wp <- c()
True.predictions.MLE.orig_nop <- c()
True.predictions.MLE.orig_wp <- c()
True.predictions.MLE_HPI.orig_nop <- c()
True.predictions.MLE_HPI.orig_wp <- c()
for(i in 1:gene_nr){
  True.predictions.MLE_nop[i] <- sum(MLE_nop[[i]]== N.matrix)/ncol(N.matrix)
  True.predictions.MLE_wp[i] <- sum(MLE_wp[[i]]== N.matrix)/ncol(N.matrix)
  True.predictions.MLE_HPI_nop[i] <- sum(HPI_l_nop[[i]][,1]<= N.matrix[,1] &  N.matrix[,1]<= HPI_u_nop[[i]][,1])   # only for 2Pop
  True.predictions.MLE_HPI_wp[i] <-sum(HPI_l_wp[[i]][,1]<= N.matrix[,1] &  N.matrix[,1]<= HPI_u_wp[[i]][,1])     # only for 2Pop
  True.predictions.MLE.orig_nop[i] <- sum(MLE.orig_nop[[i]]== N.matrix)/ncol(N.matrix)
  True.predictions.MLE.orig_wp[i] <- sum(MLE.orig_wp[[i]]== N.matrix)/ncol(N.matrix)
  True.predictions.MLE_HPI.orig_nop[i] <- sum(HPI.orig_l_nop[[i]][,1]<= N.matrix[,1] &  N.matrix[,1]<= HPI.orig_u_nop[[i]][,1])   # only for 2Pop
  True.predictions.MLE_HPI.orig_wp[i] <-sum(HPI.orig_l_wp[[i]][,1]<= N.matrix[,1] &  N.matrix[,1]<= HPI.orig_u_wp[[i]][,1])     # only for 2Pop
  
}

True.predictions.MLE.allGenes_nop <- sum(MLE.allGenes_nop== N.matrix)/ncol(N.matrix)
True.predictions.MLE.allGenes_wp <- sum(MLE.allGenes_wp== N.matrix)/ncol(N.matrix)
True.predictions.MLE_HPI.allGenes_nop <- sum(HPI_l.allGenes_nop[,1]<= N.matrix[,1] &  N.matrix[,1]<= HPI_u.allGenes_nop[,1])   # only for 2Pop
True.predictions.MLE_HPI.allGenes_wp <-sum(HPI_l.allGenes_wp[,1]<= N.matrix[,1] &  N.matrix[,1]<= HPI_u.allGenes_wp[,1])     # only for 2Pop

True.predictions.MLE.allGenes.orig_nop <- sum(MLE.allGenes.orig_nop== N.matrix)/ncol(N.matrix)
True.predictions.MLE.allGenes.orig_wp <- sum(MLE.allGenes.orig_wp== N.matrix)/ncol(N.matrix)
True.predictions.MLE_HPI.allGenes.orig_nop <- sum(HPI_l.allGenes.orig_nop[,1]<= N.matrix[,1] &  N.matrix[,1]<= HPI_u.allGenes.orig_nop[,1])   # only for 2Pop
True.predictions.MLE_HPI.allGenes.orig_wp <-sum(HPI_l.allGenes.orig_wp[,1]<= N.matrix[,1] &  N.matrix[,1]<= HPI_u.allGenes.orig_wp[,1])     # only for 2Pop


###############################################################################################
# Combine three genes for prediciton
###############################################################################################

Pred.all.genes_nop <- matrix(0,nrow = ncol(Dataset[1,,drop =FALSE]), ncol = ncol(j.vector.all))
Pred.all.genes_wp <- matrix(0,nrow = ncol(Dataset[1,,drop =FALSE]), ncol = ncol(j.vector.all))
Pred.all.genes.orig_nop <- matrix(0,nrow = ncol(Dataset[1,,drop =FALSE]), ncol = ncol(j.vector.all))
Pred.all.genes.orig_wp <- matrix(0,nrow = ncol(Dataset[1,,drop =FALSE]), ncol = ncol(j.vector.all))

for(i in 1:gene_nr){
  Pred.all.genes_nop <- Pred.all.genes_nop +  Predic_nop[[i]]
  Pred.all.genes_wp <- Pred.all.genes_wp +  Predic_wp[[i]]
  Pred.all.genes.orig_nop <- Pred.all.genes.orig_nop +  Predic_nop[[i]]
  Pred.all.genes.orig_wp <- Pred.all.genes.orig_wp +  Predic_wp[[i]]
}




Pred.all.genes_nop <-  Pred.all.genes_nop / gene_nr
Pred.all.genes_wp <-  Pred.all.genes_wp / gene_nr
Pred.all.genes.orig_nop <-  Pred.all.genes.orig_nop / gene_nr
Pred.all.genes.orig_wp <-  Pred.all.genes.orig_wp / gene_nr

True.predictions.all_nop <- sum(round(Pred.all.genes_nop,0)== N.matrix)/ncol(N.matrix)
True.predictions.all_wp <- sum(round(Pred.all.genes_wp,0)== N.matrix)/ncol(N.matrix)
True.predictions.all.orig_nop <- sum(round(Pred.all.genes.orig_nop,0)== N.matrix)/ncol(N.matrix)
True.predictions.all.orig_wp <- sum(round(Pred.all.genes.orig_wp,0)== N.matrix)/ncol(N.matrix)

############################# with median
Pred.median.all.genes_nop <- matrix(0,nrow = ncol(Dataset[1,,drop =FALSE]), ncol = ncol(j.vector.all))
Pred.median.all.genes_wp <- matrix(0,nrow = ncol(Dataset[1,,drop =FALSE]), ncol = ncol(j.vector.all))
Pred.median.all.genes.orig_nop <- matrix(0,nrow = ncol(Dataset[1,,drop =FALSE]), ncol = ncol(j.vector.all))
Pred.median.all.genes.orig_wp <- matrix(0,nrow = ncol(Dataset[1,,drop =FALSE]), ncol = ncol(j.vector.all))


for(i in 1:ncol(Dataset[1,,drop =FALSE])){
  for(j in 1:ncol(j.vector.all)){
    Pred.median.all.genes_nop[i,j]<- median( c(Predic_nop[[1]][i,j],Predic_nop[[2]][i,j],Predic_nop[[3]][i,j]))
    Pred.median.all.genes_wp[i,j]<- median(  c(Predic_wp[[1]][i,j],Predic_wp[[2]][i,j],Predic_wp[[3]][i,j]))
    Pred.median.all.genes.orig_nop[i,j]<- median(  c(Predic.orig_nop[[1]][i,j],Predic.orig_nop[[2]][i,j],Predic.orig_nop[[3]][i,j]))
    Pred.median.all.genes.orig_wp[i,j]<- median(  c(Predic.orig_wp[[1]][i,j],Predic.orig_wp[[2]][i,j],Predic.orig_wp[[3]][i,j]))
  }
}


True.predictions.all.median_nop <- sum(round(Pred.median.all.genes_nop,0)== N.matrix)/ncol(N.matrix)
True.predictions.all.median_wp <- sum(round(Pred.median.all.genes_wp,0)== N.matrix)/ncol(N.matrix)
True.predictions.all.median.orig_nop <- sum(round(Pred.median.all.genes.orig_nop,0)== N.matrix)/ncol(N.matrix)
True.predictions.all.median.orig_wp <- sum(round(Pred.median.all.genes.orig_wp,0)== N.matrix)/ncol(N.matrix)


######################################################################################################################################
save(p.vector, mu.vector, sigma.vector, k, n, gene_nr, set_seeds, N.matrix, Dataset, estimated.parameter, j.vector.all, dens.all, cond.dens.prob.all_nop, cond.dens.prob.all_wp, Predic_nop, Predic_wp, True.predictions_nop, True.predictions_wp, Pred.all.genes_nop, Pred.all.genes_wp, True.predictions.all_nop, True.predictions.all_wp, Pred.median.all.genes_nop,Pred.median.all.genes_wp, True.predictions.all.median_nop, True.predictions.all.median_wp, True.predictions.MLE_nop, True.predictions.MLE_wp, True.predictions.MLE_HPI_nop, True.predictions.MLE_HPI_wp, Predic_nop,Predic_wp, MLE_nop, MLE_wp, HPI_v_nop, HPI_v_wp, HPI_l_nop, HPI_l_wp, HPI_u_nop, HPI_u_wp, cond.dens.prob.all.orig_nop ,cond.dens.prob.all.orig_wp, cond.prob.allGenes_nop ,cond.prob.allGenes_wp,cond.prob.allGenes.orig_nop ,cond.prob.allGenes.orig_wp, Predic.allGenes_nop,Predic.allGenes_wp ,MLE.allGenes_nop ,MLE.allGenes_wp ,HPI_v.allGenes_nop ,HPI_v.allGenes_wp ,HPI_l.allGenes_nop ,HPI_l.allGenes_wp ,HPI_u.allGenes_nop ,HPI_u.allGenes_wp ,True.predictions.allGenes_nop ,True.predictions.allGenes_wp ,True.predictions.MLE.allGenes_nop, True.predictions.MLE.allGenes_wp, True.predictions.MLE_HPI.allGenes_nop ,True.predictions.MLE_HPI.allGenes_wp ,Predic.orig_nop ,Predic.orig_wp ,MLE.orig_nop,MLE.orig_wp,HPI.orig_v_nop ,HPI.orig_v_wp ,HPI.orig_l_nop,HPI.orig_l_wp ,HPI.orig_u_nop,HPI.orig_u_wp ,True.predictions.orig_nop ,True.predictions.orig_wp ,True.predictions.allGenes.orig_nop,True.predictions.allGenes.orig_wp ,True.predictions.MLE.orig_nop ,True.predictions.MLE.orig_wp ,True.predictions.MLE_HPI.orig_nop ,True.predictions.MLE_HPI.orig_wp ,True.predictions.MLE.allGenes.orig_nop ,True.predictions.MLE.allGenes.orig_wp,True.predictions.MLE_HPI.allGenes.orig_nop ,True.predictions.MLE_HPI.allGenes.orig_wp ,Predic.allGenes.orig_nop ,Predic.allGenes.orig_wp ,MLE.allGenes.orig_nop ,HPI_l.allGenes.orig_nop ,HPI_u.allGenes.orig_nop,HPI_v.allGenes.orig_nop ,MLE.allGenes.orig_wp,HPI_l.allGenes.orig_wp ,HPI_u.allGenes.orig_wp ,HPI_v.allGenes.orig_wp,True.predictions.all.orig_nop ,True.predictions.all.orig_wp ,Pred.all.genes.orig_nop,Pred.all.genes.orig_wp,True.predictions.all.median.orig_nop,True.predictions.all.median.orig_wp ,Pred.median.all.genes.orig_nop ,Pred.median.all.genes.orig_wp , file="Dataset1_3_all.rda")




############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

load("Dataset1_3_all.rda")



# I choose gene 2 to be displayed in the paper

i <- 2
pdf(paste0( "Hist_Dens_Wellprediction.pdf"),  width = 8, height = 3.5)
  hist(Dataset[i,],xlab="Measurements", main = paste("Histogram of simulated dataset "),freq=FALSE,breaks= 100)#seq(0,28,1))
  lines(dens.all[[i]][1,],dens.all[[i]][nrow(j.vector.all)+3,],col="#FF7F00",lwd=3)
  lines(dens.all[[i]][1,],dens.all[[i]][2,],col="#1F78B4",lwd=3)
  legend("topright", legend=c("Estimated density", "True density"), col=c("#1F78B4", "#FF7F00"),lwd=c(3,3), lty=c(2,2), cex=c(0.8,0.8))
dev.off()
  


#################################################################################

dt_new <- data.frame("Prob"=c(cond.dens.prob.all_nop[[2]]$Prob_j[1:6,],cond.dens.prob.all.orig_nop[[2]]$Prob_j[1:6,],cond.dens.prob.all_wp[[2]]$Prob_j[1:6,],cond.dens.prob.all.orig_wp[[2]]$Prob_j[1:6,]),"Obs"= rep(rep(c("Measurement 1","Measurement 2","Measurement 3","Measurement 4","Measurement 5","Measurement 6"),6),4),"Nr_Pop1"=rep(sort(rep(0:5,6)),4), "Method"=sort(rep(c(1,2,3,4),36)))
dt_new$fill <- c( "no p", "original no p","with p", "original with p")[dt_new$Method]

p1<-ggplot(dt_new[1:72,], aes(Nr_Pop1,Prob,fill=as.factor(Method)))+ 
  geom_bar(position="dodge",stat="identity")+ 
  scale_fill_brewer(palette="Paired", labels=c("estimated parameters", "true parameters", 
                                               "with p", "original with p"))+
  theme(legend.position="top")+
  ylab("Conditional probability")+
  xlab("# cells of population 1")+
  labs(fill="p excluded")+
  scale_x_continuous(breaks=0:5,labels = 0:5)+
  facet_wrap(~Obs, nrow = 6)

p2<-ggplot(dt_new[73:144,], aes(Nr_Pop1,Prob,fill=as.factor(Method)))+ 
  geom_bar(position="dodge",stat="identity")+ 
  scale_fill_manual(values= brewer.pal(n = 8, name = "Paired")[7:8], labels=c("estimated parameters", "true parameters"))+
  theme(legend.position="top")+
  ylab("Conditional probability")+
  xlab("# cells of population 1")+
  labs(fill="p included")+
  scale_x_continuous(breaks=0:5,labels = 0:5)+
  facet_wrap(~Obs, nrow = 6)

pdf(paste0( "ConDensHistos_Wellprediction.pdf"),  width = 10, height = 11)
grid.arrange( p1,p2, ncol=2)
dev.off()  
  

###############################################################################
i <- 2
dt <- data.frame(c("estimated Parameter without p","estimated Parameter with p","original Parameter without p","original Parameter with p"), c(True.predictions_nop[i],True.predictions_wp[i],True.predictions.orig_nop[i],True.predictions.orig_wp[i]), c(paste(bquote(.(True.predictions.MLE_nop[i]))," (",bquote(.(True.predictions.MLE_HPI_nop[i])) ,")", sep=""), paste(bquote(.(True.predictions.MLE_wp[i]))," (",bquote(.(True.predictions.MLE_HPI_wp[i])) ,")", sep=""),  paste(bquote(.(True.predictions.MLE.orig_nop[i]))," (",bquote(.(True.predictions.MLE_HPI.orig_nop[i])),")", sep=""),paste(bquote(.(True.predictions.MLE.orig_wp[i]))," (",bquote(.(True.predictions.MLE_HPI.orig_wp[i])) ,")", sep="") )
)
names(dt)<- c("Method","Exp Value", "MLE Genes (95% HPI)")
kable(dt)


dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"),  round(Predic_nop[[2]][1:6,1],2),  N.matrix[1:6,1])
names(dt)<- c("Observation", "estimated Parameter without p: Mean", "True Value")
kable(dt)

dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"),paste(bquote(.(MLE_nop[[2]][1:6,1],2) ) , " (", bquote(.(HPI_l_nop[[2]][1:6,1],2) ), ",", bquote(.(HPI_u_nop[[2]][1:6,1],2) ),")", sep =""), N.matrix[1:6,1])
names(dt)<- c("Observation", "estimated Parameter without p: MLE (95% HPI)", "True Value")
kable(dt)

dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"),  round(Predic_wp[[2]][1:6,1],2), N.matrix[1:6,1])
names(dt)<- c("Observation", "estimated Parameter wit p: Mean" ,"True Value")
kable(dt)

dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"), paste(bquote(.(MLE_wp[[2]][1:6,1],2) ) , " (", bquote(.(HPI_l_wp[[2]][1:6,1],2) ), ",", bquote(.(HPI_u_wp[[2]][1:6,1],2) ),")", sep =""), N.matrix[1:6,1])
names(dt)<- c("Observation", "estimated Parameter with p: MLE (95% HPI)" , "True Value")
kable(dt)

dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"), round(Predic.orig_nop[[2]][1:6,1],2), N.matrix[1:6,1])
names(dt)<- c("Observation",  "true Parameter without p: Mean", "True Value")
kable(dt)

dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"), paste(bquote(.(MLE.orig_nop[[2]][1:6,1],2) ) , " (", bquote(.(HPI.orig_l_nop[[2]][1:6,1],2) ), ",", bquote(.(HPI.orig_u_nop[[2]][1:6,1],2) ),")", sep =""), N.matrix[1:6,1])
names(dt)<- c("Observation",  "true Parameter without p: MLE (95% HPI)", "True Value")
kable(dt)

dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"), round(Predic.orig_wp[[2]][1:6,1],2), N.matrix[1:6,1])
names(dt)<- c("Observation", "true Parameter with p: Mean", "True Value")
kable(dt)

dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"),  paste(bquote(.(MLE.orig_wp[[2]][1:6,1],2) ) , " (", bquote(.(HPI.orig_l_wp[[2]][1:6,1],2) ), ",", bquote(.(HPI.orig_u_wp[[2]][1:6,1],2) ),")", sep =""), N.matrix[1:6,1])
names(dt)<- c("Observation", "true Parameter without p: MLE (95% HPI)", "True Value")
kable(dt)

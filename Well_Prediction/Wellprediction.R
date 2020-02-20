###################################################################################################
# load packages
library("stochprofML")
library("ggplot2") #für ggplot
library("dplyr")
library("knitr") 
library("RColorBrewer") # für palette
#library("gridExtra") # for grid
library("cowplot")

###################################################################################################
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
###################################################################################################



###################################################################################################
# Define Functions 
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




#Perform all Calculations needed
###############################################################################################
# Calculate all the densities for each gene fit
###############################################################################################

well.prediction <- function(dataset, k, n, p.vector, mu.vector, sigma.vector, N.matrix,set_a_seed, save_name){
    
  estimated.parameter <- list()
  TY <- length(p.vector)
  
    
  set.seed(set_a_seed)
  estimated.parameter[[1]] <- stochprof.loop(model=model,dataset=t(dataset[1,,drop =FALSE]),n=n,TY=TY,print.output=FALSE,show.plots=FALSE)
    
  j.vector.all<-stochprofML:::comb.summands(n=n, k=TY)
  #calculate all necessary densities
  dens.all <- list()
  i <- 1
  x= seq(round(min(dataset[i,])),round(max(dataset[i,])),length.out=200)
  dens.all[[i]] <- matrix(nrow=2*dim(j.vector.all)[1]+3,ncol = length(x))
  # value where densities are evaluated
  dens.all[[i]][1,] <- x
  # densities at all values for estimated parameters
  dens.all[[i]][2,] <- stochprofML:::d.sum.of.mixtures.LNLN(x,n,c(estimated.parameter[[i]]$mle[1],(1-estimated.parameter[[i]]$mle[1])),c(estimated.parameter[[i]]$mle[2],estimated.parameter[[i]]$mle[3]),rep(estimated.parameter[[i]]$mle[4],TY),logdens = F)
  # conditional densities at all values for estimated parameters conditioned on well composition
  for (j in 1:nrow(j.vector.all)){
    dens.all[[i]][j+2,]<- stochprofML:::d.sum.of.types.LNLN(x, j.vector=j.vector.all[j,], mu.vector=c(estimated.parameter[[i]]$mle[2],estimated.parameter[[i]]$mle[3]), sigma.vector=rep(estimated.parameter[[i]]$mle[4],TY), logdens = FALSE)
  }
  # densities at all values for true parameters
  dens.all[[i]][nrow(j.vector.all)+3,] <- stochprofML:::d.sum.of.mixtures.LNLN(x,n,p.vector,mu.vector,sigma.vector,logdens = F)
  # conditional densities at all values for true parameters conditioned on well composition
  for (j in 1:nrow(j.vector.all)){
    dens.all[[i]][j+nrow(j.vector.all)+3,]<- stochprofML:::d.sum.of.types.LNLN(x, j.vector=j.vector.all[j,], mu.vector=mu.vector, sigma.vector=sigma.vector, logdens = FALSE)
  }


###############################################################################################
# Calculate all conditional densities for all possible well compositions for estimated and true
# parameters for all measurements
###############################################################################################

  cond.dens.prob.all_wp <- list()
  cond.dens.prob.all.orig_wp <- list()

  for (j in 1:nrow(j.vector.all)){
    cond.dens.prob.all_wp[[i]]<-p.combo.sum.of.mix(y=t(dataset[i,,drop =FALSE]),n,p.vector.est=c(estimated.parameter[[i]]$mle[1],(1-estimated.parameter[[i]]$mle[1])),mu.vector.est=c(estimated.parameter[[i]]$mle[2],estimated.parameter[[i]]$mle[3]),sigma.vector.est=rep(estimated.parameter[[i]]$mle[4],TY))
    cond.dens.prob.all.orig_wp[[i]]<-p.combo.sum.of.mix(y=t(dataset[i,,drop =FALSE]),n,p.vector.est=p.vector,mu.vector.est=mu.vector,sigma.vector.est=sigma.vector)
  }


###############################################################################################
# Predict Cell Combos
###############################################################################################
# HPI done here only works properly for 2Pop


Predic_wp <- list()
MLE_wp <- list()
HPI_v_wp <- list()
HPI_l_wp <- list()
HPI_u_wp <- list()

Predic_wp[[i]] <- matrix(nrow = ncol(dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
MLE_wp[[i]] <- matrix(nrow = ncol(dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_l_wp[[i]] <- matrix(nrow = ncol(dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_u_wp[[i]] <- matrix(nrow = ncol(dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
HPI_v_wp[[i]] <- matrix(nrow = ncol(dataset[i,,drop =FALSE]), ncol = 1)
for (j in 1:ncol(dataset[i,,drop =FALSE])){
    Predic_wp[[i]][j,] <- colSums(j.vector.all * cond.dens.prob.all_wp[[i]]$Prob_j[j,])
    mle_hpi95_wp <- hpi(probs_distr = cond.dens.prob.all_wp[[i]]$Prob_j[j,], j.vec.all = j.vector.all, hpi_p = 0.95)
    MLE_wp[[i]][j,] <- j.vector.all[which.max(cond.dens.prob.all_wp[[i]]$Prob_j[j,]), ]
    HPI_l_wp[[i]][j,] <- mle_hpi95_wp$HPI_l
    HPI_u_wp[[i]][j,] <- mle_hpi95_wp$HPI_u
    HPI_v_wp[[i]][j,] <-  mle_hpi95_wp$HPI
    
}

# same for original parameters
Predic.orig_wp <- list()
MLE.orig_wp <- list()
HPI.orig_v_wp <- list()
HPI.orig_l_wp <- list()
HPI.orig_u_wp <- list()

  Predic.orig_wp[[i]] <- matrix(nrow = ncol(dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  MLE.orig_wp[[i]] <- matrix(nrow = ncol(dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI.orig_l_wp[[i]] <- matrix(nrow = ncol(dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI.orig_u_wp[[i]] <- matrix(nrow = ncol(dataset[i,,drop =FALSE]), ncol = ncol(j.vector.all))
  HPI.orig_v_wp[[i]] <- matrix(nrow = ncol(dataset[i,,drop =FALSE]), ncol = 1)
  for (j in 1:ncol(dataset[i,,drop =FALSE])){
    Predic.orig_wp[[i]][j,] <- colSums(j.vector.all * cond.dens.prob.all.orig_wp[[i]]$Prob_j[j,])
    mle_hpi95.orig_wp <- hpi(probs_distr = cond.dens.prob.all.orig_wp[[i]]$Prob_j[j,] , j.vec.all = j.vector.all, hpi_p = 0.95)
    MLE.orig_wp[[i]][j,] <- j.vector.all[which.max(cond.dens.prob.all.orig_wp[[i]]$Prob_j[j,]), ]
    HPI.orig_l_wp[[i]][j,] <- mle_hpi95.orig_wp$HPI_l
    HPI.orig_u_wp[[i]][j,] <- mle_hpi95.orig_wp$HPI_u
    HPI.orig_v_wp[[i]][j,] <-  mle_hpi95.orig_wp$HPI
    
  }







# Check how often the true combo gets predicted
True.predictions_wp <- c()
True.predictions.orig_wp <- c()

  True.predictions_wp[i] <- sum(round(Predic_wp[[i]],0)== N.matrix)/ncol(N.matrix)
  True.predictions.orig_wp[i] <- sum(round(Predic.orig_wp[[i]],0)== N.matrix)/ncol(N.matrix)  


# Check how often the true combo gets predicted via MLE and how often it is contained in the interval
True.predictions.MLE_wp <- c()
True.predictions.MLE_HPI_wp <- c()
True.predictions.MLE.orig_wp <- c()
True.predictions.MLE_HPI.orig_wp <- c()

  True.predictions.MLE_wp[i] <- sum(MLE_wp[[i]]== N.matrix)/ncol(N.matrix)
  True.predictions.MLE_HPI_wp[i] <-sum(HPI_l_wp[[i]][,1]<= N.matrix[,1] &  N.matrix[,1]<= HPI_u_wp[[i]][,1])     # only for 2Pop
  True.predictions.MLE.orig_wp[i] <- sum(MLE.orig_wp[[i]]== N.matrix)/ncol(N.matrix)
  True.predictions.MLE_HPI.orig_wp[i] <-sum(HPI.orig_l_wp[[i]][,1]<= N.matrix[,1] &  N.matrix[,1]<= HPI.orig_u_wp[[i]][,1])     # only for 2Pop
  





######################################################################################################################################
save(p.vector, mu.vector, sigma.vector, k, n, set_a_seed, N.matrix, dataset, estimated.parameter, j.vector.all, dens.all, cond.dens.prob.all_wp, Predic_wp,  True.predictions_wp, True.predictions.MLE_wp,  True.predictions.MLE_HPI_wp, Predic_wp,  MLE_wp, HPI_v_wp,  HPI_l_wp, HPI_u_wp, cond.dens.prob.all.orig_wp, Predic.orig_wp ,MLE.orig_wp ,HPI.orig_v_wp ,HPI.orig_l_wp ,HPI.orig_u_wp ,True.predictions.orig_wp, True.predictions.MLE.orig_wp  ,True.predictions.MLE_HPI.orig_wp , file=save_name)

}



###################################################################################################
# define Parameter values for first dataset

p.vector_1 <- c(0.2,0.8)
mu.vector_1 <- c(2,0)
sigma.vector_1 <- c(0.2,0.2)
k_1 <- 100 #number of measurement samples
n_1 <- 5 #the number of cells taken from each measurement sample
set_seeds_1 <- c(1234,345)

##############################################################################################
# Generate dataset 1
library(stochprofML)
model="LN-LN"

set_a_seed <- set_seeds_1[1]
set.seed(set_a_seed)
N.matrix_1 <- t(rmultinom(n = k_1, size = n_1, prob = p.vector_1))
Dataset_1 <- matrix(ncol=k_1, nrow = 1)

set.seed(set_a_seed+102)
Dataset_1[1,] <- r.sum.of.mixtures.LNLN(k_1, n_1, p.vector_1, mu.vector_1, sigma.vector_1, N.matrix_1)


##############################################################################################
# Calculate all predictions etc and save
well.prediction(dataset=Dataset_1, k=k_1, n=n_1, p.vector=  p.vector_1, mu.vector=mu.vector_1, sigma.vector=sigma.vector_1, N.matrix=N.matrix_1,set_a_seed=set_seeds_1[2]+102, save_name="Dataset_1_all.rda")






##############################################################################################
# Generate dataset 2
p.vector_2 <- c(0.2,0.8)
mu.vector_2 <- c(2,0)
sigma.vector_2 <- c(0.2,0.2)
k_2 <- 100 #number of measurement samples
n_2 <- 10 #the number of cells taken from each measurement sample
set_seeds_2 <- c(9876,6543)

##############################################################################################
# Generate dataset 1
set_a_seed <- set_seeds_2[1]
set.seed(set_a_seed)
N.matrix_2 <- t(rmultinom(n = k_2, size = n_2, prob = p.vector_2))
Dataset_2 <- matrix(ncol=k_2, nrow = 1)

set.seed(set_a_seed+102)
Dataset_2[1,] <- r.sum.of.mixtures.LNLN(k_2, n_2, p.vector_2, mu.vector_2, sigma.vector_2, N.matrix_2)

##############################################################################################
# Calculate all predictions etc and save
well.prediction(dataset=Dataset_2, k=k_2, n=n_2, p.vector=  p.vector_2, mu.vector=mu.vector_2, sigma.vector=sigma.vector_2, N.matrix=N.matrix_2,set_a_seed=set_seeds_2[2]+102, save_name="Dataset_2_all.rda")




# Generate the plots
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

load("Dataset_1_all.rda")



# I choose gene 2 to be displayed in the paper

i <- 1
pdf(paste0( "Hist_Dens_Wellprediction_1.pdf"),  width = 8, height = 3.5)
  hist(dataset[i,],xlab="Measurements", main = paste("Histogram of simulated dataset "),freq=FALSE,breaks= 100, cex=1.5)#seq(0,28,1))
  lines(dens.all[[i]][1,],dens.all[[i]][nrow(j.vector.all)+3,],col="#FF7F00",lwd=3)
  lines(dens.all[[i]][1,],dens.all[[i]][2,],col="#1F78B4",lwd=3)
  legend("topright", legend=c("Estimated density", "True density"), col=c("#1F78B4", "#FF7F00"),lwd=c(3,3),  bty ="n",cex=1)
dev.off()
  
############################################################################################################

load("Dataset_2_all.rda")



# I choose gene 2 to be displayed in the paper

i <- 1
pdf(paste0( "Hist_Dens_Wellprediction_2.pdf"),  width = 8, height = 3.5)
hist(dataset[i,],xlab="Measurements", main = paste("Histogram of simulated dataset "),freq=FALSE,breaks= 100, cex=1.5)#seq(0,28,1))
lines(dens.all[[i]][1,],dens.all[[i]][nrow(j.vector.all)+3,],col="#FF7F00",lwd=3)
lines(dens.all[[i]][1,],dens.all[[i]][2,],col="#1F78B4",lwd=3)
legend("topright", legend=c("Estimated density", "True density"), col=c("#1F78B4", "#FF7F00"),lwd=c(3,3),  bty ="n",cex=1)
dev.off()


#################################################################################
#Conditional Probabilities Prediciton PLot

load("Dataset_1_all.rda")
dt_new <- data.frame("Prob"=c(cond.dens.prob.all_wp[[1]]$Prob_j[1:6,],cond.dens.prob.all.orig_wp[[1]]$Prob_j[1:6,]),"Obs"= rep(rep(c("Measurement 1","Measurement 2","Measurement 3","Measurement 4","Measurement 5","Measurement 6"),6),2),"Nr_Pop1"=rep(sort(rep(0:5,6)),4), "Method"=sort(rep(c(1,2),36)),"True_CN"=rep( N.matrix[1:6,1],12),"Pred_Truth"=(rep(sort(rep(0:5,6)),4)==rep( N.matrix[1:6,1],12)), "Fill" =paste(sort(rep(c(1,2),36)),"_",(rep(sort(rep(0:5,6)),4)==rep( N.matrix[1:6,1],12))))
dt_new$fill <- c( "estimated parameters","true composition", "true parameters", "true composition")[dt_new$Fill]

p1<- ggplot(dt_new, aes(Nr_Pop1,Prob,fill=as.factor(Method),col=as.factor(fill)))+ 
  geom_bar(position="dodge",stat="identity")+ 
  scale_colour_manual(values= c(brewer.pal(n = 8, name = "Paired")[c(2)],"#000000",brewer.pal(n = 8, name = "Paired")[c(8)],"#000000"))+ 
  scale_fill_manual(values= brewer.pal(n = 8, name = "Paired")[c(2,8)], labels=c("estimated parameters", "true parameters"))+
  ylab("Conditional probability")+
  xlab("# cells of population 1")+
  labs(fill="Prediction with")+
  scale_x_continuous(breaks=0:5,labels = 0:5)+
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size =14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position="top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white")) +
  facet_wrap(~Obs, nrow = 6)
fig_label("A")



load("Dataset_2_all.rda")
dt_new <- data.frame("Prob"=c(cond.dens.prob.all_wp[[1]]$Prob_j[1:6,],cond.dens.prob.all.orig_wp[[1]]$Prob_j[1:6,]),"Obs"= rep(rep(c("Measurement 1","Measurement 2","Measurement 3","Measurement 4","Measurement 5","Measurement 6"),11),2),"Nr_Pop1"=rep(sort(rep(0:10,6)),2), "Method"=sort(rep(c(1,2),66)),"True_CN"=rep( N.matrix[1:6,1],22),"Pred_Truth"=(rep(sort(rep(0:10,6)),2)==rep( N.matrix[1:6,1],22)), "Fill" =paste(sort(rep(c(1,2),66)),"_",(rep(sort(rep(0:10,6)),2)==rep( N.matrix[1:6,1],22))))
dt_new$fill <- c( "estimated parameters","true composition", "true parameters", "true composition")[dt_new$Fill]


p2 <- ggplot(dt_new, aes(Nr_Pop1,Prob,fill=as.factor(Method),col=as.factor(fill)))+ 
  geom_bar(position="dodge",stat="identity")+ 
  scale_colour_manual(values= c(brewer.pal(n = 8, name = "Paired")[c(2)],"#000000",brewer.pal(n = 8, name = "Paired")[c(8)],"#000000"))+ 
  scale_fill_manual(values= brewer.pal(n = 8, name = "Paired")[c(2,8)], labels=c("estimated parameters", "true parameters"))+
  ylab("Conditional probability")+
  xlab("# cells of population 1")+
  labs(fill="Prediction with")+
  scale_x_continuous(breaks=0:10,labels = 0:10)+
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size =14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position="top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white")) +
  facet_wrap(~Obs, nrow = 6)








pdf(paste0( "ConDensHistos_Wellprediction_12.pdf"),  width = 14, height = 11)
plot_grid(p1, p2, align = "h", ncol = 2, rel_widths = c(1.2, 2), labels = c("AUTO"), label_size = 24)
dev.off()  



###############################################################################
# Values to Fill the table
load("Dataset_1_all.rda")

i <- 1

dt <- data.frame(c("estimated Parameter","original Parameter with p"), c(True.predictions_wp[i],True.predictions.orig_wp[i]), c(paste(bquote(.(True.predictions.MLE_wp[i]))," (",bquote(.(True.predictions.MLE_HPI_wp[i])) ,")", sep=""),paste(bquote(.(True.predictions.MLE.orig_wp[i]))," (",bquote(.(True.predictions.MLE_HPI.orig_wp[i])) ,")", sep="") )
)
names(dt)<- c("Method","Exp Value", "MLE Genes (95% HPI)")
kable(dt)




dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"),  round(Predic_wp[[1]][1:6,1],2), N.matrix[1:6,1])
names(dt)<- c("Observation", "estimated Parameter wit p: Mean" ,"True Value")
kable(dt)

dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"), paste(bquote(.(MLE_wp[[1]][1:6,1],2) ) , " (", bquote(.(HPI_l_wp[[1]][1:6,1],2) ), ",", bquote(.(HPI_u_wp[[1]][1:6,1],2) ),")", sep =""), N.matrix[1:6,1])
names(dt)<- c("Observation", "estimated Parameter with p: MLE (95% HPI)" , "True Value")
kable(dt)



dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"), round(Predic.orig_wp[[1]][1:6,1],2), N.matrix[1:6,1])
names(dt)<- c("Observation", "true Parameter with p: Mean", "True Value")
kable(dt)

dt <- data.frame(c("Observation 1","Observation 2","Observation 3","Observation 4","Observation 5","Observation 6"),  paste(bquote(.(MLE.orig_wp[[1]][1:6,1],2) ) , " (", bquote(.(HPI.orig_l_wp[[1]][1:6,1],2) ), ",", bquote(.(HPI.orig_u_wp[[1]][1:6,1],2) ),")", sep =""), N.matrix[1:6,1])
names(dt)<- c("Observation", "true Parameter without p: MLE (95% HPI)", "True Value")
kable(dt)

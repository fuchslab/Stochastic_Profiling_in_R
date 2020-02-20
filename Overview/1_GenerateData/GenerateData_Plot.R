# Generate Datasets
setwd("~/stochProfML/Paper/Pics/Overview")

#For A
library(stochprofML)
set.seed(4)
stochasticProfilingData()

# wich model Exp-LN
3
# how many populations
3
# number of samples
1000
# all same cellnumber
1
# cellnumber
10
# number of co-expressed genes
1
# probs for three populations
0.3 0.6 0.1
# log-means for the 2 lognormal pops
2 0.8
# log-standard deviation
0.4 0.4
# rate of exponential distribution
3
#output in file?
no



A <- .Last.value
########################################################

save(A, file= "A_1000.rda")

########################################################
########################################################
########################################################
# Dataset B

set.seed(10)
n_vektor <- c(5,15,8,sample(2:20,996,replace=TRUE),10)

set.seed(34)
stochasticProfilingData()

# wich model Exp-LN
3
# how many populations
3
# number of samples
100
# all same cellnumber
2
# enter name of variable
3
#variable name
n_vektor
# vektor is converted to matrix with 1 column and 100 rows, what are the observations
1
# number of co-expressed genes
1
# probs for three populations
0.3 0.6 0.1
# log-means for the 2 lognormal pops
2 0.8
# log-standard deviation
0.4 0.4
# rate of exponential distribution
3
#output in file?
no



B <- .Last.value

########################################################

save(B, file= "B_1000.rda")

########################################################
########################################################
########################################################

# save and plot the data in Histograms

# Plot each in a histogram

load(file="A_1000.rda")
hist(A, breaks=seq(-0.5, 100, by = 1), freq = FALSE, col = "black", main = "Sample A")

load(file="B_1000.rda")
hist(B, breaks=seq(-0.5, 150, by = 1), freq = FALSE, col = "black", main = "Sample B")


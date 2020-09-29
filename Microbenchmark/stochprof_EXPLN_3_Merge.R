library(stochprofML)

res_EXPLN_3_all <- c()

for(i in 1:5){
    filename = paste0("EXPLN_3_1times_",i,".rda")
    try(load(filename), silent = TRUE)
    try(res_EXPLN_3_all[i] <- res_EXPLN_3$time , silent = TRUE)
    try(rm(res_EXPLN_3), silent = TRUE)
}

#Results in seconds: /1000000000
# calculate median (and mean)
# select min
# select max


median(res_EXPLN_3_all)/1000000000
# Median is run number 1

mean(res_EXPLN_3_all)/1000000000
# Mean is close to median

min(res_EXPLN_3_all)/1000000000
# Minimum is Run number 3

max(res_EXPLN_3_all)/1000000000
# Maximum is run number 5

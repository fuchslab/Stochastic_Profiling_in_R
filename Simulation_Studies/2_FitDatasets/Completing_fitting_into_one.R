############## Merge all fits into one file (per Set)

#Set1
path <- "SimulationStudies/2_FitDatasets/Set1"

complete_fitting <- list()

load(file = paste(path, "Set1_1_fitting.rda", sep = "/"))
complete_fitting[[1]] <- result_fit_list

load(file = paste(path, "Set1_2_fitting.rda", sep = "/"))
complete_fitting[[2]] <- result_fit_list

load(file = paste(path, "Set1_5_fitting.rda", sep = "/"))
complete_fitting[[3]] <- result_fit_list

load(file = paste(path, "Set1_10_fitting.rda", sep = "/"))
complete_fitting[[4]] <- result_fit_list

load(file = paste(path, "Set1_15_fitting.rda", sep = "/"))
complete_fitting[[5]] <- result_fit_list

load(file = paste(path, "Set1_20_fitting.rda", sep = "/"))
complete_fitting[[6]] <- result_fit_list

load(file = paste(path, "Set1_50_fitting.rda", sep = "/"))
complete_fitting[[7]] <- result_fit_list

load(file = paste(path, "Set1_mix_fitting.rda", sep = "/"))
complete_fitting[[8]] <- result_fit_list

load(file = paste(path, "Set1_mix2_fitting.rda", sep = "/"))
complete_fitting[[9]] <- result_fit_list

save(complete_fitting, file = "Complete_fitting.rda")
#########################################################
#Set2
path <- "SimulationStudies/2_FitDatasets/Set2"

complete_fitting <- list()

load(file = paste(path, "Set2_1_fitting.rda", sep = "/"))
complete_fitting[[1]] <- result_fit_list

load(file = paste(path, "Set2_2_fitting.rda", sep = "/"))
complete_fitting[[2]] <- result_fit_list

load(file = paste(path, "Set2_5_fitting.rda", sep = "/"))
complete_fitting[[3]] <- result_fit_list

load(file = paste(path, "Set2_10_fitting.rda", sep = "/"))
complete_fitting[[4]] <- result_fit_list

load(file = paste(path, "Set2_15_fitting.rda", sep = "/"))
complete_fitting[[5]] <- result_fit_list

load(file = paste(path, "Set2_20_fitting.rda", sep = "/"))
complete_fitting[[6]] <- result_fit_list

load(file = paste(path,"Set2_50_fitting.rda", sep = "/"))
complete_fitting[[7]]<- result_fit_list

load(file = paste(path, "Set2_mix_fitting.rda", sep = "/"))
complete_fitting[[8]] <- result_fit_list

load(file = paste(path, "Set2_mix2_fitting.rda", sep = "/"))
complete_fitting[[9]] <- result_fit_list

save(complete_fitting, file = "Complete_fitting.rda")
#########################################################
#Set3
path <- "SimulationStudies/2_FitDatasets/Set3"

complete_fitting <- list()

load(file = paste(path, "Set3_1_fitting.rda", sep = "/"))
complete_fitting[[1]] <- result_fit_list

load(file = paste(path, "Set3_2_fitting.rda", sep = "/"))
complete_fitting[[2]] <- result_fit_list

load(file = paste(path, "Set3_5_fitting.rda", sep = "/"))
complete_fitting[[3]] <- result_fit_list

load(file = paste(path, "Set3_10_fitting.rda", sep = "/"))
complete_fitting[[4]] <- result_fit_list

load(file = paste(path, "Set3_15_fitting.rda", sep = "/"))
complete_fitting[[5]] <- result_fit_list

load(file = paste(path, "Set3_20_fitting.rda", sep = "/"))
complete_fitting[[6]] <- result_fit_list

load(file = paste(path, "Set3_50_fitting.rda", sep = "/"))
complete_fitting[[7]] <- result_fit_list

load(file = paste(path,"Set3_mix_fitting.rda", sep = "/"))
complete_fitting[[8]]<- result_fit_list

load(file = paste(path,"Set3_mix2_fitting.rda", sep = "/"))
complete_fitting[[9]]<- result_fit_list

save(complete_fitting, file = "Complete_fitting.rda")
#########################################################
#Set4
path <- "SimulationStudies/2_FitDatasets/Set4"

complete_fitting <- list()

load(file = paste(path, "Set4_1_fitting.rda", sep = "/"))
complete_fitting[[1]] <- result_fit_list

load(file = paste(path, "Set4_2_fitting.rda", sep = "/"))
complete_fitting[[2]] <- result_fit_list

load(file = paste(path, "Set4_5_fitting.rda", sep = "/"))
complete_fitting[[3]] <- result_fit_list

load(file = paste(path, "Set4_10_fitting.rda", sep = "/"))
complete_fitting[[4]] <- result_fit_list

load(file = paste(path, "Set4_15_fitting.rda", sep = "/"))
complete_fitting[[5]] <- result_fit_list

load(file = paste(path, "Set4_20_fitting.rda", sep = "/"))
complete_fitting[[6]] <- result_fit_list

load(file = paste(path, "Set4_50_fitting.rda", sep = "/"))
complete_fitting[[7]] <- result_fit_list

load(file = paste(path, "Set4_mix_fitting.rda", sep = "/"))
complete_fitting[[8]] <- result_fit_list

load(file = paste(path, "Set4_mix2_fitting.rda", sep = "/"))
complete_fitting[[9]] <- result_fit_list

save(complete_fitting, file = "Complete_fitting.rda")

#########################################################
#Set5
path <- "SimulationStudies/2_FitDatasets/Set5"

complete_fitting <- list()

load(file = paste(path," Set5_1_fitting.rda", sep = "/"))
complete_fitting[[1]] <- result_fit_list

load(file = paste(path, "Set5_2_fitting.rda", sep = "/"))
complete_fitting[[2]] <- result_fit_list

load(file = paste(path, "Set5_5_fitting.rda", sep = "/"))
complete_fitting[[3]] <- result_fit_list

load(file = paste(path, "Set5_10_fitting.rda", sep = "/"))
complete_fitting[[4]] <- result_fit_list

load(file = paste(path, "Set5_15_fitting.rda", sep = "/"))
complete_fitting[[5]] <- result_fit_list

load(file = paste(path, "Set5_20_fitting.rda", sep = "/"))
complete_fitting[[6]] <- result_fit_list

load(file = paste(path, "Set5_50_fitting.rda", sep = "/"))
complete_fitting[[7]] <- result_fit_list

load(file = paste(path, "Set5_mix_fitting.rda", sep = "/"))
complete_fitting[[8]] <- result_fit_list

load(file = paste(path, "Set5_mix2_fitting.rda", sep = "/"))
complete_fitting[[9]] <- result_fit_list

save(complete_fitting, file = "Complete_fitting.rda")

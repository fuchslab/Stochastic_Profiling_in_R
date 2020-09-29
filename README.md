# Stochastic Profiling in R
This repository accompanies the article: Lisa Amrhein and Christiane Fuchs (2019) "stochprofML: Stochastic Profiling using Maximum Likelihood Estimation in R".

All code used in this article can be found in this repository. 

# Overview Figure 2
All code can be found in the folder **Overview**.
- *1_generateData/GenerateData_Plot.R* generates the two datasets and the histogram in the upper part of the Figure.
- *2_OriginalDistributions/OrigPop_Plot.R* calculates the density and generates the plots shown in the Stochastic Profiling model.
- *3_ModelFitting/FitBothDatasets.R* applies to both datasets the stochprofML algorithm and saves the results.
- *4_PlotFittedModels/PlotA_results.R* generates all plots in the lower part of the Figure. 

# Microbenchmark Table 1
All code can be found in the folder **Microbenchmark**.
- *RunMicrobenchmark* contains all scripts needed to run the microbenchmarking with 5 repetitions.
- *stochprof_EXPLN_3_Merge.R* merges the separate runs of the 3 Population EXPLN model and calculates the table entries.

# Example Figure 3
All code can be found in the folder **Example_3Pop**.
- *3Pop_LNLN.R* generates the 1000 datasets and performes parameter estimation.
- *Plot_results.R* generates Figure 3.

# Usage of stochprofML: Figures 4 + 5 and Additional File 5
All code to generate these figures is directly given in the paper. Same holds for the code given in Appendix E. Nevertheless, we included everything in this repository in the folder **stochprofML_Usage**.
- *Script/Script.R* contains all functions to generate Figures 4 and 5.
- *Interactive/Interactive.R* shows all code to use the interactive functions presented in Additional File 5.


# Simulation studies
All scripts needed for the simulation studies can be found in folder **Simulation_Studies**. All simulation studies are based on the same datasets.
- *1_GenerateDatasets/GenerateDatasets.R* generates for all five parameter settings datasets for all different pool size settings. These are saved in the respective folders.
- *2_FitDatasets/FitDatasets.R* applies the stochprofML algorithm to all generated datasets. 
- *2_FitDatasets/Completing_fitting_into_one.R* merges all resulting fits in one file per parameter setting.

## Simulation study on optimal pool size: Figures 6 and S1-S4
- *3_Plots/Simstudy_PoolSizes_plots.R* generates the five Figures wehere for each parameter setting the parameter estimates are compared dependend on the pool sizes.

## Simulation study on impact of parameter values: Figures 7 and S5-S12
- *3_Plots/Simstudy_ParameterChanges_plots.R* generates the nine figures where for each pool size settings the parameter estimates are compared for small parameter changes.

## Simulation study on uncertainty in pool sizes: Figures 8 and 9
The files necessary for this simulation study can be found in the folder **Uncertainty_Cellnumber**.

### Figure 8
- *Uncertainty_10/FitDatasets_Uncertainty_10.R* is based on one dataset generated in the pervious simulation studies and thus loads *1_GenerateDatasets/Set1/Set1_10.rda*. 
  The fit with the true n, wich was already performed before and given in *2_FitDatasets/Set1/Set1_10_fitting.rda* is also used, this file is also loaded. The script performs the analogous fits with the other wrong cellnumbers.
- *Uncertainty_10/Uncertainty_plot_10.R* creates the Figure.

### Figure 9
- *Uncertainty_mix/Fitting_Set1_mix1_n1000.R* is based on one dataset generated in the pervious simulation studies and thus loads *1_GenerateDatasets/Set1/Set1_mix.rda*. 
  This script generates new n-vectors and uses them in new stochprofML runs.
- *Uncertainty_mix/Fitting_Set1_mix2_n1000.R* is based on one dataset generated in the pervious simulation studies and thus loads *1_GenerateDatasets/Set1/Set1_mix2.rda*. 
  This script generates new n-vectors and uses them in new stochprofML runs.
- *Uncertainty_mix/Uncertainty_mix_plots.R* generates the figure. As the fit with the true n-vectors are also needed (and where already performed before) the file with those results, namely 
  *2_FitDatasets/Set1/Complete_fitting.rda* is needed.

# Comparison of populations: Figures 10 and 11
All scripts needed for these figures are given in **Overlap**.
- *Overlap_LN_LN.R* contains the overlap function which is also displayed in the article.
- *OverlapTwoLN_Example.R* generates Figure 9.
- *Overlap_D.R* generates Figure 10.

# Cell pool composition: Figures 12 and 13 and Table 2
All scripts needed for these figures are given in **Well_Prediction**.
- *Wellprediction.R* contains all code needed for the whole application and generates the Figures and the Table.



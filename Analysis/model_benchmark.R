rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(dampack)  # for CEA and calculate ICERs
library(tidyverse)
library(rbenchmark)
library(microbenchmark)
library(tictoc)
library(rBeta2009)
library(parallel)
library(foreach)
library(doParallel)
library(tidyr)

# Call model setup functions
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/calibration_functions.R")
source("R/ICER_functions.R")
source("R/generate_psa_parameters.R")

# Load parameters
# Load parameters
source("Analysis/00_load_parameters.R")

# Set number of cores
n_cores <- detectCores()
#n_cores <- 4
registerDoParallel(n_cores)

# Benchmark baseline deterministic model
# df_model_benchmark <- microbenchmark("Markov Model (Base)"    = markov_model(l_params_all = l_params_all),
#                                      "Markov Model (Outputs)" = outcomes(l_params_all = l_params_all, v_params_calib = v_calib_post_map),
#                                      times = 100)
# plot <- autoplot(df_model_benchmark)
# 
# ggsave(plot, 
#        filename = "Plots/Benchmark/model_benchmark.png", 
#        width = 10, height = 7)
# 
# # Benchmark model calibration
# l_cali_targets <- list(ODF = read.csv(file = "data/cali_target_odf.csv", header = TRUE),
#                        ODN = read.csv(file = "data/cali_target_odn.csv", header = TRUE))
# 
# n_cali_max_per <- max(c(l_cali_targets$ODF$Time, l_cali_targets$ODN$Time))
# 
# df_model_benchmark_cali <- microbenchmark("Markov Model (Cali)" = markov_model(l_params_all = l_params_all, cali = TRUE), times = 100)
# plot_cali <- autoplot(df_model_benchmark_cali)
# 
# ggsave(plot_cali, 
#        filename = "Plots/Benchmark/model_benchmark_cali.png", 
#        width = 10, height = 7)

# Benchmark PSA

# Set population size for dirichlet draws
n_pop_cohort <- 29000
n_pop_trial  <- 272
n_sim <- 10 # just to test function (will be set as n_sim)

df_psa_params_MMS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_trial, scenario = "MMS",
                                         file.death_hr = "data/death_hr.csv",
                                         file.frailty = "data/frailty.csv",
                                         file.weibull = "data/Modified Model Specification/weibull.csv",
                                         file.unconditional = "data/Modified Model Specification/unconditional.csv",
                                         file.overdose = "data/overdose.csv",
                                         file.fentanyl = "data/fentanyl.csv",
                                         file.hiv = "data/hiv_sero.csv",
                                         file.hcv = "data/hcv_sero.csv",
                                         file.costs = "data/Modified Model Specification/costs.csv",
                                         file.crime_costs = "data/Modified Model Specification/crime_costs.csv",
                                         file.qalys = "data/Modified Model Specification/qalys.csv",
                                         file.imis_output = "outputs/Calibration/imis_output.RData")
df_outcomes_MET_PSA_MMS <- data.frame()
df_outcomes_BUP_PSA_MMS <- data.frame()
df_incremental_PSA_MMS <- data.frame()
df_ICER_PSA_MMS <- data.frame()

combine_custom_i <- function(LL1, LL2) {
  df_outcomes_MET_PSA_MMS <- rbind(LL1$df_outcomes_MET_PSA_MMS, LL2$df_outcomes_MET_PSA_MMS)
  df_outcomes_BUP_PSA_MMS <- rbind(LL1$df_outcomes_BUP_PSA_MMS, LL2$df_outcomes_BUP_PSA_MMS)
  df_incremental_PSA_MMS  <- rbind(LL1$df_incremental_PSA_MMS, LL2$df_incremental_PSA_MMS)
  df_ICER_PSA_MMS <- rbind(LL1$df_ICER_PSA_MMS, LL2$df_ICER_PSA_MMS)
  
  return(list(df_outcomes_MET_PSA_MMS = df_outcomes_MET_PSA_MMS, 
              df_outcomes_BUP_PSA_MMS = df_outcomes_BUP_PSA_MMS, 
              df_incremental_PSA_MMS = df_incremental_PSA_MMS,
              df_ICER_PSA_MMS = df_ICER_PSA_MMS))
}

df_model_benchmark_PSA <- microbenchmark("PSA - Single Core" = {for (i in 1:n_sim){
  # Update parameter set for each scenario with next set of PSA drawn parameters
  l_psa_input_MET_MMS <- update_param_list(l_params_all = l_params_MET_MMS, params_updated = df_psa_params_MMS[i, ])
  l_psa_input_BUP_MMS <- update_param_list(l_params_all = l_params_BUP_MMS, params_updated = df_psa_params_MMS[i, ])
  
  # Run model and generate outputs
  l_outcomes_MET_MMS <- outcomes(l_params_all = l_psa_input_MET_MMS, v_params_calib = v_calib_post_map, PSA = TRUE)
  l_outcomes_BUP_MMS <- outcomes(l_params_all = l_psa_input_BUP_MMS, v_params_calib = v_calib_post_map, PSA = TRUE)
  
  # Extract cost and QALY outputs
  df_outcomes_MET_PSA_MMS <- rbind(df_outcomes_MET_PSA_MMS, l_outcomes_MET_MMS$df_outcomes)
  df_outcomes_BUP_PSA_MMS <- rbind(df_outcomes_BUP_PSA_MMS, l_outcomes_BUP_MMS$df_outcomes)
  
  # Calculate ICER (societal and health sector perspective)
  l_ICER_MMS <- ICER(outcomes_comp = l_outcomes_MET_MMS, outcomes_int = l_outcomes_BUP_MMS)
  
  df_incremental_PSA_MMS <- rbind(df_incremental_PSA_MMS, l_ICER_MMS$df_incremental)
  
  df_ICER_PSA_MMS <- rbind(df_ICER_PSA_MMS, l_ICER_MMS$df_icer)
}}, 
"PSA - Multi Core" = {foreach(i = 1:n_sim, .combine = combine_custom_i, .packages = 'tidyr') %dopar% {
  # Update parameter set for each scenario with next set of PSA drawn parameters
  l_psa_input_MET_MMS <- update_param_list(l_params_all = l_params_MET_MMS, params_updated = df_psa_params_MMS[i, ])
  l_psa_input_BUP_MMS <- update_param_list(l_params_all = l_params_BUP_MMS, params_updated = df_psa_params_MMS[i, ])
  
  # Run model and generate outputs
  l_outcomes_MET_MMS <- outcomes(l_params_all = l_psa_input_MET_MMS, v_params_calib = v_calib_post_map, PSA = TRUE)
  l_outcomes_BUP_MMS <- outcomes(l_params_all = l_psa_input_BUP_MMS, v_params_calib = v_calib_post_map, PSA = TRUE)
  
  # Extract cost and QALY outputs
  df_outcomes_MET_PSA_MMS <- l_outcomes_MET_MMS$df_outcomes
  df_outcomes_BUP_PSA_MMS <- l_outcomes_BUP_MMS$df_outcomes
  
  # Calculate ICER (societal and health sector perspective)
  l_ICER_MMS <- ICER(outcomes_comp = l_outcomes_MET_MMS, outcomes_int = l_outcomes_BUP_MMS)
  
  df_incremental_PSA_MMS <- l_ICER_MMS$df_incremental
  
  df_ICER_PSA_MMS <- l_ICER_MMS$df_icer
  
  return(list(df_outcomes_MET_PSA_MMS = df_outcomes_MET_PSA_MMS, 
              df_outcomes_BUP_PSA_MMS = df_outcomes_BUP_PSA_MMS, 
              df_incremental_PSA_MMS = df_incremental_PSA_MMS,
              df_ICER_PSA_MMS = df_ICER_PSA_MMS))}},
times = 2)

plot_PSA <- autoplot(df_model_benchmark_PSA)

ggsave(plot_PSA, 
       filename = "Plots/Benchmark/model_benchmark_PSA.png", 
       width = 10, height = 7)
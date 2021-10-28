rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(dampack)  # for CEA and calculate ICERs
library(tidyverse)
library(rbenchmark)
library(microbenchmark)

# Call model setup functions
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/calibration_functions.R")
source("R/ICER_functions.R")

# Load parameters
l_params_all <- load_all_params(file.init = "data/init_params.csv",
                                          file.init_dist = "data/init_dist.csv",
                                          file.mort = "data/all_cause_mortality.csv",
                                          file.death_hr = "data/death_hr.csv",
                                          file.frailty = "data/frailty.csv",
                                          file.weibull_scale = "data/weibull_scale.csv",
                                          file.weibull_shape = "data/weibull_shape.csv",
                                          file.unconditional = "data/unconditional.csv",
                                          file.overdose = "data/overdose.csv",
                                          file.hiv = "data/hiv_sero.csv",
                                          file.hcv = "data/hcv_sero.csv",
                                          file.costs = "data/costs.csv",
                                          file.crime_costs = "data/crime_costs.csv",
                                          file.qalys = "data/qalys.csv")

# Calibrated parameter values
load(file = "outputs/imis_output.RData")

#### Load calibration targets ####
l_cali_targets <- list(ODF = read.csv(file = "data/cali_target_odf.csv", header = TRUE),
                       ODN = read.csv(file = "data/cali_target_odn.csv", header = TRUE))

# Benchmark baseline deterministic model
df_model_benchmark <- microbenchmark("Markov Model (Base)"    = markov_model(l_params_all = l_params_all),
                                     "Markov Model (Outputs)" = outcomes(l_params_all = l_params_all, v_params_calib = v_calib_post_map),
                                     times = 100)
plot <- autoplot(df_model_benchmark)

ggsave(plot, 
       filename = "Plots/Benchmark/model_benchmark.png", 
       width = 10, height = 7)

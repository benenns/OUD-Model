library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(dampack)  # for CEA and calculate ICERs
library(tidyverse)
library(rbenchmark)

# Call model setup functions
source("R/OPTIMA_00_input_parameter_functions.R")
source("R/OPTIMA_01_model_setup_functions.R")

# Load parameters
benchmark(l_params_all <- load_all_params(file.init = "data/init_params.csv",
                                          file.init_dist = "data/init_dist_bup.csv", # Change initial distributions (100% in BUP)
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
                                          file.qalys = "data/qalys.csv"), replications = 1)

# Run model
benchmark(l_out_markov <- markov_model(l_params_all = l_params_all, err_stop = FALSE, verbose = TRUE), replications = 1)

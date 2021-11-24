library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(rBeta2009)

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")
source("R/generate_psa_parameters.R")
source("R/model_setup_functions.R")
source("R/ICER_functions.R")

# Load parameters
# Calibrated parameter values
load(file = "outputs/imis_output.RData")

# All parameters
# BNX scenario (primary model definition)
l_params_BUP <- load_all_params(file.init = "data/init_params.csv",
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
                                file.qalys = "data/qalys.csv")

# Methadone scenario (primary model definition)
l_params_MET <- load_all_params(file.init = "data/init_params.csv",
                                file.init_dist = "data/init_dist_met.csv", # Change initial distributions (100% in MET)
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

### PSA model outputs
### Run Markov model for PSA draws and return outputs ###
# Generate PSA parameter draws
df_psa_params <- generate_psa_params(n_sim = 1000, seed = 3730687, n_samp = 250,
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
                                     file.qalys = "data/qalys.csv",
                                     file.imis_output = "outputs/imis_output.RData")


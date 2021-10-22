library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")
source("R/generate_psa_parameter_functions.R")
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

### Main deterministic model outputs ###
# Run Markov model and return outputs (using MAP point estimates from posterior distribution for calibrated params)
l_outcomes_MET <- outcomes(l_params_all = l_params_MET, v_params_calib = v_calib_post_map)
l_outcomes_BUP <- outcomes(l_params_all = l_params_BUP, v_params_calib = v_calib_post_map)

# Calculate ICERs
ICER <- ICER(outcomes_comp = l_outcomes_MET, outcomes_int = l_outcomes_BUP)

# Output to csv files
# Full model trace
write.csv(outcomes_MET$m_M_trace,"outputs/trace/trace_MET.csv", row.names = TRUE)
write.csv(outcomes_BUP$m_M_trace,"outputs/trace/trace_BUP.csv", row.names = TRUE)

# Full model costs
write.csv(outcomes_MET$m_TOTAL_costs_states,"outputs/trace/full_trace_costs_MET.csv", row.names = TRUE)
write.csv(outcomes_BUP$m_TOTAL_costs_states,"outputs/trace/full_trace_costs_BUP.csv", row.names = TRUE)

# Costs
# Total
write.csv(outcomes_MET$v_costs,"outputs/costs/costs_MET.csv", row.names = TRUE)
write.csv(outcomes_BUP$v_costs,"outputs/costs/costs_BUP.csv", row.names = TRUE)

# Treatment
write.csv(outcomes_MET$m_TX_costs,"outputs/costs/tx_costs_MET.csv", row.names = TRUE)
write.csv(outcomes_BUP$m_TX_costs,"outputs/costs/tx_costs_BUP.csv", row.names = TRUE)

# Health sector
write.csv(outcomes_MET$m_HRU_costs,"outputs/costs/hru_costs_MET.csv", row.names = TRUE)
write.csv(outcomes_BUP$m_HRU_costs,"outputs/costs/hru_costs_BUP.csv", row.names = TRUE)

# Crime
write.csv(outcomes_MET$m_crime_costs,"outputs/costs/crime_costs_MET.csv", row.names = TRUE)
write.csv(outcomes_BUP$m_crime_costs,"outputs/costs/crime_costs_BUP.csv", row.names = TRUE)

# QALYs
write.csv(outcomes_MET$v_qalys,"outputs/qalys/qalys_MET.csv", row.names = TRUE)
write.csv(outcomes_BUP$v_qalys,"outputs/qalys/qalys_BUP.csv", row.names = TRUE)

# ICER
write.csv(ICER$v_icer,"outputs/ICER/ICER.csv", row.names = TRUE)

### PSA model outputs
### Run Markov model for PSA draws and return outputs ###
# Generate PSA parameter draws
df_psa_params <- generate_psa_params(n_sim = 2000, seed = 3730687, n_samp = 250,
                                     file.death_hr = NULL,
                                     file.frailty = NULL,
                                     file.weibull_scale = NULL,
                                     file.weibull_shape = NULL,
                                     file.unconditional = NULL,
                                     file.overdose = NULL,
                                     file.hiv = NULL,
                                     file.hcv = NULL,
                                     file.costs = NULL,
                                     file.crime_costs = NULL,
                                     file.qalys = NULL,
                                     file.imis_output = NULL)


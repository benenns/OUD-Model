rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(xlsx)

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")
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

#### OUTPUT RESULTS ####
# Full model trace
write.csv(l_outcomes_MET$m_M_trace,"outputs/trace/trace_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP$m_M_trace,"outputs/trace/trace_BUP.csv", row.names = TRUE)

# Full model costs
write.csv(l_outcomes_MET$m_TOTAL_costs_states,"outputs/trace/full_trace_costs_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP$m_TOTAL_costs_states,"outputs/trace/full_trace_costs_BUP.csv", row.names = TRUE)

# Outcomes
# Disaggregated
df_outcomes <- data.frame(l_outcomes_BUP$v_outcomes, l_outcomes_MET$v_outcomes)
df_outcomes <- df_costs %>% rename("Early take-home BNX" = l_outcomes_BUP.v_outcomes, "Methadone" = l_outcomes_MET.v_outcomes)

# ICER
df_icer <- data.frame(ICER$v_icer)
df_icer <- df_icer %>% rename("Early take-home BNX vs. Methadone" = ICER.v_icer)

# Output
write.csv(df_outcomes,"outputs/main_output_det.csv", row.names = TRUE)
write.csv(df_icer,"outputs/icer_det.csv", row.names = TRUE)

#write.xlsx2(df_outcomes, file = "outputs/main_output_det.xlsx", sheetName = "Outcomes",
#            col.names = TRUE, row.names = TRUE, append = FALSE)
#write.xlsx2(df_icer, file = "outputs/main_output_det.xlsx", sheetName = "ICER",
#            col.names = TRUE, row.names = TRUE, append = TRUE)

# Raw outputs
# Costs
write.csv(l_outcomes_MET$v_costs,"outputs/costs/costs_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP$v_costs,"outputs/costs/costs_BUP.csv", row.names = TRUE)

# Treatment
write.csv(l_outcomes_MET$m_TX_costs,"outputs/costs/tx_costs_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP$m_TX_costs,"outputs/costs/tx_costs_BUP.csv", row.names = TRUE)

# Health sector
write.csv(l_outcomes_MET$m_HRU_costs,"outputs/costs/hru_costs_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP$m_HRU_costs,"outputs/costs/hru_costs_BUP.csv", row.names = TRUE)

# Crime
write.csv(l_outcomes_MET$m_crime_costs,"outputs/costs/crime_costs_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP$m_crime_costs,"outputs/costs/crime_costs_BUP.csv", row.names = TRUE)

# QALYs
write.csv(l_outcomes_MET$v_qalys,"outputs/qalys/qalys_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP$v_qalys,"outputs/qalys/qalys_BUP.csv", row.names = TRUE)

# ICER
write.csv(ICER$v_icer,"outputs/ICER/ICER.csv", row.names = TRUE)


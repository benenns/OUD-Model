rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(formattable)
#library(xlsx)

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/ICER_functions.R")

# Load parameters
source("Analysis/00_load_parameters.R") # load all model parameters for each scenario + calibrated parameters

### Main deterministic model outputs ###
# Run Markov model and return outputs (using MAP point estimates from posterior distribution for calibrated params)
#l_outcomes_MET <- outcomes(l_params_all = l_params_MET, v_params_calib = v_calib_post_map)
#l_outcomes_BUP <- outcomes(l_params_all = l_params_BUP, v_params_calib = v_calib_post_map)

# Calculate ICERs
#ICER <- ICER(outcomes_comp = l_outcomes_MET, outcomes_int = l_outcomes_BUP)

#### Produce model outputs ####
#### Modified Model Specification ####
l_outcomes_BUP_MMS  <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map)
l_outcomes_MET_MMS  <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map)
l_outcomes_validation_MMS  <- outcomes(l_params_all = l_params_all_validation_MMS, v_params_calib = v_calib_post_map)

df_outcomes_MMS <- rbind(l_outcomes_BUP_MMS$df_outcomes, l_outcomes_MET_MMS$df_outcomes)
rownames(df_outcomes_MMS) <- c("Early take-home BNX", "Methadone")

# Save output
saveRDS(l_outcomes_BUP_MMS, file = "outputs/outcomes/outcomes_BUP_MMS.RData")
saveRDS(l_outcomes_MET_MMS, file = "outputs/outcomes/outcomes_MET_MMS.RData")
saveRDS(df_outcomes_MMS, file = "outputs/outcomes/outcomes_MMS.RData")

#### Trial Specification ####
# l_outcomes_BUP_TS  <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map)
# l_outcomes_MET_TS  <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map)
# 
# df_outcomes_TS <- rbind(l_outcomes_BUP_TS$df_outcomes, l_outcomes_MET_TS$df_outcomes)
# rownames(df_outcomes_TS) <- c("Early take-home BNX", "Methadone")
# 
# # Save output
# saveRDS(l_outcomes_BUP_TS, file = "outputs/outcomes/outcomes_BUP_TS.RData")
# saveRDS(l_outcomes_MET_TS, file = "outputs/outcomes/outcomes_MET_TS.RData")
# saveRDS(df_outcomes_TS, file = "outputs/outcomes/df_outcomes_TS.RData")

# Generate ICERs
l_ICER_MMS <- ICER(outcomes_comp = l_outcomes_MET_MMS, outcomes_int = l_outcomes_BUP_MMS)
#l_ICER_TS <- ICER(outcomes_comp = l_outcomes_MET_TS, outcomes_int = l_outcomes_BUP_TS)

# Save output
saveRDS(l_ICER_MMS, file = "outputs/outcomes/ICER_MMS.RData")
#saveRDS(l_ICER_TS, file = "outputs/outcomes/ICER_TS.RData")

######################################
#### Modified Model Specification ####
######################################
# Full model trace
write.csv(l_outcomes_MET_MMS$m_M_trace,"outputs/trace/Modified Model Specification/trace_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP_MMS$m_M_trace,"outputs/trace/Modified Model Specification/trace_BUP.csv", row.names = TRUE)
write.csv(l_outcomes_validation_MMS$m_M_trace,"outputs/trace/Modified Model Specification/trace_validation.csv", row.names = TRUE)

# Aggregate trace
write.csv(l_outcomes_MET_MMS$m_M_agg_trace,"outputs/trace/Modified Model Specification/trace_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP_MMS$m_M_agg_trace,"outputs/trace/Modified Model Specification/trace_BUP.csv", row.names = TRUE)
write.csv(l_outcomes_validation_MMS$m_M_agg_trace,"outputs/trace/Modified Model Specification/trace_validation.csv", row.names = TRUE)

# Full model costs
write.csv(l_outcomes_MET_MMS$m_TOTAL_costs_states,"outputs/trace/Modified Model Specification/full_trace_costs_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP_MMS$m_TOTAL_costs_states,"outputs/trace/Modified Model Specification/full_trace_costs_BUP.csv", row.names = TRUE)
#write.csv(l_outcomes_validation_MMS$m_TOTAL_costs_states,"outputs/trace/Modified Model Specification/full_trace_costs_BUP.csv", row.names = TRUE)

# Outcomes
# Disaggregated
df_outcomes_MMS <- rbind(l_outcomes_BUP_MMS$df_outcomes, l_outcomes_MET_MMS$df_outcomes)
rownames(df_outcomes_MMS) <- c("Early take-home BNX", "Methadone")

# ICER
df_icer_MMS <- l_ICER_MMS$df_icer
rownames(df_icer_MMS) <- c("Early take-home BNX vs. Methadone")

# Incremental costs & QALYs
df_incremental_MMS <- l_ICER_MMS$df_incremental
rownames(df_incremental_MMS) <- c("Early take-home BNX vs. Methadone")

# Output
save(df_incremental_MMS, 
     file = "outputs/ICER/incremental_det_MMS.RData")

write.csv(df_outcomes_MMS,"outputs/main_output_det_MMS.csv", row.names = TRUE)
write.csv(df_icer_MMS,"outputs/ICER/icer_det_MMS.csv", row.names = TRUE)
write.csv(df_incremental_MMS,"outputs/ICER/incremental_det_MMS.csv", row.names = TRUE)

# Raw outputs
# Costs
write.csv(l_outcomes_MET_MMS$v_costs,"outputs/costs/Modified Model Specification/costs_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP_MMS$v_costs,"outputs/costs/Modified Model Specification/costs_BUP.csv", row.names = TRUE)

# Treatment
write.csv(l_outcomes_MET_MMS$m_TX_costs,"outputs/costs/Modified Model Specification/tx_costs_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP_MMS$m_TX_costs,"outputs/costs/Modified Model Specification/tx_costs_BUP.csv", row.names = TRUE)

# Health sector
write.csv(l_outcomes_MET_MMS$m_HRU_costs,"outputs/costs/Modified Model Specification/hru_costs_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP_MMS$m_HRU_costs,"outputs/costs/Modified Model Specification/hru_costs_BUP.csv", row.names = TRUE)

# Crime
write.csv(l_outcomes_MET_MMS$m_crime_costs,"outputs/costs/Modified Model Specification/crime_costs_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP_MMS$m_crime_costs,"outputs/costs/Modified Model Specification/crime_costs_BUP.csv", row.names = TRUE)

# QALYs
write.csv(l_outcomes_MET_MMS$v_qalys,"outputs/qalys/Modified Model Specification/qalys_MET.csv", row.names = TRUE)
write.csv(l_outcomes_BUP_MMS$v_qalys,"outputs/qalys/Modified Model Specification/qalys_BUP.csv", row.names = TRUE)

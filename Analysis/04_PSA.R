rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(rBeta2009)
library(parallel)
library(foreach)
library(doParallel)
library(tidyr)
#library(mail)
library(RPushbullet)

# Set number of cores
#n_cores <- detectCores()
n_cores <- 4
registerDoParallel(n_cores)

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")
source("R/generate_psa_parameters.R")
source("R/model_setup_functions.R")
source("R/ICER_functions.R")

# Load parameters
source("Analysis/00_load_parameters.R")

# Set population size for dirichlet draws
n_pop_cohort <- 29000
n_pop_trial  <- 272
n_sim <- 5000 # just to test function (will be set as n_sim)

### PSA model outputs
### Run Markov model for PSA draws and return outputs ###
# Generate PSA parameter draws

######################################
#### Modified Model Specification ####
######################################
## Base case (EQ-5D-5L)
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

write.csv(df_psa_params_MMS,"outputs/PSA/Modified Model Specification/input_PSA_MMS.csv", row.names = TRUE)

#############################
#### Trial Specification ####
#############################
# Need to draw parameters for both scenarios due to differences in allowed transitions (Dirichlet)
## Base case (EQ-5D-5L)
# BNX scenario
df_psa_params_BUP_TS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_trial, scenario = "TS_BUP",
                                            file.death_hr = "data/death_hr.csv",
                                            file.frailty = "data/frailty.csv",
                                            file.weibull = "data/Trial Specification/weibull.csv",
                                            file.unconditional = "data/Trial Specification/unconditional_bup.csv",
                                            file.overdose = "data/overdose.csv",
                                            file.fentanyl = "data/fentanyl.csv",
                                            file.hiv = "data/hiv_sero.csv",
                                            file.hcv = "data/hcv_sero.csv",
                                            file.costs = "data/Trial Specification/costs_bup.csv",
                                            file.crime_costs = "data/Trial Specification/crime_costs_bup.csv",
                                            file.qalys = "data/Trial Specification/qalys_bup.csv",
                                            file.imis_output = "outputs/Calibration/imis_output.RData")

write.csv(df_psa_params_BUP_TS,"outputs/PSA/Trial Specification/input_PSA_BUP_TS.csv", row.names = TRUE)

# Extract non-state-exit parameters from BUP
#df_psa_params_TS <- df_psa_params_BUP_TS %>% select(-c(p_BUP_BUPC_NI:p_ODN_REL_INJ))

# MET scenario
df_psa_params_MET_TS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_trial, scenario = "TS_MET",
                                            file.death_hr = "data/death_hr.csv",
                                            file.frailty = "data/frailty.csv",
                                            file.weibull = "data/Trial Specification/weibull.csv",
                                            file.unconditional = "data/Trial Specification/unconditional_met.csv",
                                            file.overdose = "data/overdose.csv",
                                            file.fentanyl = "data/fentanyl.csv",
                                            file.hiv = "data/hiv_sero.csv",
                                            file.hcv = "data/hcv_sero.csv",
                                            file.costs = "data/Trial Specification/costs_met.csv",
                                            file.crime_costs = "data/Trial Specification/crime_costs_met.csv",
                                            file.qalys = "data/Trial Specification/qalys_met.csv",
                                            file.imis_output = "outputs/Calibration/imis_output.RData")

write.csv(df_psa_params_MET_TS,"outputs/PSA/Trial Specification/input_PSA_MET_TS.csv", row.names = TRUE)

# Extract MET state-exit
#df_psa_params_MET_TS_UP <- df_psa_params_MET_TS %>% select(c(p_BUP_BUPC_NI:p_ODN_REL_INJ))

# Add MET state-exit to overall
#df_psa_params_MET_TS <- bind_cols(df_psa_params_TS, df_psa_params_MET_TS_UP)

# Initialize data frames
# # Modified Model Specification
# df_outcomes_MET_PSA_MMS <- data.frame()
# df_outcomes_BUP_PSA_MMS <- data.frame()
# df_incremental_PSA_MMS <- data.frame()
# df_ICER_PSA_MMS <- data.frame()
# 
# # Trial Specification
# df_outcomes_MET_PSA_TS <- data.frame()
# df_outcomes_BUP_PSA_TS <- data.frame()
# df_incremental_PSA_TS <- data.frame()
# df_ICER_PSA_TS <- data.frame()

######################################
#### Modified Model Specification ####
######################################
combine_custom_MMS <- function(LL1, LL2) {
  df_outcomes_MET_PSA_MMS <- rbind(LL1$df_outcomes_MET_PSA_MMS, LL2$df_outcomes_MET_PSA_MMS)
  df_outcomes_BUP_PSA_MMS <- rbind(LL1$df_outcomes_BUP_PSA_MMS, LL2$df_outcomes_BUP_PSA_MMS)
  df_incremental_PSA_MMS  <- rbind(LL1$df_incremental_PSA_MMS, LL2$df_incremental_PSA_MMS)
  df_ICER_PSA_MMS <- rbind(LL1$df_ICER_PSA_MMS, LL2$df_ICER_PSA_MMS)
  
  return(list(df_outcomes_MET_PSA_MMS = df_outcomes_MET_PSA_MMS, 
              df_outcomes_BUP_PSA_MMS = df_outcomes_BUP_PSA_MMS, 
              df_incremental_PSA_MMS = df_incremental_PSA_MMS,
              df_ICER_PSA_MMS = df_ICER_PSA_MMS))
}

# Run PSA for each block to help memory issues
n_sim <- 1500
n_block_size <- 500 # size of block for each loop
n_blocks <- n_sim/n_block_size
n_start <- 500 # set to 0 if running full PSA

# Initialize lists
l_outcomes_MET_PSA_MMS <- list()
l_outcomes_BUP_PSA_MMS <- list()
l_ICER_PSA_MMS         <- list()        
l_incremental_PSA_MMS  <- list()

#k<-0
#s1 <- n_start + 1 + k*n_block_size
#e1 <- n_start + (k + 1)*n_block_size

for (j in (0:(n_blocks - 1))){
  l_PSA_MMS <- foreach(i = (n_start + 1 + j*n_block_size):(n_start + (j + 1)*n_block_size), .combine = combine_custom_MMS, .packages = 'tidyr') %dopar% {
    # Update parameter set for each scenario with next set of PSA drawn parameters
    l_psa_input_MET_MMS <- update_param_list(l_params_all = l_params_MET_MMS, params_updated = df_psa_params_MMS[i, ])
    l_psa_input_BUP_MMS <- update_param_list(l_params_all = l_params_BUP_MMS, params_updated = df_psa_params_MMS[i, ])
  
    # Run model and generate outputs
    l_outcomes_MET_MMS <- outcomes(l_params_all = l_psa_input_MET_MMS, v_params_calib = v_calib_post_map, PSA = TRUE)
    l_outcomes_BUP_MMS <- outcomes(l_params_all = l_psa_input_BUP_MMS, v_params_calib = v_calib_post_map, PSA = TRUE)
  
    df_outcomes_MET_PSA_MMS <- l_outcomes_MET_MMS$df_outcomes
    df_outcomes_BUP_PSA_MMS <- l_outcomes_BUP_MMS$df_outcomes
  
    # Calculate ICER (societal and health sector perspective)
    l_ICER_MMS <- ICER(outcomes_comp = l_outcomes_MET_MMS, outcomes_int = l_outcomes_BUP_MMS)
  
    #df_incremental_PSA_MMS <- rbind(df_incremental_PSA_MMS, l_ICER_MMS$df_incremental)
    df_incremental_PSA_MMS <- l_ICER_MMS$df_incremental
  
    #df_ICER_PSA_MMS <- rbind(df_ICER_PSA_MMS, l_ICER_MMS$df_icer)
    df_ICER_PSA_MMS <- l_ICER_MMS$df_icer
    
    return(list(df_outcomes_MET_PSA_MMS = df_outcomes_MET_PSA_MMS, 
                df_outcomes_BUP_PSA_MMS = df_outcomes_BUP_PSA_MMS, 
                df_incremental_PSA_MMS = df_incremental_PSA_MMS,
                df_ICER_PSA_MMS = df_ICER_PSA_MMS))
  }
  
  df_outcomes_MET_PSA_MMS <- l_PSA_MMS$df_outcomes_MET_PSA_MMS
  df_outcomes_BUP_PSA_MMS <- l_PSA_MMS$df_outcomes_BUP_PSA_MMS
  df_ICER_PSA_MMS         <- l_PSA_MMS$df_ICER_PSA_MMS
  df_incremental_PSA_MMS  <- l_PSA_MMS$df_incremental_PSA_MMS
  
  l_outcomes_MET_PSA_MMS[[j + 1]] <- df_outcomes_MET_PSA_MMS
  l_outcomes_BUP_PSA_MMS[[j + 1]] <- df_outcomes_BUP_PSA_MMS
  l_ICER_PSA_MMS[[j + 1]]         <- df_ICER_PSA_MMS
  l_incremental_PSA_MMS[[j + 1]]  <- df_incremental_PSA_MMS
}

#sendmail("ben.enns@gmail.com", subject = "Notification from R", message = "PSA finished running!", password="rmail")
pbPost("note", "Notification from R", "PSA runs complete (office desktop)")

# combine output data sets
# initialize empty data frames
df_outcomes_MET_PSA_MMS_comb <- df_outcomes_BUP_PSA_MMS_comb <- df_ICER_PSA_MMS_comb <- df_incremental_PSA_MMS_comb <- data.frame()

for (i in 1:n_blocks){
  df_temp <- l_outcomes_MET_PSA_MMS[[i]]
  df_outcomes_MET_PSA_MMS_comb <- rbind(df_outcomes_MET_PSA_MMS_comb, df_temp)
}
for (i in 1:n_blocks){
  df_temp <- l_outcomes_BUP_PSA_MMS[[i]]
  df_outcomes_BUP_PSA_MMS_comb <- rbind(df_outcomes_BUP_PSA_MMS_comb, df_temp)
}
for (i in 1:n_blocks){
  df_temp <- l_ICER_PSA_MMS[[i]]
  df_ICER_PSA_MMS_comb <- rbind(df_ICER_PSA_MMS_comb, df_temp)
}
for (i in 1:n_blocks){
  df_temp <- l_incremental_PSA_MMS[[i]]
  df_incremental_PSA_MMS_comb <- rbind(df_incremental_PSA_MMS_comb, df_temp)
}

#stopImplicitCluster()

# df_outcomes_MET_PSA_MMS <- l_incremental_PSA_MMS$df_outcomes_MET_PSA_MMS
# df_outcomes_BUP_PSA_MMS <- l_incremental_PSA_MMS$df_outcomes_BUP_PSA_MMS
# df_ICER_PSA_MMS         <- l_incremental_PSA_MMS$df_ICER_PSA_MMS
# df_incremental_PSA_MMS  <- l_incremental_PSA_MMS$df_incremental_PSA_MMS

stopImplicitCluster()

### Output results
## As .RData
save(df_outcomes_MET_PSA_MMS_comb, 
     file = "outputs/PSA/Modified Model Specification/outcomes_MET_PSA_MMS.RData")
save(df_outcomes_BUP_PSA_MMS_comb, 
     file = "outputs/PSA/Modified Model Specification/outcomes_BUP_PSA_MMS.RData")
save(df_ICER_PSA_MMS_comb, 
     file = "outputs/PSA/Modified Model Specification/ICER_PSA_MMS.RData")
save(df_incremental_PSA_MMS_comb,
     file = "outputs/PSA/Modified Model Specification/incremental_PSA_MMS.RData")
## As .csv
write.csv(df_outcomes_MET_PSA_MMS_comb, 
          file = "outputs/PSA/Modified Model Specification/outcomes_MET_PSA_MMS.csv",
          row.names = FALSE)
write.csv(df_outcomes_BUP_PSA_MMS_comb, 
          file = "outputs/PSA/Modified Model Specification/outcomes_BUP_PSA_MMS.csv",
          row.names = FALSE)
write.csv(df_ICER_PSA_MMS_comb, 
          file = "outputs/PSA/Modified Model Specification/ICER_PSA_MMS.csv",
          row.names = FALSE)
write.csv(df_incremental_PSA_MMS_comb, 
          file = "outputs/PSA/Modified Model Specification/incremental_PSA_MMS.csv",
          row.names = FALSE)

#############################
#### Trial Specification ####
#############################
# combine_custom_TS <- function(LL1, LL2) {
#   df_outcomes_MET_PSA_TS <- rbind(LL1$df_outcomes_MET_PSA_TS, LL2$df_outcomes_MET_PSA_TS)
#   df_outcomes_BUP_PSA_TS <- rbind(LL1$df_outcomes_BUP_PSA_TS, LL2$df_outcomes_BUP_PSA_TS)
#   df_incremental_PSA_TS  <- rbind(LL1$df_incremental_PSA_TS, LL2$df_incremental_PSA_TS)
#   df_ICER_PSA_TS <- rbind(LL1$df_ICER_PSA_TS, LL2$df_ICER_PSA_TS)
#   
#   return(list(df_outcomes_MET_PSA_TS = df_outcomes_MET_PSA_TS, 
#               df_outcomes_BUP_PSA_TS = df_outcomes_BUP_PSA_TS, 
#               df_incremental_PSA_TS = df_incremental_PSA_TS,
#               df_ICER_PSA_TS = df_ICER_PSA_TS))
# }
# 
# l_incremental_PSA_TS <- foreach(i = 1:n_sim, .combine = combine_custom_TS, .packages = 'tidyr') %dopar% { # i <- 1
#   # Update parameter set for each scenario with next set of PSA drawn parameters
#   l_psa_input_MET_TS <- update_param_list(l_params_all = l_params_MET_TS, params_updated = df_psa_params_MET_TS[i, ])
#   l_psa_input_BUP_TS <- update_param_list(l_params_all = l_params_BUP_TS, params_updated = df_psa_params_BUP_TS[i, ])
#   
#   # Run model and generate outputs
#   l_outcomes_MET_TS <- outcomes(l_params_all = l_psa_input_MET_TS, v_params_calib = v_calib_post_map, PSA = TRUE)
#   l_outcomes_BUP_TS <- outcomes(l_params_all = l_psa_input_BUP_TS, v_params_calib = v_calib_post_map, PSA = TRUE)
#   
#   df_outcomes_MET_PSA_TS <- l_outcomes_MET_TS$df_outcomes
#   df_outcomes_BUP_PSA_TS <- l_outcomes_BUP_TS$df_outcomes
#   
#   # Calculate ICER (societal and health sector perspective)
#   l_ICER_TS <- ICER(outcomes_comp = l_outcomes_MET_TS, outcomes_int = l_outcomes_BUP_TS)
#   
#   df_incremental_PSA_TS <- l_ICER_TS$df_incremental
#   
#   df_ICER_PSA_TS <- l_ICER_TS$df_icer
#   
#   return(list(df_outcomes_MET_PSA_TS = df_outcomes_MET_PSA_TS, 
#               df_outcomes_BUP_PSA_TS = df_outcomes_BUP_PSA_TS, 
#               df_incremental_PSA_TS = df_incremental_PSA_TS,
#               df_ICER_PSA_TS = df_ICER_PSA_TS))
# }
# stopImplicitCluster()
# 
# df_outcomes_MET_PSA_TS <- l_incremental_PSA_TS$df_outcomes_MET_PSA_TS
# df_outcomes_BUP_PSA_TS <- l_incremental_PSA_TS$df_outcomes_BUP_PSA_TS
# df_ICER_PSA_TS         <- l_incremental_PSA_TS$df_ICER_PSA_TS
# df_incremental_PSA_TS  <- l_incremental_PSA_TS$df_incremental_PSA_TS
# 
# ### Output results
# ## As .RData
# save(df_outcomes_MET_PSA_TS,
#      file = "outputs/PSA/Trial Specification/outcomes_MET_PSA_TS.RData")
# save(df_outcomes_BUP_PSA_TS,
#      file = "outputs/PSA/Trial Specification/outcomes_BUP_PSA_TS.RData")
# save(df_ICER_PSA_TS,
#      file = "outputs/PSA/Trial Specification/ICER_PSA_TS.RData")
# save(df_incremental_PSA_TS,
#      file = "outputs/PSA/Trial Specification/incremental_PSA_TS.RData")
# ## As .csv
# write.csv(df_outcomes_MET_PSA_TS,
#           file = "outputs/PSA/Trial Specification/outcomes_MET_PSA_TS.csv",
#           row.names = FALSE)
# write.csv(df_outcomes_BUP_PSA_TS,
#           file = "outputs/PSA/Trial Specification/outcomes_BUP_PSA_TS.csv",
#           row.names = FALSE)
# write.csv(df_ICER_PSA_TS,
#           file = "outputs/PSA/Trial Specification/ICER_PSA_TS.csv",
#           row.names = FALSE)
# write.csv(df_incremental_PSA_TS,
#           file = "outputs/PSA/Trial Specification/incremental_PSA_TS.csv",
#           row.names = FALSE)


### Process PSA results
## Read-in saved results
## Modified Model Specification
load(file = "outputs/PSA/Modified Model Specification/outcomes_MET_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/outcomes_BUP_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/incremental_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/ICER_PSA_MMS.RData")

## Trial Specification
# load(file = "outputs/PSA/Trial Specification/outcomes_MET_PSA_TS.RData")
# load(file = "outputs/PSA/Trial Specification/outcomes_BUP_PSA_TS.RData")
# load(file = "outputs/PSA/Trial Specification/incremental_PSA_TS.RData")
# load(file = "outputs/PSA/Trial Specification/ICER_PSA_TS.RData")

# ## TEMP CODE ##
# m_outcomes_BUP_PSA_MMS <- as.matrix(df_outcomes_BUP_PSA_MMS)
# m_outcomes_MET_PSA_MMS <- as.matrix(df_outcomes_MET_PSA_MMS)
# 
# m_incremental_PSA_MMS <- m_outcomes_BUP_PSA_MMS - m_outcomes_MET_PSA_MMS
# df_incremental_PSA_MMS <- as.data.frame(m_incremental_PSA_MMS)
# names(df_incremental_PSA_MMS) <- c("n_inc_costs_TOTAL_6mo", "n_inc_costs_TOTAL_1yr", "n_inc_costs_TOTAL_5yr", "n_inc_costs_TOTAL_10yr", "n_inc_costs_TOTAL_life", 
#                                    "n_inc_costs_HEALTH_SECTOR_6mo", "n_inc_costs_HEALTH_SECTOR_1yr", "n_inc_costs_HEALTH_SECTOR_5yr", "n_inc_costs_HEALTH_SECTOR_10yr", "n_inc_costs_HEALTH_SECTOR_life",
#                                    "n_inc_costs_CRIMINAL_6mo", "n_inc_costs_CRIMINAL_1yr", "n_inc_costs_CRIMINAL_5yr", "n_inc_costs_CRIMINAL_10yr", "n_inc_costs_CRIMINAL_life",  
#                                    "n_inc_costs_TX_6mo", "n_inc_costs_TX_1yr", "n_inc_costs_TX_5yr", "n_inc_costs_TX_10yr", "n_inc_costs_TX_life", 
#                                    "n_inc_costs_HRU_6mo", "n_inc_costs_HRU_1yr", "n_inc_costs_HRU_5yr", "n_inc_costs_HRU_10yr", "n_inc_costs_HRU_life",
#                                    "n_inc_qalys_TOTAL_6mo", "n_inc_qalys_TOTAL_1yr", "n_inc_qalys_TOTAL_5yr", "n_inc_qalys_TOTAL_10yr", "n_inc_qalys_TOTAL_life")

#df_ICER_PSA_MMS <- df_ICER_PSA_MMS
#df_outcomes_BUP_PSA_MMS <- df_outcomes_BUP_PSA_MMS_comb
#df_outcomes_MET_PSA_MMS <- df_outcomes_MET_PSA_MMS_comb

### Summary stats ###
## Modified Model Specification ##
# Methadone
tbl_df_summary_MET_MMS <- df_outcomes_MET_PSA_MMS %>% as.tibble() %>% select(n_TOTAL_costs_6mo, 
                                                                             n_HEALTH_SECTOR_costs_6mo, 
                                                                             n_CRIMINAL_costs_6mo,
                                                                             n_TX_costs_6mo,
                                                                             n_HRU_costs_6mo,
                                                                             n_TOTAL_qalys_6mo,
                                                                             n_TOTAL_costs_10yr, 
                                                                             n_HEALTH_SECTOR_costs_10yr, 
                                                                             n_CRIMINAL_costs_10yr,
                                                                             n_TX_costs_10yr,
                                                                             n_HRU_costs_10yr,
                                                                             n_TOTAL_qalys_10yr,
                                                                             n_TOTAL_costs_life, 
                                                                             n_HEALTH_SECTOR_costs_life, 
                                                                             n_CRIMINAL_costs_life,
                                                                             n_TX_costs_life,
                                                                             n_HRU_costs_life,
                                                                             n_TOTAL_qalys_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))
# BNX
tbl_df_summary_BUP_MMS <- df_outcomes_BUP_PSA_MMS %>% as.tibble() %>% select(n_TOTAL_costs_6mo, 
                                                                             n_HEALTH_SECTOR_costs_6mo, 
                                                                             n_CRIMINAL_costs_6mo,
                                                                             n_TX_costs_6mo,
                                                                             n_HRU_costs_6mo,
                                                                             n_TOTAL_qalys_6mo,
                                                                             n_TOTAL_costs_10yr, 
                                                                             n_HEALTH_SECTOR_costs_10yr, 
                                                                             n_CRIMINAL_costs_10yr,
                                                                             n_TX_costs_10yr,
                                                                             n_HRU_costs_10yr,
                                                                             n_TOTAL_qalys_10yr,
                                                                             n_TOTAL_costs_life, 
                                                                             n_HEALTH_SECTOR_costs_life, 
                                                                             n_CRIMINAL_costs_life,
                                                                             n_TX_costs_life,
                                                                             n_HRU_costs_life,
                                                                             n_TOTAL_qalys_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

# Incremental
tbl_df_summary_incremental_MMS <- df_incremental_PSA_MMS %>% as.tibble() %>% select(n_inc_costs_TOTAL_6mo,
                                                                                    n_inc_costs_TOTAL_10yr,
                                                                                    n_inc_costs_TOTAL_life,
                                                                                    n_inc_costs_HEALTH_SECTOR_6mo,
                                                                                    n_inc_costs_HEALTH_SECTOR_10yr,
                                                                                    n_inc_costs_HEALTH_SECTOR_life,
                                                                                    n_inc_qalys_TOTAL_6mo,
                                                                                    n_inc_qalys_TOTAL_10yr,
                                                                                    n_inc_qalys_TOTAL_life,
                                                                                    n_inc_costs_TX_6mo,
                                                                                    n_inc_costs_TX_10yr,
                                                                                    n_inc_costs_TX_life,
                                                                                    n_inc_costs_HRU_6mo,
                                                                                    n_inc_costs_HRU_10yr,
                                                                                    n_inc_costs_HRU_life,
                                                                                    n_inc_costs_CRIMINAL_6mo,
                                                                                    n_inc_costs_CRIMINAL_10yr,
                                                                                    n_inc_costs_CRIMINAL_life) %>%
  
  # tbl_df_summary_incremental_MMS <- df_incremental_PSA_MMS %>% as.tibble() %>% select(n_TOTAL_costs_6mo,
  #                                                                                     n_TOTAL_costs_10yr,
  #                                                                                     n_TOTAL_costs_life,
  #                                                                                     n_HEALTH_SECTOR_costs_6mo,
  #                                                                                     n_HEALTH_SECTOR_costs_10yr,
  #                                                                                     n_HEALTH_SECTOR_costs_life,
  #                                                                                     n_TOTAL_qalys_6mo,
  #                                                                                     n_TOTAL_qalys_10yr,
  #                                                                                     n_TOTAL_qalys_life,
  #                                                                                     n_TX_costs_6mo,
  #                                                                                     n_TX_costs_10yr,
  #                                                                                     n_TX_costs_life,
  #                                                                                     n_HRU_costs_6mo,
  #                                                                                     n_HRU_costs_10yr,
  #                                                                                     n_HRU_costs_life,
  #                                                                                     n_CRIMINAL_costs_6mo,
  #                                                                                     n_CRIMINAL_costs_10yr,
  #                                                                                     n_CRIMINAL_costs_life) %>%
  
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

# ICER
tbl_df_summary_ICER_MMS <- df_ICER_PSA_MMS %>% as.tibble() %>% select(n_icer_TOTAL_6mo,
                                                                      n_icer_HEALTH_SECTOR_6mo,
                                                                      n_icer_TOTAL_10yr,
                                                                      n_icer_HEALTH_SECTOR_10yr,
                                                                      n_icer_TOTAL_life,
                                                                      n_icer_HEALTH_SECTOR_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

# ICER CI STABILITY
tbl_df_summary_ICER_MMS_2000 <- df_ICER_PSA_MMS[1:2000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_2500 <- df_ICER_PSA_MMS[1:2500, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_3000 <- df_ICER_PSA_MMS[1:3000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_3500 <- df_ICER_PSA_MMS[1:3500, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_4000 <- df_ICER_PSA_MMS[1:4000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_4500 <- df_ICER_PSA_MMS[1:4500, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_5000 <- df_ICER_PSA_MMS[1:5000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_CI_stability <- rbind(tbl_df_summary_ICER_MMS_2000, tbl_df_summary_ICER_MMS_2500, tbl_df_summary_ICER_MMS_3000, tbl_df_summary_ICER_MMS_3500, tbl_df_summary_ICER_MMS_4000, tbl_df_summary_ICER_MMS_4500, tbl_df_summary_ICER_MMS_5000)
#add_rownames(tbl_df_CI_stability) <- c("2000", "2500", "3000", "3500", "4000")


## As .csv
write.csv(tbl_df_summary_MET_MMS,
          file = "outputs/PSA/Modified Model Specification/summary_outcomes_MET_PSA_MMS.csv",
          row.names = FALSE)
write.csv(tbl_df_summary_BUP_MMS,
          file = "outputs/PSA/Modified Model Specification/summary_outcomes_BUP_PSA_MMS.csv",
          row.names = FALSE)
write.csv(tbl_df_summary_incremental_MMS,
          file = "outputs/PSA/Modified Model Specification/summary_incremental_PSA_MMS.csv",
          row.names = FALSE)
write.csv(tbl_df_summary_ICER_MMS,
          file = "outputs/PSA/Modified Model Specification/summary_ICER_PSA_MMS.csv",
          row.names = FALSE)

## As .RData ##
save(tbl_df_summary_incremental_MMS, 
     file = "outputs/PSA/Modified Model Specification/summary_incremental_PSA_MMS.RData")

## Trial Specification ##
# Methadone
# tbl_df_summary_MET_TS <- df_outcomes_MET_PSA_TS %>% as.tibble() %>% select(n_TOTAL_costs_6mo, 
#                                                                            n_HEALTH_SECTOR_costs_6mo, 
#                                                                            n_CRIMINAL_costs_6mo,
#                                                                            n_TX_costs_6mo,
#                                                                            n_HRU_costs_6mo,
#                                                                            n_TOTAL_qalys_6mo) %>% gather("variable", "value") %>% 
#   group_by(variable) %>% summarize(mean = mean(value),
#                                    sd = sd(value),
#                                    q50 = quantile(value, probs = .5),
#                                    q025 = quantile(value, probs = .025),
#                                    q975 = quantile(value, probs = .975),
#                                    min = min(value),
#                                    max = max(value))
# BNX
# tbl_df_summary_BUP_TS <- df_outcomes_BUP_PSA_TS %>% as.tibble() %>% select(n_TOTAL_costs_6mo, 
#                                                                            n_HEALTH_SECTOR_costs_6mo, 
#                                                                            n_CRIMINAL_costs_6mo,
#                                                                            n_TX_costs_6mo,
#                                                                            n_HRU_costs_6mo,
#                                                                            n_TOTAL_qalys_6mo) %>% gather("variable", "value") %>% group_by(variable) %>% summarize(mean = mean(value),
#                                                                                                                                        sd = sd(value),
#                                                                                                                                        q50 = quantile(value, probs = .5),
#                                                                                                                                        q025 = quantile(value, probs = .025),
#                                                                                                                                        q975 = quantile(value, probs = .975),
#                                                                                                                                        min = min(value),
#                                                                                                                                        max = max(value))
# ICER
# tbl_df_summary_ICER_TS <- df_ICER_PSA_TS %>% as.tibble() %>% select(n_icer_TOTAL_6mo,
#                                                                     n_icer_HEALTH_SECTOR_6mo) %>% gather("variable", "value") %>% group_by(variable) %>% summarize(mean = mean(value),
#                                                                                                                                 sd = sd(value),
#                                                                                                                                 q50 = quantile(value, probs = .5),
#                                                                                                                                 q025 = quantile(value, probs = .025),
#                                                                                                                                 q975 = quantile(value, probs = .975),
#                                                                                                                                 min = min(value),
#                                                                                                                                 max = max(value))

df_PSA_summary <- as.data.frame()

## Incremental ##
# MMS
tbl_df_summary_inc_MMS <- df_incremental_PSA_MMS %>% as.tibble() %>% mutate(BNX_dominate_TOTAL_6mo = ifelse(n_inc_costs_TOTAL_6mo < 0 & n_inc_qalys_TOTAL_6mo > 0, 1, 0),
                                                                            BNX_dominate_TOTAL_10yr = ifelse(n_inc_costs_TOTAL_10yr < 0 & n_inc_qalys_TOTAL_10yr > 0, 1, 0),
                                                                            BNX_dominate_TOTAL_life = ifelse(n_inc_costs_TOTAL_life < 0 & n_inc_qalys_TOTAL_life > 0, 1, 0),
                                                                            BNX_dominate_HEALTH_SECTOR_6mo = ifelse(n_inc_costs_HEALTH_SECTOR_6mo < 0 & n_inc_qalys_TOTAL_6mo > 0, 1, 0),
                                                                            BNX_dominate_HEALTH_SECTOR_10yr = ifelse(n_inc_costs_HEALTH_SECTOR_10yr < 0 & n_inc_qalys_TOTAL_10yr > 0, 1, 0),
                                                                            BNX_dominate_HEALTH_SECTOR_life = ifelse(n_inc_costs_HEALTH_SECTOR_life < 0 & n_inc_qalys_TOTAL_life > 0, 1, 0),
                                                                            MET_dominate_TOTAL_6mo = ifelse(n_inc_costs_TOTAL_6mo > 0 & n_inc_qalys_TOTAL_6mo < 0, 1, 0),
                                                                            MET_dominate_TOTAL_10yr = ifelse(n_inc_costs_TOTAL_10yr > 0 & n_inc_qalys_TOTAL_10yr < 0, 1, 0),
                                                                            MET_dominate_TOTAL_life = ifelse(n_inc_costs_TOTAL_life > 0 & n_inc_qalys_TOTAL_life < 0, 1, 0),
                                                                            MET_dominate_HEALTH_SECTOR_6mo = ifelse(n_inc_costs_HEALTH_SECTOR_6mo > 0 & n_inc_qalys_TOTAL_6mo < 0, 1, 0),
                                                                            MET_dominate_HEALTH_SECTOR_10yr = ifelse(n_inc_costs_HEALTH_SECTOR_10yr > 0 & n_inc_qalys_TOTAL_10yr < 0, 1, 0),
                                                                            MET_dominate_HEALTH_SECTOR_life = ifelse(n_inc_costs_HEALTH_SECTOR_life > 0 & n_inc_qalys_TOTAL_life < 0, 1, 0))

# TS
# tbl_df_summary_inc_TS <- df_incremental_PSA_TS %>% as.tibble() %>% mutate(BNX_dominate_TOTAL_6mo = ifelse(n_inc_costs_TOTAL_6mo < 0 & n_inc_qalys_TOTAL_6mo > 0, 1, 0),
#                                                                           BNX_dominate_HEALTH_SECTOR_6mo = ifelse(n_inc_costs_HEALTH_SECTOR_6mo < 0 & n_inc_qalys_TOTAL_6mo > 0, 1, 0),
#                                                                           MET_dominate_TOTAL_6mo = ifelse(n_inc_costs_TOTAL_6mo > 0 & n_inc_qalys_TOTAL_6mo < 0, 1, 0),
#                                                                           MET_dominate_HEALTH_SECTOR_6mo = ifelse(n_inc_costs_HEALTH_SECTOR_6mo > 0 & n_inc_qalys_TOTAL_6mo < 0, 1, 0))

# MMS
tbl_df_dominant_sim_MMS <- tbl_df_summary_inc_MMS %>% select(BNX_dominate_TOTAL_6mo, BNX_dominate_TOTAL_10yr, BNX_dominate_TOTAL_life, BNX_dominate_HEALTH_SECTOR_6mo, BNX_dominate_HEALTH_SECTOR_10yr, BNX_dominate_HEALTH_SECTOR_life,
                                                             MET_dominate_TOTAL_6mo, MET_dominate_TOTAL_10yr, MET_dominate_TOTAL_life, MET_dominate_HEALTH_SECTOR_6mo, MET_dominate_HEALTH_SECTOR_10yr, MET_dominate_HEALTH_SECTOR_life)
# TS
#tbl_df_dominant_sim_TS <- tbl_df_summary_inc_TS %>% select(BNX_dominate_TOTAL_6mo, BNX_dominate_HEALTH_SECTOR_6mo, MET_dominate_TOTAL_6mo, MET_dominate_HEALTH_SECTOR_6mo)

# MMS
tbl_df_dom_summary_MMS <- summarise_all(tbl_df_dominant_sim_MMS, mean)
## As .csv
write.csv(tbl_df_dom_summary_MMS,
          file = "outputs/PSA/Modified Model Specification/summary_dom_MMS.csv",
          row.names = FALSE)

# TS
tbl_df_dom_summary_TS <- summarise_all(tbl_df_dominant_sim_TS, mean)

######################################
### Produce scatter plot for ICERs ###
######################################
## Modified Model Specification ##
# 6-month
# Total
inc_qalys_6mo <- df_incremental_PSA_MMS[, "n_inc_qalys_TOTAL_6mo"]
inc_costs_6mo <- df_incremental_PSA_MMS[, "n_inc_costs_TOTAL_6mo"]

plot_PSA_MMS_6mo_scatter <- ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_6mo, y = inc_costs_6mo)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_MMS_6mo_scatter, 
       filename = "Plots/PSA/PSA-MMS-6mo.png", 
       width = 7, height = 10)

# Health Sector
inc_costs_6mo_health_sector <- df_incremental_PSA_MMS[, "n_inc_costs_HEALTH_SECTOR_6mo"]

plot_PSA_MMS_6mo_scatter_health_sector <- ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_6mo, y = inc_costs_6mo_health_sector)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_MMS_6mo_scatter_health_sector, 
       filename = "Plots/PSA/PSA-MMS-6mo-Health-Sector.png", 
       width = 7, height = 10)

# 10-year
inc_qalys_10yr <- df_incremental_PSA_MMS[, "n_inc_qalys_TOTAL_10yr"]
inc_costs_10yr <- df_incremental_PSA_MMS[, "n_inc_costs_TOTAL_10yr"]

plot_PSA_MMS_10yr_scatter <- ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_10yr, y = inc_costs_10yr)) +
 geom_point() +
 geom_hline(yintercept = 0) +
 geom_vline(xintercept = 0) +
 geom_abline(slope = 100000, intercept = 0)

 ggsave(plot_PSA_MMS_10yr_scatter,
        filename = "Plots/PSA/PSA-MMS-10yr.png",
        width = 7, height = 10)

# Health Sector
 inc_costs_10yr_health_sector <- df_incremental_PSA_MMS[, "n_inc_costs_HEALTH_SECTOR_10yr"]

 plot_PSA_MMS_10yr_scatter_health_sector <- ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_10yr, y = inc_costs_10yr_health_sector)) +
   geom_point() +
   geom_hline(yintercept = 0) +
   geom_vline(xintercept = 0) +
   geom_abline(slope = 100000, intercept = 0)

 ggsave(plot_PSA_MMS_10yr_scatter_health_sector,
        filename = "Plots/PSA/PSA-MMS-10yr-Health-Sector.png",
        width = 7, height = 10)

# Lifetime
inc_qalys_lifetime <- df_incremental_PSA_MMS[, "n_inc_costs_TOTAL_life"]
inc_costs_lifetime <- df_incremental_PSA_MMS[, "n_inc_qalys_TOTAL_life"]

ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_lifetime, y = inc_costs_lifetime)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
xlim(min(a), max(a)) +
ylim(min(b), max(b))

## Trial Specification ##
# 1-year
# Total
inc_qalys_6mo <- df_incremental_PSA_TS[, "n_inc_qalys_TOTAL_6mo"]
inc_costs_6mo <- df_incremental_PSA_TS[, "n_inc_costs_TOTAL_6mo"]

plot_PSA_TS_6mo_scatter <- ggplot(df_incremental_PSA_TS, aes(x = inc_qalys_6mo, y = inc_costs_6mo)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_TS_6mo_scatter,
       filename = "Plots/PSA/PSA-TS-6mo.png",
       width = 7, height = 10)

# Health Sector
inc_costs_6mo_health_sector <- df_incremental_PSA_TS[, "n_inc_costs_HEALTH_SECTOR_6mo"]

plot_PSA_TS_6mo_scatter_health_sector <- ggplot(df_incremental_PSA_TS, aes(x = inc_qalys_6mo, y = inc_costs_6mo_health_sector)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_TS_6mo_scatter_health_sector,
       filename = "Plots/PSA/PSA-TS-6mo-Health-Sector.png",
       width = 7, height = 10)

#####################
### Plot ellipses ###
#####################
### Societal perspective ###
# MMS
df_incremental_PSA_MMS_TOTAL_6mo <- df_incremental_PSA_MMS %>% as_tibble() %>% mutate(inc_qalys_MMS_6mo = n_inc_qalys_TOTAL_6mo,
                                                                                inc_costs_MMS_6mo = n_inc_costs_TOTAL_6mo) %>% select(inc_qalys_MMS_6mo, inc_costs_MMS_6mo)

df_incremental_PSA_MMS_TOTAL_10yr <- df_incremental_PSA_MMS %>% as_tibble() %>% mutate(inc_qalys_MMS_10yr = n_inc_qalys_TOTAL_10yr,
                                                                                 inc_costs_MMS_10yr = n_inc_costs_TOTAL_10yr) %>% select(inc_qalys_MMS_10yr, inc_costs_MMS_10yr)

df_incremental_PSA_MMS_TOTAL_life <- df_incremental_PSA_MMS %>% as_tibble() %>% mutate(inc_qalys_MMS_life = n_inc_qalys_TOTAL_life,
                                                                                inc_costs_MMS_life = n_inc_costs_TOTAL_life) %>% select(inc_qalys_MMS_life, inc_costs_MMS_life)

# TS
# df_incremental_PSA_TS_TOTAL_6mo <- df_incremental_PSA_TS %>% as_tibble() %>% mutate(inc_qalys_TS_6mo = n_inc_qalys_TOTAL_6mo,
#                                                                               inc_costs_TS_6mo = n_inc_costs_TOTAL_6mo) %>% select(inc_qalys_TS_6mo, inc_costs_TS_6mo)

# Combine
df_PSA_ellipse_TOTAL <- cbind(df_incremental_PSA_MMS_TOTAL_6mo, df_incremental_PSA_MMS_TOTAL_10yr, df_incremental_PSA_MMS_TOTAL_life)

### Health sector perspective ###
# MMS
df_incremental_PSA_MMS_HEALTH_SECTOR_6mo <- df_incremental_PSA_MMS %>% as_tibble() %>% mutate(inc_qalys_MMS_6mo = n_inc_qalys_TOTAL_6mo,
                                                                                              inc_costs_MMS_6mo = n_inc_costs_HEALTH_SECTOR_6mo) %>% select(inc_qalys_MMS_6mo, inc_costs_MMS_6mo)

df_incremental_PSA_MMS_HEALTH_SECTOR_10yr <- df_incremental_PSA_MMS %>% as_tibble() %>% mutate(inc_qalys_MMS_10yr = n_inc_qalys_TOTAL_10yr,
                                                                                               inc_costs_MMS_10yr = n_inc_costs_HEALTH_SECTOR_10yr) %>% select(inc_qalys_MMS_10yr, inc_costs_MMS_10yr)

df_incremental_PSA_MMS_HEALTH_SECTOR_life <- df_incremental_PSA_MMS %>% as_tibble() %>% mutate(inc_qalys_MMS_life = n_inc_qalys_TOTAL_life,
                                                                                               inc_costs_MMS_life = n_inc_costs_HEALTH_SECTOR_life) %>% select(inc_qalys_MMS_life, inc_costs_MMS_life)

# TS
# df_incremental_PSA_TS_HEALTH_SECTOR_6mo <- df_incremental_PSA_TS %>% as_tibble() %>% mutate(inc_qalys_TS_6mo = n_inc_qalys_TOTAL_6mo,
#                                                                                             inc_costs_TS_6mo = n_inc_costs_HEALTH_SECTOR_6mo) %>% select(inc_qalys_TS_6mo, inc_costs_TS_6mo)

# Combine
df_PSA_ellipse_HEALTH_SECTOR <- cbind(df_incremental_PSA_MMS_HEALTH_SECTOR_6mo, df_incremental_PSA_MMS_HEALTH_SECTOR_10yr, df_incremental_PSA_MMS_HEALTH_SECTOR_life)

#############
### Plots ###
#############
## Full plot ##
### Societal perspective ###
plot_PSA_ellipse <- ggplot() +
  # Points
  # MMS (Lifetime)
  geom_point(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life, alpha = 0.1), colour = "#08306b", size = 1) +
  #geom_raster(fill = "maroon")+
  
  # MMS (10-year)
  geom_point(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_10yr, y = inc_costs_MMS_10yr, alpha = 0.1), colour = "#4292c6", size = 1) +
  #geom_hex(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_10yr, y = inc_costs_MMS_10yr), bins = 250) +  scale_fill_viridis_c(option = "magma") +
  
  # MMS (6-month)
  geom_point(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo, alpha = 0.1), colour = "#c6dbef", size = 1) +
  #geom_point(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), colour = "#e5f5f9", size = 1.5) +
  #geom_hex(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), bins = 250) + scale_fill_viridis_c(option = "magma") +
  
  # Colours: darkblue: #08306b; med blue: #4292c6; light: #c6dbef
  # Miami dolphins: #005778; #FC4C02; #008E97
  # Balimore: #000000; #241773; #9E7C0C
  # Reds: #67000d; #ef3b2c; #fcbba1
  
  # Ellipses
  # MMS (6-month)
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = 2, color = "#000000", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +
  # MMS (10-year)
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_10yr, y = inc_costs_MMS_10yr), linetype = 2, color = "#000000", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_10yr, y = inc_costs_MMS_10yr), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +
  # MMS (Lifetime)
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = 2, color = "#000000", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +

### Health sector perspective ###
  # Points
  # MMS (Lifetime)
  geom_point(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life, alpha = 0.1), colour = "#67000d", size = 1) +
  
  # MMS (10-year)
  geom_point(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_10yr, y = inc_costs_MMS_10yr, alpha = 0.1), colour = "#ef3b2c", size = 1) +
  
  # MMS (6-month)
  geom_point(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo, alpha = 0.1), colour = "#fcbba1", size = 1) +
  
  # Ellipses
  # MMS (6-month)
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = 2, color = "#000000", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +
  # MMS (10-year)
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_10yr, y = inc_costs_MMS_10yr), linetype = 2, color = "#000000", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_10yr, y = inc_costs_MMS_10yr), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +
  # MMS (Lifetime)
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = 2, color = "#000000", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +

  # Add labels
  #annotate("text", x = 0.05, y = 15000, label = "TS (6-month) \n Societal", fontface = "bold", size = 3) +
  annotate("text", x =  0.023, y = 15000, label = "Societal Perspective \n (6-month)", fontface = "bold", size = 3) +
  annotate("text", x = -0.08, y = 50000, label = "Societal Perspective \n (10-year)", fontface = "bold", size = 3) +
  annotate("text", x = -0.24, y = 27000, label = "Societal Perspective \n (Lifetime)", fontface = "bold", size = 3) +
  
  #annotate("text", x = 0.06, y = 3000, label = "TS (6-month) \n Health Sector", fontface = "bold", size = 3, color = "royalblue") +
  annotate("text", x = 0.03, y = -7000, label = "Health-Sector Perspective \n (6-month)", fontface = "bold", size = 3, color = "royalblue") +
  annotate("text", x = -0.05, y = -10000, label = "Health-Sector Perspective \n (10-year)", fontface = "bold", size = 3, color = "royalblue") +
  annotate("text", x =  -0.2, y = -10000, label = "Health-Sector Perspective \n (Lifetime)", fontface = "bold", size = 3, color = "royalblue") +
  
  annotate("text", x =  -0.38, y = -48000, label = "ICER = $100,000/QALY", size = 3) +
  
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 1.0) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1.0) +
  geom_abline(slope = 100000, intercept = 0) +
  labs(y = "Incremental costs (2020 CAD)", x = "Incremental QALYs") +
  xlim(-0.4, 0.1) +
  #ylim(-40000, 75000) +
  scale_y_continuous(labels = scales::dollar_format(scale = .001, suffix = "K"), limits = c(-80000, 75000)) +
  
  # Try this part for adding legend
  # scale_fill_manual(values = c("red", "green", "blue"), name = "My name", 
  #                   guide = guide_legend(reverse = TRUE)) +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        plot.title = element_text(hjust=0.02, vjust=-7))#, legend.position = "none")

plot_PSA_ellipse

ggsave(plot_PSA_ellipse, 
       filename = "Plots/PSA/PSA-Ellipse.png", 
       width = 10, height = 7)

## Zoomed in plot ##

plot_PSA_ellipse_zoom <- ggplot() +
  # Points
  # MMS (6-month)
  geom_point(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), colour = "wheat", size = 1.5) +
  # TS
  geom_point(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_TS_6mo, y = inc_costs_TS_6mo), colour = "lightblue", size = 1.5) +
  # MMS (6-month)
  geom_point(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), colour = "tomato", size = 1.5) +
  # TS
  geom_point(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_TS_6mo, y = inc_costs_TS_6mo), colour = "plum1", size = 1.5) +
  
  # MMS (6-month)
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = 2, color = "wheat4", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = "solid", color = "wheat3", size = 1, alpha = 1, level = 0.5) +
  # TS
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_TS_6mo, y = inc_costs_TS_6mo), linetype = 2, color = "navyblue", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_TS_6mo, y = inc_costs_TS_6mo), linetype = "solid", color = "royalblue", size = 1, alpha = 1, level = 0.5) +
  # MMS (6-month)
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = 2, color = "firebrick1", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = "solid", color = "firebrick4", size = 1, alpha = 1, level = 0.5) +
  # TS
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_TS_6mo, y = inc_costs_TS_6mo), linetype = 2, color = "purple4", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_TS_6mo, y = inc_costs_TS_6mo), linetype = "solid", color = "mediumpurple", size = 1, alpha = 1, level = 0.5) +
  
  # Add labels
  annotate("text", x = 0.01, y = 7000, label = "TS (6-month) \n Societal", fontface = "bold", size = 3) +
  annotate("text", x =  -0.012, y = 22000, label = "MMS (6-month) \n Societal", fontface = "bold", size = 3) +
  
  annotate("text", x = 0.02, y = -3000, label = "TS (6-month) \n Health Sector", fontface = "bold", size = 3, color = "royalblue") +
  annotate("text", x = -0.02, y = -5000, label = "MMS (6-month) \n Health Sector", fontface = "bold", size = 3, color = "royalblue") +
  
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 1.0) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1.0) +
  geom_abline(slope = 100000, intercept = 0) +
  labs(y = "Incremental costs", x = "Incremental QALYs") +
  #xlim(-0.4, 0.1) +
  #ylim(-50000, 75000) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_rect(fill = "transparent", colour = "transparent"),
        plot.title=element_text(hjust=0.02, vjust=-7), legend.justification=c(1,0), legend.position=c(1,0))
  
  #plot_PSA_ellipse_zoom <- plot_PSA_ellipse + xlim(-0.025, 0.025) + ylim(-25000, 25000)

plot_PSA_ellipse_zoom

ggsave(plot_PSA_ellipse_zoom, 
       filename = "Plots/PSA/PSA-Ellipse-zoom.png", 
       width = 10, height = 7)

### Produce CEAC plot from cost-effectiveness results
## Prep data
# Set number of points (one for each quantile)
# quantiles <- seq(0,1,0.01)
# 
# tbl_df_CEAC <- df_ICER_PSA %>% as.tibble() %>% gather("variable", "value") %>% group_by(variable) %>% summarise_at(vars('value'), 
#                                                                                                                    funs(quantile = list(as.tibble(as.list(quantile(., probs = quantiles)))))) %>% unnest
# 
# ## Create plot
# #CEAC <- function(scenario = scenario){
# # Subset by ICER (e.g. societal-lifetime)
# tbl_df_CEAC_subset <- tbl_df_CEAC %>% filter(variable == "ICER (1-year)") %>% select(-variable)
# tbl_df_CEAC_subset <- gather(tbl_df_CEAC_subset)
# 
# len <- nrow(tbl_df_CEAC_subset)
# increment <- 1/(len - 1)
# p <- seq(0, 1, increment)
# 
# # Final data for plotting
# tbl_df_CEAC_plot <- cbind(tbl_df_CEAC_subset, p)
# 
# # Create plot
# ggplot(data = tbl_df_CEAC_plot, aes(x = value, y = p)) +
#   geom_line() + xlim(0, NA)
# #}
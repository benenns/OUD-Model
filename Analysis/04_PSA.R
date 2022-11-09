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
#library(RPushbullet)

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
n_sim <- 10000 # just to test function (will be set as n_sim)

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
#n_sim <- 1500
n_block_size <- 500 # size of block for each loop
n_blocks <- n_sim/n_block_size
n_start <- 0 # set to 0 if running full PSA

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
#pbPost("note", "Notification from R", "PSA runs complete (office desktop)")

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

### Process PSA results
## Read-in saved results
## Modified Model Specification
load(file = "outputs/PSA/Modified Model Specification/outcomes_MET_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/outcomes_BUP_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/incremental_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/ICER_PSA_MMS.RData")

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
tbl_df_summary_MET_MMS <- df_outcomes_MET_PSA_MMS_comb %>% as.tibble() %>% select(n_TOTAL_costs_6mo, 
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
tbl_df_summary_BUP_MMS <- df_outcomes_BUP_PSA_MMS_comb %>% as.tibble() %>% select(n_TOTAL_costs_6mo, 
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
tbl_df_summary_incremental_MMS <- df_incremental_PSA_MMS_comb %>% as.tibble() %>% select(n_inc_costs_TOTAL_6mo,
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
tbl_df_summary_ICER_MMS <- df_ICER_PSA_MMS_comb %>% as.tibble() %>% select(n_icer_TOTAL_6mo,
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
tbl_df_summary_ICER_MMS_2000 <- df_ICER_PSA_MMS_comb[1:2000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_2500 <- df_ICER_PSA_MMS_comb[1:2500, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_3000 <- df_ICER_PSA_MMS_comb[1:3000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_3500 <- df_ICER_PSA_MMS_comb[1:3500, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_4000 <- df_ICER_PSA_MMS_comb[1:4000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_4500 <- df_ICER_PSA_MMS_comb[1:4500, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_5000 <- df_ICER_PSA_MMS_comb[1:5000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_5500 <- df_ICER_PSA_MMS_comb[1:5500, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_6000 <- df_ICER_PSA_MMS_comb[1:6000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_6500 <- df_ICER_PSA_MMS_comb[1:6500, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_7000 <- df_ICER_PSA_MMS_comb[1:7000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_7500 <- df_ICER_PSA_MMS_comb[1:7500, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_8000 <- df_ICER_PSA_MMS_comb[1:8000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_8500 <- df_ICER_PSA_MMS_comb[1:8500, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_9000 <- df_ICER_PSA_MMS_comb[1:9000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_9500 <- df_ICER_PSA_MMS_comb[1:9500, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_summary_ICER_MMS_10000 <- df_ICER_PSA_MMS_comb[1:10000, ] %>% as.tibble() %>% select(n_icer_TOTAL_life) %>%
  gather("variable", "value") %>% 
  group_by(variable) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            q50 = quantile(value, probs = .5),
            q025 = quantile(value, probs = .025),
            q975 = quantile(value, probs = .975),
            min = min(value),
            max = max(value))

tbl_df_CI_stability <- rbind(tbl_df_summary_ICER_MMS_2000, tbl_df_summary_ICER_MMS_2500, tbl_df_summary_ICER_MMS_3000, tbl_df_summary_ICER_MMS_3500, tbl_df_summary_ICER_MMS_4000, tbl_df_summary_ICER_MMS_4500, tbl_df_summary_ICER_MMS_5000,
                             tbl_df_summary_ICER_MMS_5500, tbl_df_summary_ICER_MMS_6000, tbl_df_summary_ICER_MMS_6500, tbl_df_summary_ICER_MMS_7000, tbl_df_summary_ICER_MMS_7500, tbl_df_summary_ICER_MMS_8000, tbl_df_summary_ICER_MMS_8500, tbl_df_summary_ICER_MMS_9000,
                             tbl_df_summary_ICER_MMS_9500, tbl_df_summary_ICER_MMS_10000)
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

#df_PSA_summary <- as.data.frame()

## Incremental ##
# MMS
tbl_df_summary_inc_MMS <- df_incremental_PSA_MMS_comb %>% as.tibble() %>% mutate(BNX_dominate_TOTAL_6mo = ifelse(n_inc_costs_TOTAL_6mo < 0 & n_inc_qalys_TOTAL_6mo > 0, 1, 0),
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

tbl_df_dominant_sim_MMS <- tbl_df_summary_inc_MMS %>% select(BNX_dominate_TOTAL_6mo, BNX_dominate_TOTAL_10yr, BNX_dominate_TOTAL_life, BNX_dominate_HEALTH_SECTOR_6mo, BNX_dominate_HEALTH_SECTOR_10yr, BNX_dominate_HEALTH_SECTOR_life,
                                                             MET_dominate_TOTAL_6mo, MET_dominate_TOTAL_10yr, MET_dominate_TOTAL_life, MET_dominate_HEALTH_SECTOR_6mo, MET_dominate_HEALTH_SECTOR_10yr, MET_dominate_HEALTH_SECTOR_life)

tbl_df_dom_summary_MMS <- summarise_all(tbl_df_dominant_sim_MMS, mean)

write.csv(tbl_df_dom_summary_MMS,
          file = "outputs/PSA/Modified Model Specification/summary_dom_PSA_MMS.csv",
          row.names = FALSE)

######################################
### Produce scatter plot for ICERs ###
######################################
## Modified Model Specification ##
# 6-month
# Total
# inc_qalys_6mo <- df_incremental_PSA_MMS[, "n_inc_qalys_TOTAL_6mo"]
# inc_costs_6mo <- df_incremental_PSA_MMS[, "n_inc_costs_TOTAL_6mo"]
# 
# plot_PSA_MMS_6mo_scatter <- ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_6mo, y = inc_costs_6mo)) +
#   geom_point() +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   geom_abline(slope = 100000, intercept = 0)
# 
# ggsave(plot_PSA_MMS_6mo_scatter, 
#        filename = "Plots/PSA/PSA-MMS-6mo.png", 
#        width = 7, height = 10)
# 
# # Health Sector
# inc_costs_6mo_health_sector <- df_incremental_PSA_MMS[, "n_inc_costs_HEALTH_SECTOR_6mo"]
# 
# plot_PSA_MMS_6mo_scatter_health_sector <- ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_6mo, y = inc_costs_6mo_health_sector)) +
#   geom_point() +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   geom_abline(slope = 100000, intercept = 0)
# 
# ggsave(plot_PSA_MMS_6mo_scatter_health_sector, 
#        filename = "Plots/PSA/PSA-MMS-6mo-Health-Sector.png", 
#        width = 7, height = 10)
# 
# # Lifetime
# inc_qalys_lifetime <- df_incremental_PSA_MMS[, "n_inc_costs_TOTAL_life"]
# inc_costs_lifetime <- df_incremental_PSA_MMS[, "n_inc_qalys_TOTAL_life"]
# 
# ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_lifetime, y = inc_costs_lifetime)) +
#   geom_point() +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0)
# xlim(min(a), max(a)) +
# ylim(min(b), max(b))

#####################
### Plot ellipses ###
#####################
### Societal perspective ###
# MMS
df_incremental_PSA_MMS_TOTAL_6mo <- df_incremental_PSA_MMS_comb %>% as_tibble() %>% mutate(inc_qalys_MMS_6mo = n_inc_qalys_TOTAL_6mo,
                                                                                inc_costs_MMS_6mo = n_inc_costs_TOTAL_6mo) %>% select(inc_qalys_MMS_6mo, inc_costs_MMS_6mo)

# df_incremental_PSA_MMS_TOTAL_10yr <- df_incremental_PSA_MMS %>% as_tibble() %>% mutate(inc_qalys_MMS_10yr = n_inc_qalys_TOTAL_10yr,
#                                                                                  inc_costs_MMS_10yr = n_inc_costs_TOTAL_10yr) %>% select(inc_qalys_MMS_10yr, inc_costs_MMS_10yr)

df_incremental_PSA_MMS_TOTAL_life <- df_incremental_PSA_MMS_comb %>% as_tibble() %>% mutate(inc_qalys_MMS_life = n_inc_qalys_TOTAL_life,
                                                                                inc_costs_MMS_life = n_inc_costs_TOTAL_life) %>% select(inc_qalys_MMS_life, inc_costs_MMS_life)


# Combine
df_PSA_ellipse_TOTAL <- cbind(df_incremental_PSA_MMS_TOTAL_6mo, df_incremental_PSA_MMS_TOTAL_life)
df_PSA_ellipse_TOTAL <- df_PSA_ellipse_TOTAL %>% mutate(Scenario = "Societal Perspective")

### Health sector perspective ###
# MMS
df_incremental_PSA_MMS_HEALTH_SECTOR_6mo <- df_incremental_PSA_MMS_comb %>% as_tibble() %>% mutate(inc_qalys_MMS_6mo = n_inc_qalys_TOTAL_6mo,
                                                                                              inc_costs_MMS_6mo = n_inc_costs_HEALTH_SECTOR_6mo) %>% select(inc_qalys_MMS_6mo, inc_costs_MMS_6mo)

# df_incremental_PSA_MMS_HEALTH_SECTOR_10yr <- df_incremental_PSA_MMS %>% as_tibble() %>% mutate(inc_qalys_MMS_10yr = n_inc_qalys_TOTAL_10yr,
#                                                                                                inc_costs_MMS_10yr = n_inc_costs_HEALTH_SECTOR_10yr) %>% select(inc_qalys_MMS_10yr, inc_costs_MMS_10yr)

df_incremental_PSA_MMS_HEALTH_SECTOR_life <- df_incremental_PSA_MMS_comb %>% as_tibble() %>% mutate(inc_qalys_MMS_life = n_inc_qalys_TOTAL_life,
                                                                                               inc_costs_MMS_life = n_inc_costs_HEALTH_SECTOR_life) %>% select(inc_qalys_MMS_life, inc_costs_MMS_life)

# Combine
df_PSA_ellipse_HEALTH_SECTOR <- cbind(df_incremental_PSA_MMS_HEALTH_SECTOR_6mo, df_incremental_PSA_MMS_HEALTH_SECTOR_life)
df_PSA_ellipse_HEALTH_SECTOR <- df_PSA_ellipse_HEALTH_SECTOR %>% mutate(Scenario = "Health Sector Perspective")

# Combine all
df_PSA_ellipse <- rbind(df_PSA_ellipse_TOTAL, df_PSA_ellipse_HEALTH_SECTOR)
df_PSA_points_temp <- df_PSA_ellipse %>% as_tibble() %>% rename(qalys.6mo = inc_qalys_MMS_6mo,
                                                                qalys.life = inc_qalys_MMS_life,
                                                                costs.6mo = inc_costs_MMS_6mo,
                                                                costs.life = inc_costs_MMS_life) %>% mutate(ID = row_number())

df_PSA_points_qalys <- df_PSA_points_temp %>% select(ID, Scenario, qalys.6mo, qalys.life)
df_PSA_points_qalys_long <- reshape(df_PSA_points_qalys, direction = 'long', 
                         varying = c('qalys.6mo', 'qalys.life'), 
                         timevar = 'var',
                         times = c('6mo', 'life'),
                         v.names = 'qalys',
                         idvar = c('ID', 'Scenario'))

df_PSA_points_costs <- df_PSA_points_temp %>% select(ID, Scenario, costs.6mo, costs.life)
df_PSA_points_costs_long <- reshape(df_PSA_points_costs, direction = 'long', 
                                    varying = c('costs.6mo', 'costs.life'), 
                                    timevar = 'var',
                                    times = c('6mo', 'life'),
                                    v.names = 'costs',
                                    idvar = c('ID', 'Scenario'))

df_PSA_points <- bind_cols(df_PSA_points_qalys_long, df_PSA_points_costs_long, .name_repair = "minimal")

df_PSA_points <- inner_join(df_PSA_points_qalys_long, df_PSA_points_costs_long, by = c('ID', 'var', 'Scenario')) %>% 
  mutate(index = ifelse(Scenario == "Societal Perspective" & var == "6mo", "Societal (6-month)", 
                 ifelse(Scenario == "Societal Perspective" & var == "life", "Societal (Lifetime)", 
                 ifelse (Scenario == "Health Sector Perspective" & var == "6mo", "Health Sector (6-month)", "Health Sector (Lifetime)"))))
  
#############
### Plots ###
#############
## Full plot ##

plot_PSA_ellipse <- ggplot() +
  # Points (all scenarios and time horizons)
  geom_point(data = df_PSA_points, aes(x = qalys, y = costs, colour = index), alpha = 0.4, size = 1) +

  # Ellipses (Societal)
  # MMS (6-month)
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = 2, color = "white", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +
  # MMS (Lifetime)
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = 2, color = "#000000", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +

  # Ellipses (Health sector)
  # MMS (6-month)
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = 2, color = "white", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +
  # MMS (Lifetime)
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = 2, color = "#000000", size = 1, alpha = 1, level = 0.95) +
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +

  # Add labels
  #annotate("text", x =  0.05, y = 25000, label = "Six-month \n Time-horizon", fontface = "bold", size = 3) +
  #annotate("text", x = -0.15, y = 45000, label = "Lifetime \n Time-horizon", fontface = "bold", size = 3) +
  
  annotate("text", x =  -0.47, y = -60000, label = "ICER: \n $100,000/QALY", size = 4) +
  
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 1.0) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1.0) +
  geom_abline(slope = 100000, intercept = 0) +
  labs(y = "Incremental costs (BNX vs. MET)", x = "Incremental QALYs (BNX vs. MET)") +
  xlim(-0.5, 0.1) +
  scale_y_continuous(labels = scales::dollar_format(scale = .001, suffix = "K"), limits = c(-100000, 100000)) +
  
  scale_color_manual(name = '',
                     breaks = c('Societal (6-month)', 'Health Sector (6-month)', 'Societal (Lifetime)', 'Health Sector (Lifetime)'),
                     values = c('Societal (6-month)' = "#313695", 'Health Sector (6-month)' = "#f46d43", 'Societal (Lifetime)' = "#2166ac", 'Health Sector (Lifetime)' = "#d7191c")) +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        plot.title = element_text(hjust=0.02, vjust=-7), 
        legend.position = "bottom",
        text = element_text(size = 15))

plot_PSA_ellipse

# Output full plot
ggsave(plot_PSA_ellipse, 
       filename = "Plots/PSA/PSA-Ellipse.png", 
       width = 10, height = 7, dpi = 350)

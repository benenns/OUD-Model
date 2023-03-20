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
library(future)
library(doFuture)

#### START HERE TO RUN PSA ####

# Set number of cores
#n_cores <- detectCores()
n_cores <- 8 # number of cores is dictated more by memory than cpu (will max out 32GB memory at >8 cores)
makeCluster(n_cores, outfile = "checks/parallel_log.txt")
registerDoParallel(n_cores)

#registerDoFuture()
#cl <- makeCluster(4)
#lan(cluster, workers = cl)

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
# df_psa_params_MMS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_trial, scenario = "MMS",
#                                              file.death_hr = "data/death_hr.csv",
#                                              file.frailty = "data/frailty.csv",
#                                              file.weibull = "data/Modified Model Specification/weibull.csv",
#                                              file.unconditional = "data/Modified Model Specification/unconditional.csv",
#                                              file.overdose = "data/overdose.csv",
#                                              file.fentanyl = "data/fentanyl.csv",
#                                              file.naloxone = "data/naloxone.csv",### R&R MODIFICATION ###
#                                              file.hiv = "data/hiv_sero.csv",
#                                              file.hcv = "data/hcv_sero.csv",
#                                              file.costs = "data/Modified Model Specification/costs.csv",
#                                              file.crime_costs = "data/Modified Model Specification/crime_costs.csv",
#                                              file.qalys = "data/Modified Model Specification/qalys.csv",
#                                              file.imis_output = "outputs/Calibration/imis_output.RData")


# Output data
## As .RData
#save(df_psa_params_MMS, 
#     file = "outputs/PSA/Modified Model Specification/df_psa_params_MMS.RData")

## As .csv
#write.csv(df_psa_params_MMS,"outputs/PSA/Modified Model Specification/input_PSA_MMS.csv", 
#          row.names = TRUE)

# Load PSA inputs
load(file = "outputs/PSA/Modified Model Specification/df_psa_params_MMS.RData")

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

# Run PSA blockwise
n_runs <- n_sim # n_sim to run entire PSA
n_block_size <- 1000 # size of block for each loop
n_blocks <- n_runs/n_block_size #to run entire set
n_start <- 0 # set to 0 if running full PSA

# Initialize lists
l_outcomes_MET_PSA_MMS <- list()
l_outcomes_BUP_PSA_MMS <- list()
l_ICER_PSA_MMS         <- list()        
l_incremental_PSA_MMS  <- list()


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
    
    gc()
    
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
  
  out <- paste0("Block ", (j + 1), " of ", n_blocks, " complete.")  # Status
  print(out)
}

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

#### START HERE IF PSA ALREADY RUN ####
### Process PSA results
## Read-in saved results
## Modified Model Specification
load(file = "outputs/PSA/Modified Model Specification/outcomes_MET_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/outcomes_BUP_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/incremental_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/ICER_PSA_MMS.RData")

# TEMP CODE
# df_outcomes_MET_PSA_MMS_comb <- df_outcomes_MET_PSA_MMS
# df_outcomes_BUP_PSA_MMS_comb <- df_outcomes_BUP_PSA_MMS
# df_incremental_PSA_MMS_comb <- df_incremental_PSA_MMS
# df_ICER_PSA_MMS_comb <- df_ICER_PSA_MMS

### Summary stats ###
# Methadone
tbl_df_summary_MET_MMS <- df_outcomes_MET_PSA_MMS_comb %>% as.tibble() %>% 
  select(n_TOTAL_costs_6mo, 
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
tbl_df_summary_BUP_MMS <- df_outcomes_BUP_PSA_MMS_comb %>% as.tibble() %>% 
  select(n_TOTAL_costs_6mo, 
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
tbl_df_summary_incremental_MMS <- df_incremental_PSA_MMS_comb %>% as.tibble() %>% 
  select(n_inc_costs_TOTAL_6mo,
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
tbl_df_summary_ICER_MMS <- df_ICER_PSA_MMS_comb %>% as.tibble() %>% 
  select(n_icer_TOTAL_6mo,
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

## As .csv ####
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

######################
### CEAC (FOR R&R) ###
### Lifetime       ###
######################
tbl_df_CEAC <- df_incremental_PSA_MMS_comb %>% as_tibble() %>% mutate(ICER_TOTAL_life = (n_inc_costs_TOTAL_life/n_inc_qalys_TOTAL_life),
                                                                      ICER_HEALTH_SECTOR_life = (n_inc_costs_HEALTH_SECTOR_life/n_inc_qalys_TOTAL_life),
                                                                      BNX_dominate_TOTAL_life = ifelse(n_inc_costs_TOTAL_life < 0 & n_inc_qalys_TOTAL_life > 0, 1, 0),
                                                                      BNX_dominate_HEALTH_SECTOR_life = ifelse(n_inc_costs_HEALTH_SECTOR_life < 0 & n_inc_qalys_TOTAL_life > 0, 1, 0),
                                                                      MET_dominate_TOTAL_life = ifelse(n_inc_costs_TOTAL_life > 0 & n_inc_qalys_TOTAL_life < 0, 1, 0),
                                                                      MET_dominate_HEALTH_SECTOR_life = ifelse(n_inc_costs_HEALTH_SECTOR_life > 0 & n_inc_qalys_TOTAL_life < 0, 1, 0),
                                                                      BNX_health_loss = ifelse(n_inc_qalys_TOTAL_life < 0, 1, 0)) %>%
                                                               select(n_inc_costs_TOTAL_life, n_inc_costs_HEALTH_SECTOR_life, n_inc_qalys_TOTAL_life,
                                                                      ICER_TOTAL_life, ICER_HEALTH_SECTOR_life, BNX_dominate_TOTAL_life, BNX_dominate_HEALTH_SECTOR_life, 
                                                                      MET_dominate_TOTAL_life, MET_dominate_HEALTH_SECTOR_life, BNX_health_loss) %>%
                                                               mutate(CE_0K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 0) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_TOTAL_life <= 0) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_10K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 10000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                           ((BNX_health_loss == 0 & ICER_TOTAL_life <= 10000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_20K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 20000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                             ((BNX_health_loss == 0 & ICER_TOTAL_life <= 20000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_30K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 30000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                             ((BNX_health_loss == 0 & ICER_TOTAL_life <= 30000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_40K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 40000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                             ((BNX_health_loss == 0 & ICER_TOTAL_life <= 40000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_50K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 50000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_TOTAL_life <= 50000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_60K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 60000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_TOTAL_life <= 60000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_70K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 70000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_TOTAL_life <= 70000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_80K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 80000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_TOTAL_life <= 80000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_90K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 90000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_TOTAL_life <= 90000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_100K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 100000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_TOTAL_life <= 100000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_110K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 110000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_TOTAL_life <= 110000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_120K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 120000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_TOTAL_life <= 120000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_130K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 130000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_TOTAL_life <= 130000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_140K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 140000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_TOTAL_life <= 140000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_150K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 150000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_TOTAL_life <= 150000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_160K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 160000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_TOTAL_life <= 160000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_170K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 170000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_TOTAL_life <= 170000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_180K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 180000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_TOTAL_life <= 180000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_190K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 190000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_TOTAL_life <= 190000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_200K_TOTAL_life = ifelse(((BNX_health_loss == 1 & ICER_TOTAL_life >= 200000) | BNX_dominate_TOTAL_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_TOTAL_life <= 200000) | BNX_dominate_TOTAL_life == 1), 1, 0),
                                                                      CE_0K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 0) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                  ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 0) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_10K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 10000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 10000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_20K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 20000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 20000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_30K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 30000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 30000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_40K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 40000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 40000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_50K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 50000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 50000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_60K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 60000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 60000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_70K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 70000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 70000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_80K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 80000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 80000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_90K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 90000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                   ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 90000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_100K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 100000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 100000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_110K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 110000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 110000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_120K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 120000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 120000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_130K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 130000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 130000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_140K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 140000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 140000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_150K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 150000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 150000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_160K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 160000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 160000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_170K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 170000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 170000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_180K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 180000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 180000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_190K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 190000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 190000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0),
                                                                      CE_200K_HEALTH_SECTOR_life = ifelse(((BNX_health_loss == 1 & ICER_HEALTH_SECTOR_life >= 200000) | BNX_dominate_HEALTH_SECTOR_life == 1) | 
                                                                                                    ((BNX_health_loss == 0 & ICER_HEALTH_SECTOR_life <= 200000) | BNX_dominate_HEALTH_SECTOR_life == 1), 1, 0))
                                                                      
  # ADD SUMMARY POINTS AT EACH THRESHOLD
tbl_df_CEAC_summary <- tbl_df_CEAC %>% select(CE_0K_TOTAL_life, CE_10K_TOTAL_life, CE_20K_TOTAL_life, CE_30K_TOTAL_life, CE_40K_TOTAL_life, CE_50K_TOTAL_life, CE_60K_TOTAL_life, CE_70K_TOTAL_life,
                                              CE_80K_TOTAL_life, CE_90K_TOTAL_life, CE_100K_TOTAL_life, CE_110K_TOTAL_life, CE_120K_TOTAL_life, CE_130K_TOTAL_life, CE_140K_TOTAL_life,
                                              CE_150K_TOTAL_life, CE_160K_TOTAL_life, CE_170K_TOTAL_life, CE_180K_TOTAL_life, CE_190K_TOTAL_life, CE_200K_TOTAL_life,
                                              CE_0K_HEALTH_SECTOR_life, CE_10K_HEALTH_SECTOR_life, CE_20K_HEALTH_SECTOR_life, CE_30K_HEALTH_SECTOR_life, CE_40K_HEALTH_SECTOR_life, CE_50K_HEALTH_SECTOR_life, CE_60K_HEALTH_SECTOR_life, CE_70K_HEALTH_SECTOR_life,
                                              CE_80K_HEALTH_SECTOR_life, CE_90K_HEALTH_SECTOR_life, CE_100K_HEALTH_SECTOR_life, CE_110K_HEALTH_SECTOR_life, CE_120K_HEALTH_SECTOR_life, CE_130K_HEALTH_SECTOR_life, CE_140K_HEALTH_SECTOR_life,
                                              CE_150K_HEALTH_SECTOR_life, CE_160K_HEALTH_SECTOR_life, CE_170K_HEALTH_SECTOR_life, CE_180K_HEALTH_SECTOR_life, CE_190K_HEALTH_SECTOR_life, CE_200K_HEALTH_SECTOR_life) %>% summarise_all(mean) #%>%
  

tbl_df_CEAC_long <- tbl_df_CEAC_summary %>% pivot_longer(
    cols = starts_with("CE_"),
    names_to = "Scenario",
    values_to = "Proportion",
    values_drop_na = TRUE
  )

tbl_df_labels <- str_split_fixed(tbl_df_CEAC_long$Scenario, '_', 4) %>% as_tibble() %>% mutate(Threshold = ifelse(V2 == "0K", 0, 
                                                                                               ifelse(V2 == "10K", 10000,
                                                                                               ifelse(V2 == "20K", 20000,
                                                                                               ifelse(V2 == "30K", 30000,
                                                                                               ifelse(V2 == "40K", 40000,
                                                                                               ifelse(V2 == "50K", 50000,
                                                                                               ifelse(V2 == "60K", 60000,
                                                                                               ifelse(V2 == "70K", 70000,
                                                                                               ifelse(V2 == "80K", 80000,
                                                                                               ifelse(V2 == "90K", 90000,
                                                                                               ifelse(V2 == "100K", 100000,
                                                                                               ifelse(V2 == "110K", 110000,
                                                                                               ifelse(V2 == "120K", 120000,
                                                                                               ifelse(V2 == "130K", 130000,
                                                                                               ifelse(V2 == "140K", 140000,
                                                                                               ifelse(V2 == "150K", 150000,
                                                                                               ifelse(V2 == "160K", 160000,
                                                                                               ifelse(V2 == "170K", 170000,
                                                                                               ifelse(V2 == "180K", 180000,
                                                                                               ifelse(V2 == "190K", 190000,
                                                                                               ifelse(V2 == "200K", 200000, NA))))))))))))))))))))),
                                                                                   Perspective = ifelse(V3 == "TOTAL", "Societal", "Health Sector")) %>% select(Threshold, Perspective)

tbl_df_CEAC_BNX <- cbind(tbl_df_labels, tbl_df_CEAC_long) %>% as_tibble() %>% mutate(Treatment = "BNX",
                                                                                     tx_perspective = ifelse(Perspective == "Societal", "BNX (Societal)",
                                                                                                                      ifelse(Perspective == "Health Sector", "BNX (Health Sector)", NA))) %>%
                                                                              select(-Scenario) %>% relocate("Proportion", .after = "Treatment")

  
tbl_df_CEAC_MET <- cbind(tbl_df_labels, tbl_df_CEAC_long) %>% as_tibble() %>% mutate(Treatment = "MET",
                                                                                     Proportion_MET = (1 - Proportion)) %>% select(-Proportion) %>%
                                                                              rename(Proportion = Proportion_MET) %>%
                                                                              mutate(tx_perspective = ifelse(Perspective == "Societal", "MET (Societal)",
                                                                                                             ifelse(Perspective == "Health Sector", "MET (Health Sector)", NA))) %>%
                                                                              select(-Scenario)


tbl_df_CEAC_plot <- rbind(tbl_df_CEAC_BNX, tbl_df_CEAC_MET)

# tbl_df_CEAC_temp2 <- tbl_df_CEAC_temp1 %>% mutate(p_BNX_soc = ifelse(tx_perspective == "BNX (Societal)", Proportion, NA),
#                                                p_BNX_hs = ifelse(tx_perspective == "BNX (Health Sector)", Proportion, NA),
#                                                p_MET_soc = ifelse(tx_perspective == "MET (Societal)", Proportion, NA),
#                                                p_MET_hs = ifelse(tx_perspective == "MET (Health Sector)", Proportion, NA))
# 
# tbl_df_CEAC_BNX_soc <- tbl_df_CEAC_temp2 %>% select(Threshold, p_BNX_soc) %>% na.omit()
# tbl_df_CEAC_BNX_hs  <- tbl_df_CEAC_temp2 %>% select(Threshold, p_BNX_hs) %>% na.omit()
# tbl_df_CEAC_MET_soc <- tbl_df_CEAC_temp2 %>% select(Threshold, p_MET_soc) %>% na.omit()
# tbl_df_CEAC_MET_hs  <- tbl_df_CEAC_temp2 %>% select(Threshold, p_MET_hs) %>% na.omit()
# 
# tbl_df_CEAC_temp3 <- merge(x = tbl_df_CEAC_BNX_soc, y = tbl_df_CEAC_MET_soc, 
#                            by = "Threshold", all.x = TRUE)
# tbl_df_CEAC_temp4 <- merge(x = tbl_df_CEAC_temp3, y = tbl_df_CEAC_BNX_hs,
#                            by = "Threshold", all.x = TRUE)
# tbl_df_CEAC_plot  <- merge(x = tbl_df_CEAC_temp4, y = tbl_df_CEAC_MET_hs,
#                            by = "Threshold", all.x = TRUE)
  

#################
### Plot CEAC ###
#################
plot_CEAC <- ggplot(data = tbl_df_CEAC_plot, aes(x = Threshold, y = Proportion)) +
  geom_line(aes(linetype = Perspective, color = Treatment)) +
  scale_linetype_manual(values = c("twodash", "solid"))+
  scale_color_manual(values = c('#999999','#E69F00'))+
  scale_x_continuous(labels = scales::dollar_format(scale = .001, suffix = "K")) +
  labs(
    x = "Willingness to pay threshold ($/QALY)",
    y = "Probability cost-effective",
    color = "Treatment Type",
    linetype = "Payer Perspective"
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        plot.title = element_text(hjust = 0.02, vjust = -7), 
        legend.position = "right",
        text = element_text(size = 15))

# Output full plot
ggsave(plot_CEAC, 
       filename = "Plots/PSA/CEAC.png", 
       width = 8, height = 6, dpi = 350)

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

df_incremental_PSA_MMS_TOTAL_life <- df_incremental_PSA_MMS_comb %>% as_tibble() %>% mutate(inc_qalys_MMS_life = n_inc_qalys_TOTAL_life,
                                                                                inc_costs_MMS_life = n_inc_costs_TOTAL_life) %>% select(inc_qalys_MMS_life, inc_costs_MMS_life)


# Combine
df_PSA_ellipse_TOTAL <- cbind(df_incremental_PSA_MMS_TOTAL_6mo, df_incremental_PSA_MMS_TOTAL_life)
df_PSA_ellipse_TOTAL <- df_PSA_ellipse_TOTAL %>% mutate(Scenario = "Societal Perspective")

### Health sector perspective ###
# MMS
df_incremental_PSA_MMS_HEALTH_SECTOR_6mo <- df_incremental_PSA_MMS_comb %>% as_tibble() %>% mutate(inc_qalys_MMS_6mo = n_inc_qalys_TOTAL_6mo,
                                                                                              inc_costs_MMS_6mo = n_inc_costs_HEALTH_SECTOR_6mo) %>% select(inc_qalys_MMS_6mo, inc_costs_MMS_6mo)

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
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = 2, color = "gold", size = 1, alpha = 1, level = 0.95) +
  #stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +
  # MMS (Lifetime)
  stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = 2, color = "black", size = 1, alpha = 1, level = 0.95) +
  #stat_ellipse(data = df_PSA_ellipse_TOTAL, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +

  # Ellipses (Health sector)
  # MMS (6-month)
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = 2, color = "gold", size = 1, alpha = 1, level = 0.95) +
  #stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_6mo, y = inc_costs_MMS_6mo), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +
  # MMS (Lifetime)
  stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = 2, color = "black", size = 1, alpha = 1, level = 0.95) +
  #stat_ellipse(data = df_PSA_ellipse_HEALTH_SECTOR, aes(x = inc_qalys_MMS_life, y = inc_costs_MMS_life), linetype = "solid", color = "#869397", size = 1, alpha = 1, level = 0.5) +

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
                     ### UPDATE COLORS FOR R&R ###
                     values = c('Societal (6-month)' = "#f46d43", 'Health Sector (6-month)' = "#d7191c", 'Societal (Lifetime)' = "#74add1", 'Health Sector (Lifetime)' = "#2c7bb6")) +
                     #values = c('Societal (6-month)' = "#d7191c", 'Health Sector (6-month)' = "#fdae61", 'Societal (Lifetime)' = "#2c7bb6", 'Health Sector (Lifetime)' = "#abd9e9")) +
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        plot.title = element_text(hjust=0.02, vjust=-7), 
        legend.position = "none",
        text = element_text(size = 15))

plot_PSA_ellipse

# Output full plot
ggsave(plot_PSA_ellipse, 
       filename = "Plots/PSA/PSA-Ellipse.png", 
       width = 10, height = 7, dpi = 350)

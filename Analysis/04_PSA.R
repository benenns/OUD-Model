rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(rBeta2009)
#library(parallel)
library(foreach)
library(doParallel)

# Set number of cores
n_cores <- detectCores()
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
n_sim <- 50 # just to test function (will be set as n_sim)

### PSA model outputs
### Run Markov model for PSA draws and return outputs ###
# Generate PSA parameter draws

######################################
#### Modified Model Specification ####
######################################
df_psa_params_MMS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_cohort, scenario = "MMS",
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


#############################
#### Trial Specification ####
#############################
# Need to draw parameters for both scenarios due to differences in allowed transitions (Dirichlet)
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
                                            file.costs = "data/Trial Specification/costs.csv",
                                            file.crime_costs = "data/Trial Specification/crime_costs.csv",
                                            file.qalys = "data/Trial Specification/qalys.csv",
                                            file.imis_output = "outputs/Calibration/imis_output.RData")

# Extract non-state-exit parameters from BUP
df_psa_params_TS <- df_psa_params_BUP_TS %>% select(-c(p_BUP_BUPC_NI:p_ODN_REL_INJ))

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
                                            file.costs = "data/Trial Specification/costs.csv",
                                            file.crime_costs = "data/Trial Specification/crime_costs.csv",
                                            file.qalys = "data/Trial Specification/qalys.csv",
                                            file.imis_output = "outputs/Calibration/imis_output.RData")

# Extract MET state-exit
df_psa_params_MET_TS_UP <- df_psa_params_MET_TS %>% select(c(p_BUP_BUPC_NI:p_ODN_REL_INJ))

# Add MET state-exit to overall
df_psa_params_MET_TS <- bind_cols(df_psa_params_TS, df_psa_params_MET_TS_UP)

################################
#### Original Specification ####
################################
df_psa_params_OS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_cohort, scenario = "OS",
                                            file.death_hr = "data/death_hr.csv",
                                            file.frailty = "data/frailty.csv",
                                            file.weibull = "data/Original Specification/weibull.csv",
                                            file.unconditional = "data/Original Specification/unconditional.csv",
                                            file.overdose = "data/overdose.csv",
                                            file.fentanyl = "data/fentanyl.csv",
                                            file.hiv = "data/hiv_sero.csv",
                                            file.hcv = "data/hcv_sero.csv",
                                            file.costs = "data/Original Specification/costs.csv",
                                            file.crime_costs = "data/Original Specification/crime_costs.csv",
                                            file.qalys = "data/Original Specification/qalys.csv",
                                            file.imis_output = "outputs/Calibration/imis_output.RData")


### Run decision model on each parameter set of PSA input dataset to produce
### PSA outputs for cost and effects
#n_sim <- 500 # just to test function (will be set as n_sim)

# Initialize data frames
# Modified Model Specification
df_outcomes_MET_PSA_MMS <- data.frame()
df_outcomes_BUP_PSA_MMS <- data.frame()
df_incremental_PSA_MMS <- data.frame()
df_ICER_PSA_MMS <- data.frame()

# Trial Specification
df_outcomes_MET_PSA_TS <- data.frame()
df_outcomes_BUP_PSA_TS <- data.frame()
df_incremental_PSA_TS <- data.frame()
df_ICER_PSA_TS <- data.frame()

# Original Specification
df_outcomes_MET_PSA_OS <- data.frame()
df_outcomes_BUP_PSA_OS <- data.frame()
df_incremental_PSA_OS <- data.frame()
df_ICER_PSA_OS <- data.frame()

#v_outcomes_names <- c("Total Costs (1-year)", "Total Costs (5-year)", "Total Costs (10-year)", "Total Costs (Lifetime)", "Health Sector Costs (1-year)", "Health Sector Costs (5-year)", "Health Sector Costs (10-year)", "Health Sector Costs (Lifetime)",
#                      "Criminal Costs (1-year)", "Criminal Costs (5-year)", "Criminal Costs (10-year)", "Criminal Costs (Lifetime)", "Treatment Costs (1-year)", "Treatment Costs (5-year)", "Treatment Costs (10-year)", "Treatment Costs (Lifetime)",
#                      "Total QALYs (1-year)", "Total QALYs (5-year)", "Total QALYs (10-year)", "Total QALYs (Lifetime)")
#v_ICER_names <- c("ICER (1-year)", "ICER (5-year)", "ICER (10-year)", "ICER (Lifetime)",
#                  "ICER (Health Sector 1-year)", "ICER (Health Sector 5-year)", "ICER (Health Sector 10-year)", "ICER (Health Sector Lifetime)",
#                  "Incremental Costs (1-year)", "Incremental QALYs (1-year)", "Incremental Costs (5-year)", "Incremental QALYs (5-year)", "Incremental Costs (10-year)", "Incremental QALYs (10-year)", "Incremental Costs (Lifetime)", "Incremental QALYs (Lifetime)",
#                  "Incremental Costs (Health Sector 1-year)", "Incremental QALYs (Health Sector 1-year)", "Incremental Costs (Health Sector 5-year)", "Incremental QALYs (Health Sector 5-year)", "Incremental Costs (Health Sector 10-year)", "Incremental QALYs (Health Sector 10-year)", "Incremental Costs (Health Sector Lifetime)", "Incremental QALYs (Health Sector Lifetime)")

######################################
#### Modified Model Specification ####
######################################
#foreach(i = 1:n_sim, .combine = c) %dopar% { # i <- 1
for (i in 1:n_sim){
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
}
#stopImplicitCluster()

### Output results
## As .RData
save(df_outcomes_MET_PSA_MMS, 
     file = "outputs/PSA/Modified Model Specification/outcomes_MET_PSA_MMS.RData")
save(df_outcomes_BUP_PSA_MMS, 
     file = "outputs/PSA/Modified Model Specification/outcomes_BUP_PSA_MMS.RData")
save(df_ICER_PSA_MMS, 
     file = "outputs/PSA/Modified Model Specification/ICER_PSA_MMS.RData")
save(df_incremental_PSA_MMS,
     file = "outputs/PSA/Modified Model Specification/incremental_PSA_MMS.RData")
## As .csv
write.csv(df_outcomes_MET_PSA_MMS, 
          file = "outputs/PSA/Modified Model Specification/outcomes_MET_PSA_MMS.csv",
          row.names = FALSE)
write.csv(df_outcomes_BUP_PSA_MMS, 
          file = "outputs/PSA/Modified Model Specification/outcomes_BUP_PSA_MMS.csv",
          row.names = FALSE)
write.csv(df_ICER_PSA_MMS, 
          file = "outputs/PSA/ICER_PSA_MMS.csv",
          row.names = FALSE)
write.csv(df_incremental_PSA_MMS, 
          file = "outputs/PSA/incremental_PSA_MMS.csv",
          row.names = FALSE)

#############################
#### Trial Specification ####
#############################
for(i in 1:n_sim){ # i <- 1
  # Update parameter set for each scenario with next set of PSA drawn parameters
  l_psa_input_MET_TS <- update_param_list(l_params_all = l_params_MET_TS, params_updated = df_psa_params_MET_TS[i, ])
  l_psa_input_BUP_TS <- update_param_list(l_params_all = l_params_BUP_TS, params_updated = df_psa_params_BUP_TS[i, ])
  
  # Run model and generate outputs
  l_outcomes_MET_TS <- outcomes(l_params_all = l_psa_input_MET_TS, v_params_calib = v_calib_post_map, PSA = TRUE)
  l_outcomes_BUP_TS <- outcomes(l_params_all = l_psa_input_BUP_TS, v_params_calib = v_calib_post_map, PSA = TRUE)
  
  # Extract cost and QALY outputs
  df_outcomes_MET_PSA_TS <- rbind(df_outcomes_MET_PSA_TS, l_outcomes_MET_TS$df_outcomes)
  df_outcomes_BUP_PSA_TS <- rbind(df_outcomes_BUP_PSA_TS, l_outcomes_BUP_TS$df_outcomes)
  
  # Calculate ICER (societal and health sector perspective)
  l_ICER_TS <- ICER(outcomes_comp = l_outcomes_MET_TS, outcomes_int = l_outcomes_BUP_TS)
  
  df_incremental_PSA_TS <- rbind(df_incremental_PSA_TS, l_ICER_TS$df_incremental)
  
  df_ICER_PSA_TS <- rbind(df_ICER_PSA_TS, l_ICER_TS$df_icer)
}

### Output results
## As .RData
save(df_outcomes_MET_PSA_TS, 
     file = "outputs/PSA/Trial Specification/outcomes_MET_PSA_TS.RData")
save(df_outcomes_BUP_PSA_TS, 
     file = "outputs/PSA/Trial Specification/outcomes_BUP_PSA_TS.RData")
save(df_ICER_PSA_TS, 
     file = "outputs/PSA/Trial Specification/ICER_PSA_TS.RData")
save(df_incremental_PSA_TS,
     file = "outputs/PSA/Trial Specification/incremental_PSA_MMS.RData")
## As .csv
write.csv(df_outcomes_MET_PSA_TS, 
          file = "outputs/PSA/Trial Specification/outcomes_MET_PSA_TS.csv",
          row.names = FALSE)
write.csv(df_outcomes_BUP_PSA_TS, 
          file = "outputs/PSA/Trial Specification/outcomes_BUP_PSA_TS.csv",
          row.names = FALSE)
write.csv(df_ICER_PSA_TS, 
          file = "outputs/PSA/Trial Specification/ICER_PSA_TS.csv",
          row.names = FALSE)
write.csv(df_incremental_PSA_TS, 
          file = "outputs/PSA/Trial Specification/incremental_PSA_TS.csv",
          row.names = FALSE)

################################
#### Original Specification ####
################################
for(i in 1:n_sim){ # i <- 1
  # Update parameter set for each scenario with next set of PSA drawn parameters
  l_psa_input_MET_OS <- update_param_list(l_params_all = l_params_MET_OS, params_updated = df_psa_params_OS[i, ])
  l_psa_input_BUP_OS <- update_param_list(l_params_all = l_params_BUP_OS, params_updated = df_psa_params_OS[i, ])
  
  # Run model and generate outputs
  l_outcomes_MET_OS <- outcomes(l_params_all = l_psa_input_MET_OS, v_params_calib = v_calib_post_map, PSA = TRUE)
  l_outcomes_BUP_OS <- outcomes(l_params_all = l_psa_input_BUP_OS, v_params_calib = v_calib_post_map, PSA = TRUE)
  
  # Extract cost and QALY outputs
  df_outcomes_MET_PSA_OS <- rbind(df_outcomes_MET_PSA_OS, l_outcomes_MET_OS$df_outcomes)
  df_outcomes_BUP_PSA_OS <- rbind(df_outcomes_BUP_PSA_OS, l_outcomes_BUP_OS$df_outcomes)
  
  # Calculate ICER (societal and health sector perspective)
  l_ICER_OS <- ICER(outcomes_comp = l_outcomes_MET_OS, outcomes_int = l_outcomes_BUP_OS)
  
  df_ICER_PSA_OS <- rbind(df_ICER_PSA_OS, l_ICER_OS$df_icer)
}

### Output results ###
## As .RData ##
save(df_outcomes_MET_PSA_OS, 
     file = "outpuOS/PSA/Modified Model Specification/outcomes_MET_PSA_OS.RData")
save(df_outcomes_BUP_PSA_OS, 
     file = "outpuOS/PSA/Modified Model Specification/outcomes_BUP_PSA_OS.RData")
save(df_ICER_PSA_OS, 
     file = "outpuOS/PSA/Modified Model Specification/ICER_PSA_OS.RData")
## As .csv
write.csv(df_outcomes_MET_PSA_OS, 
          file = "outpuOS/PSA/Modified Model Specification/outcomes_MET_PSA_OS.csv",
          row.names = FALSE)
write.csv(df_outcomes_BUP_PSA_OS, 
          file = "outpuOS/PSA/Modified Model Specification/outcomes_BUP_PSA_OS.csv",
          row.names = FALSE)
write.csv(df_ICER_PSA_OS, 
          file = "outpuOS/PSA/ICER_PSA_OS.csv",
          row.names = FALSE)

### Process PSA results
## Read-in saved results
## Modified Model Specification
load(file = "outputs/PSA/outcomes_MET_PSA_MMS.RData")
load(file = "outputs/PSA/outcomes_BUP_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/incremental_PSA_MMS.RData")
load(file = "outputs/PSA/ICER_PSA_MMS.RData")

## Trial Specification
#load(file = "outputs/PSA/outcomes_MET_PSA_TS.RData")
#load(file = "outputs/PSA/outcomes_BUP_PSA_TS.RData")
#load(file = "outputs/PSA/Modified Model Specification/incremental_PSA_TS.RData")
#load(file = "outputs/PSA/ICER_PSA_TS.RData")

## Original Specification
#load(file = "outputs/PSA/outcomes_MET_PSA_OS.RData")
#load(file = "outputs/PSA/outcomes_BUP_PSA_OS.RData")
#load(file = "outputs/PSA/Modified Model Specification/incremental_PSA_OS.RData")
#load(file = "outputs/PSA/ICER_PSA_OS.RData")

### Summary stats ###
## Modified Model Specification ##
# Methadone
tbl_df_summary_MET_MMS <- df_outcomes_MET_PSA_MMS %>% as.tibble() %>% gather("variable", "value") %>% group_by(variable) %>% summarize(mean = mean(value),
                                                                                                                                       sd = sd(value),
                                                                                                                                       q50 = quantile(value, probs = .5),
                                                                                                                                       q025 = quantile(value, probs = .025),
                                                                                                                                       q975 = quantile(value, probs = .975),
                                                                                                                                       min = min(value),
                                                                                                                                       max = max(value))
# BNX
tbl_df_summary_BUP_MMS <- df_outcomes_BUP_PSA_MMS %>% as.tibble() %>% gather("variable", "value") %>% group_by(variable) %>% summarize(mean = mean(value),
                                                                                                                                       sd = sd(value),
                                                                                                                                       q50 = quantile(value, probs = .5),
                                                                                                                                       q025 = quantile(value, probs = .025),
                                                                                                                                       q975 = quantile(value, probs = .975),
                                                                                                                                       min = min(value),
                                                                                                                                       max = max(value))
# ICER
tbl_df_summary_ICER_MMS <- df_ICER_PSA_MMS %>% as.tibble() %>% gather("variable", "value") %>% group_by(variable) %>% summarize(mean = mean(value),
                                                                                                                                sd = sd(value),
                                                                                                                                q50 = quantile(value, probs = .5),
                                                                                                                                q025 = quantile(value, probs = .025),
                                                                                                                                q975 = quantile(value, probs = .975),
                                                                                                                                min = min(value),
                                                                                                                                max = max(value)) # want to add countif for number of CS


df_PSA_summary <- as.data.frame()


### Produce scatter plot for ICERs
## Modified Model Specification ##
# 1-year
# Total
inc_qalys_1yr <- df_incremental_PSA_MMS[, "n_inc_qalys_TOTAL_1yr"]
inc_costs_1yr <- df_incremental_PSA_MMS[, "n_inc_costs_TOTAL_1yr"]

plot_PSA_MMS_1yr_scatter <- ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_1yr, y = inc_costs_1yr)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_MMS_1yr_scatter, 
       filename = "Plots/PSA/PSA-MMS-1yr.png", 
       width = 7, height = 10)

# Health Sector
inc_costs_1yr_health_sector <- df_incremental_PSA_MMS[, "n_inc_costs_HEALTH_SECTOR_1yr"]

plot_PSA_MMS_1yr_scatter_health_sector <- ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_1yr, y = inc_costs_1yr_health_sector)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_MMS_1yr_scatter_health_sector, 
       filename = "Plots/PSA/PSA-MMS-1yr-Health-Sector.png", 
       width = 7, height = 10)

# 5-year
inc_qalys_5yr <- df_incremental_PSA_MMS[, "n_inc_qalys_TOTAL_5yr"]
inc_costs_5yr <- df_incremental_PSA_MMS[, "n_inc_costs_TOTAL_5yr"]

plot_PSA_MMS_5yr_scatter <- ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_5yr, y = inc_costs_5yr)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_MMS_5yr_scatter, 
       filename = "Plots/PSA/PSA-MMS-5yr.png", 
       width = 7, height = 10)

# Health Sector
inc_costs_5yr_health_sector <- df_incremental_PSA_MMS[, "n_inc_costs_HEALTH_SECTOR_5yr"]

plot_PSA_MMS_5yr_scatter_health_sector <- ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_5yr, y = inc_costs_5yr_health_sector)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_MMS_5yr_scatter_health_sector, 
       filename = "Plots/PSA/PSA-MMS-5yr-Health-Sector.png", 
       width = 7, height = 10)

# 10-year
#inc_qalys_10yr <- df_incremental_PSA_MMS[, "Incremental QALYs (10-year)"]
#inc_costs_10yr <- df_incremental_PSA_MMS[, "Incremental Costs (10-year)"]

#ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_10yr, y = inc_costs_10yr)) +
#  geom_point() +
#  geom_hline(yintercept = 0) +
#  geom_vline(xintercept = 0)
#xlim(min(a), max(a)) +
#ylim(min(b), max(b))

# Lifetime
#inc_qalys_lifetime <- df_incremental_PSA_MMS[, "Incremental QALYs (10-year)"]
#inc_costs_lifetime <- df_incremental_PSA_MMS[, "Incremental Costs (10-year)"]

#ggplot(df_incremental_PSA_MMS, aes(x = inc_qalys_lifetime, y = inc_costs_lifetime)) +
#  geom_point() +
#  geom_hline(yintercept = 0) +
#  geom_vline(xintercept = 0)
#xlim(min(a), max(a)) +
#ylim(min(b), max(b))

## Trial Specification ##
# 1-year
# Total
inc_qalys_1yr <- df_incremental_PSA_TS[, "n_inc_qalys_TOTAL_1yr"]
inc_costs_1yr <- df_incremental_PSA_TS[, "n_inc_costs_TOTAL_1yr"]

plot_PSA_TS_1yr_scatter <- ggplot(df_incremental_PSA_TS, aes(x = inc_qalys_1yr, y = inc_costs_1yr)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_TS_1yr_scatter, 
       filename = "Plots/PSA/PSA-TS-1yr.png", 
       width = 7, height = 10)

# Health Sector
inc_costs_1yr_health_sector <- df_incremental_PSA_TS[, "n_inc_costs_HEALTH_SECTOR_1yr"]

plot_PSA_TS_1yr_scatter_health_sector <- ggplot(df_incremental_PSA_TS, aes(x = inc_qalys_1yr, y = inc_costs_1yr_health_sector)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_TS_1yr_scatter_health_sector, 
       filename = "Plots/PSA/PSA-TS-1yr-Health-Sector.png", 
       width = 7, height = 10)

# 5-year
# Total
inc_qalys_5yr <- df_incremental_PSA_TS[, "n_inc_qalys_TOTAL_5yr"]
inc_costs_5yr <- df_incremental_PSA_TS[, "n_inc_costs_TOTAL_5yr"]

plot_PSA_TS_5yr_scatter <- ggplot(df_incremental_PSA_TS, aes(x = inc_qalys_5yr, y = inc_costs_5yr)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_TS_5yr_scatter, 
       filename = "Plots/PSA/PSA-TS-5yr.png", 
       width = 7, height = 10)

# Health Sector
inc_costs_5yr_health_sector <- df_incremental_PSA_TS[, "n_inc_costs_HEALTH_SECTOR_5yr"]

plot_PSA_TS_5yr_scatter_health_sector <- ggplot(df_incremental_PSA_TS, aes(x = inc_qalys_1yr, y = inc_costs_5yr_health_sector)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 100000, intercept = 0)

ggsave(plot_PSA_TS_5yr_scatter_health_sector, 
       filename = "Plots/PSA/PSA-TS-5yr-Health-Sector.png", 
       width = 7, height = 10)

### Produce CEAC plot from cost-effectiveness results
## Prep data
# Set number of points (one for each quantile)
quantiles <- seq(0,1,0.01)

tbl_df_CEAC <- df_ICER_PSA %>% as.tibble() %>% gather("variable", "value") %>% group_by(variable) %>% summarise_at(vars('value'), 
                                                                                                                   funs(quantile = list(as.tibble(as.list(quantile(., probs = quantiles)))))) %>% unnest

## Create plot
#CEAC <- function(scenario = scenario){
# Subset by ICER (e.g. societal-lifetime)
tbl_df_CEAC_subset <- tbl_df_CEAC %>% filter(variable == "ICER (1-year)") %>% select(-variable)
tbl_df_CEAC_subset <- gather(tbl_df_CEAC_subset)

len <- nrow(tbl_df_CEAC_subset)
increment <- 1/(len - 1)
p <- seq(0, 1, increment)

# Final data for plotting
tbl_df_CEAC_plot <- cbind(tbl_df_CEAC_subset, p)

# Create plot
ggplot(data = tbl_df_CEAC_plot, aes(x = value, y = p)) +
  geom_line() + xlim(0, NA)
#}
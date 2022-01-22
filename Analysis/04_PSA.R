rm(list = ls()) # to clean the workspace

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
source("Analysis/00_load_parameters.R")

# Set population size for dirichlet draws
n_pop_cohort <- 29000
n_pop_trial  <- 272
n_sim <- 3 # just to test function (will be set as n_sim)

### PSA model outputs
### Run Markov model for PSA draws and return outputs ###
# Generate PSA parameter draws

######################################
#### Modified Model Specification ####
######################################
# BNX scenario
df_psa_params_BUP_MMS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_cohort,
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

# MET scenario
df_psa_params_MET_MMS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_cohort,
                                             file.death_hr = "data/death_hr.csv",
                                             file.frailty = "data/frailty.csv",
                                             file.weibull = "data/Modified Model Specification/weibull.csv",
                                             file.unconditional = "data/Modified Model Specification/unconditional.csv",
                                             file.overdose = "data/overdose.csv",
                                             file.hiv = "data/hiv_sero.csv",
                                             file.hcv = "data/hcv_sero.csv",
                                             file.costs = "data/Modified Model Specification/costs.csv",
                                             file.crime_costs = "data/Modified Model Specification/crime_costs.csv",
                                             file.qalys = "data/Modified Model Specification/qalys.csv",
                                             file.imis_output = "outputs/Calibration/imis_output.RData")

#############################
#### Trial Specification ####
#############################
# BNX scenario
df_psa_params_BUP_TS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_trial,
                                            file.death_hr = "data/death_hr.csv",
                                            file.frailty = "data/frailty.csv",
                                            file.weibull = "data/Trial Specification/weibull.csv",
                                            file.unconditional = "data/Trial Specification/unconditional_bup.csv",
                                            file.overdose = "data/overdose.csv",
                                            file.hiv = "data/hiv_sero.csv",
                                            file.hcv = "data/hcv_sero.csv",
                                            file.costs = "data/Trial Specification/costs.csv",
                                            file.crime_costs = "data/Trial Specification/crime_costs.csv",
                                            file.qalys = "data/Trial Specification/qalys.csv",
                                            file.imis_output = "outputs/imis_output.RData")

# MET scenario
df_psa_params_MET_TS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_trial,
                                            file.death_hr = "data/death_hr.csv",
                                            file.frailty = "data/frailty.csv",
                                            file.weibull = "data/Trial Specification/weibull.csv",
                                            file.unconditional = "data/Trial Specification/unconditional_met.csv",
                                            file.overdose = "data/overdose.csv",
                                            file.hiv = "data/hiv_sero.csv",
                                            file.hcv = "data/hcv_sero.csv",
                                            file.costs = "data/Trial Specification/costs.csv",
                                            file.crime_costs = "data/Trial Specification/crime_costs.csv",
                                            file.qalys = "data/Trial Specification/qalys.csv",
                                            file.imis_output = "outputs/imis_output.RData")

################################
#### Original Specification ####
################################
# BNX scenario
df_psa_params_BUP_OS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_cohort,
                                            file.death_hr = "data/death_hr.csv",
                                            file.frailty = "data/frailty.csv",
                                            file.weibull = "data/Original Specification/weibull.csv",
                                            file.unconditional = "data/Original Specification/unconditional.csv",
                                            file.overdose = "data/overdose.csv",
                                            file.hiv = "data/hiv_sero.csv",
                                            file.hcv = "data/hcv_sero.csv",
                                            file.costs = "data/Original Specification/costs.csv",
                                            file.crime_costs = "data/Original Specification/crime_costs.csv",
                                            file.qalys = "data/Original Specification/qalys.csv",
                                            file.imis_output = "outputs/imis_output.RData")

# MET scenario
df_psa_params_MET_OS <- generate_psa_params(n_sim = n_sim, seed = 3730687, n_pop = n_pop_cohort,
                                            file.death_hr = "data/death_hr.csv",
                                            file.frailty = "data/frailty.csv",
                                            file.weibull = "data/Original Specification/weibull.csv",
                                            file.unconditional = "data/Original Specification/unconditional.csv",
                                            file.overdose = "data/overdose.csv",
                                            file.hiv = "data/hiv_sero.csv",
                                            file.hcv = "data/hcv_sero.csv",
                                            file.costs = "data/Original Specification/costs.csv",
                                            file.crime_costs = "data/Original Specification/crime_costs.csv",
                                            file.qalys = "data/Original Specification/qalys.csv",
                                            file.imis_output = "outputs/imis_output.RData")

### Run decision model on each parameter set of PSA input dataset to produce
### PSA outputs for cost and effects
#n_sim <- 500 # just to test function (will be set as n_sim)

# Initialize data frames
# Modified Model Specification
df_outcomes_MET_PSA_MMS <- data.frame()
df_outcomes_BUP_PSA_MMS <- data.frame()
df_ICER_PSA_MMS <- data.frame()

# Trial Specification
df_outcomes_MET_PSA_TS <- data.frame()
df_outcomes_BUP_PSA_TS <- data.frame()
df_ICER_PSA_TS <- data.frame()

# Original Specification
df_outcomes_MET_PSA_OS <- data.frame()
df_outcomes_BUP_PSA_OS <- data.frame()
df_ICER_PSA_OS <- data.frame()

v_outcomes_names <- c("Total Costs (1-year)", "Total Costs (5-year)", "Total Costs (10-year)", "Total Costs (Lifetime)", "Health Sector Costs (1-year)", "Health Sector Costs (5-year)", "Health Sector Costs (10-year)", "Health Sector Costs (Lifetime)",
                      "Criminal Costs (1-year)", "Criminal Costs (5-year)", "Criminal Costs (10-year)", "Criminal Costs (Lifetime)", "Treatment Costs (1-year)", "Treatment Costs (5-year)", "Treatment Costs (10-year)", "Treatment Costs (Lifetime)",
                      "Total QALYs (1-year)", "Total QALYs (5-year)", "Total QALYs (10-year)", "Total QALYs (Lifetime)")
v_ICER_names <- c("ICER (1-year)", "ICER (5-year)", "ICER (10-year)", "ICER (Lifetime)",
                  "ICER (Health Sector 1-year)", "ICER (Health Sector 5-year)", "ICER (Health Sector 10-year)", "ICER (Health Sector Lifetime)",
                  "Incremental Costs (1-year)", "Incremental QALYs (1-year)", "Incremental Costs (5-year)", "Incremental QALYs (5-year)", "Incremental Costs (10-year)", "Incremental QALYs (10-year)", "Incremental Costs (Lifetime)", "Incremental QALYs (Lifetime)",
                  "Incremental Costs (Health Sector 1-year)", "Incremental QALYs (Health Sector 1-year)", "Incremental Costs (Health Sector 5-year)", "Incremental QALYs (Health Sector 5-year)", "Incremental Costs (Health Sector 10-year)", "Incremental QALYs (Health Sector 10-year)", "Incremental Costs (Health Sector Lifetime)", "Incremental QALYs (Health Sector Lifetime)")

# Modified Model Specification
for(i in 1:n_sim){ # i <- 1
  # Update parameter set for each scenario with next set of PSA drawn parameters
  l_psa_input_MET_MMS <- update_param_list(l_params_all = l_params_MET_MMS, params_updated = df_psa_params_MET_MMS[i, ])
  l_psa_input_BUP_MMS <- update_param_list(l_params_all = l_params_BUP_MMS, params_updated = df_psa_params_BUP_MMS[i, ])
  
  # Run model and generate outputs
  l_outcomes_MET_MMS <- outcomes(l_params_all = l_psa_input_MET_MMS, v_params_calib = v_calib_post_map)
  l_outcomes_BUP_MMS <- outcomes(l_params_all = l_psa_input_BUP_MMS, v_params_calib = v_calib_post_map)
  
  # Extract cost and QALY outputs
  df_outcomes_MET_PSA_MMS <- rbind(df_outcomes_MET_PSA_MMS, l_outcomes_MET_MMS$df_outcomes)
  df_outcomes_BUP_PSA_MMS <- rbind(df_outcomes_BUP_PSA_MMS, l_outcomes_BUP_MMS$df_outcomes)

  # Calculate ICER (societal and health sector perspective)
  l_ICER_MMS <- ICER(outcomes_comp = l_outcomes_MET_MMS, outcomes_int = l_outcomes_BUP_MMS)
  
  df_ICER_PSA_MMS <- rbind(df_ICER_PSA_MMS, l_ICER$df_icer)
}

### Output results
## As .RData
save(df_outcomes_MET_PSA, 
     file = "outputs/PSA/outcomes_MET_PSA.RData")
save(df_outcomes_BUP_PSA, 
     file = "outputs/PSA/outcomes_BUP_PSA.RData")
save(df_ICER_PSA, 
     file = "outputs/PSA/ICER_PSA.RData")
## As .csv
write.csv(df_outcomes_MET_PSA, 
          file = "outputs/PSA/outcomes_MET_PSA.csv",
          row.names = FALSE)
write.csv(df_outcomes_BUP_PSA, 
          file = "outputs/PSA/outcomes_BUP_PSA.csv",
          row.names = FALSE)
write.csv(df_ICER_PSA, 
          file = "outputs/PSA/ICER_PSA.csv",
          row.names = FALSE)

### Process PSA results
## Read-in saved results
load(file = "outputs/PSA/outcomes_MET_PSA.RData")
load(file = "outputs/PSA/outcomes_BUP_PSA.RData")
load(file = "outputs/PSA/ICER_PSA.RData")

## Summary stats
# Methadone
tbl_df_summary_MET <- df_outcomes_MET_PSA %>% as.tibble() %>% gather("variable", "value") %>% group_by(variable) %>% summarize(mean = mean(value),
                                                                                                                               sd = sd(value),
                                                                                                                               q50 = quantile(value, probs = .5),
                                                                                                                               q025 = quantile(value, probs = .025),
                                                                                                                               q975 = quantile(value, probs = .975),
                                                                                                                               min = min(value),
                                                                                                                               max = max(value))
# BNX
tbl_df_summary_BUP <- df_outcomes_BUP_PSA %>% as.tibble() %>% gather("variable", "value") %>% group_by(variable) %>% summarize(mean = mean(value),
                                                                                                                               sd = sd(value),
                                                                                                                               q50 = quantile(value, probs = .5),
                                                                                                                               q025 = quantile(value, probs = .025),
                                                                                                                               q975 = quantile(value, probs = .975),
                                                                                                                               min = min(value),
                                                                                                                               max = max(value))
# ICER
tbl_df_summary_ICER <- df_ICER_PSA %>% as.tibble() %>% gather("variable", "value") %>% group_by(variable) %>% summarize(mean = mean(value),
                                                                                                                        sd = sd(value),
                                                                                                                        q50 = quantile(value, probs = .5),
                                                                                                                        q025 = quantile(value, probs = .025),
                                                                                                                        q975 = quantile(value, probs = .975),
                                                                                                                        min = min(value),
                                                                                                                        max = max(value)) # want to add countif for number of CS


df_PSA_summary <- as.data.frame()


### Produce scatter plot for ICERs


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
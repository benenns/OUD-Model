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

### Run decision model on each parameter set of PSA input dataset to produce
### PSA outputs for cost and effects
n_sim <- 250 # just to test function (will be set as n_sim)

# Initialize data frames
df_outcomes_MET_PSA <- data.frame()
df_outcomes_BUP_PSA <- data.frame()
df_ICER_PSA <- data.frame()

v_outcomes_names <- c("Total Costs (1-year)", "Total Costs (5-year)", "Total Costs (10-year)", "Total Costs (Lifetime)", "Health Sector Costs (1-year)", "Health Sector Costs (5-year)", "Health Sector Costs (10-year)", "Health Sector Costs (Lifetime)",
                      "Criminal Costs (1-year)", "Criminal Costs (5-year)", "Criminal Costs (10-year)", "Criminal Costs (Lifetime)", "Treatment Costs (1-year)", "Treatment Costs (5-year)", "Treatment Costs (10-year)", "Treatment Costs (Lifetime)",
                      "Total QALYs (1-year)", "Total QALYs (5-year)", "Total QALYs (10-year)", "Total QALYs (Lifetime)")
v_ICER_names <- c("ICER (1-year)", "ICER (5-year)", "ICER (10-year)", "ICER (Lifetime)", "ICER (Health Sector 1-year)", "ICER (Health Sector 5-year)", "ICER (Health Sector 10-year)", "ICER (Health Sector Lifetime)")

for(i in 1:n_sim){ # i <- 1
  # Update parameter set for each scenario with next set of PSA drawn parameters
  l_psa_input_MET <- update_param_list(l_params_all = l_params_MET, params_updated = df_psa_params[i,])
  l_psa_input_BUP <- update_param_list(l_params_all = l_params_BUP, params_updated = df_psa_params[i,])
  
  # Run model and generate outputs
  l_outcomes_MET <- outcomes(l_params_all = l_psa_input_MET, v_params_calib = v_calib_post_map)
  l_outcomes_BUP <- outcomes(l_params_all = l_psa_input_BUP, v_params_calib = v_calib_post_map)
  
  # Extract cost and QALY outputs
  df_outcomes_MET_PSA <- rbind(df_outcomes_MET_PSA, l_outcomes_MET$df_outcomes)
  df_outcomes_BUP_PSA <- rbind(df_outcomes_BUP_PSA, l_outcomes_BUP$df_outcomes)

  # Calculate ICER (societal and health sector perspective)
  l_ICER <- ICER(outcomes_comp = l_outcomes_MET, outcomes_int = l_outcomes_BUP)
  
  df_ICER_PSA <- rbind(df_ICER_PSA, l_ICER$df_icer)
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
                                                                                                                               max = max(value),
                                                                                                                               percent_cs = length(value[variable<=0])) # want to add countif for number of CS
# ICER
tbl_df_summary_ICER <- df_ICER_PSA %>% as.tibble() %>% gather("variable", "value") %>% group_by(variable) %>% summarize(mean = mean(value),
                                                                                                                        sd = sd(value),
                                                                                                                        q50 = quantile(value, probs = .5),
                                                                                                                        q025 = quantile(value, probs = .025),
                                                                                                                        q975 = quantile(value, probs = .975),
                                                                                                                        min = min(value),
                                                                                                                        max = max(value))


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
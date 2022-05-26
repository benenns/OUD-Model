rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(data.table)
library(formattable)
library(tidyr)
library(RColorBrewer)

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/ICER_functions.R")

# Load parameters
source("Analysis/00_load_parameters.R")

# Load DSA parameters
#################
### HRU costs ###
#################
# DSA data
df_dsa_HRU_costs_MMS <- read.csv(file = "data/DSA/Modified Model Specification/HRU_costs.csv", row.names = 1, header = TRUE)

# MMS
v_dsa_HRU_costs_low_MMS <- unlist(df_dsa_HRU_costs_MMS["pe_low",])
v_dsa_HRU_costs_high_MMS <- unlist(df_dsa_HRU_costs_MMS["pe_high",])

###################
### Crime costs ###
###################
# DSA data
df_dsa_crime_costs_MMS <- read.csv(file = "data/DSA/Modified Model Specification/crime_costs.csv", row.names = 1, header = TRUE)

# MMS
v_dsa_crime_costs_low_MMS <- unlist(df_dsa_crime_costs_MMS["pe_low",])
v_dsa_crime_costs_high_MMS <- unlist(df_dsa_crime_costs_MMS["pe_high",])
v_dsa_crime_costs_alt_MMS <- unlist(df_dsa_crime_costs_MMS["pe_alt",])
#v_dsa_crime_costs_reduced_MMS <- unlist(df_dsa_crime_costs_MMS["pe_reduced_trial",])

############################################
#### Deterministic sensitivity analysis ####
############################################

################
### Baseline ###
################
# MMS
#l_outcomes_MET_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map)
#l_outcomes_BUP_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map)
#ICER_MMS <- ICER(outcomes_comp = l_outcomes_MET_MMS, outcomes_int = l_outcomes_BUP_MMS)
#################
### HRU Costs ###
#################
# Low
# MMS
l_outcomes_MET_HRU_costs_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_HRU_costs_low_MMS)
l_outcomes_BUP_HRU_costs_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_HRU_costs_low_MMS)
ICER_HRU_costs_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_HRU_costs_low_MMS, outcomes_int = l_outcomes_BUP_HRU_costs_low_MMS)

# High
# MMS
l_outcomes_MET_HRU_costs_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_HRU_costs_high_MMS)
l_outcomes_BUP_HRU_costs_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_HRU_costs_high_MMS)
ICER_HRU_costs_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_HRU_costs_high_MMS, outcomes_int = l_outcomes_BUP_HRU_costs_high_MMS)

###################
### Crime Costs ###
###################
# Low
# MMS
l_outcomes_MET_crime_costs_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_low_MMS)
l_outcomes_BUP_crime_costs_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_low_MMS)
ICER_crime_costs_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_low_MMS, outcomes_int = l_outcomes_BUP_crime_costs_low_MMS)

# High
# MMS
l_outcomes_MET_crime_costs_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_high_MMS)
l_outcomes_BUP_crime_costs_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_high_MMS)
ICER_crime_costs_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_high_MMS, outcomes_int = l_outcomes_BUP_crime_costs_high_MMS)

# Alternative (Krebs 2014 estimates)
# MMS
l_outcomes_MET_crime_costs_alt_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_MMS)
l_outcomes_BUP_crime_costs_alt_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_MMS)
ICER_crime_costs_alt_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_alt_MMS, outcomes_int = l_outcomes_BUP_crime_costs_alt_MMS)

#############
### Costs ###
#############
# Baseline
# HRU costs
v_HRU_costs_MMS <- c(ICER_HRU_costs_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_HRU_costs_high_MMS$df_incremental$n_inc_costs_TOTAL_life)

# Crime costs
v_crime_costs_MMS <- c(ICER_crime_costs_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_crime_costs_high_MMS$df_incremental$n_inc_costs_TOTAL_life)
#v_crime_costs_high_MMS <- c(ICER_crime_costs_high_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_crime_costs_high_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_crime_costs_alt_MMS <- c(ICER_crime_costs_alt_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_crime_costs_alt_MMS$df_incremental$n_inc_costs_TOTAL_life)
#v_crime_costs_reduced_MMS <- c(ICER_crime_costs_reduced_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_crime_costs_reduced_MMS$df_incremental$n_inc_costs_TOTAL_life)

m_costs_MMS <- rbind(v_HRU_costs_MMS, v_crime_costs_MMS, v_crime_costs_alt_MMS)

df_costs_MMS <- as.data.frame(m_costs_MMS)
colnames(df_costs_MMS) <- c("Lower", "Upper")
df_costs_MMS <- as_data_frame(df_costs_MMS) %>% #mutate(base = ICER_MMS$df_incremental$n_inc_costs_TOTAL_life,
                                                #       diff = ifelse(abs(Upper - Lower) > 0, abs(Upper - Lower), abs(base - Upper))) %>%
  add_column(var_name = c("HRU Costs (Range)", "Crime Costs (Range)", "Crime Costs (Alternative Estimates)"))

# Save outputs
## As .RData ##
save(df_costs_MMS, 
     file = "outputs/DSA/Modified Model Specification/df_costs_MMS.RData")

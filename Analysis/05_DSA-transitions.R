rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(data.table)
library(formattable)
library(tidyr)
library(RColorBrewer)
library(Rmisc)
library(grid)
library(gridExtra)
library(lattice)

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/ICER_functions.R")

# Load parameters
source("Analysis/00_load_parameters.R")

# Load DSA parameters
###################
### Transitions ###
###################
df_dsa_frailty <- read.csv(file = "data/DSA/frailty.csv", row.names = 1, header = TRUE)
#df_dsa_frailty_TS <- read.csv(file = "data/DSA/Trial Specification/frailty.csv", row.names = 1, header = TRUE)

v_dsa_frailty_episode <- unlist(df_dsa_frailty["pe_episode_frailty",])
v_dsa_frailty_concurrent <- unlist(df_dsa_frailty["pe_concurrent_frailty",])
v_dsa_frailty_inj <- unlist(df_dsa_frailty["pe_inj_frailty",])
# v_dsa_frailty_kurz <- unlist(df_dsa_frailty["pe_kurz_frailty",])

############################################
#### Deterministic sensitivity analysis ####
############################################

################
### Baseline ###
################
# MMS
# l_outcomes_MET_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map)
# l_outcomes_BUP_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map)
# ICER_MMS <- ICER(outcomes_comp = l_outcomes_MET_MMS, outcomes_int = l_outcomes_BUP_MMS)

###################
### Transitions ###
###################
## Frailty ##
# Episode
l_outcomes_MET_frailty_episode_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_episode)
l_outcomes_BUP_frailty_episode_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_episode)
ICER_frailty_episode_MMS <- ICER(outcomes_comp = l_outcomes_MET_frailty_episode_MMS, outcomes_int = l_outcomes_BUP_frailty_episode_MMS)

# Concurrent opioid use
l_outcomes_MET_frailty_concurrent_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_concurrent)
l_outcomes_BUP_frailty_concurrent_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_concurrent)
ICER_frailty_concurrent_MMS <- ICER(outcomes_comp = l_outcomes_MET_frailty_concurrent_MMS, outcomes_int = l_outcomes_BUP_frailty_concurrent_MMS)

# Injection multiplier
l_outcomes_MET_frailty_inj_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_inj)
l_outcomes_BUP_frailty_inj_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_inj)
ICER_frailty_inj_MMS <- ICER(outcomes_comp = l_outcomes_MET_frailty_inj_MMS, outcomes_int = l_outcomes_BUP_frailty_inj_MMS)

# Kurz results (combined for treatment states)
# l_outcomes_MET_frailty_kurz_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_kurz)
# l_outcomes_BUP_frailty_kurz_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_kurz)
# ICER_frailty_kurz_MMS <- ICER(outcomes_comp = l_outcomes_MET_frailty_kurz_MMS, outcomes_int = l_outcomes_BUP_frailty_kurz_MMS)

################
### Baseline ###
################
#df_baseline_MMS <- data.frame(ICER_MMS$df_incremental, ICER_MMS$df_icer)

###################
### Transitions ###
###################
# Baseline
# df_transitions_baseline_MMS <- data.frame(ICER_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_MMS$df_incremental$n_inc_qalys_TOTAL_life, 
#                                           ICER_MMS$df_icer$n_icer_TOTAL_life)

# Costs
v_transitions_frailty_episode_MMS <- c(ICER_frailty_episode_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_frailty_episode_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_transitions_frailty_concurrent_MMS <- c(ICER_frailty_concurrent_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_frailty_concurrent_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_transitions_frailty_inj_MMS <- c(ICER_frailty_inj_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_frailty_inj_MMS$df_incremental$n_inc_costs_TOTAL_life)
#v_transitions_frailty_kurz_MMS <- c(ICER_frailty_kurz_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_frailty_kurz_MMS$df_incremental$n_inc_costs_TOTAL_life)

m_transitions_costs_MMS <- rbind(v_transitions_frailty_episode_MMS, v_transitions_frailty_concurrent_MMS, v_transitions_frailty_inj_MMS)

df_transitions_costs_MMS <- as.data.frame(m_transitions_costs_MMS)
colnames(df_transitions_costs_MMS) <- c("Lower", "Upper")
df_transitions_costs_MMS <- as_data_frame(df_transitions_costs_MMS) %>% #mutate(diff = abs(Upper - Lower),
                                                                         #base = ICER_MMS$df_incremental$n_inc_costs_TOTAL_life) %>%
  add_column(var_name = c("No retention difference by episode", "No retention difference for concurrent use", "No retention difference for injection"))

# QALYs
v_transitions_frailty_episode_MMS <- c(ICER_frailty_episode_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_frailty_episode_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_transitions_frailty_concurrent_MMS <- c(ICER_frailty_concurrent_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_frailty_concurrent_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_transitions_frailty_inj_MMS <- c(ICER_frailty_inj_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_frailty_inj_MMS$df_incremental$n_inc_qalys_TOTAL_life)
#v_transitions_frailty_kurz_MMS <- c(ICER_frailty_kurz_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_frailty_kurz_MMS$df_incremental$n_inc_qalys_TOTAL_life)

m_transitions_qalys_MMS <- rbind(v_transitions_frailty_episode_MMS, v_transitions_frailty_concurrent_MMS, v_transitions_frailty_inj_MMS)

df_transitions_qalys_MMS <- as.data.frame(m_transitions_qalys_MMS)
colnames(df_transitions_qalys_MMS) <- c("Lower", "Upper")
df_transitions_qalys_MMS <- as_data_frame(df_transitions_qalys_MMS) %>% #mutate(diff = abs(Upper - Lower),
                                                                         #      base = ICER_MMS$df_incremental$n_inc_qalys_TOTAL_life) %>%
  add_column(var_name = c("No retention difference by episode", "No retention difference for concurrent use", "No retention difference for injection"))

# Save output
## As .RData ##
save(df_transitions_costs_MMS, 
     file = "outputs/DSA/Modified Model Specification/df_transitions_costs_MMS.RData")
save(df_transitions_qalys_MMS, 
     file = "outputs/DSA/Modified Model Specification/df_transitions_qalys_MMS.RData")

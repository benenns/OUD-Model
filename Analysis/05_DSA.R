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

################
### Overdose ###
################
# DSA data
df_dsa_overdose <- read.csv(file = "data/DSA/overdose.csv", row.names = 1, header = TRUE)
df_dsa_fentanyl <- read.csv(file = "data/DSA/fentanyl.csv", row.names = 1, header = TRUE)
# Load posterior
load(file = "outputs/Calibration/summary_posterior.RData")

## Fatal overdose ##
# Witnessed OD
v_dsa_witness_low <- unlist(df_dsa_overdose["low", "p_witness"])
v_dsa_witness_high <- unlist(df_dsa_overdose["high", "p_witness"])
names(v_dsa_witness_low) <- c("p_witness")
names(v_dsa_witness_high) <- c("p_witness")

# Naloxone prevalence
v_dsa_NX_prob_low <- unlist(df_dsa_overdose["low", "p_NX_used"])
v_dsa_NX_prob_high <- unlist(df_dsa_overdose["high", "p_NX_used"])
names(v_dsa_NX_prob_low) <- c("p_NX_used")
names(v_dsa_NX_prob_high) <- c("p_NX_used")

# Naloxone effectiveness
v_dsa_NX_success_low <- unlist(df_dsa_overdose["low", "p_NX_success"])
v_dsa_NX_success_high <- unlist(df_dsa_overdose["high", "p_NX_success"])
names(v_dsa_NX_success_low) <- c("p_NX_success")
names(v_dsa_NX_success_high) <- c("p_NX_success")

# Fatal overdose risk
v_dsa_fatal_OD_low <- unlist(df_posterior_summ["n_fatal_OD", "2.5%"])
v_dsa_fatal_OD_high <- unlist(df_posterior_summ["n_fatal_OD", "97.5%"])
names(v_dsa_fatal_OD_low) <- c("n_fatal_OD")
names(v_dsa_fatal_OD_high) <- c("n_fatal_OD")

## Non-fatal overdose ##
# Fentanyl prevalence
v_dsa_fent_exp_2020_low <- unlist(df_dsa_fentanyl["low", "BC"])
v_dsa_fent_exp_2020_high <- unlist(df_dsa_fentanyl["high", "BC"])
names(v_dsa_fent_exp_2020_low) <- c("p_fent_exp_2020")
names(v_dsa_fent_exp_2020_high) <- c("p_fent_exp_2020")

# Fent OD multiplier
v_dsa_fent_OD_mult_low <- unlist(df_posterior_summ["n_fent_OD_mult", "2.5%"])
v_dsa_fent_OD_mult_high <- unlist(df_posterior_summ["n_fent_OD_mult", "97.5%"])
names(v_dsa_fent_OD_mult_low) <- c("n_fent_OD_mult")
names(v_dsa_fent_OD_mult_high) <- c("n_fent_OD_mult")

# Reduction in fentanyl exposure for non-injection v. injection
v_dsa_ni_fent_reduction_low <- unlist(df_dsa_overdose["low", "p_ni_fent_reduction"])
v_dsa_ni_fent_reduction_high <- unlist(df_dsa_overdose["high", "p_ni_fent_reduction"])
names(v_dsa_ni_fent_reduction_low) <- c("p_ni_fent_reduction")
names(v_dsa_ni_fent_reduction_high) <- c("p_ni_fent_reduction")

# BUP OD multiplier
v_dsa_BUP_OD_mult_low <- unlist(df_dsa_overdose["low", "n_BUP_OD_mult"])
v_dsa_BUP_OD_mult_high <- unlist(df_dsa_overdose["high", "n_BUP_OD_mult"])

# MET OD multiplier
v_dsa_MET_OD_mult_low <- unlist(df_dsa_overdose["low", "n_MET_OD_mult"])
v_dsa_MET_OD_mult_high <- unlist(df_dsa_overdose["high", "n_MET_OD_mult"])

# REL OD multiplier
v_dsa_REL_OD_mult_low <- unlist(df_dsa_overdose["low", "n_REL_OD_mult"])
v_dsa_REL_OD_mult_high <- unlist(df_dsa_overdose["high", "n_REL_OD_mult"])

#################
### HRU costs ###
#################
# DSA data
#df_dsa_HRU_costs_MMS <- read.csv(file = "data/DSA/Modified Model Specification/HRU_costs.csv", row.names = 1, header = TRUE)
#df_dsa_HRU_costs_TS <- read.csv(file = "data/DSA/Trial Specification/HRU_costs.csv", row.names = 1, header = TRUE)
# MMS
#v_dsa_HRU_costs_alt_MMS <- unlist(df_dsa_HRU_costs_MMS["pe_alt",])
# TS
#v_dsa_HRU_costs_alt_TS <- unlist(df_dsa_HRU_costs_TS["pe_alt",])

###################
### Crime costs ###
###################
# DSA data
df_dsa_crime_costs_MMS <- read.csv(file = "data/DSA/Modified Model Specification/crime_costs.csv", row.names = 1, header = TRUE)
df_dsa_crime_costs_TS <- read.csv(file = "data/DSA/Trial Specification/crime_costs.csv", row.names = 1, header = TRUE)
# MMS
v_dsa_crime_costs_low_MMS <- unlist(df_dsa_crime_costs_MMS["pe_low",])
v_dsa_crime_costs_high_MMS <- unlist(df_dsa_crime_costs_MMS["pe_high",])
v_dsa_crime_costs_alt_MMS <- unlist(df_dsa_crime_costs_MMS["pe_alt",])
v_dsa_crime_costs_reduced_MMS <- unlist(df_dsa_crime_costs_MMS["pe_reduced",])
# TS
v_dsa_crime_costs_low_TS <- unlist(df_dsa_crime_costs_TS["pe_low",])
v_dsa_crime_costs_high_TS <- unlist(df_dsa_crime_costs_TS["pe_high",])
v_dsa_crime_costs_alt_TS <- unlist(df_dsa_crime_costs_TS["pe_alt",])
v_dsa_crime_costs_reduced_TS <- unlist(df_dsa_crime_costs_TS["pe_reduced",])

#############
### QALYs ###
#############
df_dsa_qalys_MMS <- read.csv(file = "data/DSA/Modified Model Specification/qalys.csv", row.names = 1, header = TRUE)
df_dsa_qalys_TS <- read.csv(file = "data/DSA/Trial Specification/qalys.csv", row.names = 1, header = TRUE)
# MMS
v_dsa_qalys_reduced_eq_5d_5l_MMS <- unlist(df_dsa_qalys_MMS["pe_reduced_eq_5d_5l",])
v_dsa_qalys_low_MMS <- unlist(df_dsa_qalys_MMS["pe_low",])
v_dsa_qalys_high_MMS <- unlist(df_dsa_qalys_MMS["pe_high",])
v_dsa_qalys_eq_5d_3l_MMS <- unlist(df_dsa_qalys_MMS["pe_eq_5d_3l",])
v_dsa_qalys_hui_3_MMS <- unlist(df_dsa_qalys_MMS["pe_hui_3",])
v_dsa_qalys_odn_low_MMS <- unlist(df_dsa_qalys_MMS["pe_odn_low",])
# TS
v_dsa_qalys_reduced_eq_5d_5l_TS <- unlist(df_dsa_qalys_TS["pe_reduced_eq_5d_5l",])
v_dsa_qalys_low_TS <- unlist(df_dsa_qalys_TS["pe_low",])
v_dsa_qalys_high_TS <- unlist(df_dsa_qalys_TS["pe_high",])
v_dsa_qalys_eq_5d_3l_TS <- unlist(df_dsa_qalys_TS["pe_eq_5d_3l",])
v_dsa_qalys_hui_3_TS <- unlist(df_dsa_qalys_TS["pe_hui_3",])
v_dsa_qalys_odn_low_TS <- unlist(df_dsa_qalys_TS["pe_odn_low",])

# Transitions

# Province-specific

############################################
#### Deterministic sensitivity analysis ####
############################################

################
### Baseline ###
################
# MMS
l_outcomes_MET_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map)
l_outcomes_BUP_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map)
ICER_MMS <- ICER(outcomes_comp = l_outcomes_MET_MMS, outcomes_int = l_outcomes_BUP_MMS)
# TS
l_outcomes_MET_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map)
l_outcomes_BUP_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map)
ICER_TS <- ICER(outcomes_comp = l_outcomes_MET_TS, outcomes_int = l_outcomes_BUP_TS)

###################
### Crime Costs ###
###################
# Low
# MMS
l_outcomes_MET_crime_costs_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_low_MMS)
l_outcomes_BUP_crime_costs_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_low_MMS)
ICER_crime_costs_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_low_MMS, outcomes_int = l_outcomes_BUP_crime_costs_low_MMS)
# TS
l_outcomes_MET_crime_costs_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_low_TS)
l_outcomes_BUP_crime_costs_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_low_TS)
ICER_crime_costs_low_TS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_low_TS, outcomes_int = l_outcomes_BUP_crime_costs_low_TS)

# High
# MMS
l_outcomes_MET_crime_costs_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_high_MMS)
l_outcomes_BUP_crime_costs_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_high_MMS)
ICER_crime_costs_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_high_MMS, outcomes_int = l_outcomes_BUP_crime_costs_high_MMS)
# TS
l_outcomes_MET_crime_costs_high_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_high_TS)
l_outcomes_BUP_crime_costs_high_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_high_TS)
ICER_crime_costs_high_TS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_high_TS, outcomes_int = l_outcomes_BUP_crime_costs_high_TS)

# Alternative (Krebs et al. 2014)
# MMS
l_outcomes_MET_crime_costs_alt_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_MMS)
l_outcomes_BUP_crime_costs_alt_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_MMS)
ICER_crime_costs_alt_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_alt_MMS, outcomes_int = l_outcomes_BUP_crime_costs_alt_MMS)
# TS
l_outcomes_MET_crime_costs_alt_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_TS)
l_outcomes_BUP_crime_costs_alt_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_TS)
ICER_crime_costs_alt_TS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_alt_TS, outcomes_int = l_outcomes_BUP_crime_costs_alt_TS)

# Reduced (equal across treatments)
# MMS
l_outcomes_MET_crime_costs_reduced_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_reduced_MMS)
l_outcomes_BUP_crime_costs_reduced_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_reduced_MMS)
ICER_crime_costs_reduced_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_reduced_MMS, outcomes_int = l_outcomes_BUP_crime_costs_reduced_MMS)
# TS
l_outcomes_MET_crime_costs_reduced_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_reduced_TS)
l_outcomes_BUP_crime_costs_reduced_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_reduced_TS)
ICER_crime_costs_reduced_TS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_reduced_TS, outcomes_int = l_outcomes_BUP_crime_costs_reduced_TS)

#############
### QALYs ###
#############
# Low
# MMS
l_outcomes_MET_qalys_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_low_MMS)
l_outcomes_BUP_qalys_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_low_MMS)
ICER_qalys_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_low_MMS, outcomes_int = l_outcomes_BUP_qalys_low_MMS)
# TS
l_outcomes_MET_qalys_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_low_TS)
l_outcomes_BUP_qalys_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_low_TS)
ICER_qalys_low_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_low_TS, outcomes_int = l_outcomes_BUP_qalys_low_TS)

# High
# MMS
l_outcomes_MET_qalys_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_high_MMS)
l_outcomes_BUP_qalys_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_high_MMS)
ICER_qalys_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_high_MMS, outcomes_int = l_outcomes_BUP_qalys_high_MMS)
# TS
l_outcomes_MET_qalys_high_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_high_TS)
l_outcomes_BUP_qalys_high_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_high_TS)
ICER_qalys_high_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_high_TS, outcomes_int = l_outcomes_BUP_qalys_high_TS)

## Reduced (EQ-5D-5L) ##
# MMS
l_outcomes_MET_qalys_reduced_eq_5d_5l_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_reduced_eq_5d_5l_MMS)
l_outcomes_BUP_qalys_reduced_eq_5d_5l_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_reduced_eq_5d_5l_MMS)
ICER_qalys_reduced_eq_5d_5l_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_reduced_eq_5d_5l_MMS, outcomes_int = l_outcomes_BUP_qalys_reduced_eq_5d_5l_MMS)
# TS
l_outcomes_MET_qalys_reduced_eq_5d_5l_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_reduced_eq_5d_5l_TS)
l_outcomes_BUP_qalys_reduced_eq_5d_5l_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_reduced_eq_5d_5l_TS)
ICER_qalys_reduced_eq_5d_5l_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_reduced_eq_5d_5l_TS, outcomes_int = l_outcomes_BUP_qalys_reduced_eq_5d_5l_TS)

## Alternative (EQ-5D-3L) ##
# MMS
l_outcomes_MET_qalys_eq_5d_3l_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_eq_5d_3l_MMS)
l_outcomes_BUP_qalys_eq_5d_3l_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_eq_5d_3l_MMS)
ICER_qalys_eq_5d_3l_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_eq_5d_3l_MMS, outcomes_int = l_outcomes_BUP_qalys_eq_5d_3l_MMS)
# TS
l_outcomes_MET_qalys_eq_5d_3l_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_eq_5d_3l_TS)
l_outcomes_BUP_qalys_eq_5d_3l_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_eq_5d_3l_TS)
ICER_qalys_eq_5d_3l_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_eq_5d_3l_TS, outcomes_int = l_outcomes_BUP_qalys_eq_5d_3l_TS)

## Alternative (HUI-3) ##
# MMS
l_outcomes_MET_qalys_hui_3_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_hui_3_MMS)
l_outcomes_BUP_qalys_hui_3_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_hui_3_MMS)
ICER_qalys_hui_3_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_hui_3_MMS, outcomes_int = l_outcomes_BUP_qalys_hui_3_MMS)
# TS
l_outcomes_MET_qalys_hui_3_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_hui_3_TS)
l_outcomes_BUP_qalys_hui_3_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_hui_3_TS)
ICER_qalys_hui_3_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_hui_3_TS, outcomes_int = l_outcomes_BUP_qalys_hui_3_TS)

## Alternative overdose ##
# MMS
l_outcomes_MET_qalys_odn_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_odn_low_MMS)
l_outcomes_BUP_qalys_odn_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_odn_low_MMS)
ICER_qalys_odn_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_odn_low_MMS, outcomes_int = l_outcomes_BUP_qalys_odn_low_MMS)
# TS
l_outcomes_MET_qalys_odn_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_odn_low_TS)
l_outcomes_BUP_qalys_odn_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_odn_low_TS)
ICER_qalys_odn_low_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_odn_low_TS, outcomes_int = l_outcomes_BUP_qalys_odn_low_TS)

################
### Overdose ###
################

### Overdose fatality ###
## Witnessed OD ##
# Low
l_outcomes_MET_witness_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_low)
l_outcomes_BUP_witness_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_low)
ICER_witness_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_witness_low_MMS, outcomes_int = l_outcomes_BUP_witness_low_MMS)

l_outcomes_MET_witness_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_low)
l_outcomes_BUP_witness_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_low)
ICER_witness_low_TS <- ICER(outcomes_comp = l_outcomes_MET_witness_low_TS, outcomes_int = l_outcomes_BUP_witness_low_TS)

# High
l_outcomes_MET_witness_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_high)
l_outcomes_BUP_witness_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_high)
ICER_witness_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_witness_high_MMS, outcomes_int = l_outcomes_BUP_witness_high_MMS)

l_outcomes_MET_witness_high_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_high)
l_outcomes_BUP_witness_high_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_high)
ICER_witness_high_TS <- ICER(outcomes_comp = l_outcomes_MET_witness_high_TS, outcomes_int = l_outcomes_BUP_witness_high_TS)

## Naloxone prevalence ##
# Low
l_outcomes_MET_NX_prob_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_prob_low)
l_outcomes_BUP_NX_prob_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_prob_low)
ICER_NX_prob_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_NX_prob_low_MMS, outcomes_int = l_outcomes_BUP_NX_prob_low_MMS)

l_outcomes_MET_NX_prob_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_prob_low)
l_outcomes_BUP_NX_prob_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_prob_low)
ICER_NX_prob_low_TS <- ICER(outcomes_comp = l_outcomes_MET_NX_prob_low_TS, outcomes_int = l_outcomes_BUP_NX_prob_low_TS)

# High
l_outcomes_MET_NX_prob_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_prob_high)
l_outcomes_BUP_NX_prob_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_prob_high)
ICER_NX_prob_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_NX_prob_high_MMS, outcomes_int = l_outcomes_BUP_NX_prob_high_MMS)

l_outcomes_MET_NX_prob_high_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_prob_high)
l_outcomes_BUP_NX_prob_high_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_prob_high)
ICER_NX_prob_high_TS <- ICER(outcomes_comp = l_outcomes_MET_NX_prob_high_TS, outcomes_int = l_outcomes_BUP_NX_prob_high_TS)

## Naloxone effectiveness ##
# Low
l_outcomes_MET_NX_success_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_success_low)
l_outcomes_BUP_NX_success_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_success_low)
ICER_NX_success_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_NX_success_low_MMS, outcomes_int = l_outcomes_BUP_NX_success_low_MMS)

l_outcomes_MET_NX_success_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_success_low)
l_outcomes_BUP_NX_success_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_success_low)
ICER_NX_success_low_TS <- ICER(outcomes_comp = l_outcomes_MET_NX_success_low_TS, outcomes_int = l_outcomes_BUP_NX_success_low_TS)

# High
l_outcomes_MET_NX_success_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_success_high)
l_outcomes_BUP_NX_success_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_success_high)
ICER_NX_success_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_NX_success_high_MMS, outcomes_int = l_outcomes_BUP_NX_success_high_MMS)

l_outcomes_MET_NX_success_high_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_success_high)
l_outcomes_BUP_NX_success_high_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_NX_success_high)
ICER_NX_success_high_TS <- ICER(outcomes_comp = l_outcomes_MET_NX_success_high_TS, outcomes_int = l_outcomes_BUP_NX_success_high_TS)

## Fatal overdose risk ##
# Low
l_outcomes_MET_fatal_OD_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fatal_OD_low)
l_outcomes_BUP_fatal_OD_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fatal_OD_low)
ICER_fatal_OD_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_fatal_OD_low_MMS, outcomes_int = l_outcomes_BUP_fatal_OD_low_MMS)

l_outcomes_MET_fatal_OD_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fatal_OD_low)
l_outcomes_BUP_fatal_OD_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fatal_OD_low)
ICER_fatal_OD_low_TS <- ICER(outcomes_comp = l_outcomes_MET_fatal_OD_low_TS, outcomes_int = l_outcomes_BUP_fatal_OD_low_TS)

# High
l_outcomes_MET_fatal_OD_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fatal_OD_high)
l_outcomes_BUP_fatal_OD_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fatal_OD_high)
ICER_fatal_OD_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_fatal_OD_high_MMS, outcomes_int = l_outcomes_BUP_fatal_OD_high_MMS)

l_outcomes_MET_fatal_OD_high_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fatal_OD_high)
l_outcomes_BUP_fatal_OD_high_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fatal_OD_high)
ICER_fatal_OD_high_TS <- ICER(outcomes_comp = l_outcomes_MET_fatal_OD_high_TS, outcomes_int = l_outcomes_BUP_fatal_OD_high_TS)


### Non-fatal overdose ###
## Fentanyl prevalence ##
# Low
l_outcomes_MET_fent_exp_2020_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_exp_2020_low)
l_outcomes_BUP_fent_exp_2020_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_exp_2020_low)
ICER_fent_exp_2020_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_fent_exp_2020_low_MMS, outcomes_int = l_outcomes_BUP_fent_exp_2020_low_MMS)

l_outcomes_MET_fent_exp_2020_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_exp_2020_low)
l_outcomes_BUP_fent_exp_2020_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_exp_2020_low)
ICER_fent_exp_2020_low_TS <- ICER(outcomes_comp = l_outcomes_MET_fent_exp_2020_low_TS, outcomes_int = l_outcomes_BUP_fent_exp_2020_low_TS)

# High
l_outcomes_MET_fent_exp_2020_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_exp_2020_high)
l_outcomes_BUP_fent_exp_2020_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_exp_2020_high)
ICER_fent_exp_2020_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_fent_exp_2020_high_MMS, outcomes_int = l_outcomes_BUP_fent_exp_2020_high_MMS)

l_outcomes_MET_fent_exp_2020_high_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_exp_2020_high)
l_outcomes_BUP_fent_exp_2020_high_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_exp_2020_high)
ICER_fent_exp_2020_high_TS <- ICER(outcomes_comp = l_outcomes_MET_fent_exp_2020_high_TS, outcomes_int = l_outcomes_BUP_fent_exp_2020_high_TS)

## Reduction in fentanyl exposure for non-injection v. injection ##
# Low
l_outcomes_MET_ni_fent_reduction_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_ni_fent_reduction_low)
l_outcomes_BUP_ni_fent_reduction_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_ni_fent_reduction_low)
ICER_ni_fent_reduction_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_ni_fent_reduction_low_MMS, outcomes_int = l_outcomes_BUP_ni_fent_reduction_low_MMS)

l_outcomes_MET_ni_fent_reduction_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_ni_fent_reduction_low)
l_outcomes_BUP_ni_fent_reduction_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_ni_fent_reduction_low)
ICER_ni_fent_reduction_low_TS <- ICER(outcomes_comp = l_outcomes_MET_ni_fent_reduction_low_TS, outcomes_int = l_outcomes_BUP_ni_fent_reduction_low_TS)

# High
l_outcomes_MET_ni_fent_reduction_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_ni_fent_reduction_high)
l_outcomes_BUP_ni_fent_reduction_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_ni_fent_reduction_high)
ICER_ni_fent_reduction_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_ni_fent_reduction_high_MMS, outcomes_int = l_outcomes_BUP_ni_fent_reduction_high_MMS)

l_outcomes_MET_ni_fent_reduction_high_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_ni_fent_reduction_high)
l_outcomes_BUP_ni_fent_reduction_high_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_ni_fent_reduction_high)
ICER_ni_fent_reduction_high_TS <- ICER(outcomes_comp = l_outcomes_MET_ni_fent_reduction_high_TS, outcomes_int = l_outcomes_BUP_ni_fent_reduction_high_TS)

## Fent OD multiplier ##
# Low
l_outcomes_MET_fent_OD_mult_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_OD_mult_low)
l_outcomes_BUP_fent_OD_mult_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_OD_mult_low)
ICER_fent_OD_mult_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_fent_OD_mult_low_MMS, outcomes_int = l_outcomes_BUP_fent_OD_mult_low_MMS)

l_outcomes_MET_fent_OD_mult_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_OD_mult_low)
l_outcomes_BUP_fent_OD_mult_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_OD_mult_low)
ICER_fent_OD_mult_low_TS <- ICER(outcomes_comp = l_outcomes_MET_fent_OD_mult_low_TS, outcomes_int = l_outcomes_BUP_fent_OD_mult_low_TS)

# High
l_outcomes_MET_fent_OD_mult_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_OD_mult_high)
l_outcomes_BUP_fent_OD_mult_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_OD_mult_high)
ICER_fent_OD_mult_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_fent_OD_mult_high_MMS, outcomes_int = l_outcomes_BUP_fent_OD_mult_high_MMS)

l_outcomes_MET_fent_OD_mult_high_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_OD_mult_high)
l_outcomes_BUP_fent_OD_mult_high_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fent_OD_mult_high)
ICER_fent_OD_mult_high_TS <- ICER(outcomes_comp = l_outcomes_MET_fent_OD_mult_high_TS, outcomes_int = l_outcomes_BUP_fent_OD_mult_high_TS)

# BUP OD multiplier

# MET OD multiplier

# REL OD multiplier

##############################
### Cohort characteristics ###
##############################
# Starting age

# Male %

## Outcomes ##
# Costs



# QALYs

###################
#### Data Prep ####
###################
# Plotting
# Function to compare against baseline
improvement_formatter <- formatter("span", 
                                   style = x ~ style(font.weight = "bold", 
                                                     color = ifelse(x > 0, customGreen, ifelse(x < 0, customRed, "black"))), 
                                   x ~ icontext(ifelse(x>0, "arrow-up", "arrow-down"), x)
)

# set colours
customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customRed = "#ff7f7f"

################
### Baseline ###
################
df_baseline_MMS <- data.frame(ICER_MMS$df_incremental, ICER_MMS$df_icer)
df_baseline_TS  <- data.frame(ICER_TS$df_incremental, ICER_TS$df_icer)

###################
### Crime costs ###
###################
df_crime_costs_baseline_MMS <- data.frame(ICER_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                          ICER_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_MMS$df_icer$n_icer_TOTAL_1yr, 
                                          ICER_MMS$df_icer$n_icer_TOTAL_5yr, ICER_MMS$df_icer$n_icer_TOTAL_10yr)
df_crime_costs_baseline_TS <- data.frame(ICER_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_TS$df_icer$n_icer_TOTAL_1yr)

df_crime_costs_reduced_MMS <- data.frame(ICER_crime_costs_reduced_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_reduced_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                         ICER_crime_costs_reduced_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_crime_costs_reduced_MMS$df_icer$n_icer_TOTAL_1yr, 
                                         ICER_crime_costs_reduced_MMS$df_icer$n_icer_TOTAL_5yr, ICER_crime_costs_reduced_MMS$df_icer$n_icer_TOTAL_10yr)
df_crime_costs_reduced_TS <- data.frame(ICER_crime_costs_reduced_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_reduced_TS$df_icer$n_icer_TOTAL_1yr)

df_crime_costs_low_MMS <- data.frame(ICER_crime_costs_low_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_low_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                     ICER_crime_costs_low_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_crime_costs_low_MMS$df_icer$n_icer_TOTAL_1yr, 
                                     ICER_crime_costs_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_crime_costs_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_crime_costs_low_TS <- data.frame(ICER_crime_costs_low_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_low_TS$df_icer$n_icer_TOTAL_1yr)

df_crime_costs_high_MMS <- data.frame(ICER_crime_costs_high_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_high_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                      ICER_crime_costs_high_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_crime_costs_high_MMS$df_icer$n_icer_TOTAL_1yr, 
                                      ICER_crime_costs_high_MMS$df_icer$n_icer_TOTAL_5yr, ICER_crime_costs_high_MMS$df_icer$n_icer_TOTAL_10yr)
df_crime_costs_high_TS <- data.frame(ICER_crime_costs_high_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_high_TS$df_icer$n_icer_TOTAL_1yr)

df_crime_costs_alt_MMS <- data.frame(ICER_crime_costs_alt_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_alt_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                     ICER_crime_costs_alt_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_crime_costs_alt_MMS$df_icer$n_icer_TOTAL_1yr, 
                                     ICER_crime_costs_alt_MMS$df_icer$n_icer_TOTAL_5yr, ICER_crime_costs_alt_MMS$df_icer$n_icer_TOTAL_10yr)
df_crime_costs_alt_TS <- data.frame(ICER_crime_costs_alt_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_alt_TS$df_icer$n_icer_TOTAL_1yr)

df_crime_costs_reduced_MMS <- data.frame(ICER_crime_costs_reduced_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_reduced_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                     ICER_crime_costs_reduced_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_crime_costs_reduced_MMS$df_icer$n_icer_TOTAL_1yr, 
                                     ICER_crime_costs_reduced_MMS$df_icer$n_icer_TOTAL_5yr, ICER_crime_costs_reduced_MMS$df_icer$n_icer_TOTAL_10yr)
df_crime_costs_reduced_TS <- data.frame(ICER_crime_costs_reduced_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_reduced_TS$df_icer$n_icer_TOTAL_1yr)

colnames(df_crime_costs_baseline_MMS) <- colnames(df_crime_costs_low_MMS) <- colnames(df_crime_costs_high_MMS) <- colnames(df_crime_costs_alt_MMS) <- colnames(df_crime_costs_reduced_MMS) <- c("inc_costs_1yr", "inc_costs_5yr", "inc_costs_10yr", "icer_1yr", "icer_5yr", "icer_10yr")
colnames(df_crime_costs_baseline_TS) <- colnames(df_crime_costs_low_TS) <- colnames(df_crime_costs_high_TS) <- colnames(df_crime_costs_alt_TS) <- colnames(df_crime_costs_reduced_TS) <- c("inc_costs_1yr", "icer_1yr")

df_crime_costs <- rbind(df_crime_costs_baseline_MMS, df_crime_costs_low_MMS, df_crime_costs_high_MMS, df_crime_costs_alt_MMS, df_crime_costs_reduced_MMS)
df_crime_costs_TS <- rbind(df_crime_costs_baseline_TS, df_crime_costs_low_TS, df_crime_costs_high_TS, df_crime_costs_alt_TS, df_crime_costs_reduced_TS)

df_crime_costs <- data.frame("Scenario" = c("Baseline", "Low", "High", "Alternative", "Combined Tx"), df_crime_costs)
df_crime_costs_TS <- data.frame("Scenario" = c("Baseline", "Low", "High", "Alternative", "Combined Tx"), df_crime_costs_TS)

# Custom table output
crime_costs_palette <- brewer.pal(3,"BrBG")

table_crime_costs <- df_crime_costs %>%
  mutate(`Incremental Costs (1-year)` = accounting(inc_costs_1yr, 0),
         `Incremental Costs (5-year)` = accounting(inc_costs_5yr, 0),
         `Incremental Costs (10-year)` = accounting(inc_costs_10yr, 0),
         `Difference v. Baseline (1-year)` = `Incremental Costs (1-year)` - `Incremental Costs (1-year)`[`Scenario` == "Baseline"],
         `Difference v. Baseline (5-year)` = `Incremental Costs (5-year)` - `Incremental Costs (5-year)`[`Scenario` == "Baseline"],
         `Difference v. Baseline (10-year)` = `Incremental Costs (10-year)` - `Incremental Costs (10-year)`[`Scenario` == "Baseline"],
         `ICER (1-year)` = accounting(icer_1yr, 0),
         `ICER (5-year)` = accounting(icer_5yr, 0),
         `ICER (10-year)` = accounting(icer_10yr, 0)) %>%
  select(c(`Scenario`, `Incremental Costs (1-year)`, `Incremental Costs (5-year)`, `Incremental Costs (10-year)`, 
           `Difference v. Baseline (1-year)`, `Difference v. Baseline (5-year)`, `Difference v. Baseline (10-year)`, `ICER (1-year)`, `ICER (5-year)`, `ICER (10-year)`))

ftable_crime_costs_out <- formattable(table_crime_costs, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental Costs (1-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[1]),
  `Incremental Costs (5-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[1]),
  `Incremental Costs (10-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[1]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `Difference v. Baseline (5-year)` = improvement_formatter,
  `Difference v. Baseline (10-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[3]),
  `ICER (5-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[3]),
  `ICER (10-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[3])))

table_crime_costs_TS <- df_crime_costs_TS %>%
  mutate(`Incremental Costs (1-year)` = accounting(inc_costs_1yr, 0),
         `Difference v. Baseline (1-year)` = `Incremental Costs (1-year)` - `Incremental Costs (1-year)`[`Scenario` == "Baseline"],
         `ICER (1-year)` = accounting(icer_1yr, 0)) %>%
  select(c(`Scenario`, `Incremental Costs (1-year)`, `Difference v. Baseline (1-year)`, `ICER (1-year)`))

ftable_crime_costs_TS_out <- formattable(table_crime_costs_TS, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental Costs (1-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[1]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[3])))

# Output
save(ftable_crime_costs_out, 
     file = "outputs/DSA/Modified Model Specification/ftable_crime_costs_out.RData")
save(ftable_crime_costs_TS_out, 
     file = "outputs/DSA/Trial Specification/ftable_crime_costs_TS_out.RData")



#############
### QALYs ###
#############
# Baseline
df_qalys_baseline_MMS <- data.frame(ICER_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                                    ICER_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_MMS$df_icer$n_icer_TOTAL_1yr, 
                                    ICER_MMS$df_icer$n_icer_TOTAL_5yr, ICER_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_baseline_TS <- data.frame(ICER_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_TS$df_icer$n_icer_TOTAL_1yr)

# Combined Tx
df_qalys_reduced_eq_5d_5l_MMS <- data.frame(ICER_qalys_reduced_eq_5d_5l_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_reduced_eq_5d_5l_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                                            ICER_qalys_reduced_eq_5d_5l_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_qalys_reduced_eq_5d_5l_MMS$df_icer$n_icer_TOTAL_1yr, 
                                            ICER_qalys_reduced_eq_5d_5l_MMS$df_icer$n_icer_TOTAL_5yr, ICER_qalys_reduced_eq_5d_5l_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_reduced_eq_5d_5l_TS <- data.frame(ICER_qalys_reduced_eq_5d_5l_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_reduced_eq_5d_5l_TS$df_icer$n_icer_TOTAL_1yr)

# Low
df_qalys_low_MMS <- data.frame(ICER_qalys_low_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_low_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                               ICER_qalys_low_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_qalys_low_MMS$df_icer$n_icer_TOTAL_1yr, 
                               ICER_qalys_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_qalys_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_low_TS <- data.frame(ICER_qalys_low_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_low_TS$df_icer$n_icer_TOTAL_1yr)

# High
df_qalys_high_MMS <- data.frame(ICER_qalys_high_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_high_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                                ICER_qalys_high_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_qalys_high_MMS$df_icer$n_icer_TOTAL_1yr, 
                                ICER_qalys_high_MMS$df_icer$n_icer_TOTAL_5yr, ICER_qalys_high_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_high_TS <- data.frame(ICER_qalys_high_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_high_TS$df_icer$n_icer_TOTAL_1yr)

# Alt (EQ-5d-3L)
df_qalys_eq_5d_3l_MMS <- data.frame(ICER_qalys_eq_5d_3l_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_eq_5d_3l_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                                    ICER_qalys_eq_5d_3l_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_qalys_eq_5d_3l_MMS$df_icer$n_icer_TOTAL_1yr, 
                                    ICER_qalys_eq_5d_3l_MMS$df_icer$n_icer_TOTAL_5yr, ICER_qalys_eq_5d_3l_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_eq_5d_3l_TS <- data.frame(ICER_qalys_eq_5d_3l_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_eq_5d_3l_TS$df_icer$n_icer_TOTAL_1yr)

# Alt (HUI-3)
df_qalys_hui_3_MMS <- data.frame(ICER_qalys_hui_3_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_hui_3_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                                 ICER_qalys_hui_3_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_qalys_hui_3_MMS$df_icer$n_icer_TOTAL_1yr, 
                                 ICER_qalys_hui_3_MMS$df_icer$n_icer_TOTAL_5yr, ICER_qalys_hui_3_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_hui_3_TS <- data.frame(ICER_qalys_hui_3_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_hui_3_TS$df_icer$n_icer_TOTAL_1yr)

# ODN low
df_qalys_odn_low_MMS <- data.frame(ICER_qalys_odn_low_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_odn_low_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                                 ICER_qalys_odn_low_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_qalys_odn_low_MMS$df_icer$n_icer_TOTAL_1yr, 
                                 ICER_qalys_odn_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_qalys_odn_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_odn_low_TS <- data.frame(ICER_qalys_odn_low_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_odn_low_TS$df_icer$n_icer_TOTAL_1yr)

colnames(df_qalys_baseline_MMS) <- colnames(df_qalys_reduced_eq_5d_5l_MMS) <- colnames(df_qalys_low_MMS) <- colnames(df_qalys_high_MMS) <- colnames(df_qalys_eq_5d_3l_MMS) <- colnames(df_qalys_hui_3_MMS) <- colnames(df_qalys_odn_low_MMS) <- c("inc_qalys_1yr", "inc_qalys_5yr", "inc_qalys_10yr", "icer_1yr", "icer_5yr", "icer_10yr")
colnames(df_qalys_baseline_TS) <- colnames(df_qalys_reduced_eq_5d_5l_TS) <- colnames(df_qalys_low_TS) <- colnames(df_qalys_high_TS) <- colnames(df_qalys_eq_5d_3l_TS) <- colnames(df_qalys_hui_3_TS) <- colnames(df_qalys_odn_low_TS) <- c("inc_qalys_1yr", "icer_1yr")

df_qalys <- rbind(df_qalys_baseline_MMS, df_qalys_reduced_eq_5d_5l_MMS, df_qalys_low_MMS, df_qalys_high_MMS, df_qalys_eq_5d_3l_MMS, df_qalys_hui_3_MMS, df_qalys_odn_low_MMS)
df_qalys_TS <- rbind(df_qalys_baseline_TS, df_qalys_reduced_eq_5d_5l_TS, df_qalys_low_TS, df_qalys_high_TS, df_qalys_eq_5d_3l_TS, df_qalys_hui_3_TS, df_qalys_odn_low_TS)

df_qalys <- data.frame("Scenario" = c("Baseline", "Combined Tx", "Low", "High", "EQ-5D-3L", "HUI-3", "ODN (low)"), df_qalys)
df_qalys_TS <- data.frame("Scenario" = c("Baseline", "Combined Tx", "Low", "High", "EQ-5D-3L", "HUI-3", "ODN (low)"), df_qalys_TS)

# Custom table output
qaly_palette <- brewer.pal(3,"PuOr")

table_qalys_MMS <- df_qalys_MMS %>%
  mutate(`Incremental QALYs (1-year)` = round(inc_qalys_1yr, 3),
         `Incremental QALYs (5-year)` = round(inc_qalys_5yr, 3),
         `Incremental QALYs (10-year)` = round(inc_qalys_10yr, 3),
         `Difference v. Baseline (1-year)` = round(`Incremental QALYs (1-year)` - `Incremental QALYs (1-year)`[`Scenario` == "Baseline"], 3),
         `Difference v. Baseline (5-year)` = round(`Incremental QALYs (5-year)` - `Incremental QALYs (5-year)`[`Scenario` == "Baseline"], 3),
         `Difference v. Baseline (10-year)` = round(`Incremental QALYs (10-year)` - `Incremental QALYs (10-year)`[`Scenario` == "Baseline"], 3),
         `ICER (1-year)` = accounting(icer_1yr, 0),
         `ICER (5-year)` = accounting(icer_5yr, 0),
         `ICER (10-year)` = accounting(icer_10yr, 0)) %>%
  select(c(`Scenario`, `Incremental QALYs (1-year)`, `Incremental QALYs (5-year)`, `Incremental QALYs (10-year)`, 
           `Difference v. Baseline (1-year)`, `Difference v. Baseline (5-year)`, `Difference v. Baseline (10-year)`,
           `ICER (1-year)`, `ICER (5-year)`, `ICER (10-year)`)) #%>%

ftable_qalys_MMS_out <- formattable(table_qalys_MMS, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
    `Incremental QALYs (1-year)` = color_tile(qaly_palette[1], qaly_palette[2]),
  `Incremental QALYs (5-year)` = color_tile(qaly_palette[1], qaly_palette[2]),
  `Incremental QALYs (10-year)` = color_tile(qaly_palette[1], qaly_palette[2]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `Difference v. Baseline (5-year)` = improvement_formatter,
  `Difference v. Baseline (10-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(qaly_palette[2], qaly_palette[3]),
  `ICER (5-year)` = color_tile(qaly_palette[2], qaly_palette[3]),
  `ICER (10-year)` = color_tile(qaly_palette[2], qaly_palette[3])))

table_qalys_TS <- df_qalys_TS %>%
  mutate(`Incremental QALYs (1-year)` = round(inc_qalys_1yr, 3),
         `Difference v. Baseline (1-year)` = round(`Incremental QALYs (1-year)` - `Incremental QALYs (1-year)`[`Scenario` == "Baseline"], 3),
         `ICER (1-year)` = accounting(icer_1yr, 0)) %>%
  select(c(`Scenario`, `Incremental QALYs (1-year)`, `Difference v. Baseline (1-year)`, `ICER (1-year)`)) #%>%

ftable_qalys_TS_out <- formattable(table_qalys_TS, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental QALYs (1-year)` = color_tile(qaly_palette[1], qaly_palette[2]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(qaly_palette[2], qaly_palette[3])))

# Output
save(ftable_qalys_MMS_out, 
     file = "outputs/DSA/Modified Model Specification/ftable_qalys_MMS_out.RData")
save(ftable_qalys_TS_out, 
     file = "outputs/DSA/Trial Specification/ftable_qalys_TS_out.RData")

################
### Overdose ###
################
# Baseline
df_overdose_baseline_MMS <- data.frame(ICER_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                       ICER_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                       ICER_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                       ICER_MMS$df_icer$n_icer_TOTAL_1yr, ICER_MMS$df_icer$n_icer_TOTAL_5yr, ICER_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_baseline_TS <- data.frame(ICER_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_TS$df_icer$n_icer_TOTAL_1yr)

## Fatal overdose ##
# Witness probability
df_overdose_witness_low_MMS <- data.frame(ICER_witness_low_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_witness_low_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                          ICER_witness_low_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_witness_low_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                          ICER_witness_low_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_witness_low_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                          ICER_witness_low_MMS$df_icer$n_icer_TOTAL_1yr, ICER_witness_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_witness_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_witness_low_TS <- data.frame(ICER_witness_low_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_witness_low_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_witness_low_TS$df_icer$n_icer_TOTAL_1yr)


df_overdose_witness_high_MMS <- data.frame(ICER_witness_high_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_witness_high_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                           ICER_witness_high_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_witness_high_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                           ICER_witness_high_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_witness_high_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                           ICER_witness_high_MMS$df_icer$n_icer_TOTAL_1yr, ICER_witness_high_MMS$df_icer$n_icer_TOTAL_5yr, ICER_witness_high_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_witness_high_TS <- data.frame(ICER_witness_high_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_witness_high_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_witness_high_TS$df_icer$n_icer_TOTAL_1yr)

# Naloxone probability
df_overdose_NX_prob_low_MMS <- data.frame(ICER_NX_prob_low_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_NX_prob_low_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                          ICER_NX_prob_low_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_NX_prob_low_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                          ICER_NX_prob_low_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_NX_prob_low_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                          ICER_NX_prob_low_MMS$df_icer$n_icer_TOTAL_1yr, ICER_NX_prob_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_NX_prob_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_NX_prob_low_TS <- data.frame(ICER_NX_prob_low_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_NX_prob_low_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_NX_prob_low_TS$df_icer$n_icer_TOTAL_1yr)


df_overdose_NX_prob_high_MMS <- data.frame(ICER_NX_prob_high_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_NX_prob_high_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                           ICER_NX_prob_high_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_NX_prob_high_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                           ICER_NX_prob_high_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_NX_prob_high_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                           ICER_NX_prob_high_MMS$df_icer$n_icer_TOTAL_1yr, ICER_NX_prob_high_MMS$df_icer$n_icer_TOTAL_5yr, ICER_NX_prob_high_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_NX_prob_high_TS <- data.frame(ICER_NX_prob_high_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_NX_prob_high_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_NX_prob_high_TS$df_icer$n_icer_TOTAL_1yr)

# Naloxone success
df_overdose_NX_success_low_MMS <- data.frame(ICER_NX_success_low_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_NX_success_low_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                          ICER_NX_success_low_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_NX_success_low_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                          ICER_NX_success_low_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_NX_success_low_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                          ICER_NX_success_low_MMS$df_icer$n_icer_TOTAL_1yr, ICER_NX_success_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_NX_success_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_NX_success_low_TS <- data.frame(ICER_NX_success_low_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_NX_success_low_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_NX_success_low_TS$df_icer$n_icer_TOTAL_1yr)


df_overdose_NX_success_high_MMS <- data.frame(ICER_NX_success_high_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_NX_success_high_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                           ICER_NX_success_high_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_NX_success_high_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                           ICER_NX_success_high_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_NX_success_high_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                           ICER_NX_success_high_MMS$df_icer$n_icer_TOTAL_1yr, ICER_NX_success_high_MMS$df_icer$n_icer_TOTAL_5yr, ICER_NX_success_high_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_NX_success_high_TS <- data.frame(ICER_NX_success_high_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_NX_success_high_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_NX_success_high_TS$df_icer$n_icer_TOTAL_1yr)

# Fatal overdose rate
df_overdose_fatal_OD_low_MMS <- data.frame(ICER_fatal_OD_low_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fatal_OD_low_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                             ICER_fatal_OD_low_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_fatal_OD_low_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                             ICER_fatal_OD_low_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_fatal_OD_low_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                             ICER_fatal_OD_low_MMS$df_icer$n_icer_TOTAL_1yr, ICER_fatal_OD_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_fatal_OD_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_fatal_OD_low_TS <- data.frame(ICER_fatal_OD_low_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fatal_OD_low_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_fatal_OD_low_TS$df_icer$n_icer_TOTAL_1yr)

df_overdose_fatal_OD_high_MMS <- data.frame(ICER_fatal_OD_high_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fatal_OD_high_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                              ICER_fatal_OD_high_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_fatal_OD_high_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                              ICER_fatal_OD_high_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_fatal_OD_high_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                              ICER_fatal_OD_high_MMS$df_icer$n_icer_TOTAL_1yr, ICER_fatal_OD_high_MMS$df_icer$n_icer_TOTAL_5yr, ICER_fatal_OD_high_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_fatal_OD_high_TS <- data.frame(ICER_fatal_OD_high_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fatal_OD_high_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_fatal_OD_high_TS$df_icer$n_icer_TOTAL_1yr)

colnames(df_overdose_baseline_MMS) <- colnames(df_overdose_witness_low_MMS) <- colnames(df_overdose_witness_high_MMS) <- colnames(df_overdose_NX_prob_low_MMS) <- colnames(df_overdose_NX_prob_high_MMS) <- colnames(df_overdose_NX_success_low_MMS) <- colnames(df_overdose_NX_success_high_MMS) <- colnames(df_overdose_fatal_OD_low_MMS) <- colnames(df_overdose_fatal_OD_high_MMS) <- c("inc_costs_1yr", "inc_costs_5yr", "inc_costs_10yr", "inc_qalys_1yr", 
                                                                                                                                                                                                                                                                                                                                                                                            "inc_qalys_5yr", "inc_qalys_10yr", "icer_1yr", "icer_5yr", "icer_10yr")

colnames(df_overdose_baseline_TS) <- colnames(df_overdose_witness_low_TS) <- colnames(df_overdose_witness_high_TS) <- colnames(df_overdose_NX_prob_low_TS) <- colnames(df_overdose_NX_prob_high_TS) <- colnames(df_overdose_NX_success_low_TS) <- colnames(df_overdose_NX_success_high_TS) <- colnames(df_overdose_fatal_OD_low_TS) <- colnames(df_overdose_fatal_OD_high_TS) <- c("inc_costs_1yr", "inc_qalys_1yr", "icer_1yr")

df_overdose_fatal_MMS <- rbind(df_overdose_baseline_MMS, df_overdose_witness_low_MMS, df_overdose_witness_high_MMS, df_overdose_NX_prob_low_MMS, df_overdose_NX_prob_high_MMS, df_overdose_NX_success_low_MMS, df_overdose_NX_success_high_MMS, df_overdose_fatal_OD_low_MMS, df_overdose_fatal_OD_high_MMS)
df_overdose_fatal_TS <- rbind(df_overdose_baseline_TS, df_overdose_witness_low_TS, df_overdose_witness_high_TS, df_overdose_NX_prob_low_TS, df_overdose_NX_prob_high_TS, df_overdose_NX_success_low_TS, df_overdose_NX_success_high_TS, df_overdose_fatal_OD_low_TS, df_overdose_fatal_OD_high_TS)

df_overdose_fatal_MMS <- data.frame("Scenario" = c("Baseline", "Witness (low)", "Witness (high)", "NX prob (low)", "NX prob (high)", "NX success (low)", "NX success (high)", "Fatal OD rate (low)", "Fatal OD rate (high)"), df_overdose_fatal_MMS)
df_overdose_fatal_TS <- data.frame("Scenario" = c("Baseline", "Witness (low)", "Witness (high)", "NX prob (low)", "NX prob (high)", "NX success (low)", "NX success (high)", "Fatal OD rate (low)", "Fatal OD rate (high)"), df_overdose_fatal_TS)

# Custom table output
overdose_fatal_palette <- brewer.pal(3,"PuOr")

table_overdose_fatal_costs_MMS <- df_overdose_fatal_MMS %>%
  mutate(`Incremental Costs (1-year)` = accounting(inc_costs_1yr, 0),
         `Incremental Costs (5-year)` = accounting(inc_costs_5yr, 0),
         `Incremental Costs (10-year)` = accounting(inc_costs_10yr, 0),
         `Difference v. Baseline (1-year)` = `Incremental Costs (1-year)` - `Incremental Costs (1-year)`[`Scenario` == "Baseline"],
         `Difference v. Baseline (5-year)` = `Incremental Costs (5-year)` - `Incremental Costs (5-year)`[`Scenario` == "Baseline"],
         `Difference v. Baseline (10-year)` = `Incremental Costs (10-year)` - `Incremental Costs (10-year)`[`Scenario` == "Baseline"],
         `ICER (1-year)` = accounting(icer_1yr, 0),
         `ICER (5-year)` = accounting(icer_5yr, 0),
         `ICER (10-year)` = accounting(icer_10yr, 0)) %>%
  select(c(`Scenario`, `Incremental Costs (1-year)`, `Incremental Costs (5-year)`, `Incremental Costs (10-year)`, 
           `Difference v. Baseline (1-year)`, `Difference v. Baseline (5-year)`, `Difference v. Baseline (10-year)`,
           `ICER (1-year)`, `ICER (5-year)`, `ICER (10-year)`))

table_overdose_fatal_costs_TS <- df_overdose_fatal_TS %>%
  mutate(`Incremental Costs (1-year)` = accounting(inc_costs_1yr, 0),
         `Difference v. Baseline (1-year)` = `Incremental Costs (1-year)` - `Incremental Costs (1-year)`[`Scenario` == "Baseline"],
         `ICER (1-year)` = accounting(icer_1yr, 0)) %>%
  select(c(`Scenario`, `Incremental Costs (1-year)`, `Difference v. Baseline (1-year)`, `ICER (1-year)`))

table_overdose_fatal_qalys_MMS <- df_overdose_fatal_MMS %>%
  mutate(`Incremental QALYs (1-year)` = round(inc_qalys_1yr, 3),
         `Incremental QALYs (5-year)` = round(inc_qalys_5yr, 3),
         `Incremental QALYs (10-year)` = round(inc_qalys_10yr, 3),
         `Difference v. Baseline (1-year)` = round(`Incremental QALYs (1-year)` - `Incremental QALYs (1-year)`[`Scenario` == "Baseline"], 3),
         `Difference v. Baseline (5-year)` = round(`Incremental QALYs (5-year)` - `Incremental QALYs (5-year)`[`Scenario` == "Baseline"], 3),
         `Difference v. Baseline (10-year)` = round(`Incremental QALYs (10-year)` - `Incremental QALYs (10-year)`[`Scenario` == "Baseline"], 3),
         `ICER (1-year)` = accounting(icer_1yr, 0),
         `ICER (5-year)` = accounting(icer_5yr, 0),
         `ICER (10-year)` = accounting(icer_10yr, 0)) %>%
  select(c(`Scenario`, `Incremental QALYs (1-year)`, `Incremental QALYs (5-year)`, `Incremental QALYs (10-year)`, 
           `Difference v. Baseline (1-year)`, `Difference v. Baseline (5-year)`, `Difference v. Baseline (10-year)`,
           `ICER (1-year)`, `ICER (5-year)`, `ICER (10-year)`))

table_overdose_fatal_qalys_TS <- df_overdose_fatal_TS %>%
  mutate(`Incremental QALYs (1-year)` = round(inc_qalys_1yr, 3),
         `Difference v. Baseline (1-year)` = round(`Incremental QALYs (1-year)` - `Incremental QALYs (1-year)`[`Scenario` == "Baseline"], 3),
         `ICER (1-year)` = accounting(icer_1yr, 0)) %>%
  select(c(`Scenario`, `Incremental QALYs (1-year)`, `Difference v. Baseline (1-year)`, `ICER (1-year)`))

ftable_overdose_fatal_qalys_MMS_out <- formattable(table_overdose_fatal_qalys_MMS, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental QALYs (1-year)` = color_tile(overdose_fatal_palette[1], overdose_fatal_palette[2]),
  `Incremental QALYs (5-year)` = color_tile(overdose_fatal_palette[1], overdose_fatal_palette[2]),
  `Incremental QALYs (10-year)` = color_tile(overdose_fatal_palette[1], overdose_fatal_palette[2]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `Difference v. Baseline (5-year)` = improvement_formatter,
  `Difference v. Baseline (10-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(overdose_fatal_palette[2], overdose_fatal_palette[3]),
  `ICER (5-year)` = color_tile(overdose_fatal_palette[2], overdose_fatal_palette[3]),
  `ICER (10-year)` = color_tile(overdose_fatal_palette[2], overdose_fatal_palette[3])))

ftable_overdose_fatal_qalys_TS_out <- formattable(table_overdose_fatal_qalys_TS, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental QALYs (1-year)` = color_tile(overdose_fatal_palette[1], overdose_fatal_palette[2]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(overdose_fatal_palette[2], overdose_fatal_palette[3])))

ftable_overdose_fatal_costs_MMS_out <- formattable(table_overdose_fatal_costs_MMS, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental Costs (1-year)` = color_tile(overdose_fatal_palette[1], overdose_fatal_palette[2]),
  `Incremental Costs (5-year)` = color_tile(overdose_fatal_palette[1], overdose_fatal_palette[2]),
  `Incremental Costs (10-year)` = color_tile(overdose_fatal_palette[1], overdose_fatal_palette[2]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `Difference v. Baseline (5-year)` = improvement_formatter,
  `Difference v. Baseline (10-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(overdose_fatal_palette[2], overdose_fatal_palette[3]),
  `ICER (5-year)` = color_tile(overdose_fatal_palette[2], overdose_fatal_palette[3]),
  `ICER (10-year)` = color_tile(overdose_fatal_palette[2], overdose_fatal_palette[3])))

ftable_overdose_fatal_costs_TS_out <- formattable(table_overdose_fatal_costs_TS, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental Costs (1-year)` = color_tile(overdose_fatal_palette[1], overdose_fatal_palette[2]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(overdose_fatal_palette[2], overdose_fatal_palette[3])))

# Output
save(ftable_overdose_fatal_qalys_MMS_out, 
     file = "outputs/DSA/Modified Model Specification/ftable_overdose_fatal_qalys_MMS_out.RData")
save(ftable_overdose_fatal_qalys_TS_out, 
     file = "outputs/DSA/Trial Specification/ftable_overdose_fatal_qalys_TS_out.RData")
save(ftable_overdose_fatal_costs_MMS_out, 
     file = "outputs/DSA/Modified Model Specification/ftable_overdose_fatal_costs_MMS_out.RData")
save(ftable_overdose_fatal_costs_TS_out, 
     file = "outputs/DSA/Trial Specification/ftable_overdose_fatal_costs_TS_out.RData")

## Non-fatal overdose ##
# Fentanyl prevalence
df_overdose_fent_exp_2020_low_MMS <- data.frame(ICER_fent_exp_2020_low_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fent_exp_2020_low_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                                ICER_fent_exp_2020_low_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_fent_exp_2020_low_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                                ICER_fent_exp_2020_low_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_fent_exp_2020_low_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                                ICER_fent_exp_2020_low_MMS$df_icer$n_icer_TOTAL_1yr, ICER_fent_exp_2020_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_fent_exp_2020_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_fent_exp_2020_low_TS <- data.frame(ICER_fent_exp_2020_low_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fent_exp_2020_low_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_fent_exp_2020_low_TS$df_icer$n_icer_TOTAL_1yr)

df_overdose_fent_exp_2020_high_MMS <- data.frame(ICER_fent_exp_2020_high_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fent_exp_2020_high_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                                 ICER_fent_exp_2020_high_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_fent_exp_2020_high_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                                 ICER_fent_exp_2020_high_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_fent_exp_2020_high_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                                 ICER_fent_exp_2020_high_MMS$df_icer$n_icer_TOTAL_1yr, ICER_fent_exp_2020_high_MMS$df_icer$n_icer_TOTAL_5yr, ICER_fent_exp_2020_high_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_fent_exp_2020_high_TS <- data.frame(ICER_fent_exp_2020_high_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fent_exp_2020_high_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_fent_exp_2020_high_TS$df_icer$n_icer_TOTAL_1yr)

# Reduction in fentanyl exposure for non-injection
df_overdose_ni_fent_reduction_low_MMS <- data.frame(ICER_ni_fent_reduction_low_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_ni_fent_reduction_low_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                                ICER_ni_fent_reduction_low_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_ni_fent_reduction_low_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                                ICER_ni_fent_reduction_low_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_ni_fent_reduction_low_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                                ICER_ni_fent_reduction_low_MMS$df_icer$n_icer_TOTAL_1yr, ICER_ni_fent_reduction_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_ni_fent_reduction_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_ni_fent_reduction_low_TS <- data.frame(ICER_ni_fent_reduction_low_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_ni_fent_reduction_low_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_ni_fent_reduction_low_TS$df_icer$n_icer_TOTAL_1yr)

df_overdose_ni_fent_reduction_high_MMS <- data.frame(ICER_ni_fent_reduction_high_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_ni_fent_reduction_high_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                                 ICER_ni_fent_reduction_high_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_ni_fent_reduction_high_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                                 ICER_ni_fent_reduction_high_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_ni_fent_reduction_high_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                                 ICER_ni_fent_reduction_high_MMS$df_icer$n_icer_TOTAL_1yr, ICER_ni_fent_reduction_high_MMS$df_icer$n_icer_TOTAL_5yr, ICER_ni_fent_reduction_high_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_ni_fent_reduction_high_TS <- data.frame(ICER_ni_fent_reduction_high_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_ni_fent_reduction_high_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_ni_fent_reduction_high_TS$df_icer$n_icer_TOTAL_1yr)

# Fentanyl overdose multiplier
df_overdose_fent_OD_mult_low_MMS <- data.frame(ICER_fent_OD_mult_low_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fent_OD_mult_low_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                                    ICER_fent_OD_mult_low_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_fent_OD_mult_low_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                                    ICER_fent_OD_mult_low_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_fent_OD_mult_low_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                                    ICER_fent_OD_mult_low_MMS$df_icer$n_icer_TOTAL_1yr, ICER_fent_OD_mult_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_fent_OD_mult_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_fent_OD_mult_low_TS <- data.frame(ICER_fent_OD_mult_low_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fent_OD_mult_low_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_fent_OD_mult_low_TS$df_icer$n_icer_TOTAL_1yr)

df_overdose_fent_OD_mult_high_MMS <- data.frame(ICER_fent_OD_mult_high_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fent_OD_mult_high_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                                     ICER_fent_OD_mult_high_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_fent_OD_mult_high_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, 
                                                     ICER_fent_OD_mult_high_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, ICER_fent_OD_mult_high_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, 
                                                     ICER_fent_OD_mult_high_MMS$df_icer$n_icer_TOTAL_1yr, ICER_fent_OD_mult_high_MMS$df_icer$n_icer_TOTAL_5yr, ICER_fent_OD_mult_high_MMS$df_icer$n_icer_TOTAL_10yr)
df_overdose_fent_OD_mult_high_TS <- data.frame(ICER_fent_OD_mult_high_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_fent_OD_mult_high_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_fent_OD_mult_high_TS$df_icer$n_icer_TOTAL_1yr)

colnames(df_overdose_baseline_MMS) <- colnames(df_overdose_fent_exp_2020_low_MMS) <- colnames(df_overdose_fent_exp_2020_high_MMS) <- colnames(df_overdose_ni_fent_reduction_low_MMS) <- colnames(df_overdose_ni_fent_reduction_high_MMS) <- colnames(df_overdose_fent_OD_mult_low_MMS) <- colnames(df_overdose_fent_OD_mult_high_MMS) <- c("inc_costs_1yr", "inc_costs_5yr", "inc_costs_10yr", "inc_qalys_1yr", 
                                                                                                                                                                                                                                                                                                                                           "inc_qalys_5yr", "inc_qalys_10yr", "icer_1yr", "icer_5yr", "icer_10yr")
colnames(df_overdose_baseline_TS) <- colnames(df_overdose_fent_exp_2020_low_TS) <- colnames(df_overdose_fent_exp_2020_high_TS) <- colnames(df_overdose_ni_fent_reduction_low_TS) <- colnames(df_overdose_ni_fent_reduction_high_TS) <- colnames(df_overdose_fent_OD_mult_low_TS) <- colnames(df_overdose_fent_OD_mult_high_TS) <- c("inc_costs_1yr", "inc_qalys_1yr", "icer_1yr")

df_overdose_nonfatal_MMS <- rbind(df_overdose_baseline_MMS, df_overdose_fent_exp_2020_low_MMS, df_overdose_fent_exp_2020_high_MMS, df_overdose_ni_fent_reduction_low_MMS, df_overdose_ni_fent_reduction_high_MMS, df_overdose_fent_OD_mult_low_MMS, df_overdose_fent_OD_mult_high_MMS)
df_overdose_nonfatal_TS <- rbind(df_overdose_baseline_TS, df_overdose_fent_exp_2020_low_TS, df_overdose_fent_exp_2020_high_TS, df_overdose_ni_fent_reduction_low_TS, df_overdose_ni_fent_reduction_high_TS, df_overdose_fent_OD_mult_low_TS, df_overdose_fent_OD_mult_high_TS)

df_overdose_nonfatal_MMS <- data.frame("Scenario" = c("Baseline", "Fentanyl Exp (low)", "Fentanyl Exp (high)", "Fent Reduction (NI) (low)", "Fent Reduction (NI) (high)", "Fent OD Mult (low)", "Fent OD Mult (high)"), df_overdose_nonfatal_MMS)
df_overdose_nonfatal_TS <- data.frame("Scenario" = c("Baseline", "Fentanyl Exp (low)", "Fentanyl Exp (high)", "Fent Reduction (NI) (low)", "Fent Reduction (NI) (high)", "Fent OD Mult (low)", "Fent OD Mult (high)"), df_overdose_nonfatal_TS)

# Custom table output
overdose_nonfatal_palette <- brewer.pal(3,"PuOr")

table_overdose_nonfatal_costs_MMS <- df_overdose_nonfatal_MMS %>%
  mutate(`Incremental Costs (1-year)` = accounting(inc_costs_1yr, 0),
         `Incremental Costs (5-year)` = accounting(inc_costs_5yr, 0),
         `Incremental Costs (10-year)` = accounting(inc_costs_10yr, 0),
         `Difference v. Baseline (1-year)` = `Incremental Costs (1-year)` - `Incremental Costs (1-year)`[`Scenario` == "Baseline"],
         `Difference v. Baseline (5-year)` = `Incremental Costs (5-year)` - `Incremental Costs (5-year)`[`Scenario` == "Baseline"],
         `Difference v. Baseline (10-year)` = `Incremental Costs (10-year)` - `Incremental Costs (10-year)`[`Scenario` == "Baseline"],
         `ICER (1-year)` = accounting(icer_1yr, 0),
         `ICER (5-year)` = accounting(icer_5yr, 0),
         `ICER (10-year)` = accounting(icer_10yr, 0)) %>%
  select(c(`Scenario`, `Incremental Costs (1-year)`, `Incremental Costs (5-year)`, `Incremental Costs (10-year)`, 
           `Difference v. Baseline (1-year)`, `Difference v. Baseline (5-year)`, `Difference v. Baseline (10-year)`,
           `ICER (1-year)`, `ICER (5-year)`, `ICER (10-year)`))

table_overdose_nonfatal_costs_TS <- df_overdose_nonfatal_TS %>%
  mutate(`Incremental Costs (1-year)` = accounting(inc_costs_1yr, 0),
         `Difference v. Baseline (1-year)` = `Incremental Costs (1-year)` - `Incremental Costs (1-year)`[`Scenario` == "Baseline"],
         `ICER (1-year)` = accounting(icer_1yr, 0)) %>%
  select(c(`Scenario`, `Incremental Costs (1-year)`, `Difference v. Baseline (1-year)`, `ICER (1-year)`))

table_overdose_nonfatal_qalys_MMS <- df_overdose_nonfatal_MMS %>%
  mutate(`Incremental QALYs (1-year)` = round(inc_qalys_1yr, 3),
         `Incremental QALYs (5-year)` = round(inc_qalys_5yr, 3),
         `Incremental QALYs (10-year)` = round(inc_qalys_10yr, 3),
         `Difference v. Baseline (1-year)` = round(`Incremental QALYs (1-year)` - `Incremental QALYs (1-year)`[`Scenario` == "Baseline"], 3),
         `Difference v. Baseline (5-year)` = round(`Incremental QALYs (5-year)` - `Incremental QALYs (5-year)`[`Scenario` == "Baseline"], 3),
         `Difference v. Baseline (10-year)` = round(`Incremental QALYs (10-year)` - `Incremental QALYs (10-year)`[`Scenario` == "Baseline"], 3),
         `ICER (1-year)` = accounting(icer_1yr, 0),
         `ICER (5-year)` = accounting(icer_5yr, 0),
         `ICER (10-year)` = accounting(icer_10yr, 0)) %>%
  select(c(`Scenario`, `Incremental QALYs (1-year)`, `Incremental QALYs (5-year)`, `Incremental QALYs (10-year)`, 
           `Difference v. Baseline (1-year)`, `Difference v. Baseline (5-year)`, `Difference v. Baseline (10-year)`,
           `ICER (1-year)`, `ICER (5-year)`, `ICER (10-year)`))

table_overdose_nonfatal_qalys_TS <- df_overdose_nonfatal_TS %>%
  mutate(`Incremental QALYs (1-year)` = round(inc_qalys_1yr, 3),
         `Difference v. Baseline (1-year)` = round(`Incremental QALYs (1-year)` - `Incremental QALYs (1-year)`[`Scenario` == "Baseline"], 3),
         `ICER (1-year)` = accounting(icer_1yr, 0)) %>%
  select(c(`Scenario`, `Incremental QALYs (1-year)`, `Difference v. Baseline (1-year)`, `ICER (1-year)`))

ftable_overdose_nonfatal_qalys_MMS_out <- formattable(table_overdose_nonfatal_qalys_MMS, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental QALYs (1-year)` = color_tile(overdose_nonfatal_palette[1], overdose_nonfatal_palette[2]),
  `Incremental QALYs (5-year)` = color_tile(overdose_nonfatal_palette[1], overdose_nonfatal_palette[2]),
  `Incremental QALYs (10-year)` = color_tile(overdose_nonfatal_palette[1], overdose_nonfatal_palette[2]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `Difference v. Baseline (5-year)` = improvement_formatter,
  `Difference v. Baseline (10-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(overdose_nonfatal_palette[2], overdose_nonfatal_palette[3]),
  `ICER (5-year)` = color_tile(overdose_nonfatal_palette[2], overdose_nonfatal_palette[3]),
  `ICER (10-year)` = color_tile(overdose_nonfatal_palette[2], overdose_nonfatal_palette[3])))

ftable_overdose_nonfatal_qalys_TS_out <- formattable(table_overdose_nonfatal_qalys_TS, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental QALYs (1-year)` = color_tile(overdose_nonfatal_palette[1], overdose_nonfatal_palette[2]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(overdose_nonfatal_palette[2], overdose_nonfatal_palette[3])))

ftable_overdose_nonfatal_costs_MMS_out <- formattable(table_overdose_nonfatal_costs_MMS, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental Costs (1-year)` = color_tile(overdose_nonfatal_palette[1], overdose_nonfatal_palette[2]),
  `Incremental Costs (5-year)` = color_tile(overdose_nonfatal_palette[1], overdose_nonfatal_palette[2]),
  `Incremental Costs (10-year)` = color_tile(overdose_nonfatal_palette[1], overdose_nonfatal_palette[2]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `Difference v. Baseline (5-year)` = improvement_formatter,
  `Difference v. Baseline (10-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(overdose_nonfatal_palette[2], overdose_nonfatal_palette[3]),
  `ICER (5-year)` = color_tile(overdose_nonfatal_palette[2], overdose_nonfatal_palette[3]),
  `ICER (10-year)` = color_tile(overdose_nonfatal_palette[2], overdose_nonfatal_palette[3])))

ftable_overdose_nonfatal_costs_TS_out <- formattable(table_overdose_nonfatal_costs_TS, align =c("l","c","c","c","c","c","c","c","c","c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental Costs (1-year)` = color_tile(overdose_nonfatal_palette[1], overdose_nonfatal_palette[2]),
  `Difference v. Baseline (1-year)` = improvement_formatter,
  `ICER (1-year)` = color_tile(overdose_nonfatal_palette[2], overdose_nonfatal_palette[3])))

# Output
save(ftable_overdose_nonfatal_qalys_MMS_out, 
     file = "outputs/DSA/Modified Model Specification/ftable_overdose_nonfatal_qalys_MMS_out.RData")
save(ftable_overdose_nonfatal_qalys_TS_out, 
     file = "outputs/DSA/Trial Specification/ftable_overdose_nonfatal_qalys_TS_out.RData")
save(ftable_overdose_nonfatal_costs_MMS_out, 
     file = "outputs/DSA/Modified Model Specification/ftable_overdose_nonfatal_costs_MMS_out.RData")
save(ftable_overdose_nonfatal_costs_TS_out, 
     file = "outputs/DSA/Trial Specification/ftable_overdose_nonfatal_costs_TS_out.RData")


# Prepare data for plotting
# Costs
#df_low_costs <- rbind(df_crime_costs_low_MMS$n_inc_costs_TOTAL_10yr, df_qalys_low_MMS$n_inc_costs_TOTAL_10yr)
#df_high_costs <- rbind(df_crime_costs_high_MMS$n_inc_costs_TOTAL_10yr, df_qalys_high_MMS$n_inc_costs_TOTAL_10yr)
#df_tornado_costs <- data.frame(df_low, df_high)

#colnames(df_tornado_costs) <- c("Low", "High")
#row.names(df_tornado_costs) <- c("Crime Costs", "QALYs")

# ICER
#df_low_icer <- rbind(df_crime_costs_low_MMS$n_icer_TOTAL_10yr, df_qalys_low_MMS$n_icer_TOTAL_10yr)
#df_high_icer <- rbind(df_crime_costs_high_MMS$n_icer_TOTAL_10yr, df_qalys_high_MMS$n_icer_TOTAL_10yr)
#df_tornado_icer <- data.frame(df_low, df_high)

#colnames(df_tornado_icer) <- c("Low", "High")
#row.names(df_tornado_icer) <- c("Crime Costs", "QALYs")

#########################
#### Tornado Diagram ####
#########################
# this is throwing some warnings in my computer, but it is reading the data frame correctly
df <- '
Parameter Lower_Bound Upper_Bound UL_Difference
Parameter01 8074 11181 3108 
Parameter02 8177 11007 2831 
Parameter03 8879 10188 1308 
Parameter04 4358 18697 14339 
Parameter05 9073 10087 1013 
Parameter06 12034 7572 4462 
Parameter07 11357 7933 3423 
Parameter08 9769 9202 567 
Parameter09 8833 10403 1570 
Parameter10 13450 4219 9231 
Parameter11 10691 7915 2776 
Parameter12 10036 8792 1244
' %>% read_table2()

# original value of output
baseline_icer <- df_baseline_MMS$n_icer_TOTAL_10yr
baseline_costs <- df_baseline_MMS$n_inc_costs_TOTAL_10yr
baseline_qalys <- df_baseline_MMS$n_inc_qalys_TOTAL_10yr

# get order of parameters according to size of intervals
# (I use this to define the ordering of the factors which I then use to define the positions in the plot)
order.parameters <- df %>% arrange(UL_Difference) %>%
  mutate(Parameter=factor(x=Parameter, levels=Parameter)) %>%
  select(Parameter) %>% unlist() %>% levels()

# width of columns in plot (value between 0 and 1)
width <- 0.95

# get data frame in shape for ggplot and geom_rect
df.2 <- df %>% 
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key='type', value='output.value', Lower_Bound:Upper_Bound) %>%
  # just reordering columns
  select(Parameter, type, output.value, UL_Difference) %>%
  # create the columns for geom_rect
  mutate(Parameter=factor(Parameter, levels=order.parameters),
         ymin=pmin(output.value, base.value),
         ymax=pmax(output.value, base.value),
         xmin=as.numeric(Parameter)-width/2,
         xmax=as.numeric(Parameter)+width/2)

# create plot
# (use scale_x_continuous to change labels in y axis to name of parameters)
png(width = 960, height = 540)
ggplot() + 
  geom_rect(data = df.2, 
            aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, fill=type)) +
  theme_bw() + 
  theme(axis.title.y=element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  geom_hline(yintercept = base.value) +
  scale_x_continuous(breaks = c(1:length(order.parameters)), 
                     labels = order.parameters) +
  coord_flip()
dev.off()

## Trial health state versus model health state definition ##

## Province-specific analysis ##
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
#df_dsa_crime_costs_TS <- read.csv(file = "data/DSA/Trial Specification/crime_costs.csv", row.names = 1, header = TRUE)
# MMS
v_dsa_crime_costs_low_MMS <- unlist(df_dsa_crime_costs_MMS["pe_low",])
v_dsa_crime_costs_high_MMS <- unlist(df_dsa_crime_costs_MMS["pe_high",])
v_dsa_crime_costs_alt_MMS <- unlist(df_dsa_crime_costs_MMS["pe_alt",])
v_dsa_crime_costs_reduced_MMS <- unlist(df_dsa_crime_costs_MMS["pe_reduced_trial",])

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

# High
# MMS
l_outcomes_MET_crime_costs_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_high_MMS)
l_outcomes_BUP_crime_costs_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_high_MMS)
ICER_crime_costs_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_high_MMS, outcomes_int = l_outcomes_BUP_crime_costs_high_MMS)

# Alternative (OPTIMA trial estimates)
# MMS
l_outcomes_MET_crime_costs_alt_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_MMS)
l_outcomes_BUP_crime_costs_alt_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_MMS)
ICER_crime_costs_alt_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_alt_MMS, outcomes_int = l_outcomes_BUP_crime_costs_alt_MMS)

# Trial estimates (equal across treatments, REL = concurrent)
# MMS
l_outcomes_MET_crime_costs_reduced_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_reduced_MMS)
l_outcomes_BUP_crime_costs_reduced_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_reduced_MMS)
ICER_crime_costs_reduced_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_reduced_MMS, outcomes_int = l_outcomes_BUP_crime_costs_reduced_MMS)

################
### Baseline ###
################
df_baseline_MMS <- data.frame(ICER_MMS$df_incremental, ICER_MMS$df_icer)
df_baseline_TS  <- data.frame(ICER_TS$df_incremental, ICER_TS$df_icer)

################
### Overdose ###
################
# Baseline
df_overdose_baseline_MMS <- data.frame(ICER_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_MMS$df_incremental$n_inc_qalys_TOTAL_life, 
                                       ICER_MMS$df_icer$n_icer_TOTAL_life)
#df_overdose_baseline_TS <- data.frame(ICER_TS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_TS$df_icer$n_icer_TOTAL_1yr)

# Costs
# Non-fatal OD
v_overdose_fent_exp_costs_MMS <- c(ICER_fent_exp_2020_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_fent_exp_2020_high_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_overdose_ni_fent_reduction_costs_MMS <- c(ICER_ni_fent_reduction_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_ni_fent_reduction_high_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_overdose_fent_OD_mult_costs_MMS <- c(ICER_fent_OD_mult_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_fent_OD_mult_high_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_overdose_BUP_OD_mult_costs_MMS <- c(ICER_BUP_OD_mult_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_BUP_OD_mult_high_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_overdose_MET_OD_mult_costs_MMS <- c(ICER_MET_OD_mult_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_MET_OD_mult_high_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_overdose_REL_OD_mult_costs_MMS <- c(ICER_REL_OD_mult_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_REL_OD_mult_high_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_overdose_INJ_OD_mult_costs_MMS <- c(ICER_INJ_OD_mult_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_INJ_OD_mult_high_MMS$df_incremental$n_inc_costs_TOTAL_life)

# Fatal OD
v_overdose_witness_costs_MMS <- c(ICER_witness_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_witness_high_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_overdose_NX_prob_costs_MMS <- c(ICER_NX_prob_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_NX_prob_high_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_overdose_NX_success_costs_MMS <- c(ICER_NX_success_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_NX_success_high_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_overdose_fatal_OD_costs_MMS <- c(ICER_fatal_OD_low_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_fatal_OD_high_MMS$df_incremental$n_inc_costs_TOTAL_life)

m_overdose_costs_MMS <- rbind(v_overdose_fent_exp_costs_MMS, v_overdose_ni_fent_reduction_costs_MMS, v_overdose_fent_OD_mult_costs_MMS,
                              v_overdose_BUP_OD_mult_costs_MMS, v_overdose_MET_OD_mult_costs_MMS, v_overdose_REL_OD_mult_costs_MMS,
                              v_overdose_INJ_OD_mult_costs_MMS, v_overdose_witness_costs_MMS, v_overdose_NX_prob_costs_MMS, 
                              v_overdose_NX_success_costs_MMS, v_overdose_fatal_OD_costs_MMS)

df_overdose_costs_MMS <- as.data.frame(m_overdose_costs_MMS)
colnames(df_overdose_costs_MMS) <- c("Lower", "Upper")
df_overdose_costs_MMS <- as_data_frame(df_overdose_costs_MMS) %>% mutate(diff = abs(Upper - Lower),
                                                                         base = ICER_MMS$df_incremental$n_inc_costs_TOTAL_life) %>%
  add_column(var_name = c("Fentanyl prevalence", "Fentanyl reduction (non-injection)", "Fentanyl OD mult", "BUP OD first month mult", "MET OD first month mult", 
                          "REL OD first month mult", "INJ OD mult", "Probability witnessed", "Probability NX", "NX success", "Fatal OD rate"))

# Save outputs
## As .RData ##
save(df_overdose_costs_MMS, 
     file = "outputs/DSA/Modified Model Specification/df_overdose_costs_MMS.RData")

# QALYs
# Non-fatal OD
v_overdose_fent_exp_qalys_MMS <- c(ICER_fent_exp_2020_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_fent_exp_2020_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_overdose_ni_fent_reduction_qalys_MMS <- c(ICER_ni_fent_reduction_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_ni_fent_reduction_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_overdose_fent_OD_mult_qalys_MMS <- c(ICER_fent_OD_mult_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_fent_OD_mult_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_overdose_BUP_OD_mult_qalys_MMS <- c(ICER_BUP_OD_mult_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_BUP_OD_mult_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_overdose_MET_OD_mult_qalys_MMS <- c(ICER_MET_OD_mult_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_MET_OD_mult_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_overdose_REL_OD_mult_qalys_MMS <- c(ICER_REL_OD_mult_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_REL_OD_mult_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_overdose_INJ_OD_mult_qalys_MMS <- c(ICER_INJ_OD_mult_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_INJ_OD_mult_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)

# Fatal OD
v_overdose_witness_qalys_MMS <- c(ICER_witness_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_witness_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_overdose_NX_prob_qalys_MMS <- c(ICER_NX_prob_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_NX_prob_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_overdose_NX_success_qalys_MMS <- c(ICER_NX_success_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_NX_success_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_overdose_fatal_OD_qalys_MMS <- c(ICER_fatal_OD_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_fatal_OD_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)

m_overdose_qalys_MMS <- rbind(v_overdose_fent_exp_qalys_MMS, v_overdose_ni_fent_reduction_qalys_MMS, v_overdose_fent_OD_mult_qalys_MMS,
                              v_overdose_BUP_OD_mult_qalys_MMS, v_overdose_MET_OD_mult_qalys_MMS, v_overdose_REL_OD_mult_qalys_MMS,
                              v_overdose_INJ_OD_mult_qalys_MMS, v_overdose_witness_qalys_MMS, v_overdose_NX_prob_qalys_MMS, 
                              v_overdose_NX_success_qalys_MMS, v_overdose_fatal_OD_qalys_MMS)

df_overdose_qalys_MMS <- as.data.frame(m_overdose_qalys_MMS)
colnames(df_overdose_qalys_MMS) <- c("Lower", "Upper")
df_overdose_qalys_MMS <- as_data_frame(df_overdose_qalys_MMS) %>% mutate(diff = abs(Upper - Lower),
                                                                         base = ICER_MMS$df_incremental$n_inc_qalys_TOTAL_life) %>%
  add_column(var_name = c("Fentanyl prevalence", "Fentanyl reduction (non-injection)", "Fentanyl OD mult", "BUP OD first month mult", "MET OD first month mult", 
                          "REL OD first month mult", "INJ OD mult", "Probability witnessed", "Probability NX", "NX success", "Fatal OD rate"))

# Save outputs
## As .RData ##
save(df_overdose_qalys_MMS, 
     file = "outputs/DSA/Modified Model Specification/df_overdose_qalys_MMS.RData")

# ICER
# v_overdose_witness_icer_MMS <- c(ICER_witness_low_MMS$df_incremental$n_inc_icer_TOTAL_life, ICER_witness_high_MMS$df_incremental$n_inc_icer_TOTAL_life)
# v_overdose_NX_prob_icer_MMS <- c(ICER_NX_prob_low_MMS$df_incremental$n_inc_icer_TOTAL_life, ICER_NX_prob_high_MMS$df_incremental$n_inc_icer_TOTAL_life)
# v_overdose_NX_success_icer_MMS <- c(ICER_NX_success_low_MMS$df_incremental$n_inc_icer_TOTAL_life, ICER_NX_success_high_MMS$df_incremental$n_inc_icer_TOTAL_life)
# v_overdose_fatal_OD_icer_MMS <- c(ICER_fatal_OD_low_MMS$df_incremental$n_inc_icer_TOTAL_life, ICER_fatal_OD_high_MMS$df_incremental$n_inc_icer_TOTAL_life)


#########################
#### Tornado Diagram ####
#########################
# Costs
v_order_parameters <- df_overdose_costs_MMS %>% arrange(diff) %>%
  mutate(var_name = factor(x = var_name, levels = var_name)) %>%
  select(var_name) %>% unlist() %>% levels()

# width of columns in plot (value between 0 and 1)
width <- 0.75
# get data frame in shape for ggplot and geom_rect
df.2 <- df_overdose_costs_MMS %>% 
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key = 'type', value = 'output.value', Lower:Upper) %>%
  # just reordering columns
  select(var_name, type, output.value, diff, base) %>%
  # create the columns for geom_rect
  mutate(var_name = factor(var_name, levels = v_order_parameters),
         ymin = pmin(output.value, base),
         ymax = pmax(output.value, base),
         xmin = as.numeric(var_name) - width/2,
         xmax = as.numeric(var_name) + width/2)

# create plot
# (use scale_x_continuous to change labels in y axis to name of parameters)
p_tornado_overdose_costs <- ggplot() + 
  geom_rect(data = df.2, 
            aes(ymax = ymax, ymin = ymin, xmax = xmax, xmin = xmin, fill = type)) +
  theme_bw() + 
  scale_fill_manual(values = c("Upper" = "#5ab4ac",
                               "Lower" = "#d8b365")) +
  theme(axis.title.y = element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  geom_hline(yintercept = df.2$base) +
  scale_x_continuous(breaks = c(1:length(v_order_parameters)), 
                     labels = v_order_parameters) +
  xlab("Parameter") + ylab("Incremental Cost") +
  coord_flip()

png(file = "Plots/DSA/Modified Model Spec/tornado_overdose_costs.png", width = 600, height = 600)
p_tornado_overdose_costs
dev.off()

# QALYs
v_order_parameters <- df_overdose_qalys_MMS %>% arrange(diff) %>%
  mutate(var_name = factor(x = var_name, levels = var_name)) %>%
  select(var_name) %>% unlist() %>% levels()

# width of columns in plot (value between 0 and 1)
width <- 0.75
# get data frame in shape for ggplot and geom_rect
df.2 <- df_overdose_qalys_MMS %>% 
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key = 'type', value = 'output.value', Lower:Upper) %>%
  # just reordering columns
  select(var_name, type, output.value, diff, base) %>%
  # create the columns for geom_rect
  mutate(var_name = factor(var_name, levels = v_order_parameters),
         ymin = pmin(output.value, base),
         ymax = pmax(output.value, base),
         xmin = as.numeric(var_name) - width/2,
         xmax = as.numeric(var_name) + width/2)

# create plot
# (use scale_x_continuous to change labels in y axis to name of parameters)
#png(file = "Plots/DSA/Modified Model Spec/tornado_overdose_qalys.png", width = 960, height = 540)
p_tornado_overdose_qalys <- ggplot() + 
  geom_rect(data = df.2, 
            aes(ymax = ymax, ymin = ymin, xmax = xmax, xmin = xmin, fill = type)) +
  theme_bw() + 
  scale_fill_manual(values = c("Upper" = "#5ab4ac",
                               "Lower" = "#d8b365")) +
  theme(axis.title.y=element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  geom_hline(yintercept = df.2$base) +
  scale_x_continuous(breaks = c(1:length(v_order_parameters)), 
                     labels = v_order_parameters) +
  xlab("Parameter") + ylab("Incremental Cost") +
  coord_flip()
#dev.off()

png(file = "Plots/DSA/Modified Model Spec/tornado_overdose_qalys.png", width = 600, height = 600)
p_tornado_overdose_qalys
dev.off()

#p_tornado_overdose_costs
#p_tornado_overdose_qalys

## Trial health state versus model health state definition ##

## Province-specific analysis ##
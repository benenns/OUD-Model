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
# v_dsa_qalys_reduced_eq_5d_5l_TS <- unlist(df_dsa_qalys_TS["pe_reduced_eq_5d_5l",])
# v_dsa_qalys_low_TS <- unlist(df_dsa_qalys_TS["pe_low",])
# v_dsa_qalys_high_TS <- unlist(df_dsa_qalys_TS["pe_high",])
# v_dsa_qalys_eq_5d_3l_TS <- unlist(df_dsa_qalys_TS["pe_eq_5d_3l",])
# v_dsa_qalys_hui_3_TS <- unlist(df_dsa_qalys_TS["pe_hui_3",])
# v_dsa_qalys_odn_low_TS <- unlist(df_dsa_qalys_TS["pe_odn_low",])

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
#l_outcomes_MET_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map)
#l_outcomes_BUP_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map)
#ICER_TS <- ICER(outcomes_comp = l_outcomes_MET_TS, outcomes_int = l_outcomes_BUP_TS)

#############
### QALYs ###
#############
# Low
# MMS
l_outcomes_MET_qalys_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_low_MMS)
l_outcomes_BUP_qalys_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_low_MMS)
ICER_qalys_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_low_MMS, outcomes_int = l_outcomes_BUP_qalys_low_MMS)
# TS
# l_outcomes_MET_qalys_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_low_TS)
# l_outcomes_BUP_qalys_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_low_TS)
# ICER_qalys_low_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_low_TS, outcomes_int = l_outcomes_BUP_qalys_low_TS)

# High
# MMS
l_outcomes_MET_qalys_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_high_MMS)
l_outcomes_BUP_qalys_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_high_MMS)
ICER_qalys_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_high_MMS, outcomes_int = l_outcomes_BUP_qalys_high_MMS)
# TS
# l_outcomes_MET_qalys_high_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_high_TS)
# l_outcomes_BUP_qalys_high_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_high_TS)
# ICER_qalys_high_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_high_TS, outcomes_int = l_outcomes_BUP_qalys_high_TS)

## Reduced (EQ-5D-5L) ##
# MMS
l_outcomes_MET_qalys_reduced_eq_5d_5l_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_reduced_eq_5d_5l_MMS)
l_outcomes_BUP_qalys_reduced_eq_5d_5l_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_reduced_eq_5d_5l_MMS)
ICER_qalys_reduced_eq_5d_5l_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_reduced_eq_5d_5l_MMS, outcomes_int = l_outcomes_BUP_qalys_reduced_eq_5d_5l_MMS)
# TS
# l_outcomes_MET_qalys_reduced_eq_5d_5l_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_reduced_eq_5d_5l_TS)
# l_outcomes_BUP_qalys_reduced_eq_5d_5l_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_reduced_eq_5d_5l_TS)
# ICER_qalys_reduced_eq_5d_5l_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_reduced_eq_5d_5l_TS, outcomes_int = l_outcomes_BUP_qalys_reduced_eq_5d_5l_TS)

## Alternative (EQ-5D-3L) ##
# MMS
l_outcomes_MET_qalys_eq_5d_3l_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_eq_5d_3l_MMS)
l_outcomes_BUP_qalys_eq_5d_3l_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_eq_5d_3l_MMS)
ICER_qalys_eq_5d_3l_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_eq_5d_3l_MMS, outcomes_int = l_outcomes_BUP_qalys_eq_5d_3l_MMS)
# TS
# l_outcomes_MET_qalys_eq_5d_3l_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_eq_5d_3l_TS)
# l_outcomes_BUP_qalys_eq_5d_3l_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_eq_5d_3l_TS)
# ICER_qalys_eq_5d_3l_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_eq_5d_3l_TS, outcomes_int = l_outcomes_BUP_qalys_eq_5d_3l_TS)

## Alternative (HUI-3) ##
# MMS
l_outcomes_MET_qalys_hui_3_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_hui_3_MMS)
l_outcomes_BUP_qalys_hui_3_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_hui_3_MMS)
ICER_qalys_hui_3_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_hui_3_MMS, outcomes_int = l_outcomes_BUP_qalys_hui_3_MMS)
# TS
# l_outcomes_MET_qalys_hui_3_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_hui_3_TS)
# l_outcomes_BUP_qalys_hui_3_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_hui_3_TS)
# ICER_qalys_hui_3_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_hui_3_TS, outcomes_int = l_outcomes_BUP_qalys_hui_3_TS)

## Alternative overdose ##
# MMS
l_outcomes_MET_qalys_odn_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_odn_low_MMS)
l_outcomes_BUP_qalys_odn_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_odn_low_MMS)
ICER_qalys_odn_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_odn_low_MMS, outcomes_int = l_outcomes_BUP_qalys_odn_low_MMS)
# TS
# l_outcomes_MET_qalys_odn_low_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_odn_low_TS)
# l_outcomes_BUP_qalys_odn_low_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_odn_low_TS)
# ICER_qalys_odn_low_TS <- ICER(outcomes_comp = l_outcomes_MET_qalys_odn_low_TS, outcomes_int = l_outcomes_BUP_qalys_odn_low_TS)

################
### Baseline ###
################
df_baseline_MMS <- data.frame(ICER_MMS$df_incremental, ICER_MMS$df_icer)
#df_baseline_TS  <- data.frame(ICER_TS$df_incremental, ICER_TS$df_icer)

#############
### QALYs ###
#############
# Baseline
v_qalys_baseline_MMS <- c(ICER_MMS$df_incremental$n_inc_qalys_TOTAL_life)
#df_qalys_baseline_TS <- data.frame(ICER_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_TS$df_icer$n_icer_TOTAL_1yr)
# Combined Tx
v_qalys_reduced_eq_5d_5l_MMS <- c(ICER_qalys_reduced_eq_5d_5l_MMS$df_incremental$n_inc_qalys_TOTAL_life)
#df_qalys_reduced_eq_5d_5l_TS <- data.frame(ICER_qalys_reduced_eq_5d_5l_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_reduced_eq_5d_5l_TS$df_icer$n_icer_TOTAL_1yr)
# Low
v_qalys_range_MMS <- c(ICER_qalys_low_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_qalys_high_MMS$df_incremental$n_inc_qalys_TOTAL_life)
#df_qalys_low_TS <- data.frame(ICER_qalys_low_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_low_TS$df_icer$n_icer_TOTAL_1yr)
# High
#v_qalys_high_MMS <- c(ICER_qalys_high_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_qalys_high_MMS$df_icer$n_icer_TOTAL_life)
#df_qalys_high_TS <- data.frame(ICER_qalys_high_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_high_TS$df_icer$n_icer_TOTAL_1yr)
# Alt (EQ-5d-3L)
v_qalys_eq_5d_3l_MMS <- c(ICER_qalys_eq_5d_3l_MMS$df_incremental$n_inc_qalys_TOTAL_life)
#df_qalys_eq_5d_3l_TS <- data.frame(ICER_qalys_eq_5d_3l_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_eq_5d_3l_TS$df_icer$n_icer_TOTAL_1yr)
# Alt (HUI-3)
v_qalys_hui_3_MMS <- c(ICER_qalys_hui_3_MMS$df_incremental$n_inc_qalys_TOTAL_life)
#df_qalys_hui_3_TS <- data.frame(ICER_qalys_hui_3_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_hui_3_TS$df_icer$n_icer_TOTAL_1yr)
# ODN low
v_qalys_odn_low_MMS <- c(ICER_qalys_odn_low_MMS$df_incremental$n_inc_qalys_TOTAL_life)
#df_qalys_odn_low_TS <- data.frame(ICER_qalys_odn_low_TS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_odn_low_TS$df_icer$n_icer_TOTAL_1yr)

m_qalys_MMS <- rbind(v_qalys_reduced_eq_5d_5l_MMS, v_qalys_range_MMS, v_qalys_eq_5d_3l_MMS, v_qalys_hui_3_MMS, v_qalys_odn_low_MMS)

df_qalys_MMS <- as.data.frame(m_qalys_MMS)
colnames(df_qalys_MMS) <- c("Lower", "Upper")
df_qalys_MMS <- as_data_frame(df_qalys_MMS) %>% mutate(diff = abs(Upper - Lower),
                                                       base = ICER_MMS$df_incremental$n_inc_qalys_TOTAL_life) %>%
  add_column(var_name = c("Combined Tx", "Range", "EQ-5D-3L", "HUI-3", "ODN (low)"))

## As .RData ##
save(df_qalys_MMS, 
     file = "outputs/DSA/Modified Model Specification/df_qalys_MMS.RData")

#########################
#### Tornado Diagram ####
#########################
# QALYs
v_order_parameters <- df_qalys_MMS %>% arrange(diff) %>%
  mutate(var_name = factor(x = var_name, levels = var_name)) %>%
  select(var_name) %>% unlist() %>% levels()

# width of columns in plot (value between 0 and 1)
width <- 0.75
# get data frame in shape for ggplot and geom_rect
df.2 <- df_qalys_MMS %>% 
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
p_tornado_qalys <- ggplot() + 
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

png(file = "Plots/DSA/Modified Model Spec/tornado_qalys.png", width = 600, height = 600)
p_tornado_qalys
dev.off()

#p_tornado_overdose_costs
#p_tornado_overdose_qalys

## Trial health state versus model health state definition ##

## Province-specific analysis ##
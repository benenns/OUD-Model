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
#source("R/generate_psa_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/ICER_functions.R")

# Load parameters
source("Analysis/00_load_parameters.R")

# Load DSA parameters


# Overdose

## Costs
# HRU costs
df_dsa_HRU_costs_MMS <- read.csv(file = "data/DSA/Modified Model Specification/HRU_costs.csv", row.names = 1, header = TRUE)
df_dsa_HRU_costs_TS <- read.csv(file = "data/DSA/Trial Specification/HRU_costs.csv", row.names = 1, header = TRUE)
# MMS
v_dsa_HRU_costs_alt_MMS <- unlist(df_dsa_HRU_costs_MMS["pe_alt",]) #as.numeric(df_dsa_costs["pe_alt",])
# TS
v_dsa_HRU_costs_alt_TS <- unlist(df_dsa_HRU_costs_TS["pe_alt",])

# Crime costs
# DSA data
df_dsa_crime_costs_MMS <- read.csv(file = "data/DSA/Modified Model Specification/crime_costs.csv", row.names = 1, header = TRUE)
df_dsa_crime_costs_TS <- read.csv(file = "data/DSA/Trial Specification/crime_costs.csv", row.names = 1, header = TRUE)
# MMS
v_dsa_crime_costs_low_MMS <- unlist(df_dsa_crime_costs_MMS["pe_low",])
v_dsa_crime_costs_high_MMS <- unlist(df_dsa_crime_costs_MMS["pe_high",])
v_dsa_crime_costs_alt_MMS <- unlist(df_dsa_crime_costs_MMS["pe_alt",])
# TS
v_dsa_crime_costs_alt_TS <- unlist(df_dsa_crime_costs_TS["pe_alt",])

# QALYs
df_dsa_qalys_MMS <- read.csv(file = "data/DSA/Modified Model Specification/qalys.csv", row.names = 1, header = TRUE)
#df_dsa_qalys_TS <- read.csv(file = "data/DSA/Trial Specification/qalys.csv", row.names = 1, header = TRUE)
# MMS
v_dsa_qalys_reduced_eq_5d_5l_MMS <- unlist(df_dsa_qalys_MMS["pe_reduced_eq_5d_5l",])
v_dsa_qalys_low_MMS <- unlist(df_dsa_qalys_MMS["pe_low",])
v_dsa_qalys_high_MMS <- unlist(df_dsa_qalys_MMS["pe_high",])
v_dsa_qalys_eq_5d_3l_MMS <- unlist(df_dsa_qalys_MMS["pe_eq_5d_3l",])
v_dsa_qalys_hui_3_MMS <- unlist(df_dsa_qalys_MMS["pe_hui_3",])
# TS

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

###################
### Crime Costs ###
###################
# Low
# MMS
l_outcomes_MET_crime_costs_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_low_MMS)
l_outcomes_BUP_crime_costs_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_low_MMS)
ICER_crime_costs_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_low_MMS, outcomes_int = l_outcomes_BUP_crime_costs_low_MMS)
# TS

# High
# MMS
l_outcomes_MET_crime_costs_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_high_MMS)
l_outcomes_BUP_crime_costs_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_high_MMS)
ICER_crime_costs_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_high_MMS, outcomes_int = l_outcomes_BUP_crime_costs_high_MMS)
# TS

# Alternative (Krebs et al. 2014)
# MMS
l_outcomes_MET_crime_costs_alt_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_MMS)
l_outcomes_BUP_crime_costs_alt_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_MMS)
ICER_crime_costs_alt_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_alt_MMS, outcomes_int = l_outcomes_BUP_crime_costs_alt_MMS)
# TS
l_outcomes_MET_crime_costs_alt_TS <- outcomes(l_params_all = l_params_MET_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_TS)
l_outcomes_BUP_crime_costs_alt_TS <- outcomes(l_params_all = l_params_BUP_TS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_TS)
ICER_crime_costs_alt_TS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_alt_TS, outcomes_int = l_outcomes_BUP_crime_costs_alt_TS)

# Reduced (all treatment costs equal)

#############
### QALYs ###
#############
# Low
# MMS
l_outcomes_MET_qalys_low_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_low_MMS)
l_outcomes_BUP_qalys_low_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_low_MMS)
ICER_qalys_low_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_low_MMS, outcomes_int = l_outcomes_BUP_qalys_low_MMS)
# TS

# High
# MMS
l_outcomes_MET_qalys_high_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_high_MMS)
l_outcomes_BUP_qalys_high_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_high_MMS)
ICER_qalys_high_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_high_MMS, outcomes_int = l_outcomes_BUP_qalys_high_MMS)
# TS

## Reduced (EQ-5D-5L) ##
# MMS
l_outcomes_MET_qalys_reduced_eq_5d_5l_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_reduced_eq_5d_5l_MMS)
l_outcomes_BUP_qalys_reduced_eq_5d_5l_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_reduced_eq_5d_5l_MMS)
ICER_qalys_reduced_eq_5d_5l_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_reduced_eq_5d_5l_MMS, outcomes_int = l_outcomes_BUP_qalys_reduced_eq_5d_5l_MMS)

## Alternative (EQ-5D-3L) ##
# MMS
l_outcomes_MET_qalys_eq_5d_3l_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_eq_5d_3l_MMS)
l_outcomes_BUP_qalys_eq_5d_3l_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_eq_5d_3l_MMS)
ICER_qalys_eq_5d_3l_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_eq_5d_3l_MMS, outcomes_int = l_outcomes_BUP_qalys_eq_5d_3l_MMS)

## Alternative (HUI-3) ##
# MMS
l_outcomes_MET_qalys_hui_3_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_hui_3_MMS)
l_outcomes_BUP_qalys_hui_3_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_qalys_hui_3_MMS)
ICER_qalys_hui_3_MMS <- ICER(outcomes_comp = l_outcomes_MET_qalys_hui_3_MMS, outcomes_int = l_outcomes_BUP_qalys_hui_3_MMS)

################
### Overdose ###
################
# Witnessed OD
# Low
l_outcomes_MET_witness_low <- outcomes(l_params_all = l_params_MET, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_low)
l_outcomes_BUP_witness_low <- outcomes(l_params_all = l_params_BUP, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_low)
ICER_witness_low <- ICER(outcomes_comp = l_outcomes_MET_witness_low, outcomes_int = l_outcomes_BUP_witness_low)
# High
l_outcomes_MET_witness_high <- outcomes(l_params_all = l_params_MET, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_high)
l_outcomes_BUP_witness_high <- outcomes(l_params_all = l_params_BUP, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_high)
ICER_witness_high <- ICER(outcomes_comp = l_outcomes_MET_witness_high, outcomes_int = l_outcomes_BUP_witness_high)

# Fentanyl prevalence
# Low
l_outcomes_MET_fentanyl_low <- outcomes(l_params_all = l_params_MET, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fentanyl_low)
l_outcomes_BUP_fentanyl_low <- outcomes(l_params_all = l_params_BUP, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fentanyl_low)
ICER_fentanyl_low <- ICER(outcomes_comp = l_outcomes_MET_fentanyl_low, outcomes_int = l_outcomes_BUP_fentanyl_low)

# Naloxone prevalence
l_outcomes_MET_naloxone <- outcomes(l_params_all = l_params_MET, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_naloxone)
l_outcomes_BUP_naloxone <- outcomes(l_params_all = l_params_BUP, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_naloxone)
ICER_naloxone <- ICER(outcomes_comp = l_outcomes_MET_naloxone, outcomes_int = l_outcomes_BUP_naloxone)

## Cohort characteristics ##
# Starting age

# Male %

## Outcomes ##
# Costs



# QALYs

###################
#### Data Prep ####
###################

# Baseline
df_baseline_MMS <- data.frame(ICER_MMS$df_incremental, ICER_MMS$df_icer)

###################
### Crime costs ###
###################
df_crime_costs_baseline_MMS <- data.frame(ICER_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                          ICER_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_MMS$df_icer$n_icer_TOTAL_1yr, 
                                          ICER_MMS$df_icer$n_icer_TOTAL_5yr, ICER_MMS$df_icer$n_icer_TOTAL_10yr)
df_crime_costs_reduced_MMS <- data.frame()
df_crime_costs_low_MMS <- data.frame(ICER_crime_costs_low_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_low_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                     ICER_crime_costs_low_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_crime_costs_low_MMS$df_icer$n_icer_TOTAL_1yr, 
                                     ICER_crime_costs_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_crime_costs_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_crime_costs_high_MMS <- data.frame(ICER_crime_costs_high_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_high_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                      ICER_crime_costs_high_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_crime_costs_high_MMS$df_icer$n_icer_TOTAL_1yr, 
                                      ICER_crime_costs_high_MMS$df_icer$n_icer_TOTAL_5yr, ICER_crime_costs_high_MMS$df_icer$n_icer_TOTAL_10yr)
df_crime_costs_alt_MMS <- data.frame(ICER_crime_costs_alt_MMS$df_incremental$n_inc_costs_TOTAL_1yr, ICER_crime_costs_alt_MMS$df_incremental$n_inc_costs_TOTAL_5yr, 
                                     ICER_crime_costs_alt_MMS$df_incremental$n_inc_costs_TOTAL_10yr, ICER_crime_costs_alt_MMS$df_icer$n_icer_TOTAL_1yr, 
                                     ICER_crime_costs_alt_MMS$df_icer$n_icer_TOTAL_5yr, ICER_crime_costs_alt_MMS$df_icer$n_icer_TOTAL_10yr)

colnames(df_crime_costs_baseline_MMS) <- colnames(df_crime_costs_low_MMS) <- colnames(df_crime_costs_high_MMS) <- colnames(df_crime_costs_alt_MMS) <- c("inc_costs_1yr", "inc_costs_5yr", "inc_costs_10yr", "icer_1yr", "icer_5yr", "icer_10yr")

df_crime_costs <- rbind(df_crime_costs_baseline_MMS, df_crime_costs_low_MMS, df_crime_costs_high_MMS, df_crime_costs_alt_MMS)
df_crime_costs <- data.frame("Scenario" = c("Baseline", "Low", "High", "Alternative"), df_crime_costs)

# Custom table output
# set colours
crime_costs_palette <- brewer.pal(3,"BrBG")

# table
table_crime_costs <- df_crime_costs %>%
  mutate(`Incremental Costs (1-year)` = round(inc_costs_1yr, 3),
         `Incremental Costs (5-year)` = round(inc_costs_5yr, 3),
         `Incremental Costs (10-year)` = round(inc_costs_10yr, 3),
         `ICER (1-year)` = accounting(icer_1yr, 0),
         `ICER (5-year)` = accounting(icer_5yr, 0),
         `ICER (10-year)` = accounting(icer_10yr, 0)) %>%
  select(c(`Scenario`, `Incremental Costs (1-year)`, `Incremental Costs (5-year)`, `Incremental Costs (10-year)`, `ICER (1-year)`, `ICER (5-year)`, `ICER (10-year)`)) #%>%

formattable(table_crime_costs, align =c("l","c","c","c","c", "c", "c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `Incremental Costs (1-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[1]),
  `Incremental Costs (5-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[1]),
  `Incremental Costs (10-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[1]),
  `ICER (1-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[3]),
  `ICER (5-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[3]),
  `ICER (10-year)` = color_tile(crime_costs_palette[2], crime_costs_palette[3])))

#############
### QALYs ###
#############
df_qalys_baseline_MMS <- data.frame(ICER_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                                    ICER_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_MMS$df_icer$n_icer_TOTAL_1yr, 
                                    ICER_MMS$df_icer$n_icer_TOTAL_5yr, ICER_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_reduced_eq_5d_5l_MMS <- data.frame(ICER_qalys_reduced_eq_5d_5l_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_reduced_eq_5d_5l_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                                            ICER_qalys_reduced_eq_5d_5l_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_qalys_reduced_eq_5d_5l_MMS$df_icer$n_icer_TOTAL_1yr, 
                                            ICER_qalys_reduced_eq_5d_5l_MMS$df_icer$n_icer_TOTAL_5yr, ICER_qalys_reduced_eq_5d_5l_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_low_MMS <- data.frame(ICER_qalys_low_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_low_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                               ICER_qalys_low_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_qalys_low_MMS$df_icer$n_icer_TOTAL_1yr, 
                               ICER_qalys_low_MMS$df_icer$n_icer_TOTAL_5yr, ICER_qalys_low_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_high_MMS <- data.frame(ICER_qalys_high_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_high_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                                ICER_qalys_high_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_qalys_high_MMS$df_icer$n_icer_TOTAL_1yr, 
                                ICER_qalys_high_MMS$df_icer$n_icer_TOTAL_5yr, ICER_qalys_high_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_eq_5d_3l_MMS <- data.frame(ICER_qalys_eq_5d_3l_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_eq_5d_3l_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                                    ICER_qalys_eq_5d_3l_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_qalys_eq_5d_3l_MMS$df_icer$n_icer_TOTAL_1yr, 
                                    ICER_qalys_eq_5d_3l_MMS$df_icer$n_icer_TOTAL_5yr, ICER_qalys_eq_5d_3l_MMS$df_icer$n_icer_TOTAL_10yr)
df_qalys_hui_3_MMS <- data.frame(ICER_qalys_hui_3_MMS$df_incremental$n_inc_qalys_TOTAL_1yr, ICER_qalys_hui_3_MMS$df_incremental$n_inc_qalys_TOTAL_5yr, 
                                 ICER_qalys_hui_3_MMS$df_incremental$n_inc_qalys_TOTAL_10yr, ICER_qalys_hui_3_MMS$df_icer$n_icer_TOTAL_1yr, 
                                 ICER_qalys_hui_3_MMS$df_icer$n_icer_TOTAL_5yr, ICER_qalys_hui_3_MMS$df_icer$n_icer_TOTAL_10yr)

colnames(df_qalys_baseline_MMS) <- colnames(df_qalys_reduced_eq_5d_5l_MMS) <- colnames(df_qalys_low_MMS) <- colnames(df_qalys_high_MMS) <- colnames(df_qalys_eq_5d_3l_MMS) <- colnames(df_qalys_hui_3_MMS) <- c("inc_qalys_1yr", "inc_qalys_5yr", "inc_qalys_10yr", "icer_1yr", "icer_5yr", "icer_10yr")

df_qalys <- rbind(df_qalys_baseline_MMS, df_qalys_reduced_eq_5d_5l_MMS, df_qalys_low_MMS, df_qalys_high_MMS, df_qalys_eq_5d_3l_MMS, df_qalys_hui_3_MMS)
df_qalys <- data.frame("Scenario" = c("Baseline", "Reduced", "Low", "High", "EQ-5D-3L", "HUI-3"), df_qalys)
#row.names(df_qalys) <- c("Baseline", "Reduced", "Low", "High", "EQ-5D-3L", "HUI-3")

# Custom table output
# set colours
#customGreen0 <- "#DeF7E9"
#customGreen <- "#71CA97"
#customRed <- "#ff7f7f"
qaly_palette <- brewer.pal(3,"PuOr")
#qaly_palette <- brewer.pal(3,"RdYlBu")
#qaly_palette[1]

# table
table_qalys <- df_qalys %>%
  mutate(`Incremental QALYs (1-year)` = round(inc_qalys_1yr, 3),
         `Incremental QALYs (5-year)` = round(inc_qalys_5yr, 3),
         `Incremental QALYs (10-year)` = round(inc_qalys_10yr, 3),
         `ICER (1-year)` = accounting(icer_1yr, 0),
         `ICER (5-year)` = accounting(icer_5yr, 0),
         `ICER (10-year)` = accounting(icer_10yr, 0)) %>%
  select(c(`Scenario`, `Incremental QALYs (1-year)`, `Incremental QALYs (5-year)`, `Incremental QALYs (10-year)`, `ICER (1-year)`, `ICER (5-year)`, `ICER (10-year)`)) #%>%

formattable(table_qalys, align =c("l","c","c","c","c", "c", "c"), list(
  `Scenario` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
    `Incremental QALYs (1-year)` = color_tile(qaly_palette[1], qaly_palette[2]),
  `Incremental QALYs (5-year)` = color_tile(qaly_palette[1], qaly_palette[2]),
  `Incremental QALYs (10-year)` = color_tile(qaly_palette[1], qaly_palette[2]),
  `ICER (1-year)` = color_tile(qaly_palette[2], qaly_palette[3]),
  `ICER (5-year)` = color_tile(qaly_palette[2], qaly_palette[3]),
  `ICER (10-year)` = color_tile(qaly_palette[2], qaly_palette[3])))

# Prepare data for plotting
# Costs
df_low_costs <- rbind(df_crime_costs_low_MMS$n_inc_costs_TOTAL_10yr, df_qalys_low_MMS$n_inc_costs_TOTAL_10yr)
df_high_costs <- rbind(df_crime_costs_high_MMS$n_inc_costs_TOTAL_10yr, df_qalys_high_MMS$n_inc_costs_TOTAL_10yr)
df_tornado_costs <- data.frame(df_low, df_high)

colnames(df_tornado_costs) <- c("Low", "High")
row.names(df_tornado_costs) <- c("Crime Costs", "QALYs")

# ICER
df_low_icer <- rbind(df_crime_costs_low_MMS$n_icer_TOTAL_10yr, df_qalys_low_MMS$n_icer_TOTAL_10yr)
df_high_icer <- rbind(df_crime_costs_high_MMS$n_icer_TOTAL_10yr, df_qalys_high_MMS$n_icer_TOTAL_10yr)
df_tornado_icer <- data.frame(df_low, df_high)

colnames(df_tornado_icer) <- c("Low", "High")
row.names(df_tornado_icer) <- c("Crime Costs", "QALYs")

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
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

v_dsa_qalys_low_MMS <- unlist(df_dsa_qalys_MMS["pe_low",])
v_dsa_qalys_high_MMS <- unlist(df_dsa_qalys_MMS["pe_high",])
v_dsa_qalys_eq_5d_3l_MMS <- unlist(df_dsa_qalys_MMS["pe_eq_5d_3l",])
v_dsa_qalys_hui_3_MMS <- unlist(df_dsa_qalys_MMS["pe_hui_3",])

############################################
#### Deterministic sensitivity analysis ####
############################################
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

# Alt (EQ-5d-3L)
v_qalys_eq_5d_3l_MMS <- c(ICER_qalys_eq_5d_3l_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_qalys_eq_5d_3l_MMS$df_incremental$n_inc_qalys_TOTAL_life)
# Alt (HUI-3)
v_qalys_hui_3_MMS <- c(ICER_qalys_hui_3_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_qalys_hui_3_MMS$df_incremental$n_inc_qalys_TOTAL_life)

m_qalys_MMS <- rbind(v_qalys_eq_5d_3l_MMS, v_qalys_hui_3_MMS)

df_qalys_MMS <- as.data.frame(m_qalys_MMS)
colnames(df_qalys_MMS) <- c("Lower", "Upper")
df_qalys_MMS <- as_data_frame(df_qalys_MMS) %>%
  add_column(var_name = c("Alternative instrument (EQ-5D-3L)", "Alternative instrument (HUI3)"))

## As .RData ##
save(df_qalys_MMS, 
     file = "outputs/DSA/Modified Model Specification/df_qalys_MMS.RData")

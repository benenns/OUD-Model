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

###################
### Crime costs ###
###################
# DSA data
df_dsa_crime_costs_MMS <- read.csv(file = "data/DSA/Modified Model Specification/crime_costs.csv", row.names = 1, header = TRUE)

v_dsa_crime_costs_alt_MMS <- unlist(df_dsa_crime_costs_MMS["pe_alt",])

############################################
#### Deterministic sensitivity analysis ####
############################################

# Alternative (Krebs 2014 estimates)
# MMS
l_outcomes_MET_crime_costs_alt_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_MMS)
l_outcomes_BUP_crime_costs_alt_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_crime_costs_alt_MMS)
ICER_crime_costs_alt_MMS <- ICER(outcomes_comp = l_outcomes_MET_crime_costs_alt_MMS, outcomes_int = l_outcomes_BUP_crime_costs_alt_MMS)

#############
### Costs ###
#############
v_crime_costs_alt_MMS <- c(ICER_crime_costs_alt_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_crime_costs_alt_MMS$df_incremental$n_inc_costs_TOTAL_life)

m_costs_MMS <- rbind(v_crime_costs_alt_MMS)

df_costs_MMS <- as.data.frame(m_costs_MMS)
colnames(df_costs_MMS) <- c("Lower", "Upper")
df_costs_MMS <- as_data_frame(df_costs_MMS) %>%
  add_column(var_name = c("Alternative crime cost estimates (Krebs et al. 2014)"))

# Save outputs
## As .RData ##
save(df_costs_MMS, 
     file = "outputs/DSA/Modified Model Specification/df_costs_MMS.RData")

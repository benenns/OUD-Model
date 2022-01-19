rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
#library(scales)   # for dollar signs and commas
library(tidyverse)
library(Rmisc)
library(rlist)

# Call model setup functions
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
#source("R/plot_functions.R")

# Load parameters
l_params_all <- load_all_params(file.init = "data/init_params.csv",
                                file.init_dist = "data/init_dist.csv", 
                                file.mort = "data/all_cause_mortality.csv",
                                file.death_hr = "data/death_hr.csv",
                                file.frailty = "data/frailty.csv",
                                file.weibull = "data/weibull.csv",
                                file.unconditional = "data/unconditional.csv",
                                file.overdose = "data/overdose.csv",
                                file.fentanyl = "data/fentanyl.csv",
                                file.hiv = "data/hiv_sero.csv",
                                file.hcv = "data/hcv_sero.csv",
                                file.costs = "data/costs.csv",
                                file.crime_costs = "data/crime_costs.csv",
                                file.qalys = "data/qalys.csv")

#### Load calibration targets ####
l_cali_targets <- list(ODF = read.csv(file = "data/cali_target_odf.csv", header = TRUE),
                       ODN = read.csv(file = "data/cali_target_odn.csv", header = TRUE))

# Max calibration periods
n_cali_max_per <- max(c(l_cali_targets$ODF$Time, l_cali_targets$ODN$Time))

# Run model with updated calibrated parameters
l_out_markov <- markov_model(l_params_all = l_params_all, cali = TRUE)

#### Epidemiological Output ####
### Overdose deaths ###
v_ODF <- l_out_markov$m_M_agg_trace[, "ODF"] # cumulative deaths at time i

### Non-fatal overdoses ###
v_ODN <- l_out_markov$m_M_agg_trace[, "ODN"] # cumulative non-fatal overdoses at time i

### Select time-points ###
### Overdose deaths ###
n_ODF_t1 <- l_cali_targets$ODF$Time[1]
n_ODF_t2 <- l_cali_targets$ODF$Time[2]
n_ODF_t3 <- l_cali_targets$ODF$Time[3]

### Non-fatal overdose ###
n_ODN_t1 <- l_cali_targets$ODN$Time[1]
n_ODN_t2 <- l_cali_targets$ODN$Time[2]
n_ODN_t3 <- l_cali_targets$ODN$Time[3]

### Subset output by time-points ###
### Overdose deaths ###
# Fatal overdoses already cumulative at each time point
#n_ODF1 <- v_ODF[n_ODF_t1]
#n_ODF2 <- v_ODF[n_ODF_t2]
#n_ODF3 <- v_ODF[n_ODF_t3]

# Yearly fatal overdoses (disaggregated)
n_ODF1 <- v_ODF[n_ODF_t1]
n_ODF2 <- v_ODF[n_ODF_t2] - v_ODF[n_ODF_t1]
n_ODF3 <- v_ODF[n_ODF_t3] - v_ODF[n_ODF_t2]

### Non-fatal overdose
# Non-fatal overdoses need to be summed across time points to generate cumulative estimates
#n_ODN1 <- sum(v_ODN[c(1:n_ODN_t1)])
#n_ODN2 <- sum(v_ODN[c(1:n_ODN_t2)])
#n_ODN3 <- sum(v_ODN[c(1:n_ODN_t3)])

# Yearly non-fatal overdose (disaggregated)
n_ODN1 <- sum(v_ODN[c(1:n_ODN_t1)])
n_ODN2 <- sum(v_ODN[c(13:n_ODN_t2)])
n_ODN3 <- sum(v_ODN[c(25:n_ODN_t3)])

#### Return Output ####
l_out <- list(fatal_overdose = c(n_ODF1, n_ODF2, n_ODF3), # deaths at t1, t2 , t3 time periods (for yearly deaths: (i + 12) - i where i = first month of year, 1 + 12 = last month)
              overdose = c(n_ODN1, n_ODN2, n_ODN3))




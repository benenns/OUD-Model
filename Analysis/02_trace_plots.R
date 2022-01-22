rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
#library(scales)   # for dollar signs and commas
library(tidyverse)
library(Rmisc)

# Call model setup functions
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/plot_functions.R")

# Load parameters
# Calibrated parameter values
load(file = "outputs/Calibration/imis_output.RData")

# Modified Model Specification
# BNX scenario
l_params_BUP_MMS <- load_all_params(file.init = "data/init_params.csv",
                                      file.init_dist = "data/init_dist_bup.csv",
                                      file.mort = "data/all_cause_mortality.csv",
                                      file.death_hr = "data/death_hr.csv",
                                      file.frailty = "data/frailty.csv",
                                      file.weibull = "data/Modified Model Specification/weibull.csv",
                                      file.unconditional = "data/Modified Model Specification/unconditional.csv",
                                      file.overdose = "data/overdose.csv",
                                      file.fentanyl = "data/fentanyl.csv",
                                      file.hiv = "data/hiv_sero.csv",
                                      file.hcv = "data/hcv_sero.csv",
                                      file.costs = "data/Modified Model Specification/costs.csv",
                                      file.crime_costs = "data/Modified Model Specification/crime_costs.csv",
                                      file.qalys = "data/Modified Model Specification/qalys.csv")
# Methadone scenario
l_params_MET_MMS <- load_all_params(file.init = "data/init_params.csv",
                                      file.init_dist = "data/init_dist_met.csv",
                                      file.mort = "data/all_cause_mortality.csv",
                                      file.death_hr = "data/death_hr.csv",
                                      file.frailty = "data/frailty.csv",
                                      file.weibull = "data/Modified Model Specification/weibull.csv",
                                      file.unconditional = "data/Modified Model Specification/unconditional.csv",
                                      file.overdose = "data/overdose.csv",
                                      file.fentanyl = "data/fentanyl.csv",
                                      file.hiv = "data/hiv_sero.csv",
                                      file.hcv = "data/hcv_sero.csv",
                                      file.costs = "data/Modified Model Specification/costs.csv",
                                      file.crime_costs = "data/Modified Model Specification/crime_costs.csv",
                                      file.qalys = "data/Modified Model Specification/qalys.csv")

# Trial Specification
# BNX scenario
l_params_BUP_TS <- load_all_params(file.init = "data/init_params.csv",
                                     file.init_dist = "data/init_dist_bup.csv",
                                     file.mort = "data/all_cause_mortality.csv",
                                     file.death_hr = "data/death_hr.csv",
                                     file.frailty = "data/frailty.csv",
                                     file.weibull = "data/Trial Specification/weibull.csv",
                                     file.unconditional = "data/Trial Specification/unconditional_bup.csv",
                                     file.overdose = "data/overdose.csv",
                                     file.fentanyl = "data/fentanyl.csv",
                                     file.hiv = "data/hiv_sero.csv",
                                     file.hcv = "data/hcv_sero.csv",
                                     file.costs = "data/Trial Specification/costs.csv",
                                     file.crime_costs = "data/Trial Specification/crime_costs.csv",
                                     file.qalys = "data/Trial Specification/qalys.csv")

# Methadone scenario
l_params_MET_TS <- load_all_params(file.init = "data/init_params.csv",
                                     file.init_dist = "data/init_dist_met.csv",
                                     file.mort = "data/all_cause_mortality.csv",
                                     file.death_hr = "data/death_hr.csv",
                                     file.frailty = "data/frailty.csv",
                                     file.weibull = "data/Trial Specification/weibull.csv",
                                     file.unconditional = "data/Trial Specification/unconditional_met.csv",
                                     file.overdose = "data/overdose.csv",
                                     file.fentanyl = "data/fentanyl.csv",
                                     file.hiv = "data/hiv_sero.csv",
                                     file.hcv = "data/hcv_sero.csv",
                                     file.costs = "data/Trial Specification/costs.csv",
                                     file.crime_costs = "data/Trial Specification/crime_costs.csv",
                                     file.qalys = "data/Trial Specification/qalys.csv")

# Update parameter list with calibrated params
#l_params_all <- update_param_list(l_params_all = l_params_all, params_updated = v_calib_post_map)
l_params_BUP_MMS <- update_param_list(l_params_all = l_params_BUP_MMS, params_updated = v_calib_post_map)
l_params_MET_MMS <- update_param_list(l_params_all = l_params_MET_MMS, params_updated = v_calib_post_map)

l_params_BUP_TS <- update_param_list(l_params_all = l_params_BUP_TS, params_updated = v_calib_post_map)
l_params_MET_TS <- update_param_list(l_params_all = l_params_MET_TS, params_updated = v_calib_post_map)

# Run model
#_out_markov_base <- markov_model(l_params_all = l_params_all, err_stop = FALSE, verbose = TRUE, checks = FALSE)
l_out_markov_BUP_MMS  <- markov_model(l_params_all = l_params_BUP_MMS, err_stop = FALSE, verbose = TRUE, checks = FALSE)
l_out_markov_MET_MMS  <- markov_model(l_params_all = l_params_MET_MMS, err_stop = FALSE, verbose = TRUE, checks = FALSE)

l_out_markov_BUP_TS  <- markov_model(l_params_all = l_params_BUP_TS, err_stop = FALSE, verbose = TRUE, checks = FALSE)
l_out_markov_MET_TS  <- markov_model(l_params_all = l_params_MET_TS, err_stop = FALSE, verbose = TRUE, checks = FALSE)

#### Create plots ####
#l_trace_base <- trace_plots(outcomes = l_out_markov_base)
l_trace_BUP_MMS  <- trace_plots(outcomes = l_out_markov_BUP_MMS)
l_trace_MET_MMS  <- trace_plots(outcomes = l_out_markov_MET_MMS)

l_trace_BUP_TS  <- trace_plots(outcomes = l_out_markov_BUP_TS)
l_trace_MET_TS  <- trace_plots(outcomes = l_out_markov_MET_TS)

### Outputs ###
# Model 1 (Modified model specification)
# Health state trace
# Base
#pdf("Plots/Markov Trace/Modified Model Spec/trace_states_base.pdf", width = 8, height = 6)
#l_trace_base[[1]]
#dev.off()

# Serostatus trace
# Base
#pdf("Plots/Markov Trace/Modified Model Spec/trace_sero_base.pdf", width = 8, height = 6)
#l_trace_base[[2]]
#dev.off()
# BUP
pdf("Plots/Markov Trace/Modified Model Spec/trace_sero_BUP.pdf", width = 8, height = 6)
l_trace_BUP_MMS[[2]]
dev.off()
# MET
pdf("Plots/Markov Trace/Modified Model Spec/trace_sero_MET.pdf", width = 8, height = 6)
l_trace_MET_MMS[[2]]
dev.off()

# Health state time
# Base
#pdf(file = "Plots/Markov Trace/Modified Model Spec/time_states_base.pdf", width = 8, height = 3)
#l_trace_base[[3]]
#dev.off()

# Health state trace + time
# Base
#pdf(file = "Plots/Markov Trace/Modified Model Spec/full_trace_base.pdf", width = 8, height = 9)
#multiplot(plotlist = l_trace_base[[4]], layout = l_trace_base[[5]])
#dev.off()
# BUP
pdf(file = "Plots/Markov Trace/Modified Model Spec/full_trace_BUP_MMS.pdf", width = 8, height = 9)
multiplot(plotlist = l_trace_BUP_MMS[[4]], layout = l_trace_BUP_MMS[[5]])
dev.off()
# MET
pdf(file = "Plots/Markov Trace/Modified Model Spec/full_trace_MET_MMS.pdf", width = 8, height = 9)
multiplot(plotlist = l_trace_MET_MMS[[4]], layout = l_trace_MET_MMS[[5]])
dev.off()

# Model 2 (Trial health state definition)
# Health state trace + time
# BUP
pdf(file = "Plots/Markov Trace/Trial Spec/full_trace_BUP_TS.pdf", width = 8, height = 9)
multiplot(plotlist = l_trace_BUP_TS[[4]], layout = l_trace_BUP_TS[[5]])
dev.off()
# MET
pdf(file = "Plots/Markov Trace/Trial Spec/full_trace_MET_TS.pdf", width = 8, height = 9)
multiplot(plotlist = l_trace_MET_TS[[4]], layout = l_trace_MET_TS[[5]])
dev.off()

# Model 3 (Original health state definition)
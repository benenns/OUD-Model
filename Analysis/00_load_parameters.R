library(dplyr)    # to manipulate data

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")

# Load parameters
# Calibrated parameter values
load(file = "outputs/Calibration/imis_output.RData")

######################################
#### Modified Model Specification ####
######################################
# BNX scenario
l_params_BUP_MMS <- load_all_params(file.init = "data/Modified Model Specification/init_params.csv",
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
l_params_MET_MMS <- load_all_params(file.init = "data/Modified Model Specification/init_params.csv",
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

#############################
#### Trial Specification ####
#############################
# BNX scenario
l_params_BUP_TS <- load_all_params(file.init = "data/Trial Specification/init_params.csv",
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
                                   file.costs = "data/Trial Specification/costs_bup.csv",
                                   file.crime_costs = "data/Trial Specification/crime_costs_bup.csv",
                                   file.qalys = "data/Trial Specification/qalys_bup.csv")

# Methadone scenario
l_params_MET_TS <- load_all_params(file.init = "data/Trial Specification/init_params.csv",
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
                                   file.costs = "data/Trial Specification/costs_met.csv",
                                   file.crime_costs = "data/Trial Specification/crime_costs_met.csv",
                                   file.qalys = "data/Trial Specification/qalys_met.csv")

################################
#### Original Specification ####
################################
# BNX scenario
l_params_BUP_OS <- load_all_params(file.init = "data/Original Specification/init_params.csv",
                                   file.init_dist = "data/init_dist_bup.csv",
                                   file.mort = "data/all_cause_mortality.csv",
                                   file.death_hr = "data/death_hr.csv",
                                   file.frailty = "data/frailty.csv",
                                   file.weibull = "data/Original Specification/weibull.csv",
                                   file.unconditional = "data/Original Specification/unconditional.csv",
                                   file.overdose = "data/overdose.csv",
                                   file.fentanyl = "data/fentanyl.csv",
                                   file.hiv = "data/hiv_sero.csv",
                                   file.hcv = "data/hcv_sero.csv",
                                   file.costs = "data/Original Specification/costs.csv",
                                   file.crime_costs = "data/Original Specification/crime_costs.csv",
                                   file.qalys = "data/Original Specification/qalys.csv")

# Methadone scenario
l_params_MET_OS <- load_all_params(file.init = "data/Original Specification/init_params.csv",
                                   file.init_dist = "data/init_dist_met.csv",
                                   file.mort = "data/all_cause_mortality.csv",
                                   file.death_hr = "data/death_hr.csv",
                                   file.frailty = "data/frailty.csv",
                                   file.weibull = "data/Original Specification/weibull.csv",
                                   file.unconditional = "data/Original Specification/unconditional.csv",
                                   file.overdose = "data/overdose.csv",
                                   file.fentanyl = "data/fentanyl.csv",
                                   file.hiv = "data/hiv_sero.csv",
                                   file.hcv = "data/hcv_sero.csv",
                                   file.costs = "data/Original Specification/costs.csv",
                                   file.crime_costs = "data/Original Specification/crime_costs.csv",
                                   file.qalys = "data/Original Specification/qalys.csv")
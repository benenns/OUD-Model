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
                                    file.naloxone = "data/naloxone.csv",
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
                                    file.naloxone = "data/naloxone.csv",
                                    file.hiv = "data/hiv_sero.csv",
                                    file.hcv = "data/hcv_sero.csv",
                                    file.costs = "data/Modified Model Specification/costs.csv",
                                    file.crime_costs = "data/Modified Model Specification/crime_costs.csv",
                                    file.qalys = "data/Modified Model Specification/qalys.csv")
# Mixed proportions for validation
l_params_all_validation_MMS <- load_all_params(file.init = "data/Calibration/init_params.csv",
                                               file.init_dist = "data/init_dist.csv",
                                               file.mort = "data/all_cause_mortality.csv",
                                               file.death_hr = "data/death_hr.csv",
                                               file.frailty = "data/frailty.csv",
                                               file.weibull = "data/Modified Model Specification/weibull.csv",
                                               file.unconditional = "data/Modified Model Specification/unconditional.csv",
                                               file.overdose = "data/overdose.csv",
                                               file.fentanyl = "data/Calibration/fentanyl.csv",
                                               file.naloxone = "data/naloxone.csv",
                                               file.hiv = "data/hiv_sero.csv",
                                               file.hcv = "data/hcv_sero.csv",
                                               file.costs = "data/Modified Model Specification/costs.csv",
                                               file.crime_costs = "data/Modified Model Specification/crime_costs.csv",
                                               file.qalys = "data/Modified Model Specification/qalys.csv")


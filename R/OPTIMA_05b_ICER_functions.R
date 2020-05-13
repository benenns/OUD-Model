library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(dampack)  # for CEA and calculate ICERs
library(tidyverse)

# Load parameters
l_params_all_CA <- load_all_params(file.init = "data/init_params_CA.csv",
                                   file.mort = "data/all_cause_mortality_CA.csv",
                                   file.death_hr = "data/death_hr_CA.csv",
                                   file.frailty = "data/frailty_CA.csv",
                                   file.weibull_scale = "data/weibull_scale_CA.csv",
                                   file.weibull_shape = "data/weibull_shape_CA.csv",
                                   file.unconditional = "data/unconditional_TP_CA.csv",
                                   file.sero = "data/hiv_sero_CA.csv",
                                   file.costs = "data/costs_CA.csv",
                                   file.crime_costs = "data/crime_costs_CA.csv",
                                   file.qalys = "data/qalys_CA.csv")
# Run model
l_out_markov <- markov_model(l_params_all = l_params_all_CA, err_stop = FALSE, verbose = TRUE)

v_index <- l_out_markov$v_index

with(as.list(l_params_all, l_out_markov), {
#### Cost-effectiveness analysis ####
### Matrix of individual costs
m_TX_costs <- m_HRU_costs <- m_HIV_costs <- m_crime_costs <- array(0, dim = c((n_t + 1), n_states),
                                                                   dimnames = list(0:n_t, v_n_states))
# Treatment
for (i in 1:nrow(l_out_markov$m_M_trace)){
m_TX_costs[i, all_BUP] <- l_out_markov$m_M_trace[i, all_BUP] * c_BUP_TX
m_TX_costs[i, all_MET] <- m_TX_costs[i, all_MET] * c_MET_TX
m_TX_costs[i, all_REL] <- m_TX_costs[i, all_REL] * 0
m_TX_costs[i, OD]  <- m_TX_costs[i, OD] * c_OD_TX
m_TX_costs[i, ABS]     <- m_TX_costs[i, ABS] * 0

# Health resource use
m_HRU_costs[i, all_BUP & NI] <- m_TX_costs[i, all_BUP & NI] * c_BUP_NI_HRU
m_HRU_costs[i, all_MET & NI] <- m_TX_costs[i, all_MET & NI] * c_MET_NI_HRU
m_HRU_costs[i, all_REL & NI] <- m_TX_costs[i, all_REL & NI] * c_REL_NI_HRU
m_HRU_costs[i, OD & NI]      <- m_TX_costs[i, OD & NI] * c_OD_NI_HRU
m_HRU_costs[i, ABS & NI]     <- m_TX_costs[i, ABS & NI] * c_ABS_NI_HRU

m_HRU_costs[i, all_BUP & INJ] <- m_TX_costs[i, all_BUP & INJ] * c_BUP_INJ_HRU
m_HRU_costs[i, all_MET & INJ] <- m_TX_costs[i, all_MET & INJ] * c_MET_INJ_HRU
m_HRU_costs[i, all_REL & INJ] <- m_TX_costs[i, all_REL & INJ] * c_REL_INJ_HRU
m_HRU_costs[i, OD & INJ]      <- m_TX_costs[i, OD & INJ] * c_OD_INJ_HRU
m_HRU_costs[i, ABS & INJ]     <- m_TX_costs[i, ABS & INJ] * c_ABS_INJ_HRU

# HIV
m_HIV_costs[i, POS] <- m_HIV_costs[i, POS] * (c_HIV_HRU + (c_HIV_ART * p_HIV_ART)) # can disaggregate ART costs if we want
m_HIV_costs[i, NEG] <- m_HIV_costs[i, NEG] * 0

# Crime costs
# Create vectors of crime costs for every time period
v_BUP_NI_crime <- rep(v_c_BUP_NI_crime[l_params_all$n_age_init:(l_params_all$n_age_max - 1), ], each = 12)
v_MET_NI_crime <- rep(v_c_MET_NI_crime[l_params_all$n_age_init:(l_params_all$n_age_max - 1), ], each = 12)
v_REL_NI_crime <- rep(v_c_REL_NI_crime[l_params_all$n_age_init:(l_params_all$n_age_max - 1), ], each = 12)
v_OD_NI_crime  <- rep(v_c_OD_NI_crime[l_params_all$n_age_init:(l_params_all$n_age_max - 1), ], each = 12)
v_ABS_NI_crime <- rep(v_c_ABS_NI_crime[l_params_all$n_age_init:(l_params_all$n_age_max - 1), ], each = 12)
v_BUP_INJ_crime <- rep(v_c_BUP_INJ_crime[l_params_all$n_age_init:(l_params_all$n_age_max - 1), ], each = 12)
v_MET_INJ_crime <- rep(v_c_MET_INJ_crime[l_params_all$n_age_init:(l_params_all$n_age_max - 1), ], each = 12)
v_REL_INJ_crime <- rep(v_c_REL_INJ_crime[l_params_all$n_age_init:(l_params_all$n_age_max - 1), ], each = 12)
v_OD_INJ_crime  <- rep(v_c_OD_INJ_crime[l_params_all$n_age_init:(l_params_all$n_age_max - 1), ], each = 12)
v_ABS_INJ_crime <- rep(v_c_ABS_INJ_crime[l_params_all$n_age_init:(l_params_all$n_age_max - 1), ], each = 12)

m_crime_costs[i, all_BUP & NI] <- m_crime_costs[i, all_BUP & NI] * v_BUP_NI_crime[i, ]
m_crime_costs[i, all_MET & NI] <- m_crime_costs[i, all_MET & NI] * v_MET_NI_crime[i, ]
m_crime_costs[i, all_REL & NI] <- m_crime_costs[i, all_REL & NI] * v_REL_NI_crime[i, ]
m_crime_costs[i, OD & NI]      <- m_crime_costs[i, OD & NI] * v_OD_NI_crime[i, ]
m_crime_costs[i, ABS & NI]     <- m_crime_costs[i, ABS & NI] * v_ABS_NI_crime[i, ]

m_crime_costs[i, all_BUP & INJ] <- m_crime_costs[i, all_BUP & INJ] * v_BUP_INJ_crime[i, ]
m_crime_costs[i, all_MET & INJ] <- m_crime_costs[i, all_MET & INJ] * v_MET_INJ_crime[i, ]
m_crime_costs[i, all_REL & INJ] <- m_crime_costs[i, all_REL & INJ] * v_REL_INJ_crime[i, ]
m_crime_costs[i, OD & INJ]      <- m_crime_costs[i, OD & INJ] * v_OD_INJ_crime[i, ]
m_crime_costs[i, ABS & INJ]     <- m_crime_costs[i, ABS & INJ] * v_ABS_INJ_crime[i, ]
}
}
)
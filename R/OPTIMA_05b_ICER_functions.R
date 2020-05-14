library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
#library(dampack)  # for CEA and calculate ICERs
library(tidyverse)

# Load parameters
l_params_all <- load_all_params(file.init = "data/init_params.csv",
                                file.init_dist = "data/init_dist.csv",
                                file.mort = "data/all_cause_mortality.csv",
                                file.death_hr = "data/death_hr.csv",
                                file.frailty = "data/frailty.csv",
                                file.weibull_scale = "data/weibull_scale.csv",
                                file.weibull_shape = "data/weibull_shape.csv",
                                file.unconditional = "data/unconditional_TP.csv",
                                file.sero = "data/hiv_sero.csv",
                                file.costs = "data/costs.csv",
                                file.crime_costs = "data/crime_costs.csv",
                                file.qalys = "data/qalys.csv")
# Run model
l_out_markov <- markov_model(l_params_all = l_params_all, err_stop = FALSE, verbose = TRUE)

outcomes <- function(l_params_all, l_out_markov){ # Move above "# Run model" so that model runs within function after verifying cost matrices
                                                  # Function needs to be modified so that ICERs can be calculated looping through PSA params
  
  # State indices
  l_index_s <-l_out_markov$l_index_s # Load state/strata indices to assign costs/QALYs
  all_BUP <- l_index_s$all_BUP
  BUP     <- l_index_s$BUP
  all_MET <- l_index_s$all_MET
  MET     <- l_index_s$MET
  all_REL <- l_index_s$all_REL
  REL     <- l_index_s$REL
  OD      <- l_index_s$OD
  ABS     <- l_index_s$ABS
  NI      <- l_index_s$NI
  INJ     <- l_index_s$INJ
  NEG     <- l_index_s$NEG
  POS     <- l_index_s$POS

  #### Cost-effectiveness analysis ####
  with(as.list(l_params_all), {
    ### Matrix of individual costs ###
    m_M_trace <- l_out_markov$m_M_trace # Load full markov trace output
    m_TX_costs <- m_HRU_costs <- m_HIV_costs <- m_crime_costs <- m_M_trace # All state occupancy to apply costs
    # Treatment
    for (i in 1:nrow(l_out_markov$m_M_trace)){
      m_TX_costs[i, all_BUP] <- m_TX_costs[i, all_BUP] * c_BUP_TX
      m_TX_costs[i, all_MET] <- m_TX_costs[i, all_MET] * c_MET_TX
      m_TX_costs[i, all_REL] <- m_TX_costs[i, all_REL] * 0
      m_TX_costs[i, OD]      <- m_TX_costs[i, OD] * c_OD_TX # If adding detox costs to OD
      m_TX_costs[i, ABS]     <- m_TX_costs[i, ABS] * 0

      # Health resource use
      m_HRU_costs[i, all_BUP & NI] <- m_HRU_costs[i, all_BUP & NI] * c_BUP_NI_HRU
      m_HRU_costs[i, all_MET & NI] <- m_HRU_costs[i, all_MET & NI] * c_MET_NI_HRU
      m_HRU_costs[i, all_REL & NI] <- m_HRU_costs[i, all_REL & NI] * c_REL_NI_HRU
      m_HRU_costs[i, OD & NI]      <- m_HRU_costs[i, OD & NI] * c_OD_NI_HRU
      m_HRU_costs[i, ABS & NI]     <- m_HRU_costs[i, ABS & NI] * c_ABS_NI_HRU
      m_HRU_costs[i, all_BUP & INJ] <- m_HRU_costs[i, all_BUP & INJ] * c_BUP_INJ_HRU
      m_HRU_costs[i, all_MET & INJ] <- m_HRU_costs[i, all_MET & INJ] * c_MET_INJ_HRU
      m_HRU_costs[i, all_REL & INJ] <- m_HRU_costs[i, all_REL & INJ] * c_REL_INJ_HRU
      m_HRU_costs[i, OD & INJ]      <- m_HRU_costs[i, OD & INJ] * c_OD_INJ_HRU
      m_HRU_costs[i, ABS & INJ]     <- m_HRU_costs[i, ABS & INJ] * c_ABS_INJ_HRU

      # HIV
      m_HIV_costs[i, POS] <- m_HIV_costs[i, POS] * (c_HIV_HRU + (c_HIV_ART * p_HIV_ART)) # Can disaggregate ART costs if we want
      m_HIV_costs[i, NEG] <- m_HIV_costs[i, NEG] * 0
      
      # Crime costs
      # Create vectors of crime costs for every time period
      v_BUP_NI_crime <- rep(v_c_BUP_NI_crime[n_age_init:(n_age_max - 1), ], each = 12)
      v_MET_NI_crime <- rep(v_c_MET_NI_crime[n_age_init:(n_age_max - 1), ], each = 12)
      v_REL_NI_crime <- rep(v_c_REL_NI_crime[n_age_init:(n_age_max - 1), ], each = 12)
      v_OD_NI_crime  <- rep(v_c_OD_NI_crime[n_age_init:(n_age_max - 1), ], each = 12)
      v_ABS_NI_crime <- rep(v_c_ABS_NI_crime[n_age_init:(n_age_max - 1), ], each = 12)
      v_BUP_INJ_crime <- rep(v_c_BUP_INJ_crime[n_age_init:(n_age_max - 1), ], each = 12)
      v_MET_INJ_crime <- rep(v_c_MET_INJ_crime[n_age_init:(n_age_max - 1), ], each = 12)
      v_REL_INJ_crime <- rep(v_c_REL_INJ_crime[n_age_init:(n_age_max - 1), ], each = 12)
      v_OD_INJ_crime  <- rep(v_c_OD_INJ_crime[n_age_init:(n_age_max - 1), ], each = 12)
      v_ABS_INJ_crime <- rep(v_c_ABS_INJ_crime[n_age_init:(n_age_max - 1), ], each = 12)

      m_crime_costs[i, all_BUP & NI] <- m_crime_costs[i, all_BUP & NI] * v_BUP_NI_crime[i]
      m_crime_costs[i, all_MET & NI] <- m_crime_costs[i, all_MET & NI] * v_MET_NI_crime[i]
      m_crime_costs[i, all_REL & NI] <- m_crime_costs[i, all_REL & NI] * v_REL_NI_crime[i]
      m_crime_costs[i, OD & NI]      <- m_crime_costs[i, OD & NI] * v_OD_NI_crime[i]
      m_crime_costs[i, ABS & NI]     <- m_crime_costs[i, ABS & NI] * v_ABS_NI_crime[i]
      m_crime_costs[i, all_BUP & INJ] <- m_crime_costs[i, all_BUP & INJ] * v_BUP_INJ_crime[i]
      m_crime_costs[i, all_MET & INJ] <- m_crime_costs[i, all_MET & INJ] * v_MET_INJ_crime[i]
      m_crime_costs[i, all_REL & INJ] <- m_crime_costs[i, all_REL & INJ] * v_REL_INJ_crime[i]
      m_crime_costs[i, OD & INJ]      <- m_crime_costs[i, OD & INJ] * v_OD_INJ_crime[i]
      m_crime_costs[i, ABS & INJ]     <- m_crime_costs[i, ABS & INJ] * v_ABS_INJ_crime[i]
    }
    
    m_TOTAL_costs_states = (m_TX_costs + m_HRU_costs + m_HIV_costs + m_crime_costs)
    m_HEALTH_SECTOR_costs_states = (m_TX_costs + m_HRU_costs + m_HIV_costs)
    
    v_TOTAL_costs = rowSums(m_TOTAL_costs_states)
    v_HEALTH_SECTOR_costs = rowSums(m_HEALTH_SECTOR_costs_states)
    
    ### QALY weights ###
    
    
    
    

    return(list(m_TX_costs = m_TX_costs,
                m_HRU_costs = m_HRU_costs,
                m_HIV_costs = m_HIV_costs,
                m_crime_costs = m_crime_costs,
                m_TOTAL_costs_states = m_TOTAL_costs_states,
                m_HEALTH_SECTOR_costs_states = m_HEALTH_SECTOR_costs_states))
  }
 )
}


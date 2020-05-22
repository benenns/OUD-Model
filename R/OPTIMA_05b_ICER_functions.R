#' Outcomes
#'
#' \code{outcomes} implements functions to apply costs & QALYs to Markov trace.
#'
#' @param l_params_all List with all parameters
#' @param markov_model Run Markov model with l_params_all as inputs
#' @return 
#' m_TX_costs: Matrix of treatment costs.
#' m_HRU_costs: Matrix of health resource use costs.
#' m_HIV_costs: Matrix of HIV-related costs.
#' m_crime_costs: Matrix of crime costs.
#' m_TOTAL_costs_states: Matrix of all costs.
#' m_HEALTH_SECTOR_costs_states: Matrix of costs excluding crime costs.
#' v_TOTAL_costs: Total costs summed across health states.
#' v_HEALTH_SECTOR_costs: Health sector costs summed across health states.
#' v_TOTAL_qalys: Total QALYs summed across health states.
#' n_TOTAL_costs_1yr:
#' n_TOTAL_costs_5yr:
#' n_TOTAL_costs_10yr:
#' n_TOTAL_costs_life:
#' n_HEALTH_SECTOR_costs_1yr:
#' n_HEALTH_SECTOR_costs_5yr:
#' n_HEALTH_SECTOR_costs_10yr:
#' n_HEALTH_SECTOR_costs_life:
#' n_TOTAL_qalys_1yr:
#' n_TOTAL_qalys_5yr:
#' n_TOTAL_qalys_10yr:
#' n_TOTAL_qalys_life:
#' @export
outcomes <- function(l_params_all){
  
  # Run model
  l_out_markov <- markov_model(l_params_all = l_params_all, err_stop = FALSE, verbose = TRUE)
  
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
    ### Matrices of individual costs/QALYs ###
    m_M_trace <- l_out_markov$m_M_trace # Load full markov trace output
    m_TX_costs <- m_HRU_costs <- m_HIV_costs <- m_crime_costs <- m_qalys <- m_M_trace # All state occupancy to apply costs
    
    # Calculate monthly discount rate
    n_mo_discount <- (1 + l_params_all$n_discount)^(1/12) - 1
    
    ### Costs ###
    # Treatment
    for (i in 1:nrow(l_out_markov$m_M_trace)){
      m_TX_costs[i, all_BUP] <- m_TX_costs[i, all_BUP] * c_BUP_TX
      m_TX_costs[i, all_MET] <- m_TX_costs[i, all_MET] * c_MET_TX
      m_TX_costs[i, all_REL] <- m_TX_costs[i, all_REL] * 0
      m_TX_costs[i, OD]      <- m_TX_costs[i, OD] * c_OD_TX # If adding detox costs to OD
      m_TX_costs[i, ABS]     <- m_TX_costs[i, ABS] * 0
      
      m_TX_costs[i, ] <- m_TX_costs[i, ] / ((1 + n_mo_discount)^i) # apply monthly discount

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
      
      m_HRU_costs[i, ] <- m_HRU_costs[i, ] / ((1 + n_mo_discount)^i) # apply monthly discount

      # HIV
      m_HIV_costs[i, POS] <- m_HIV_costs[i, POS] * (c_HIV_HRU + (c_HIV_ART * n_HIV_ART)) # Can disaggregate ART costs
      m_HIV_costs[i, NEG] <- m_HIV_costs[i, NEG] * 0
      
      m_HIV_costs[i, ] <- m_HIV_costs[i, ] / ((1 + n_mo_discount)^i) # apply monthly discount
      
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
      
      m_crime_costs[i, ] <- m_crime_costs[i, ] / ((1 + n_mo_discount)^i) # apply monthly discount
      
      ### QALY weights ###
      # HIV negative
      m_qalys[i, all_BUP & NI & NEG] <- m_qalys[i, all_BUP & NI & NEG] * u_BUP_NI_NEG
      m_qalys[i, all_MET & NI & NEG] <- m_qalys[i, all_MET & NI & NEG] * u_MET_NI_NEG
      m_qalys[i, all_REL & NI & NEG] <- m_qalys[i, all_REL & NI & NEG] * u_REL_NI_NEG
      m_qalys[i, OD & NI & NEG]  <- m_qalys[i, OD & NI & NEG] * u_OD_NI_NEG
      m_qalys[i, ABS & NI & NEG] <- m_qalys[i, ABS & NI & NEG] * u_BUP_NI_NEG
      
      m_qalys[i, all_BUP & INJ & NEG] <- m_qalys[i, all_BUP & INJ & NEG] * u_BUP_INJ_NEG
      m_qalys[i, all_MET & INJ & NEG] <- m_qalys[i, all_MET & INJ & NEG] * u_MET_INJ_NEG
      m_qalys[i, all_REL & INJ & NEG] <- m_qalys[i, all_REL & INJ & NEG] * u_REL_INJ_NEG
      m_qalys[i, OD & INJ & NEG]  <- m_qalys[i, OD & INJ & NEG] * u_OD_INJ_NEG
      m_qalys[i, ABS & INJ & NEG] <- m_qalys[i, ABS & INJ & NEG] * u_BUP_INJ_NEG
      
      # HIV positive
      m_qalys[i, all_BUP & NI & POS] <- m_qalys[i, all_BUP & NI & POS] * u_BUP_NI_POS
      m_qalys[i, all_MET & NI & POS] <- m_qalys[i, all_MET & NI & POS] * u_MET_NI_POS
      m_qalys[i, all_REL & NI & POS] <- m_qalys[i, all_REL & NI & POS] * u_REL_NI_POS
      m_qalys[i, OD & NI & POS]  <- m_qalys[i, OD & NI & POS] * u_OD_NI_POS
      m_qalys[i, ABS & NI & POS] <- m_qalys[i, ABS & NI & POS] * u_BUP_NI_POS
      
      m_qalys[i, all_BUP & INJ & POS] <- m_qalys[i, all_BUP & INJ & POS] * u_BUP_INJ_POS
      m_qalys[i, all_MET & INJ & POS] <- m_qalys[i, all_MET & INJ & POS] * u_MET_INJ_POS
      m_qalys[i, all_REL & INJ & POS] <- m_qalys[i, all_REL & INJ & POS] * u_REL_INJ_POS
      m_qalys[i, OD & INJ & POS]  <- m_qalys[i, OD & INJ & POS] * u_OD_INJ_POS
      m_qalys[i, ABS & INJ & POS] <- m_qalys[i, ABS & INJ & POS] * u_BUP_INJ_POS
      
      m_qalys[i, ] <- m_qalys[i, ] / ((1 + n_mo_discount)^i) # apply monthly discount
    }
    
    m_TOTAL_costs_states = (m_TX_costs + m_HRU_costs + m_HIV_costs + m_crime_costs)
    m_HEALTH_SECTOR_costs_states = (m_TX_costs + m_HRU_costs + m_HIV_costs)
    
    v_TOTAL_costs = rowSums(m_TOTAL_costs_states)
    v_HEALTH_SECTOR_costs = rowSums(m_HEALTH_SECTOR_costs_states)
    
    v_TOTAL_qalys = rowSums(m_qalys)
    
    ### Sum costs ###
    ## Societal ##
    # 1-year
    n_TOTAL_costs_1yr <- sum(v_TOTAL_costs[1:12])
    # 5-year
    n_TOTAL_costs_5yr <- sum(v_TOTAL_costs[1:60])
    # 10-year
    n_TOTAL_costs_10yr <- sum(v_TOTAL_costs[1:120])
    # Lifetime
    n_TOTAL_costs_life <- sum(v_TOTAL_costs[1:n_t])
    ## Health sector ##
    # 1-year
    n_HEALTH_SECTOR_costs_1yr <- sum(v_HEALTH_SECTOR_costs[1:12])
    # 5-year
    n_HEALTH_SECTOR_costs_5yr <- sum(v_HEALTH_SECTOR_costs[1:60])
    # 10-year
    n_HEALTH_SECTOR_costs_10yr <- sum(v_HEALTH_SECTOR_costs[1:120])
    # Lifetime
    n_HEALTH_SECTOR_costs_life <- sum(v_HEALTH_SECTOR_costs[1:n_t])
    
    ### Sum QALYs ###
    # 1-year
    n_TOTAL_qalys_1yr <- sum(v_TOTAL_qalys[1:12])
    # 5-year
    n_TOTAL_qalys_5yr <- sum(v_TOTAL_qalys[1:60])
    # 10-year
    n_TOTAL_qalys_10yr <- sum(v_TOTAL_qalys[1:120])
    # Lifetime
    n_TOTAL_qalys_life <- sum(v_TOTAL_qalys[1:n_t])

    return(list(m_TX_costs = m_TX_costs,
                m_HRU_costs = m_HRU_costs,
                m_HIV_costs = m_HIV_costs,
                m_crime_costs = m_crime_costs,
                m_TOTAL_costs_states = m_TOTAL_costs_states,
                m_HEALTH_SECTOR_costs_states = m_HEALTH_SECTOR_costs_states,
                v_TOTAL_costs = v_TOTAL_costs,
                v_HEALTH_SECTOR_costs = v_HEALTH_SECTOR_costs,
                v_TOTAL_qalys = v_TOTAL_qalys,
                n_TOTAL_costs_1yr = n_TOTAL_costs_1yr,
                n_TOTAL_costs_5yr = n_TOTAL_costs_5yr,
                n_TOTAL_costs_10yr = n_TOTAL_costs_10yr,
                n_TOTAL_costs_life = n_TOTAL_costs_life,
                n_HEALTH_SECTOR_costs_1yr = n_HEALTH_SECTOR_costs_1yr,
                n_HEALTH_SECTOR_costs_5yr = n_HEALTH_SECTOR_costs_5yr,
                n_HEALTH_SECTOR_costs_10yr = n_HEALTH_SECTOR_costs_10yr,
                n_HEALTH_SECTOR_costs_life = n_HEALTH_SECTOR_costs_life,
                n_TOTAL_qalys_1yr = n_TOTAL_qalys_1yr,
                n_TOTAL_qalys_5yr = n_TOTAL_qalys_5yr,
                n_TOTAL_qalys_10yr = n_TOTAL_qalys_10yr,
                n_TOTAL_qalys_life = n_TOTAL_qalys_life))
  }
 )
}

#' ICER
#'
#' \code{ICER} calculate ICERs.
#'
#' @param outcomes_comp List with outputs from comparator
#' @param outcomes_int List with outputs from intervention
#' @return 
#' m_TX_costs: Matrix of treatment costs.
#' m_HRU_costs: Matrix of health resource use costs.
#' m_HIV_costs: Matrix of HIV-related costs.
#' m_crime_costs: Matrix of crime costs.
#' m_TOTAL_costs_states: Matrix of all costs.
#' m_HEALTH_SECTOR_costs_states: Matrix of costs excluding crime costs.
#' v_TOTAL_costs: Total costs summed across health states.
#' v_HEALTH_SECTOR_costs: Health sector costs summed across health states.
#' v_TOTAL_qalys: Total QALYs summed across health states.
#' n_TOTAL_costs_1yr:
#' n_TOTAL_costs_5yr:
#' n_TOTAL_costs_10yr:
#' n_TOTAL_costs_life:
#' n_HEALTH_SECTOR_costs_1yr:
#' n_HEALTH_SECTOR_costs_5yr:
#' n_HEALTH_SECTOR_costs_10yr:
#' n_HEALTH_SECTOR_costs_life:
#' n_TOTAL_qalys_1yr:
#' n_TOTAL_qalys_5yr:
#' n_TOTAL_qalys_10yr:
#' n_TOTAL_qalys_life:
#' @export
ICER <- function(outcomes_comp, outcomes_int){
  ### Calculate ICERs ###
  ## Societal ##
  # 1-year
  n_icer_TOTAL_1yr <- (outcomes_int$n_TOTAL_costs_1yr - outcomes_comp$n_TOTAL_costs_1yr)/(outcomes_int$n_TOTAL_qalys_1yr - outcomes_comp$n_TOTAL_qalys_1yr)
  # 5-year
  n_icer_TOTAL_5yr <- (outcomes_int$n_TOTAL_costs_5yr - outcomes_comp$n_TOTAL_costs_5yr)/(outcomes_int$n_TOTAL_qalys_5yr - outcomes_comp$n_TOTAL_qalys_5yr)
  # 10-year
  n_icer_TOTAL_10yr <- (outcomes_int$n_TOTAL_costs_10yr - outcomes_comp$n_TOTAL_costs_10yr)/(outcomes_int$n_TOTAL_qalys_10yr - outcomes_comp$n_TOTAL_qalys_10yr)
  # Lifetime
  n_icer_TOTAL_life <- (outcomes_int$n_TOTAL_costs_life - outcomes_comp$n_TOTAL_costs_life)/(outcomes_int$n_TOTAL_qalys_10yr - outcomes_comp$n_TOTAL_qalys_10yr)
  ## Health sector ##
  # 1-year
  n_icer_HEALTH_SECTOR_1yr <- (outcomes_int$n_HEALTH_SECTOR_costs_1yr - outcomes_comp$n_HEALTH_SECTOR_costs_1yr)/(outcomes_int$n_TOTAL_qalys_1yr - outcomes_comp$n_TOTAL_qalys_1yr)
  # 5-year
  n_icer_HEALTH_SECTOR_5yr <- (outcomes_int$n_HEALTH_SECTOR_costs_5yr - outcomes_comp$n_HEALTH_SECTOR_costs_5yr)/(outcomes_int$n_TOTAL_qalys_5yr - outcomes_comp$n_TOTAL_qalys_5yr)
  # 10-year
  n_icer_HEALTH_SECTOR_10yr <- (outcomes_int$n_HEALTH_SECTOR_costs_10yr - outcomes_comp$n_HEALTH_SECTOR_costs_10yr)/(outcomes_int$n_TOTAL_qalys_10yr - outcomes_comp$n_TOTAL_qalys_10yr)
  # Lifetime
  n_icer_HEALTH_SECTOR_life <- (outcomes_int$n_HEALTH_SECTOR_costs_life - outcomes_comp$n_HEALTH_SECTOR_costs_life)/(outcomes_int$n_TOTAL_qalys_10yr - outcomes_comp$n_TOTAL_qalys_10yr)

  return(list(n_icer_TOTAL_1yr = n_icer_TOTAL_1yr,
              n_icer_TOTAL_5yr = n_icer_TOTAL_5yr,
              n_icer_TOTAL_10yr = n_icer_TOTAL_10yr,
              n_icer_TOTAL_life = n_icer_TOTAL_life,
              n_icer_HEALTH_SECTOR_1yr = n_icer_HEALTH_SECTOR_1yr,
              n_icer_HEALTH_SECTOR_5yr = n_icer_HEALTH_SECTOR_5yr,
              n_icer_HEALTH_SECTOR_10yr = n_icer_HEALTH_SECTOR_10yr,
              n_icer_HEALTH_SECTOR_life = n_icer_HEALTH_SECTOR_life))  
}

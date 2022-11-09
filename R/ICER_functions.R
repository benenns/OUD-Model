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
outcomes <- function(l_params_all, 
                     v_params_calib,
                     v_params_dsa = NULL,
                     PSA = FALSE){
  
  # Substitute values of calibrated parameters in base-case with calibrated values
  # Combine with DSA ranges
  dsa_check <- !is.null(v_params_dsa)
  
  if(dsa_check){
  v_params_updated <- c(v_params_dsa, v_params_calib)
  } else{
    v_params_updated <- v_params_calib
  }
  
  #l_params_all <- update_param_list(l_params_all = l_params_all, params_updated = v_params_updated)
  # PSA update already includes draws from posterior for cali parameters
  if(PSA){
    l_params_all <- l_params_all
  } else{
    l_params_all <- update_param_list(l_params_all = l_params_all, params_updated = v_params_updated)
  }
  
  # Substitute values of calibrated parameters in base-case with calibrated values 
  #l_params_all <- update_param_list(l_params_all = l_params_all, params_updated = v_params_calib)

  # Run model
  l_out_markov <- markov_model(l_params_all = l_params_all, err_stop = FALSE, verbose = TRUE)
  
  # State indices
  l_index_s <-l_out_markov$l_index_s # Load state/strata indices to assign costs/QALYs
  BUP     <- l_index_s$BUP
  BUPC    <- l_index_s$BUPC
  all_BUP <- l_index_s$all_BUP
  MET     <- l_index_s$MET
  METC    <- l_index_s$METC
  all_MET <- l_index_s$all_MET
  REL     <- l_index_s$REL
  ODN     <- l_index_s$ODN
  ODF     <- l_index_s$ODF
  ABS     <- l_index_s$ABS
  NI      <- l_index_s$NI
  INJ     <- l_index_s$INJ
  NEG     <- l_index_s$NEG
  HIV     <- l_index_s$HIV
  HCV     <- l_index_s$HCV
  COI     <- l_index_s$COI
  all_HIV <- l_index_s$all_HIV
  all_HCV <- l_index_s$all_HCV

  #### Cost-effectiveness analysis ####
  with(as.list(l_params_all), {
    ### Matrices of individual costs/QALYs ###
    a_M_trace <- l_out_markov$a_M_trace
    m_M_trace <- l_out_markov$m_M_trace # Load full markov trace output
    m_M_agg_trace <- l_out_markov$m_M_agg_trace # Load aggregated trace output
    m_M_agg_trace_sero <- l_out_markov$m_M_agg_trace_sero # Load serostatus trace output
    m_TX_costs <- m_HRU_costs <- m_HIV_costs <- m_HCV_costs <- m_crime_costs <- m_qalys <- m_M_trace # All state occupancy to apply costs
    
    # Calculate periodic discount rate
    n_per_discount <- (1 + l_params_all$n_discount)^(1/n_per) - 1 # Convert yearly rate to model periods
    
    ### Costs ###
    # Treatment
    for (i in 1:nrow(l_out_markov$m_M_trace)){
      m_TX_costs[i, BUP]     <- m_TX_costs[i, BUP] * c_BUP_TX
      m_TX_costs[i, BUPC]    <- m_TX_costs[i, BUPC] * c_BUP_TX
      m_TX_costs[i, MET]     <- m_TX_costs[i, MET] * c_MET_TX
      m_TX_costs[i, METC]    <- m_TX_costs[i, METC] * c_MET_TX
      m_TX_costs[i, REL]     <- m_TX_costs[i, REL] * 0
      m_TX_costs[i, ODN]     <- m_TX_costs[i, ODN] * 0 # Change if adding detox costs to OD
      m_TX_costs[i, ODF]     <- m_TX_costs[i, ODF] * 0
      m_TX_costs[i, ABS]     <- m_TX_costs[i, ABS] * 0
      
      m_TX_costs[i, ] <- m_TX_costs[i, ] / ((1 + n_per_discount)^i) # apply monthly discount

      # Health resource use
      m_HRU_costs[i, BUP & NI]     <- m_HRU_costs[i, BUP & NI] * c_BUP_NI_HRU
      m_HRU_costs[i, BUPC & NI]    <- m_HRU_costs[i, BUPC & NI] * c_BUPC_NI_HRU
      m_HRU_costs[i, MET & NI]     <- m_HRU_costs[i, MET & NI] * c_MET_NI_HRU
      m_HRU_costs[i, METC & NI]    <- m_HRU_costs[i, METC & NI] * c_METC_NI_HRU
      m_HRU_costs[i, REL & NI]     <- m_HRU_costs[i, REL & NI] * c_REL_NI_HRU
      m_HRU_costs[i, ODN & NI]     <- m_HRU_costs[i, ODN & NI] * (c_ODN_NI_HRU + (c_OD_NX * p_witness * p_NX_used) + (c_OD_AMB * p_attended) + (c_OD_TX * p_attended))
      m_HRU_costs[i, ODF & NI]     <- m_HRU_costs[i, ODF & NI] * (c_ODN_NI_HRU + (c_OD_NX * p_witness * p_NX_used) + (c_OD_AMB * p_attended))
      m_HRU_costs[i, ABS & NI]     <- m_HRU_costs[i, ABS & NI] * c_ABS_NI_HRU
      m_HRU_costs[i, BUP & INJ]    <- m_HRU_costs[i, BUP & INJ] * c_BUP_INJ_HRU
      m_HRU_costs[i, BUPC & INJ]   <- m_HRU_costs[i, BUPC & INJ] * c_BUPC_INJ_HRU
      m_HRU_costs[i, MET & INJ]    <- m_HRU_costs[i, MET & INJ] * c_MET_INJ_HRU
      m_HRU_costs[i, METC & INJ]   <- m_HRU_costs[i, METC & INJ] * c_METC_INJ_HRU
      m_HRU_costs[i, REL & INJ]    <- m_HRU_costs[i, REL & INJ] * c_REL_INJ_HRU
      m_HRU_costs[i, ODN & INJ]    <- m_HRU_costs[i, ODN & INJ] * (c_ODN_INJ_HRU + (c_OD_NX * p_witness * p_NX_used) + (c_OD_AMB * p_attended) + (c_OD_TX * p_attended))
      m_HRU_costs[i, ODF & INJ]    <- m_HRU_costs[i, ODF & INJ] * (c_ODN_INJ_HRU + (c_OD_NX * p_witness * p_NX_used) + (c_OD_AMB * p_attended))
      m_HRU_costs[i, ABS & INJ]    <- m_HRU_costs[i, ABS & INJ] * c_ABS_INJ_HRU
      
      ## ADD COSTS FOR FATAL OVERDOSE AND DEBUG vvv
      # For fatal overdose costs (need to use i - 1 since trace is cumulative)
      if(i == 1){
        m_HRU_costs[i, ODF & NI]  <- m_HRU_costs[i, ODF & NI] #* (c_OD_NI_HRU + (c_OD_NX * p_witness * p_NX_used) + (c_OD_AMB * p_attended))
        m_HRU_costs[i, ODF & INJ] <- m_HRU_costs[i, ODF & INJ] #* (c_OD_INJ_HRU + (c_OD_NX * p_witness * p_NX_used) + (c_OD_AMB * p_attended))
      } else{
        m_HRU_costs[i, ODF & NI]  <- (m_HRU_costs[i, ODF & NI] - m_HRU_costs[i - 1, ODF & NI]) #* (c_OD_NI_HRU + (c_OD_NX * p_witness * p_NX_used) + (c_OD_AMB * p_attended))
        m_HRU_costs[i, ODF & INJ] <- (m_HRU_costs[i, ODF & INJ] - m_HRU_costs[i - 1, ODF & INJ]) #* (c_OD_INJ_HRU + (c_OD_NX * p_witness * p_NX_used) + (c_OD_AMB * p_attended))
      }
      
      m_HRU_costs[i, ] <- m_HRU_costs[i, ] / ((1 + n_per_discount)^i) # apply monthly discount

      # HIV and HCV costs (applied to all HIV and HCV, i.e. including co-infection for both)
      # HIV
      m_HIV_costs[i, all_HIV] <- m_HIV_costs[i, all_HIV] * (c_HIV_HRU + (c_HIV_ART * n_HIV_ART)) # Can disaggregate ART costs
      m_HIV_costs[i, NEG] <- m_HIV_costs[i, NEG] * 0
      
      m_HIV_costs[i, ] <- m_HIV_costs[i, ] / ((1 + n_per_discount)^i) # apply monthly discount
      
      # HCV
      m_HCV_costs[i, all_HCV] <- m_HCV_costs[i, all_HCV] * (c_HCV_HRU + (c_HCV_DAA * n_HCV_DAA))
      m_HCV_costs[i, NEG] <- m_HCV_costs[i, NEG] * 0
      
      m_HCV_costs[i, ] <- m_HCV_costs[i, ] / ((1 + n_per_discount)^i) # apply monthly discount
      
      # Crime costs
      # Create vectors of crime costs for every time period
      m_crime_costs[i, BUP & NI]  <- m_crime_costs[i, BUP & NI] * c_BUP_NI_crime
      m_crime_costs[i, BUPC & NI] <- m_crime_costs[i, BUPC & NI] * c_BUPC_NI_crime
      m_crime_costs[i, MET & NI]  <- m_crime_costs[i, MET & NI] * c_MET_NI_crime
      m_crime_costs[i, METC & NI] <- m_crime_costs[i, METC & NI] * c_METC_NI_crime
      m_crime_costs[i, REL & NI]  <- m_crime_costs[i, REL & NI] * c_REL_NI_crime
      m_crime_costs[i, ODN & NI]  <- m_crime_costs[i, ODN & NI] * c_ODN_NI_crime
      m_crime_costs[i, ODF & NI]  <- m_crime_costs[i, ODF & NI] * 0
      m_crime_costs[i, ABS & NI]  <- m_crime_costs[i, ABS & NI] * c_ABS_NI_crime
      m_crime_costs[i, BUP & INJ] <- m_crime_costs[i, BUP & INJ] * c_BUP_INJ_crime
      m_crime_costs[i, BUPC & INJ] <- m_crime_costs[i, BUPC & INJ] * c_BUPC_INJ_crime
      m_crime_costs[i, MET & INJ]  <- m_crime_costs[i, MET & INJ] * c_MET_INJ_crime
      m_crime_costs[i, METC & INJ] <- m_crime_costs[i, METC & INJ] * c_METC_INJ_crime
      m_crime_costs[i, REL & INJ]  <- m_crime_costs[i, REL & INJ] * c_REL_INJ_crime
      m_crime_costs[i, ODN & INJ]  <- m_crime_costs[i, ODN & INJ] * c_ODN_INJ_crime
      m_crime_costs[i, ODF & INJ]  <- m_crime_costs[i, ODF & INJ] * 0
      m_crime_costs[i, ABS & INJ]  <- m_crime_costs[i, ABS & INJ] * c_ABS_INJ_crime
      
      m_crime_costs[i, ] <- m_crime_costs[i, ] / ((1 + n_per_discount)^i) # apply monthly discount
      
      ### QALY weights ###
      # HIV/HCV negative
      # Apply state-specific weights to all states regardless of serostatus
      m_qalys[i, BUP & NI] <- m_qalys[i, BUP & NI] * u_BUP_NI_NEG # different QALYs for BUP and BUPC?
      m_qalys[i, BUPC & NI] <- m_qalys[i, BUPC & NI] * u_BUPC_NI_NEG
      m_qalys[i, MET & NI] <- m_qalys[i, MET & NI] * u_MET_NI_NEG # different QALYs for MET and METC?
      m_qalys[i, METC & NI] <- m_qalys[i, METC & NI] * u_METC_NI_NEG # different QALYs for MET and METC?
      m_qalys[i, REL & NI] <- m_qalys[i, REL & NI] * u_REL_NI_NEG
      m_qalys[i, ODN & NI] <- m_qalys[i, ODN & NI] * u_ODN_NI_NEG
      m_qalys[i, ODF & NI] <- m_qalys[i, ODF & NI] * u_ODF_NI_NEG
      m_qalys[i, ABS & NI] <- m_qalys[i, ABS & NI] * u_BUP_NI_NEG
      
      m_qalys[i, BUP & INJ] <- m_qalys[i, BUP & INJ] * u_BUP_INJ_NEG
      m_qalys[i, BUPC & INJ] <- m_qalys[i, BUPC & INJ] * u_BUPC_INJ_NEG
      m_qalys[i, MET & INJ] <- m_qalys[i, MET & INJ] * u_MET_INJ_NEG
      m_qalys[i, METC & INJ] <- m_qalys[i, METC & INJ] * u_METC_INJ_NEG
      m_qalys[i, REL & INJ] <- m_qalys[i, REL & INJ] * u_REL_INJ_NEG
      m_qalys[i, ODN & INJ] <- m_qalys[i, ODN & INJ] * u_ODN_INJ_NEG
      m_qalys[i, ODF & INJ] <- m_qalys[i, ODF & INJ] * u_ODF_INJ_NEG
      m_qalys[i, ABS & INJ] <- m_qalys[i, ABS & INJ] * u_BUP_INJ_NEG
      
      # HIV/HCV positive
      # HIV only
      m_qalys[i, HIV] <- m_qalys[i, HIV] * u_HIV_mult
      # HCV only
      m_qalys[i, HCV] <- m_qalys[i, HCV] * u_HCV_mult
      # HIV-HCV co-infection
      m_qalys[i, COI] <- m_qalys[i, COI] * u_COI_mult
      
      # Adjust yearly QALY weights to monthly
      m_qalys[i, ] <- m_qalys[i, ] / 12
      
      # Apply monthly discount
      m_qalys[i, ] <- m_qalys[i, ] / ((1 + n_per_discount)^i) 
    }
    
    m_TOTAL_costs_states = (m_TX_costs + m_HRU_costs + m_HIV_costs + m_HCV_costs + m_crime_costs)
    m_HEALTH_SECTOR_costs_states = (m_TX_costs + m_HRU_costs + m_HIV_costs + m_HCV_costs)
    
    v_TOTAL_costs = rowSums(m_TOTAL_costs_states)
    v_HEALTH_SECTOR_costs = rowSums(m_HEALTH_SECTOR_costs_states)
    v_CRIMINAL_costs = rowSums(m_crime_costs)
    v_TX_costs = rowSums(m_TX_costs)
    v_HRU_costs = rowSums(m_HRU_costs)
    
    v_TOTAL_qalys = rowSums(m_qalys)
    
    ### Sum costs ###
    # Max model periods
    n_t <- (n_age_max - n_age_init) * n_per
    
    ## Societal ##
    # 6-month
    n_TOTAL_costs_6mo <- sum(v_TOTAL_costs[1:6])
    # 1-year
    n_TOTAL_costs_1yr <- sum(v_TOTAL_costs[1:12])
    # 5-year
    n_TOTAL_costs_5yr <- sum(v_TOTAL_costs[1:60])
    # 10-year
    n_TOTAL_costs_10yr <- sum(v_TOTAL_costs[1:120])
    # Lifetime
    n_TOTAL_costs_life <- sum(v_TOTAL_costs[1:n_t])
    
    ## Health sector ##
    # 6-month
    n_HEALTH_SECTOR_costs_6mo <- sum(v_HEALTH_SECTOR_costs[1:6])
    # 1-year
    n_HEALTH_SECTOR_costs_1yr <- sum(v_HEALTH_SECTOR_costs[1:12])
    # 5-year
    n_HEALTH_SECTOR_costs_5yr <- sum(v_HEALTH_SECTOR_costs[1:60])
    # 10-year
    n_HEALTH_SECTOR_costs_10yr <- sum(v_HEALTH_SECTOR_costs[1:120])
    # Lifetime
    n_HEALTH_SECTOR_costs_life <- sum(v_HEALTH_SECTOR_costs[1:n_t])
    
    ## Criminal justice costs ##
    # 6-month
    n_CRIMINAL_costs_6mo <- sum(v_CRIMINAL_costs[1:6])
    # 1-year
    n_CRIMINAL_costs_1yr <- sum(v_CRIMINAL_costs[1:12])
    # 5-year
    n_CRIMINAL_costs_5yr <- sum(v_CRIMINAL_costs[1:60])
    # 10-year
    n_CRIMINAL_costs_10yr <- sum(v_CRIMINAL_costs[1:120])
    # Lifetime
    n_CRIMINAL_costs_life <- sum(v_CRIMINAL_costs[1:n_t])
    
    ## Treatment costs ##
    # 6-month
    n_TX_costs_6mo <- sum(v_TX_costs[1:6])
    # 1-year
    n_TX_costs_1yr <- sum(v_TX_costs[1:12])
    # 5-year
    n_TX_costs_5yr <- sum(v_TX_costs[1:60])
    # 10-year
    n_TX_costs_10yr <- sum(v_TX_costs[1:120])
    # Lifetime
    n_TX_costs_life <- sum(v_TX_costs[1:n_t])
    
    ## HRU costs ##
    # 6-month
    n_HRU_costs_6mo <- sum(v_HRU_costs[1:6])
    # 1-year
    n_HRU_costs_1yr <- sum(v_HRU_costs[1:12])
    # 5-year
    n_HRU_costs_5yr <- sum(v_HRU_costs[1:60])
    # 10-year
    n_HRU_costs_10yr <- sum(v_HRU_costs[1:120])
    # Lifetime
    n_HRU_costs_life <- sum(v_HRU_costs[1:n_t])
    
    ## Combined costs ##
    v_costs <- c(n_TOTAL_costs_6mo, n_TOTAL_costs_1yr, n_TOTAL_costs_5yr, n_TOTAL_costs_10yr, n_TOTAL_costs_life, 
                 n_HEALTH_SECTOR_costs_6mo, n_HEALTH_SECTOR_costs_1yr, n_HEALTH_SECTOR_costs_5yr, n_HEALTH_SECTOR_costs_10yr, n_HEALTH_SECTOR_costs_life,
                 n_CRIMINAL_costs_6mo, n_CRIMINAL_costs_1yr, n_CRIMINAL_costs_5yr, n_CRIMINAL_costs_10yr, n_CRIMINAL_costs_life, 
                 n_TX_costs_6mo, n_TX_costs_1yr, n_TX_costs_5yr, n_TX_costs_10yr, n_TX_costs_life,
                 n_HRU_costs_6mo, n_HRU_costs_1yr, n_HRU_costs_5yr, n_HRU_costs_10yr, n_HRU_costs_life)
    #names(v_costs) <- c("Total Costs (1-year)", "Total Costs (5-year)", "Total Costs (10-year)", "Total Costs (Lifetime)", "Health Sector Costs (1-year)", "Health Sector Costs (5-year)", "Health Sector Costs (10-year)", "Health Sector Costs (Lifetime)",
    #                    "Criminal Costs (1-year)", "Criminal Costs (5-year)", "Criminal Costs (10-year)", "Criminal Costs (Lifetime)", "Treatment Costs (1-year)", "Treatment Costs (5-year)", "Treatment Costs (10-year)", "Treatment Costs (Lifetime)")
    
    ### Sum QALYs ###
    # 6-month
    n_TOTAL_qalys_6mo <- sum(v_TOTAL_qalys[1:6])
    # 1-year
    n_TOTAL_qalys_1yr <- sum(v_TOTAL_qalys[1:12])
    # 5-year
    n_TOTAL_qalys_5yr <- sum(v_TOTAL_qalys[1:60])
    # 10-year
    n_TOTAL_qalys_10yr <- sum(v_TOTAL_qalys[1:120])
    # Lifetime
    n_TOTAL_qalys_life <- sum(v_TOTAL_qalys[1:n_t])
    
    ## Combined QALYs ##
    v_qalys <- c(n_TOTAL_qalys_6mo, n_TOTAL_qalys_1yr, n_TOTAL_qalys_5yr, n_TOTAL_qalys_10yr, n_TOTAL_qalys_life)
    #names(v_qalys) <- c("Total QALYs (1-year)", "Total QALYs (5-year)", "Total QALYs (10-year)", "Total QALYs (Lifetime)")
    
    ## Combined outcomes ##
    df_outcomes <- data.frame(n_TOTAL_costs_6mo, n_TOTAL_costs_1yr, n_TOTAL_costs_5yr, n_TOTAL_costs_10yr, n_TOTAL_costs_life, 
                              n_HEALTH_SECTOR_costs_6mo, n_HEALTH_SECTOR_costs_1yr, n_HEALTH_SECTOR_costs_5yr, n_HEALTH_SECTOR_costs_10yr, n_HEALTH_SECTOR_costs_life,
                              n_CRIMINAL_costs_6mo, n_CRIMINAL_costs_1yr, n_CRIMINAL_costs_5yr, n_CRIMINAL_costs_10yr, n_CRIMINAL_costs_life, 
                              n_TX_costs_6mo, n_TX_costs_1yr, n_TX_costs_5yr, n_TX_costs_10yr, n_TX_costs_life,
                              n_HRU_costs_6mo, n_HRU_costs_1yr, n_HRU_costs_5yr, n_HRU_costs_10yr, n_HRU_costs_life,
                              n_TOTAL_qalys_6mo, n_TOTAL_qalys_1yr, n_TOTAL_qalys_5yr, n_TOTAL_qalys_10yr, n_TOTAL_qalys_life)
    #names(df_outcomes) <- c("Total Costs (1-year)", "Total Costs (5-year)", "Total Costs (10-year)", "Total Costs (Lifetime)", "Health Sector Costs (1-year)", "Health Sector Costs (5-year)", "Health Sector Costs (10-year)", "Health Sector Costs (Lifetime)",
    #                        "Criminal Costs (1-year)", "Criminal Costs (5-year)", "Criminal Costs (10-year)", "Criminal Costs (Lifetime)", "Treatment Costs (1-year)", "Treatment Costs (5-year)", "Treatment Costs (10-year)", "Treatment Costs (Lifetime)",
    #                        "Total QALYs (1-year)", "Total QALYs (5-year)", "Total QALYs (10-year)", "Total QALYs (Lifetime)")

    return(list(a_M_trace = a_M_trace,
                m_M_trace = m_M_trace,
                m_M_agg_trace = m_M_agg_trace,
                m_M_agg_trace_sero = m_M_agg_trace_sero,
                m_TX_costs = m_TX_costs,
                m_HRU_costs = m_HRU_costs,
                m_HIV_costs = m_HIV_costs,
                m_crime_costs = m_crime_costs,
                m_TOTAL_costs_states = m_TOTAL_costs_states,
                m_HEALTH_SECTOR_costs_states = m_HEALTH_SECTOR_costs_states,
                v_TOTAL_costs = v_TOTAL_costs,
                v_HEALTH_SECTOR_costs = v_HEALTH_SECTOR_costs,
                v_CRIMINAL_costs = v_CRIMINAL_costs,
                v_TX_costs = v_TX_costs,
                v_TOTAL_qalys = v_TOTAL_qalys,
                n_TOTAL_costs_6mo = n_TOTAL_costs_6mo,
                n_TOTAL_costs_1yr = n_TOTAL_costs_1yr,
                n_TOTAL_costs_5yr = n_TOTAL_costs_5yr,
                n_TOTAL_costs_10yr = n_TOTAL_costs_10yr,
                n_TOTAL_costs_life = n_TOTAL_costs_life,
                n_HEALTH_SECTOR_costs_6mo = n_HEALTH_SECTOR_costs_6mo,
                n_HEALTH_SECTOR_costs_1yr = n_HEALTH_SECTOR_costs_1yr,
                n_HEALTH_SECTOR_costs_5yr = n_HEALTH_SECTOR_costs_5yr,
                n_HEALTH_SECTOR_costs_10yr = n_HEALTH_SECTOR_costs_10yr,
                n_HEALTH_SECTOR_costs_life = n_HEALTH_SECTOR_costs_life,
                n_CRIMINAL_costs_6mo = n_CRIMINAL_costs_6mo,
                n_CRIMINAL_costs_1yr = n_CRIMINAL_costs_1yr,
                n_CRIMINAL_costs_5yr = n_CRIMINAL_costs_5yr,
                n_CRIMINAL_costs_10yr = n_CRIMINAL_costs_10yr,
                n_CRIMINAL_costs_life = n_CRIMINAL_costs_life,
                n_TX_costs_6mo = n_TX_costs_6mo,
                n_TX_costs_1yr = n_TX_costs_1yr,
                n_TX_costs_5yr = n_TX_costs_5yr,
                n_TX_costs_10yr = n_TX_costs_10yr,
                n_TX_costs_life = n_TX_costs_life,
                n_HRU_costs_6mo = n_HRU_costs_6mo,
                n_HRU_costs_1yr = n_HRU_costs_1yr,
                n_HRU_costs_5yr = n_HRU_costs_5yr,
                n_HRU_costs_10yr = n_HRU_costs_10yr,
                n_HRU_costs_life = n_HRU_costs_life,
                v_costs = v_costs,
                n_TOTAL_qalys_6mo = n_TOTAL_qalys_6mo,
                n_TOTAL_qalys_1yr = n_TOTAL_qalys_1yr,
                n_TOTAL_qalys_5yr = n_TOTAL_qalys_5yr,
                n_TOTAL_qalys_10yr = n_TOTAL_qalys_10yr,
                n_TOTAL_qalys_life = n_TOTAL_qalys_life,
                v_qalys = v_qalys,
                df_outcomes = df_outcomes))
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
  # 6-month
  n_icer_TOTAL_6mo <- (outcomes_int$n_TOTAL_costs_6mo - outcomes_comp$n_TOTAL_costs_6mo)/(outcomes_int$n_TOTAL_qalys_6mo - outcomes_comp$n_TOTAL_qalys_6mo)
  # 1-year
  n_icer_TOTAL_1yr <- (outcomes_int$n_TOTAL_costs_1yr - outcomes_comp$n_TOTAL_costs_1yr)/(outcomes_int$n_TOTAL_qalys_1yr - outcomes_comp$n_TOTAL_qalys_1yr)
  # 5-year
  n_icer_TOTAL_5yr <- (outcomes_int$n_TOTAL_costs_5yr - outcomes_comp$n_TOTAL_costs_5yr)/(outcomes_int$n_TOTAL_qalys_5yr - outcomes_comp$n_TOTAL_qalys_5yr)
  # 10-year
  n_icer_TOTAL_10yr <- (outcomes_int$n_TOTAL_costs_10yr - outcomes_comp$n_TOTAL_costs_10yr)/(outcomes_int$n_TOTAL_qalys_10yr - outcomes_comp$n_TOTAL_qalys_10yr)
  # Lifetime
  n_icer_TOTAL_life <- (outcomes_int$n_TOTAL_costs_life - outcomes_comp$n_TOTAL_costs_life)/(outcomes_int$n_TOTAL_qalys_life - outcomes_comp$n_TOTAL_qalys_life)
  
  ## Health sector ##
  # 1-year
  n_icer_HEALTH_SECTOR_6mo <- (outcomes_int$n_HEALTH_SECTOR_costs_6mo - outcomes_comp$n_HEALTH_SECTOR_costs_6mo)/(outcomes_int$n_TOTAL_qalys_6mo - outcomes_comp$n_TOTAL_qalys_6mo)
  # 1-year
  n_icer_HEALTH_SECTOR_1yr <- (outcomes_int$n_HEALTH_SECTOR_costs_1yr - outcomes_comp$n_HEALTH_SECTOR_costs_1yr)/(outcomes_int$n_TOTAL_qalys_1yr - outcomes_comp$n_TOTAL_qalys_1yr)
  # 5-year
  n_icer_HEALTH_SECTOR_5yr <- (outcomes_int$n_HEALTH_SECTOR_costs_5yr - outcomes_comp$n_HEALTH_SECTOR_costs_5yr)/(outcomes_int$n_TOTAL_qalys_5yr - outcomes_comp$n_TOTAL_qalys_5yr)
  # 10-year
  n_icer_HEALTH_SECTOR_10yr <- (outcomes_int$n_HEALTH_SECTOR_costs_10yr - outcomes_comp$n_HEALTH_SECTOR_costs_10yr)/(outcomes_int$n_TOTAL_qalys_10yr - outcomes_comp$n_TOTAL_qalys_10yr)
  # Lifetime
  n_icer_HEALTH_SECTOR_life <- (outcomes_int$n_HEALTH_SECTOR_costs_life - outcomes_comp$n_HEALTH_SECTOR_costs_life)/(outcomes_int$n_TOTAL_qalys_life - outcomes_comp$n_TOTAL_qalys_life)
  
  ## Incremental outcomes - societal
  # 6-month
  n_inc_costs_TOTAL_6mo <- outcomes_int$n_TOTAL_costs_6mo - outcomes_comp$n_TOTAL_costs_6mo
  n_inc_qalys_TOTAL_6mo <- outcomes_int$n_TOTAL_qalys_6mo - outcomes_comp$n_TOTAL_qalys_6mo
  # 1-year
  n_inc_costs_TOTAL_1yr <- outcomes_int$n_TOTAL_costs_1yr - outcomes_comp$n_TOTAL_costs_1yr
  n_inc_qalys_TOTAL_1yr <- outcomes_int$n_TOTAL_qalys_1yr - outcomes_comp$n_TOTAL_qalys_1yr
  # 5-year
  n_inc_costs_TOTAL_5yr <- outcomes_int$n_TOTAL_costs_5yr - outcomes_comp$n_TOTAL_costs_5yr
  n_inc_qalys_TOTAL_5yr <- outcomes_int$n_TOTAL_qalys_5yr - outcomes_comp$n_TOTAL_qalys_5yr
  # 10-year
  n_inc_costs_TOTAL_10yr <- outcomes_int$n_TOTAL_costs_10yr - outcomes_comp$n_TOTAL_costs_10yr
  n_inc_qalys_TOTAL_10yr <- outcomes_int$n_TOTAL_qalys_10yr - outcomes_comp$n_TOTAL_qalys_10yr
  # Lifetime
  n_inc_costs_TOTAL_life <- outcomes_int$n_TOTAL_costs_life - outcomes_comp$n_TOTAL_costs_life
  n_inc_qalys_TOTAL_life <- outcomes_int$n_TOTAL_qalys_life - outcomes_comp$n_TOTAL_qalys_life
  
  ## Incremental outcomes - Health sector (all costs excluding crime)
  # 6-month
  n_inc_costs_HEALTH_SECTOR_6mo <- outcomes_int$n_HEALTH_SECTOR_costs_6mo - outcomes_comp$n_HEALTH_SECTOR_costs_6mo
  # 1-year
  n_inc_costs_HEALTH_SECTOR_1yr <- outcomes_int$n_HEALTH_SECTOR_costs_1yr - outcomes_comp$n_HEALTH_SECTOR_costs_1yr
  # 5-year
  n_inc_costs_HEALTH_SECTOR_5yr <- outcomes_int$n_HEALTH_SECTOR_costs_5yr - outcomes_comp$n_HEALTH_SECTOR_costs_5yr
  # 10-year
  n_inc_costs_HEALTH_SECTOR_10yr <- outcomes_int$n_HEALTH_SECTOR_costs_10yr - outcomes_comp$n_HEALTH_SECTOR_costs_10yr
  # Lifetime
  n_inc_costs_HEALTH_SECTOR_life <- outcomes_int$n_HEALTH_SECTOR_costs_life - outcomes_comp$n_HEALTH_SECTOR_costs_life
  
  ## Incremental outcomes - OUD treatment costs
  # 6-month
  n_inc_costs_TX_6mo <- outcomes_int$n_TX_costs_6mo - outcomes_comp$n_TX_costs_6mo
  # 1-year
  n_inc_costs_TX_1yr <- outcomes_int$n_TX_costs_1yr - outcomes_comp$n_TX_costs_1yr
  # 5-year
  n_inc_costs_TX_5yr <- outcomes_int$n_TX_costs_5yr - outcomes_comp$n_TX_costs_5yr
  # 10-year
  n_inc_costs_TX_10yr <- outcomes_int$n_TX_costs_10yr - outcomes_comp$n_TX_costs_10yr
  # Lifetime
  n_inc_costs_TX_life <- outcomes_int$n_TX_costs_life - outcomes_comp$n_TX_costs_life
  
  ## Incremental outcomes - HRU costs
  # 6-month
  n_inc_costs_HRU_6mo <- outcomes_int$n_HRU_costs_6mo - outcomes_comp$n_HRU_costs_6mo
  # 1-year
  n_inc_costs_HRU_1yr <- outcomes_int$n_HRU_costs_1yr - outcomes_comp$n_HRU_costs_1yr
  # 5-year
  n_inc_costs_HRU_5yr <- outcomes_int$n_HRU_costs_5yr - outcomes_comp$n_HRU_costs_5yr
  # 10-year
  n_inc_costs_HRU_10yr <- outcomes_int$n_HRU_costs_10yr - outcomes_comp$n_HRU_costs_10yr
  # Lifetime
  n_inc_costs_HRU_life <- outcomes_int$n_HRU_costs_life - outcomes_comp$n_HRU_costs_life
  
  ## Incremental outcomes - crime costs
  # 6-month
  n_inc_costs_CRIMINAL_6mo <- outcomes_int$n_CRIMINAL_costs_6mo - outcomes_comp$n_CRIMINAL_costs_6mo
  # 1-year
  n_inc_costs_CRIMINAL_1yr <- outcomes_int$n_CRIMINAL_costs_1yr - outcomes_comp$n_CRIMINAL_costs_1yr
  # 5-year
  n_inc_costs_CRIMINAL_5yr <- outcomes_int$n_CRIMINAL_costs_5yr - outcomes_comp$n_CRIMINAL_costs_5yr
  # 10-year
  n_inc_costs_CRIMINAL_10yr <- outcomes_int$n_CRIMINAL_costs_10yr - outcomes_comp$n_CRIMINAL_costs_10yr
  # Lifetime
  n_inc_costs_CRIMINAL_life <- outcomes_int$n_CRIMINAL_costs_life - outcomes_comp$n_CRIMINAL_costs_life

  ## Combined dataframes ##
  df_incremental <- data.frame(n_inc_costs_TOTAL_6mo, n_inc_costs_TOTAL_1yr, n_inc_costs_TOTAL_5yr, n_inc_costs_TOTAL_10yr, n_inc_costs_TOTAL_life,
                               n_inc_costs_HEALTH_SECTOR_6mo, n_inc_costs_HEALTH_SECTOR_1yr, n_inc_costs_HEALTH_SECTOR_5yr, n_inc_costs_HEALTH_SECTOR_10yr, n_inc_costs_HEALTH_SECTOR_life,
                               n_inc_qalys_TOTAL_6mo, n_inc_qalys_TOTAL_1yr, n_inc_qalys_TOTAL_5yr, n_inc_qalys_TOTAL_10yr, n_inc_qalys_TOTAL_life,
                               n_inc_costs_TX_6mo, n_inc_costs_TX_1yr, n_inc_costs_TX_5yr, n_inc_costs_TX_10yr, n_inc_costs_TX_life,
                               n_inc_costs_HRU_6mo, n_inc_costs_HRU_1yr, n_inc_costs_HRU_5yr, n_inc_costs_HRU_10yr, n_inc_costs_HRU_life,
                               n_inc_costs_CRIMINAL_6mo, n_inc_costs_CRIMINAL_1yr, n_inc_costs_CRIMINAL_5yr, n_inc_costs_CRIMINAL_10yr, n_inc_costs_CRIMINAL_life)
  #names(df_incremental) <- c("Incremental Costs (1-year)", "Incremental Costs (5-year)", "Incremental Costs (10-year)", "Incremental Costs (Lifetime)", 
  #                           "Incremental Costs (Health Sector 1-year)", "Incremental Costs (Health Sector 5-year)", "Incremental Costs (Health Sector 10-year)", "Incremental Costs (Health Sector Lifetime)",
  #                           "Incremental QALYs (1-year)", "Incremental QALYs (5-year)", "Incremental QALYs (10-year)", "Incremental QALYs (Lifetime)")
  
  df_icer <- data.frame(n_icer_TOTAL_6mo, n_icer_TOTAL_1yr, n_icer_TOTAL_5yr, n_icer_TOTAL_10yr, n_icer_TOTAL_life,
                        n_icer_HEALTH_SECTOR_6mo, n_icer_HEALTH_SECTOR_1yr, n_icer_HEALTH_SECTOR_5yr, n_icer_HEALTH_SECTOR_10yr, n_icer_HEALTH_SECTOR_life)
  #names(df_icer) <- c("ICER (1-year)", "ICER (5-year)", "ICER (10-year)", "ICER (Lifetime)",
  #                    "ICER (Health Sector 1-year)", "ICER (Health Sector 5-year)", "ICER (Health Sector 10-year)", "ICER (Health Sector Lifetime)")
                      #"Incremental Costs (1-year)", "Incremental QALYs (1-year)", "Incremental Costs (5-year)", "Incremental QALYs (5-year)", "Incremental Costs (10-year)", "Incremental QALYs (10-year)", "Incremental Costs (Lifetime)", "Incremental QALYs (Lifetime)",
                      #"Incremental Costs (Health Sector 1-year)", "Incremental QALYs (Health Sector 1-year)", "Incremental Costs (Health Sector 5-year)", "Incremental QALYs (Health Sector 5-year)", "Incremental Costs (Health Sector 10-year)", "Incremental QALYs (Health Sector 10-year)", "Incremental Costs (Health Sector Lifetime)", "Incremental QALYs (Health Sector Lifetime)")
  
  return(list(df_incremental = df_incremental,
              df_icer = df_icer))  
}

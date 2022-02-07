#' Generate PSA dataset of CEA parameters
#'
#' \code{generate_psa_params} generates PSA input dataset by sampling decision 
#' model parameters from their distributions. The sample of the calibrated
#' parameters is a draw from their posterior distribution obtained with the
#' IMIS algorithm.
#' 
#' Parameters that are not sampled in PSA do not need to be defined here, they will
#' default to their original input value for each simulation.
#' @param n_sim Number of PSA samples.
#' @param seed Seed for reproducibility of Monte Carlo sampling.
#' @param n_pop Sample size to determine dirichlet distribution variance.
#' @return 
#' A data frame with \code{n_sim} rows and {n_states} columns of parameters for PSA. 
#' Each row is a parameter set sampled from distributions that characterize 
#' their uncertainty
#' @examples 
#' generate_psa_params()
#' @export
generate_psa_params <- function(n_sim = n_sim, seed = seed, n_pop = n_pop, scenario = scenario,
                                file.death_hr = NULL,
                                file.frailty = NULL,
                                file.weibull = NULL,
                                file.unconditional = NULL,
                                file.overdose = NULL,
                                file.fentanyl = NULL,
                                file.hiv = NULL,
                                file.hcv = NULL,
                                file.costs = NULL,
                                file.crime_costs = NULL,
                                file.qalys = NULL,
                                file.imis_output = NULL){
  
  #Load files with parameter distribution values
  df_death_hr <- read.csv(file = file.death_hr, row.names = 1, header = TRUE) # Mortality hazard ratios
  df_frailty <- read.csv(file = file.frailty, row.names = 1, header = TRUE) # Episode frailty params
  df_weibull <- read.csv(file = file.weibull, row.names = 1, header = TRUE) # Weibull shape and scale
  df_UP <- read.csv(file = file.unconditional, row.names = 1, header = TRUE) # Unconditional transition probs
  df_overdose <- read.csv(file = file.overdose, row.names = 1, header = TRUE) # Overdose params
  df_fentanyl <- read.csv(file = file.fentanyl, row.names = 1, header = TRUE) # Fentanyl params
  df_hiv <- read.csv(file = file.hiv, row.names = 1, header = TRUE) # HIV seroconversion probs
  df_hcv <- read.csv(file = file.hcv, row.names = 1, header = TRUE) # HCV seroconversion probs
  df_costs <- read.csv(file = file.costs, row.names = 1, header = TRUE) # All costs excluding crime
  df_crime_costs <- read.csv(file = file.crime_costs, row.names = 1, header = TRUE) # Crime costs
  df_qalys <- read.csv(file = file.qalys, row.names = 1, header = TRUE) # QALYs
  
  ## Load calibrated parameters
  load(file = file.imis_output)
  #load(file = "outputs/imis_output.RData")
  df_calib_post <- as.data.frame(m_calib_post)
  
  # Set scenario
  scenario <- scenario
  
  # Number of simulations
  n_sim <- n_sim
  if(n_sim != nrow(df_calib_post)){
    warning("Number of PSA simulations and posterior draws not equal")
    
    # Truncate to match simulations
    df_calib_post <- df_calib_post[1:n_sim, ]
  }
  
  # Set seed for random number generator
  seed <- seed
  if (!missing(seed)) 
    set.seed(seed)
  #set_seed <- seed
  
  # Set sample-size for trial-based parameter uncertainty
  n_pop <- n_pop
  #write.csv(n_pop, file = "checks/n_pop.csv", row.names = TRUE)
  # Function to generate lognormal parameter
  location <- function(m = m, s = s){
    log(m^2 / sqrt(s^2 + m^2))
  }
  shape <- function(m = m, s = s){
    shape <- sqrt(log(1 + (s^2 / m^2)))
  }
  
  ########################################
  #### Set up dirichlet random sample ####
  ########################################
  # Need to adjust Dirichlet draws to account for zero transition probabilities from different scenarios
  df_dirichlet_UP = df_UP * n_pop
  write.csv(df_dirichlet_UP, file = "checks/PSA/df_dirichlet_UP.csv", row.names = TRUE)
  m_UP_zeros <- matrix(rep(0, n_sim), nrow = n_sim, ncol = 1) # empty zero matrix to populate impossible transitions
  
  ######################################
  #### Modified Model Specification ####
  ######################################
  if (scenario == "MMS"){
    ## Non-injection ##
    # From BUP
    v_dirichlet_UP_BUP_NI = df_dirichlet_UP["BUP_NI",]
    m_BUP_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_BUP_NI["BUP_NI", "ABS_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "REL_NI"])))
    p_BUP_BUPC_NI = p_BUP_MET_NI = p_BUP_METC_NI = m_UP_zeros[,1]
    p_BUP_ABS_NI = m_BUP_UP_NI[,1]
    p_BUP_REL_NI = 1 - m_BUP_UP_NI[,1]
    # From BUPC
    v_dirichlet_UP_BUPC_NI = df_dirichlet_UP["BUPC_NI",]
    m_BUPC_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_BUPC_NI["BUPC_NI", "ABS_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "REL_NI"])))
    p_BUPC_BUP_NI = p_BUPC_MET_NI = p_BUPC_METC_NI = m_UP_zeros[,1]
    p_BUPC_ABS_NI = m_BUPC_UP_NI[,1]
    p_BUPC_REL_NI = 1 - m_BUPC_UP_NI[,1]
    # From MET
    v_dirichlet_UP_MET_NI = df_dirichlet_UP["MET_NI",]
    m_MET_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_MET_NI["MET_NI", "ABS_NI"], v_dirichlet_UP_MET_NI["MET_NI", "REL_NI"])))
    p_MET_BUPC_NI = p_MET_BUP_NI = p_MET_METC_NI = m_UP_zeros[,1]
    p_MET_ABS_NI = m_MET_UP_NI[,1]
    p_MET_REL_NI = 1 - m_MET_UP_NI[,1]
    # From METC
    v_dirichlet_UP_METC_NI = df_dirichlet_UP["METC_NI",]
    m_METC_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_METC_NI["METC_NI", "ABS_NI"], v_dirichlet_UP_METC_NI["METC_NI", "REL_NI"])))
    p_METC_BUPC_NI = p_METC_BUP_NI = p_METC_MET_NI = m_UP_zeros[,1]
    p_METC_ABS_NI = m_METC_UP_NI[,1]
    p_METC_REL_NI = 1 - m_METC_UP_NI[,1]
    # From ABS (not sampled, all return to relapse)
    p_ABS_BUP_NI = p_ABS_BUPC_NI = p_ABS_MET_NI = p_ABS_METC_NI = p_ABS_REL_NI = m_UP_zeros[,1]
    # From REL
    v_dirichlet_UP_REL_NI = df_dirichlet_UP["REL_NI",]
    m_REL_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_REL_NI["REL_NI", "MET_NI"], v_dirichlet_UP_REL_NI["REL_NI", "METC_NI"], v_dirichlet_UP_REL_NI["REL_NI", "BUP_NI"], v_dirichlet_UP_REL_NI["REL_NI", "BUPC_NI"])))
    p_REL_ABS_NI = m_UP_zeros[,1]
    p_REL_MET_NI  = m_REL_UP_NI[,1]
    p_REL_METC_NI = m_REL_UP_NI[,2]
    p_REL_BUP_NI  = m_REL_UP_NI[,3]
    p_REL_BUPC_NI = m_REL_UP_NI[,4]
    # From OD
    v_dirichlet_UP_OD_NI = df_dirichlet_UP["ODN_NI",]
    m_OD_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_OD_NI["ODN_NI", "MET_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "METC_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "BUP_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "BUPC_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "REL_NI"])))
    p_ODN_ABS_NI = m_UP_zeros[,1]
    p_ODN_MET_NI  = m_OD_UP_NI[,1]
    p_ODN_METC_NI = m_OD_UP_NI[,2]
    p_ODN_BUP_NI  = m_OD_UP_NI[,3]
    p_ODN_BUPC_NI = m_OD_UP_NI[,4]
    p_ODN_REL_NI  = m_OD_UP_NI[,5]
    
    ## Injection ##
    # From BUP
    v_dirichlet_UP_BUP_INJ = df_dirichlet_UP["BUP_INJ",]
    m_BUP_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_BUP_INJ["BUP_INJ", "ABS_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "REL_INJ"])))
    p_BUP_BUPC_INJ = p_BUP_MET_INJ = p_BUP_METC_INJ = m_UP_zeros[,1]
    p_BUP_ABS_INJ = m_BUP_UP_INJ[,1]
    p_BUP_REL_INJ = 1 - m_BUP_UP_INJ[,1]
    # From BUPC
    v_dirichlet_UP_BUPC_INJ = df_dirichlet_UP["BUPC_INJ",]
    m_BUPC_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "ABS_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "REL_INJ"])))
    p_BUPC_BUP_INJ = p_BUPC_MET_INJ = p_BUPC_METC_INJ = m_UP_zeros[,1]
    p_BUPC_ABS_INJ = m_BUPC_UP_INJ[,1]
    p_BUPC_REL_INJ = 1 - m_BUPC_UP_INJ[,1]
    # From MET
    v_dirichlet_UP_MET_INJ = df_dirichlet_UP["MET_INJ",]
    m_MET_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_MET_INJ["MET_INJ", "ABS_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "REL_INJ"])))
    p_MET_BUPC_INJ = p_MET_BUP_INJ = p_MET_METC_INJ = m_UP_zeros[,1]
    p_MET_ABS_INJ = m_MET_UP_INJ[,1]
    p_MET_REL_INJ = 1 - m_MET_UP_INJ[,1]
    # From METC
    v_dirichlet_UP_METC_INJ = df_dirichlet_UP["METC_INJ",]
    m_METC_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_METC_INJ["METC_INJ", "ABS_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "REL_INJ"])))
    p_METC_BUPC_INJ = p_METC_BUP_INJ = p_METC_MET_INJ = m_UP_zeros[,1]
    p_METC_ABS_INJ = m_METC_UP_INJ[,1]
    p_METC_REL_INJ = 1 - m_METC_UP_INJ[,1]
    # From ABS
    p_ABS_BUP_INJ = p_ABS_BUPC_INJ = p_ABS_MET_INJ = p_ABS_METC_INJ = p_ABS_REL_INJ = 0
    # From REL
    v_dirichlet_UP_REL_INJ = df_dirichlet_UP["REL_INJ",]
    m_REL_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_REL_INJ["REL_INJ", "MET_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "METC_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "BUP_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "BUPC_INJ"])))
    p_REL_ABS_INJ = m_UP_zeros[,1]
    p_REL_MET_INJ  = m_REL_UP_INJ[,1]
    p_REL_METC_INJ = m_REL_UP_INJ[,2]
    p_REL_BUP_INJ  = m_REL_UP_INJ[,3]
    p_REL_BUPC_INJ = m_REL_UP_INJ[,4]
    # From OD
    v_dirichlet_UP_OD_INJ = df_dirichlet_UP["ODN_INJ",]
    m_OD_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_OD_INJ["ODN_INJ", "MET_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "METC_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "BUP_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "BUPC_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "REL_INJ"])))
    p_ODN_ABS_INJ = m_UP_zeros[,1]
    p_ODN_MET_INJ  = m_OD_UP_INJ[,1]
    p_ODN_METC_INJ = m_OD_UP_INJ[,2]
    p_ODN_BUP_INJ  = m_OD_UP_INJ[,3]
    p_ODN_BUPC_INJ = m_OD_UP_INJ[,4]
    p_ODN_REL_INJ  = m_OD_UP_INJ[,5]
    
  #############################  
  #### Trial Specification ####
  #############################
  } else if(scenario == "TS_BUP"){
    ### BNX Scenario ###
    ## Non-injection ##
    # From BUP
    v_dirichlet_UP_BUP_NI = df_dirichlet_UP["BUP_NI",]
    m_BUP_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_BUP_NI["BUP_NI", "ABS_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "REL_NI"])))
    p_BUP_BUPC_NI = p_BUP_MET_NI = p_BUP_METC_NI = m_UP_zeros[,1]
    p_BUP_ABS_NI = m_BUP_UP_NI[,1]
    p_BUP_REL_NI = 1 - m_BUP_UP_NI[,1]
    # From BUPC
    p_BUPC_BUP_NI = p_BUPC_MET_NI = p_BUPC_METC_NI = p_BUPC_ABS_NI = p_BUPC_REL_NI = m_UP_zeros[,1]
    # From MET
    p_MET_BUP_NI = p_MET_BUPC_NI = p_MET_METC_NI = p_MET_ABS_NI = p_MET_REL_NI = m_UP_zeros[,1]
    # From METC
    p_METC_BUP_NI = p_METC_BUPC_NI = p_METC_MET_NI = p_METC_ABS_NI = p_METC_REL_NI = m_UP_zeros[,1]
    # From ABS
    v_dirichlet_UP_ABS_NI = df_dirichlet_UP["ABS_NI",]
    m_ABS_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_ABS_NI["ABS_NI", "BUP_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "REL_NI"])))
    p_ABS_BUP_NI = m_ABS_UP_NI[,1]
    p_ABS_REL_NI = 1 - m_ABS_UP_NI[,1]
    p_ABS_BUPC_NI = p_ABS_MET_NI = p_ABS_METC_NI = m_UP_zeros[,1]
    # From REL
    v_dirichlet_UP_REL_NI = df_dirichlet_UP["REL_NI",]
    m_REL_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_REL_NI["REL_NI", "BUP_NI"], v_dirichlet_UP_REL_NI["REL_NI", "ABS_NI"])))
    p_REL_BUP_NI = m_REL_UP_NI[,1]
    p_REL_ABS_NI = 1 - m_REL_UP_NI[,1]
    p_REL_BUP_NI = p_REL_BUPC_NI = p_REL_MET_NI = p_REL_METC_NI = m_UP_zeros[,1]
    # From OD
    v_dirichlet_UP_OD_NI = df_dirichlet_UP["ODN_NI",]
    m_OD_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_OD_NI["ODN_NI", "BUP_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "ABS_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "REL_NI"])))
    p_ODN_BUP_NI = m_OD_UP_NI[,1]
    p_ODN_ABS_NI = m_OD_UP_NI[,2]
    p_ODN_REL_NI = m_OD_UP_NI[,3]
    p_ODN_BUPC_NI = p_ODN_MET_NI = p_ODN_METC_NI = m_UP_zeros[,1]
    
    ## Injection ##
    # From BUP
    v_dirichlet_UP_BUP_INJ = df_dirichlet_UP["BUP_INJ",]
    m_BUP_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_BUP_INJ["BUP_INJ", "ABS_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "REL_INJ"])))
    p_BUP_BUPC_INJ = p_BUP_MET_INJ = p_BUP_METC_INJ = m_UP_zeros[,1]
    p_BUP_ABS_INJ = m_BUP_UP_INJ[,1]
    p_BUP_REL_INJ = 1 - m_BUP_UP_INJ[,1]
    # From BUPC
    p_BUPC_BUP_INJ = p_BUPC_MET_INJ = p_BUPC_METC_INJ = p_BUPC_ABS_INJ = p_BUPC_REL_INJ = m_UP_zeros[,1]
    # From MET
    p_MET_BUP_INJ = p_MET_BUPC_INJ = p_MET_METC_INJ = p_MET_ABS_INJ = p_MET_REL_INJ = m_UP_zeros[,1]
    # From METC
    p_METC_BUP_INJ = p_METC_BUPC_INJ = p_METC_MET_INJ = p_METC_ABS_INJ = p_METC_REL_INJ = m_UP_zeros[,1]
    # From ABS
    v_dirichlet_UP_ABS_INJ = df_dirichlet_UP["ABS_INJ",]
    m_ABS_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_ABS_INJ["ABS_INJ", "BUP_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "REL_INJ"])))
    p_ABS_BUP_INJ = m_ABS_UP_INJ[,1]
    p_ABS_REL_INJ = 1 - m_ABS_UP_INJ[,1]
    p_ABS_BUPC_INJ = p_ABS_MET_INJ = p_ABS_METC_INJ = m_UP_zeros[,1]
    # From REL
    v_dirichlet_UP_REL_INJ = df_dirichlet_UP["REL_INJ",]
    m_REL_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_REL_INJ["REL_INJ", "BUP_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "ABS_INJ"])))
    p_REL_BUP_INJ = m_REL_UP_INJ[,1]
    p_REL_ABS_INJ = 1 - m_REL_UP_INJ[,1]
    p_REL_BUP_INJ = p_REL_BUPC_INJ = p_REL_MET_INJ = p_REL_METC_INJ = m_UP_zeros[,1]
    # From OD
    v_dirichlet_UP_OD_INJ = df_dirichlet_UP["ODN_INJ",]
    m_OD_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_OD_INJ["ODN_INJ", "BUP_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "ABS_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "REL_INJ"])))
    p_ODN_BUP_INJ = m_OD_UP_INJ[,1]
    p_ODN_ABS_INJ = m_OD_UP_INJ[,2]
    p_ODN_REL_INJ = m_OD_UP_INJ[,3]
    p_ODN_BUPC_INJ = p_ODN_MET_INJ = p_ODN_METC_INJ = m_UP_zeros[,1]
    
  } else if (scenario == "TS_MET"){
    ### MET Scenario ###
    # Non-injection
    # From MET
    v_dirichlet_UP_MET_NI = df_dirichlet_UP["MET_NI",]
    m_MET_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_MET_NI["MET_NI", "ABS_NI"], v_dirichlet_UP_MET_NI["MET_NI", "REL_NI"])))
    p_MET_BUPC_NI = p_MET_BUP_NI = p_MET_METC_NI = m_UP_zeros[,1]
    p_MET_ABS_NI = m_MET_UP_NI[,1]
    p_MET_REL_NI = 1 - m_MET_UP_NI[,1]
    # From METC
    p_METC_MET_NI = p_METC_BUP_NI = p_METC_BUPC_NI = p_METC_ABS_NI = p_METC_REL_NI = m_UP_zeros[,1]
    # From BUP
    p_BUP_BUPC_NI = p_BUP_MET_NI = p_BUP_METC_NI = p_BUP_ABS_NI = p_BUP_REL_NI = m_UP_zeros[,1]
    # From BUPC
    p_BUPC_BUP_NI = p_BUPC_MET_NI = p_BUPC_METC_NI = p_BUPC_ABS_NI = p_BUPC_REL_NI = m_UP_zeros[,1]
    # From ABS
    v_dirichlet_UP_ABS_NI = df_dirichlet_UP["ABS_NI",]
    m_ABS_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_ABS_NI["ABS_NI", "MET_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "REL_NI"])))
    p_ABS_MET_NI = m_ABS_UP_NI[,1]
    p_ABS_REL_NI = 1 - m_ABS_UP_NI[,1]
    p_ABS_BUP_NI = p_ABS_BUPC_NI = p_ABS_METC_NI = m_UP_zeros[,1]
    # From REL
    v_dirichlet_UP_REL_NI = df_dirichlet_UP["REL_NI",]
    m_REL_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_REL_NI["REL_NI", "MET_NI"], v_dirichlet_UP_REL_NI["REL_NI", "ABS_NI"])))
    p_REL_MET_NI = m_REL_UP_NI[,1]
    p_REL_ABS_NI = 1 - m_REL_UP_NI[,1]
    p_REL_BUP_NI = p_REL_BUPC_NI = p_REL_METC_NI = m_UP_zeros[,1]
    # From OD
    v_dirichlet_UP_OD_NI = df_dirichlet_UP["ODN_NI",]
    m_OD_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_OD_NI["ODN_NI", "MET_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "ABS_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "REL_NI"])))
    p_ODN_MET_NI = m_OD_UP_NI[,1]
    p_ODN_ABS_NI = m_OD_UP_NI[,2]
    p_ODN_REL_NI = m_OD_UP_NI[,3]
    p_ODN_BUP_NI = p_ODN_BUPC_NI = p_ODN_METC_NI = m_UP_zeros[,1]
    
    # Injection
    # From MET
    v_dirichlet_UP_MET_INJ = df_dirichlet_UP["MET_INJ",]
    m_MET_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_MET_INJ["MET_INJ", "ABS_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "REL_INJ"])))
    p_MET_BUPC_INJ = p_MET_BUP_INJ = p_MET_METC_INJ = m_UP_zeros[,1]
    p_MET_ABS_INJ = m_MET_UP_INJ[,1]
    p_MET_REL_INJ = 1 - m_MET_UP_INJ[,1]
    # From METC
    p_METC_MET_INJ = p_METC_BUP_INJ = p_METC_BUPC_INJ = p_METC_ABS_INJ = p_METC_REL_INJ = m_UP_zeros[,1]
    # From BUP
    p_BUP_BUPC_INJ = p_BUP_MET_INJ = p_BUP_METC_INJ = p_BUP_ABS_INJ = p_BUP_REL_INJ = m_UP_zeros[,1]
    # From BUPC
    p_BUPC_BUP_INJ = p_BUPC_MET_INJ = p_BUPC_METC_INJ = p_BUPC_ABS_INJ = p_BUPC_REL_INJ = m_UP_zeros[,1]
    # From ABS
    v_dirichlet_UP_ABS_INJ = df_dirichlet_UP["ABS_INJ",]
    m_ABS_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_ABS_INJ["ABS_INJ", "MET_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "REL_INJ"])))
    p_ABS_MET_INJ = m_ABS_UP_INJ[,1]
    p_ABS_REL_INJ = 1 - m_ABS_UP_INJ[,1]
    p_ABS_BUP_INJ = p_ABS_BUPC_INJ = p_ABS_METC_INJ = m_UP_zeros[,1]
    # From REL
    v_dirichlet_UP_REL_INJ = df_dirichlet_UP["REL_INJ",]
    m_REL_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_REL_INJ["REL_INJ", "MET_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "ABS_INJ"])))
    p_REL_MET_INJ = m_REL_UP_INJ[,1]
    p_REL_ABS_INJ = 1 - m_REL_UP_INJ[,1]
    p_REL_BUP_INJ = p_REL_BUPC_INJ = p_REL_METC_INJ = m_UP_zeros[,1]
    # From OD
    v_dirichlet_UP_OD_INJ = df_dirichlet_UP["ODN_INJ",]
    m_OD_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_OD_INJ["ODN_INJ", "MET_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "ABS_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "REL_INJ"])))
    p_ODN_MET_INJ = m_OD_UP_INJ[,1]
    p_ODN_ABS_INJ = m_OD_UP_INJ[,2]
    p_ODN_REL_INJ = m_OD_UP_INJ[,3]
    p_ODN_BUP_INJ = p_ODN_BUPC_INJ = p_ODN_METC_INJ = m_UP_zeros[,1]
    
  ################################  
  #### Original Specification ####
  ################################ 
    
  } else if(scenario == "OS"){
    ## Non-injection ##
    # From BUP
    v_dirichlet_UP_BUP_NI = df_dirichlet_UP["BUP_NI",]
    m_BUP_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_BUP_NI["BUP_NI", "ABS_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "REL_NI"])))
    p_BUP_BUPC_NI = p_BUP_MET_NI = p_BUP_METC_NI = m_UP_zeros[,1]
    p_BUP_ABS_NI = m_BUP_UP_NI[,1]
    p_BUP_REL_NI = 1 - m_BUP_UP_NI[,1]
    # From BUPC
    p_BUPC_BUP_NI = p_BUPC_MET_NI = p_BUPC_METC_NI = p_BUPC_ABS_NI = p_BUPC_REL_NI = m_UP_zeros[,1]
    # From MET
    v_dirichlet_UP_MET_NI = df_dirichlet_UP["MET_NI",]
    m_MET_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_MET_NI["MET_NI", "ABS_NI"], v_dirichlet_UP_MET_NI["MET_NI", "REL_NI"])))
    p_MET_BUPC_NI = p_MET_BUP_NI = p_MET_METC_NI = m_UP_zeros[,1]
    p_MET_ABS_NI = m_MET_UP_NI[,1]
    p_MET_REL_NI = 1 - m_MET_UP_NI[,1]
    # From METC
    p_METC_BUPC_NI = p_METC_BUP_NI = p_METC_MET_NI = p_METC_ABS_NI = p_METC_REL_NI = m_UP_zeros[,1]
    # From ABS (not sampled, all return to relapse)
    p_ABS_BUP_NI = p_ABS_BUPC_NI = p_ABS_MET_NI = p_ABS_METC_NI = p_ABS_REL_NI = m_UP_zeros[,1]
    # From REL
    v_dirichlet_UP_REL_NI = df_dirichlet_UP["REL_NI",]
    m_REL_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_REL_NI["REL_NI", "MET_NI"], v_dirichlet_UP_REL_NI["REL_NI", "BUP_NI"])))
    p_REL_BUPC_NI = p_REL_METC_NI = p_REL_ABS_NI = m_UP_zeros[,1]
    p_REL_MET_NI  = m_REL_UP_NI[,1]
    p_REL_BUP_NI  = 1 - m_REL_UP_NI[,1]
    # From OD
    v_dirichlet_UP_OD_NI = df_dirichlet_UP["ODN_NI",]
    m_OD_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_OD_NI["ODN_NI", "MET_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "BUP_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "REL_NI"])))
    p_ODN_BUPC_NI = p_ODN_METC_NI = p_ODN_ABS_NI = m_UP_zeros[,1]
    p_ODN_MET_NI  = m_OD_UP_NI[,1]
    p_ODN_BUP_NI  = m_OD_UP_NI[,2]
    p_ODN_REL_NI  = m_OD_UP_NI[,3]
    
    ## Injection ##
    # From BUP
    v_dirichlet_UP_BUP_INJ = df_dirichlet_UP["BUP_INJ",]
    m_BUP_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_BUP_INJ["BUP_INJ", "ABS_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "REL_INJ"])))
    p_BUP_BUPC_INJ = p_BUP_MET_INJ = p_BUP_METC_INJ = m_UP_zeros[,1]
    p_BUP_ABS_INJ = m_BUP_UP_INJ[,1]
    p_BUP_REL_INJ = 1 - m_BUP_UP_INJ[,1]
    # From BUPC
    p_BUPC_BUP_INJ = p_BUPC_MET_INJ = p_BUPC_METC_INJ = p_BUPC_ABS_INJ = p_BUPC_REL_INJ = m_UP_zeros[,1]
    # From MET
    v_dirichlet_UP_MET_INJ = df_dirichlet_UP["MET_INJ",]
    m_MET_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_MET_INJ["MET_INJ", "ABS_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "REL_INJ"])))
    p_MET_BUPC_INJ = p_MET_BUP_INJ = p_MET_METC_INJ = m_UP_zeros[,1]
    p_MET_ABS_INJ = m_MET_UP_INJ[,1]
    p_MET_REL_INJ = 1 - m_MET_UP_INJ[,1]
    # From METC
    p_METC_BUPC_INJ = p_METC_BUP_INJ = p_METC_MET_INJ = p_METC_ABS_INJ = p_METC_REL_INJ = m_UP_zeros[,1]
    # From ABS (not sampled, all return to relapse)
    p_ABS_BUP_INJ = p_ABS_BUPC_INJ = p_ABS_MET_INJ = p_ABS_METC_INJ = p_ABS_REL_INJ = m_UP_zeros[,1]
    # From REL
    v_dirichlet_UP_REL_INJ = df_dirichlet_UP["REL_INJ",]
    m_REL_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_REL_INJ["REL_INJ", "MET_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "BUP_INJ"])))
    p_REL_BUPC_INJ = p_REL_METC_INJ = p_REL_ABS_INJ = m_UP_zeros[,1]
    p_REL_MET_INJ  = m_REL_UP_INJ[,1]
    p_REL_BUP_INJ  = 1 - m_REL_UP_INJ[,1]
    # From OD
    v_dirichlet_UP_OD_INJ = df_dirichlet_UP["ODN_INJ",]
    m_OD_UP_INJ = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_OD_INJ["ODN_INJ", "MET_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "BUP_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "REL_INJ"])))
    p_ODN_BUPC_INJ = p_ODN_METC_INJ = p_ODN_ABS_INJ = m_UP_zeros[,1]
    p_ODN_MET_INJ  = m_OD_UP_INJ[,1]
    p_ODN_BUP_INJ  = m_OD_UP_INJ[,2]
    p_ODN_REL_INJ  = m_OD_UP_INJ[,3]
  } else{
    print("No scenario selected")
  }
  
  #### Hazard ratios for death probability ####
  # ***NEW ESTIMATES*** #
  # Non-injection
  hr_TX  = rlnorm(n_sim, location(m = df_death_hr["pe", "TX"], s = df_death_hr["sd", "TX"]), shape(m = df_death_hr["pe", "TX"], s = df_death_hr["sd", "TX"]))
  hr_REL = rlnorm(n_sim, location(m = df_death_hr["pe", "REL"], s = df_death_hr["sd", "REL"]), shape(m = df_death_hr["pe", "REL"], s = df_death_hr["sd", "REL"]))
  hr_ABS = rlnorm(n_sim, location(m = df_death_hr["pe", "ABS"], s = df_death_hr["sd", "ABS"]), shape(m = df_death_hr["pe", "ABS"], s = df_death_hr["sd", "ABS"]))
  hr_HIV = rlnorm(n_sim, location(m = df_death_hr["pe", "HIV"], s = df_death_hr["sd", "HIV"]), shape(m = df_death_hr["pe", "HIV"], s = df_death_hr["sd", "HIV"]))
  hr_HCV = rlnorm(n_sim, location(m = df_death_hr["pe", "HCV"], s = df_death_hr["sd", "HCV"]), shape(m = df_death_hr["pe", "HCV"], s = df_death_hr["sd", "HCV"]))
  hr_COI = rlnorm(n_sim, location(m = df_death_hr["pe", "COI"], s = df_death_hr["sd", "COI"]), shape(m = df_death_hr["pe", "COI"], s = df_death_hr["sd", "COI"]))

  # Set other parameters that are structurally equivalent for identical draws (e.g. probability of HIV seroconversion among all NI states)
  # HIV
  p_HIV_NI = rbeta(n_sim, shape1 = df_hiv["shape1", "HIV_NI"], shape2 = df_hiv["shape2", "HIV_NI"])
  p_HIV_TX_INJ  = rbeta(n_sim, shape1 = df_hiv["shape1", "HIV_TX_INJ"], shape2 = df_hiv["shape2", "HIV_TX_INJ"])
  p_HIV_TXC_INJ = rbeta(n_sim, shape1 = df_hiv["shape1", "HIV_TXC_INJ"], shape2 = df_hiv["shape2", "HIV_TXC_INJ"])
  p_HIV_REL_INJ  = rbeta(n_sim, shape1 = df_hiv["shape1", "HIV_REL_INJ"], shape2 = df_hiv["shape2", "HIV_REL_INJ"])
  
  # HCV
  p_HCV_NI = rbeta(n_sim, shape1 = df_hcv["shape1", "HCV_NI"], shape2 = df_hcv["shape2", "HCV_NI"])
  p_HCV_TX_INJ  = rbeta(n_sim, shape1 = df_hcv["shape1", "HCV_TX_INJ"], shape2 = df_hcv["shape2", "HCV_TX_INJ"])
  p_HCV_TXC_INJ = rbeta(n_sim, shape1 = df_hcv["shape1", "HCV_TXC_INJ"], shape2 = df_hcv["shape2", "HCV_TXC_INJ"])
  p_HCV_REL_INJ  = rbeta(n_sim, shape1 = df_hcv["shape1", "HCV_REL_INJ"], shape2 = df_hcv["shape2", "HCV_REL_INJ"])
  
  # HRU Costs
  c_BUP_HRU = rgamma(n_sim, shape = df_costs["shape", "BUP_NI_HRU"], scale = df_costs["scale", "BUP_NI_HRU"])
  c_BUPC_HRU = rgamma(n_sim, shape = df_costs["shape", "BUPC_NI_HRU"], scale = df_costs["scale", "BUPC_NI_HRU"])
  c_MET_HRU = rgamma(n_sim, shape = df_costs["shape", "MET_NI_HRU"], scale = df_costs["scale", "MET_NI_HRU"])
  c_METC_HRU = rgamma(n_sim, shape = df_costs["shape", "METC_NI_HRU"], scale = df_costs["scale", "METC_NI_HRU"])
  c_ABS_HRU = rgamma(n_sim, shape = df_costs["shape", "ABS_NI_HRU"], scale = df_costs["scale", "ABS_NI_HRU"]) 
  c_REL_HRU = rgamma(n_sim, shape = df_costs["shape", "REL_NI_HRU"], scale = df_costs["scale", "REL_NI_HRU"]) 
  c_ODN_HRU  = rgamma(n_sim, shape = df_costs["shape", "ODN_NI_HRU"] , scale = df_costs["scale", "ODN_NI_HRU"])
  c_ODF_HRU  = rgamma(n_sim, shape = df_costs["shape", "ODF_NI_HRU"] , scale = df_costs["scale", "ODF_NI_HRU"])
  
  # Crime Costs
  c_BUP_crime = rgamma(n_sim, shape = df_crime_costs["shape", "BUP"], scale = df_crime_costs["scale", "BUP"])
  c_BUPC_crime = rgamma(n_sim, shape = df_crime_costs["shape", "BUPC"], scale = df_crime_costs["scale", "BUPC"])
  c_MET_crime = rgamma(n_sim, shape = df_crime_costs["shape", "MET"], scale = df_crime_costs["scale", "MET"])
  c_METC_crime = rgamma(n_sim, shape = df_crime_costs["shape", "METC"], scale = df_crime_costs["scale", "METC"])
  c_REL_crime = rgamma(n_sim, shape = df_crime_costs["shape", "REL"], scale = df_crime_costs["scale", "REL"])
  c_ODN_crime = rgamma(n_sim, shape = df_crime_costs["shape", "REL"], scale = df_crime_costs["scale", "REL"])
  
  # QALYs
  u_BUP_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "BUP"], sd = df_qalys["sd", "BUP"], b = 1)
  u_BUPC_NEG = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "BUPC"], sd = df_qalys["sd", "BUPC"], b = 1)
  u_MET_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "MET"], sd = df_qalys["sd", "MET"], b = 1)
  u_METC_NEG = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "METC"], sd = df_qalys["sd", "METC"], b = 1)
  u_REL_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "REL"], sd = df_qalys["sd", "REL"], b = 1)
  u_ODN_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "ODN"], sd = df_qalys["sd", "ODN"], b = 1)
  u_ABS_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "ABS"], sd = df_qalys["sd", "ABS"], b = 1)
  
  u_HIV_mult = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "HIV_mult"], sd = df_qalys["sd", "HIV_mult"], b = 1)
  u_HCV_mult = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "HCV_mult"], sd = df_qalys["sd", "HCV_mult"], b = 1)
  u_COI_mult = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "COI_mult"], sd = df_qalys["sd", "COI_mult"], b = 1)
  
  df_psa_params <- data.frame(
    ### Calibrated parameters
    df_calib_post, # Matrix of calibration parameters drawn from posterior distribution
    
    # Hazard ratios for death probability
    # Log-normal distribution
    # Non-injection
    hr_BUP_NI  = hr_TX,
    hr_BUPC_NI = hr_TX,
    hr_MET_NI  = hr_TX,
    hr_METC_NI = hr_TX,
    hr_REL_NI  = hr_REL,
    hr_ODN_NI  = hr_REL,
    hr_HIV_NI  = hr_HIV,
    hr_HCV_NI  = hr_HCV,
    hr_COI_NI  = hr_COI,
    
    # Injection
    hr_BUP_INJ  = hr_TX,
    hr_BUPC_INJ = hr_TX,
    hr_MET_INJ  = hr_TX,
    hr_METC_INJ = hr_TX,
    hr_REL_INJ  = hr_REL,
    hr_ODN_INJ  = hr_REL,
    hr_HIV_INJ  = hr_HIV,
    hr_HCV_INJ  = hr_HCV,
    hr_COI_INJ  = hr_COI,

    # Frailty terms for successive episodes
    # Log-normal distribution (mean, sd)
    # Episodes
    p_frailty_BUP_2 = rlnorm(n_sim, location(m = df_frailty["pe", "BUP_2"], s = df_frailty["sd", "BUP_2"]), shape(m = df_frailty["pe", "BUP_2"], s = df_frailty["sd", "BUP_2"])),
    p_frailty_BUP_3 = rlnorm(n_sim, location(m = df_frailty["pe", "BUP_3"], s = df_frailty["sd", "BUP_3"]), shape(m = df_frailty["pe", "BUP_3"], s = df_frailty["sd", "BUP_3"])),

    p_frailty_MET_2 = rlnorm(n_sim, location(m = df_frailty["pe", "MET_2"], s = df_frailty["sd", "MET_2"]), shape(m = df_frailty["pe", "MET_2"], s = df_frailty["sd", "MET_2"])),
    p_frailty_MET_3 = rlnorm(n_sim, location(m = df_frailty["pe", "MET_3"], s = df_frailty["sd", "MET_3"]), shape(m = df_frailty["pe", "MET_3"], s = df_frailty["sd", "MET_3"])),

    p_frailty_REL_2 = rlnorm(n_sim, location(m = df_frailty["pe", "REL_2"], s = df_frailty["sd", "REL_2"]), shape(m = df_frailty["pe", "REL_2"], s = df_frailty["sd", "REL_2"])),
    p_frailty_REL_3 = rlnorm(n_sim, location(m = df_frailty["pe", "REL_3"], s = df_frailty["sd", "REL_3"]), shape(m = df_frailty["pe", "REL_3"], s = df_frailty["sd", "REL_3"])),

    p_frailty_ABS_2 = rlnorm(n_sim, location(m = df_frailty["pe", "ABS_2"], s = df_frailty["sd", "ABS_2"]), shape(m = df_frailty["pe", "ABS_2"], s = df_frailty["sd", "ABS_2"])),
    p_frailty_ABS_3 = rlnorm(n_sim, location(m = df_frailty["pe", "ABS_3"], s = df_frailty["sd", "ABS_3"]), shape(m = df_frailty["pe", "ABS_3"], s = df_frailty["sd", "ABS_3"])),
    
    # Injection vs. non-injection
    p_frailty_BUP_INJ = rlnorm(n_sim, location(m = df_frailty["pe", "BUP_INJ"], s = df_frailty["sd", "BUP_INJ"]), shape(m = df_frailty["pe", "BUP_INJ"], s = df_frailty["sd", "BUP_INJ"])),
    p_frailty_MET_INJ = rlnorm(n_sim, location(m = df_frailty["pe", "MET_INJ"], s = df_frailty["sd", "MET_INJ"]), shape(m = df_frailty["pe", "MET_INJ"], s = df_frailty["sd", "MET_INJ"])),
    p_frailty_REL_INJ = rlnorm(n_sim, location(m = df_frailty["pe", "REL_INJ"], s = df_frailty["sd", "REL_INJ"]), shape(m = df_frailty["pe", "REL_INJ"], s = df_frailty["sd", "REL_INJ"])),
    p_frailty_ABS_INJ = rlnorm(n_sim, location(m = df_frailty["pe", "ABS_INJ"], s = df_frailty["sd", "ABS_INJ"]), shape(m = df_frailty["pe", "ABS_INJ"], s = df_frailty["sd", "ABS_INJ"])),
    
    # Concurrent opioid use
    p_frailty_BUPC = rlnorm(n_sim, location(m = df_frailty["pe", "BUPC"], s = df_frailty["sd", "BUPC"]), shape(m = df_frailty["pe", "BUPC"], s = df_frailty["sd", "BUPC"])),
    p_frailty_METC = rlnorm(n_sim, location(m = df_frailty["pe", "METC"], s = df_frailty["sd", "METC"]), shape(m = df_frailty["pe", "METC"], s = df_frailty["sd", "METC"])),
    
    # Weibull
    # Shape
    # EP1
    p_weibull_shape_BUP = rnorm(n_sim, mean = df_weibull["pe", "BUP_shape_1"], sd = df_weibull["sd", "BUP_shape_1"]),
    p_weibull_shape_MET = rnorm(n_sim, mean = df_weibull["pe", "MET_shape_1"], sd = df_weibull["sd", "MET_shape_1"]),
    p_weibull_shape_REL = rnorm(n_sim, mean = df_weibull["pe", "REL_shape_1"], sd = df_weibull["sd", "REL_shape_1"]),
    p_weibull_shape_ABS = rnorm(n_sim, mean = df_weibull["pe", "ABS_shape_1"], sd = df_weibull["sd", "ABS_shape_1"]),
    
    # scale
    # EP1
    p_weibull_scale_BUP = rnorm(n_sim, mean = df_weibull["pe", "BUP_scale_1"], sd = df_weibull["sd", "BUP_scale_1"]),
    p_weibull_scale_MET = rnorm(n_sim, mean = df_weibull["pe", "MET_scale_1"], sd = df_weibull["sd", "MET_scale_1"]),
    p_weibull_scale_REL = rnorm(n_sim, mean = df_weibull["pe", "REL_scale_1"], sd = df_weibull["sd", "REL_scale_1"]),
    p_weibull_scale_ABS = rnorm(n_sim, mean = df_weibull["pe", "ABS_scale_1"], sd = df_weibull["sd", "ABS_scale_1"]),
    
    ### Transition probabilities conditional on leaving (derived from Dirichlet above)
    ## Non-Injection ##
    # From BUP
    p_BUP_BUPC_NI, p_BUP_MET_NI, p_BUP_METC_NI, p_BUP_ABS_NI, p_BUP_REL_NI,
    # From BUPC
    p_BUPC_BUP_NI, p_BUPC_MET_NI, p_BUPC_METC_NI, p_BUPC_ABS_NI, p_BUPC_REL_NI,
    # From MET
    p_MET_BUP_NI, p_MET_BUPC_NI, p_MET_METC_NI, p_MET_ABS_NI, p_MET_REL_NI,
    # From METC
    p_METC_BUP_NI, p_METC_BUPC_NI, p_METC_MET_NI, p_METC_ABS_NI, p_METC_REL_NI,
    # From ABS
    p_ABS_BUP_NI, p_ABS_BUPC_NI, p_ABS_MET_NI, p_ABS_METC_NI, p_ABS_REL_NI,
    # From REL
    p_REL_BUP_NI, p_REL_BUPC_NI, p_REL_MET_NI, p_REL_METC_NI, p_REL_ABS_NI,
    # From OD
    p_ODN_BUP_NI, p_ODN_BUPC_NI, p_ODN_MET_NI, p_ODN_METC_NI, p_ODN_ABS_NI, p_ODN_REL_NI,

    ## Injection ##
    # From BUP
    p_BUP_BUPC_INJ, p_BUP_MET_INJ, p_BUP_METC_INJ, p_BUP_ABS_INJ, p_BUP_REL_INJ,
    # From BUPC
    p_BUPC_BUP_INJ, p_BUPC_MET_INJ, p_BUPC_METC_INJ, p_BUPC_ABS_INJ, p_BUPC_REL_INJ,
    # From MET
    p_MET_BUP_INJ, p_MET_BUPC_INJ, p_MET_METC_INJ, p_MET_ABS_INJ, p_MET_REL_INJ,
    # From METC
    p_METC_BUP_INJ, p_METC_BUPC_INJ, p_METC_MET_INJ, p_METC_ABS_INJ, p_METC_REL_INJ,
    # From ABS
    p_ABS_BUP_INJ, p_ABS_BUPC_INJ, p_ABS_MET_INJ, p_ABS_METC_INJ, p_ABS_REL_INJ,
    # From REL
    p_REL_BUP_INJ, p_REL_BUPC_INJ, p_REL_MET_INJ, p_REL_METC_INJ, p_REL_ABS_INJ,
    # From OD
    p_ODN_BUP_INJ, p_ODN_BUPC_INJ, p_ODN_MET_INJ, p_ODN_METC_INJ, p_ODN_ABS_INJ, p_ODN_REL_INJ,

    ### Overdose ###
    # Overdose transition multipliers
    # First month
    # BUP
    n_BUP_OD_mult = rgamma(n_sim, shape = df_overdose["shape", "BUP_OD_mult"], scale = df_overdose["scale", "BUP_OD_mult"]),
    #n_BUP_OD_mult = df_overdose["pe", "BUP_OD_mult"],
    #n_BUP_OD_mult_shape = df_overdose["shape", "BUP_OD_mult"],
    #n_BUP_OD_mult_scale = df_overdose["scale", "BUP_OD_mult"],
    
    # MET
    n_MET_OD_mult = rgamma(n_sim, shape = df_overdose["shape", "MET_OD_mult"], scale = df_overdose["scale", "MET_OD_mult"]),
    #n_MET_OD_mult = df_overdose["pe", "MET_OD_mult"],
    #n_MET_OD_mult_shape = df_overdose["shape", "MET_OD_mult"],
    #n_MET_OD_mult_scale = df_overdose["scale", "MET_OD_mult"],

    # Relapse
    n_REL_OD_mult = rgamma(n_sim, shape = df_overdose["shape", "REL_OD_mult"], scale = df_overdose["scale", "REL_OD_mult"]),
    #n_REL_OD_mult  = df_overdose["pe", "REL_OD_mult"],
    #n_REL_OD_mult_shape  = df_overdose["shape", "REL_OD_mult"],
    #n_REL_OD_mult_scale  = df_overdose["scale", "REL_OD_mult"],
    
    # Abstinence
    #n_ABS_OD_mult  = df_overdose["pe", "ABS_OD_mult"],
    #n_ABS_OD_mult_shape  = df_overdose["shape", "ABS_OD_mult"],
    #n_ABS_OD_mult_scale  = df_overdose["scale", "ABS_OD_mult"],
    
    # Injection (vs. non-injection)
    #n_INJ_OD_mult = df_overdose["pe", "INJ_OD_mult"],
    #n_INJ_OD_mult_shape = df_overdose["shape", "INJ_OD_mult"],
    #n_INJ_OD_mult_scale = df_overdose["scale", "INJ_OD_mult"],

    #n_fent_OD = rgamma(n_sim, shape = df_overdose["shape", "fent_OD_rate"], scale = df_overdose["scale", "fent_OD_rate"]),
    #p_fent_exp = rbeta(n_sim, shape1 = df_overdose["shape1", "fent_exp_prob"], shape2 = df_overdose["shape2", "fent_exp_prob"]),
    # Fentanyl
    p_ni_fent_reduction = rbeta(n_sim, shape1 = df_overdose["shape1", "ni_fent_reduction"], shape2 = df_overdose["shape2", "ni_fent_reduction"]),
    #p_fent_exp_2020 = df_fentanyl["2020", "CAN"], # add distribution parameters (uniform between range of multiple fentanyl % estimates)
    p_fent_exp_2020 = runif(n_sim, min = df_fentanyl["low", "BC"], max = df_fentanyl["high", "BC"]),
    
    # Naloxone (OD reversal)
    p_witness = rbeta(n_sim, shape1 = df_overdose["shape1", "witness_prob"], shape2 = df_overdose["shape2", "witness_prob"]),
    p_attended = rbeta(n_sim, shape1 = df_overdose["shape1", "attended_prob"], shape2 = df_overdose["shape2", "attended_prob"]),
    p_NX_used = rbeta(n_sim, shape1 = df_overdose["shape1", "NX_prob"], shape2 = df_overdose["shape2", "NX_prob"]),
    p_NX_success = rbeta(n_sim, shape1 = df_overdose["shape1", "NX_success_prob"], shape2 = df_overdose["shape2", "NX_success_prob"]),

    ### HIV seroconversion ###
    # Ensure that seed is set and produces identical draws for parameters that are set to be equal by assumption (e.g. all non-injection HIV seroconversion)
    # From negative
    # Non-injection
    #p_HIV_NI = rbeta(n_sim, shape1 = df_hiv["shape1", "HIV_NI"], shape2 = df_hiv["shape2", "HIV_NI"]),
    p_HIV_BUP_NI  = p_HIV_NI,
    p_HIV_BUPC_NI = p_HIV_NI,
    p_HIV_MET_NI  = p_HIV_NI,
    p_HIV_METC_NI = p_HIV_NI,
    p_HIV_REL_NI  = p_HIV_NI,
    p_HIV_ODN_NI  = p_HIV_NI,
    p_HIV_ABS_NI  = p_HIV_NI,
    # Injection
    p_HIV_BUP_INJ  = p_HIV_TX_INJ,
    p_HIV_BUPC_INJ = p_HIV_TXC_INJ,
    p_HIV_MET_INJ  = p_HIV_TX_INJ,
    p_HIV_METC_INJ = p_HIV_TXC_INJ,
    p_HIV_REL_INJ  = p_HIV_REL_INJ,
    p_HIV_ODN_INJ  = p_HIV_REL_INJ,
    p_HIV_ABS_INJ  = p_HIV_NI,
    
    # Co-infection conditional on HCV
    # Same as HIV from negative by assumption
    # Non-injection
    p_HCV_HIV_BUP_NI  = p_HIV_NI,
    p_HCV_HIV_BUPC_NI = p_HIV_NI,
    p_HCV_HIV_MET_NI  = p_HIV_NI,
    p_HCV_HIV_METC_NI = p_HIV_NI,
    p_HCV_HIV_REL_NI  = p_HIV_NI,
    p_HCV_HIV_ODN_NI  = p_HIV_NI,
    p_HCV_HIV_ABS_NI  = p_HIV_NI,
    # Injection
    p_HCV_HIV_BUP_INJ  = p_HIV_TX_INJ,
    p_HCV_HIV_BUPC_INJ = p_HIV_TXC_INJ,
    p_HCV_HIV_MET_INJ  = p_HIV_TX_INJ,
    p_HCV_HIV_METC_INJ = p_HIV_TXC_INJ,
    p_HCV_HIV_REL_INJ  = p_HIV_REL_INJ,
    p_HCV_HIV_ODN_INJ  = p_HIV_REL_INJ,
    p_HCV_HIV_ABS_INJ  = p_HIV_NI,
    
    ### HCV seroconversion ###
    # HCV Seroconversion
    # From negative
    # Non-injection
    #p_HCV_NI = rbeta(n_sim, shape1 = df_hcv["shape1", "HCV_NI"], shape2 = df_hcv["shape2", "HCV_NI"]),
    p_HCV_BUP_NI  = p_HCV_NI,
    p_HCV_BUPC_NI = p_HCV_NI,
    p_HCV_MET_NI  = p_HCV_NI,
    p_HCV_METC_NI = p_HCV_NI,
    p_HCV_REL_NI  = p_HCV_NI,
    p_HCV_ODN_NI  = p_HCV_NI,
    p_HCV_ABS_NI  = p_HCV_NI,
    # Injection
    p_HCV_BUP_INJ  = p_HCV_TX_INJ,
    p_HCV_BUPC_INJ = p_HCV_TXC_INJ,
    p_HCV_MET_INJ  = p_HCV_TX_INJ,
    p_HCV_METC_INJ = p_HCV_TXC_INJ,
    p_HCV_REL_INJ  = p_HCV_REL_INJ,
    p_HCV_ODN_INJ  = p_HCV_REL_INJ,
    p_HCV_ABS_INJ  = p_HCV_NI,
    
    # Co-infection conditional on HIV
    # Same as HCV from negative by assumption
    # Non-injection
    p_HIV_HCV_BUP_NI  = p_HCV_NI,
    p_HIV_HCV_BUPC_NI = p_HCV_NI,
    p_HIV_HCV_MET_NI  = p_HCV_NI,
    p_HIV_HCV_METC_NI = p_HCV_NI,
    p_HIV_HCV_REL_NI  = p_HCV_NI,
    p_HIV_HCV_ODN_NI  = p_HCV_NI,
    p_HIV_HCV_ABS_NI  = p_HCV_NI,
    # Injection
    p_HIV_HCV_BUP_INJ  = p_HCV_TX_INJ,
    p_HIV_HCV_BUPC_INJ = p_HCV_TXC_INJ,
    p_HIV_HCV_MET_INJ  = p_HCV_TX_INJ,
    p_HIV_HCV_METC_INJ = p_HCV_TXC_INJ,
    p_HIV_HCV_REL_INJ  = p_HCV_REL_INJ,
    p_HIV_HCV_ODN_INJ  = p_HCV_REL_INJ,
    p_HIV_HCV_ABS_INJ  = p_HCV_NI,
    
    ### Costs ###
    # Treatment Costs
    c_BUP_TX  = rgamma(n_sim, shape = df_costs["shape", "BUP_TX"], scale = df_costs["scale", "BUP_TX"]), # BUP treatment costs - change to normal
    c_MET_TX  = rgamma(n_sim, shape = df_costs["shape", "MET_TX"], scale = df_costs["scale", "MET_TX"]), # MET treatment costs - change to normal
    
    # HRU Costs
    c_BUP_NI_HRU  = c_BUP_HRU, 
    c_BUPC_NI_HRU = c_BUPC_HRU, 
    c_MET_NI_HRU = c_MET_HRU,
    c_METC_NI_HRU = c_METC_HRU,
    c_ABS_NI_HRU = c_ABS_HRU, 
    c_REL_NI_HRU = c_REL_HRU, 
    c_ODN_NI_HRU  = c_ODN_HRU,
    c_ODF_NI_HRU  = c_ODF_HRU, 
    c_BUP_INJ_HRU  = c_BUP_HRU, 
    c_BUPC_INJ_HRU = c_BUPC_HRU, 
    c_MET_INJ_HRU = c_MET_HRU,
    c_METC_INJ_HRU = c_METC_HRU,
    c_ABS_INJ_HRU = c_ABS_HRU, 
    c_REL_INJ_HRU = c_REL_HRU, 
    c_ODN_INJ_HRU  = c_ODN_HRU,
    c_ODF_INJ_HRU  = c_ODF_HRU,
    
    # HIV Costs
    c_HIV_HRU = rgamma(n_sim, shape = df_costs["shape", "HIV_HRU"], scale = df_costs["scale", "HIV_HRU"]),
    c_HIV_ART = rgamma(n_sim, shape = df_costs["shape", "HIV_ART"], scale = df_costs["scale", "HIV_ART"]),
    
    # HCV Costs
    c_HCV_HRU = rgamma(n_sim, shape = df_costs["shape", "HCV_HRU"], scale = df_costs["scale", "HCV_HRU"]),
    c_HCV_DAA = rgamma(n_sim, shape = df_costs["shape", "HCV_DAA"], scale = df_costs["scale", "HCV_DAA"]),

    # Crime Costs
    c_BUP_NI_crime = c_BUP_crime,
    c_BUP_INJ_crime = c_BUP_crime,
    c_BUPC_NI_crime = c_BUPC_crime,
    c_BUPC_INJ_crime = c_BUPC_crime,
    c_MET_NI_crime = c_MET_crime,
    c_MET_INJ_crime = c_MET_crime,
    c_METC_NI_crime = c_METC_crime,
    c_METC_INJ_crime = c_METC_crime,
    c_REL_NI_crime = c_REL_crime,
    c_REL_INJ_crime = c_REL_crime,
    c_ODN_NI_crime = c_ODN_crime,
    c_ODN_INJ_crime = c_ODN_crime,
    
    ### Utilities ###
    # HIV/HCV negative
    u_BUP_NI_NEG  = u_BUP_NEG,
    u_BUPC_NI_NEG = u_BUPC_NEG,
    u_MET_NI_NEG  = u_MET_NEG,
    u_METC_NI_NEG = u_METC_NEG,
    u_REL_NI_NEG  = u_REL_NEG,
    u_ODN_NI_NEG  = u_ODN_NEG,
    u_ABS_NI_NEG  = u_ABS_NEG,
    
    u_BUP_INJ_NEG  = u_BUP_NEG,
    u_BUPC_INJ_NEG = u_BUPC_NEG,
    u_MET_INJ_NEG  = u_MET_NEG,
    u_METC_INJ_NEG = u_METC_NEG,
    u_REL_INJ_NEG  = u_REL_NEG,
    u_ODN_INJ_NEG  = u_ODN_NEG,
    u_ABS_INJ_NEG  = u_ABS_NEG,

    u_HIV_mult = u_HIV_mult, # HIV multiplier for negative states
    u_HCV_mult = u_HCV_mult, # HCV multiplier for negative states
    U_COI_mult = u_COI_mult)
  
  return(df_psa_params)
}

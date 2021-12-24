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
#' @param n_samp Sample size to determine dirichlet distribution variance.
#' @return 
#' A data frame with \code{n_sim} rows and {n_states} columns of parameters for PSA. 
#' Each row is a parameter set sampled from distributions that characterize 
#' their uncertainty
#' @examples 
#' generate_psa_params()
#' @export
generate_psa_params <- function(n_sim = n_sim, seed = seed, n_samp = n_samp,
                                file.death_hr = NULL,
                                file.frailty = NULL,
                                file.weibull_scale = NULL,
                                file.weibull_shape = NULL,
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
  df_weibull_scale <- read.csv(file = file.weibull_shape, row.names = 1, header = TRUE) # Weibull scale params
  df_weibull_shape <- read.csv(file = file.weibull_scale, row.names = 1, header = TRUE) # Weibull shape params
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
  
  # Number of simulations
  n_sim <- n_sim
  if(n_sim != nrow(df_calib_post)){
    warning("Number of PSA simulations and posterior draws not equal")
  }
  
  # Set seed for random number generator
  seed <- seed
  if (!missing(seed)) 
    set.seed(seed)
  #set_seed <- seed
  
  # Set sample-size for trial-based parameter uncertainty
  n_samp <- n_samp
  #write.csv(n_samp, file = "checks/n_samp.csv", row.names = TRUE)
  # Function to generate lognormal parameter
  location <- function(m = m, s = s){
    log(m^2 / sqrt(s^2 + m^2))
  }
  shape <- function(m = m, s = s){
    shape <- sqrt(log(1 + (s^2 / m^2)))
  }
  
  # Set up dirichlet random sample
  df_dirichlet_UP = df_UP * n_samp
  write.csv(df_dirichlet_UP, file = "checks/df_dirichlet_UP.csv", row.names = TRUE)
  
  # Non-injection
  # From BUP
  v_dirichlet_UP_BUP_NI = df_dirichlet_UP["BUP_NI",]
  m_BUP_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_BUP_NI["BUP_NI", "BUPC_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "MET_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "METC_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "ABS_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "REL_NI"])))

  # From BUPC
  v_dirichlet_UP_BUPC_NI = df_dirichlet_UP["BUPC_NI",]
  m_BUPC_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_BUPC_NI["BUPC_NI", "BUP_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "MET_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "METC_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "ABS_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "REL_NI"])) 
  
  # From MET
  v_dirichlet_UP_MET_NI = df_dirichlet_UP["MET_NI",]
  m_MET_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_MET_NI["MET_NI", "METC_NI"], v_dirichlet_UP_MET_NI["MET_NI", "BUP_NI"], v_dirichlet_UP_MET_NI["MET_NI", "BUPC_NI"], v_dirichlet_UP_MET_NI["MET_NI", "ABS_NI"], v_dirichlet_UP_MET_NI["MET_NI", "REL_NI"]))
  
  # From METC
  v_dirichlet_UP_METC_NI = df_dirichlet_UP["METC_NI",]
  m_METC_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_METC_NI["METC_NI", "MET_NI"], v_dirichlet_UP_METC_NI["METC_NI", "BUP_NI"], v_dirichlet_UP_METC_NI["METC_NI", "BUPC_NI"], v_dirichlet_UP_METC_NI["METC_NI", "ABS_NI"], v_dirichlet_UP_METC_NI["METC_NI", "REL_NI"]))
  
  # From ABS
  v_dirichlet_UP_ABS_NI = df_dirichlet_UP["ABS_NI",]
  m_ABS_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_ABS_NI["ABS_NI", "MET_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "METC_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "BUP_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "BUPC_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "REL_NI"]))
  
  # From REL
  v_dirichlet_UP_REL_NI = df_dirichlet_UP["REL_NI",]
  m_REL_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_REL_NI["REL_NI", "MET_NI"], v_dirichlet_UP_REL_NI["REL_NI", "METC_NI"], v_dirichlet_UP_REL_NI["REL_NI", "BUP_NI"], v_dirichlet_UP_REL_NI["REL_NI", "BUPC_NI"], v_dirichlet_UP_REL_NI["REL_NI", "ABS_NI"]))
  
  # From OD
  v_dirichlet_UP_OD_NI = df_dirichlet_UP["ODN_NI",]
  m_OD_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_OD_NI["ODN_NI", "MET_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "METC_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "BUP_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "BUPC_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "ABS_NI"], v_dirichlet_UP_OD_NI["ODN_NI", "REL_NI"]))
  
  # Injection
  # From BUP
  v_dirichlet_UP_BUP_INJ = df_dirichlet_UP["BUP_INJ",]
  m_BUP_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_BUP_INJ["BUP_INJ", "BUPC_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "MET_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "METC_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "ABS_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "REL_INJ"]))
  
  # From BUPC
  v_dirichlet_UP_BUPC_INJ = df_dirichlet_UP["BUPC_INJ",]
  m_BUPC_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "BUP_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "MET_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "METC_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "ABS_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "REL_INJ"]))
  
  # From MET
  v_dirichlet_UP_MET_INJ = df_dirichlet_UP["MET_INJ",]
  m_MET_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_MET_INJ["MET_INJ", "METC_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "BUP_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "BUPC_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "ABS_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "REL_INJ"]))
  
  # From METC
  v_dirichlet_UP_METC_INJ = df_dirichlet_UP["METC_INJ",]
  m_METC_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_METC_INJ["METC_INJ", "MET_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "BUP_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "BUPC_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "ABS_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "REL_INJ"]))
  
  # From ABS
  v_dirichlet_UP_ABS_INJ = df_dirichlet_UP["ABS_INJ",]
  m_ABS_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_ABS_INJ["ABS_INJ", "MET_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "METC_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "BUP_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "BUPC_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "REL_INJ"]))
  
  # From REL
  v_dirichlet_UP_REL_INJ = df_dirichlet_UP["REL_INJ",]
  m_REL_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_REL_INJ["REL_INJ", "MET_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "METC_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "BUP_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "BUPC_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "ABS_INJ"]))
  
  # From OD
  v_dirichlet_UP_OD_INJ = df_dirichlet_UP["ODN_INJ",]
  m_OD_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_OD_INJ["ODN_INJ", "MET_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "METC_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "BUP_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "BUPC_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "ABS_INJ"], v_dirichlet_UP_OD_INJ["ODN_INJ", "REL_INJ"]))
  
  
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
  c_OD_HRU  = rgamma(n_sim, shape = df_costs["shape", "OD_NI_HRU"] , scale = df_costs["scale", "OD_NI_HRU"] )
  
  # Crime Costs
  c_BUP_crime = rgamma(n_sim, shape = df_crime_costs["shape", "BUP"], scale = df_crime_costs["scale", "BUP"])
  c_BUPC_crime = rgamma(n_sim, shape = df_crime_costs["shape", "BUPC"], scale = df_crime_costs["scale", "BUPC"])
  c_MET_crime = rgamma(n_sim, shape = df_crime_costs["shape", "MET"], scale = df_crime_costs["scale", "MET"])
  c_METC_crime = rgamma(n_sim, shape = df_crime_costs["shape", "METC"], scale = df_crime_costs["scale", "METC"])
  c_REL_crime = rgamma(n_sim, shape = df_crime_costs["shape", "REL"], scale = df_crime_costs["scale", "REL"])
  c_ODN_crime = rgamma(n_sim, shape = df_crime_costs["shape", "REL"], scale = df_crime_costs["scale", "REL"])
  
  df_psa_params <- data.frame(
    ### Calibrated parameters
    df_calib_post, # Matrix of calibration parameters drawn from posterior distribution
    
    # Hazard ratios for death probability
    # Log-normal distribution
    # Non-injection
    hr_BUP_NI  = rlnorm(n_sim, location(m = df_death_hr["pe", "BUP_NI"], s = df_death_hr["sd", "BUP_NI"]), shape(m = df_death_hr["pe", "BUP_NI"], s = df_death_hr["sd", "BUP_NI"])),
    hr_BUPC_NI = rlnorm(n_sim, location(m = df_death_hr["pe", "BUPC_NI"], s = df_death_hr["sd", "BUPC_NI"]), shape(m = df_death_hr["pe", "BUPC_NI"], s = df_death_hr["sd", "BUPC_NI"])),
    hr_MET_NI  = rlnorm(n_sim, location(m = df_death_hr["pe", "MET_NI"], s = df_death_hr["sd", "MET_NI"]), shape(m = df_death_hr["pe", "MET_NI"], s = df_death_hr["sd", "MET_NI"])),
    hr_METC_NI = rlnorm(n_sim, location(m = df_death_hr["pe", "METC_NI"], s = df_death_hr["sd", "METC_NI"]), shape(m = df_death_hr["pe", "METC_NI"], s = df_death_hr["sd", "METC_NI"])),
    hr_REL_NI  = rlnorm(n_sim, location(m = df_death_hr["pe", "REL_NI"], s = df_death_hr["sd", "REL_NI"]), shape(m = df_death_hr["pe", "REL_NI"], s = df_death_hr["sd", "REL_NI"])),
    hr_ODN_NI  = rlnorm(n_sim, location(m = df_death_hr["pe", "REL_NI"], s = df_death_hr["sd", "REL_NI"]), shape(m = df_death_hr["pe", "REL_NI"], s = df_death_hr["sd", "REL_NI"])),
    #hr_ODN_NI  = hr_REL_NI,
    #hr_ABS_NI  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "ABS_NI"] ), sd = log(df_death_hr["sd", "ABS_NI"]))),
    hr_HIV_NI  = rlnorm(n_sim, location(m = df_death_hr["pe", "HIV_NI"], s = df_death_hr["sd", "HIV_NI"]), shape(m = df_death_hr["pe", "HIV_NI"], s = df_death_hr["sd", "HIV_NI"])),
    hr_HCV_NI  = rlnorm(n_sim, location(m = df_death_hr["pe", "HCV_NI"], s = df_death_hr["sd", "HCV_NI"]), shape(m = df_death_hr["pe", "HCV_NI"], s = df_death_hr["sd", "HCV_NI"])),
    hr_COI_NI  = rlnorm(n_sim, location(m = df_death_hr["pe", "COI_NI"], s = df_death_hr["sd", "COI_NI"]), shape(m = df_death_hr["pe", "COI_NI"], s = df_death_hr["sd", "COI_NI"])),
    
    # Injection
    hr_BUP_INJ  = rlnorm(n_sim, location(m = df_death_hr["pe", "BUP_INJ"], s = df_death_hr["sd", "BUP_INJ"]), shape(m = df_death_hr["pe", "BUP_INJ"], s = df_death_hr["sd", "BUP_INJ"])),
    hr_BUPC_INJ = rlnorm(n_sim, location(m = df_death_hr["pe", "BUPC_INJ"], s = df_death_hr["sd", "BUPC_INJ"]), shape(m = df_death_hr["pe", "BUPC_INJ"], s = df_death_hr["sd", "BUPC_INJ"])),
    hr_MET_INJ  = rlnorm(n_sim, location(m = df_death_hr["pe", "MET_INJ"], s = df_death_hr["sd", "MET_INJ"]), shape(m = df_death_hr["pe", "MET_INJ"], s = df_death_hr["sd", "MET_INJ"])),
    hr_METC_INJ = rlnorm(n_sim, location(m = df_death_hr["pe", "METC_INJ"], s = df_death_hr["sd", "METC_INJ"]), shape(m = df_death_hr["pe", "METC_INJ"], s = df_death_hr["sd", "METC_INJ"])),
    hr_REL_INJ  = rlnorm(n_sim, location(m = df_death_hr["pe", "REL_INJ"], s = df_death_hr["sd", "REL_INJ"]), shape(m = df_death_hr["pe", "REL_INJ"], s = df_death_hr["sd", "REL_INJ"])),
    hr_ODN_INJ  = rlnorm(n_sim, location(m = df_death_hr["pe", "REL_INJ"], s = df_death_hr["sd", "REL_INJ"]), shape(m = df_death_hr["pe", "REL_INJ"], s = df_death_hr["sd", "REL_INJ"])),
    #hr_ODN_INJ  = hr_REL_INJ,
    #hr_ABS_INJ  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "ABS_INJ"] ), sd = log(df_death_hr["sd", "ABS_INJ"]))),
    hr_HIV_INJ  = rlnorm(n_sim, location(m = df_death_hr["pe", "HIV_INJ"], s = df_death_hr["sd", "HIV_INJ"]), shape(m = df_death_hr["pe", "HIV_INJ"], s = df_death_hr["sd", "HIV_INJ"])),
    hr_HCV_INJ  = rlnorm(n_sim, location(m = df_death_hr["pe", "HCV_INJ"], s = df_death_hr["sd", "HCV_INJ"]), shape(m = df_death_hr["pe", "HCV_INJ"], s = df_death_hr["sd", "HCV_INJ"])),
    hr_COI_INJ  = rlnorm(n_sim, location(m = df_death_hr["pe", "COI_INJ"], s = df_death_hr["sd", "COI_INJ"]), shape(m = df_death_hr["pe", "COI_INJ"], s = df_death_hr["sd", "COI_INJ"])),

    # Frailty terms for successive episodes
    # Log-normal distribution (mean, sd)
    # Non-injection
    # BUP
    p_frailty_BUP_NI_2 = rlnorm(n_sim, location(m = df_frailty["pe", "BUP_NI_2"], s = df_frailty["sd", "BUP_NI_2"]), shape(m = df_frailty["pe", "BUP_NI_2"], s = df_frailty["sd", "BUP_NI_2"])),
    p_frailty_BUP_NI_3 = rlnorm(n_sim, location(m = df_frailty["pe", "BUP_NI_3"], s = df_frailty["sd", "BUP_NI_3"]), shape(m = df_frailty["pe", "BUP_NI_3"], s = df_frailty["sd", "BUP_NI_3"])),
    # BUPC
    p_frailty_BUPC_NI_2 = rlnorm(n_sim, location(m = df_frailty["pe", "BUPC_NI_2"], s = df_frailty["sd", "BUPC_NI_2"]), shape(m = df_frailty["pe", "BUPC_NI_2"], s = df_frailty["sd", "BUPC_NI_2"])),
    p_frailty_BUPC_NI_3 = rlnorm(n_sim, location(m = df_frailty["pe", "BUPC_NI_3"], s = df_frailty["sd", "BUPC_NI_3"]), shape(m = df_frailty["pe", "BUPC_NI_3"], s = df_frailty["sd", "BUPC_NI_3"])),
    # MET
    p_frailty_MET_NI_2 = rlnorm(n_sim, location(m = df_frailty["pe", "MET_NI_2"], s = df_frailty["sd", "MET_NI_2"]), shape(m = df_frailty["pe", "MET_NI_2"], s = df_frailty["sd", "MET_NI_2"])),
    p_frailty_MET_NI_3 = rlnorm(n_sim, location(m = df_frailty["pe", "MET_NI_3"], s = df_frailty["sd", "MET_NI_3"]), shape(m = df_frailty["pe", "MET_NI_3"], s = df_frailty["sd", "MET_NI_3"])),
    # METC
    p_frailty_METC_NI_2 = rlnorm(n_sim, location(m = df_frailty["pe", "METC_NI_2"], s = df_frailty["sd", "METC_NI_2"]), shape(m = df_frailty["pe", "METC_NI_2"], s = df_frailty["sd", "METC_NI_2"])),
    p_frailty_METC_NI_3 = rlnorm(n_sim, location(m = df_frailty["pe", "METC_NI_3"], s = df_frailty["sd", "METC_NI_3"]), shape(m = df_frailty["pe", "METC_NI_3"], s = df_frailty["sd", "METC_NI_3"])),
    # ABS
    p_frailty_ABS_NI_2 = rlnorm(n_sim, location(m = df_frailty["pe", "ABS_NI_2"], s = df_frailty["sd", "ABS_NI_2"]), shape(m = df_frailty["pe", "ABS_NI_2"], s = df_frailty["sd", "ABS_NI_2"])),
    p_frailty_ABS_NI_3 = rlnorm(n_sim, location(m = df_frailty["pe", "ABS_NI_3"], s = df_frailty["sd", "ABS_NI_3"]), shape(m = df_frailty["pe", "ABS_NI_3"], s = df_frailty["sd", "ABS_NI_3"])),
    # REL
    p_frailty_REL_NI_2 = rlnorm(n_sim, location(m = df_frailty["pe", "REL_NI_2"], s = df_frailty["sd", "REL_NI_2"]), shape(m = df_frailty["pe", "REL_NI_2"], s = df_frailty["sd", "REL_NI_2"])),
    p_frailty_REL_NI_3 = rlnorm(n_sim, location(m = df_frailty["pe", "REL_NI_3"], s = df_frailty["sd", "REL_NI_3"]), shape(m = df_frailty["pe", "REL_NI_3"], s = df_frailty["sd", "REL_NI_3"])),
    # OD
    p_frailty_OD_NI_2 = rlnorm(n_sim, location(m = df_frailty["pe", "REL_NI_2"], s = df_frailty["sd", "REL_NI_2"]), shape(m = df_frailty["pe", "REL_NI_2"], s = df_frailty["sd", "REL_NI_2"])),
    p_frailty_OD_NI_3 = rlnorm(n_sim, location(m = df_frailty["pe", "REL_NI_3"], s = df_frailty["sd", "REL_NI_3"]), shape(m = df_frailty["pe", "REL_NI_3"], s = df_frailty["sd", "REL_NI_3"])),
    
    # Injection
    # BUP
    p_frailty_BUP_INJ_2 = rlnorm(n_sim, location(m = df_frailty["pe", "BUP_INJ_2"], s = df_frailty["sd", "BUP_INJ_2"]), shape(m = df_frailty["pe", "BUP_INJ_2"], s = df_frailty["sd", "BUP_INJ_2"])),
    p_frailty_BUP_INJ_3 = rlnorm(n_sim, location(m = df_frailty["pe", "BUP_INJ_3"], s = df_frailty["sd", "BUP_INJ_3"]), shape(m = df_frailty["pe", "BUP_INJ_3"], s = df_frailty["sd", "BUP_INJ_3"])),
    # BUPC
    p_frailty_BUPC_INJ_2 = rlnorm(n_sim, location(m = df_frailty["pe", "BUPC_INJ_2"], s = df_frailty["sd", "BUPC_INJ_2"]), shape(m = df_frailty["pe", "BUPC_INJ_2"], s = df_frailty["sd", "BUPC_INJ_2"])),
    p_frailty_BUPC_INJ_3 = rlnorm(n_sim, location(m = df_frailty["pe", "BUPC_INJ_3"], s = df_frailty["sd", "BUPC_INJ_3"]), shape(m = df_frailty["pe", "BUPC_INJ_3"], s = df_frailty["sd", "BUPC_INJ_3"])),
    # MET
    p_frailty_MET_INJ_2 = rlnorm(n_sim, location(m = df_frailty["pe", "MET_INJ_2"], s = df_frailty["sd", "MET_INJ_2"]), shape(m = df_frailty["pe", "MET_INJ_2"], s = df_frailty["sd", "MET_INJ_2"])),
    p_frailty_MET_INJ_3 = rlnorm(n_sim, location(m = df_frailty["pe", "MET_INJ_3"], s = df_frailty["sd", "MET_INJ_3"]), shape(m = df_frailty["pe", "MET_INJ_3"], s = df_frailty["sd", "MET_INJ_3"])),
    # METC
    p_frailty_METC_INJ_2 = rlnorm(n_sim, location(m = df_frailty["pe", "METC_INJ_2"], s = df_frailty["sd", "METC_INJ_2"]), shape(m = df_frailty["pe", "METC_INJ_2"], s = df_frailty["sd", "METC_INJ_2"])),
    p_frailty_METC_INJ_3 = rlnorm(n_sim, location(m = df_frailty["pe", "METC_INJ_3"], s = df_frailty["sd", "METC_INJ_3"]), shape(m = df_frailty["pe", "METC_INJ_3"], s = df_frailty["sd", "METC_INJ_3"])),
    # ABS
    p_frailty_ABS_INJ_2 = rlnorm(n_sim, location(m = df_frailty["pe", "ABS_INJ_2"], s = df_frailty["sd", "ABS_INJ_2"]), shape(m = df_frailty["pe", "ABS_INJ_2"], s = df_frailty["sd", "ABS_INJ_2"])),
    p_frailty_ABS_INJ_3 = rlnorm(n_sim, location(m = df_frailty["pe", "ABS_INJ_3"], s = df_frailty["sd", "ABS_INJ_3"]), shape(m = df_frailty["pe", "ABS_INJ_3"], s = df_frailty["sd", "ABS_INJ_3"])),
    # REL
    p_frailty_REL_INJ_2 = rlnorm(n_sim, location(m = df_frailty["pe", "REL_INJ_2"], s = df_frailty["sd", "REL_INJ_2"]), shape(m = df_frailty["pe", "REL_INJ_2"], s = df_frailty["sd", "REL_INJ_2"])),
    p_frailty_REL_INJ_3 = rlnorm(n_sim, location(m = df_frailty["pe", "REL_INJ_3"], s = df_frailty["sd", "REL_INJ_3"]), shape(m = df_frailty["pe", "REL_INJ_3"], s = df_frailty["sd", "REL_INJ_3"])),
    # OD
    p_frailty_OD_INJ_2 = rlnorm(n_sim, location(m = df_frailty["pe", "REL_INJ_2"], s = df_frailty["sd", "REL_INJ_2"]), shape(m = df_frailty["pe", "REL_INJ_2"], s = df_frailty["sd", "REL_INJ_2"])),
    p_frailty_OD_INJ_3 = rlnorm(n_sim, location(m = df_frailty["pe", "REL_INJ_3"], s = df_frailty["sd", "REL_INJ_3"]), shape(m = df_frailty["pe", "REL_INJ_3"], s = df_frailty["sd", "REL_INJ_3"])),
    
    # Weibull scale
    # Non-injection
    p_weibull_scale_BUP_NI  = rnorm(n_sim, mean = df_weibull_scale["pe", "BUP_NI"], sd = df_weibull_scale["sd", "BUP_NI"]),
    p_weibull_scale_BUPC_NI = rnorm(n_sim, mean = df_weibull_scale["pe", "BUPC_NI"], sd = df_weibull_scale["sd", "BUPC_NI"]),
    p_weibull_scale_MET_NI = rnorm(n_sim, mean = df_weibull_scale["pe", "MET_NI"], sd = df_weibull_scale["sd", "MET_NI"]),
    p_weibull_scale_METC_NI = rnorm(n_sim, mean = df_weibull_scale["pe", "METC_NI"], sd = df_weibull_scale["sd", "METC_NI"]),
    p_weibull_scale_ABS_NI = rnorm(n_sim, mean = df_weibull_scale["pe", "ABS_NI"], sd = df_weibull_scale["sd", "ABS_NI"]),
    p_weibull_scale_REL_NI = rnorm(n_sim, mean = df_weibull_scale["pe", "REL_NI"], sd = df_weibull_scale["sd", "REL_NI"]),
        
    # Injection
    p_weibull_scale_BUP_INJ = rnorm(n_sim, mean = df_weibull_scale["pe", "BUP_INJ"], sd = df_weibull_scale["sd", "BUP_INJ"]),
    p_weibull_scale_BUPC_INJ = rnorm(n_sim, mean = df_weibull_scale["pe", "BUPC_INJ"], sd = df_weibull_scale["sd", "BUPC_INJ"]),
    p_weibull_scale_MET_INJ = rnorm(n_sim, mean = df_weibull_scale["pe", "MET_INJ"], sd = df_weibull_scale["sd", "MET_INJ"]),
    p_weibull_scale_METC_INJ = rnorm(n_sim, mean = df_weibull_scale["pe", "METC_INJ"], sd = df_weibull_scale["sd", "METC_INJ"]),
    p_weibull_scale_ABS_INJ = rnorm(n_sim, mean = df_weibull_scale["pe", "ABS_INJ"], sd = df_weibull_scale["sd", "ABS_INJ"]),
    p_weibull_scale_REL_INJ = rnorm(n_sim, mean = df_weibull_scale["pe", "REL_INJ"], sd = df_weibull_scale["sd", "REL_INJ"]),
    
    # Weibull shape
    # Non-injection
    p_weibull_shape_BUP_NI = rnorm(n_sim, mean = df_weibull_shape["pe", "BUP_NI"], sd = df_weibull_shape["sd", "BUP_NI"]),
    p_weibull_shape_BUPC_NI = rnorm(n_sim, mean = df_weibull_shape["pe", "BUPC_NI"], sd = df_weibull_shape["sd", "BUPC_NI"]),
    p_weibull_shape_MET_NI = rnorm(n_sim, mean = df_weibull_shape["pe", "MET_NI"], sd = df_weibull_shape["sd", "MET_NI"]),
    p_weibull_shape_METC_NI = rnorm(n_sim, mean = df_weibull_shape["pe", "METC_NI"], sd = df_weibull_shape["sd", "METC_NI"]),
    p_weibull_shape_ABS_NI = rnorm(n_sim, mean = df_weibull_shape["pe", "ABS_NI"], sd = df_weibull_shape["sd", "ABS_NI"]),
    p_weibull_shape_REL_NI = rnorm(n_sim, mean = df_weibull_shape["pe", "REL_NI"], sd = df_weibull_shape["sd", "REL_NI"]),
    
    # Injection
    p_weibull_shape_BUP_INJ = rnorm(n_sim, mean = df_weibull_shape["pe", "BUP_INJ"], sd = df_weibull_shape["sd", "BUP_INJ"]),
    p_weibull_shape_BUPC_INJ = rnorm(n_sim, mean = df_weibull_shape["pe", "BUPC_INJ"], sd = df_weibull_shape["sd", "BUPC_INJ"]),
    p_weibull_shape_MET_INJ = rnorm(n_sim, mean = df_weibull_shape["pe", "MET_INJ"], sd = df_weibull_shape["sd", "MET_INJ"]),
    p_weibull_shape_METC_INJ = rnorm(n_sim, mean = df_weibull_shape["pe", "METC_INJ"], sd = df_weibull_shape["sd", "METC_INJ"]),
    p_weibull_shape_ABS_INJ = rnorm(n_sim, mean = df_weibull_shape["pe", "ABS_INJ"], sd = df_weibull_shape["sd", "ABS_INJ"]),
    p_weibull_shape_REL_INJ = rnorm(n_sim, mean = df_weibull_shape["pe", "REL_INJ"], sd = df_weibull_shape["sd", "REL_INJ"]),
    
    ### Transition probabilities conditional on leaving (use Dirichlet)
    #write.csv(df_UP, file = "checks/df_UP.csv", row.names = TRUE),
    #df_dirichlet_UP = mutate_if(df_UP, is.numeric, ~ . * n_samp), # weight unconditional matrix by sample size to generate dirichlet PSA
    #check = df_UP * n_samp,
    #df_dirichlet_UP = df_UP * n_samp,
    
    #write.csv(df_dirichlet_UP, file = "checks/df_dirichlet_UP.csv", row.names = TRUE),
    
    # Non-Injection
    # From BUP
    #v_dirichlet_UP_BUP_NI = df_dirichlet_UP["BUP_NI",],
    #m_BUP_UP_NI = matrix(0, nrow = n_sim, ncol = ncol(v_dirichlet_UP_BUP_NI)),
    #m_BUP_UP_NI = as.matrix(rdirichlet(n_sim, c(v_dirichlet_UP_BUP_NI["BUP_NI", "BUPC_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "MET_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "METC_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "ABS_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "REL_NI"]))),
    p_BUP_BUPC_NI = m_BUP_UP_NI[,1], # assign probabilities by place
    p_BUP_MET_NI  = m_BUP_UP_NI[,2],
    p_BUP_METC_NI = m_BUP_UP_NI[,3],
    p_BUP_ABS_NI  = m_BUP_UP_NI[,4],
    p_BUP_REL_NI  = m_BUP_UP_NI[,5],
    
    # From BUPC
    #v_dirichlet_UP_BUPC_NI = df_dirichlet_UP["BUPC_NI",],
    #m_BUPC_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_BUPC_NI["BUPC_NI", "BUP_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "MET_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "METC_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "ABS_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "REL_NI"])),
    p_BUPC_BUP_NI  = m_BUPC_UP_NI[,1],
    p_BUPC_MET_NI  = m_BUPC_UP_NI[,2],
    p_BUPC_METC_NI = m_BUPC_UP_NI[,3],
    p_BUPC_ABS_NI  = m_BUPC_UP_NI[,4],
    p_BUPC_REL_NI  = m_BUPC_UP_NI[,5],
    
    # From MET
    #v_dirichlet_UP_MET_NI = df_dirichlet_UP["MET_NI",],
    #m_MET_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_MET_NI["MET_NI", "METC_NI"], v_dirichlet_UP_MET_NI["MET_NI", "BUP_NI"], v_dirichlet_UP_MET_NI["MET_NI", "BUPC_NI"], v_dirichlet_UP_MET_NI["MET_NI", "ABS_NI"], v_dirichlet_UP_MET_NI["MET_NI", "REL_NI"])),
    p_MET_METC_NI = m_MET_UP_NI[,1],
    p_MET_BUP_NI  = m_MET_UP_NI[,2],
    p_MET_BUPC_NI = m_MET_UP_NI[,3],
    p_MET_ABS_NI  = m_MET_UP_NI[,4],
    p_MET_REL_NI  = m_MET_UP_NI[,5],
    
    # From METC
    #v_dirichlet_UP_METC_NI = df_dirichlet_UP["METC_NI",],
    #m_METC_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_METC_NI["METC_NI", "MET_NI"], v_dirichlet_UP_METC_NI["METC_NI", "BUP_NI"], v_dirichlet_UP_METC_NI["METC_NI", "BUPC_NI"], v_dirichlet_UP_METC_NI["METC_NI", "ABS_NI"], v_dirichlet_UP_METC_NI["METC_NI", "REL_NI"])),
    p_METC_MET_NI  = m_METC_UP_NI[,1],
    p_METC_BUP_NI  = m_METC_UP_NI[,2],
    p_METC_BUPC_NI = m_METC_UP_NI[,3],
    p_METC_ABS_NI  = m_METC_UP_NI[,4],
    p_METC_REL_NI  = m_METC_UP_NI[,5],
    
    # From ABS
    #v_dirichlet_UP_ABS_NI = df_dirichlet_UP["ABS_NI",],
    #m_ABS_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_ABS_NI["ABS_NI", "MET_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "METC_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "BUP_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "BUPC_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "REL_NI"])),
    p_ABS_MET_NI  = m_ABS_UP_NI[,1],
    p_ABS_METC_NI = m_ABS_UP_NI[,2],
    p_ABS_BUP_NI  = m_ABS_UP_NI[,3],
    p_ABS_BUPC_NI = m_ABS_UP_NI[,4],
    p_ABS_REL_NI  = m_ABS_UP_NI[,5],
    
    # From REL
    #v_dirichlet_UP_REL_NI = df_dirichlet_UP["REL_NI",],
    #m_REL_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_REL_NI["REL_NI", "MET_NI"], v_dirichlet_UP_REL_NI["REL_NI", "METC_NI"], v_dirichlet_UP_REL_NI["REL_NI", "BUP_NI"], v_dirichlet_UP_REL_NI["REL_NI", "BUPC_NI"], v_dirichlet_UP_REL_NI["REL_NI", "ABS_NI"])),
    p_REL_MET_NI  = m_REL_UP_NI[,1],
    p_REL_METC_NI = m_REL_UP_NI[,2],
    p_REL_BUP_NI  = m_REL_UP_NI[,3],
    p_REL_BUPC_NI = m_REL_UP_NI[,4],
    p_REL_ABS_NI  = m_REL_UP_NI[,5],
    
    # From OD
    #v_dirichlet_UP_OD_NI = df_dirichlet_UP["ODN_NI",],
    #m_OD_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_OD_NI["OD_NI", "MET_NI"], v_dirichlet_UP_OD_NI["OD_NI", "METC_NI"], v_dirichlet_UP_OD_NI["OD_NI", "BUP_NI"], v_dirichlet_UP_OD_NI["OD_NI", "BUPC_NI"], v_dirichlet_UP_OD_NI["OD_NI", "ABS_NI"], v_dirichlet_UP_OD_NI["OD_NI", "REL_NI"])),
    p_OD_MET_NI  = m_OD_UP_NI[,1],
    p_OD_METC_NI = m_OD_UP_NI[,2],
    p_OD_BUP_NI  = m_OD_UP_NI[,3],
    p_OD_BUPC_NI = m_OD_UP_NI[,4],
    p_OD_ABS_NI  = m_OD_UP_NI[,5],
    p_OD_REL_NI  = m_OD_UP_NI[,6],
    
    # Injection
    # From BUP
    #v_dirichlet_UP_BUP_INJ = df_dirichlet_UP["BUP_INJ",],
    #m_BUP_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_BUP_INJ["BUP_INJ", "BUPC_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "MET_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "METC_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "ABS_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "REL_INJ"])),
    p_BUP_BUPC_INJ = m_BUP_UP_INJ[,1], # assign probabilities by place
    p_BUP_MET_INJ  = m_BUP_UP_INJ[,2],
    p_BUP_METC_INJ = m_BUP_UP_INJ[,3],
    p_BUP_ABS_INJ  = m_BUP_UP_INJ[,4],
    p_BUP_REL_INJ  = m_BUP_UP_INJ[,5],
    
    # From BUPC
    #v_dirichlet_UP_BUPC_INJ = df_dirichlet_UP["BUPC_INJ",],
    #m_BUPC_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "BUP_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "MET_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "METC_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "ABS_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "REL_INJ"])),
    p_BUPC_BUP_INJ  = m_BUPC_UP_INJ[,1],
    p_BUPC_MET_INJ  = m_BUPC_UP_INJ[,2],
    p_BUPC_METC_INJ = m_BUPC_UP_INJ[,3],
    p_BUPC_ABS_INJ  = m_BUPC_UP_INJ[,4],
    p_BUPC_REL_INJ  = m_BUPC_UP_INJ[,5],
    
    # From MET
    #v_dirichlet_UP_MET_INJ = df_dirichlet_UP["MET_INJ",],
    #m_MET_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_MET_INJ["MET_INJ", "METC_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "BUP_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "BUPC_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "ABS_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "REL_INJ"])),
    p_MET_METC_INJ = m_MET_UP_INJ[,1],
    p_MET_BUP_INJ  = m_MET_UP_INJ[,2],
    p_MET_BUPC_INJ = m_MET_UP_INJ[,3],
    p_MET_ABS_INJ  = m_MET_UP_INJ[,4],
    p_MET_REL_INJ  = m_MET_UP_INJ[,5],
    
    # From METC
    #v_dirichlet_UP_METC_INJ = df_dirichlet_UP["METC_INJ",],
    #m_METC_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_METC_INJ["METC_INJ", "MET_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "BUP_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "BUPC_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "ABS_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "REL_INJ"])),
    p_METC_MET_INJ  = m_METC_UP_INJ[,1],
    p_METC_BUP_INJ  = m_METC_UP_INJ[,2],
    p_METC_BUPC_INJ = m_METC_UP_INJ[,3],
    p_METC_ABS_INJ  = m_METC_UP_INJ[,4],
    p_METC_REL_INJ  = m_METC_UP_INJ[,5],
    
    # From ABS
    #v_dirichlet_UP_ABS_INJ = df_dirichlet_UP["ABS_INJ",],
    #m_ABS_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_ABS_INJ["ABS_INJ", "MET_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "METC_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "BUP_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "BUPC_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "REL_INJ"])),
    p_ABS_MET_INJ  = m_ABS_UP_INJ[,1],
    p_ABS_METC_INJ = m_ABS_UP_INJ[,2],
    p_ABS_BUP_INJ  = m_ABS_UP_INJ[,3],
    p_ABS_BUPC_INJ = m_ABS_UP_INJ[,4],
    p_ABS_REL_INJ  = m_ABS_UP_INJ[,5],
    
    # From REL
    #v_dirichlet_UP_REL_INJ = df_dirichlet_UP["REL_INJ",],
    #m_REL_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_REL_INJ["REL_INJ", "MET_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "METC_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "BUP_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "BUPC_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "ABS_INJ"])),
    p_REL_MET_INJ  = m_REL_UP_INJ[,1],
    p_REL_METC_INJ = m_REL_UP_INJ[,2],
    p_REL_BUP_INJ  = m_REL_UP_INJ[,3],
    p_REL_BUPC_INJ = m_REL_UP_INJ[,4],
    p_REL_ABS_INJ  = m_REL_UP_INJ[,5],
    
    # From OD
    #v_dirichlet_UP_OD_INJ = df_dirichlet_UP["ODN_INJ",],
    #m_OD_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_OD_INJ["OD_INJ", "MET_INJ"], v_dirichlet_UP_OD_INJ["OD_INJ", "METC_INJ"], v_dirichlet_UP_OD_INJ["OD_INJ", "BUP_INJ"], v_dirichlet_UP_OD_INJ["OD_INJ", "BUPC_INJ"], v_dirichlet_UP_OD_INJ["OD_INJ", "ABS_INJ"], v_dirichlet_UP_OD_INJ["OD_INJ", "REL_INJ"])),
    p_OD_MET_INJ  = m_OD_UP_INJ[,1],
    p_OD_METC_INJ = m_OD_UP_INJ[,2],
    p_OD_BUP_INJ  = m_OD_UP_INJ[,3],
    p_OD_BUPC_INJ = m_OD_UP_INJ[,4],
    p_OD_ABS_INJ  = m_OD_UP_INJ[,5],
    p_OD_REL_INJ  = m_OD_UP_INJ[,6],
    
    ### Overdose ###
    #n_fent_OD = rgamma(n_sim, shape = df_overdose["shape", "fent_OD_rate"], scale = df_overdose["scale", "fent_OD_rate"]),
    #p_fent_exp = rbeta(n_sim, shape1 = df_overdose["shape1", "fent_exp_prob"], shape2 = df_overdose["shape2", "fent_exp_prob"]),
    p_ni_fent_reduction = rbeta(n_sim, shape1 = df_overdose["shape1", "ni_fent_reduction"], shape2 = df_overdose["shape2", "ni_fent_reduction"]),
    p_witness = rbeta(n_sim, shape1 = df_overdose["shape1", "witness_prob"], shape2 = df_overdose["shape2", "witness_prob"]),
    p_attended = rbeta(n_sim, shape1 = df_overdose["shape1", "attended_prob"], shape2 = df_overdose["shape2", "attended_prob"]),
    p_NX_used = rbeta(n_sim, shape1 = df_overdose["shape1", "NX_prob"], shape2 = df_overdose["shape2", "NX_prob"]),
    p_NX_success = rbeta(n_sim, shape1 = df_overdose["shape1", "NX_success_prob"], shape2 = df_overdose["shape2", "NX_success_prob"]),
    
    ### Fentanyl ###
    # PLACEHOLDER
    
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
    c_OD_NI_HRU  = c_OD_HRU, 
    c_BUP_INJ_HRU  = c_BUP_HRU, 
    c_BUPC_INJ_HRU = c_BUPC_HRU, 
    c_MET_INJ_HRU = c_MET_HRU,
    c_METC_INJ_HRU = c_METC_HRU,
    c_ABS_INJ_HRU = c_ABS_HRU, 
    c_REL_INJ_HRU = c_REL_HRU, 
    c_OD_INJ_HRU  = c_OD_HRU, 
    
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
    # Consider having equivalent QALYs between injection and non-injection
    # HIV/HCV negative
    u_BUP_NI_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "BUP_NI_NEG"], sd = df_qalys["sd", "BUP_NI_NEG"], b = 1),
    u_BUPC_NI_NEG = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "BUPC_NI_NEG"], sd = df_qalys["sd", "BUPC_NI_NEG"], b = 1),
    u_MET_NI_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "MET_NI_NEG"], sd = df_qalys["sd", "MET_NI_NEG"], b = 1),
    u_METC_NI_NEG = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "METC_NI_NEG"], sd = df_qalys["sd", "METC_NI_NEG"], b = 1),
    u_REL_NI_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "REL_NI_NEG"], sd = df_qalys["sd", "REL_NI_NEG"], b = 1),
    u_ODN_NI_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "ODN_NI_NEG"], sd = df_qalys["sd", "ODN_NI_NEG"], b = 1),
    u_ABS_NI_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "ABS_NI_NEG"], sd = df_qalys["sd", "ABS_NI_NEG"], b = 1),
    
    u_BUP_INJ_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "BUP_INJ_NEG"], sd = df_qalys["sd", "BUP_INJ_NEG"], b = 1),
    u_BUPC_INJ_NEG = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "BUPC_INJ_NEG"], sd = df_qalys["sd", "BUPC_INJ_NEG"], b = 1),
    u_MET_INJ_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "MET_INJ_NEG"], sd = df_qalys["sd", "MET_INJ_NEG"], b = 1),
    u_METC_INJ_NEG = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "METC_INJ_NEG"], sd = df_qalys["sd", "METC_INJ_NEG"], b = 1),
    u_REL_INJ_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "REL_INJ_NEG"], sd = df_qalys["sd", "REL_INJ_NEG"], b = 1),
    u_ODN_INJ_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "ODN_INJ_NEG"], sd = df_qalys["sd", "ODN_INJ_NEG"], b = 1),
    u_ABS_INJ_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "ABS_INJ_NEG"], sd = df_qalys["sd", "ABS_INJ_NEG"], b = 1),

    # Consider beta (other?) distributions for these
    u_HIV_mult = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "HIV_mult"], sd = df_qalys["sd", "HIV_mult"], b = 1), # HIV multiplier for negative states
    u_HCV_mult = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "HCV_mult"], sd = df_qalys["sd", "HCV_mult"], b = 1), # HCV multiplier for negative states
    U_COI_mult = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "COI_mult"], sd = df_qalys["sd", "COI_mult"], b = 1))
  
  return(df_psa_params)
}

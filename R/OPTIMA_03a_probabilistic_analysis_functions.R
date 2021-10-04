#' Generate PSA dataset of CEA parameters
#'
#' \code{generate_psa_params} generates PSA input dataset by sampling decision 
#' model parameters from their distributions. The sample of the calibrated
#' parameters is a draw from their posterior distribution obtained with the
#' IMIS algorithm.
#' @param n_sim Number of PSA samples.
#' @param seed Seed for reproducibility of Monte Carlo sampling.
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
                                file.hiv = NULL,
                                file.hcv = NULL,
                                file.costs = NULL,
                                file.crime_costs = NULL,
                                file.qalys = NULL){
  
  #Load files with parameter distribution values
  df_death_hr <- read.csv(file = file.death_hr, row.names = 1, header = TRUE) # Mortality hazard ratios
  df_frailty <- read.csv(file = file.frailty, row.names = 1, header = TRUE) # Episode frailty params
  df_weibull_scale <- read.csv(file = file.weibull_shape, row.names = 1, header = TRUE) # Weibull scale params
  df_weibull_shape <- read.csv(file = file.weibull_scale, row.names = 1, header = TRUE) # Weibull shape params
  df_UP <- read.csv(file = file.unconditional, row.names = 1, header = TRUE) # Unconditional transition probs
  df_overdose <- read.csv(file = file.overdose, row.names = 1, header = TRUE) # Overdose params
  df_hiv <- read.csv(file = file.hiv, row.names = 1, header = TRUE) # HIV seroconversion probs
  df_hcv <- read.csv(file = file.hcv, row.names = 1, header = TRUE) # HCV seroconversion probs
  df_costs <- read.csv(file = file.costs, row.names = 1, header = TRUE) # All costs excluding crime
  df_crime_costs <- read.csv(file = file.crime_costs, row.names = 1, header = TRUE) # Age-dependent crime costs
  df_qalys <- read.csv(file = file.qalys, row.names = 1, header = TRUE) # QALYs
  
  ## Load calibrated parameters
  n_sim <- nrow(m_calib_post)
  set_seed <- seed
  df_psa_params <- data.frame(
    ### Calibrated parameters
    m_calib_post, # Matrix of calibration parameters drawn from posterior distribution
    
    # Hazard ratios for death probability
    # Log-normal distribution
    # Non-injection
    hr_BUP_NI  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "BUP_NI"] ), sd = ((log(df_death_hr["high", "BUP_NI"])  - log(df_death_hr["low", "BUP_NI"])) / 2 / 1.96))),
    hr_BUPC_NI = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "BUPC_NI"]), sd = ((log(df_death_hr["high", "BUPC_NI"]) - log(df_death_hr["low", "BUPC_NI"])) / 2 / 1.96))),
    hr_MET_NI  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "MET_NI"] ), sd = ((log(df_death_hr["high", "MET_NI"])  - log(df_death_hr["low", "MET_NI"])) / 2 / 1.96))),
    hr_METC_NI = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "METC_NI"]), sd = ((log(df_death_hr["high", "METC_NI"]) - log(df_death_hr["low", "METC_NI"])) / 2 / 1.96))),
    hr_REL_NI  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "REL_NI"] ), sd = ((log(df_death_hr["high", "REL_NI"])  - log(df_death_hr["low", "REL_NI"])) / 2 / 1.96))),
    hr_ODN_NI  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "ODN_NI"] ), sd = ((log(df_death_hr["high", "ODN_NI"])  - log(df_death_hr["low", "ODN_NI"])) / 2 / 1.96))),
    hr_ABS_NI  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "ABS_NI"] ), sd = ((log(df_death_hr["high", "ABS_NI"])  - log(df_death_hr["low", "ABS_NI"])) / 2 / 1.96))),
    hr_HIV_NI  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "HIV_NI"] ), sd = ((log(df_death_hr["high", "HIV_NI"])  - log(df_death_hr["low", "HIV_NI"])) / 2 / 1.96))),
    hr_HCV_NI  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "HCV_NI"] ), sd = ((log(df_death_hr["high", "HCV_NI"])  - log(df_death_hr["low", "HCV_NI"])) / 2 / 1.96))),
    hr_COI_NI  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "COI_NI"] ), sd = ((log(df_death_hr["high", "COI_NI"])  - log(df_death_hr["low", "COI_NI"])) / 2 / 1.96))),
    
    # Injection
    hr_BUP_INJ  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "BUP_INJ"] ), sd = ((log(df_death_hr["high", "BUP_INJ"])  - log(df_death_hr["low", "BUP_INJ"])) / 2 / 1.96))),
    hr_BUPC_INJ = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "BUPC_INJ"]), sd = ((log(df_death_hr["high", "BUPC_INJ"]) - log(df_death_hr["low", "BUPC_INJ"])) / 2 / 1.96))),
    hr_MET_INJ  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "MET_INJ"] ), sd = ((log(df_death_hr["high", "MET_INJ"])  - log(df_death_hr["low", "MET_INJ"])) / 2 / 1.96))),
    hr_METC_INJ = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "METC_INJ"]), sd = ((log(df_death_hr["high", "METC_INJ"]) - log(df_death_hr["low", "METC_INJ"])) / 2 / 1.96))),
    hr_REL_INJ  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "REL_INJ"] ), sd = ((log(df_death_hr["high", "REL_INJ"])  - log(df_death_hr["low", "REL_INJ"])) / 2 / 1.96))),
    hr_ODN_INJ  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "ODN_INJ"] ), sd = ((log(df_death_hr["high", "ODN_INJ"])  - log(df_death_hr["low", "ODN_INJ"])) / 2 / 1.96))),
    hr_ABS_INJ  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "ABS_INJ"] ), sd = ((log(df_death_hr["high", "ABS_INJ"])  - log(df_death_hr["low", "ABS_INJ"])) / 2 / 1.96))),
    hr_HIV_INJ  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "HIV_INJ"] ), sd = ((log(df_death_hr["high", "HIV_INJ"])  - log(df_death_hr["low", "HIV_INJ"])) / 2 / 1.96))),
    hr_HCV_INJ  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "HCV_INJ"] ), sd = ((log(df_death_hr["high", "HCV_INJ"])  - log(df_death_hr["low", "HCV_INJ"])) / 2 / 1.96))),
    hr_COI_INJ  = exp(rnorm(n_sim, mean = log(df_death_hr["pe", "COI_INJ"] ), sd = ((log(df_death_hr["high", "COI_INJ"])  - log(df_death_hr["low", "COI_INJ"])) / 2 / 1.96))),

    # Frailty terms for successive episodes
    # Log-normal distribution (mean, sd)
    # Non-injection
    # BUP
    #p_frailty_BUP_NI_1 = rep(1, n_sim),
    p_frailty_BUP_NI_2 = rnorm(n_sim, mean = log(df_frailty["pe", "BUP_NI_2"]), sd = ((log(df_frailty["high", "BUP_NI_2"]) - log(df_frailty["low", "BUP_NI_2"])) / 2 / 1.96)),
    p_frailty_BUP_NI_3 = rnorm(n_sim, mean = log(df_frailty["pe", "BUP_NI_3"]), sd = ((log(df_frailty["high", "BUP_NI_3"]) - log(df_frailty["low", "BUP_NI_3"])) / 2 / 1.96)),
    # BUP + concurrent opioid
    #p_frailty_BUPC_NI_1 = rep(1, n_sim),
    p_frailty_BUPC_NI_2 = rnorm(n_sim, mean = df_frailty["pe", "BUPC_NI_2"], sd = df_frailty["sd", "BUPC_NI_2"]),
    p_frailty_BUPC_NI_3 = rnorm(n_sim, mean = df_frailty["pe", "BUPC_NI_3"], sd = df_frailty["sd", "BUPC_NI_3"]),
    # MET
    p_frailty_MET_NI_1 = rep(1, n_sim),
    p_frailty_MET_NI_2 = rnorm(n_sim, mean = df_frailty["pe", "MET_NI_2"], sd = df_frailty["sd", "MET_NI_2"]),
    p_frailty_MET_NI_3 = rnorm(n_sim, mean = df_frailty["pe", "MET_NI_3"], sd = df_frailty["sd", "MET_NI_3"]),
    # MET + concurrent opioid
    p_frailty_METC_NI_1 = rep(1, n_sim),
    p_frailty_METC_NI_2 = rnorm(n_sim, mean = df_frailty["pe", "METC_NI_2"], sd = df_frailty["sd", "METC_NI_2"]),
    p_frailty_METC_NI_3 = rnorm(n_sim, mean = df_frailty["pe", "METC_NI_3"], sd = df_frailty["sd", "METC_NI_3"]),
    # ABS
    p_frailty_ABS_NI_1 = rep(1, n_sim),
    p_frailty_ABS_NI_2 = rnorm(n_sim, mean = df_frailty["pe", "ABS_NI_2"], sd = df_frailty["sd", "ABS_NI_2"]),
    p_frailty_ABS_NI_3 = rnorm(n_sim, mean = df_frailty["pe", "ABS_NI_3"], sd = df_frailty["sd", "ABS_NI_3"]),
    # REL
    p_frailty_REL_NI_1 = rep(1, n_sim),
    p_frailty_REL_NI_2 = rnorm(n_sim, mean = df_frailty["pe", "REL_NI_2"], sd = df_frailty["sd", "REL_NI_2"]),
    p_frailty_REL_NI_3 = rnorm(n_sim, mean = df_frailty["pe", "REL_NI_3"], sd = df_frailty["sd", "REL_NI_3"]),
    # OD
    p_frailty_OD_NI_1  = rep(1, n_sim),
    p_frailty_OD_NI_2  = rnorm(n_sim, mean = df_frailty["pe", "OD_NI_2"], sd = df_frailty["sd", "OD_NI_2"]),
    p_frailty_OD_NI_3  = rnorm(n_sim, mean = df_frailty["pe", "OD_NI_3"], sd = df_frailty["sd", "OD_NI_3"]),
    
    # Injection
    # BUP
    p_frailty_BUP_INJ_1 = rep(1, n_sim),
    p_frailty_BUP_INJ_2 = rnorm(n_sim, mean = df_frailty["pe", "BUP_INJ_2"], sd = df_frailty["sd", "BUP_INJ_2"]),
    p_frailty_BUP_INJ_3 = rnorm(n_sim, mean = df_frailty["pe", "BUP_INJ_3"], sd = df_frailty["sd", "BUP_INJ_3"]),
    # BUP + concurrent opioid
    p_frailty_BUPC_INJ_1 = rep(1, n_sim),
    p_frailty_BUPC_INJ_2 = rnorm(n_sim, mean = df_frailty["pe", "BUPC_INJ_2"], sd = df_frailty["sd", "BUPC_INJ_2"]),
    p_frailty_BUPC_INJ_3 = rnorm(n_sim, mean = df_frailty["pe", "BUPC_INJ_3"], sd = df_frailty["sd", "BUPC_INJ_3"]),
    # MET
    p_frailty_MET_INJ_1 = rep(1, n_sim),
    p_frailty_MET_INJ_2 = rnorm(n_sim, mean = df_frailty["pe", "MET_INJ_2"], sd = df_frailty["sd", "MET_INJ_2"]),
    p_frailty_MET_INJ_3 = rnorm(n_sim, mean = df_frailty["pe", "MET_INJ_3"], sd = df_frailty["sd", "MET_INJ_3"]),
    # MET + concurrent opioid
    p_frailty_METC_INJ_1 = rep(1, n_sim),
    p_frailty_METC_INJ_2 = rnorm(n_sim, mean = df_frailty["pe", "METC_INJ_2"], sd = df_frailty["sd", "METC_INJ_2"]),
    p_frailty_METC_INJ_3 = rnorm(n_sim, mean = df_frailty["pe", "METC_INJ_3"], sd = df_frailty["sd", "METC_INJ_3"]),
    # ABS
    p_frailty_ABS_INJ_1 = rep(1, n_sim),
    p_frailty_ABS_INJ_2 = rnorm(n_sim, mean = df_frailty["pe", "ABS_INJ_2"], sd = df_frailty["sd", "ABS_INJ_2"]),
    p_frailty_ABS_INJ_3 = rnorm(n_sim, mean = df_frailty["pe", "ABS_INJ_3"], sd = df_frailty["sd", "ABS_INJ_3"]),
    # REL
    p_frailty_REL_INJ_1 = rep(1, n_sim),
    p_frailty_REL_INJ_2 = rnorm(n_sim, mean = df_frailty["pe", "REL_INJ_2"], sd = df_frailty["sd", "REL_INJ_2"]),
    p_frailty_REL_INJ_3 = rnorm(n_sim, mean = df_frailty["pe", "REL_INJ_3"], sd = df_frailty["sd", "REL_INJ_3"]),
    # OD
    p_frailty_OD_INJ_1  = rep(1, n_sim),
    p_frailty_OD_INJ_2  = rnorm(n_sim, mean = df_frailty["pe", "OD_INJ_2"], sd = df_frailty["sd", "OD_INJ_2"]),
    p_frailty_OD_INJ_3  = rnorm(n_sim, mean = df_frailty["pe", "OD_INJ_3"], sd = df_frailty["sd", "OD_INJ_3"]),
    
    # Weibull scale
    # Uniform distribution (low, high)
    # Non-injection
    p_weibull_scale_BUP_NI = runif(n_sim, min =, max = ),
    p_weibull_scale_BUPC_NI = runif(),
    p_weibull_scale_MET_NI = runif(),
    p_weibull_scale_METC_NI = runif(),
    p_weibull_scale_ABS_NI = runif(),
    p_weibull_scale_REL_NI = runif(),
    #p_weibull_scale_OD_NI  = runif(),
    
    # Injection
    p_weibull_scale_BUP_INJ = runif(),
    p_weibull_scale_BUPC_INJ = runif(),
    p_weibull_scale_MET_INJ = runif(),
    p_weibull_scale_METC_INJ = runif(),
    p_weibull_scale_ABS_INJ = runif(),
    p_weibull_scale_REL_INJ = runif(),
    #p_weibull_scale_OD_INJ  = runif(),
    
    # Weibull shape
    # Uniform distribution
    # Non-injection
    p_weibull_shape_BUP_NI = runif(),
    p_weibull_shape_BUPC_NI = runif(),
    p_weibull_shape_MET_NI = runif(),
    p_weibull_shape_METC_NI = runif(),
    p_weibull_shape_ABS_NI = runif(),
    p_weibull_shape_REL_NI = runif(),
    #p_weibull_shape_OD_NI  = runif(),
    
    # Injection
    p_weibull_shape_BUP_INJ = runif(),
    p_weibull_shape_BUPC_INJ = runif(),
    p_weibull_shape_MET_INJ = runif(),
    p_weibull_shape_METC_INJ = runif(),
    p_weibull_shape_ABS_INJ = runif(),
    p_weibull_shape_REL_INJ = runif(),
    #p_weibull_shape_OD_INJ  = runif(),
    
    ### Transition probabilities conditional on leaving (use Dirichlet)
    m_dirichlet_UP = df_UP * n_samp, # weight unconditional matrix by sample size to generate dirichlet PSA
    
    # Non-Injection
    # From BUP
    v_dirichlet_UP_BUP_NI = m_dirichlet_UP["BUP_NI",],
    m_BUP_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_BUP_NI["BUP_NI", "BUPC_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "MET_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "METC_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "ABS_NI"], v_dirichlet_UP_BUP_NI["BUP_NI", "REL_NI"])),
    p_BUP_BUPC_NI = m_BUP_UP_NI[,1], # assign probabilities by place
    p_BUP_MET_NI  = m_BUP_UP_NI[,2],
    p_BUP_METC_NI = m_BUP_UP_NI[,3],
    p_BUP_ABS_NI  = m_BUP_UP_NI[,4],
    p_BUP_REL_NI  = m_BUP_UP_NI[,5],
    
    # From BUPC
    v_dirichlet_UP_BUPC_NI = m_dirichlet_UP["BUPC_NI",],
    m_BUPC_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_BUPC_NI["BUPC_NI", "BUP_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "MET_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "METC_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "ABS_NI"], v_dirichlet_UP_BUPC_NI["BUPC_NI", "REL_NI"])),
    p_BUPC_BUP_NI  = m_BUPC_UP_NI[,1],
    p_BUPC_MET_NI  = m_BUPC_UP_NI[,2],
    p_BUPC_METC_NI = m_BUPC_UP_NI[,3],
    p_BUPC_ABS_NI  = m_BUPC_UP_NI[,4],
    p_BUPC_REL_NI  = m_BUPC_UP_NI[,5],
    
    # From MET
    v_dirichlet_UP_MET_NI = m_dirichlet_UP["MET_NI",],
    m_MET_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_MET_NI["MET_NI", "METC_NI"], v_dirichlet_UP_MET_NI["MET_NI", "BUP_NI"], v_dirichlet_UP_MET_NI["MET_NI", "BUPC_NI"], v_dirichlet_UP_MET_NI["MET_NI", "ABS_NI"], v_dirichlet_UP_MET_NI["MET_NI", "REL_NI"])),
    p_MET_METC_NI = m_MET_UP_NI[,1],
    p_MET_BUP_NI  = m_MET_UP_NI[,2],
    p_MET_BUPC_NI = m_MET_UP_NI[,3],
    p_MET_ABS_NI  = m_MET_UP_NI[,4],
    p_MET_REL_NI  = m_MET_UP_NI[,5],
    
    # From METC
    v_dirichlet_UP_METC_NI = m_dirichlet_UP["METC_NI",],
    m_METC_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_METC_NI["METC_NI", "MET_NI"], v_dirichlet_UP_METC_NI["METC_NI", "BUP_NI"], v_dirichlet_UP_METC_NI["METC_NI", "BUPC_NI"], v_dirichlet_UP_METC_NI["METC_NI", "ABS_NI"], v_dirichlet_UP_METC_NI["METC_NI", "REL_NI"])),
    p_METC_MET_NI  = m_METC_UP_NI[,1],
    p_METC_BUP_NI  = m_METC_UP_NI[,2],
    p_METC_BUPC_NI = m_METC_UP_NI[,3],
    p_METC_ABS_NI  = m_METC_UP_NI[,4],
    p_METC_REL_NI  = m_METC_UP_NI[,5],
    
    # From ABS
    v_dirichlet_UP_ABS_NI = m_dirichlet_UP["ABS_NI",],
    m_ABS_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_ABS_NI["ABS_NI", "MET_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "METC_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "BUP_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "BUPC_NI"], v_dirichlet_UP_ABS_NI["ABS_NI", "REL_NI"])),
    p_ABS_MET_NI  = m_ABS_UP_NI[1,],
    p_ABS_METC_NI = m_ABS_UP_NI[2,],
    p_ABS_BUP_NI  = m_ABS_UP_NI[3,],
    p_ABS_BUPC_NI = m_ABS_UP_NI[4,],
    p_ABS_REL_NI  = m_ABS_UP_NI[5,],
    
    # From REL
    v_dirichlet_UP_REL_NI = m_dirichlet_UP["REL_NI",],
    m_REL_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_REL_NI["REL_NI", "MET_NI"], v_dirichlet_UP_REL_NI["REL_NI", "METC_NI"], v_dirichlet_UP_REL_NI["REL_NI", "BUP_NI"], v_dirichlet_UP_REL_NI["REL_NI", "BUPC_NI"], v_dirichlet_UP_REL_NI["REL_NI", "ABS_NI"])),
    p_REL_MET_NI  = m_REL_UP_NI[1,],
    p_REL_METC_NI = m_REL_UP_NI[2,],
    p_REL_BUP_NI  = m_REL_UP_NI[3,],
    p_REL_BUPC_NI = m_REL_UP_NI[4,],
    p_REL_ABS_NI  = m_REL_UP_NI[5,],
    
    # From OD
    v_dirichlet_UP_OD_NI = m_dirichlet_UP["OD_NI",],
    m_OD_UP_NI = rdirichlet(n_sim, c(v_dirichlet_UP_OD_NI["OD_NI", "MET_NI"], v_dirichlet_UP_OD_NI["OD_NI", "METC_NI"], v_dirichlet_UP_OD_NI["OD_NI", "BUP_NI"], v_dirichlet_UP_OD_NI["OD_NI", "BUPC_NI"], v_dirichlet_UP_OD_NI["OD_NI", "ABS_NI"], v_dirichlet_UP_OD_NI["OD_NI", "REL_NI"])),
    p_OD_MET_NI  = m_OD_UP_NI[1,],
    p_OD_METC_NI = m_OD_UP_NI[2,],
    p_OD_BUP_NI  = m_OD_UP_NI[3,],
    p_OD_BUPC_NI = m_OD_UP_NI[4,],
    p_OD_ABS_NI  = m_OD_UP_NI[5,],
    p_OD_REL_NI  = m_OD_UP_NI[6,],
    
    # Injection
    # From BUP
    v_dirichlet_UP_BUP_INJ = m_dirichlet_UP["BUP_INJ",],
    m_BUP_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_BUP_INJ["BUP_INJ", "BUPC_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "MET_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "METC_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "ABS_INJ"], v_dirichlet_UP_BUP_INJ["BUP_INJ", "REL_INJ"])),
    p_BUP_BUPC_INJ = m_BUP_UP_INJ[,1], # assign probabilities by place
    p_BUP_MET_INJ  = m_BUP_UP_INJ[,2],
    p_BUP_METC_INJ = m_BUP_UP_INJ[,3],
    p_BUP_ABS_INJ  = m_BUP_UP_INJ[,4],
    p_BUP_REL_INJ  = m_BUP_UP_INJ[,5],
    
    # From BUPC
    v_dirichlet_UP_BUPC_INJ = m_dirichlet_UP["BUPC_INJ",],
    m_BUPC_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "BUP_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "MET_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "METC_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "ABS_INJ"], v_dirichlet_UP_BUPC_INJ["BUPC_INJ", "REL_INJ"])),
    p_BUPC_BUP_INJ  = m_BUPC_UP_INJ[,1],
    p_BUPC_MET_INJ  = m_BUPC_UP_INJ[,2],
    p_BUPC_METC_INJ = m_BUPC_UP_INJ[,3],
    p_BUPC_ABS_INJ  = m_BUPC_UP_INJ[,4],
    p_BUPC_REL_INJ  = m_BUPC_UP_INJ[,5],
    
    # From MET
    v_dirichlet_UP_MET_INJ = m_dirichlet_UP["MET_INJ",],
    m_MET_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_MET_INJ["MET_INJ", "METC_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "BUP_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "BUPC_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "ABS_INJ"], v_dirichlet_UP_MET_INJ["MET_INJ", "REL_INJ"])),
    p_MET_METC_INJ = m_MET_UP_INJ[,1],
    p_MET_BUP_INJ  = m_MET_UP_INJ[,2],
    p_MET_BUPC_INJ = m_MET_UP_INJ[,3],
    p_MET_ABS_INJ  = m_MET_UP_INJ[,4],
    p_MET_REL_INJ  = m_MET_UP_INJ[,5],
    
    # From METC
    v_dirichlet_UP_METC_INJ = m_dirichlet_UP["METC_INJ",],
    m_METC_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_METC_INJ["METC_INJ", "MET_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "BUP_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "BUPC_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "ABS_INJ"], v_dirichlet_UP_METC_INJ["METC_INJ", "REL_INJ"])),
    p_METC_MET_INJ  = m_METC_UP_INJ[,1],
    p_METC_BUP_INJ  = m_METC_UP_INJ[,2],
    p_METC_BUPC_INJ = m_METC_UP_INJ[,3],
    p_METC_ABS_INJ  = m_METC_UP_INJ[,4],
    p_METC_REL_INJ  = m_METC_UP_INJ[,5],
    
    # From ABS
    v_dirichlet_UP_ABS_INJ = m_dirichlet_UP["ABS_INJ",],
    m_ABS_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_ABS_INJ["ABS_INJ", "MET_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "METC_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "BUP_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "BUPC_INJ"], v_dirichlet_UP_ABS_INJ["ABS_INJ", "REL_INJ"])),
    p_ABS_MET_INJ  = m_ABS_UP_INJ[1,],
    p_ABS_METC_INJ = m_ABS_UP_INJ[2,],
    p_ABS_BUP_INJ  = m_ABS_UP_INJ[3,],
    p_ABS_BUPC_INJ = m_ABS_UP_INJ[4,],
    p_ABS_REL_INJ  = m_ABS_UP_INJ[5,],
    
    # From REL
    v_dirichlet_UP_REL_INJ = m_dirichlet_UP["REL_INJ",],
    m_REL_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_REL_INJ["REL_INJ", "MET_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "METC_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "BUP_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "BUPC_INJ"], v_dirichlet_UP_REL_INJ["REL_INJ", "ABS_INJ"])),
    p_REL_MET_INJ  = m_REL_UP_INJ[1,],
    p_REL_METC_INJ = m_REL_UP_INJ[2,],
    p_REL_BUP_INJ  = m_REL_UP_INJ[3,],
    p_REL_BUPC_INJ = m_REL_UP_INJ[4,],
    p_REL_ABS_INJ  = m_REL_UP_INJ[5,],
    
    # From OD
    v_dirichlet_UP_OD_INJ = m_dirichlet_UP["OD_INJ",],
    m_OD_UP_INJ = rdirichlet(n_sim, c(v_dirichlet_UP_OD_INJ["OD_INJ", "MET_INJ"], v_dirichlet_UP_OD_INJ["OD_INJ", "METC_INJ"], v_dirichlet_UP_OD_INJ["OD_INJ", "BUP_INJ"], v_dirichlet_UP_OD_INJ["OD_INJ", "BUPC_INJ"], v_dirichlet_UP_OD_INJ["OD_INJ", "ABS_INJ"], v_dirichlet_UP_OD_INJ["OD_INJ", "REL_INJ"])),
    p_OD_MET_INJ  = m_OD_UP_INJ[1,],
    p_OD_METC_INJ = m_OD_UP_INJ[2,],
    p_OD_BUP_INJ  = m_OD_UP_INJ[3,],
    p_OD_BUPC_INJ = m_OD_UP_INJ[4,],
    p_OD_ABS_INJ  = m_OD_UP_INJ[5,],
    p_OD_REL_INJ  = m_OD_UP_INJ[6,],
    
    ### Overdose ###
    n_fent_OD = runif(n_sim, df_overdose["low", "fent_OD_rate"], df_overdose["high", "fent_OD_rate"]),
    p_fent_exp = runif(n_sim, df_overdose["low", "fent_exp_prob"], df_overdose["high", "fent_exp_prob"]),
    p_witness = runif(n_sim, df_overdose["low", "witness_prob"], df_overdose["high", "witness_prob"]),
    p_attended = runif(n_sim, df_overdose["low", "attended_prob"], df_overdose["low", "attended_prob"]),
    p_NX_used = runif(n_sim, df_overdose["low", "NX_prob"], df_overdose["high", "NX_prob"]),
    p_NX_success = runif(n_sim, df_overdose["low", "NX_success_prob"], df_overdose["high", "NX_success_prob"]),
    
    ### HIV seroconversion ###
    # ** DOUBLE-CHECK SEROCONVERSION PSA PROBS IF WE WANT TO HAVE THE EXACT SAME VALUES FOR BUP/BUPC, etc. **
    # CHANGE TO BETA PROBS
    # HIV Seroconversion
    # From negative
    # Non-injection
    p_HIV_NI = rbeta(),
    p_HIV_BUP_NI  = df_hiv["pe", "HIV_BUP_NI"],
    p_HIV_BUPC_NI = df_hiv["pe", "HIV_BUPC_NI"],
    p_HIV_MET_NI  = df_hiv["pe", "HIV_MET_NI"],
    p_HIV_METC_NI = df_hiv["pe", "HIV_METC_NI"],
    p_HIV_REL_NI  = df_hiv["pe", "HIV_REL_NI"],
    p_HIV_ODN_NI  = df_hiv["pe", "HIV_REL_NI"],
    p_HIV_ABS_NI  = df_hiv["pe", "HIV_ABS_NI"],
    # Injection
    p_HIV_BUP_INJ  = df_hiv["pe", "HIV_BUP_INJ"],
    p_HIV_BUPC_INJ = df_hiv["pe", "HIV_BUPC_INJ"],
    p_HIV_MET_INJ  = df_hiv["pe", "HIV_MET_INJ"],
    p_HIV_METC_INJ = df_hiv["pe", "HIV_METC_INJ"],
    p_HIV_REL_INJ  = df_hiv["pe", "HIV_REL_INJ"],
    p_HIV_ODN_INJ  = df_hiv["pe", "HIV_REL_INJ"],
    p_HIV_ABS_INJ  = df_hiv["pe", "HIV_ABS_INJ"],
    
    # Co-infection conditional on HCV
    # Non-injection
    p_HCV_HIV_BUP_NI  = df_hiv["pe", "COI_BUP_NI"],
    p_HCV_HIV_BUPC_NI = df_hiv["pe", "COI_BUPC_NI"],
    p_HCV_HIV_MET_NI  = df_hiv["pe", "COI_MET_NI"],
    p_HCV_HIV_METC_NI = df_hiv["pe", "COI_METC_NI"],
    p_HCV_HIV_REL_NI  = df_hiv["pe", "COI_REL_NI"],
    p_HCV_HIV_ODN_NI  = df_hiv["pe", "COI_REL_NI"],
    p_HCV_HIV_ABS_NI  = df_hiv["pe", "COI_ABS_NI"],
    # Injection
    p_HCV_HIV_BUP_INJ  = df_hiv["pe", "COI_BUP_INJ"],
    p_HCV_HIV_BUPC_INJ = df_hiv["pe", "COI_BUPC_INJ"],
    p_HCV_HIV_MET_INJ  = df_hiv["pe", "COI_MET_INJ"],
    p_HCV_HIV_METC_INJ = df_hiv["pe", "COI_METC_INJ"],
    p_HCV_HIV_REL_INJ  = df_hiv["pe", "COI_REL_INJ"],
    p_HCV_HIV_ODN_INJ  = df_hiv["pe", "COI_REL_INJ"],
    p_HCV_HIV_ABS_INJ  = df_hiv["pe", "COI_ABS_INJ"],
    
    ### Costs ###
    # Treatment Costs
    c_BUP_TX  = rgamma(n_sim, shape = df_costs["shape", "BUP_TX"], scale = df_costs["scale", "BUP_TX"]), # BUP treatment costs - change to normal
    c_MET_TX  = rgamma(n_sim, shape = df_costs["shape", "MET_TX"], scale = df_costs["scale", "MET_TX"]), # MET treatment costs - change to normal
    
    # HRU Costs
    c_BUP_NI_HRU = rgamma(n_sim, shape = df_costs["shape", "BUP_NI_HRU"], scale = df_costs["scale", "BUP_NI_HRU"]), 
    c_MET_NI_HRU = rgamma(n_sim, shape = df_costs["shape", "MET_NI_HRU"], scale = df_costs["scale", "MET_NI_HRU"]), 
    c_ABS_NI_HRU = rgamma(n_sim, shape = df_costs["shape", "ABS_NI_HRU"], scale = df_costs["scale", "ABS_NI_HRU"]), 
    c_REL_NI_HRU = rgamma(n_sim, shape = df_costs["shape", "REL_NI_HRU"], scale = df_costs["scale", "REL_NI_HRU"]), 
    c_OD_NI_HRU  = rgamma(n_sim, shape = df_costs["shape", "OD_NI_HRU"] , scale = df_costs["scale", "OD_NI_HRU"] ), 
    c_BUP_INJ_HRU = rgamma(n_sim, shape = df_costs["shape", "BUP_INJ_HRU"], scale = df_costs["scale", "BUP_INJ_HRU"]), 
    c_MET_INJ_HRU = rgamma(n_sim, shape = df_costs["shape", "MET_INJ_HRU"], scale = df_costs["scale", "MET_INJ_HRU"]), 
    c_ABS_INJ_HRU = rgamma(n_sim, shape = df_costs["shape", "ABS_INJ_HRU"], scale = df_costs["scale", "ABS_INJ_HRU"]), 
    c_REL_INJ_HRU = rgamma(n_sim, shape = df_costs["shape", "REL_INJ_HRU"], scale = df_costs["scale", "REL_INJ_HRU"]), 
    c_OD_INJ_HRU  = rgamma(n_sim, shape = df_costs["shape", "OD_INJ_HRU"] , scale = df_costs["scale", "OD_INJ_HRU"] ),
    
    # HIV Costs
    c_HIV_HRU = rgamma(n_sim, shape = ((df_costs["pe", "HIV_HRU"])^2 / ((df_costs["high", "HIV_HRU"] - df_costs["low", "HIV_HRU"]) / 2 / 1.96)^2), scale = ((df_costs["high", "HIV_HRU"] - df_costs["low", "HIV_HRU"]) / 2 / 1.96)^2 / df_costs["pe", "HIV_HRU"]),
    c_HIV_ART = rgamma(n_sim, shape = ((df_costs["pe", "HIV_ART"])^2 / ((df_costs["high", "HIV_ART"] - df_costs["low", "HIV_ART"]) / 2 / 1.96)^2), scale = ((df_costs["high", "HIV_ART"] - df_costs["low", "HIV_ART"]) / 2 / 1.96)^2 / df_costs["pe", "HIV_ART"]),
    
    # HCV Costs
    c_HCV_HRU = rgamma(n_sim, shape = ((df_costs["pe", "HCV_HRU"])^2 / ((df_costs["high", "HCV_HRU"] - df_costs["low", "HCV_HRU"]) / 2 / 1.96)^2), scale = ((df_costs["high", "HCV_HRU"] - df_costs["low", "HCV_HRU"]) / 2 / 1.96)^2 / df_costs["pe", "HCV_HRU"]),
    c_HCV_DAA = rgamma(n_sim, shape = ((df_costs["pe", "HCV_DAA"])^2 / ((df_costs["high", "HCV_DAA"] - df_costs["low", "HCV_DAA"]) / 2 / 1.96)^2), scale = ((df_costs["high", "HCV_DAA"] - df_costs["low", "HCV_DAA"]) / 2 / 1.96)^2 / df_costs["pe", "HCV_DAA"]),

    # Crime Costs
    c_BUP_NI_crime = rgamma(n_sim, shape = ((df_costs["pe", "BUP"])^2 / ((df_costs["high", "BUP"] - df_costs["low", "BUP"]) / 2 / 1.96)^2), scale = ((df_costs["high", "BUP"] - df_costs["low", "BUP"]) / 2 / 1.96)^2 / df_costs["pe", "BUP"]),
    c_BUP_INJ_crime = c_BUP_NI_crime,
    c_BUPC_NI_crime = rgamma(n_sim, shape = ((df_costs["pe", "BUPC"])^2 / ((df_costs["high", "BUPC"] - df_costs["low", "BUPC"]) / 2 / 1.96)^2), scale = ((df_costs["high", "BUPC"] - df_costs["low", "BUPC"]) / 2 / 1.96)^2 / df_costs["pe", "BUPC"]),
    c_BUPC_INJ_crime = c_BUPC_NI_crime,
    c_MET_NI_crime = rgamma(n_sim, shape = ((df_costs["pe", "MET"])^2 / ((df_costs["high", "MET"] - df_costs["low", "MET"]) / 2 / 1.96)^2), scale = ((df_costs["high", "MET"] - df_costs["low", "MET"]) / 2 / 1.96)^2 / df_costs["pe", "MET"]),
    c_MET_INJ_crime = c_MET_NI_crime,
    c_METC_NI_crime = rgamma(n_sim, shape = ((df_costs["pe", "METC"])^2 / ((df_costs["high", "METC"] - df_costs["low", "METC"]) / 2 / 1.96)^2), scale = ((df_costs["high", "METC"] - df_costs["low", "METC"]) / 2 / 1.96)^2 / df_costs["pe", "METC"]),
    c_METC_INJ_crime = c_METC_NI_crime,
    c_REL_NI_crime = rgamma(n_sim, shape = ((df_costs["pe", "REL"])^2 / ((df_costs["high", "REL"] - df_costs["low", "REL"]) / 2 / 1.96)^2), scale = ((df_costs["high", "REL"] - df_costs["low", "REL"]) / 2 / 1.96)^2 / df_costs["pe", "REL"]),
    c_REL_INJ_crime = c_REL_NI_crime,
    c_ODN_NI_crime = c_REL_NI_crime,
    c_ODN_INJ_crime = c_REL_NI_crime,
    #c_ODF_NI_crime = rep(0, n_sim),
    
    ### Utilities ###
    # Consider having equivalent QALYs between injection and non-injection
    # HIV/HCV negative
    u_BUP_NI_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "BUP_NI_NEG"], sd = df_qalys["sd", "BUP_NI_NEG"], b = 1),
    u_BUPC_NI_NEG = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "BUPC_NI_NEG"], sd = df_qalys["sd", "BUPC_NI_NEG"], b = 1),
    u_MET_NI_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "MET_NI_NEG"], sd = df_qalys["sd", "MET_NI_NEG"], b = 1),
    u_METC_NI_NEG = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "METC_NI_NEG"], sd = df_qalys["sd", "METC_NI_NEG"], b = 1),
    u_REL_NI_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "REL_NI_NEG"], sd = df_qalys["sd", "REL_NI_NEG"], b = 1),
    u_ODN_NI_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "ODN_NI_NEG"], sd = df_qalys["sd", "ODN_NI_NEG"], b = 1),
    U_ODF_NI_NEG  = rep(0, n_sim),
    u_ABS_NI_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "ABS_NI_NEG"], sd = df_qalys["sd", "ABS_NI_NEG"], b = 1),
    
    u_BUP_INJ_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "BUP_INJ_NEG"], sd = df_qalys["sd", "BUP_INJ_NEG"], b = 1),
    u_BUPC_INJ_NEG = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "BUPC_INJ_NEG"], sd = df_qalys["sd", "BUPC_INJ_NEG"], b = 1),
    u_MET_INJ_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "MET_INJ_NEG"], sd = df_qalys["sd", "MET_INJ_NEG"], b = 1),
    u_METC_INJ_NEG = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "METC_INJ_NEG"], sd = df_qalys["sd", "METC_INJ_NEG"], b = 1),
    u_REL_INJ_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "REL_INJ_NEG"], sd = df_qalys["sd", "REL_INJ_NEG"], b = 1),
    u_ODN_INJ_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "ODN_INJ_NEG"], sd = df_qalys["sd", "ODN_INJ_NEG"], b = 1),
    U_ODF_INJ_NEG  = rep(0, n_sim),
    u_ABS_INJ_NEG  = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "ABS_INJ_NEG"], sd = df_qalys["sd", "ABS_INJ_NEG"], b = 1),

    # Consider beta (other?) distributions for these
    u_HIV_mult = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "HIV_mult"], sd = df_qalys["sd", "HIV_mult"], b = 1), # HIV multiplier for negative states
    u_HCV_mult = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "HCV_mult"], sd = df_qalys["sd", "HCV_mult"], b = 1), # HCV multiplier for negative states
    U_COI_mult = truncnorm::rtruncnorm(n_sim, mean = df_qalys["pe", "COI_mult"], sd = df_qalys["sd", "COI_mult"], b = 1))
  
  return(df_psa_params)
}

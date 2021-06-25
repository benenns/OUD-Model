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
generate_psa_params <- function(n_sim = n_sim, seed = seed,
                                file.death_hr = NULL,
                                file.frailty = NULL,
                                file.weibull_scale = NULL,
                                file.weibull_shape = NULL,
                                file.unconditional = NULL,
                                file.sero = NULL,
                                file.costs = NULL,
                                file.crime_costs = NULL,
                                file.qalys = NULL){
  
  #Load files with parameter distribution values
  df_death_hr <- read.csv(file = file.death_hr, row.names = 1, header = TRUE) # Mortality hazard ratios
  df_frailty <- read.csv(file = file.frailty, row.names = 1, header = TRUE) # Episode frailty params
  df_weibull_scale <- read.csv(file = file.weibull_shape, row.names = 1, header = TRUE) # Weibull scale params
  df_weibull_shape <- read.csv(file = file.weibull_scale, row.names = 1, header = TRUE) # Weibull shape params
  df_UP <- read.csv(file = file.unconditional, row.names = 1, header = TRUE) # Unconditional transition probs
  df_sero <- read.csv(file = file.sero, row.names = 1, header = TRUE) # Seroconversion probs
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
    # HIV-negative
    hr_BUP_NI  = rnorm(n_sim, mean = df_death_hr["pe", "BUP_NI"] , sd = df_death_hr["sd", "BUP_NI"] ),
    hr_BUPC_NI  = rnorm(n_sim, mean = df_death_hr["pe", "BUPC_NI"] , sd = df_death_hr["sd", "BUPC_NI"] ),
    
    hr_MET_NI  = rnorm(n_sim, mean = df_death_hr["pe", "MET_NI"] , sd = df_death_hr["sd", "MET_NI"] ),
    hr_METC_NI  = rnorm(n_sim, mean = df_death_hr["pe", "METC_NI"] , sd = df_death_hr["sd", "METC_NI"] ),
    
    hr_REL_NI  = rnorm(n_sim, mean = df_death_hr["pe", "REL_NI"] , sd = df_death_hr["sd", "REL_NI"] ),
    hr_ODN_NI  = rnorm(n_sim, mean = df_death_hr["pe", "ODN_NI"] , sd = df_death_hr["sd", "ODN_NI"] ),
    
    hr_ABS_NI  = rnorm(n_sim, mean = df_death_hr["pe", "ABS_NI"] , sd = df_death_hr["sd", "ABS_NI"] ),
    hr_HIV_NI  = rnorm(n_sim, mean = df_death_hr["pe", "HIV_NI"] , sd = df_death_hr["sd", "HIV_NI"] ),
    hr_HCV_NI  = rnorm(n_sim, mean = df_death_hr["pe", "HCV_NI"] , sd = df_death_hr["sd", "HCV_NI"] ),
    hr_COI_NI  = rnorm(n_sim, mean = df_death_hr["pe", "COI_NI"] , sd = df_death_hr["sd", "COI_NI"] ),
    
    hr_BUP_INJ  = rnorm(n_sim, mean = df_death_hr["pe", "BUP_INJ"] , sd = df_death_hr["sd", "BUP_INJ"] ),
    hr_MET_INJ  = rnorm(n_sim, mean = df_death_hr["pe", "MET_INJ"] , sd = df_death_hr["sd", "MET_INJ"] ),
    hr_REL_INJ  = rnorm(n_sim, mean = df_death_hr["pe", "REL_INJ"] , sd = df_death_hr["sd", "REL_INJ"] ),
    hr_ODN_INJ  = rnorm(n_sim, mean = df_death_hr["pe", "ODN_INJ"] , sd = df_death_hr["sd", "ODN_INJ"] ),
    hr_ABS_INJ  = rnorm(n_sim, mean = df_death_hr["pe", "ABS_INJ"] , sd = df_death_hr["sd", "ABS_INJ"] ),
    hr_HIV_INJ  = rnorm(n_sim, mean = df_death_hr["pe", "HIV_INJ"] , sd = df_death_hr["sd", "HIV_INJ"] ),
  
    # Frailty terms for successive episodes
    # Non-injection
    p_frailty_BUP_NI_1 = rep(1, n_sim),
    p_frailty_BUP_NI_2 = rnorm(n_sim, mean = par_frailty_BUP_NI_2_mean, sd = par_frailty_BUP_NI_2_sd),
    p_frailty_BUP_NI_3 = rnorm(n_sim, mean = par_frailty_BUP_NI_3_mean, sd = par_frailty_BUP_NI_3_sd),
    p_frailty_MET_NI_1 = rep(1, n_sim),
    p_frailty_MET_NI_2 = rnorm(n_sim, mean = par_frailty_MET_NI_2_mean, sd = par_frailty_MET_NI_2_sd),
    p_frailty_MET_NI_3 = rnorm(n_sim, mean = par_frailty_MET_NI_3_mean, sd = par_frailty_MET_NI_3_sd),
    p_frailty_ABS_NI_1 = rep(1, n_sim),
    p_frailty_ABS_NI_2 = rnorm(n_sim, mean = par_frailty_ABS_NI_2_mean, sd = par_frailty_ABS_NI_2_sd),
    p_frailty_ABS_NI_3 = rnorm(n_sim, mean = par_frailty_ABS_NI_3_mean, sd = par_frailty_ABS_NI_3_sd),
    p_frailty_REL_NI_1 = rep(1, n_sim),
    p_frailty_REL_NI_2 = rnorm(n_sim, mean = par_frailty_REL_NI_2_mean, sd = par_frailty_REL_NI_2_sd),
    p_frailty_REL_NI_3 = rnorm(n_sim, mean = par_frailty_REL_NI_3_mean, sd = par_frailty_REL_NI_3_sd),
    p_frailty_OD_NI_1  = rep(1, n_sim),
    p_frailty_OD_NI_2  = rnorm(n_sim, mean = par_frailty_OD_NI_2_mean, sd = par_frailty_OD_NI_2_sd),
    p_frailty_OD_NI_3  = rnorm(n_sim, mean = par_frailty_OD_NI_3_mean, sd = par_frailty_OD_NI_3_sd),
    
    # Injection
    p_frailty_BUP_INJ_1 = rep(1, n_sim),
    p_frailty_BUP_INJ_2 = rnorm(n_sim, mean = par_frailty_BUP_INJ_2_mean, sd = par_frailty_BUP_INJ_2_sd),
    p_frailty_BUP_INJ_3 = rnorm(n_sim, mean = par_frailty_BUP_INJ_3_mean, sd = par_frailty_BUP_INJ_3_sd),
    p_frailty_MET_INJ_1 = rep(1, n_sim),
    p_frailty_MET_INJ_2 = rnorm(n_sim, mean = par_frailty_MET_INJ_2_mean, sd = par_frailty_MET_INJ_2_sd),
    p_frailty_MET_INJ_3 = rnorm(n_sim, mean = par_frailty_MET_INJ_3_mean, sd = par_frailty_MET_INJ_3_sd),
    p_frailty_ABS_INJ_1 = rep(1, n_sim),
    p_frailty_ABS_INJ_2 = rnorm(n_sim, mean = par_frailty_ABS_INJ_2_mean, sd = par_frailty_ABS_INJ_2_sd),
    p_frailty_ABS_INJ_3 = rnorm(n_sim, mean = par_frailty_ABS_INJ_3_mean, sd = par_frailty_ABS_INJ_3_sd),
    p_frailty_REL_INJ_1 = rep(1, n_sim),
    p_frailty_REL_INJ_2 = rnorm(n_sim, mean = par_frailty_REL_INJ_2_mean, sd = par_frailty_REL_INJ_2_sd), 
    p_frailty_REL_INJ_3 = rnorm(n_sim, mean = par_frailty_REL_INJ_3_mean, sd = par_frailty_REL_INJ_3_sd),
    p_frailty_OD_INJ_1  = rep(1, n_sim),
    p_frailty_OD_INJ_2  = rnorm(n_sim, mean = par_frailty_OD_INJ_2_mean, sd = par_frailty_OD_INJ_2_sd),
    p_frailty_OD_INJ_3  = rnorm(n_sim, mean = par_frailty_OD_INJ_3_mean, sd = par_frailty_OD_INJ_3_sd),
    
    # Weibull scale  
    p_weibull_scale_BUP_NI = runif(),
    p_weibull_scale_MET_NI = runif(),
    p_weibull_scale_ABS_NI = runif(),
    p_weibull_scale_REL_NI = runif(),
    #p_weibull_scale_OD_NI  = runif(),
    
    p_weibull_scale_BUP_INJ = runif(),
    p_weibull_scale_MET_INJ = runif(),
    p_weibull_scale_ABS_INJ = runif(),
    p_weibull_scale_REL_INJ = runif(),
    #p_weibull_scale_OD_INJ  = runif(),
    
    # Weibull shape
    p_weibull_shape_BUP_NI = runif(),
    p_weibull_shape_MET_NI = runif(),
    p_weibull_shape_ABS_NI = runif(),
    p_weibull_shape_REL_NI = runif(),
    #p_weibull_shape_OD_NI  = runif(),
    
    p_weibull_shape_BUP_INJ = runif(),
    p_weibull_shape_MET_INJ = runif(),
    p_weibull_shape_ABS_INJ = runif(),
    p_weibull_shape_REL_INJ = runif(),
    #p_weibull_shape_OD_INJ  = runif(),
    
    ### Transition probabilities conditional on leaving (use Dirichlet)
    # Non-Injection
    # From BUP
    p_BUP_MET_NI  = rdirichlet(),
    p_BUP_ABS_NI  = rdirichlet(),
    p_BUP_OD_NI   = rdirichlet(),
    # From MET
    p_MET_ABS_NI  = rdirichlet(),
    p_MET_OD_NI   = rdirichlet(),
    # From ABS
    p_ABS_REL1_NI = rdirichlet(),
    p_ABS_OD_NI   = rdirichlet(),
    # From REL
    p_REL_MET1_NI = rdirichlet(),
    p_REL_BUP1_NI = rdirichlet(),
    p_REL_ABS_NI  = rdirichlet(),
    p_REL_OD_NI   = rdirichlet(),
    # From OD
    p_OD_MET1_NI = rdirichlet(),
    p_OD_BUP1_NI = rdirichlet(),
    p_OD_ABS_NI  = rdirichlet(),
    p_OD_REL1_NI = rdirichlet(),
    
    # Injection
    # From BUP1
    p_BUP1_MET1_INJ = rdirichlet(),
    p_BUP1_ABS_INJ  = rdirichlet(),
    p_BUP1_REL1_INJ = rdirichlet(),
    p_BUP1_OD_INJ   = rdirichlet(),
    # From BUP
    p_BUP_MET1_INJ = rdirichlet(),
    p_BUP_ABS_INJ  = rdirichlet(),
    p_BUP_REL1_INJ = rdirichlet(),
    p_BUP_OD_INJ   = rdirichlet(),
    # From MET1
    p_MET1_BUP1_INJ = rdirichlet(),
    p_MET1_ABS_INJ  = rdirichlet(),
    p_MET1_REL1_INJ = rdirichlet(),
    p_MET1_OD_INJ   = rdirichlet(),
    # From MET
    p_MET_BUP1_INJ = rdirichlet(),
    p_MET_ABS_INJ  = rdirichlet(),
    p_MET_REL1_INJ = rdirichlet(),
    p_MET_OD_INJ   = rdirichlet(),
    # From ABS
    p_ABS_REL1_INJ = rdirichlet(),
    p_ABS_OD_INJ   = rdirichlet(),
    # From REL1
    p_REL1_MET1_INJ = rdirichlet(),
    p_REL1_BUP1_INJ = rdirichlet(),
    p_REL1_ABS_INJ  = rdirichlet(),
    p_REL1_OD_INJ   = rdirichlet(),
    # From REL
    p_REL_MET1_INJ = rdirichlet(),
    p_REL_BUP1_INJ = rdirichlet(),
    p_REL_ABS_INJ  = rdirichlet(),
    p_REL_OD_INJ   = rdirichlet(),
    # From OD
    p_OD_MET1_INJ = rdirichlet(),
    p_OD_BUP1_INJ = rdirichlet(),
    p_OD_ABS_INJ  = rdirichlet(),
    p_OD_REL1_INJ = rdirichlet(),
    

    
    ### HIV seroconversion ###
    # Non-injection
    p_sero_BUP1_NI = runif(n_sim, df_sero["low", "BUP_NI"], df_sero["high", "BUP_NI"]),
    p_sero_BUP_NI  = runif(n_sim, df_sero["low", "BUP_NI"], df_sero["high", "BUP_NI"]),
    p_sero_MET1_NI = runif(n_sim, df_sero["low", "MET_NI"], df_sero["high", "MET_NI"]),
    p_sero_MET_NI  = runif(n_sim, df_sero["low", "MET_NI"], df_sero["high", "MET_NI"]),
    p_sero_REL1_NI = runif(n_sim, df_sero["low", "REL_NI"], df_sero["high", "REL_NI"]),
    p_sero_REL_NI  = runif(n_sim, df_sero["low", "REL_NI"], df_sero["high", "REL_NI"]),
    p_sero_OD_NI   = runif(n_sim, df_sero["low", "REL_NI"], df_sero["high", "REL_NI"]),
    p_sero_ABS_NI  = runif(n_sim, df_sero["low", "ABS_NI"], df_sero["high", "ABS_NI"]),
    # Injection
    p_sero_BUP1_INJ = runif(n_sim, df_sero["low", "BUP_INJ"], df_sero["high", "BUP_INJ"]),
    p_sero_BUP_INJ  = runif(n_sim, df_sero["low", "BUP_INJ"], df_sero["high", "BUP_INJ"]),
    p_sero_MET1_INJ = runif(n_sim, df_sero["low", "MET_INJ"], df_sero["high", "MET_INJ"]),
    p_sero_MET_INJ  = runif(n_sim, df_sero["low", "MET_INJ"], df_sero["high", "MET_INJ"]),
    p_sero_REL1_INJ = runif(n_sim, df_sero["low", "REL_INJ"], df_sero["high", "REL_INJ"]),
    p_sero_REL_INJ  = runif(n_sim, df_sero["low", "REL_INJ"], df_sero["high", "REL_INJ"]),
    p_sero_OD_INJ   = runif(n_sim, df_sero["low", "REL_INJ"], df_sero["high", "REL_INJ"]),
    p_sero_ABS_INJ  = runif(n_sim, df_sero["low", "ABS_INJ"], df_sero["high", "ABS_INJ"]),
    
    ### Costs ###
    # Treatment Costs
    c_BUP_TX  = rgamma(n_sim, shape = df_costs["shape", "BUP_TX"], scale = df_costs["scale", "BUP_TX"]), # BUP treatment costs
    c_MET_TX  = rgamma(n_sim, shape = df_costs["shape", "MET_TX"], scale = df_costs["scale", "MET_TX"]), # MET treatment costs
    
    # HRU Costs
    # Modify if age-specific
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
    c_HIV = rgamma(n_sim, shape = df_costs["shape", "HIV_HRU"], scale = df_costs["scale", "HIV_HRU"]),
    c_ART = rgamma(n_sim, shape = df_costs["shape", "HIV_ART"], scale = df_costs["scale", "HIV_ART"]),
    
    # Crime Costs
    # Age-specific
    

    
    
    
    ### Utilities ###
    u_H   = truncnorm::rtruncnorm(n_sim, mean =    1, sd = 0.01, b = 1), # utility when healthy
    u_S1  = truncnorm::rtruncnorm(n_sim, mean = 0.75, sd = 0.02, b = 1), # utility when sick
    u_S2  = truncnorm::rtruncnorm(n_sim, mean = 0.50, sd = 0.03, b = 1), # utility when sicker
    u_D   = 0                                               , # utility when dead
    u_Trt = truncnorm::rtruncnorm(n_sim, mean = 0.95, sd = 0.02, b = 1)  # utility when being treated
  )
  return(df_psa_params)
}

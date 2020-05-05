#' Generate PSA dataset of CEA parameters
#'
#' \code{generate_psa_params} generates PSA input dataset by sampling decision 
#' model parameters from their distributions. The sample of the calibrated
#' parameters is a draw from their posterior distribution obtained with the
#' IMIS algorithm.
#' @param n_sim Number of PSA samples.
#' @param seed Seed for reproducibility of Monte Carlo sampling.
#' @return 
#' A data frame with \code{n_sim} rows and 15 columns of parameters for PSA. 
#' Each row is a parameter set sampled from distributions that characterize 
#' their uncertainty
#' @examples 
#' generate_psa_params()
#' @export
generate_psa_params <- function(n_sim = 2000, seed = 3730687){ # User defined
  ## Load calibrated parameters
  n_sim <- nrow(m_calib_post)
  set_seed <- seed
  df_psa_params <- data.frame(
    ### Calibrated parameters
    m_calib_post,
  
    # Frailty terms for successive episodes
    # Consider "truncnorm" if wanting to restrict frailty terms to >1 for EP1
    # Non-injection
    p_frailty_BUP_NI_1 <- rep(1, n_sim),
    #p_frailty_BUP_NI_2 <- truncnorm::rtruncnorm(n_sim, mean = par_frailty_BUP_NI_2$mean, sd = par_frailty_BUP_NI_2$sd, a = 1), # Disallow frailty <1
    p_frailty_BUP_NI_2 <- rnorm(n_sim, mean = par_frailty_BUP_NI_2_mean, sd = par_frailty_BUP_NI_2_sd),
    p_frailty_BUP_NI_3 <- rnorm(n_sim, mean = par_frailty_BUP_NI_3_mean, sd = par_frailty_BUP_NI_3_sd),
    p_frailty_MET_NI_1 <- rep(1, n_sim),
    p_frailty_MET_NI_2 <- rnorm(n_sim, mean = par_frailty_MET_NI_2_mean, sd = par_frailty_MET_NI_2_sd),
    p_frailty_MET_NI_3 <- rnorm(n_sim, mean = par_frailty_MET_NI_3_mean, sd = par_frailty_MET_NI_3_sd),
    p_frailty_ABS_NI_1 <- rep(1, n_sim),
    p_frailty_ABS_NI_2 <- rnorm(n_sim, mean = par_frailty_ABS_NI_2_mean, sd = par_frailty_ABS_NI_2_sd),
    p_frailty_ABS_NI_3 <- rnorm(n_sim, mean = par_frailty_ABS_NI_3_mean, sd = par_frailty_ABS_NI_3_sd),
    p_frailty_REL_NI_1 <- rep(1, n_sim),
    p_frailty_REL_NI_2 <- rnorm(n_sim, mean = par_frailty_REL_NI_2_mean, sd = par_frailty_REL_NI_2_sd),
    p_frailty_REL_NI_3 <- rnorm(n_sim, mean = par_frailty_REL_NI_3_mean, sd = par_frailty_REL_NI_3_sd),
    p_frailty_OD_NI_1  <- rep(1, n_sim),
    p_frailty_OD_NI_2  <- rnorm(n_sim, mean = par_frailty_OD_NI_2_mean, sd = par_frailty_OD_NI_2_sd),
    p_frailty_OD_NI_3  <- rnorm(n_sim, mean = par_frailty_OD_NI_3_mean, sd = par_frailty_OD_NI_3_sd),
    
    # Injection
    p_frailty_BUP_INJ_1 <- rep(1, n_sim),
    p_frailty_BUP_INJ_2 <- rnorm(n_sim, mean = par_frailty_BUP_INJ_2_mean, sd = par_frailty_BUP_INJ_2_sd),
    p_frailty_BUP_INJ_3 <- rnorm(n_sim, mean = par_frailty_BUP_INJ_3_mean, sd = par_frailty_BUP_INJ_3_sd),
    p_frailty_MET_INJ_1 <- rep(1, n_sim),
    p_frailty_MET_INJ_2 <- rnorm(n_sim, mean = par_frailty_MET_INJ_2_mean, sd = par_frailty_MET_INJ_2_sd),
    p_frailty_MET_INJ_3 <- rnorm(n_sim, mean = par_frailty_MET_INJ_3_mean, sd = par_frailty_MET_INJ_3_sd),
    p_frailty_ABS_INJ_1 <- rep(1, n_sim),
    p_frailty_ABS_INJ_2 <- rnorm(n_sim, mean = par_frailty_ABS_INJ_2_mean, sd = par_frailty_ABS_INJ_2_sd),
    p_frailty_ABS_INJ_3 <- rnorm(n_sim, mean = par_frailty_ABS_INJ_3_mean, sd = par_frailty_ABS_INJ_3_sd),
    p_frailty_REL_INJ_1 <- rep(1, n_sim),
    p_frailty_REL_INJ_2 <- rnorm(n_sim, mean = par_frailty_REL_INJ_2_mean, sd = par_frailty_REL_INJ_2_sd), 
    p_frailty_REL_INJ_3 <- rnorm(n_sim, mean = par_frailty_REL_INJ_3_mean, sd = par_frailty_REL_INJ_3_sd),
    p_frailty_OD_INJ_1  <- rep(1, n_sim),
    p_frailty_OD_INJ_2  <- rnorm(n_sim, mean = par_frailty_OD_INJ_2_mean, sd = par_frailty_OD_INJ_2_sd),
    p_frailty_OD_INJ_3  <- rnorm(n_sim, mean = par_frailty_OD_INJ_3_mean, sd = par_frailty_OD_INJ_3_sd),
    
    # Weibull scale  
    p_weibull_scale_BUP_NI <- runif(),
    p_weibull_scale_MET_NI <- runif(),
    p_weibull_scale_ABS_NI <- runif(),
    p_weibull_scale_REL_NI <- runif(),
    p_weibull_scale_OD_NI  <- runif(),
    
    p_weibull_scale_BUP_INJ <- runif(),
    p_weibull_scale_MET_INJ <- runif(),
    p_weibull_scale_ABS_INJ <- runif(),
    p_weibull_scale_REL_INJ <- runif(),
    p_weibull_scale_OD_INJ  <- runif(),
    
    # Weibull shape
    p_weibull_shape_BUP_NI <- runif(),
    p_weibull_shape_MET_NI <- runif(),
    p_weibull_shape_ABS_NI <- runif(),
    p_weibull_shape_REL_NI <- runif(),
    p_weibull_shape_OD_NI  <- runif(),
    
    p_weibull_shape_BUP_INJ <- runif(),
    p_weibull_shape_MET_INJ <- runif(),
    p_weibull_shape_ABS_INJ <- runif(),
    p_weibull_shape_REL_INJ <- runif(),
    p_weibull_shape_OD_INJ  <- runif(),
    
    ### Transition probabilities conditional on leaving (use Dirichlet)
    # Sample code
    #p_HS1   = rbeta(n_sim, 30, 170),        # probability to become sick when healthy
    #p_S1H   = rbeta(n_sim, 60, 60) ,        # probability to become healthy when sick
    
    ######## Non-Injection #########
    # From BUP1
    p_BUP1_MET1_NI <- rdirichlet(),
    p_BUP1_ABS_NI  <- rdirichlet(),
    p_BUP1_REL1_NI <- rdirichlet(),
    p_BUP1_OD_NI   <- rdirichlet(),
    # From BUP
    p_BUP_MET1_NI <- rdirichlet(),
    p_BUP_ABS_NI  <- rdirichlet(),
    p_BUP_REL1_NI <- rdirichlet(),
    p_BUP_OD_NI   <- rdirichlet(),
    # From MET1
    p_MET1_BUP1_NI <- rdirichlet(),
    p_MET1_ABS_NI  <- rdirichlet(),
    p_MET1_REL1_NI <- rdirichlet(),
    p_MET1_OD_NI   <- rdirichlet(),
    # From MET
    p_MET_BUP1_NI <- rdirichlet(),
    p_MET_ABS_NI  <- rdirichlet(),
    p_MET_REL1_NI <- rdirichlet(),
    p_MET_OD_NI   <- rdirichlet(),
    # From ABS
    p_ABS_REL1_NI <- rdirichlet(),
    p_ABS_OD_NI   <- rdirichlet(),
    # From REL1
    p_REL1_MET1_NI <- rdirichlet(),
    p_REL1_BUP1_NI <- rdirichlet(),
    p_REL1_ABS_NI  <- rdirichlet(),
    p_REL1_OD_NI   <- rdirichlet(),
    # From REL
    p_REL_MET1_NI <- rdirichlet(),
    p_REL_BUP1_NI <- rdirichlet(),
    p_REL_ABS_NI  <- rdirichlet(),
    p_REL_OD_NI   <- rdirichlet(),
    # From OD
    p_OD_MET1_NI <- rdirichlet(),
    p_OD_BUP1_NI <- rdirichlet(),
    p_OD_ABS_NI  <- rdirichlet(),
    p_OD_REL1_NI <- rdirichlet(),
    
    ######## Injection ##########
    # From BUP1
    p_BUP1_MET1_INJ <- rdirichlet(),
    p_BUP1_ABS_INJ  <- rdirichlet(),
    p_BUP1_REL1_INJ <- rdirichlet(),
    p_BUP1_OD_INJ   <- rdirichlet(),
    # From BUP
    p_BUP_MET1_INJ <- rdirichlet(),
    p_BUP_ABS_INJ  <- rdirichlet(),
    p_BUP_REL1_INJ <- rdirichlet(),
    p_BUP_OD_INJ   <- rdirichlet(),
    # From MET1
    p_MET1_BUP1_INJ <- rdirichlet(),
    p_MET1_ABS_INJ  <- rdirichlet(),
    p_MET1_REL1_INJ <- rdirichlet(),
    p_MET1_OD_INJ   <- rdirichlet(),
    # From MET
    p_MET_BUP1_INJ <- rdirichlet(),
    p_MET_ABS_INJ  <- rdirichlet(),
    p_MET_REL1_INJ <- rdirichlet(),
    p_MET_OD_INJ   <- rdirichlet(),
    # From ABS
    p_ABS_REL1_INJ <- rdirichlet(),
    p_ABS_OD_INJ   <- rdirichlet(),
    # From REL1
    p_REL1_MET1_INJ <- rdirichlet(),
    p_REL1_BUP1_INJ <- rdirichlet(),
    p_REL1_ABS_INJ  <- rdirichlet(),
    p_REL1_OD_INJ   <- rdirichlet(),
    # From REL
    p_REL_MET1_INJ <- rdirichlet(),
    p_REL_BUP1_INJ <- rdirichlet(),
    p_REL_ABS_INJ  <- rdirichlet(),
    p_REL_OD_INJ   <- rdirichlet(),
    # From OD
    p_OD_MET1_INJ <- rdirichlet(),
    p_OD_BUP1_INJ <- rdirichlet(),
    p_OD_ABS_INJ  <- rdirichlet(),
    p_OD_REL1_INJ <- rdirichlet(),
    
    
    
    
    
    ### State rewards
    ## Costs
    c_BUP_INJ_TX  = rgamma(n_sim, shape = , scale = ), #
    c_BUP_INJ_HRU = rgamma(n_sim, shape = , scale = ), #
    
    
    c_H   = rgamma(n_sim, shape = 100, scale = 20)    , # cost of remaining one cycle in state H
    c_S1  = rgamma(n_sim, shape = 177.8, scale = 22.5), # cost of remaining one cycle in state S1
    c_S2  = rgamma(n_sim, shape = 225, scale = 66.7)  , # cost of remaining one cycle in state S2
    c_Trt = rgamma(n_sim, shape = 73.5, scale = 163.3), # cost of treatment (per cycle)
    c_D   = 0                                         , # cost of being in the death state
    ## Utilities
    u_H   = truncnorm::rtruncnorm(n_sim, mean =    1, sd = 0.01, b = 1), # utility when healthy
    u_S1  = truncnorm::rtruncnorm(n_sim, mean = 0.75, sd = 0.02, b = 1), # utility when sick
    u_S2  = truncnorm::rtruncnorm(n_sim, mean = 0.50, sd = 0.03, b = 1), # utility when sicker
    u_D   = 0                                               , # utility when dead
    u_Trt = truncnorm::rtruncnorm(n_sim, mean = 0.95, sd = 0.02, b = 1)  # utility when being treated
  )
  return(df_psa_params)
}

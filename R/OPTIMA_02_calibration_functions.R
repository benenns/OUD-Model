#' Generate model outputs for calibration from a parameter set
#'
#' \code{calibration_out} computes model outputs to be used for calibration 
#' routines.
#'
#' @param v_params_calib is a vector of parameters that need to be calibrated.   
#' @param l_params_all is a list with all parameters of the decision model.
#' @return 
#' A list with all cause deaths, and non-fatal overdoses.
#' @export
calibration_out <- function(v_params_calib, 
                            l_params_all){
  
  #Load file with calibration parameters and priors
  #df_overdose <- read.csv(file = file.overdose, row.names = 1, header = TRUE) # Overdose-fentanyl parameters
  #df_targets <- read.csv(file = file.targets, row.names = 1, header = TRUE) # Overdose and fatal overdose targets
  
  # Substitute values of calibrated parameters in base-case with 
  # calibrated values
  l_params_all <- update_param_list(l_params_all = l_params_all, params_updated = v_params_calib)
  
  # Run model with updated calibrated parameters
  l_out_markov <- markov_model(l_params_all = l_params_all)
  
  #### Epidemiological Output ####
  ### Overdose deaths ###
  v_ODF <- l_out_markov$m_M_agg_trace[, "ODF"] # cumulative deaths at time i

  ### Non-fatal overdoses ###
  v_ODN <- l_out_markov$m_M_agg_trace[, "ODN"] # cumulative non-fatal overdoses at time i
  
  #### Return Output ####
  l_out <- list(fatal_overdose = v_ODF[c(l_cali_targets$ODF$Time[1], l_cali_targets$ODF$Time[2], l_cali_targets$ODF$Time[3])], # cumulative deaths at t1, t2 , t3 time periods (for yearly deaths: (i + 12) - i where i = first month of year, 1 + 12 = last month)
                overdose = v_ODN[c(l_cali_targets$ODN$Time[1], l_cali_targets$ODN$Time[2], l_cali_targets$ODN$Time[3])])
  return(l_out)
}

#' Sample from prior distributions of calibrated parameters (THIS ASSUMES UNIFORM PRIORS)
sample.prior <- function(n_samp,
                         v_param_names = c("p_BUP_OD_NI", 
                                           "p_MET_OD_NI", 
                                           "p_REL_OD_NI", 
                                           "p_ABS_OD_NI", 
                                           "p_BUP_OD_INJ", 
                                           "p_MET_OD_INJ", 
                                           "p_REL_OD_INJ", 
                                           "p_ABS_OD_INJ"),
                         v_lb = c(p_BUP_OD_NI_lb = l_params_all$p_BUP_OD_NI_lb, 
                                  p_MET_OD_NI_lb = l_params_all$p_MET_OD_NI_lb, 
                                  p_REL_OD_NI_lb = l_params_all$p_REL_OD_NI_lb, 
                                  p_ABS_OD_NI_lb = l_params_all$p_ABS_OD_NI_lb, 
                                  p_BUP_OD_INJ_lb = l_params_all$p_BUP_OD_INJ_lb, 
                                  p_MET_OD_INJ_lb = l_params_all$p_MET_OD_INJ_lb, 
                                  p_REL_OD_INJ_lb = l_params_all$p_REL_OD_INJ_lb, 
                                  p_ABS_OD_INJ_lb = l_params_all$p_ABS_OD_INJ_lb), # lower bound estimate for each param
                         v_ub = c(p_BUP_OD_NI_ub = l_params_all$p_BUP_OD_NI_ub, 
                                  p_MET_OD_NI_ub = l_params_all$p_MET_OD_NI_ub, 
                                  p_REL_OD_NI_ub = l_params_all$p_REL_OD_NI_ub, 
                                  p_ABS_OD_NI_ub = l_params_all$p_ABS_OD_NI_ub, 
                                  p_BUP_OD_INJ_ub = l_params_all$p_BUP_OD_INJ_ub, 
                                  p_MET_OD_INJ_ub = l_params_all$p_MET_OD_INJ_ub, 
                                  p_REL_OD_INJ_ub = l_params_all$p_REL_OD_INJ_ub, 
                                  p_ABS_OD_INJ_ub = l_params_all$p_ABS_OD_INJ_ub)){ # higher bound estimate for each param
                         
                         
                         #v_ub = c(p_S1S2 = 0.50, hr_S1 = 4.5, hr_S2 = 15)){
  n_param <- length(v_param_names)
  m_lhs_unit   <- lhs::randomLHS(n = n_samp, k = n_param) # random latin hypercube sampling
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  for (i in 1:n_param){ # draw parameters
    m_param_samp[, i] <- qunif(m_lhs_unit[,i],
                               min = v_lb[i],
                               max = v_ub[i])
    # ALTERNATIVE prior using beta (or other) distributions
    # m_param_samp[, i] <- qbeta(m_lhs_unit[,i],
    #                            min = 1,
    #                            max = 1)
  }
  return(m_param_samp)
}

#' Evaluate log-prior of calibrated parameters
log_prior <- function(v_params, 
                      v_param_names = c("p_BUP_OD_NI", 
                                        "p_MET_OD_NI", 
                                        "p_REL_OD_NI", 
                                        "p_ABS_OD_NI", 
                                        "p_BUP_OD_INJ", 
                                        "p_MET_OD_INJ", 
                                        "p_REL_OD_INJ", 
                                        "p_ABS_OD_INJ"), # Calibrated params: baseline overdose probability from each health state
                      v_lb = c(p_BUP_OD_NI_lb = l_params_all$p_BUP_OD_NI_lb, 
                               p_MET_OD_NI_lb = l_params_all$p_MET_OD_NI_lb, 
                               p_REL_OD_NI_lb = l_params_all$p_REL_OD_NI_lb, 
                               p_ABS_OD_NI_lb = l_params_all$p_ABS_OD_NI_lb, 
                               p_BUP_OD_INJ_lb = l_params_all$p_BUP_OD_INJ_lb, 
                               p_MET_OD_INJ_lb = l_params_all$p_MET_OD_INJ_lb, 
                               p_REL_OD_INJ_lb = l_params_all$p_REL_OD_INJ_lb, 
                               p_ABS_OD_INJ_lb = l_params_all$p_ABS_OD_INJ_lb), # lower bound estimate for each param
                      v_ub = c(p_BUP_OD_NI_ub = l_params_all$p_BUP_OD_NI_ub, 
                               p_MET_OD_NI_ub = l_params_all$p_MET_OD_NI_ub, 
                               p_REL_OD_NI_ub = l_params_all$p_REL_OD_NI_ub, 
                               p_ABS_OD_NI_ub = l_params_all$p_ABS_OD_NI_ub, 
                               p_BUP_OD_INJ_ub = l_params_all$p_BUP_OD_INJ_ub, 
                               p_MET_OD_INJ_ub = l_params_all$p_MET_OD_INJ_ub, 
                               p_REL_OD_INJ_ub = l_params_all$p_REL_OD_INJ_ub, 
                               p_ABS_OD_INJ_ub = l_params_all$p_ABS_OD_INJ_ub)){ # higher bound estimate for each param
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_param <- length(v_param_names)
  n_samp <- nrow(v_params)
  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)
  for (i in 1:n_param){
    lprior <- lprior + dunif(v_params[, i],
                             min = v_lb[i],
                             max = v_ub[i], 
                             log = T)
    # ALTERNATIVE prior using beta distributions
    # lprior <- lprior + dbeta(v_params[, i],
    #                          min = 1,
    #                          max = 1, 
    #                          log = T)
  }
  return(lprior)
}

#' Evaluate prior of calibrated parameters
prior <- function(v_params) { 
  v_prior <- exp(log_prior(v_params)) 
  return(v_prior)
}

#' Log-likelihood function for a parameter set
log_lik <- function(v_params,
                    l_params_all = load_all_params()){ # User defined
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  
  n_samp <- nrow(v_params)
  v_target_names <- c("Overdoses", "Fatal Overdoses")
  n_target       <- length(v_target_names)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_target) 
  colnames(v_llik) <- v_target_names
  v_llik_overall <- numeric(n_samp)
  for(j in 1:n_samp) { # j=1
    jj <- tryCatch( { 
      
      ###   Run model for parameter set "v_params" ###
      l_model_res <- calibration_out(v_params_calib = v_params[j, ], 
                                     l_params_all = l_params_all)
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      ## TARGET 1: Fatal overdoses ("fatal_overdose")
      ## Normal log-likelihood  
      v_llik[j, "fatal_overdose"] <- sum(dnorm(x = l_cali_targets$ODF$pe,
                                               mean = l_model_res$fatal_overdose,
                                               sd = l_cali_targets$ODF$se,
                                               log = T))
      
      ## TARGET 2: Non-fatal overdoses ("overdose")
      ## Normal log-likelihood
      v_llik[j, "overdose"] <- sum(dnorm(x = l_cali_targets$ODN$pe,
                                         mean = l_model_res$overdose,
                                         sd = l_cali_targets$ODN$se,
                                         log = T))
      

      ## OVERALL
      ## can give different targets different weights (user must change this)
      # To-do: Confirm this calculation
      v_weights <- rep(1, n_target) # weight fatal overdoses 1:1 to overall overdoses
      
      ## weighted sum
      v_llik_overall[j] <- v_llik[j, ] %*% v_weights
    }, error = function(e) NA) 
    if(is.na(jj)) { v_llik_overall <- -Inf }
  } ## End loop over sampled parameter sets
  
  ## return GOF
  return(v_llik_overall)
}

#' Likelihood
likelihood <- function(v_params){ 
  v_like <- exp(log_lik(v_params)) 
  return(v_like)
}

#' Evaluate log-posterior of calibrated parameters
log_post <- function(v_params) { 
  v_lpost <- log_prior(v_params) + log_lik(v_params)
  return(v_lpost) 
}

#' Evaluate posterior of calibrated parameters
posterior <- function(v_params) { 
  v_posterior <- exp(log_post(v_params)) 
  return(v_posterior)
}
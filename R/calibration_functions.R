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

  # Substitute values of calibrated parameters in base-case with calibrated values 
  l_params_all <- update_param_list(l_params_all = l_params_all, params_updated = v_params_calib)
  
  # Run model with updated calibrated parameters
  l_out_markov <- markov_model(l_params_all = l_params_all)
  
  #### Epidemiological Output ####
  ### Overdose deaths ###
  v_ODF <- l_out_markov$m_M_agg_trace[, "ODF"] # cumulative deaths at time i

  ### Non-fatal overdoses ###
  v_ODN <- l_out_markov$m_M_agg_trace[, "ODN"] # cumulative non-fatal overdoses at time i
  
  ### Select time-points ###
  ### Overdose deaths ###
  n_ODF_t1 <- l_cali_targets$ODF$Time[1]
  n_ODF_t2 <- l_cali_targets$ODF$Time[2]
  n_ODF_t3 <- l_cali_targets$ODF$Time[3]
  
  ### Non-fatal overdose ###
  n_ODN_t1 <- l_cali_targets$ODN$Time[1]
  n_ODN_t2 <- l_cali_targets$ODN$Time[2]
  n_ODN_t3 <- l_cali_targets$ODN$Time[3]
  
  ### Subset output by time-points ###
  ### Overdose deaths ###
  # Fatal overdoses already cumulative at each time point
  n_ODF1 <- v_ODF[n_ODF_t1]
  n_ODF2 <- v_ODF[n_ODF_t2]
  n_ODF3 <- v_ODF[n_ODF_t3]
  
  ### Non-fatal overdose
  # Non-fatal overdoses need to be summed across time points to generate cumulative estimates
  n_ODN1 <- sum(v_ODN[c(1:n_ODN_t1)])
  n_ODN2 <- sum(v_ODN[c(1:n_ODN_t2)])
  n_ODN3 <- sum(v_ODN[c(1:n_ODN_t3)])
  
  #### Return Output ####
  l_out <- list(fatal_overdose = c(n_ODF1, n_ODF2, n_ODF3), # cumulative deaths at t1, t2 , t3 time periods (for yearly deaths: (i + 12) - i where i = first month of year, 1 + 12 = last month)
                overdose = c(n_ODN1, n_ODN2, n_ODN3))
  return(l_out)
}

## Sample prior distribution ##
sample.prior <- function(n_samp, 
                         v_param_names = v_cali_param_names, 
                         v_alpha = v_shape, 
                         v_beta = v_scale){
  n_param <- length(v_param_names)
  # random latin hypercube sampling
  m_lhs_unit   <- lhs::randomLHS(n = n_samp, k = n_param) 
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  
  # draw parameters
  draws  <- data.frame(n_TX_OD  = qgamma(m_lhs_unit[,1], shape = v_alpha[1], scale = v_beta[1]),  # n_TX_OD
                       n_TXC_OD  = qgamma(m_lhs_unit[,2], shape = v_alpha[2], scale = v_beta[2]),  # n_TXC_OD
                       n_REL_OD  = qgamma(m_lhs_unit[,3], shape = v_alpha[3], scale = v_beta[3]), # n_REL_OD
                       #n_ABS_OD  = qgamma(m_lhs_unit[,4], shape = n_ABS_OD_shape, scale = n_ABS_OD_scale),   # n_ABS_OD
                       n_TX_OD_mult  = qgamma(m_lhs_unit[,4], shape = v_alpha[4], scale = v_beta[4]),   # n_TX_OD_mult
                       n_TXC_OD_mult  = qgamma(m_lhs_unit[,5], shape = v_alpha[5], scale = v_beta[5]),   # n_TXC_OD_mult
                       n_REL_OD_mult  = qgamma(m_lhs_unit[,6], shape = v_alpha[6], scale = v_beta[6]),                      # n_REL_OD_mult
                       n_INJ_OD_mult  = qgamma(m_lhs_unit[,7], shape = v_alpha[7], scale = v_beta[7]),  # n_INJ_OD_mult
                       n_fatal_OD  = qgamma(m_lhs_unit[,8], shape = v_alpha[8], scale = v_beta[8]))  # n_fatal_OD
  
  # draw parameters (uniform distribution)
  #for (i in 1:n_param){ 
  #  m_param_samp[, i] <- qunif(m_lhs_unit[,i],
  #                             min = v_lb[i],
  #                             max = v_ub[i])
    # ALTERNATIVE prior using beta (or other) distributions
    # m_param_samp[, i] <- qbeta(m_lhs_unit[,i],
    #                            min = 1,
    #                            max = 1)
  #}
  #return(m_param_samp)
  return(as.matrix(draws))
}

#### Log prior ####
log_prior <- function(v_params, 
                      v_param_names = v_cali_param_names, 
                      v_alpha = v_shape, 
                      v_beta = v_scale){
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_param <- length(v_param_names)
  n_samp <- nrow(v_params)
  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)
  
  lprior <- lprior + dgamma(v_params[, 1], shape = v_alpha[1], scale = v_beta[1], log = TRUE) # n_TX_OD
  lprior <- lprior + dgamma(v_params[, 2], shape = v_alpha[2], scale = v_beta[2], log = TRUE) # n_TXC_OD
  lprior <- lprior + dgamma(v_params[, 3], shape = v_alpha[3], scale = v_beta[3], log = TRUE) # n_REL_OD
  lprior <- lprior + dgamma(v_params[, 4], shape = v_alpha[4], scale = v_beta[4], log = TRUE) # n_TX_OD_mult
  lprior <- lprior + dgamma(v_params[, 5], shape = v_alpha[5], scale = v_beta[5], log = TRUE) # n_TXC_OD_mult
  lprior <- lprior + dgamma(v_params[, 6], shape = v_alpha[6], scale = v_beta[6], log = TRUE) # n_REL_OD_mult
  lprior <- lprior + dgamma(v_params[, 7], shape = v_alpha[7], scale = v_beta[7], log = TRUE) # n_INJ_OD_mult
  lprior <- lprior + dgamma(v_params[, 8], shape = v_alpha[8], scale = v_beta[8], log = TRUE) # n_fatal_OD
  
  
  #for (i in 1:n_param){
  #  lprior <- lprior + dunif(v_params[, i],
  #                           min = v_lb[i],
  #                           max = v_ub[i], 
  #                           log = T)
    # ALTERNATIVE prior using beta distributions
    # lprior <- lprior + dbeta(v_params[, i],
    #                          min = 1,
    #                          max = 1, 
    #                          log = T)
  #}
  return(lprior)
}

#' Evaluate prior of calibrated parameters
prior <- function(v_params) { 
  v_prior <- exp(log_prior(v_params)) 
  return(v_prior)
}

#' Log-likelihood function for a parameter set
log_lik <- function(v_params){ # User defined
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  
  n_samp <- nrow(v_params)
  v_target_names <- c("Fatal Overdoses", "Overdoses")
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
      v_llik[j, "Fatal Overdoses"] <- sum(dnorm(x = l_cali_targets$ODF$pe,
                                                mean = l_model_res$fatal_overdose,
                                                sd = l_cali_targets$ODF$se,
                                                log = T))
      ## TARGET 2: Non-fatal overdoses ("overdose")
      ## Normal log-likelihood
      v_llik[j, "Overdoses"] <- sum(dnorm(x = l_cali_targets$ODN$pe,
                                          mean = l_model_res$overdose,
                                          sd = l_cali_targets$ODN$se,
                                          log = T))
      ## can give different targets different weights
      v_weights <- rep(1, n_target) # currently weight fatal overdoses 1:1 to overall overdoses
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
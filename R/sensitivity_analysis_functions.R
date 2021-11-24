#' Generate model outputs for deterministic sensitivity analysis
#'
#' \code{dsa_out} generates model outputs from DSA parameter ranges
#' 
#' @param v_params_dsa is a vector (or single) parameter to change in deterministic analysis.
#' @param v_params_calib is a vector of parameters that need to be calibrated.   
#' @param l_params_all is a list with all parameters of the decision model.
#' @return 
#' A list with all cause deaths, and non-fatal overdoses.
#' @export
dsa_out <- function(v_params_dsa,
                    v_params_calib,
                    l_params_all){
  
  # Substitute values of calibrated parameters in base-case with calibrated values
  # Combine with DSA ranges
  v_params_updated <- c(v_params_dsa, v_params_calib)
  
  l_params_all <- update_param_list(l_params_all = l_params_all, params_updated = v_params_updated)
  
  # Run model with updated calibrated parameters
  l_out_markov <- markov_model(l_params_all = l_params_all)
  
  #### Return Output ####
  l_out <- list(fatal_overdose = c(n_ODF1, n_ODF2, n_ODF3), # cumulative deaths at t1, t2 , t3 time periods (for yearly deaths: (i + 12) - i where i = first month of year, 1 + 12 = last month)
                overdose = c(n_ODN1, n_ODN2, n_ODN3))
  return(l_out)
}
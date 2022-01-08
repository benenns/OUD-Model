#' Markov model
#'
#' \code{markov_model} implements the main model functions to calculate Markov trace.
#'
#' @param l_params_all List with all parameters
#' @param err_stop Logical variable to stop model run if transition array is invalid, if TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE
#' @param checks Logical variable to indicate output of visual checks (e.g. slices of transition array)
#' @param cali Logical variable to adjust model cutoff to only calibration period
#' @return 
#' a_TDP: Transition probability array
#' m_M_trace: Fully stratified markov cohort trace
#' m_M_agg_trace: Aggregated markov trace over base health states
#' m_M_agg_trace_death: State-specific mortality from each health state
#' m_M_agg_trace_sero: HIV seroconversions from each health state
#' @export
markov_model <- function(l_params_all, err_stop = FALSE, verbose = FALSE, checks = FALSE, cali = FALSE){
  ### Definition:
  ##   Markov model implementation function
  ### Prefixes:
  ##   l_* denotes list
  ##   n_* denotes number
  ##   a_* denotes 3-D array
  ##   m_* denotes 2-D matrix
  ##   v_* denotes vector
  ##   df_* denotes data frame
  ##   p_* denotes transition parameters
  ##   c_* denotes costs
  ##   u_* denotes utilities
  ### Arguments:  
  ##   l_params_all: List with all parameters
  ##   verbose: Logical variable to indicate print out of messages
  ### Returns:
  ##   a_TDP: Transition probability array.
  ##   m_M_trace: Fully disaggregated matrix cohort trace.
  ##   m_M_agg_trace: Aggregated trace over selected base health states.
  ##   m_M_agg_trace_death: State-specific mortality from each health state.
  ##   m_M_agg_trace_sero: HIV seroconversions from each health state.
  ##
  with(as.list(l_params_all), {

  #### Set up model states ####
  l_dim_s  <- list() # list of health states
  
  # Base health states
  BASE <- l_dim_s[[1]] <- c("MET", "METC", "BUP", "BUPC", "ABS", "REL", "ODN", "ODF")
  # Injection/non-injection stratification
  INJECT <- l_dim_s[[2]] <- c("NI", "INJ")
  # Episodes (1-3)
  EP <-  l_dim_s[[3]] <- c("1", "2", "3")
  # HIV/HCV status
  SERO <- l_dim_s[[4]] <- c("NEG", "HIV", "HCV", "COI")
  # Set model periods
  if(cali == TRUE){
    # Calibration periods
    n_t <- (n_cali_max_per + 1) # if calibrating, cut model off at max calibration output (e.g. 36 months for three-years)
  } else{
    # Maximum model periods(regular)
    n_t <- (n_age_max - n_age_init) * n_per # convert years
  }
  
  df_flat <- expand.grid(l_dim_s) #combine all elements together into vector of health states
  df_flat <- dplyr::rename(df_flat, BASE    = Var1, 
                                    INJECT  = Var2, 
                                    EP      = Var3, 
                                    SERO    = Var4)

  # Create index of states to populate transition matrices
  # All treatment
  TX <- df_flat$BASE == "BUP" | df_flat$BASE == "MET"
  TXC <- df_flat$BASE == "BUPC" | df_flat$BASE == "METC"
  all_TX <- df_flat$BASE == "BUP" | df_flat$BASE == "BUPC" | df_flat$BASE == "MET" | df_flat$BASE == "METC"
  
  # All out-of-treatment (incl ABS)
  OOT <- df_flat$BASE == "REL" | df_flat$BASE == "ABS" | df_flat$BASE == "ODN" | df_flat$BASE == "ODF" 
  
  # Buprenorphine
  BUP  <- df_flat$BASE == "BUP" # treatment only 
  BUPC <- df_flat$BASE == "BUPC" # concurrent opioid use
  all_BUP <- df_flat$BASE == "BUP" | df_flat$BASE == "BUPC"
  
  # Methadone
  MET  <- df_flat$BASE == "MET" # treatment only
  METC <- df_flat$BASE == "METC" # concurrent opioid use
  all_MET <- df_flat$BASE == "MET" | df_flat$BASE == "METC"
  
  # Relapse
  REL <- df_flat$BASE == "REL"
  
  # Overdose
  all_OD <- df_flat$BASE == "ODN" | df_flat$BASE == "ODF"
  ODN <- df_flat$BASE == "ODN" # non-fatal overdose
  ODF <- df_flat$BASE == "ODF" # fatal overdose
  
  # Abstinence
  ABS <- df_flat$BASE == "ABS"
  
  # Serostatus
  NEG <- df_flat$SERO == "NEG"
  HIV <- df_flat$SERO == "HIV"
  HCV <- df_flat$SERO == "HCV"
  COI <- df_flat$SERO == "COI"
  
  all_HIV <- df_flat$SERO == "HIV" | df_flat$SERO == "COI"
  all_HCV <- df_flat$SERO == "HCV" | df_flat$SERO == "COI"
  
  # Injection
  INJ <- df_flat$INJECT == "INJ"
  NI <- df_flat$INJECT == "NI"
  
  # Episodes
  EP1 <- df_flat$EP == "1"
  EP2 <- df_flat$EP == "2"
  EP3 <- df_flat$EP == "3"

  df_n <- unite(df_flat, newCol) # combine columns into one data frame of all health states
  v_n_states <- df_n[,1] # convert df into vector
  n_states <- length(v_n_states) # total number of health states
  l_index_s  <- list(TX = TX, OOT = OOT, 
                     BUP = BUP, BUPC = BUPC,
                     MET = MET, METC = METC,
                     REL = REL, 
                     all_OD = all_OD, ODN = ODN, ODF = ODF, 
                     ABS = ABS, 
                     NEG = NEG, HIV = HIV, HCV = HCV, COI = COI, 
                     INJ = INJ, NI = NI, 
                     EP1 = EP1, EP2 = EP2, EP3 = EP3)
  
  #### Overdose probability ####
  #' Probability of non-fatal and fatal overdose
  #'
  #' \code{p_OD} is used to calculate overdose probabilities from health states. This function also requires additional overdose/fentanyl/naloxone parameters included in `l_params_all`
  #'
  #' @param rate Baseline overdose rate for health states
  #' @param rate_fatal Fatal overdose rate
  #' @param rate_fent Fentanyl overdose rate
  #' @param multiplier Multiplier for elevated overdose in first month of health state
  #' @param fent_mult Multiplier for overdose rate in health states when exposed to fentanyl
  #' @param first_month Logical parameter to switch between month 1 and month 2+ for parameter estimation
  #' @param fatal Logical parameter to switch between fatal/non-fatal overdose
  #' @param injection Logical parameter to adjust rate calculation for injection/non-injection use
  #' @param time Time period for time-varying parameters
  #' 
  #' @return 
  #' `p_OD` monthly probability of fatal or non-fatal overdose from a given health state
  #' @export
  p_OD <- function(rate,# = rate,
                   rate_fatal,# = rate_fatal,
                   #rate_fent = n_fent_OD,
                   multiplier,# = multiplier,
                   fent_mult,# = fent_mult,
                   #fent_reduction_state = fent_reduction_state,
                   time,# = time,
                   first_month = FALSE,
                   fatal = FALSE,
                   injection = FALSE){
    
    # Probability of successful naloxone use
    p_NX_rev <- (p_witness * p_NX_used * p_NX_success)
    
    # Probability of mortality from overdose accounting for baseline overdose fatality and effectiveness of naloxone
    # Subsets overdose into fatal and non-fatal, conditional on different parameters
    
    # Convert fatal overdose rate into probability of death following overdose
    p_fatal_OD <- 1 - exp(-(rate_fatal))
    
    # Convert fentanyl overdose rate into probability
    #p_fent_OD <- 1 - exp(-(rate_fent))
    
    # Probability of fentanyl exposure 
    # Generate time-varying probability of fentanyl exposure
    # Currently modeling logarithmic growth based on 2018-2020 period (consider cutting off at some point, e.g., 5-years, 10-years)
    #v_fent_exp_rate <- rep(0, n_t)
    
    #for(i in 2:n_t){
    #  v_fent_exp_rate[1] <- -log(1- p_fent_exp_base)
    #  v_fent_exp_rate[i] <- v_fent_exp_rate[i-1] + n_fent_growth_rate
    #}
    #v_fent_exp_prob <- 1 - exp(-v_fent_exp_rate) # create vector of monthly fentanyl exposure probabilities (generated by growth rates from 2018-2020)
    
    v_fent_exp_prob <- c(p_fent_exp_2018, p_fent_exp_2019, p_fent_exp_2020)
    
    # Adjustment for injection/non-injection
    #if (injection){
    #  v_fent_exp_prob <- v_fent_exp_prob
    #}  else{
    #  v_fent_exp_prob <- v_fent_exp_prob * p_ni_fent_reduction
    #}
    
    # Convert input monthly rates to monthly probabilities - multiply rates by first month multiplier before converting
    if (injection == TRUE && first_month == TRUE){
      p_base_OD <- 1 - exp(-(rate * n_INJ_OD_mult * multiplier * (v_fent_exp_prob[time] * fent_mult)))
      #p_fent_OD <- 1 - exp(-(rate_fent * multiplier))
    }
    else if (injection == TRUE && first_month == FALSE){
      p_base_OD <- 1 - exp(-(rate * n_INJ_OD_mult * (v_fent_exp_prob[time] * fent_mult)))
      #p_fent_OD <- 1 - exp(-(rate_fent))
    }
    else if (injection == FALSE && first_month == TRUE){
      p_base_OD <- 1 - exp(-(rate * multiplier * (v_fent_exp_prob[time] * p_ni_fent_reduction * fent_mult)))
      #p_fent_OD <- 1 - exp(-(rate_fent * multiplier))
    }
    else if (injection == FALSE && first_month == FALSE){
      p_base_OD <- 1 - exp(-(rate * (v_fent_exp_prob[time] * p_ni_fent_reduction * fent_mult)))
      #p_fent_OD <- 1 - exp(-(rate_fent))
    }

    # Naloxone effect on fatal overdose
    p_fatal_OD_NX <- p_fatal_OD * (1 - p_NX_rev)
    
    # Probability of fentanyl exposure (adjusted for injection/non-injection)
    #if (injection){
    #p_fent_exp <- p_fent_exp
    #}  else{
    #  p_fent_exp <- p_fent_exp * p_ni_fent_reduction
    #}
    
    # Calculate fatal and non-fatal overdose probabilities
    if (fatal == TRUE){
      p_OD <- p_base_OD * p_fatal_OD_NX #((p_base_OD * (1 - p_fent_exp)) + (p_fent_OD * (p_fent_exp))) * p_fatal_OD_NX
    } else{
      p_OD <- p_base_OD * (1 - p_fatal_OD_NX) #((p_base_OD * (1 - p_fent_exp)) + (p_fent_OD * (p_fent_exp))) * (1 - p_fatal_OD_NX)
    }
    return(p_OD)
  }

  # Module to calculate probability of overdose from states
  # Four separate matrices to account for state-time (first month vs. second+), and model-time (changing fentanyl prevalence, etc.)
  #### Time-dependent overdose probabilities ####
  # Time periods
  time_periods <- n_cali_per
  
  # Empty 2-D matrix
  m_ODN <- m_ODN_first <- m_ODF <- m_ODF_first <- array(0, dim = c(n_states, time_periods),
                                                  dimnames = list(v_n_states, 1:time_periods))
  
  #test<- p_OD(rate = n_TX_OD, rate_fatal = n_fatal_OD, multiplier = n_TX_OD_mult, fent_mult = n_fent_OD_mult, time = 1, first_month = TRUE, fatal = FALSE, injection = FALSE)
  
  for(i in 1:time_periods){
  # Probability of overdose
  # Non-fatal (first month)
  m_ODN_first[TX & NI, i]   <- p_OD(rate = n_TX_OD,    rate_fatal = n_fatal_OD, multiplier = n_TX_OD_mult,  fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = FALSE, injection = FALSE)
  m_ODN_first[TXC & NI, i]  <- p_OD(rate = n_TXC_OD,  rate_fatal = n_fatal_OD, multiplier = n_TXC_OD_mult,  fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = FALSE, injection = FALSE)
  m_ODN_first[REL & NI, i]  <- p_OD(rate = n_REL_OD,  rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult,  fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = FALSE, injection = FALSE)
  m_ODN_first[ABS & NI, i]  <- p_OD(rate = n_ABS_OD,  rate_fatal = n_fatal_OD, multiplier = n_ABS_OD_mult,  fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = FALSE, injection = FALSE)
  m_ODN_first[TX & INJ, i]  <- p_OD(rate = n_TX_OD,   rate_fatal = n_fatal_OD, multiplier = n_TX_OD_mult,  fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = FALSE, injection = TRUE)
  m_ODN_first[TXC & INJ, i] <- p_OD(rate = n_TXC_OD, rate_fatal = n_fatal_OD, multiplier = n_TXC_OD_mult,  fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = FALSE, injection = TRUE)
  m_ODN_first[REL & INJ, i] <- p_OD(rate = n_REL_OD, rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult,  fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = FALSE, injection = TRUE)
  m_ODN_first[ABS & INJ, i] <- p_OD(rate = n_ABS_OD, rate_fatal = n_fatal_OD, multiplier = n_ABS_OD_mult,  fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = FALSE, injection = TRUE)
    
  # Fatal (first month)
  m_ODF_first[TX & NI, i] <- p_OD(rate = n_TX_OD,    rate_fatal = n_fatal_OD, multiplier = n_TX_OD_mult,   fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = TRUE, injection = FALSE)
  m_ODF_first[TXC & NI, i] <- p_OD(rate = n_TXC_OD,  rate_fatal = n_fatal_OD, multiplier = n_TXC_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = TRUE, injection = FALSE)
  m_ODF_first[REL & NI, i] <- p_OD(rate = n_REL_OD,  rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = TRUE, injection = FALSE)
  m_ODF_first[ODN & NI, i] <- p_OD(rate = n_REL_OD,  rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = TRUE, injection = FALSE)
  m_ODF_first[ABS & NI, i] <- p_OD(rate = n_ABS_OD,  rate_fatal = n_fatal_OD, multiplier = n_ABS_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = TRUE, injection = FALSE)
  m_ODF_first[TX & INJ, i] <- p_OD(rate = n_TX_OD,   rate_fatal = n_fatal_OD, multiplier = n_TX_OD_mult,   fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = TRUE, injection = TRUE)
  m_ODF_first[TXC & INJ, i] <- p_OD(rate = n_TXC_OD, rate_fatal = n_fatal_OD, multiplier = n_TXC_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = TRUE, injection = TRUE)
  m_ODF_first[REL & INJ, i] <- p_OD(rate = n_REL_OD, rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = TRUE, injection = TRUE)
  m_ODF_first[ODN & INJ, i] <- p_OD(rate = n_REL_OD, rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = TRUE, injection = TRUE)
  m_ODF_first[ABS & INJ, i] <- p_OD(rate = n_ABS_OD, rate_fatal = n_fatal_OD, multiplier = n_ABS_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = TRUE, fatal = TRUE, injection = TRUE)
    
  # Non-fatal (month 2+)
  m_ODN[TX & NI, i] <- p_OD(rate = n_TX_OD,    rate_fatal = n_fatal_OD, multiplier = n_TX_OD_mult,   fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = FALSE, injection = FALSE)
  m_ODN[TXC & NI, i] <- p_OD(rate = n_TXC_OD,  rate_fatal = n_fatal_OD, multiplier = n_TXC_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = FALSE, injection = FALSE)
  m_ODN[REL & NI, i] <- p_OD(rate = n_REL_OD,  rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = FALSE, injection = FALSE)
  m_ODN[ABS & NI, i] <- p_OD(rate = n_ABS_OD,  rate_fatal = n_fatal_OD, multiplier = n_ABS_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = FALSE, injection = FALSE)
  m_ODN[TX & INJ, i] <- p_OD(rate = n_TX_OD,   rate_fatal = n_fatal_OD, multiplier = n_TX_OD_mult,   fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = FALSE, injection = TRUE)
  m_ODN[TXC & INJ, i] <- p_OD(rate = n_TXC_OD, rate_fatal = n_fatal_OD, multiplier = n_TXC_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = FALSE, injection = TRUE)
  m_ODN[REL & INJ, i] <- p_OD(rate = n_REL_OD, rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = FALSE, injection = TRUE)
  m_ODN[ABS & INJ, i] <- p_OD(rate = n_ABS_OD, rate_fatal = n_fatal_OD, multiplier = n_ABS_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = FALSE, injection = TRUE)
  
  # Fatal (month 2+)
  m_ODF[TX & NI, i] <- p_OD(rate = n_TX_OD,    rate_fatal = n_fatal_OD, multiplier = n_TX_OD_mult,   fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = TRUE, injection = FALSE)
  m_ODF[TXC & NI, i] <- p_OD(rate = n_TXC_OD,  rate_fatal = n_fatal_OD, multiplier = n_TXC_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = TRUE, injection = FALSE)
  m_ODF[REL & NI, i] <- p_OD(rate = n_REL_OD,  rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = TRUE, injection = FALSE)
  m_ODF[ODN & NI, i] <- p_OD(rate = n_REL_OD,  rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = TRUE, injection = FALSE)
  m_ODF[ABS & NI, i] <- p_OD(rate = n_ABS_OD,  rate_fatal = n_fatal_OD, multiplier = n_ABS_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = TRUE, injection = FALSE)
  m_ODF[TX & INJ, i] <- p_OD(rate = n_TX_OD,   rate_fatal = n_fatal_OD, multiplier = n_TX_OD_mult,   fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = TRUE, injection = TRUE)
  m_ODF[TXC & INJ, i] <- p_OD(rate = n_TXC_OD, rate_fatal = n_fatal_OD, multiplier = n_TXC_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = TRUE, injection = TRUE)
  m_ODF[REL & INJ, i] <- p_OD(rate = n_REL_OD, rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = TRUE, injection = TRUE)
  m_ODF[ODN & INJ, i] <- p_OD(rate = n_REL_OD, rate_fatal = n_fatal_OD, multiplier = n_REL_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = TRUE, injection = TRUE)
  m_ODF[ABS & INJ, i] <- p_OD(rate = n_ABS_OD, rate_fatal = n_fatal_OD, multiplier = n_ABS_OD_mult, fent_mult = n_fent_OD_mult, time = i, first_month = FALSE, fatal = TRUE, injection = TRUE)
  }
  
  if (checks){
    # Overdose
    write.csv(m_ODN_first,"checks/overdose/m_ODN_first.csv", row.names = TRUE)
    write.csv(m_ODF_first,"checks/overdose/m_ODF_first.csv", row.names = TRUE)
    write.csv(m_ODN,"checks/overdose/m_ODN.csv", row.names = TRUE)
    write.csv(m_ODF,"checks/overdose/m_ODF.csv", row.names = TRUE)
  } else{}
  
  #### Time-dependent remain probabilities ####
  # Empty 2-D matrix
  # Three matrices to account for overdose probabilities in 2018, 2019, and 2020+
  m_TDP_1 <- m_TDP_2 <- m_TDP_3 <- array(0, dim = c(n_states, n_t),
                                     dimnames = list(v_n_states, 1:n_t))
  
  # Generate state-specific weibull parameters
  # Shape
  p_weibull_shape_BUP_NI <- p_weibull_shape_BUP_INJ <- p_weibull_shape_BUP
  p_weibull_shape_MET_NI <- p_weibull_shape_MET_INJ <- p_weibull_shape_MET
  p_weibull_shape_ABS_NI <- p_weibull_shape_ABS_INJ <- p_weibull_shape_ABS
  p_weibull_shape_REL_NI <- p_weibull_shape_REL_INJ <- p_weibull_shape_REL
  
  # Scale
  p_weibull_scale_BUP_NI <- p_weibull_scale_BUP_INJ <- p_weibull_scale_BUP
  p_weibull_scale_MET_NI <- p_weibull_scale_MET_INJ <- p_weibull_scale_MET
  p_weibull_scale_ABS_NI <- p_weibull_scale_ABS_INJ <- p_weibull_scale_ABS
  p_weibull_scale_REL_NI <- p_weibull_scale_REL_INJ <- p_weibull_scale_REL

  # Generate state-specific frailty terms
  # Non-injection
  p_frailty_BUP_NI_1 <- p_frailty_BUP_1
  p_frailty_BUP_NI_2 <- p_frailty_BUP_2
  p_frailty_BUP_NI_3 <- p_frailty_BUP_3
  p_frailty_BUPC_NI_1 <- p_frailty_BUP_1 * p_frailty_BUPC
  p_frailty_BUPC_NI_2 <- p_frailty_BUP_2 * p_frailty_BUPC
  p_frailty_BUPC_NI_3 <- p_frailty_BUP_3 * p_frailty_BUPC
  p_frailty_MET_NI_1 <- p_frailty_MET_1
  p_frailty_MET_NI_2 <- p_frailty_MET_2
  p_frailty_MET_NI_3 <- p_frailty_MET_3
  p_frailty_METC_NI_1 <- p_frailty_MET_1 * p_frailty_METC
  p_frailty_METC_NI_2 <- p_frailty_MET_2 * p_frailty_METC
  p_frailty_METC_NI_3 <- p_frailty_MET_3 * p_frailty_METC
  p_frailty_ABS_NI_1 <- p_frailty_ABS_1
  p_frailty_ABS_NI_2 <- p_frailty_ABS_2
  p_frailty_ABS_NI_3 <- p_frailty_ABS_3
  p_frailty_REL_NI_1 <- p_frailty_REL_1
  p_frailty_REL_NI_2 <- p_frailty_REL_2
  p_frailty_REL_NI_3 <- p_frailty_REL_3
  
  # Injection
  p_frailty_BUP_INJ_1 <- p_frailty_BUP_1 * p_frailty_BUP_INJ
  p_frailty_BUP_INJ_2 <- p_frailty_BUP_2 * p_frailty_BUP_INJ
  p_frailty_BUP_INJ_3 <- p_frailty_BUP_3 * p_frailty_BUP_INJ
  p_frailty_BUPC_INJ_1 <- p_frailty_BUP_1 * p_frailty_BUP_INJ * p_frailty_BUPC
  p_frailty_BUPC_INJ_2 <- p_frailty_BUP_2 * p_frailty_BUP_INJ * p_frailty_BUPC
  p_frailty_BUPC_INJ_3 <- p_frailty_BUP_3 * p_frailty_BUP_INJ * p_frailty_BUPC
  p_frailty_MET_INJ_1 <- p_frailty_MET_1 * p_frailty_MET_INJ
  p_frailty_MET_INJ_2 <- p_frailty_MET_2 * p_frailty_MET_INJ
  p_frailty_MET_INJ_3 <- p_frailty_MET_3 * p_frailty_MET_INJ
  p_frailty_METC_INJ_1 <- p_frailty_MET_1 * p_frailty_MET_INJ * p_frailty_METC
  p_frailty_METC_INJ_2 <- p_frailty_MET_2 * p_frailty_MET_INJ * p_frailty_METC
  p_frailty_METC_INJ_3 <- p_frailty_MET_3 * p_frailty_MET_INJ * p_frailty_METC
  p_frailty_ABS_INJ_1 <- p_frailty_ABS_1 * p_frailty_ABS_INJ
  p_frailty_ABS_INJ_2 <- p_frailty_ABS_2 * p_frailty_ABS_INJ
  p_frailty_ABS_INJ_3 <- p_frailty_ABS_3 * p_frailty_ABS_INJ
  p_frailty_REL_INJ_1 <- p_frailty_REL_1 * p_frailty_REL_INJ
  p_frailty_REL_INJ_2 <- p_frailty_REL_2 * p_frailty_REL_INJ
  p_frailty_REL_INJ_3 <- p_frailty_REL_3 * p_frailty_REL_INJ

  # Probability of remaining in health state
  # All remain in fatal overdose, remain probability = 1
  for(i in 1:n_t){
    # Non-injection
      # Episode 1
      m_TDP_1[EP1 & BUP & NI, i]  <- m_TDP_2[EP1 & BUP & NI, i]  <- m_TDP_3[EP1 & BUP & NI, i]  <- as.vector(exp(p_frailty_BUP_NI_1 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI)))) # vector of remain probabilities
      m_TDP_1[EP1 & BUPC & NI, i] <- m_TDP_2[EP1 & BUPC & NI, i] <- m_TDP_3[EP1 & BUPC & NI, i] <- as.vector(exp(p_frailty_BUPC_NI_1 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI))))
      m_TDP_1[EP1 & MET & NI, i]  <- m_TDP_2[EP1 & MET & NI, i]  <- m_TDP_3[EP1 & MET & NI, i]  <- as.vector(exp(p_frailty_MET_NI_1 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
      m_TDP_1[EP1 & METC & NI, i] <- m_TDP_2[EP1 & METC & NI, i] <- m_TDP_3[EP1 & METC & NI, i] <- as.vector(exp(p_frailty_METC_NI_1 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
      m_TDP_1[EP1 & ABS & NI, i]  <- m_TDP_2[EP1 & ABS & NI, i]  <- m_TDP_3[EP1 & ABS & NI, i]  <- as.vector(exp(p_frailty_ABS_NI_1 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
      m_TDP_1[EP1 & REL & NI, i]  <- m_TDP_2[EP1 & REL & NI, i]  <- m_TDP_3[EP1 & REL & NI, i]  <- as.vector(exp(p_frailty_REL_NI_1 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
      m_TDP_1[EP1 & ODN & NI, i]  <- m_ODN_first[EP1 & REL & NI, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
      m_TDP_2[EP1 & ODN & NI, i]  <- m_ODN_first[EP1 & REL & NI, 2]
      m_TDP_3[EP1 & ODN & NI, i]  <- m_ODN_first[EP1 & REL & NI, 3]
      m_TDP_1[EP1 & ODF & NI, i]  <- m_TDP_2[EP1 & ODF & NI, i] <- m_TDP_3[EP1 & ODF & NI, i]  <- 1 # all remain in ODF
      # Episode 2
      m_TDP_1[EP2 & BUP & NI, i]  <- m_TDP_2[EP2 & BUP & NI, i]  <- m_TDP_3[EP2 & BUP & NI, i]  <- as.vector(exp(p_frailty_BUP_NI_2 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI)))) # vector of remain probabilities
      m_TDP_1[EP2 & BUPC & NI, i] <- m_TDP_2[EP2 & BUPC & NI, i] <- m_TDP_3[EP2 & BUPC & NI, i] <- as.vector(exp(p_frailty_BUPC_NI_2 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI))))
      m_TDP_1[EP2 & MET & NI, i]  <- m_TDP_2[EP2 & MET & NI, i]  <- m_TDP_3[EP2 & MET & NI, i]  <- as.vector(exp(p_frailty_MET_NI_2 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
      m_TDP_1[EP2 & METC & NI, i] <- m_TDP_2[EP2 & METC & NI, i] <- m_TDP_3[EP2 & METC & NI, i] <- as.vector(exp(p_frailty_METC_NI_2 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
      m_TDP_1[EP2 & ABS & NI, i]  <- m_TDP_2[EP2 & ABS & NI, i]  <- m_TDP_3[EP2 & ABS & NI, i]  <- as.vector(exp(p_frailty_ABS_NI_2 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
      m_TDP_1[EP2 & REL & NI, i]  <- m_TDP_2[EP2 & REL & NI, i]  <- m_TDP_3[EP2 & REL & NI, i]  <- as.vector(exp(p_frailty_REL_NI_2 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
      m_TDP_1[EP2 & ODN & NI, i]  <- m_ODN_first[EP2 & REL & NI, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
      m_TDP_2[EP2 & ODN & NI, i]  <- m_ODN_first[EP2 & REL & NI, 2]
      m_TDP_3[EP2 & ODN & NI, i]  <- m_ODN_first[EP2 & REL & NI, 3]
      m_TDP_1[EP2 & ODF & NI, i]  <- m_TDP_2[EP2 & ODF & NI, i] <- m_TDP_3[EP2 & ODF & NI, i]  <- 1 # all remain in ODF
      # Episode 3
      m_TDP_1[EP3 & BUP & NI, i]  <- m_TDP_2[EP3 & BUP & NI, i]  <- m_TDP_3[EP3 & BUP & NI, i]  <- as.vector(exp(p_frailty_BUP_NI_3 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI)))) # vector of remain probabilities
      m_TDP_1[EP3 & BUPC & NI, i] <- m_TDP_2[EP3 & BUPC & NI, i] <- m_TDP_3[EP3 & BUPC & NI, i] <- as.vector(exp(p_frailty_BUPC_NI_3 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI))))
      m_TDP_1[EP3 & MET & NI, i]  <- m_TDP_2[EP3 & MET & NI, i]  <- m_TDP_3[EP3 & MET & NI, i]  <- as.vector(exp(p_frailty_MET_NI_3 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
      m_TDP_1[EP3 & METC & NI, i] <- m_TDP_2[EP3 & METC & NI, i] <- m_TDP_3[EP3 & METC & NI, i] <- as.vector(exp(p_frailty_METC_NI_3 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
      m_TDP_1[EP3 & ABS & NI, i]  <- m_TDP_2[EP3 & ABS & NI, i]  <- m_TDP_3[EP3 & ABS & NI, i]  <- as.vector(exp(p_frailty_ABS_NI_3 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
      m_TDP_1[EP3 & REL & NI, i]  <- m_TDP_2[EP3 & REL & NI, i]  <- m_TDP_3[EP3 & REL & NI, i]  <- as.vector(exp(p_frailty_REL_NI_3 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
      m_TDP_1[EP3 & ODN & NI, i]  <- m_ODN_first[EP3 & REL & NI, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
      m_TDP_2[EP3 & ODN & NI, i]  <- m_ODN_first[EP3 & REL & NI, 2]
      m_TDP_3[EP3 & ODN & NI, i]  <- m_ODN_first[EP3 & REL & NI, 3]
      m_TDP_1[EP3 & ODF & NI, i]  <- m_TDP_2[EP3 & ODF & NI, i] <- m_TDP_3[EP3 & ODF & NI, i]  <- 1 # all remain in ODF
    # Injection
      # Episode 1
      m_TDP_1[EP1 & BUP & INJ, i]  <- m_TDP_2[EP1 & BUP & INJ, i]  <- m_TDP_3[EP1 & BUP & INJ, i]  <- as.vector(exp(p_frailty_BUP_INJ_1 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ)))) # vector of remain probabilities
      m_TDP_1[EP1 & BUPC & INJ, i] <- m_TDP_2[EP1 & BUPC & INJ, i] <- m_TDP_3[EP1 & BUPC & INJ, i] <- as.vector(exp(p_frailty_BUPC_INJ_1 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ))))
      m_TDP_1[EP1 & MET & INJ, i]  <- m_TDP_2[EP1 & MET & INJ, i]  <- m_TDP_3[EP1 & MET & INJ, i]  <- as.vector(exp(p_frailty_MET_INJ_1 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
      m_TDP_1[EP1 & METC & INJ, i] <- m_TDP_2[EP1 & METC & INJ, i] <- m_TDP_3[EP1 & METC & INJ, i] <- as.vector(exp(p_frailty_METC_INJ_1 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
      m_TDP_1[EP1 & ABS & INJ, i]  <- m_TDP_2[EP1 & ABS & INJ, i]  <- m_TDP_3[EP1 & ABS & INJ, i]  <- as.vector(exp(p_frailty_ABS_INJ_1 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
      m_TDP_1[EP1 & REL & INJ, i]  <- m_TDP_2[EP1 & REL & INJ, i]  <- m_TDP_3[EP1 & REL & INJ, i]  <- as.vector(exp(p_frailty_REL_INJ_1 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
      m_TDP_1[EP1 & ODN & INJ, i]  <- m_ODN_first[EP1 & REL & INJ, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
      m_TDP_2[EP1 & ODN & INJ, i]  <- m_ODN_first[EP1 & REL & INJ, 2]
      m_TDP_3[EP1 & ODN & INJ, i]  <- m_ODN_first[EP1 & REL & INJ, 3]
      m_TDP_1[EP1 & ODF & INJ, i]  <- m_TDP_2[EP1 & ODF & INJ, i] <- m_TDP_3[EP1 & ODF & INJ, i]  <- 1 # all remain in ODF
      # Episode 2
      m_TDP_1[EP2 & BUP & INJ, i]  <- m_TDP_2[EP2 & BUP & INJ, i]  <- m_TDP_3[EP2 & BUP & INJ, i]  <- as.vector(exp(p_frailty_BUP_INJ_2 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ)))) # vector of remain probabilities
      m_TDP_1[EP2 & BUPC & INJ, i] <- m_TDP_2[EP2 & BUPC & INJ, i] <- m_TDP_3[EP2 & BUPC & INJ, i] <- as.vector(exp(p_frailty_BUPC_INJ_2 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ))))
      m_TDP_1[EP2 & MET & INJ, i]  <- m_TDP_2[EP2 & MET & INJ, i]  <- m_TDP_3[EP2 & MET & INJ, i]  <- as.vector(exp(p_frailty_MET_INJ_2 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
      m_TDP_1[EP2 & METC & INJ, i] <- m_TDP_2[EP2 & METC & INJ, i] <- m_TDP_3[EP2 & METC & INJ, i] <- as.vector(exp(p_frailty_METC_INJ_2 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
      m_TDP_1[EP2 & ABS & INJ, i]  <- m_TDP_2[EP2 & ABS & INJ, i]  <- m_TDP_3[EP2 & ABS & INJ, i]  <- as.vector(exp(p_frailty_ABS_INJ_2 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
      m_TDP_1[EP2 & REL & INJ, i]  <- m_TDP_2[EP2 & REL & INJ, i]  <- m_TDP_3[EP2 & REL & INJ, i]  <- as.vector(exp(p_frailty_REL_INJ_2 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
      m_TDP_1[EP2 & ODN & INJ, i]  <- m_ODN_first[EP2 & REL & INJ, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
      m_TDP_2[EP2 & ODN & INJ, i]  <- m_ODN_first[EP2 & REL & INJ, 2]
      m_TDP_3[EP2 & ODN & INJ, i]  <- m_ODN_first[EP2 & REL & INJ, 3]
      m_TDP_1[EP2 & ODF & INJ, i]  <- m_TDP_2[EP2 & ODF & INJ, i] <- m_TDP_3[EP2 & ODF & INJ, i]  <- 1 # all remain in ODF
      # Episode 3
      m_TDP_1[EP3 & BUP & INJ, i]  <- m_TDP_2[EP3 & BUP & INJ, i]  <- m_TDP_3[EP3 & BUP & INJ, i]  <- as.vector(exp(p_frailty_BUP_INJ_3 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ)))) # vector of remain probabilities
      m_TDP_1[EP3 & BUPC & INJ, i] <- m_TDP_2[EP3 & BUPC & INJ, i] <- m_TDP_3[EP3 & BUPC & INJ, i] <- as.vector(exp(p_frailty_BUPC_INJ_3 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ))))
      m_TDP_1[EP3 & MET & INJ, i]  <- m_TDP_2[EP3 & MET & INJ, i]  <- m_TDP_3[EP3 & MET & INJ, i]  <- as.vector(exp(p_frailty_MET_INJ_3 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
      m_TDP_1[EP3 & METC & INJ, i] <- m_TDP_2[EP3 & METC & INJ, i] <- m_TDP_3[EP3 & METC & INJ, i] <- as.vector(exp(p_frailty_METC_INJ_3 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
      m_TDP_1[EP3 & ABS & INJ, i]  <- m_TDP_2[EP3 & ABS & INJ, i]  <- m_TDP_3[EP3 & ABS & INJ, i]  <- as.vector(exp(p_frailty_ABS_INJ_3 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
      m_TDP_1[EP3 & REL & INJ, i]  <- m_TDP_2[EP3 & REL & INJ, i]  <- m_TDP_3[EP3 & REL & INJ, i]  <- as.vector(exp(p_frailty_REL_INJ_3 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
      m_TDP_1[EP3 & ODN & INJ, i]  <- m_ODN_first[EP3 & REL & INJ, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
      m_TDP_2[EP3 & ODN & INJ, i]  <- m_ODN_first[EP3 & REL & INJ, 2]
      m_TDP_3[EP3 & ODN & INJ, i]  <- m_ODN_first[EP3 & REL & INJ, 3]
      m_TDP_1[EP3 & ODF & INJ, i]  <- m_TDP_2[EP3 & ODF & INJ, i] <- m_TDP_3[EP3 & ODF & INJ, i]  <- 1 # all remain in ODF
    }
    # Month 2+ (state-time)
    #for(i in 2:n_t){
      # Non-injection
        # Episode 1
        #m_TDP_1[EP1 & BUP & NI, i]  <- m_TDP_2[EP1 & BUP & NI, i]  <- m_TDP_3[EP1 & BUP & NI, i]  <- as.vector(exp(p_frailty_BUP_NI_1 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI)))) # vector of remain probabilities
        #m_TDP_1[EP1 & BUPC & NI, i] <- m_TDP_2[EP1 & BUPC & NI, i] <- m_TDP_3[EP1 & BUPC & NI, i] <- as.vector(exp(p_frailty_BUPC_NI_1 * p_weibull_scale_BUPC_NI * (((i - 1)^p_weibull_shape_BUPC_NI) - (i^p_weibull_shape_BUPC_NI))))
        #m_TDP_1[EP1 & MET & NI, i]  <- m_TDP_2[EP1 & MET & NI, i]  <- m_TDP_3[EP1 & MET & NI, i]  <- as.vector(exp(p_frailty_MET_NI_1 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        #m_TDP_1[EP1 & METC & NI, i] <- m_TDP_2[EP1 & METC & NI, i] <- m_TDP_3[EP1 & METC & NI, i] <- as.vector(exp(p_frailty_METC_NI_1 * p_weibull_scale_METC_NI * (((i - 1)^p_weibull_shape_METC_NI) - (i^p_weibull_shape_METC_NI))))
        #m_TDP_1[EP1 & ABS & NI, i]  <- m_TDP_2[EP1 & ABS & NI, i]  <- m_TDP_3[EP1 & ABS & NI, i]  <- as.vector(exp(p_frailty_ABS_NI_1 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        #m_TDP_1[EP1 & REL & NI, i]  <- m_TDP_2[EP1 & REL & NI, i]  <- m_TDP_3[EP1 & REL & NI, i]  <- as.vector(exp(p_frailty_REL_NI_1 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
        #m_TDP_1[EP1 & ODN & NI, i]  <- m_ODN[EP1 & REL & NI, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
        #m_TDP_2[EP1 & ODN & NI, i]  <- m_ODN[EP1 & REL & NI, 2]
        #m_TDP_3[EP1 & ODN & NI, i]  <- m_ODN[EP1 & REL & NI, 3]
        #m_TDP_1[EP1 & ODF & NI, i]  <- m_TDP_2[EP1 & ODF & NI, i] <- m_TDP_3[EP1 & ODF & NI, i]  <- 1 # all remain in ODF
        # Episode 2
        #m_TDP_1[EP2 & BUP & NI, i] <- m_TDP_2[EP2 & BUP & NI, i] <- m_TDP_3[EP2 & BUP & NI, i] <- as.vector(exp(p_frailty_BUP_NI_1 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI)))) # vector of remain probabilities
        #m_TDP_1[EP2 & BUPC & NI, i] <- m_TDP_2[EP2 & BUPC & NI, i] <- m_TDP_3[EP2 & BUPC & NI, i] <- as.vector(exp(p_frailty_BUPC_NI_1 * p_weibull_scale_BUPC_NI * (((i - 1)^p_weibull_shape_BUPC_NI) - (i^p_weibull_shape_BUPC_NI))))
        #m_TDP_1[EP2 & MET & NI, i] <- m_TDP_2[EP2 & MET & NI, i] <- m_TDP_3[EP2 & MET & NI, i]  <- as.vector(exp(p_frailty_MET_NI_1 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        #m_TDP_1[EP2 & METC & NI, i] <- m_TDP_2[EP2 & METC & NI, i] <- m_TDP_3[EP2 & METC & NI, i] <- as.vector(exp(p_frailty_METC_NI_1 * p_weibull_scale_METC_NI * (((i - 1)^p_weibull_shape_METC_NI) - (i^p_weibull_shape_METC_NI))))
        #m_TDP_1[EP2 & ABS & NI, i] <- m_TDP_2[EP2 & ABS & NI, i] <- m_TDP_3[EP2 & ABS & NI, i]  <- as.vector(exp(p_frailty_ABS_NI_1 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        #m_TDP_1[EP2 & REL & NI, i] <- m_TDP_2[EP2 & REL & NI, i] <- m_TDP_3[EP2 & REL & NI, i]  <- as.vector(exp(p_frailty_REL_NI_1 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
        #m_TDP_1[EP2 & ODN & NI, i]  <- m_ODN[EP2 & REL & NI, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
        #m_TDP_2[EP2 & ODN & NI, i]  <- m_ODN[EP2 & REL & NI, 2]
        #m_TDP_3[EP2 & ODN & NI, i]  <- m_ODN[EP2 & REL & NI, 3]
        #m_TDP_1[EP2 & ODF & NI, i] <- m_TDP_2[EP2 & ODF & NI, i] <- m_TDP_3[EP2 & ODF & NI, i]  <- 1 # all remain in ODF
        # Episode 3
        #m_TDP_1[EP3 & BUP & NI, i] <- m_TDP_2[EP3 & BUP & NI, i] <- m_TDP_3[EP3 & BUP & NI, i] <- as.vector(exp(p_frailty_BUP_NI_1 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI)))) # vector of remain probabilities
        #m_TDP_1[EP3 & BUPC & NI, i] <- m_TDP_2[EP3 & BUPC & NI, i] <- m_TDP_3[EP3 & BUPC & NI, i] <- as.vector(exp(p_frailty_BUPC_NI_1 * p_weibull_scale_BUPC_NI * (((i - 1)^p_weibull_shape_BUPC_NI) - (i^p_weibull_shape_BUPC_NI))))
        #m_TDP_1[EP3 & MET & NI, i] <- m_TDP_2[EP3 & MET & NI, i] <- m_TDP_3[EP3 & MET & NI, i]  <- as.vector(exp(p_frailty_MET_NI_1 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        #m_TDP_1[EP3 & METC & NI, i] <- m_TDP_2[EP3 & METC & NI, i] <- m_TDP_3[EP3 & METC & NI, i] <- as.vector(exp(p_frailty_METC_NI_1 * p_weibull_scale_METC_NI * (((i - 1)^p_weibull_shape_METC_NI) - (i^p_weibull_shape_METC_NI))))
        #m_TDP_1[EP3 & ABS & NI, i] <- m_TDP_2[EP3 & ABS & NI, i] <- m_TDP_3[EP3 & ABS & NI, i]  <- as.vector(exp(p_frailty_ABS_NI_1 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        #m_TDP_1[EP3 & REL & NI, i] <- m_TDP_2[EP3 & REL & NI, i] <- m_TDP_3[EP3 & REL & NI, i]  <- as.vector(exp(p_frailty_REL_NI_1 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
        #m_TDP_1[EP3 & ODN & NI, i]  <- m_ODN[EP3 & REL & NI, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
        #m_TDP_2[EP3 & ODN & NI, i]  <- m_ODN[EP3 & REL & NI, 2]
        #m_TDP_3[EP3 & ODN & NI, i]  <- m_ODN[EP3 & REL & NI, 3]
        #m_TDP_1[EP3 & ODF & NI, i] <- m_TDP_2[EP3 & ODF & NI, i] <- m_TDP_3[EP3 & ODF & NI, i]  <- 1 # all remain in ODF
      # Injection
        # Episode 1
        #m_TDP_1[EP1 & BUP & INJ, i] <- m_TDP_2[EP1 & BUP & INJ, i] <- m_TDP_3[EP1 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_1 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ)))) # vector of remain probabilities
        #m_TDP_1[EP1 & BUPC & INJ, i] <- m_TDP_2[EP1 & BUPC & INJ, i] <- m_TDP_3[EP1 & BUPC & INJ, i] <- as.vector(exp(p_frailty_BUPC_INJ_1 * p_weibull_scale_BUPC_INJ * (((i - 1)^p_weibull_shape_BUPC_INJ) - (i^p_weibull_shape_BUPC_INJ))))
        #m_TDP_1[EP1 & MET & INJ, i] <- m_TDP_2[EP1 & MET & INJ, i] <- m_TDP_3[EP1 & MET & INJ, i]  <- as.vector(exp(p_frailty_MET_INJ_1 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        #m_TDP_1[EP1 & METC & INJ, i] <- m_TDP_2[EP1 & METC & INJ, i] <- m_TDP_3[EP1 & METC & INJ, i] <- as.vector(exp(p_frailty_METC_INJ_1 * p_weibull_scale_METC_INJ * (((i - 1)^p_weibull_shape_METC_INJ) - (i^p_weibull_shape_METC_INJ))))
        #m_TDP_1[EP1 & ABS & INJ, i] <- m_TDP_2[EP1 & ABS & INJ, i] <- m_TDP_3[EP1 & ABS & INJ, i]  <- as.vector(exp(p_frailty_ABS_INJ_1 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        #m_TDP_1[EP1 & REL & INJ, i] <- m_TDP_2[EP1 & REL & INJ, i] <- m_TDP_3[EP1 & REL & INJ, i]  <- as.vector(exp(p_frailty_REL_INJ_1 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
        #m_TDP_1[EP1 & ODN & INJ, i]  <- m_ODN[EP1 & REL & INJ, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
        #m_TDP_2[EP1 & ODN & INJ, i]  <- m_ODN[EP1 & REL & INJ, 2]
        #m_TDP_3[EP1 & ODN & INJ, i]  <- m_ODN[EP1 & REL & INJ, 3]
        #m_TDP_1[EP1 & ODF & INJ, i] <- m_TDP_2[EP1 & ODF & INJ, i] <- m_TDP_3[EP1 & ODF & INJ, i]  <- 1 # all remain in ODF
        # Episode 2
        #m_TDP_1[EP2 & BUP & INJ, i] <- m_TDP_2[EP2 & BUP & INJ, i] <- m_TDP_3[EP2 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_1 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ)))) # vector of remain probabilities
        #m_TDP_1[EP2 & BUPC & INJ, i] <- m_TDP_2[EP2 & BUPC & INJ, i] <- m_TDP_3[EP2 & BUPC & INJ, i] <- as.vector(exp(p_frailty_BUPC_INJ_1 * p_weibull_scale_BUPC_INJ * (((i - 1)^p_weibull_shape_BUPC_INJ) - (i^p_weibull_shape_BUPC_INJ))))
        #m_TDP_1[EP2 & MET & INJ, i] <- m_TDP_2[EP2 & MET & INJ, i] <- m_TDP_3[EP2 & MET & INJ, i]  <- as.vector(exp(p_frailty_MET_INJ_1 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        #m_TDP_1[EP2 & METC & INJ, i] <- m_TDP_2[EP2 & METC & INJ, i] <- m_TDP_3[EP2 & METC & INJ, i] <- as.vector(exp(p_frailty_METC_INJ_1 * p_weibull_scale_METC_INJ * (((i - 1)^p_weibull_shape_METC_INJ) - (i^p_weibull_shape_METC_INJ))))
        #m_TDP_1[EP2 & ABS & INJ, i] <- m_TDP_2[EP2 & ABS & INJ, i] <- m_TDP_3[EP2 & ABS & INJ, i]  <- as.vector(exp(p_frailty_ABS_INJ_1 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        #m_TDP_1[EP2 & REL & INJ, i] <- m_TDP_2[EP2 & REL & INJ, i] <- m_TDP_3[EP2 & REL & INJ, i]  <- as.vector(exp(p_frailty_REL_INJ_1 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
        #m_TDP_1[EP2 & ODN & INJ, i]  <- m_ODN[EP2 & REL & INJ, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
        #m_TDP_2[EP2 & ODN & INJ, i]  <- m_ODN[EP2 & REL & INJ, 2]
        #m_TDP_3[EP2 & ODN & INJ, i]  <- m_ODN[EP2 & REL & INJ, 3]
        #m_TDP_1[EP2 & ODF & INJ, i] <- m_TDP_2[EP2 & ODF & INJ, i] <- m_TDP_3[EP2 & ODF & INJ, i]  <- 1 # all remain in ODF
        # Episode 3
        #m_TDP_1[EP3 & BUP & INJ, i] <- m_TDP_2[EP3 & BUP & INJ, i] <- m_TDP_3[EP3 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_1 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ)))) # vector of remain probabilities
        #m_TDP_1[EP3 & BUPC & INJ, i] <- m_TDP_2[EP3 & BUPC & INJ, i] <- m_TDP_3[EP3 & BUPC & INJ, i] <- as.vector(exp(p_frailty_BUPC_INJ_1 * p_weibull_scale_BUPC_INJ * (((i - 1)^p_weibull_shape_BUPC_INJ) - (i^p_weibull_shape_BUPC_INJ))))
        #m_TDP_1[EP3 & MET & INJ, i] <- m_TDP_2[EP3 & MET & INJ, i] <- m_TDP_3[EP3 & MET & INJ, i]  <- as.vector(exp(p_frailty_MET_INJ_1 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        #m_TDP_1[EP3 & METC & INJ, i] <- m_TDP_2[EP3 & METC & INJ, i] <- m_TDP_3[EP3 & METC & INJ, i] <- as.vector(exp(p_frailty_METC_INJ_1 * p_weibull_scale_METC_INJ * (((i - 1)^p_weibull_shape_METC_INJ) - (i^p_weibull_shape_METC_INJ))))
        #m_TDP_1[EP3 & ABS & INJ, i] <- m_TDP_2[EP3 & ABS & INJ, i] <- m_TDP_3[EP3 & ABS & INJ, i]  <- as.vector(exp(p_frailty_ABS_INJ_1 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        #m_TDP_1[EP3 & REL & INJ, i] <- m_TDP_2[EP3 & REL & INJ, i] <- m_TDP_3[EP3 & REL & INJ, i]  <- as.vector(exp(p_frailty_REL_INJ_1 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
        #m_TDP_1[EP3 & ODN & INJ, i]  <- m_ODN[EP3 & REL & INJ, 1] #use same probability as overdose from relapse (no episode multipliers at this point)
        #m_TDP_2[EP3 & ODN & INJ, i]  <- m_ODN[EP3 & REL & INJ, 2]
        #m_TDP_3[EP3 & ODN & INJ, i]  <- m_ODN[EP3 & REL & INJ, 3]
        #m_TDP_1[EP3 & ODF & INJ, i] <- m_TDP_2[EP3 & ODF & INJ, i] <- m_TDP_3[EP3 & ODF & INJ, i]  <- 1 # all remain in ODF
    #}
  # Probability of state-exit
  m_leave_1 <- 1 - m_TDP_1
  m_leave_2 <- 1 - m_TDP_2
  m_leave_3 <- 1 - m_TDP_3
  
  if(checks){
    # Time dependent state-exit probabilities (from weibull estimates)
    write.csv(m_TDP_1,"checks/state-time dependent transitions/m_TDP_1.csv", row.names = TRUE)
    write.csv(m_TDP_2,"checks/state-time dependent transitions/m_TDP_2.csv", row.names = TRUE)
    write.csv(m_TDP_3,"checks/state-time dependent transitions/m_TDP_3.csv", row.names = TRUE)
    write.csv(m_leave_1,"checks/state-time dependent transitions/m_leave_1.csv", row.names = TRUE)
    write.csv(m_leave_2,"checks/state-time dependent transitions/m_leave_2.csv", row.names = TRUE)
    write.csv(m_leave_3,"checks/state-time dependent transitions/m_leave_3.csv", row.names = TRUE)
  } else{}
  
  #### Mortality ####
  #' Mortality probability estimates
  #'
  #' \code{v_mort} is used to populate mortality probability vectors.
  #'
  #' @param hr
  #' @param per
  #' @return 
  #' Mortality vectors for each age applied to model periods (months or weeks), includes state-specific hr.
  #' Overdose deaths tracked as "ODF"
  #' @export
  v_mort <- function(hr = hr, per = per){
    v_mort <- rep((1 - exp(-v_r_mort_by_age[n_age_init:(n_age_max - 1), ] * (1/per) * hr)), each = per) #currently working in months
    return(v_mort)
  }
  # Non-injection
  v_mort_BUP_NI     <- v_mort(hr = hr_BUP_NI, per = n_per)
  v_mort_BUPC_NI    <- v_mort(hr = hr_BUPC_NI, per = n_per)
  v_mort_MET_NI     <- v_mort(hr = hr_MET_NI, per = n_per)
  v_mort_METC_NI    <- v_mort(hr = hr_METC_NI, per = n_per)
  v_mort_REL_NI     <- v_mort(hr = hr_REL_NI, per = n_per)
  v_mort_ODN_NI     <- v_mort(hr = hr_ODN_NI, per = n_per) # Mortality equal for REL/ODN (fatal overdoses counted separately)
  v_mort_ODF_NI     <- rep(0, n_t) # stay in ODF, not tracked in "death"
  v_mort_ABS_NEG_NI <- v_mort(hr = hr_ABS_NI, per = n_per)
  v_mort_ABS_HIV_NI <- v_mort(hr = hr_HIV_NI, per = n_per)
  v_mort_ABS_HCV_NI <- v_mort(hr = hr_HCV_NI, per = n_per)
  v_mort_ABS_COI_NI <- v_mort(hr = hr_COI_NI, per = n_per)
  
  # Injection
  v_mort_BUP_INJ     <- v_mort(hr = hr_BUP_INJ, per = n_per)
  v_mort_BUPC_INJ    <- v_mort(hr = hr_BUPC_INJ, per = n_per)
  v_mort_MET_INJ     <- v_mort(hr = hr_MET_INJ, per = n_per)
  v_mort_METC_INJ    <- v_mort(hr = hr_METC_INJ, per = n_per)
  v_mort_REL_INJ     <- v_mort(hr = hr_REL_INJ, per = n_per)
  v_mort_ODN_INJ     <- v_mort(hr = hr_ODN_INJ, per = n_per) # Mortality equal for REL/ODN (fatal overdoses counted separately)
  v_mort_ODF_INJ     <- rep(0, n_t) # stay in ODF, not tracked in "death"
  v_mort_ABS_NEG_INJ <- v_mort(hr = hr_ABS_INJ, per = n_per)
  v_mort_ABS_HIV_INJ <- v_mort(hr = hr_HIV_INJ, per = n_per)
  v_mort_ABS_HCV_INJ <- v_mort(hr = hr_HCV_INJ, per = n_per)
  v_mort_ABS_COI_INJ <- v_mort(hr = hr_COI_INJ, per = n_per)

  # Create empty mortality matrix
  m_mort <- array(0, dim = c(n_states, n_t),
                  dimnames = list(v_n_states, 1:n_t))
  # Populate mortality matrix (death probability from each state)
  for (i in 1:n_t){
    # Non-injection
    m_mort[BUP & NI, i]       <- v_mort_BUP_NI[i]
    m_mort[BUPC & NI, i]      <- v_mort_BUPC_NI[i]
    m_mort[MET & NI, i]       <- v_mort_MET_NI[i]
    m_mort[METC & NI, i]      <- v_mort_METC_NI[i]
    m_mort[REL & NI, i]       <- v_mort_REL_NI[i]
    m_mort[ODN & NI, i]       <- v_mort_ODN_NI[i] # using background excess mortality for relapse in non-fatal overdose
    m_mort[ODF & NI, i]       <- v_mort_ODF_NI[i] # transition to death = 0
    m_mort[ABS & NI & NEG, i] <- v_mort_ABS_NEG_NI[i]
    m_mort[ABS & NI & HIV, i] <- v_mort_ABS_HIV_NI[i]
    m_mort[ABS & NI & HCV, i] <- v_mort_ABS_HCV_NI[i]
    m_mort[ABS & NI & COI, i] <- v_mort_ABS_COI_NI[i]
    # Injection
    m_mort[BUP & INJ, i]       <- v_mort_BUP_INJ[i]
    m_mort[BUPC & INJ, i]      <- v_mort_BUPC_INJ[i]
    m_mort[MET & INJ, i]       <- v_mort_MET_INJ[i]
    m_mort[METC & INJ, i]      <- v_mort_METC_INJ[i]
    m_mort[REL & INJ, i]       <- v_mort_REL_INJ[i]
    m_mort[ODN & INJ, i]       <- v_mort_ODN_INJ[i] # using background excess mortality for relapse in non-fatal overdose
    m_mort[ODF & INJ, i]       <- v_mort_ODF_INJ[i] # transition to death = 0
    m_mort[ABS & INJ & NEG, i] <- v_mort_ABS_NEG_INJ[i]
    m_mort[ABS & INJ & HIV, i] <- v_mort_ABS_HIV_INJ[i]
    m_mort[ABS & INJ & HCV, i] <- v_mort_ABS_HCV_INJ[i]
    m_mort[ABS & INJ & COI, i] <- v_mort_ABS_COI_INJ[i]
  }
  
  if(checks){
    # Mortality matrix
    write.csv(m_mort,"checks/mortality/m_mort.csv", row.names = TRUE)
  } else{}

  # Alive probability in each period
  m_alive <- 1 - m_mort

  #### Unconditional transition probabilities ####
  # Empty 2-D unconditional transition matrix (from states, to states)
  # Create 2 matrices: 1 - First four weeks (higher OD prob)
  #m_UP <- m_UP_4wk <- array(0, dim = c(n_states, n_states),
  #                          dimnames = list(v_n_states, v_n_states))
  
  # Create as array (different probabilities for model-time varying overdose)
  a_UP <- a_UP_first <- array(0, dim = c(n_states, n_states, time_periods),
                              dimnames = list(v_n_states, v_n_states, 1:time_periods))
  # Populate unconditional transition matrix
  # Overdose probability populated first, accounting for higher probability of overdose transition in first month
  for (i in 1:time_periods){
  # Non-Injection
  # From BUP
  # First month
  a_UP_first[BUP & NI, BUPC & NI, i] <- p_BUP_BUPC_NI * (1 - m_ODN_first[BUP & NI, i] - m_ODF_first[BUP & NI, i])
  a_UP_first[BUP & NI, MET & NI, i]  <- p_BUP_MET_NI * (1 - m_ODN_first[BUP & NI, i] - m_ODF_first[BUP & NI, i])
  a_UP_first[BUP & NI, METC & NI, i] <- p_BUP_METC_NI * (1 - m_ODN_first[BUP & NI, i] - m_ODF_first[BUP & NI, i])
  a_UP_first[BUP & NI, ABS & NI, i]  <- p_BUP_ABS_NI * (1 - m_ODN_first[BUP & NI, i] - m_ODF_first[BUP & NI, i])
  a_UP_first[BUP & NI, REL & NI, i]  <- p_BUP_REL_NI * (1 - m_ODN_first[BUP & NI, i] - m_ODF_first[BUP & NI, i])
  a_UP_first[BUP & NI, ODN & NI, i]  <- m_ODN_first[BUP & NI, i]
  a_UP_first[BUP & NI, ODF & NI, i]  <- m_ODF_first[BUP & NI, i]
  
  # Month 2+
  a_UP[BUP & NI, BUPC & NI, i] <- p_BUP_BUPC_NI * (1 - m_ODN[BUP & NI, i] - m_ODF[BUP & NI, i])
  a_UP[BUP & NI, MET & NI, i]  <- p_BUP_MET_NI * (1 - m_ODN[BUP & NI, i] - m_ODF[BUP & NI, i])
  a_UP[BUP & NI, METC & NI, i] <- p_BUP_METC_NI * (1 - m_ODN[BUP & NI, i] - m_ODF[BUP & NI, i])
  a_UP[BUP & NI, ABS & NI, i]  <- p_BUP_ABS_NI * (1 - m_ODN[BUP & NI, i] - m_ODF[BUP & NI, i])
  a_UP[BUP & NI, REL & NI, i]  <- p_BUP_REL_NI * (1 - m_ODN[BUP & NI, i] - m_ODF[BUP & NI, i])
  a_UP[BUP & NI, ODN & NI, i]  <- m_ODN[BUP & NI, i]
  a_UP[BUP & NI, ODF & NI, i]  <- m_ODF[BUP & NI, i]

  # From BUPC
  # First month
  a_UP_first[BUPC & NI, BUP & NI, i]  <- p_BUPC_BUP_NI * (1 - m_ODN_first[BUPC & NI, i] - m_ODF_first[BUPC & NI, i])
  a_UP_first[BUPC & NI, MET & NI, i]  <- p_BUPC_MET_NI * (1 - m_ODN_first[BUPC & NI, i] - m_ODF_first[BUPC & NI, i])
  a_UP_first[BUPC & NI, METC & NI, i] <- p_BUPC_METC_NI * (1 - m_ODN_first[BUPC & NI, i] - m_ODF_first[BUPC & NI, i])
  a_UP_first[BUPC & NI, ABS & NI, i]  <- p_BUPC_ABS_NI * (1 - m_ODN_first[BUPC & NI, i] - m_ODF_first[BUPC & NI, i])
  a_UP_first[BUPC & NI, REL & NI, i]  <- p_BUPC_REL_NI * (1 - m_ODN_first[BUPC & NI, i] - m_ODF_first[BUPC & NI, i])
  a_UP_first[BUPC & NI, ODN & NI, i]  <- m_ODN_first[BUPC & NI, i]
  a_UP_first[BUPC & NI, ODF & NI, i]  <- m_ODF_first[BUPC & NI, i]
  
  # Month 2+
  a_UP[BUPC & NI, BUP & NI, i]  <- p_BUPC_BUP_NI * (1 - m_ODN[BUPC & NI, i] - m_ODF[BUPC & NI, i])
  a_UP[BUPC & NI, MET & NI, i]  <- p_BUPC_MET_NI * (1 - m_ODN[BUPC & NI, i] - m_ODF[BUPC & NI, i])
  a_UP[BUPC & NI, METC & NI, i] <- p_BUPC_METC_NI * (1 - m_ODN[BUPC & NI, i] - m_ODF[BUPC & NI, i])
  a_UP[BUPC & NI, ABS & NI, i]  <- p_BUPC_ABS_NI * (1 - m_ODN[BUPC & NI, i] - m_ODF[BUPC & NI, i])
  a_UP[BUPC & NI, REL & NI, i]  <- p_BUPC_REL_NI * (1 - m_ODN[BUPC & NI, i] - m_ODF[BUPC & NI, i])
  a_UP[BUPC & NI, ODN & NI, i]  <- m_ODN[BUPC & NI, i]
  a_UP[BUPC & NI, ODF & NI, i]  <- m_ODF[BUPC & NI, i]
  
  # From MET
  # First month
  a_UP_first[MET & NI, METC & NI, i] <- p_MET_METC_NI * (1 - m_ODN_first[MET & NI, i] - m_ODF_first[MET & NI, i])
  a_UP_first[MET & NI, BUP & NI, i]  <- p_MET_BUP_NI * (1 - m_ODN_first[MET & NI, i] - m_ODF_first[MET & NI, i])
  a_UP_first[MET & NI, BUPC & NI, i]  <- p_MET_BUPC_NI * (1 - m_ODN_first[MET & NI, i] - m_ODF_first[MET & NI, i])
  a_UP_first[MET & NI, ABS & NI, i]  <- p_MET_ABS_NI * (1 - m_ODN_first[MET & NI, i] - m_ODF_first[MET & NI, i])
  a_UP_first[MET & NI, REL & NI, i]  <- p_MET_REL_NI * (1 - m_ODN_first[MET & NI, i] - m_ODF_first[MET & NI, i])
  a_UP_first[MET & NI, ODN & NI, i]  <- m_ODN_first[MET & NI, i]
  a_UP_first[MET & NI, ODF & NI, i]  <- m_ODF_first[MET & NI, i]
  
  # Month 2+
  a_UP[MET & NI, METC & NI, i] <- p_MET_METC_NI * (1 - m_ODN[MET & NI, i] - m_ODF[MET & NI, i])
  a_UP[MET & NI, BUP & NI, i]  <- p_MET_BUP_NI * (1 - m_ODN[MET & NI, i] - m_ODF[MET & NI, i])
  a_UP[MET & NI, BUPC & NI, i]  <- p_MET_BUPC_NI * (1 - m_ODN[MET & NI, i] - m_ODF[MET & NI, i])
  a_UP[MET & NI, ABS & NI, i]  <- p_MET_ABS_NI * (1 - m_ODN[MET & NI, i] - m_ODF[MET & NI, i])
  a_UP[MET & NI, REL & NI, i]  <- p_MET_REL_NI * (1 - m_ODN[MET & NI, i] - m_ODF[MET & NI, i])
  a_UP[MET & NI, ODN & NI, i]  <- m_ODN[MET & NI, i]
  a_UP[MET & NI, ODF & NI, i]  <- m_ODF[MET & NI, i]

  # From METC
  # First month
  a_UP_first[METC & NI, MET & NI, i]  <- p_METC_MET_NI * (1 - m_ODN_first[METC & NI, i] - m_ODF_first[METC & NI, i])
  a_UP_first[METC & NI, BUP & NI, i]  <- p_METC_BUP_NI * (1 - m_ODN_first[METC & NI, i] - m_ODF_first[METC & NI, i])
  a_UP_first[METC & NI, BUPC & NI, i] <- p_METC_BUPC_NI * (1 - m_ODN_first[METC & NI, i] - m_ODF_first[METC & NI, i])
  a_UP_first[METC & NI, ABS & NI, i]  <- p_METC_ABS_NI * (1 - m_ODN_first[METC & NI, i] - m_ODF_first[METC & NI, i])
  a_UP_first[METC & NI, REL & NI, i]  <- p_METC_REL_NI * (1 - m_ODN_first[METC & NI, i] - m_ODF_first[METC & NI, i])
  a_UP_first[METC & NI, ODN & NI, i]  <- m_ODN_first[METC & NI, i]
  a_UP_first[METC & NI, ODF & NI, i]  <- m_ODF_first[METC & NI, i]
  
  # Month 2+
  a_UP[METC & NI, MET & NI, i]  <- p_METC_MET_NI * (1 - m_ODN[METC & NI, i] - m_ODF[METC & NI, i])
  a_UP[METC & NI, BUP & NI, i]  <- p_METC_BUP_NI * (1 - m_ODN[METC & NI, i] - m_ODF[METC & NI, i])
  a_UP[METC & NI, BUPC & NI, i] <- p_METC_BUPC_NI * (1 - m_ODN[METC & NI, i] - m_ODF[METC & NI, i])
  a_UP[METC & NI, ABS & NI, i]  <- p_METC_ABS_NI * (1 - m_ODN[METC & NI, i] - m_ODF[METC & NI, i])
  a_UP[METC & NI, REL & NI, i]  <- p_METC_REL_NI * (1 - m_ODN[METC & NI, i] - m_ODF[METC & NI, i])
  a_UP[METC & NI, ODN & NI, i]  <- m_ODN[METC & NI, i]
  a_UP[METC & NI, ODF & NI, i]  <- m_ODF[METC & NI, i]
  
  # From ABS
  # First month
  a_UP_first[ABS & NI, REL & NI, i]  <- p_ABS_REL_NI * (1 - m_ODN_first[ABS & NI, i] - m_ODF_first[ABS & NI, i])
  a_UP_first[ABS & NI, MET & NI, i]  <- p_ABS_MET_NI * (1 - m_ODN_first[ABS & NI, i] - m_ODF_first[ABS & NI, i])
  a_UP_first[ABS & NI, METC & NI, i] <- p_ABS_METC_NI * (1 - m_ODN_first[ABS & NI, i] - m_ODF_first[ABS & NI, i])
  a_UP_first[ABS & NI, BUP & NI, i]  <- p_ABS_BUP_NI * (1 - m_ODN_first[ABS & NI, i] - m_ODF_first[ABS & NI, i])
  a_UP_first[ABS & NI, BUPC & NI, i] <- p_ABS_BUPC_NI * (1 - m_ODN_first[ABS & NI, i] - m_ODF_first[ABS & NI, i])
  a_UP_first[ABS & NI, ODN & NI, i]  <- m_ODN_first[ABS & NI, i]
  a_UP_first[ABS & NI, ODF & NI, i]  <- m_ODF_first[ABS & NI, i]
  
  # Month 2+
  a_UP[ABS & NI, REL & NI, i]  <- p_ABS_REL_NI * (1 - m_ODN[ABS & NI, i] - m_ODF[ABS & NI, i])
  a_UP[ABS & NI, MET & NI, i]  <- p_ABS_MET_NI * (1 - m_ODN[ABS & NI, i] - m_ODF[ABS & NI, i])
  a_UP[ABS & NI, METC & NI, i] <- p_ABS_METC_NI * (1 - m_ODN[ABS & NI, i] - m_ODF[ABS & NI, i])
  a_UP[ABS & NI, BUP & NI, i]  <- p_ABS_BUP_NI * (1 - m_ODN[ABS & NI, i] - m_ODF[ABS & NI, i])
  a_UP[ABS & NI, BUPC & NI, i] <- p_ABS_BUPC_NI * (1 - m_ODN[ABS & NI, i] - m_ODF[ABS & NI, i])
  a_UP[ABS & NI, ODN & NI, i]  <- m_ODN[ABS & NI, i]
  a_UP[ABS & NI, ODF & NI, i]  <- m_ODF[ABS & NI, i]
  
  # From REL
  # First month
  a_UP_first[REL & NI, MET & NI, i]  <- p_REL_MET_NI * (1 - m_ODN_first[REL & NI, i] - m_ODF_first[REL & NI, i])
  a_UP_first[REL & NI, METC & NI, i] <- p_REL_METC_NI * (1 - m_ODN_first[REL & NI, i] - m_ODF_first[REL & NI, i])
  a_UP_first[REL & NI, BUP & NI, i]  <- p_REL_BUP_NI * (1 - m_ODN_first[REL & NI, i] - m_ODF_first[REL & NI, i])
  a_UP_first[REL & NI, BUPC & NI, i] <- p_REL_BUPC_NI * (1 - m_ODN_first[REL & NI, i] - m_ODF_first[REL & NI, i])
  a_UP_first[REL & NI, ABS & NI, i]  <- p_REL_ABS_NI * (1 - m_ODN_first[REL & NI, i] - m_ODF_first[REL & NI, i])
  a_UP_first[REL & NI, ODN & NI, i]  <- m_ODN_first[REL & NI, i]
  a_UP_first[REL & NI, ODF & NI, i]  <- m_ODF_first[REL & NI, i]
  
  # Month 2+
  a_UP[REL & NI, MET & NI, i]  <- p_REL_MET_NI * (1 - m_ODN[REL & NI, i] - m_ODF[REL & NI, i])
  a_UP[REL & NI, METC & NI, i] <- p_REL_METC_NI * (1 - m_ODN[REL & NI, i] - m_ODF[REL & NI, i])
  a_UP[REL & NI, BUP & NI, i]  <- p_REL_BUP_NI * (1 - m_ODN[REL & NI, i] - m_ODF[REL & NI, i])
  a_UP[REL & NI, BUPC & NI, i] <- p_REL_BUPC_NI * (1 - m_ODN[REL & NI, i] - m_ODF[REL & NI, i])
  a_UP[REL & NI, ABS & NI, i]  <- p_REL_ABS_NI * (1 - m_ODN[REL & NI, i] - m_ODF[REL & NI, i])
  a_UP[REL & NI, ODN & NI, i]  <- m_ODN[REL & NI, i]
  a_UP[REL & NI, ODF & NI, i]  <- m_ODF[REL & NI, i]

  # From OD (first month same)
  a_UP[ODN & NI, MET & NI, i]  <- a_UP_first[ODN & NI, MET & NI, i] <- p_ODN_MET_NI
  a_UP[ODN & NI, METC & NI, i] <- a_UP_first[ODN & NI, METC & NI, i] <- p_ODN_METC_NI
  a_UP[ODN & NI, BUP & NI, i]  <- a_UP_first[ODN & NI, BUP & NI, i] <- p_ODN_BUP_NI
  a_UP[ODN & NI, BUPC & NI, i] <- a_UP_first[ODN & NI, BUPC & NI, i] <- p_ODN_BUPC_NI
  a_UP[ODN & NI, ABS & NI, i]  <- a_UP_first[ODN & NI, ABS & NI, i] <- p_ODN_ABS_NI
  a_UP[ODN & NI, REL & NI, i]  <- a_UP_first[ODN & NI, REL & NI, i] <- p_ODN_REL_NI

  # Injection
  # From BUP
  # First month
  a_UP_first[BUP & INJ, BUPC & INJ, i] <- p_BUP_BUPC_INJ * (1 - m_ODN_first[BUP & INJ, i] - m_ODF_first[BUP & INJ, i])
  a_UP_first[BUP & INJ, MET & INJ, i]  <- p_BUP_MET_INJ * (1 - m_ODN_first[BUP & INJ, i] - m_ODF_first[BUP & INJ, i])
  a_UP_first[BUP & INJ, METC & INJ, i] <- p_BUP_METC_INJ * (1 - m_ODN_first[BUP & INJ, i] - m_ODF_first[BUP & INJ, i])
  a_UP_first[BUP & INJ, ABS & INJ, i]  <- p_BUP_ABS_INJ * (1 - m_ODN_first[BUP & INJ, i] - m_ODF_first[BUP & INJ, i])
  a_UP_first[BUP & INJ, REL & INJ, i]  <- p_BUP_REL_INJ * (1 - m_ODN_first[BUP & INJ, i] - m_ODF_first[BUP & INJ, i])
  a_UP_first[BUP & INJ, ODN & INJ, i]  <- m_ODN_first[BUP & INJ, i]
  a_UP_first[BUP & INJ, ODF & INJ, i]  <- m_ODF_first[BUP & INJ, i]
  
  # Month 2+
  a_UP[BUP & INJ, BUPC & INJ, i] <- p_BUP_BUPC_INJ * (1 - m_ODN[BUP & INJ, i] - m_ODF[BUP & INJ, i])
  a_UP[BUP & INJ, MET & INJ, i]  <- p_BUP_MET_INJ * (1 - m_ODN[BUP & INJ, i] - m_ODF[BUP & INJ, i])
  a_UP[BUP & INJ, METC & INJ, i] <- p_BUP_METC_INJ * (1 - m_ODN[BUP & INJ, i] - m_ODF[BUP & INJ, i])
  a_UP[BUP & INJ, ABS & INJ, i]  <- p_BUP_ABS_INJ * (1 - m_ODN[BUP & INJ, i] - m_ODF[BUP & INJ, i])
  a_UP[BUP & INJ, REL & INJ, i]  <- p_BUP_REL_INJ * (1 - m_ODN[BUP & INJ, i] - m_ODF[BUP & INJ, i])
  a_UP[BUP & INJ, ODN & INJ, i]  <- m_ODN[BUP & INJ, i]
  a_UP[BUP & INJ, ODF & INJ, i]  <- m_ODF[BUP & INJ, i]
  
  # From BUPC
  # First month
  a_UP_first[BUPC & INJ, BUP & INJ, i]  <- p_BUPC_BUP_INJ * (1 - m_ODN_first[BUPC & INJ, i] - m_ODF_first[BUPC & INJ, i])
  a_UP_first[BUPC & INJ, MET & INJ, i]  <- p_BUPC_MET_INJ * (1 - m_ODN_first[BUPC & INJ, i] - m_ODF_first[BUPC & INJ, i])
  a_UP_first[BUPC & INJ, METC & INJ, i] <- p_BUPC_METC_INJ * (1 - m_ODN_first[BUPC & INJ, i] - m_ODF_first[BUPC & INJ, i])
  a_UP_first[BUPC & INJ, ABS & INJ, i]  <- p_BUPC_ABS_INJ * (1 - m_ODN_first[BUPC & INJ, i] - m_ODF_first[BUPC & INJ, i])
  a_UP_first[BUPC & INJ, REL & INJ, i]  <- p_BUPC_REL_INJ * (1 - m_ODN_first[BUPC & INJ, i] - m_ODF_first[BUPC & INJ, i])
  a_UP_first[BUPC & INJ, ODN & INJ, i]  <- m_ODN_first[BUPC & INJ, i]
  a_UP_first[BUPC & INJ, ODF & INJ, i]  <- m_ODF_first[BUPC & INJ, i]
  
  # Month 2+
  a_UP[BUPC & INJ, BUP & INJ, i]  <- p_BUPC_BUP_INJ * (1 - m_ODN[BUPC & INJ, i] - m_ODF[BUPC & INJ, i])
  a_UP[BUPC & INJ, MET & INJ, i]  <- p_BUPC_MET_INJ * (1 - m_ODN[BUPC & INJ, i] - m_ODF[BUPC & INJ, i])
  a_UP[BUPC & INJ, METC & INJ, i] <- p_BUPC_METC_INJ * (1 - m_ODN[BUPC & INJ, i] - m_ODF[BUPC & INJ, i])
  a_UP[BUPC & INJ, ABS & INJ, i]  <- p_BUPC_ABS_INJ * (1 - m_ODN[BUPC & INJ, i] - m_ODF[BUPC & INJ, i])
  a_UP[BUPC & INJ, REL & INJ, i]  <- p_BUPC_REL_INJ * (1 - m_ODN[BUPC & INJ, i] - m_ODF[BUPC & INJ, i])
  a_UP[BUPC & INJ, ODN & INJ, i]  <- m_ODN[BUPC & INJ, i]
  a_UP[BUPC & INJ, ODF & INJ, i]  <- m_ODF[BUPC & INJ, i]
  
  # From MET
  # First month
  a_UP_first[MET & INJ, METC & INJ, i] <- p_MET_METC_INJ * (1 - m_ODN_first[MET & INJ, i] - m_ODF_first[MET & INJ, i])
  a_UP_first[MET & INJ, BUP & INJ, i]  <- p_MET_BUP_INJ * (1 - m_ODN_first[MET & INJ, i] - m_ODF_first[MET & INJ, i])
  a_UP_first[MET & INJ, BUPC & INJ, i]  <- p_MET_BUPC_INJ * (1 - m_ODN_first[MET & INJ, i] - m_ODF_first[MET & INJ, i])
  a_UP_first[MET & INJ, ABS & INJ, i]  <- p_MET_ABS_INJ * (1 - m_ODN_first[MET & INJ, i] - m_ODF_first[MET & INJ, i])
  a_UP_first[MET & INJ, REL & INJ, i]  <- p_MET_REL_INJ * (1 - m_ODN_first[MET & INJ, i] - m_ODF_first[MET & INJ, i])
  a_UP_first[MET & INJ, ODN & INJ, i]  <- m_ODN_first[MET & INJ, i]
  a_UP_first[MET & INJ, ODF & INJ, i]  <- m_ODF_first[MET & INJ, i]
  
  # Month 2+
  a_UP[MET & INJ, METC & INJ, i] <- p_MET_METC_INJ * (1 - m_ODN[MET & INJ, i] - m_ODF[MET & INJ, i])
  a_UP[MET & INJ, BUP & INJ, i]  <- p_MET_BUP_INJ * (1 - m_ODN[MET & INJ, i] - m_ODF[MET & INJ, i])
  a_UP[MET & INJ, BUPC & INJ, i]  <- p_MET_BUPC_INJ * (1 - m_ODN[MET & INJ, i] - m_ODF[MET & INJ, i])
  a_UP[MET & INJ, ABS & INJ, i]  <- p_MET_ABS_INJ * (1 - m_ODN[MET & INJ, i] - m_ODF[MET & INJ, i])
  a_UP[MET & INJ, REL & INJ, i]  <- p_MET_REL_INJ * (1 - m_ODN[MET & INJ, i] - m_ODF[MET & INJ, i])
  a_UP[MET & INJ, ODN & INJ, i]  <- m_ODN[MET & INJ, i]
  a_UP[MET & INJ, ODF & INJ, i]  <- m_ODF[MET & INJ, i]
  
  # From METC
  # First month
  a_UP_first[METC & INJ, MET & INJ, i]  <- p_METC_MET_INJ * (1 - m_ODN_first[METC & INJ, i] - m_ODF_first[METC & INJ, i])
  a_UP_first[METC & INJ, BUP & INJ, i]  <- p_METC_BUP_INJ * (1 - m_ODN_first[METC & INJ, i] - m_ODF_first[METC & INJ, i])
  a_UP_first[METC & INJ, BUPC & INJ, i] <- p_METC_BUPC_INJ * (1 - m_ODN_first[METC & INJ, i] - m_ODF_first[METC & INJ, i])
  a_UP_first[METC & INJ, ABS & INJ, i]  <- p_METC_ABS_INJ * (1 - m_ODN_first[METC & INJ, i] - m_ODF_first[METC & INJ, i])
  a_UP_first[METC & INJ, REL & INJ, i]  <- p_METC_REL_INJ * (1 - m_ODN_first[METC & INJ, i] - m_ODF_first[METC & INJ, i])
  a_UP_first[METC & INJ, ODN & INJ, i]  <- m_ODN_first[METC & INJ, i]
  a_UP_first[METC & INJ, ODF & INJ, i]  <- m_ODF_first[METC & INJ, i]
  
  # Month 2+
  a_UP[METC & INJ, MET & INJ, i]  <- p_METC_MET_INJ * (1 - m_ODN[METC & INJ, i] - m_ODF[METC & INJ, i])
  a_UP[METC & INJ, BUP & INJ, i]  <- p_METC_BUP_INJ * (1 - m_ODN[METC & INJ, i] - m_ODF[METC & INJ, i])
  a_UP[METC & INJ, BUPC & INJ, i] <- p_METC_BUPC_INJ * (1 - m_ODN[METC & INJ, i] - m_ODF[METC & INJ, i])
  a_UP[METC & INJ, ABS & INJ, i]  <- p_METC_ABS_INJ * (1 - m_ODN[METC & INJ, i] - m_ODF[METC & INJ, i])
  a_UP[METC & INJ, REL & INJ, i]  <- p_METC_REL_INJ * (1 - m_ODN[METC & INJ, i] - m_ODF[METC & INJ, i])
  a_UP[METC & INJ, ODN & INJ, i]  <- m_ODN[METC & INJ, i]
  a_UP[METC & INJ, ODF & INJ, i]  <- m_ODF[METC & INJ, i]
  
  # From ABS
  # First month
  a_UP_first[ABS & INJ, REL & INJ, i]  <- p_ABS_REL_INJ * (1 - m_ODN_first[ABS & INJ, i] - m_ODF_first[ABS & INJ, i])
  a_UP_first[ABS & INJ, MET & INJ, i]  <- p_ABS_MET_INJ * (1 - m_ODN_first[ABS & INJ, i] - m_ODF_first[ABS & INJ, i])
  a_UP_first[ABS & INJ, METC & INJ, i] <- p_ABS_METC_INJ * (1 - m_ODN_first[ABS & INJ, i] - m_ODF_first[ABS & INJ, i])
  a_UP_first[ABS & INJ, BUP & INJ, i]  <- p_ABS_BUP_INJ * (1 - m_ODN_first[ABS & INJ, i] - m_ODF_first[ABS & INJ, i])
  a_UP_first[ABS & INJ, BUPC & INJ, i] <- p_ABS_BUPC_INJ * (1 - m_ODN_first[ABS & INJ, i] - m_ODF_first[ABS & INJ, i])
  a_UP_first[ABS & INJ, ODN & INJ, i]  <- m_ODN_first[ABS & INJ, i]
  a_UP_first[ABS & INJ, ODF & INJ, i]  <- m_ODF_first[ABS & INJ, i]
  
  # Month 2+
  a_UP[ABS & INJ, REL & INJ, i]  <- p_ABS_REL_INJ * (1 - m_ODN[ABS & INJ, i] - m_ODF[ABS & INJ, i])
  a_UP[ABS & INJ, MET & INJ, i]  <- p_ABS_MET_INJ * (1 - m_ODN[ABS & INJ, i] - m_ODF[ABS & INJ, i])
  a_UP[ABS & INJ, METC & INJ, i] <- p_ABS_METC_INJ * (1 - m_ODN[ABS & INJ, i] - m_ODF[ABS & INJ, i])
  a_UP[ABS & INJ, BUP & INJ, i]  <- p_ABS_BUP_INJ * (1 - m_ODN[ABS & INJ, i] - m_ODF[ABS & INJ, i])
  a_UP[ABS & INJ, BUPC & INJ, i] <- p_ABS_BUPC_INJ * (1 - m_ODN[ABS & INJ, i] - m_ODF[ABS & INJ, i])
  a_UP[ABS & INJ, ODN & INJ, i]  <- m_ODN[ABS & INJ, i]
  a_UP[ABS & INJ, ODF & INJ, i]  <- m_ODF[ABS & INJ, i]
  
  # From REL
  # First month
  a_UP_first[REL & INJ, MET & INJ, i]  <- p_REL_MET_INJ * (1 - m_ODN_first[REL & INJ, i] - m_ODF_first[REL & INJ, i])
  a_UP_first[REL & INJ, METC & INJ, i] <- p_REL_METC_INJ * (1 - m_ODN_first[REL & INJ, i] - m_ODF_first[REL & INJ, i])
  a_UP_first[REL & INJ, BUP & INJ, i]  <- p_REL_BUP_INJ * (1 - m_ODN_first[REL & INJ, i] - m_ODF_first[REL & INJ, i])
  a_UP_first[REL & INJ, BUPC & INJ, i] <- p_REL_BUPC_INJ * (1 - m_ODN_first[REL & INJ, i] - m_ODF_first[REL & INJ, i])
  a_UP_first[REL & INJ, ABS & INJ, i]  <- p_REL_ABS_INJ * (1 - m_ODN_first[REL & INJ, i] - m_ODF_first[REL & INJ, i])
  a_UP_first[REL & INJ, ODN & INJ, i]  <- m_ODN_first[REL & INJ, i]
  a_UP_first[REL & INJ, ODF & INJ, i]  <- m_ODF_first[REL & INJ, i]
  
  # Month 2+
  a_UP[REL & INJ, MET & INJ, i]  <- p_REL_MET_INJ * (1 - m_ODN[REL & INJ, i] - m_ODF[REL & INJ, i])
  a_UP[REL & INJ, METC & INJ, i] <- p_REL_METC_INJ * (1 - m_ODN[REL & INJ, i] - m_ODF[REL & INJ, i])
  a_UP[REL & INJ, BUP & INJ, i]  <- p_REL_BUP_INJ * (1 - m_ODN[REL & INJ, i] - m_ODF[REL & INJ, i])
  a_UP[REL & INJ, BUPC & INJ, i] <- p_REL_BUPC_INJ * (1 - m_ODN[REL & INJ, i] - m_ODF[REL & INJ, i])
  a_UP[REL & INJ, ABS & INJ, i]  <- p_REL_ABS_INJ * (1 - m_ODN[REL & INJ, i] - m_ODF[REL & INJ, i])
  a_UP[REL & INJ, ODN & INJ, i]  <- m_ODN[REL & INJ, i]
  a_UP[REL & INJ, ODF & INJ, i]  <- m_ODF[REL & INJ, i]
  
  # From OD (first month same)
  a_UP[ODN & INJ, MET & INJ, i]  <- a_UP_first[ODN & INJ, MET & INJ, i] <- p_ODN_MET_INJ
  a_UP[ODN & INJ, METC & INJ, i] <- a_UP_first[ODN & INJ, METC & INJ, i] <- p_ODN_METC_INJ
  a_UP[ODN & INJ, BUP & INJ, i]  <- a_UP_first[ODN & INJ, BUP & INJ, i] <- p_ODN_BUP_INJ
  a_UP[ODN & INJ, BUPC & INJ, i] <- a_UP_first[ODN & INJ, BUPC & INJ, i] <- p_ODN_BUPC_INJ
  a_UP[ODN & INJ, ABS & INJ, i]  <- a_UP_first[ODN & INJ, ABS & INJ, i] <- p_ODN_ABS_INJ
  a_UP[ODN & INJ, REL & INJ, i]  <- a_UP_first[ODN & INJ, REL & INJ, i] <- p_ODN_REL_INJ
  }

  #### Create full time-dependent transition array ####
  # Empty 3-D array
  a_TDP_1 <- a_TDP_2 <- a_TDP_3 <- array(0, dim = c(n_states, n_states, n_t),
                                   dimnames = list(v_n_states, v_n_states, 1:n_t))
  
  # Add transitions conditional on state-exit (m_leave = 1 - remain)
  # Three transition arrays to account for model-time-varying overdose (2018, 2019, 2020+)
  # Modified transitions for first month (state-time)
  for (i in 1){
    a_TDP_1[, , i] <- a_UP_first[, , 1] * m_leave_1[, i]
    a_TDP_2[, , i] <- a_UP_first[, , 2] * m_leave_2[, i]
    a_TDP_3[, , i] <- a_UP_first[, , 3] * m_leave_3[, i]
  }
  # All transitions 2+ months
  for (i in 2:n_t){
    #a_TDP[, , i] <- m_UP * m_leave[, i]
    a_TDP_1[, , i] <- a_UP[, , 1] * m_leave_1[, i]
    a_TDP_2[, , i] <- a_UP[, , 2] * m_leave_2[, i]
    a_TDP_3[, , i] <- a_UP[, , 3] * m_leave_3[, i]
  }

  # Add time-dependent remain probabilities
  for (i in 1:n_t){
    # Non-injection
    # Episode 1
    a_TDP_1[EP1 & BUP & NI, EP1 & BUP & NI, i] <- m_TDP_1[EP1 & BUP & NI, i]
    a_TDP_2[EP1 & BUP & NI, EP1 & BUP & NI, i] <- m_TDP_2[EP1 & BUP & NI, i]
    a_TDP_3[EP1 & BUP & NI, EP1 & BUP & NI, i] <- m_TDP_3[EP1 & BUP & NI, i]
    a_TDP_1[EP1 & BUPC & NI, EP1 & BUPC & NI, i] <- m_TDP_1[EP1 & BUPC & NI, i]
    a_TDP_2[EP1 & BUPC & NI, EP1 & BUPC & NI, i] <- m_TDP_2[EP1 & BUPC & NI, i]
    a_TDP_3[EP1 & BUPC & NI, EP1 & BUPC & NI, i] <- m_TDP_3[EP1 & BUPC & NI, i]
    a_TDP_1[EP1 & MET & NI, EP1 & MET & NI, i] <- m_TDP_1[EP1 & MET & NI, i]
    a_TDP_2[EP1 & MET & NI, EP1 & MET & NI, i] <- m_TDP_2[EP1 & MET & NI, i]
    a_TDP_3[EP1 & MET & NI, EP1 & MET & NI, i] <- m_TDP_3[EP1 & MET & NI, i]
    a_TDP_1[EP1 & METC & NI, EP1 & METC & NI, i] <- m_TDP_1[EP1 & METC & NI, i]
    a_TDP_2[EP1 & METC & NI, EP1 & METC & NI, i] <- m_TDP_2[EP1 & METC & NI, i]
    a_TDP_3[EP1 & METC & NI, EP1 & METC & NI, i] <- m_TDP_3[EP1 & METC & NI, i]
    a_TDP_1[EP1 & ABS & NI, EP1 & ABS & NI, i] <- m_TDP_1[EP1 & ABS & NI, i]
    a_TDP_2[EP1 & ABS & NI, EP1 & ABS & NI, i] <- m_TDP_2[EP1 & ABS & NI, i]
    a_TDP_3[EP1 & ABS & NI, EP1 & ABS & NI, i] <- m_TDP_3[EP1 & ABS & NI, i]
    a_TDP_1[EP1 & REL & NI, EP1 & REL & NI, i] <- m_TDP_1[EP1 & REL & NI, i]
    a_TDP_2[EP1 & REL & NI, EP1 & REL & NI, i] <- m_TDP_2[EP1 & REL & NI, i]
    a_TDP_3[EP1 & REL & NI, EP1 & REL & NI, i] <- m_TDP_3[EP1 & REL & NI, i]
    a_TDP_1[EP1 & ODN & NI, EP1 & ODN & NI, i] <- m_TDP_1[EP1 & ODN & NI, i]
    a_TDP_2[EP1 & ODN & NI, EP1 & ODN & NI, i] <- m_TDP_2[EP1 & ODN & NI, i]
    a_TDP_3[EP1 & ODN & NI, EP1 & ODN & NI, i] <- m_TDP_3[EP1 & ODN & NI, i]
    a_TDP_1[EP1 & ODF & NI, EP1 & ODF & NI, i] <- m_TDP_1[EP1 & ODF & NI, i]
    a_TDP_2[EP1 & ODF & NI, EP1 & ODF & NI, i] <- m_TDP_2[EP1 & ODF & NI, i]
    a_TDP_3[EP1 & ODF & NI, EP1 & ODF & NI, i] <- m_TDP_3[EP1 & ODF & NI, i]
    # Episode 2
    a_TDP_1[EP2 & BUP & NI, EP2 & BUP & NI, i] <- m_TDP_1[EP2 & BUP & NI, i]
    a_TDP_2[EP2 & BUP & NI, EP2 & BUP & NI, i] <- m_TDP_2[EP2 & BUP & NI, i]
    a_TDP_3[EP2 & BUP & NI, EP2 & BUP & NI, i] <- m_TDP_3[EP2 & BUP & NI, i]
    a_TDP_1[EP2 & BUPC & NI, EP2 & BUPC & NI, i] <- m_TDP_1[EP2 & BUPC & NI, i]
    a_TDP_2[EP2 & BUPC & NI, EP2 & BUPC & NI, i] <- m_TDP_2[EP2 & BUPC & NI, i]
    a_TDP_3[EP2 & BUPC & NI, EP2 & BUPC & NI, i] <- m_TDP_3[EP2 & BUPC & NI, i]
    a_TDP_1[EP2 & MET & NI, EP2 & MET & NI, i] <- m_TDP_1[EP2 & MET & NI, i]
    a_TDP_2[EP2 & MET & NI, EP2 & MET & NI, i] <- m_TDP_2[EP2 & MET & NI, i]
    a_TDP_3[EP2 & MET & NI, EP2 & MET & NI, i] <- m_TDP_3[EP2 & MET & NI, i]
    a_TDP_1[EP2 & METC & NI, EP2 & METC & NI, i] <- m_TDP_1[EP2 & METC & NI, i]
    a_TDP_2[EP2 & METC & NI, EP2 & METC & NI, i] <- m_TDP_2[EP2 & METC & NI, i]
    a_TDP_3[EP2 & METC & NI, EP2 & METC & NI, i] <- m_TDP_3[EP2 & METC & NI, i]
    a_TDP_1[EP2 & ABS & NI, EP2 & ABS & NI, i] <- m_TDP_1[EP2 & ABS & NI, i]
    a_TDP_2[EP2 & ABS & NI, EP2 & ABS & NI, i] <- m_TDP_2[EP2 & ABS & NI, i]
    a_TDP_3[EP2 & ABS & NI, EP2 & ABS & NI, i] <- m_TDP_3[EP2 & ABS & NI, i]
    a_TDP_1[EP2 & REL & NI, EP2 & REL & NI, i] <- m_TDP_1[EP2 & REL & NI, i]
    a_TDP_2[EP2 & REL & NI, EP2 & REL & NI, i] <- m_TDP_2[EP2 & REL & NI, i]
    a_TDP_3[EP2 & REL & NI, EP2 & REL & NI, i] <- m_TDP_3[EP2 & REL & NI, i]
    a_TDP_1[EP2 & ODN & NI, EP2 & ODN & NI, i] <- m_TDP_1[EP2 & ODN & NI, i]
    a_TDP_2[EP2 & ODN & NI, EP2 & ODN & NI, i] <- m_TDP_2[EP2 & ODN & NI, i]
    a_TDP_3[EP2 & ODN & NI, EP2 & ODN & NI, i] <- m_TDP_3[EP2 & ODN & NI, i]
    a_TDP_1[EP2 & ODF & NI, EP2 & ODF & NI, i] <- m_TDP_1[EP2 & ODF & NI, i]
    a_TDP_2[EP2 & ODF & NI, EP2 & ODF & NI, i] <- m_TDP_2[EP2 & ODF & NI, i]
    a_TDP_3[EP2 & ODF & NI, EP2 & ODF & NI, i] <- m_TDP_3[EP2 & ODF & NI, i]
    # Episode 3
    a_TDP_1[EP3 & BUP & NI, EP3 & BUP & NI, i] <- m_TDP_1[EP3 & BUP & NI, i]
    a_TDP_2[EP3 & BUP & NI, EP3 & BUP & NI, i] <- m_TDP_2[EP3 & BUP & NI, i]
    a_TDP_3[EP3 & BUP & NI, EP3 & BUP & NI, i] <- m_TDP_3[EP3 & BUP & NI, i]
    a_TDP_1[EP3 & BUPC & NI, EP3 & BUPC & NI, i] <- m_TDP_1[EP3 & BUPC & NI, i]
    a_TDP_2[EP3 & BUPC & NI, EP3 & BUPC & NI, i] <- m_TDP_2[EP3 & BUPC & NI, i]
    a_TDP_3[EP3 & BUPC & NI, EP3 & BUPC & NI, i] <- m_TDP_3[EP3 & BUPC & NI, i]
    a_TDP_1[EP3 & MET & NI, EP3 & MET & NI, i] <- m_TDP_1[EP3 & MET & NI, i]
    a_TDP_2[EP3 & MET & NI, EP3 & MET & NI, i] <- m_TDP_2[EP3 & MET & NI, i]
    a_TDP_3[EP3 & MET & NI, EP3 & MET & NI, i] <- m_TDP_3[EP3 & MET & NI, i]
    a_TDP_1[EP3 & METC & NI, EP3 & METC & NI, i] <- m_TDP_1[EP3 & METC & NI, i]
    a_TDP_2[EP3 & METC & NI, EP3 & METC & NI, i] <- m_TDP_2[EP3 & METC & NI, i]
    a_TDP_3[EP3 & METC & NI, EP3 & METC & NI, i] <- m_TDP_3[EP3 & METC & NI, i]
    a_TDP_1[EP3 & ABS & NI, EP3 & ABS & NI, i] <- m_TDP_1[EP3 & ABS & NI, i]
    a_TDP_2[EP3 & ABS & NI, EP3 & ABS & NI, i] <- m_TDP_2[EP3 & ABS & NI, i]
    a_TDP_3[EP3 & ABS & NI, EP3 & ABS & NI, i] <- m_TDP_3[EP3 & ABS & NI, i]
    a_TDP_1[EP3 & REL & NI, EP3 & REL & NI, i] <- m_TDP_1[EP3 & REL & NI, i]
    a_TDP_2[EP3 & REL & NI, EP3 & REL & NI, i] <- m_TDP_2[EP3 & REL & NI, i]
    a_TDP_3[EP3 & REL & NI, EP3 & REL & NI, i] <- m_TDP_3[EP3 & REL & NI, i]
    a_TDP_1[EP3 & ODN & NI, EP3 & ODN & NI, i] <- m_TDP_1[EP3 & ODN & NI, i]
    a_TDP_2[EP3 & ODN & NI, EP3 & ODN & NI, i] <- m_TDP_2[EP3 & ODN & NI, i]
    a_TDP_3[EP3 & ODN & NI, EP3 & ODN & NI, i] <- m_TDP_3[EP3 & ODN & NI, i]
    a_TDP_1[EP3 & ODF & NI, EP3 & ODF & NI, i] <- m_TDP_1[EP3 & ODF & NI, i]
    a_TDP_2[EP3 & ODF & NI, EP3 & ODF & NI, i] <- m_TDP_2[EP3 & ODF & NI, i]
    a_TDP_3[EP3 & ODF & NI, EP3 & ODF & NI, i] <- m_TDP_3[EP3 & ODF & NI, i]

    # Injection
    # Episode 1
    a_TDP_1[EP1 & BUP & INJ, EP1 & BUP & INJ, i] <- m_TDP_1[EP1 & BUP & INJ, i]
    a_TDP_2[EP1 & BUP & INJ, EP1 & BUP & INJ, i] <- m_TDP_2[EP1 & BUP & INJ, i]
    a_TDP_3[EP1 & BUP & INJ, EP1 & BUP & INJ, i] <- m_TDP_3[EP1 & BUP & INJ, i]
    a_TDP_1[EP1 & BUPC & INJ, EP1 & BUPC & INJ, i] <- m_TDP_1[EP1 & BUPC & INJ, i]
    a_TDP_2[EP1 & BUPC & INJ, EP1 & BUPC & INJ, i] <- m_TDP_2[EP1 & BUPC & INJ, i]
    a_TDP_3[EP1 & BUPC & INJ, EP1 & BUPC & INJ, i] <- m_TDP_3[EP1 & BUPC & INJ, i]
    a_TDP_1[EP1 & MET & INJ, EP1 & MET & INJ, i] <- m_TDP_1[EP1 & MET & INJ, i]
    a_TDP_2[EP1 & MET & INJ, EP1 & MET & INJ, i] <- m_TDP_2[EP1 & MET & INJ, i]
    a_TDP_3[EP1 & MET & INJ, EP1 & MET & INJ, i] <- m_TDP_3[EP1 & MET & INJ, i]
    a_TDP_1[EP1 & METC & INJ, EP1 & METC & INJ, i] <- m_TDP_1[EP1 & METC & INJ, i]
    a_TDP_2[EP1 & METC & INJ, EP1 & METC & INJ, i] <- m_TDP_2[EP1 & METC & INJ, i]
    a_TDP_3[EP1 & METC & INJ, EP1 & METC & INJ, i] <- m_TDP_3[EP1 & METC & INJ, i]
    a_TDP_1[EP1 & ABS & INJ, EP1 & ABS & INJ, i] <- m_TDP_1[EP1 & ABS & INJ, i]
    a_TDP_2[EP1 & ABS & INJ, EP1 & ABS & INJ, i] <- m_TDP_2[EP1 & ABS & INJ, i]
    a_TDP_3[EP1 & ABS & INJ, EP1 & ABS & INJ, i] <- m_TDP_3[EP1 & ABS & INJ, i]
    a_TDP_1[EP1 & REL & INJ, EP1 & REL & INJ, i] <- m_TDP_1[EP1 & REL & INJ, i]
    a_TDP_2[EP1 & REL & INJ, EP1 & REL & INJ, i] <- m_TDP_2[EP1 & REL & INJ, i]
    a_TDP_3[EP1 & REL & INJ, EP1 & REL & INJ, i] <- m_TDP_3[EP1 & REL & INJ, i]
    a_TDP_1[EP1 & ODN & INJ, EP1 & ODN & INJ, i] <- m_TDP_1[EP1 & ODN & INJ, i]
    a_TDP_2[EP1 & ODN & INJ, EP1 & ODN & INJ, i] <- m_TDP_2[EP1 & ODN & INJ, i]
    a_TDP_3[EP1 & ODN & INJ, EP1 & ODN & INJ, i] <- m_TDP_3[EP1 & ODN & INJ, i]
    a_TDP_1[EP1 & ODF & INJ, EP1 & ODF & INJ, i] <- m_TDP_1[EP1 & ODF & INJ, i]
    a_TDP_2[EP1 & ODF & INJ, EP1 & ODF & INJ, i] <- m_TDP_2[EP1 & ODF & INJ, i]
    a_TDP_3[EP1 & ODF & INJ, EP1 & ODF & INJ, i] <- m_TDP_3[EP1 & ODF & INJ, i]
    # Episode 2
    a_TDP_1[EP2 & BUP & INJ, EP2 & BUP & INJ, i] <- m_TDP_1[EP2 & BUP & INJ, i]
    a_TDP_2[EP2 & BUP & INJ, EP2 & BUP & INJ, i] <- m_TDP_2[EP2 & BUP & INJ, i]
    a_TDP_3[EP2 & BUP & INJ, EP2 & BUP & INJ, i] <- m_TDP_3[EP2 & BUP & INJ, i]
    a_TDP_1[EP2 & BUPC & INJ, EP2 & BUPC & INJ, i] <- m_TDP_1[EP2 & BUPC & INJ, i]
    a_TDP_2[EP2 & BUPC & INJ, EP2 & BUPC & INJ, i] <- m_TDP_2[EP2 & BUPC & INJ, i]
    a_TDP_3[EP2 & BUPC & INJ, EP2 & BUPC & INJ, i] <- m_TDP_3[EP2 & BUPC & INJ, i]
    a_TDP_1[EP2 & MET & INJ, EP2 & MET & INJ, i] <- m_TDP_1[EP2 & MET & INJ, i]
    a_TDP_2[EP2 & MET & INJ, EP2 & MET & INJ, i] <- m_TDP_2[EP2 & MET & INJ, i]
    a_TDP_3[EP2 & MET & INJ, EP2 & MET & INJ, i] <- m_TDP_3[EP2 & MET & INJ, i]
    a_TDP_1[EP2 & METC & INJ, EP2 & METC & INJ, i] <- m_TDP_1[EP2 & METC & INJ, i]
    a_TDP_2[EP2 & METC & INJ, EP2 & METC & INJ, i] <- m_TDP_2[EP2 & METC & INJ, i]
    a_TDP_3[EP2 & METC & INJ, EP2 & METC & INJ, i] <- m_TDP_3[EP2 & METC & INJ, i]
    a_TDP_1[EP2 & ABS & INJ, EP2 & ABS & INJ, i] <- m_TDP_1[EP2 & ABS & INJ, i]
    a_TDP_2[EP2 & ABS & INJ, EP2 & ABS & INJ, i] <- m_TDP_2[EP2 & ABS & INJ, i]
    a_TDP_3[EP2 & ABS & INJ, EP2 & ABS & INJ, i] <- m_TDP_3[EP2 & ABS & INJ, i]
    a_TDP_1[EP2 & REL & INJ, EP2 & REL & INJ, i] <- m_TDP_1[EP2 & REL & INJ, i]
    a_TDP_2[EP2 & REL & INJ, EP2 & REL & INJ, i] <- m_TDP_2[EP2 & REL & INJ, i]
    a_TDP_3[EP2 & REL & INJ, EP2 & REL & INJ, i] <- m_TDP_3[EP2 & REL & INJ, i]
    a_TDP_1[EP2 & ODN & INJ, EP2 & ODN & INJ, i] <- m_TDP_1[EP2 & ODN & INJ, i]
    a_TDP_2[EP2 & ODN & INJ, EP2 & ODN & INJ, i] <- m_TDP_2[EP2 & ODN & INJ, i]
    a_TDP_3[EP2 & ODN & INJ, EP2 & ODN & INJ, i] <- m_TDP_3[EP2 & ODN & INJ, i]
    a_TDP_1[EP2 & ODF & INJ, EP2 & ODF & INJ, i] <- m_TDP_1[EP2 & ODF & INJ, i]
    a_TDP_2[EP2 & ODF & INJ, EP2 & ODF & INJ, i] <- m_TDP_2[EP2 & ODF & INJ, i]
    a_TDP_3[EP2 & ODF & INJ, EP2 & ODF & INJ, i] <- m_TDP_3[EP2 & ODF & INJ, i]
    # Episode 3
    a_TDP_1[EP3 & BUP & INJ, EP3 & BUP & INJ, i] <- m_TDP_1[EP3 & BUP & INJ, i]
    a_TDP_2[EP3 & BUP & INJ, EP3 & BUP & INJ, i] <- m_TDP_2[EP3 & BUP & INJ, i]
    a_TDP_3[EP3 & BUP & INJ, EP3 & BUP & INJ, i] <- m_TDP_3[EP3 & BUP & INJ, i]
    a_TDP_1[EP3 & BUPC & INJ, EP3 & BUPC & INJ, i] <- m_TDP_1[EP3 & BUPC & INJ, i]
    a_TDP_2[EP3 & BUPC & INJ, EP3 & BUPC & INJ, i] <- m_TDP_2[EP3 & BUPC & INJ, i]
    a_TDP_3[EP3 & BUPC & INJ, EP3 & BUPC & INJ, i] <- m_TDP_3[EP3 & BUPC & INJ, i]
    a_TDP_1[EP3 & MET & INJ, EP3 & MET & INJ, i] <- m_TDP_1[EP3 & MET & INJ, i]
    a_TDP_2[EP3 & MET & INJ, EP3 & MET & INJ, i] <- m_TDP_2[EP3 & MET & INJ, i]
    a_TDP_3[EP3 & MET & INJ, EP3 & MET & INJ, i] <- m_TDP_3[EP3 & MET & INJ, i]
    a_TDP_1[EP3 & METC & INJ, EP3 & METC & INJ, i] <- m_TDP_1[EP3 & METC & INJ, i]
    a_TDP_2[EP3 & METC & INJ, EP3 & METC & INJ, i] <- m_TDP_2[EP3 & METC & INJ, i]
    a_TDP_3[EP3 & METC & INJ, EP3 & METC & INJ, i] <- m_TDP_3[EP3 & METC & INJ, i]
    a_TDP_1[EP3 & ABS & INJ, EP3 & ABS & INJ, i] <- m_TDP_1[EP3 & ABS & INJ, i]
    a_TDP_2[EP3 & ABS & INJ, EP3 & ABS & INJ, i] <- m_TDP_2[EP3 & ABS & INJ, i]
    a_TDP_3[EP3 & ABS & INJ, EP3 & ABS & INJ, i] <- m_TDP_3[EP3 & ABS & INJ, i]
    a_TDP_1[EP3 & REL & INJ, EP3 & REL & INJ, i] <- m_TDP_1[EP3 & REL & INJ, i]
    a_TDP_2[EP3 & REL & INJ, EP3 & REL & INJ, i] <- m_TDP_2[EP3 & REL & INJ, i]
    a_TDP_3[EP3 & REL & INJ, EP3 & REL & INJ, i] <- m_TDP_3[EP3 & REL & INJ, i]
    a_TDP_1[EP3 & ODN & INJ, EP3 & ODN & INJ, i] <- m_TDP_1[EP3 & ODN & INJ, i]
    a_TDP_2[EP3 & ODN & INJ, EP3 & ODN & INJ, i] <- m_TDP_2[EP3 & ODN & INJ, i]
    a_TDP_3[EP3 & ODN & INJ, EP3 & ODN & INJ, i] <- m_TDP_3[EP3 & ODN & INJ, i]
    a_TDP_1[EP3 & ODF & INJ, EP3 & ODF & INJ, i] <- m_TDP_1[EP3 & ODF & INJ, i]
    a_TDP_2[EP3 & ODF & INJ, EP3 & ODF & INJ, i] <- m_TDP_2[EP3 & ODF & INJ, i]
    a_TDP_3[EP3 & ODF & INJ, EP3 & ODF & INJ, i] <- m_TDP_3[EP3 & ODF & INJ, i]
  }
  
  #### Seroconversion ####
  # Apply seroconversion probability to re-weight NEG -> POS for to-states each time period
  # Probabilities applied equally across POS/NEG initially, re-weight by sero prob
  # Non-injection
  # BUP
  # From NEG
  a_TDP_1[NEG & NI, BUP & NI & NEG, ]  <- a_TDP_1[NEG & NI, BUP & NI & NEG, ] * (1 - p_HIV_BUP_NI - p_HCV_BUP_NI)
  a_TDP_2[NEG & NI, BUP & NI & NEG, ]  <- a_TDP_2[NEG & NI, BUP & NI & NEG, ] * (1 - p_HIV_BUP_NI - p_HCV_BUP_NI)
  a_TDP_3[NEG & NI, BUP & NI & NEG, ]  <- a_TDP_3[NEG & NI, BUP & NI & NEG, ] * (1 - p_HIV_BUP_NI - p_HCV_BUP_NI)
  a_TDP_1[NEG & NI, BUP & NI & HIV, ]  <- a_TDP_1[NEG & NI, BUP & NI & HIV, ] * p_HIV_BUP_NI
  a_TDP_2[NEG & NI, BUP & NI & HIV, ]  <- a_TDP_2[NEG & NI, BUP & NI & HIV, ] * p_HIV_BUP_NI
  a_TDP_3[NEG & NI, BUP & NI & HIV, ]  <- a_TDP_3[NEG & NI, BUP & NI & HIV, ] * p_HIV_BUP_NI
  a_TDP_1[NEG & NI, BUP & NI & HCV, ]  <- a_TDP_1[NEG & NI, BUP & NI & HCV, ] * p_HCV_BUP_NI
  a_TDP_2[NEG & NI, BUP & NI & HCV, ]  <- a_TDP_2[NEG & NI, BUP & NI & HCV, ] * p_HCV_BUP_NI
  a_TDP_3[NEG & NI, BUP & NI & HCV, ]  <- a_TDP_3[NEG & NI, BUP & NI & HCV, ] * p_HCV_BUP_NI
  # From HIV
  a_TDP_1[HIV & NI, BUP & NI & HIV, ]  <- a_TDP_1[HIV & NI, BUP & NI & HIV, ] * (1 - p_HIV_HCV_BUP_NI)
  a_TDP_2[HIV & NI, BUP & NI & HIV, ]  <- a_TDP_2[HIV & NI, BUP & NI & HIV, ] * (1 - p_HIV_HCV_BUP_NI)
  a_TDP_3[HIV & NI, BUP & NI & HIV, ]  <- a_TDP_3[HIV & NI, BUP & NI & HIV, ] * (1 - p_HIV_HCV_BUP_NI)
  a_TDP_1[HIV & NI, BUP & NI & COI, ]  <- a_TDP_1[HIV & NI, BUP & NI & COI, ] * p_HIV_HCV_BUP_NI # Probability of HCV conditional on HIV
  a_TDP_2[HIV & NI, BUP & NI & COI, ]  <- a_TDP_2[HIV & NI, BUP & NI & COI, ] * p_HIV_HCV_BUP_NI
  a_TDP_3[HIV & NI, BUP & NI & COI, ]  <- a_TDP_3[HIV & NI, BUP & NI & COI, ] * p_HIV_HCV_BUP_NI
  # From HCV
  a_TDP_1[HCV & NI, BUP & NI & HCV, ]  <- a_TDP_1[HCV & NI, BUP & NI & HCV, ] * (1 - p_HCV_HIV_BUP_NI)
  a_TDP_2[HCV & NI, BUP & NI & HCV, ]  <- a_TDP_2[HCV & NI, BUP & NI & HCV, ] * (1 - p_HCV_HIV_BUP_NI)
  a_TDP_3[HCV & NI, BUP & NI & HCV, ]  <- a_TDP_3[HCV & NI, BUP & NI & HCV, ] * (1 - p_HCV_HIV_BUP_NI)
  a_TDP_1[HCV & NI, BUP & NI & COI, ]  <- a_TDP_1[HCV & NI, BUP & NI & COI, ] * p_HCV_HIV_BUP_NI # Probability of HIV conditional on HCV
  a_TDP_2[HCV & NI, BUP & NI & COI, ]  <- a_TDP_2[HCV & NI, BUP & NI & COI, ] * p_HCV_HIV_BUP_NI
  a_TDP_3[HCV & NI, BUP & NI & COI, ]  <- a_TDP_3[HCV & NI, BUP & NI & COI, ] * p_HCV_HIV_BUP_NI
  
  # BUPC
  # From NEG
  a_TDP_1[NEG & NI, BUPC & NI & NEG, ]  <- a_TDP_1[NEG & NI, BUPC & NI & NEG, ] * (1 - p_HIV_BUPC_NI - p_HCV_BUPC_NI)
  a_TDP_2[NEG & NI, BUPC & NI & NEG, ]  <- a_TDP_2[NEG & NI, BUPC & NI & NEG, ] * (1 - p_HIV_BUPC_NI - p_HCV_BUPC_NI)
  a_TDP_3[NEG & NI, BUPC & NI & NEG, ]  <- a_TDP_3[NEG & NI, BUPC & NI & NEG, ] * (1 - p_HIV_BUPC_NI - p_HCV_BUPC_NI)
  a_TDP_1[NEG & NI, BUPC & NI & HIV, ]  <- a_TDP_1[NEG & NI, BUPC & NI & HIV, ] * p_HIV_BUPC_NI
  a_TDP_2[NEG & NI, BUPC & NI & HIV, ]  <- a_TDP_2[NEG & NI, BUPC & NI & HIV, ] * p_HIV_BUPC_NI
  a_TDP_3[NEG & NI, BUPC & NI & HIV, ]  <- a_TDP_3[NEG & NI, BUPC & NI & HIV, ] * p_HIV_BUPC_NI
  a_TDP_1[NEG & NI, BUPC & NI & HCV, ]  <- a_TDP_1[NEG & NI, BUPC & NI & HCV, ] * p_HCV_BUPC_NI
  a_TDP_2[NEG & NI, BUPC & NI & HCV, ]  <- a_TDP_2[NEG & NI, BUPC & NI & HCV, ] * p_HCV_BUPC_NI
  a_TDP_3[NEG & NI, BUPC & NI & HCV, ]  <- a_TDP_3[NEG & NI, BUPC & NI & HCV, ] * p_HCV_BUPC_NI
  # From HIV
  a_TDP_1[HIV & NI, BUPC & NI & HIV, ]  <- a_TDP_1[HIV & NI, BUPC & NI & HIV, ] * (1 - p_HIV_HCV_BUPC_NI)
  a_TDP_2[HIV & NI, BUPC & NI & HIV, ]  <- a_TDP_2[HIV & NI, BUPC & NI & HIV, ] * (1 - p_HIV_HCV_BUPC_NI)
  a_TDP_3[HIV & NI, BUPC & NI & HIV, ]  <- a_TDP_3[HIV & NI, BUPC & NI & HIV, ] * (1 - p_HIV_HCV_BUPC_NI)
  a_TDP_1[HIV & NI, BUPC & NI & COI, ]  <- a_TDP_1[HIV & NI, BUPC & NI & COI, ] * p_HIV_HCV_BUPC_NI # Probability of HCV conditional on HIV
  a_TDP_2[HIV & NI, BUPC & NI & COI, ]  <- a_TDP_2[HIV & NI, BUPC & NI & COI, ] * p_HIV_HCV_BUPC_NI
  a_TDP_3[HIV & NI, BUPC & NI & COI, ]  <- a_TDP_3[HIV & NI, BUPC & NI & COI, ] * p_HIV_HCV_BUPC_NI
  # From HCV
  a_TDP_1[HCV & NI, BUPC & NI & HCV, ]  <- a_TDP_1[HCV & NI, BUPC & NI & HCV, ] * (1 - p_HCV_HIV_BUPC_NI)
  a_TDP_2[HCV & NI, BUPC & NI & HCV, ]  <- a_TDP_2[HCV & NI, BUPC & NI & HCV, ] * (1 - p_HCV_HIV_BUPC_NI)
  a_TDP_3[HCV & NI, BUPC & NI & HCV, ]  <- a_TDP_3[HCV & NI, BUPC & NI & HCV, ] * (1 - p_HCV_HIV_BUPC_NI)
  a_TDP_1[HCV & NI, BUPC & NI & COI, ]  <- a_TDP_1[HCV & NI, BUPC & NI & COI, ] * p_HCV_HIV_BUPC_NI # Probability of HIV conditional on HCV
  a_TDP_2[HCV & NI, BUPC & NI & COI, ]  <- a_TDP_2[HCV & NI, BUPC & NI & COI, ] * p_HCV_HIV_BUPC_NI
  a_TDP_3[HCV & NI, BUPC & NI & COI, ]  <- a_TDP_3[HCV & NI, BUPC & NI & COI, ] * p_HCV_HIV_BUPC_NI
  
  # MET
  # From NEG
  a_TDP_1[NEG & NI, MET & NI & NEG, ]  <- a_TDP_1[NEG & NI, MET & NI & NEG, ] * (1 - p_HIV_MET_NI - p_HCV_MET_NI)
  a_TDP_2[NEG & NI, MET & NI & NEG, ]  <- a_TDP_2[NEG & NI, MET & NI & NEG, ] * (1 - p_HIV_MET_NI - p_HCV_MET_NI)
  a_TDP_3[NEG & NI, MET & NI & NEG, ]  <- a_TDP_3[NEG & NI, MET & NI & NEG, ] * (1 - p_HIV_MET_NI - p_HCV_MET_NI)
  a_TDP_1[NEG & NI, MET & NI & HIV, ]  <- a_TDP_1[NEG & NI, MET & NI & HIV, ] * p_HIV_MET_NI
  a_TDP_2[NEG & NI, MET & NI & HIV, ]  <- a_TDP_2[NEG & NI, MET & NI & HIV, ] * p_HIV_MET_NI
  a_TDP_3[NEG & NI, MET & NI & HIV, ]  <- a_TDP_3[NEG & NI, MET & NI & HIV, ] * p_HIV_MET_NI
  a_TDP_1[NEG & NI, MET & NI & HCV, ]  <- a_TDP_1[NEG & NI, MET & NI & HCV, ] * p_HCV_MET_NI
  a_TDP_2[NEG & NI, MET & NI & HCV, ]  <- a_TDP_2[NEG & NI, MET & NI & HCV, ] * p_HCV_MET_NI
  a_TDP_3[NEG & NI, MET & NI & HCV, ]  <- a_TDP_3[NEG & NI, MET & NI & HCV, ] * p_HCV_MET_NI
  # From HIV
  a_TDP_1[HIV & NI, MET & NI & HIV, ]  <- a_TDP_1[HIV & NI, MET & NI & HIV, ] * (1 - p_HIV_HCV_MET_NI)
  a_TDP_2[HIV & NI, MET & NI & HIV, ]  <- a_TDP_2[HIV & NI, MET & NI & HIV, ] * (1 - p_HIV_HCV_MET_NI)
  a_TDP_3[HIV & NI, MET & NI & HIV, ]  <- a_TDP_3[HIV & NI, MET & NI & HIV, ] * (1 - p_HIV_HCV_MET_NI)
  a_TDP_1[HIV & NI, MET & NI & COI, ]  <- a_TDP_1[HIV & NI, MET & NI & COI, ] * p_HIV_HCV_MET_NI # Probability of HCV conditional on HIV
  a_TDP_2[HIV & NI, MET & NI & COI, ]  <- a_TDP_2[HIV & NI, MET & NI & COI, ] * p_HIV_HCV_MET_NI
  a_TDP_3[HIV & NI, MET & NI & COI, ]  <- a_TDP_3[HIV & NI, MET & NI & COI, ] * p_HIV_HCV_MET_NI
  # From HCV
  a_TDP_1[HCV & NI, MET & NI & HCV, ]  <- a_TDP_1[HCV & NI, MET & NI & HCV, ] * (1 - p_HCV_HIV_MET_NI)
  a_TDP_2[HCV & NI, MET & NI & HCV, ]  <- a_TDP_2[HCV & NI, MET & NI & HCV, ] * (1 - p_HCV_HIV_MET_NI)
  a_TDP_3[HCV & NI, MET & NI & HCV, ]  <- a_TDP_3[HCV & NI, MET & NI & HCV, ] * (1 - p_HCV_HIV_MET_NI)
  a_TDP_1[HCV & NI, MET & NI & COI, ]  <- a_TDP_1[HCV & NI, MET & NI & COI, ] * p_HCV_HIV_MET_NI # Probability of HIV conditional on HCV
  a_TDP_2[HCV & NI, MET & NI & COI, ]  <- a_TDP_2[HCV & NI, MET & NI & COI, ] * p_HCV_HIV_MET_NI
  a_TDP_3[HCV & NI, MET & NI & COI, ]  <- a_TDP_3[HCV & NI, MET & NI & COI, ] * p_HCV_HIV_MET_NI
  
  # METC
  # From NEG
  a_TDP_1[NEG & NI, METC & NI & NEG, ]  <- a_TDP_1[NEG & NI, METC & NI & NEG, ] * (1 - p_HIV_METC_NI - p_HCV_METC_NI)
  a_TDP_2[NEG & NI, METC & NI & NEG, ]  <- a_TDP_2[NEG & NI, METC & NI & NEG, ] * (1 - p_HIV_METC_NI - p_HCV_METC_NI)
  a_TDP_3[NEG & NI, METC & NI & NEG, ]  <- a_TDP_3[NEG & NI, METC & NI & NEG, ] * (1 - p_HIV_METC_NI - p_HCV_METC_NI)
  a_TDP_1[NEG & NI, METC & NI & HIV, ]  <- a_TDP_1[NEG & NI, METC & NI & HIV, ] * p_HIV_METC_NI
  a_TDP_2[NEG & NI, METC & NI & HIV, ]  <- a_TDP_2[NEG & NI, METC & NI & HIV, ] * p_HIV_METC_NI
  a_TDP_3[NEG & NI, METC & NI & HIV, ]  <- a_TDP_3[NEG & NI, METC & NI & HIV, ] * p_HIV_METC_NI
  a_TDP_1[NEG & NI, METC & NI & HCV, ]  <- a_TDP_1[NEG & NI, METC & NI & HCV, ] * p_HCV_METC_NI
  a_TDP_2[NEG & NI, METC & NI & HCV, ]  <- a_TDP_2[NEG & NI, METC & NI & HCV, ] * p_HCV_METC_NI
  a_TDP_3[NEG & NI, METC & NI & HCV, ]  <- a_TDP_3[NEG & NI, METC & NI & HCV, ] * p_HCV_METC_NI
  # From HIV
  a_TDP_1[HIV & NI, METC & NI & HIV, ]  <- a_TDP_1[HIV & NI, METC & NI & HIV, ] * (1 - p_HIV_HCV_METC_NI)
  a_TDP_2[HIV & NI, METC & NI & HIV, ]  <- a_TDP_2[HIV & NI, METC & NI & HIV, ] * (1 - p_HIV_HCV_METC_NI)
  a_TDP_3[HIV & NI, METC & NI & HIV, ]  <- a_TDP_3[HIV & NI, METC & NI & HIV, ] * (1 - p_HIV_HCV_METC_NI)
  a_TDP_1[HIV & NI, METC & NI & COI, ]  <- a_TDP_1[HIV & NI, METC & NI & COI, ] * p_HIV_HCV_METC_NI # Probability of HCV conditional on HIV
  a_TDP_2[HIV & NI, METC & NI & COI, ]  <- a_TDP_2[HIV & NI, METC & NI & COI, ] * p_HIV_HCV_METC_NI
  a_TDP_3[HIV & NI, METC & NI & COI, ]  <- a_TDP_3[HIV & NI, METC & NI & COI, ] * p_HIV_HCV_METC_NI
  # From HCV
  a_TDP_1[HCV & NI, METC & NI & HCV, ]  <- a_TDP_1[HCV & NI, METC & NI & HCV, ] * (1 - p_HCV_HIV_METC_NI)
  a_TDP_2[HCV & NI, METC & NI & HCV, ]  <- a_TDP_2[HCV & NI, METC & NI & HCV, ] * (1 - p_HCV_HIV_METC_NI)
  a_TDP_3[HCV & NI, METC & NI & HCV, ]  <- a_TDP_3[HCV & NI, METC & NI & HCV, ] * (1 - p_HCV_HIV_METC_NI)
  a_TDP_1[HCV & NI, METC & NI & COI, ]  <- a_TDP_1[HCV & NI, METC & NI & COI, ] * p_HCV_HIV_METC_NI # Probability of HIV conditional on HCV
  a_TDP_2[HCV & NI, METC & NI & COI, ]  <- a_TDP_2[HCV & NI, METC & NI & COI, ] * p_HCV_HIV_METC_NI
  a_TDP_3[HCV & NI, METC & NI & COI, ]  <- a_TDP_3[HCV & NI, METC & NI & COI, ] * p_HCV_HIV_METC_NI
  
  # REL
  # From NEG
  a_TDP_1[NEG & NI, REL & NI & NEG, ]  <- a_TDP_1[NEG & NI, REL & NI & NEG, ] * (1 - p_HIV_REL_NI - p_HCV_REL_NI)
  a_TDP_2[NEG & NI, REL & NI & NEG, ]  <- a_TDP_2[NEG & NI, REL & NI & NEG, ] * (1 - p_HIV_REL_NI - p_HCV_REL_NI)
  a_TDP_3[NEG & NI, REL & NI & NEG, ]  <- a_TDP_3[NEG & NI, REL & NI & NEG, ] * (1 - p_HIV_REL_NI - p_HCV_REL_NI)
  a_TDP_1[NEG & NI, REL & NI & HIV, ]  <- a_TDP_1[NEG & NI, REL & NI & HIV, ] * p_HIV_REL_NI
  a_TDP_2[NEG & NI, REL & NI & HIV, ]  <- a_TDP_2[NEG & NI, REL & NI & HIV, ] * p_HIV_REL_NI
  a_TDP_3[NEG & NI, REL & NI & HIV, ]  <- a_TDP_3[NEG & NI, REL & NI & HIV, ] * p_HIV_REL_NI
  a_TDP_1[NEG & NI, REL & NI & HCV, ]  <- a_TDP_1[NEG & NI, REL & NI & HCV, ] * p_HCV_REL_NI
  a_TDP_2[NEG & NI, REL & NI & HCV, ]  <- a_TDP_2[NEG & NI, REL & NI & HCV, ] * p_HCV_REL_NI
  a_TDP_3[NEG & NI, REL & NI & HCV, ]  <- a_TDP_3[NEG & NI, REL & NI & HCV, ] * p_HCV_REL_NI
  # From HIV
  a_TDP_1[HIV & NI, REL & NI & HIV, ]  <- a_TDP_1[HIV & NI, REL & NI & HIV, ] * (1 - p_HIV_HCV_REL_NI)
  a_TDP_2[HIV & NI, REL & NI & HIV, ]  <- a_TDP_2[HIV & NI, REL & NI & HIV, ] * (1 - p_HIV_HCV_REL_NI)
  a_TDP_3[HIV & NI, REL & NI & HIV, ]  <- a_TDP_3[HIV & NI, REL & NI & HIV, ] * (1 - p_HIV_HCV_REL_NI)
  a_TDP_1[HIV & NI, REL & NI & COI, ]  <- a_TDP_1[HIV & NI, REL & NI & COI, ] * p_HIV_HCV_REL_NI # Probability of HCV conditional on HIV
  a_TDP_2[HIV & NI, REL & NI & COI, ]  <- a_TDP_2[HIV & NI, REL & NI & COI, ] * p_HIV_HCV_REL_NI
  a_TDP_3[HIV & NI, REL & NI & COI, ]  <- a_TDP_3[HIV & NI, REL & NI & COI, ] * p_HIV_HCV_REL_NI
  # From HCV
  a_TDP_1[HCV & NI, REL & NI & HCV, ]  <- a_TDP_1[HCV & NI, REL & NI & HCV, ] * (1 - p_HCV_HIV_REL_NI)
  a_TDP_2[HCV & NI, REL & NI & HCV, ]  <- a_TDP_2[HCV & NI, REL & NI & HCV, ] * (1 - p_HCV_HIV_REL_NI)
  a_TDP_3[HCV & NI, REL & NI & HCV, ]  <- a_TDP_3[HCV & NI, REL & NI & HCV, ] * (1 - p_HCV_HIV_REL_NI)
  a_TDP_1[HCV & NI, REL & NI & COI, ]  <- a_TDP_1[HCV & NI, REL & NI & COI, ] * p_HCV_HIV_REL_NI # Probability of HIV conditional on HCV
  a_TDP_2[HCV & NI, REL & NI & COI, ]  <- a_TDP_2[HCV & NI, REL & NI & COI, ] * p_HCV_HIV_REL_NI
  a_TDP_3[HCV & NI, REL & NI & COI, ]  <- a_TDP_3[HCV & NI, REL & NI & COI, ] * p_HCV_HIV_REL_NI
  
  # ODN
  # From NEG
  a_TDP_1[NEG & NI, ODN & NI & NEG, ]  <- a_TDP_1[NEG & NI, ODN & NI & NEG, ] * (1 - p_HIV_ODN_NI - p_HCV_ODN_NI)
  a_TDP_2[NEG & NI, ODN & NI & NEG, ]  <- a_TDP_2[NEG & NI, ODN & NI & NEG, ] * (1 - p_HIV_ODN_NI - p_HCV_ODN_NI)
  a_TDP_3[NEG & NI, ODN & NI & NEG, ]  <- a_TDP_3[NEG & NI, ODN & NI & NEG, ] * (1 - p_HIV_ODN_NI - p_HCV_ODN_NI)
  a_TDP_1[NEG & NI, ODN & NI & HIV, ]  <- a_TDP_1[NEG & NI, ODN & NI & HIV, ] * p_HIV_ODN_NI
  a_TDP_2[NEG & NI, ODN & NI & HIV, ]  <- a_TDP_2[NEG & NI, ODN & NI & HIV, ] * p_HIV_ODN_NI
  a_TDP_3[NEG & NI, ODN & NI & HIV, ]  <- a_TDP_3[NEG & NI, ODN & NI & HIV, ] * p_HIV_ODN_NI
  a_TDP_1[NEG & NI, ODN & NI & HCV, ]  <- a_TDP_1[NEG & NI, ODN & NI & HCV, ] * p_HCV_ODN_NI
  a_TDP_2[NEG & NI, ODN & NI & HCV, ]  <- a_TDP_2[NEG & NI, ODN & NI & HCV, ] * p_HCV_ODN_NI
  a_TDP_3[NEG & NI, ODN & NI & HCV, ]  <- a_TDP_3[NEG & NI, ODN & NI & HCV, ] * p_HCV_ODN_NI
  # From HIV
  a_TDP_1[HIV & NI, ODN & NI & HIV, ]  <- a_TDP_1[HIV & NI, ODN & NI & HIV, ] * (1 - p_HIV_HCV_ODN_NI)
  a_TDP_2[HIV & NI, ODN & NI & HIV, ]  <- a_TDP_2[HIV & NI, ODN & NI & HIV, ] * (1 - p_HIV_HCV_ODN_NI)
  a_TDP_3[HIV & NI, ODN & NI & HIV, ]  <- a_TDP_3[HIV & NI, ODN & NI & HIV, ] * (1 - p_HIV_HCV_ODN_NI)
  a_TDP_1[HIV & NI, ODN & NI & COI, ]  <- a_TDP_1[HIV & NI, ODN & NI & COI, ] * p_HIV_HCV_ODN_NI # Probability of HCV conditional on HIV
  a_TDP_2[HIV & NI, ODN & NI & COI, ]  <- a_TDP_2[HIV & NI, ODN & NI & COI, ] * p_HIV_HCV_ODN_NI
  a_TDP_3[HIV & NI, ODN & NI & COI, ]  <- a_TDP_3[HIV & NI, ODN & NI & COI, ] * p_HIV_HCV_ODN_NI
  # From HCV
  a_TDP_1[HCV & NI, ODN & NI & HCV, ]  <- a_TDP_1[HCV & NI, ODN & NI & HCV, ] * (1 - p_HCV_HIV_ODN_NI)
  a_TDP_2[HCV & NI, ODN & NI & HCV, ]  <- a_TDP_2[HCV & NI, ODN & NI & HCV, ] * (1 - p_HCV_HIV_ODN_NI)
  a_TDP_3[HCV & NI, ODN & NI & HCV, ]  <- a_TDP_3[HCV & NI, ODN & NI & HCV, ] * (1 - p_HCV_HIV_ODN_NI)
  a_TDP_1[HCV & NI, ODN & NI & COI, ]  <- a_TDP_1[HCV & NI, ODN & NI & COI, ] * p_HCV_HIV_ODN_NI # Probability of HIV conditional on HCV
  a_TDP_2[HCV & NI, ODN & NI & COI, ]  <- a_TDP_2[HCV & NI, ODN & NI & COI, ] * p_HCV_HIV_ODN_NI
  a_TDP_3[HCV & NI, ODN & NI & COI, ]  <- a_TDP_3[HCV & NI, ODN & NI & COI, ] * p_HCV_HIV_ODN_NI

  # ODF
  # From NEG
  a_TDP_1[NEG & NI, ODF & NI & NEG, ]  <- a_TDP_1[NEG & NI, ODF & NI & NEG, ] * 1
  a_TDP_2[NEG & NI, ODF & NI & NEG, ]  <- a_TDP_2[NEG & NI, ODF & NI & NEG, ] * 1
  a_TDP_3[NEG & NI, ODF & NI & NEG, ]  <- a_TDP_3[NEG & NI, ODF & NI & NEG, ] * 1
  a_TDP_1[NEG & NI, ODF & NI & HIV, ]  <- a_TDP_1[NEG & NI, ODF & NI & HIV, ] * 0
  a_TDP_2[NEG & NI, ODF & NI & HIV, ]  <- a_TDP_2[NEG & NI, ODF & NI & HIV, ] * 0
  a_TDP_3[NEG & NI, ODF & NI & HIV, ]  <- a_TDP_3[NEG & NI, ODF & NI & HIV, ] * 0
  a_TDP_1[NEG & NI, ODF & NI & HCV, ]  <- a_TDP_1[NEG & NI, ODF & NI & HCV, ] * 0
  a_TDP_2[NEG & NI, ODF & NI & HCV, ]  <- a_TDP_2[NEG & NI, ODF & NI & HCV, ] * 0
  a_TDP_3[NEG & NI, ODF & NI & HCV, ]  <- a_TDP_3[NEG & NI, ODF & NI & HCV, ] * 0
  # From HIV
  a_TDP_1[HIV & NI, ODF & NI & HIV, ]  <- a_TDP_1[HIV & NI, ODF & NI & HIV, ] * 1
  a_TDP_2[HIV & NI, ODF & NI & HIV, ]  <- a_TDP_2[HIV & NI, ODF & NI & HIV, ] * 1
  a_TDP_3[HIV & NI, ODF & NI & HIV, ]  <- a_TDP_3[HIV & NI, ODF & NI & HIV, ] * 1
  a_TDP_1[HIV & NI, ODF & NI & COI, ]  <- a_TDP_1[HIV & NI, ODF & NI & COI, ] * 0
  a_TDP_2[HIV & NI, ODF & NI & COI, ]  <- a_TDP_2[HIV & NI, ODF & NI & COI, ] * 0
  a_TDP_3[HIV & NI, ODF & NI & COI, ]  <- a_TDP_3[HIV & NI, ODF & NI & COI, ] * 0
  # From HCV
  a_TDP_1[HCV & NI, ODF & NI & HCV, ]  <- a_TDP_1[HCV & NI, ODF & NI & HCV, ] * 1
  a_TDP_2[HCV & NI, ODF & NI & HCV, ]  <- a_TDP_2[HCV & NI, ODF & NI & HCV, ] * 1
  a_TDP_3[HCV & NI, ODF & NI & HCV, ]  <- a_TDP_3[HCV & NI, ODF & NI & HCV, ] * 1
  a_TDP_1[HCV & NI, ODF & NI & COI, ]  <- a_TDP_1[HCV & NI, ODF & NI & COI, ] * 0
  a_TDP_2[HCV & NI, ODF & NI & COI, ]  <- a_TDP_2[HCV & NI, ODF & NI & COI, ] * 0
  a_TDP_3[HCV & NI, ODF & NI & COI, ]  <- a_TDP_3[HCV & NI, ODF & NI & COI, ] * 0
  
  # ABS
  # From NEG
  a_TDP_1[NEG & NI, ABS & NI & NEG, ]  <- a_TDP_1[NEG & NI, ABS & NI & NEG, ] * (1 - p_HIV_ABS_NI - p_HCV_ABS_NI)
  a_TDP_2[NEG & NI, ABS & NI & NEG, ]  <- a_TDP_2[NEG & NI, ABS & NI & NEG, ] * (1 - p_HIV_ABS_NI - p_HCV_ABS_NI)
  a_TDP_3[NEG & NI, ABS & NI & NEG, ]  <- a_TDP_3[NEG & NI, ABS & NI & NEG, ] * (1 - p_HIV_ABS_NI - p_HCV_ABS_NI)
  a_TDP_1[NEG & NI, ABS & NI & HIV, ]  <- a_TDP_1[NEG & NI, ABS & NI & HIV, ] * p_HIV_ABS_NI
  a_TDP_2[NEG & NI, ABS & NI & HIV, ]  <- a_TDP_2[NEG & NI, ABS & NI & HIV, ] * p_HIV_ABS_NI
  a_TDP_3[NEG & NI, ABS & NI & HIV, ]  <- a_TDP_3[NEG & NI, ABS & NI & HIV, ] * p_HIV_ABS_NI
  a_TDP_1[NEG & NI, ABS & NI & HCV, ]  <- a_TDP_1[NEG & NI, ABS & NI & HCV, ] * p_HCV_ABS_NI
  a_TDP_2[NEG & NI, ABS & NI & HCV, ]  <- a_TDP_2[NEG & NI, ABS & NI & HCV, ] * p_HCV_ABS_NI
  a_TDP_3[NEG & NI, ABS & NI & HCV, ]  <- a_TDP_3[NEG & NI, ABS & NI & HCV, ] * p_HCV_ABS_NI
  # From HIV
  a_TDP_1[HIV & NI, ABS & NI & HIV, ]  <- a_TDP_1[HIV & NI, ABS & NI & HIV, ] * (1 - p_HIV_HCV_ABS_NI)
  a_TDP_2[HIV & NI, ABS & NI & HIV, ]  <- a_TDP_2[HIV & NI, ABS & NI & HIV, ] * (1 - p_HIV_HCV_ABS_NI)
  a_TDP_3[HIV & NI, ABS & NI & HIV, ]  <- a_TDP_3[HIV & NI, ABS & NI & HIV, ] * (1 - p_HIV_HCV_ABS_NI)
  a_TDP_1[HIV & NI, ABS & NI & COI, ]  <- a_TDP_1[HIV & NI, ABS & NI & COI, ] * p_HIV_HCV_ABS_NI # Probability of HCV conditional on HIV
  a_TDP_2[HIV & NI, ABS & NI & COI, ]  <- a_TDP_2[HIV & NI, ABS & NI & COI, ] * p_HIV_HCV_ABS_NI
  a_TDP_3[HIV & NI, ABS & NI & COI, ]  <- a_TDP_3[HIV & NI, ABS & NI & COI, ] * p_HIV_HCV_ABS_NI
  # From HCV
  a_TDP_1[HCV & NI, ABS & NI & HCV, ]  <- a_TDP_1[HCV & NI, ABS & NI & HCV, ] * (1 - p_HCV_HIV_ABS_NI)
  a_TDP_2[HCV & NI, ABS & NI & HCV, ]  <- a_TDP_2[HCV & NI, ABS & NI & HCV, ] * (1 - p_HCV_HIV_ABS_NI)
  a_TDP_3[HCV & NI, ABS & NI & HCV, ]  <- a_TDP_3[HCV & NI, ABS & NI & HCV, ] * (1 - p_HCV_HIV_ABS_NI)
  a_TDP_1[HCV & NI, ABS & NI & COI, ]  <- a_TDP_1[HCV & NI, ABS & NI & COI, ] * p_HCV_HIV_ABS_NI # Probability of HIV conditional on HCV
  a_TDP_2[HCV & NI, ABS & NI & COI, ]  <- a_TDP_2[HCV & NI, ABS & NI & COI, ] * p_HCV_HIV_ABS_NI
  a_TDP_3[HCV & NI, ABS & NI & COI, ]  <- a_TDP_3[HCV & NI, ABS & NI & COI, ] * p_HCV_HIV_ABS_NI
  
  # Injection
  # BUP
  # From NEG
  a_TDP_1[NEG & INJ, BUP & INJ & NEG, ]  <- a_TDP_1[NEG & INJ, BUP & INJ & NEG, ] * (1 - p_HIV_BUP_INJ - p_HCV_BUP_INJ)
  a_TDP_2[NEG & INJ, BUP & INJ & NEG, ]  <- a_TDP_2[NEG & INJ, BUP & INJ & NEG, ] * (1 - p_HIV_BUP_INJ - p_HCV_BUP_INJ)
  a_TDP_3[NEG & INJ, BUP & INJ & NEG, ]  <- a_TDP_3[NEG & INJ, BUP & INJ & NEG, ] * (1 - p_HIV_BUP_INJ - p_HCV_BUP_INJ)
  a_TDP_1[NEG & INJ, BUP & INJ & HIV, ]  <- a_TDP_1[NEG & INJ, BUP & INJ & HIV, ] * p_HIV_BUP_INJ
  a_TDP_2[NEG & INJ, BUP & INJ & HIV, ]  <- a_TDP_2[NEG & INJ, BUP & INJ & HIV, ] * p_HIV_BUP_INJ
  a_TDP_3[NEG & INJ, BUP & INJ & HIV, ]  <- a_TDP_3[NEG & INJ, BUP & INJ & HIV, ] * p_HIV_BUP_INJ
  a_TDP_1[NEG & INJ, BUP & INJ & HCV, ]  <- a_TDP_1[NEG & INJ, BUP & INJ & HCV, ] * p_HCV_BUP_INJ
  a_TDP_2[NEG & INJ, BUP & INJ & HCV, ]  <- a_TDP_2[NEG & INJ, BUP & INJ & HCV, ] * p_HCV_BUP_INJ
  a_TDP_3[NEG & INJ, BUP & INJ & HCV, ]  <- a_TDP_3[NEG & INJ, BUP & INJ & HCV, ] * p_HCV_BUP_INJ
  # From HIV
  a_TDP_1[HIV & INJ, BUP & INJ & HIV, ]  <- a_TDP_1[HIV & INJ, BUP & INJ & HIV, ] * (1 - p_HIV_HCV_BUP_INJ)
  a_TDP_2[HIV & INJ, BUP & INJ & HIV, ]  <- a_TDP_2[HIV & INJ, BUP & INJ & HIV, ] * (1 - p_HIV_HCV_BUP_INJ)
  a_TDP_3[HIV & INJ, BUP & INJ & HIV, ]  <- a_TDP_3[HIV & INJ, BUP & INJ & HIV, ] * (1 - p_HIV_HCV_BUP_INJ)
  a_TDP_1[HIV & INJ, BUP & INJ & COI, ]  <- a_TDP_1[HIV & INJ, BUP & INJ & COI, ] * p_HIV_HCV_BUP_INJ # Probability of HCV conditional on HIV
  a_TDP_2[HIV & INJ, BUP & INJ & COI, ]  <- a_TDP_2[HIV & INJ, BUP & INJ & COI, ] * p_HIV_HCV_BUP_INJ
  a_TDP_3[HIV & INJ, BUP & INJ & COI, ]  <- a_TDP_3[HIV & INJ, BUP & INJ & COI, ] * p_HIV_HCV_BUP_INJ
  # From HCV
  a_TDP_1[HCV & INJ, BUP & INJ & HCV, ]  <- a_TDP_1[HCV & INJ, BUP & INJ & HCV, ] * (1 - p_HCV_HIV_BUP_INJ)
  a_TDP_2[HCV & INJ, BUP & INJ & HCV, ]  <- a_TDP_2[HCV & INJ, BUP & INJ & HCV, ] * (1 - p_HCV_HIV_BUP_INJ)
  a_TDP_3[HCV & INJ, BUP & INJ & HCV, ]  <- a_TDP_3[HCV & INJ, BUP & INJ & HCV, ] * (1 - p_HCV_HIV_BUP_INJ)
  a_TDP_1[HCV & INJ, BUP & INJ & COI, ]  <- a_TDP_1[HCV & INJ, BUP & INJ & COI, ] * p_HCV_HIV_BUP_INJ # Probability of HIV conditional on HCV
  a_TDP_2[HCV & INJ, BUP & INJ & COI, ]  <- a_TDP_2[HCV & INJ, BUP & INJ & COI, ] * p_HCV_HIV_BUP_INJ
  a_TDP_3[HCV & INJ, BUP & INJ & COI, ]  <- a_TDP_3[HCV & INJ, BUP & INJ & COI, ] * p_HCV_HIV_BUP_INJ
  
  # BUPC
  # From NEG
  a_TDP_1[NEG & INJ, BUPC & INJ & NEG, ]  <- a_TDP_1[NEG & INJ, BUPC & INJ & NEG, ] * (1 - p_HIV_BUPC_INJ - p_HCV_BUPC_INJ)
  a_TDP_2[NEG & INJ, BUPC & INJ & NEG, ]  <- a_TDP_2[NEG & INJ, BUPC & INJ & NEG, ] * (1 - p_HIV_BUPC_INJ - p_HCV_BUPC_INJ)
  a_TDP_3[NEG & INJ, BUPC & INJ & NEG, ]  <- a_TDP_3[NEG & INJ, BUPC & INJ & NEG, ] * (1 - p_HIV_BUPC_INJ - p_HCV_BUPC_INJ)
  a_TDP_1[NEG & INJ, BUPC & INJ & HIV, ]  <- a_TDP_1[NEG & INJ, BUPC & INJ & HIV, ] * p_HIV_BUPC_INJ
  a_TDP_2[NEG & INJ, BUPC & INJ & HIV, ]  <- a_TDP_2[NEG & INJ, BUPC & INJ & HIV, ] * p_HIV_BUPC_INJ
  a_TDP_3[NEG & INJ, BUPC & INJ & HIV, ]  <- a_TDP_3[NEG & INJ, BUPC & INJ & HIV, ] * p_HIV_BUPC_INJ
  a_TDP_1[NEG & INJ, BUPC & INJ & HCV, ]  <- a_TDP_1[NEG & INJ, BUPC & INJ & HCV, ] * p_HCV_BUPC_INJ
  a_TDP_2[NEG & INJ, BUPC & INJ & HCV, ]  <- a_TDP_2[NEG & INJ, BUPC & INJ & HCV, ] * p_HCV_BUPC_INJ
  a_TDP_3[NEG & INJ, BUPC & INJ & HCV, ]  <- a_TDP_3[NEG & INJ, BUPC & INJ & HCV, ] * p_HCV_BUPC_INJ
  # From HIV
  a_TDP_1[HIV & INJ, BUPC & INJ & HIV, ]  <- a_TDP_1[HIV & INJ, BUPC & INJ & HIV, ] * (1 - p_HIV_HCV_BUPC_INJ)
  a_TDP_2[HIV & INJ, BUPC & INJ & HIV, ]  <- a_TDP_2[HIV & INJ, BUPC & INJ & HIV, ] * (1 - p_HIV_HCV_BUPC_INJ)
  a_TDP_3[HIV & INJ, BUPC & INJ & HIV, ]  <- a_TDP_3[HIV & INJ, BUPC & INJ & HIV, ] * (1 - p_HIV_HCV_BUPC_INJ)
  a_TDP_1[HIV & INJ, BUPC & INJ & COI, ]  <- a_TDP_1[HIV & INJ, BUPC & INJ & COI, ] * p_HIV_HCV_BUPC_INJ # Probability of HCV conditional on HIV
  a_TDP_2[HIV & INJ, BUPC & INJ & COI, ]  <- a_TDP_2[HIV & INJ, BUPC & INJ & COI, ] * p_HIV_HCV_BUPC_INJ
  a_TDP_3[HIV & INJ, BUPC & INJ & COI, ]  <- a_TDP_3[HIV & INJ, BUPC & INJ & COI, ] * p_HIV_HCV_BUPC_INJ
  # From HCV
  a_TDP_1[HCV & INJ, BUPC & INJ & HCV, ]  <- a_TDP_1[HCV & INJ, BUPC & INJ & HCV, ] * (1 - p_HCV_HIV_BUPC_INJ)
  a_TDP_2[HCV & INJ, BUPC & INJ & HCV, ]  <- a_TDP_2[HCV & INJ, BUPC & INJ & HCV, ] * (1 - p_HCV_HIV_BUPC_INJ)
  a_TDP_3[HCV & INJ, BUPC & INJ & HCV, ]  <- a_TDP_3[HCV & INJ, BUPC & INJ & HCV, ] * (1 - p_HCV_HIV_BUPC_INJ)
  a_TDP_1[HCV & INJ, BUPC & INJ & COI, ]  <- a_TDP_1[HCV & INJ, BUPC & INJ & COI, ] * p_HCV_HIV_BUPC_INJ # Probability of HIV conditional on HCV
  a_TDP_2[HCV & INJ, BUPC & INJ & COI, ]  <- a_TDP_2[HCV & INJ, BUPC & INJ & COI, ] * p_HCV_HIV_BUPC_INJ
  a_TDP_3[HCV & INJ, BUPC & INJ & COI, ]  <- a_TDP_3[HCV & INJ, BUPC & INJ & COI, ] * p_HCV_HIV_BUPC_INJ
  
  # MET
  # From NEG
  a_TDP_1[NEG & INJ, MET & INJ & NEG, ]  <- a_TDP_1[NEG & INJ, MET & INJ & NEG, ] * (1 - p_HIV_MET_INJ - p_HCV_MET_INJ)
  a_TDP_2[NEG & INJ, MET & INJ & NEG, ]  <- a_TDP_2[NEG & INJ, MET & INJ & NEG, ] * (1 - p_HIV_MET_INJ - p_HCV_MET_INJ)
  a_TDP_3[NEG & INJ, MET & INJ & NEG, ]  <- a_TDP_3[NEG & INJ, MET & INJ & NEG, ] * (1 - p_HIV_MET_INJ - p_HCV_MET_INJ)
  a_TDP_1[NEG & INJ, MET & INJ & HIV, ]  <- a_TDP_1[NEG & INJ, MET & INJ & HIV, ] * p_HIV_MET_INJ
  a_TDP_2[NEG & INJ, MET & INJ & HIV, ]  <- a_TDP_2[NEG & INJ, MET & INJ & HIV, ] * p_HIV_MET_INJ
  a_TDP_3[NEG & INJ, MET & INJ & HIV, ]  <- a_TDP_3[NEG & INJ, MET & INJ & HIV, ] * p_HIV_MET_INJ
  a_TDP_1[NEG & INJ, MET & INJ & HCV, ]  <- a_TDP_1[NEG & INJ, MET & INJ & HCV, ] * p_HCV_MET_INJ
  a_TDP_2[NEG & INJ, MET & INJ & HCV, ]  <- a_TDP_2[NEG & INJ, MET & INJ & HCV, ] * p_HCV_MET_INJ
  a_TDP_3[NEG & INJ, MET & INJ & HCV, ]  <- a_TDP_3[NEG & INJ, MET & INJ & HCV, ] * p_HCV_MET_INJ
  # From HIV
  a_TDP_1[HIV & INJ, MET & INJ & HIV, ]  <- a_TDP_1[HIV & INJ, MET & INJ & HIV, ] * (1 - p_HIV_HCV_MET_INJ)
  a_TDP_2[HIV & INJ, MET & INJ & HIV, ]  <- a_TDP_2[HIV & INJ, MET & INJ & HIV, ] * (1 - p_HIV_HCV_MET_INJ)
  a_TDP_3[HIV & INJ, MET & INJ & HIV, ]  <- a_TDP_3[HIV & INJ, MET & INJ & HIV, ] * (1 - p_HIV_HCV_MET_INJ)
  a_TDP_1[HIV & INJ, MET & INJ & COI, ]  <- a_TDP_1[HIV & INJ, MET & INJ & COI, ] * p_HIV_HCV_MET_INJ # Probability of HCV conditional on HIV
  a_TDP_2[HIV & INJ, MET & INJ & COI, ]  <- a_TDP_2[HIV & INJ, MET & INJ & COI, ] * p_HIV_HCV_MET_INJ
  a_TDP_3[HIV & INJ, MET & INJ & COI, ]  <- a_TDP_3[HIV & INJ, MET & INJ & COI, ] * p_HIV_HCV_MET_INJ
  # From HCV
  a_TDP_1[HCV & INJ, MET & INJ & HCV, ]  <- a_TDP_1[HCV & INJ, MET & INJ & HCV, ] * (1 - p_HCV_HIV_MET_INJ)
  a_TDP_2[HCV & INJ, MET & INJ & HCV, ]  <- a_TDP_2[HCV & INJ, MET & INJ & HCV, ] * (1 - p_HCV_HIV_MET_INJ)
  a_TDP_3[HCV & INJ, MET & INJ & HCV, ]  <- a_TDP_3[HCV & INJ, MET & INJ & HCV, ] * (1 - p_HCV_HIV_MET_INJ)
  a_TDP_1[HCV & INJ, MET & INJ & COI, ]  <- a_TDP_1[HCV & INJ, MET & INJ & COI, ] * p_HCV_HIV_MET_INJ # Probability of HIV conditional on HCV
  a_TDP_2[HCV & INJ, MET & INJ & COI, ]  <- a_TDP_2[HCV & INJ, MET & INJ & COI, ] * p_HCV_HIV_MET_INJ
  a_TDP_3[HCV & INJ, MET & INJ & COI, ]  <- a_TDP_3[HCV & INJ, MET & INJ & COI, ] * p_HCV_HIV_MET_INJ
  
  # METC
  # From NEG
  a_TDP_1[NEG & INJ, METC & INJ & NEG, ]  <- a_TDP_1[NEG & INJ, METC & INJ & NEG, ] * (1 - p_HIV_METC_INJ - p_HCV_METC_INJ)
  a_TDP_2[NEG & INJ, METC & INJ & NEG, ]  <- a_TDP_2[NEG & INJ, METC & INJ & NEG, ] * (1 - p_HIV_METC_INJ - p_HCV_METC_INJ)
  a_TDP_3[NEG & INJ, METC & INJ & NEG, ]  <- a_TDP_3[NEG & INJ, METC & INJ & NEG, ] * (1 - p_HIV_METC_INJ - p_HCV_METC_INJ)
  a_TDP_1[NEG & INJ, METC & INJ & HIV, ]  <- a_TDP_1[NEG & INJ, METC & INJ & HIV, ] * p_HIV_METC_INJ
  a_TDP_2[NEG & INJ, METC & INJ & HIV, ]  <- a_TDP_2[NEG & INJ, METC & INJ & HIV, ] * p_HIV_METC_INJ
  a_TDP_3[NEG & INJ, METC & INJ & HIV, ]  <- a_TDP_3[NEG & INJ, METC & INJ & HIV, ] * p_HIV_METC_INJ
  a_TDP_1[NEG & INJ, METC & INJ & HCV, ]  <- a_TDP_1[NEG & INJ, METC & INJ & HCV, ] * p_HCV_METC_INJ
  a_TDP_2[NEG & INJ, METC & INJ & HCV, ]  <- a_TDP_2[NEG & INJ, METC & INJ & HCV, ] * p_HCV_METC_INJ
  a_TDP_3[NEG & INJ, METC & INJ & HCV, ]  <- a_TDP_3[NEG & INJ, METC & INJ & HCV, ] * p_HCV_METC_INJ
  # From HIV
  a_TDP_1[HIV & INJ, METC & INJ & HIV, ]  <- a_TDP_1[HIV & INJ, METC & INJ & HIV, ] * (1 - p_HIV_HCV_METC_INJ)
  a_TDP_2[HIV & INJ, METC & INJ & HIV, ]  <- a_TDP_2[HIV & INJ, METC & INJ & HIV, ] * (1 - p_HIV_HCV_METC_INJ)
  a_TDP_3[HIV & INJ, METC & INJ & HIV, ]  <- a_TDP_3[HIV & INJ, METC & INJ & HIV, ] * (1 - p_HIV_HCV_METC_INJ)
  a_TDP_1[HIV & INJ, METC & INJ & COI, ]  <- a_TDP_1[HIV & INJ, METC & INJ & COI, ] * p_HIV_HCV_METC_INJ # Probability of HCV conditional on HIV
  a_TDP_2[HIV & INJ, METC & INJ & COI, ]  <- a_TDP_2[HIV & INJ, METC & INJ & COI, ] * p_HIV_HCV_METC_INJ
  a_TDP_3[HIV & INJ, METC & INJ & COI, ]  <- a_TDP_3[HIV & INJ, METC & INJ & COI, ] * p_HIV_HCV_METC_INJ
  # From HCV
  a_TDP_1[HCV & INJ, METC & INJ & HCV, ]  <- a_TDP_1[HCV & INJ, METC & INJ & HCV, ] * (1 - p_HCV_HIV_METC_INJ)
  a_TDP_2[HCV & INJ, METC & INJ & HCV, ]  <- a_TDP_2[HCV & INJ, METC & INJ & HCV, ] * (1 - p_HCV_HIV_METC_INJ)
  a_TDP_3[HCV & INJ, METC & INJ & HCV, ]  <- a_TDP_3[HCV & INJ, METC & INJ & HCV, ] * (1 - p_HCV_HIV_METC_INJ)
  a_TDP_1[HCV & INJ, METC & INJ & COI, ]  <- a_TDP_1[HCV & INJ, METC & INJ & COI, ] * p_HCV_HIV_METC_INJ # Probability of HIV conditional on HCV
  a_TDP_2[HCV & INJ, METC & INJ & COI, ]  <- a_TDP_2[HCV & INJ, METC & INJ & COI, ] * p_HCV_HIV_METC_INJ
  a_TDP_3[HCV & INJ, METC & INJ & COI, ]  <- a_TDP_3[HCV & INJ, METC & INJ & COI, ] * p_HCV_HIV_METC_INJ
  
  # REL
  # From NEG
  a_TDP_1[NEG & INJ, REL & INJ & NEG, ]  <- a_TDP_1[NEG & INJ, REL & INJ & NEG, ] * (1 - p_HIV_REL_INJ - p_HCV_REL_INJ)
  a_TDP_2[NEG & INJ, REL & INJ & NEG, ]  <- a_TDP_2[NEG & INJ, REL & INJ & NEG, ] * (1 - p_HIV_REL_INJ - p_HCV_REL_INJ)
  a_TDP_3[NEG & INJ, REL & INJ & NEG, ]  <- a_TDP_3[NEG & INJ, REL & INJ & NEG, ] * (1 - p_HIV_REL_INJ - p_HCV_REL_INJ)
  a_TDP_1[NEG & INJ, REL & INJ & HIV, ]  <- a_TDP_1[NEG & INJ, REL & INJ & HIV, ] * p_HIV_REL_INJ
  a_TDP_2[NEG & INJ, REL & INJ & HIV, ]  <- a_TDP_2[NEG & INJ, REL & INJ & HIV, ] * p_HIV_REL_INJ
  a_TDP_3[NEG & INJ, REL & INJ & HIV, ]  <- a_TDP_3[NEG & INJ, REL & INJ & HIV, ] * p_HIV_REL_INJ
  a_TDP_1[NEG & INJ, REL & INJ & HCV, ]  <- a_TDP_1[NEG & INJ, REL & INJ & HCV, ] * p_HCV_REL_INJ
  a_TDP_2[NEG & INJ, REL & INJ & HCV, ]  <- a_TDP_2[NEG & INJ, REL & INJ & HCV, ] * p_HCV_REL_INJ
  a_TDP_3[NEG & INJ, REL & INJ & HCV, ]  <- a_TDP_3[NEG & INJ, REL & INJ & HCV, ] * p_HCV_REL_INJ
  # From HIV
  a_TDP_1[HIV & INJ, REL & INJ & HIV, ]  <- a_TDP_1[HIV & INJ, REL & INJ & HIV, ] * (1 - p_HIV_HCV_REL_INJ)
  a_TDP_2[HIV & INJ, REL & INJ & HIV, ]  <- a_TDP_2[HIV & INJ, REL & INJ & HIV, ] * (1 - p_HIV_HCV_REL_INJ)
  a_TDP_3[HIV & INJ, REL & INJ & HIV, ]  <- a_TDP_3[HIV & INJ, REL & INJ & HIV, ] * (1 - p_HIV_HCV_REL_INJ)
  a_TDP_1[HIV & INJ, REL & INJ & COI, ]  <- a_TDP_1[HIV & INJ, REL & INJ & COI, ] * p_HIV_HCV_REL_INJ # Probability of HCV conditional on HIV
  a_TDP_2[HIV & INJ, REL & INJ & COI, ]  <- a_TDP_2[HIV & INJ, REL & INJ & COI, ] * p_HIV_HCV_REL_INJ
  a_TDP_3[HIV & INJ, REL & INJ & COI, ]  <- a_TDP_3[HIV & INJ, REL & INJ & COI, ] * p_HIV_HCV_REL_INJ
  # From HCV
  a_TDP_1[HCV & INJ, REL & INJ & HCV, ]  <- a_TDP_1[HCV & INJ, REL & INJ & HCV, ] * (1 - p_HCV_HIV_REL_INJ)
  a_TDP_2[HCV & INJ, REL & INJ & HCV, ]  <- a_TDP_2[HCV & INJ, REL & INJ & HCV, ] * (1 - p_HCV_HIV_REL_INJ)
  a_TDP_3[HCV & INJ, REL & INJ & HCV, ]  <- a_TDP_3[HCV & INJ, REL & INJ & HCV, ] * (1 - p_HCV_HIV_REL_INJ)
  a_TDP_1[HCV & INJ, REL & INJ & COI, ]  <- a_TDP_1[HCV & INJ, REL & INJ & COI, ] * p_HCV_HIV_REL_INJ # Probability of HIV conditional on HCV
  a_TDP_2[HCV & INJ, REL & INJ & COI, ]  <- a_TDP_2[HCV & INJ, REL & INJ & COI, ] * p_HCV_HIV_REL_INJ
  a_TDP_3[HCV & INJ, REL & INJ & COI, ]  <- a_TDP_3[HCV & INJ, REL & INJ & COI, ] * p_HCV_HIV_REL_INJ
  
  # ODN
  # From NEG
  a_TDP_1[NEG & INJ, ODN & INJ & NEG, ]  <- a_TDP_1[NEG & INJ, ODN & INJ & NEG, ] * (1 - p_HIV_ODN_INJ - p_HCV_ODN_INJ)
  a_TDP_2[NEG & INJ, ODN & INJ & NEG, ]  <- a_TDP_2[NEG & INJ, ODN & INJ & NEG, ] * (1 - p_HIV_ODN_INJ - p_HCV_ODN_INJ)
  a_TDP_3[NEG & INJ, ODN & INJ & NEG, ]  <- a_TDP_3[NEG & INJ, ODN & INJ & NEG, ] * (1 - p_HIV_ODN_INJ - p_HCV_ODN_INJ)
  a_TDP_1[NEG & INJ, ODN & INJ & HIV, ]  <- a_TDP_1[NEG & INJ, ODN & INJ & HIV, ] * p_HIV_ODN_INJ
  a_TDP_2[NEG & INJ, ODN & INJ & HIV, ]  <- a_TDP_2[NEG & INJ, ODN & INJ & HIV, ] * p_HIV_ODN_INJ
  a_TDP_3[NEG & INJ, ODN & INJ & HIV, ]  <- a_TDP_3[NEG & INJ, ODN & INJ & HIV, ] * p_HIV_ODN_INJ
  a_TDP_1[NEG & INJ, ODN & INJ & HCV, ]  <- a_TDP_1[NEG & INJ, ODN & INJ & HCV, ] * p_HCV_ODN_INJ
  a_TDP_2[NEG & INJ, ODN & INJ & HCV, ]  <- a_TDP_2[NEG & INJ, ODN & INJ & HCV, ] * p_HCV_ODN_INJ
  a_TDP_3[NEG & INJ, ODN & INJ & HCV, ]  <- a_TDP_3[NEG & INJ, ODN & INJ & HCV, ] * p_HCV_ODN_INJ
  # From HIV
  a_TDP_1[HIV & INJ, ODN & INJ & HIV, ]  <- a_TDP_1[HIV & INJ, ODN & INJ & HIV, ] * (1 - p_HIV_HCV_ODN_INJ)
  a_TDP_2[HIV & INJ, ODN & INJ & HIV, ]  <- a_TDP_2[HIV & INJ, ODN & INJ & HIV, ] * (1 - p_HIV_HCV_ODN_INJ)
  a_TDP_3[HIV & INJ, ODN & INJ & HIV, ]  <- a_TDP_3[HIV & INJ, ODN & INJ & HIV, ] * (1 - p_HIV_HCV_ODN_INJ)
  a_TDP_1[HIV & INJ, ODN & INJ & COI, ]  <- a_TDP_1[HIV & INJ, ODN & INJ & COI, ] * p_HIV_HCV_ODN_INJ # Probability of HCV conditional on HIV
  a_TDP_2[HIV & INJ, ODN & INJ & COI, ]  <- a_TDP_2[HIV & INJ, ODN & INJ & COI, ] * p_HIV_HCV_ODN_INJ
  a_TDP_3[HIV & INJ, ODN & INJ & COI, ]  <- a_TDP_3[HIV & INJ, ODN & INJ & COI, ] * p_HIV_HCV_ODN_INJ
  # From HCV
  a_TDP_1[HCV & INJ, ODN & INJ & HCV, ]  <- a_TDP_1[HCV & INJ, ODN & INJ & HCV, ] * (1 - p_HCV_HIV_ODN_INJ)
  a_TDP_2[HCV & INJ, ODN & INJ & HCV, ]  <- a_TDP_2[HCV & INJ, ODN & INJ & HCV, ] * (1 - p_HCV_HIV_ODN_INJ)
  a_TDP_3[HCV & INJ, ODN & INJ & HCV, ]  <- a_TDP_3[HCV & INJ, ODN & INJ & HCV, ] * (1 - p_HCV_HIV_ODN_INJ)
  a_TDP_1[HCV & INJ, ODN & INJ & COI, ]  <- a_TDP_1[HCV & INJ, ODN & INJ & COI, ] * p_HCV_HIV_ODN_INJ # Probability of HIV conditional on HCV
  a_TDP_2[HCV & INJ, ODN & INJ & COI, ]  <- a_TDP_2[HCV & INJ, ODN & INJ & COI, ] * p_HCV_HIV_ODN_INJ
  a_TDP_3[HCV & INJ, ODN & INJ & COI, ]  <- a_TDP_3[HCV & INJ, ODN & INJ & COI, ] * p_HCV_HIV_ODN_INJ
  
  # ODF
  # From NEG
  a_TDP_1[NEG & INJ, ODF & INJ & NEG, ]  <- a_TDP_1[NEG & INJ, ODF & INJ & NEG, ] * 1
  a_TDP_2[NEG & INJ, ODF & INJ & NEG, ]  <- a_TDP_2[NEG & INJ, ODF & INJ & NEG, ] * 1
  a_TDP_3[NEG & INJ, ODF & INJ & NEG, ]  <- a_TDP_3[NEG & INJ, ODF & INJ & NEG, ] * 1
  a_TDP_1[NEG & INJ, ODF & INJ & HIV, ]  <- a_TDP_1[NEG & INJ, ODF & INJ & HIV, ] * 0
  a_TDP_2[NEG & INJ, ODF & INJ & HIV, ]  <- a_TDP_2[NEG & INJ, ODF & INJ & HIV, ] * 0
  a_TDP_3[NEG & INJ, ODF & INJ & HIV, ]  <- a_TDP_3[NEG & INJ, ODF & INJ & HIV, ] * 0
  a_TDP_1[NEG & INJ, ODF & INJ & HCV, ]  <- a_TDP_1[NEG & INJ, ODF & INJ & HCV, ] * 0
  a_TDP_2[NEG & INJ, ODF & INJ & HCV, ]  <- a_TDP_2[NEG & INJ, ODF & INJ & HCV, ] * 0
  a_TDP_3[NEG & INJ, ODF & INJ & HCV, ]  <- a_TDP_3[NEG & INJ, ODF & INJ & HCV, ] * 0
  # From HIV
  a_TDP_1[HIV & INJ, ODF & INJ & HIV, ]  <- a_TDP_1[HIV & INJ, ODF & INJ & HIV, ] * 1
  a_TDP_2[HIV & INJ, ODF & INJ & HIV, ]  <- a_TDP_2[HIV & INJ, ODF & INJ & HIV, ] * 1
  a_TDP_3[HIV & INJ, ODF & INJ & HIV, ]  <- a_TDP_3[HIV & INJ, ODF & INJ & HIV, ] * 1
  a_TDP_1[HIV & INJ, ODF & INJ & COI, ]  <- a_TDP_1[HIV & INJ, ODF & INJ & COI, ] * 0
  a_TDP_2[HIV & INJ, ODF & INJ & COI, ]  <- a_TDP_2[HIV & INJ, ODF & INJ & COI, ] * 0
  a_TDP_3[HIV & INJ, ODF & INJ & COI, ]  <- a_TDP_3[HIV & INJ, ODF & INJ & COI, ] * 0
  # From HCV
  a_TDP_1[HCV & INJ, ODF & INJ & HCV, ]  <- a_TDP_1[HCV & INJ, ODF & INJ & HCV, ] * 1
  a_TDP_2[HCV & INJ, ODF & INJ & HCV, ]  <- a_TDP_2[HCV & INJ, ODF & INJ & HCV, ] * 1
  a_TDP_3[HCV & INJ, ODF & INJ & HCV, ]  <- a_TDP_3[HCV & INJ, ODF & INJ & HCV, ] * 1
  a_TDP_1[HCV & INJ, ODF & INJ & COI, ]  <- a_TDP_1[HCV & INJ, ODF & INJ & COI, ] * 0
  a_TDP_2[HCV & INJ, ODF & INJ & COI, ]  <- a_TDP_2[HCV & INJ, ODF & INJ & COI, ] * 0
  a_TDP_3[HCV & INJ, ODF & INJ & COI, ]  <- a_TDP_3[HCV & INJ, ODF & INJ & COI, ] * 0
  
  # ABS
  # From NEG
  a_TDP_1[NEG & INJ, ABS & INJ & NEG, ]  <- a_TDP_1[NEG & INJ, ABS & INJ & NEG, ] * (1 - p_HIV_ABS_INJ - p_HCV_ABS_INJ)
  a_TDP_2[NEG & INJ, ABS & INJ & NEG, ]  <- a_TDP_2[NEG & INJ, ABS & INJ & NEG, ] * (1 - p_HIV_ABS_INJ - p_HCV_ABS_INJ)
  a_TDP_3[NEG & INJ, ABS & INJ & NEG, ]  <- a_TDP_3[NEG & INJ, ABS & INJ & NEG, ] * (1 - p_HIV_ABS_INJ - p_HCV_ABS_INJ)
  a_TDP_1[NEG & INJ, ABS & INJ & HIV, ]  <- a_TDP_1[NEG & INJ, ABS & INJ & HIV, ] * p_HIV_ABS_INJ
  a_TDP_2[NEG & INJ, ABS & INJ & HIV, ]  <- a_TDP_2[NEG & INJ, ABS & INJ & HIV, ] * p_HIV_ABS_INJ
  a_TDP_3[NEG & INJ, ABS & INJ & HIV, ]  <- a_TDP_3[NEG & INJ, ABS & INJ & HIV, ] * p_HIV_ABS_INJ
  a_TDP_1[NEG & INJ, ABS & INJ & HCV, ]  <- a_TDP_1[NEG & INJ, ABS & INJ & HCV, ] * p_HCV_ABS_INJ
  a_TDP_2[NEG & INJ, ABS & INJ & HCV, ]  <- a_TDP_2[NEG & INJ, ABS & INJ & HCV, ] * p_HCV_ABS_INJ
  a_TDP_3[NEG & INJ, ABS & INJ & HCV, ]  <- a_TDP_3[NEG & INJ, ABS & INJ & HCV, ] * p_HCV_ABS_INJ
  # From HIV
  a_TDP_1[HIV & INJ, ABS & INJ & HIV, ]  <- a_TDP_1[HIV & INJ, ABS & INJ & HIV, ] * (1 - p_HIV_HCV_ABS_INJ)
  a_TDP_2[HIV & INJ, ABS & INJ & HIV, ]  <- a_TDP_2[HIV & INJ, ABS & INJ & HIV, ] * (1 - p_HIV_HCV_ABS_INJ)
  a_TDP_3[HIV & INJ, ABS & INJ & HIV, ]  <- a_TDP_3[HIV & INJ, ABS & INJ & HIV, ] * (1 - p_HIV_HCV_ABS_INJ)
  a_TDP_1[HIV & INJ, ABS & INJ & COI, ]  <- a_TDP_1[HIV & INJ, ABS & INJ & COI, ] * p_HIV_HCV_ABS_INJ # Probability of HCV conditional on HIV
  a_TDP_2[HIV & INJ, ABS & INJ & COI, ]  <- a_TDP_2[HIV & INJ, ABS & INJ & COI, ] * p_HIV_HCV_ABS_INJ
  a_TDP_3[HIV & INJ, ABS & INJ & COI, ]  <- a_TDP_3[HIV & INJ, ABS & INJ & COI, ] * p_HIV_HCV_ABS_INJ
  # From HCV
  a_TDP_1[HCV & INJ, ABS & INJ & HCV, ]  <- a_TDP_1[HCV & INJ, ABS & INJ & HCV, ] * (1 - p_HCV_HIV_ABS_INJ)
  a_TDP_2[HCV & INJ, ABS & INJ & HCV, ]  <- a_TDP_2[HCV & INJ, ABS & INJ & HCV, ] * (1 - p_HCV_HIV_ABS_INJ)
  a_TDP_3[HCV & INJ, ABS & INJ & HCV, ]  <- a_TDP_3[HCV & INJ, ABS & INJ & HCV, ] * (1 - p_HCV_HIV_ABS_INJ)
  a_TDP_1[HCV & INJ, ABS & INJ & COI, ]  <- a_TDP_1[HCV & INJ, ABS & INJ & COI, ] * p_HCV_HIV_ABS_INJ # Probability of HIV conditional on HCV
  a_TDP_2[HCV & INJ, ABS & INJ & COI, ]  <- a_TDP_2[HCV & INJ, ABS & INJ & COI, ] * p_HCV_HIV_ABS_INJ
  a_TDP_3[HCV & INJ, ABS & INJ & COI, ]  <- a_TDP_3[HCV & INJ, ABS & INJ & COI, ] * p_HCV_HIV_ABS_INJ

  #### Disallowed transitions (ensure that impossible transitions are set to zero) ####
  # Episode rules
  a_TDP_1[EP1, EP3, ] <- a_TDP_2[EP1, EP3, ] <- a_TDP_3[EP1, EP3, ] <- 0
  a_TDP_1[EP2, EP1, ] <- a_TDP_2[EP2, EP1, ] <- a_TDP_3[EP2, EP1, ] <- 0
  a_TDP_1[EP3, EP1, ] <- a_TDP_2[EP3, EP1, ] <- a_TDP_3[EP3, EP1, ] <- 0
  a_TDP_1[EP3, EP2, ] <- a_TDP_2[EP3, EP2, ] <- a_TDP_3[EP3, EP2, ] <- 0
  
  # Seroconversions
  a_TDP_1[HIV, NEG, ] <- a_TDP_2[HIV, NEG, ] <- a_TDP_3[HIV, NEG, ] <- 0
  a_TDP_1[HCV, NEG, ] <- a_TDP_2[HCV, NEG, ] <- a_TDP_3[HCV, NEG, ] <- 0 
  a_TDP_1[COI, NEG, ] <- a_TDP_2[COI, NEG, ] <- a_TDP_3[COI, NEG, ] <- 0
  a_TDP_1[HIV, HCV, ] <- a_TDP_2[HIV, HCV, ] <- a_TDP_3[HIV, HCV, ] <- 0
  a_TDP_1[HCV, HIV, ] <- a_TDP_2[HCV, HIV, ] <- a_TDP_3[HCV, HIV, ] <- 0
  a_TDP_1[COI, HIV, ] <- a_TDP_2[COI, HIV, ] <- a_TDP_3[COI, HIV, ] <- 0 
  a_TDP_1[COI, HCV, ] <- a_TDP_2[COI, HCV, ] <- a_TDP_3[COI, HCV, ] <- 0
  a_TDP_1[COI, NEG, ] <- a_TDP_2[COI, NEG, ] <- a_TDP_3[COI, NEG, ] <- 0
  a_TDP_1[NEG, COI, ] <- a_TDP_2[NEG, COI, ] <- a_TDP_3[NEG, COI, ] <- 0
  
  # Conditional transitions
  # Next episode with out-of-treatment(OOT) EPi -> treatment(TX) EP(i+1)
  a_TDP_1[all_TX & EP1, all_TX & EP2, ] <- a_TDP_2[all_TX & EP1, all_TX & EP2, ] <- a_TDP_3[all_TX & EP1, all_TX & EP2, ] <- 0
  a_TDP_1[all_TX & EP1, all_TX & EP3, ] <- a_TDP_2[all_TX & EP1, all_TX & EP3, ] <- a_TDP_3[all_TX & EP1, all_TX & EP3, ] <- 0
  a_TDP_1[all_TX & EP2, all_TX & EP3, ] <- a_TDP_2[all_TX & EP2, all_TX & EP3, ] <- a_TDP_3[all_TX & EP2, all_TX & EP3, ] <- 0
  a_TDP_1[OOT & EP1, OOT & EP2, ] <- a_TDP_2[OOT & EP1, OOT & EP2, ] <- a_TDP_3[OOT & EP1, OOT & EP2, ] <- 0
  a_TDP_1[OOT & EP1, OOT & EP3, ] <- a_TDP_2[OOT & EP1, OOT & EP3, ] <- a_TDP_3[OOT & EP1, OOT & EP3, ] <- 0
  a_TDP_1[OOT & EP2, OOT & EP3, ] <- a_TDP_2[OOT & EP2, OOT & EP3, ] <- a_TDP_3[OOT & EP2, OOT & EP3, ] <- 0
  a_TDP_1[all_TX & EP1, OOT & EP2, ] <- a_TDP_2[all_TX & EP1, OOT & EP2, ] <- a_TDP_3[all_TX & EP1, OOT & EP2, ] <- 0
  a_TDP_1[all_TX & EP2, OOT & EP3, ] <- a_TDP_2[all_TX & EP2, OOT & EP3, ] <- a_TDP_3[all_TX & EP2, OOT & EP3, ] <- 0
  a_TDP_1[OOT & EP1, all_TX & EP1, ] <- a_TDP_2[OOT & EP1, all_TX & EP1, ] <- a_TDP_3[OOT & EP1, all_TX & EP1, ] <- 0
  a_TDP_1[OOT & EP2, all_TX & EP2, ] <- a_TDP_2[OOT & EP2, all_TX & EP2, ] <- a_TDP_3[OOT & EP2, all_TX & EP2, ] <- 0
  
  if(checks){
    # State transitions
    # First month (state-time)
    write.csv(a_UP_first[, , 1], "checks/state transitions/a_UP_first_2018.csv")
    write.csv(a_UP_first[, , 2], "checks/state transitions/a_UP_first_2019.csv")
    write.csv(a_UP_first[, , 3], "checks/state transitions/a_UP_first_2020.csv")
    
    # Month 2+ (state-time)
    write.csv(a_UP[, , 1], "checks/state transitions/a_UP_2018.csv")
    write.csv(a_UP[, , 2], "checks/state transitions/a_UP_2019.csv")
    write.csv(a_UP[, , 3], "checks/state transitions/a_UP_2020.csv")
    
    # Full array at time = 1
    array_2018_1m <- a_TDP_1[, , 1]
    array_2019_1m <- a_TDP_2[, , 1]
    array_2020_1m <- a_TDP_3[, , 1]
    write.csv(array_2018_1m,"checks/full array/array_2018_1m.csv", row.names = TRUE)
    write.csv(array_2019_1m,"checks/full array/array_2019_1m.csv", row.names = TRUE)
    write.csv(array_2020_1m,"checks/full array/array_2020_1m.csv", row.names = TRUE)
    
    # Full array at time = 24 months
    array_2018_24m <- a_TDP_1[, , 24]
    array_2019_24m <- a_TDP_2[, , 24]
    array_2020_24m <- a_TDP_3[, , 24]
    write.csv(array_2018_24m,"checks/full array/array_2018_24m.csv", row.names = TRUE)
    write.csv(array_2019_24m,"checks/full array/array_2019_24m.csv", row.names = TRUE)
    write.csv(array_2020_24m,"checks/full array/array_2020_24m.csv", row.names = TRUE)
    
    # Full array at time = max
    array_2018_last <- a_TDP_1[, , n_t]
    array_2019_last <- a_TDP_2[, , n_t]
    array_2020_last <- a_TDP_3[, , n_t]
    write.csv(array_2018_last,"checks/full array/array_2018_last.csv", row.names = TRUE)
    write.csv(array_2019_last,"checks/full array/array_2019_last.csv", row.names = TRUE)
    write.csv(array_2020_last,"checks/full array/array_2020_last.csv", row.names = TRUE)
  } else{}

  #### Check transition array ####
  #check_transition_probability(a_P = a_TDP, err_stop = err_stop, verbose = verbose) # check all probs [0, 1]
  #check_sum_of_transition_array(a_P = a_TDP, n_states = n_states, n_t = n_t, err_stop = err_stop, verbose = verbose) # check prob sums = 1

  #### Run Markov model ####
  # Create empty initial state vectors
  v_s_init <- rep(0, n_states)
  names(v_s_init) <- v_n_states

  #### Set initial state vector ####
  # Baseline
  # Populate baseline states
  # Think about this e.g. (all start in EP1, or according to observed %'s in each episode)
  # Episode 1
  v_s_init[BUP & EP1]  <- v_init_dist["pe", "BUP"] # Empirically observed proportions from base states
  v_s_init[BUPC & EP1] <- v_init_dist["pe", "BUPC"]
  v_s_init[MET & EP1]  <- v_init_dist["pe", "MET"]
  v_s_init[METC & EP1] <- v_init_dist["pe", "METC"]
  v_s_init[REL & EP1]  <- v_init_dist["pe", "REL"]
  v_s_init[ODN & EP1]  <- v_init_dist["pe", "ODN"]
  v_s_init[ODF & EP1]  <- v_init_dist["pe", "ODF"]
  v_s_init[ABS & EP1]  <- v_init_dist["pe", "ABS"]
  
  # Distribute by injection/non-injection
  v_s_init[NI]  <- v_s_init[NI] * (1 - n_INJ)
  v_s_init[INJ] <- v_s_init[INJ] * n_INJ
  
  # Distribute HIV/HCV/COI
  # Injection
  v_s_init[NEG & INJ] <- v_s_init[NEG & INJ] * (1 - n_HIV_INJ - n_HCV_INJ - n_COI_INJ)
  v_s_init[HIV & INJ] <- v_s_init[HIV & INJ] * n_HIV_INJ
  v_s_init[HCV & INJ] <- v_s_init[HCV & INJ] * n_HCV_INJ
  v_s_init[COI & INJ] <- v_s_init[COI & INJ] * n_COI_INJ
  
  # Non-injection
  v_s_init[NEG & NI] <- v_s_init[NEG & NI] * (1 - n_HIV_NI - n_HCV_NI - n_COI_NI)
  v_s_init[HIV & NI] <- v_s_init[HIV & NI] * n_HIV_NI
  v_s_init[HCV & NI] <- v_s_init[HCV & NI] * n_HCV_NI
  v_s_init[COI & NI] <- v_s_init[COI & NI] * n_COI_NI
  
  if(checks){
    write.csv(v_s_init,"checks/initial/v_s_init.csv", row.names = TRUE)
  }
  
  # Create Markov Trace
  # Initialize population
  a_M_trace <- array(0, dim = c((n_t + 1), n_states, (n_t + 1)),
                     dimnames = list(0:n_t, v_n_states, 0:n_t))
  a_M_trace[1, , 1] <- v_s_init
    
  # Calbrate over three time points (2018, 2019, 2020)
  if(cali == TRUE){
    # All model time periods
    # Split into three periods 0-12; 12-24; 24+ to account for overall time-varying overdose parameters
    # Months 0-12 (model-time)
      for(i in 2:11){
        # Time spent in given health state
        # First month (state-time)
        for(j in 1:(i - 1)){
          # state-time-dependent transition probability (j) * age (model-time)-specific mortality (i) * model-time-specific overdose (track in separate matrix)
          m_sojourn <- a_TDP_1[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point
          v_current_state <- as.vector(a_M_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_M_trace[i, ,j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_M_trace[i, ,1] <- v_new_state + a_M_trace[i, ,1] # add new state %'s to array
        }
      }
    # Months 12-23 (model-time)
      for(i in 12:23){
        # Time spent in given health state
        # First month (state-time)
        for(j in 1:(i - 1)){
          # state-time-dependent transition probability (j) * age (model-time)-specific mortality (i) * model-time-specific overdose (track in separate matrix)
          m_sojourn <- a_TDP_2[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point
          v_current_state <- as.vector(a_M_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_M_trace[i, ,j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_M_trace[i, ,1] <- v_new_state + a_M_trace[i, ,1] # add new state %'s to array
        }
      }
    # Months 24+ (model-time)
      for(i in 24:(n_t)){
        # Time spent in given health state
        # First month (state-time)
        for(j in 1:(i - 1)){
          m_sojourn <- a_TDP_3[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_M_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_M_trace[i, ,j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_M_trace[i, ,1] <- v_new_state + a_M_trace[i, ,1] # add new state %'s to array
        }
      }
    } else{ # For non-calibration, run model for 2020+
      # Full periods (model-time)
      for(i in 2:(n_t)){
        # Time spent in given health state
        # First month (state-time)
        for(j in 1:(i - 1)){
          m_sojourn <- a_TDP_3[, , j] * m_alive[, i - 1] # state-time transition matrix at state-time j, re-weighted for model-time (age) varying mortality at each time point, time-varying overdose (t=3)
          v_current_state <- as.vector(a_M_trace[i - 1, , j]) # all in current state
          v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining in state next period
          a_M_trace[i, ,j + 1] <- v_same_state # add remain to next period
          diag(m_sojourn) <- 0 # reset remain to 0 once counted
          v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states post-transition (excluding remaining)
          a_M_trace[i, ,1] <- v_new_state + a_M_trace[i, ,1] # add new state %'s to array
        }
      }
    }
    
  # Collect trace for time-periods across all model states  
  m_M_trace <- array(0, dim = c((n_t + 1), n_states),
                     dimnames = list(0:n_t, v_n_states))
  for (i in 1:n_t){
    m_M_trace[i, ] <- rowSums(a_M_trace[i, ,])
  }

  # Count cumulative state-specific deaths
  m_M_trace_death <- array(0, dim = c((n_t + 1), n_states),
                           dimnames = list(0:n_t, v_n_states))
  for (i in 2:n_t){
    m_M_trace_death[i, ] <- m_M_trace[i - 1, ] * m_mort[, i - 1] # State-specific deaths at each time point as function of state-occupancy in t-1
  }
  m_M_trace_cumsum_death <- apply(m_M_trace_death, 2, cumsum) # Cumulative non-overdose deaths at each time point (use m_M_trace_death for individual period deaths)
  
  #### Create aggregated trace matrices ####
  v_agg_trace_states <- c("Alive", "Death", "ODN", "ODF", "REL", "BUP", "BUPC", "MET", "METC", "ABS") # states to aggregate
  v_agg_trace_death_states <- c("Total", "ODN", "ODF", "REL", "BUP", "BUPC", "MET", "METC", "ABS") # states to aggregate
  v_agg_trace_sero_states <- c("NEG-Alive", "HIV-Alive", "HCV-Alive", "COI-Alive", "NEG-Dead", "HIV-Dead", "HCV-Dead", "COI-Dead") # states to aggregate
  
  n_agg_trace_states <- length(v_agg_trace_states)
  n_agg_trace_death_states <- length(v_agg_trace_death_states)
  n_agg_trace_sero_states <- length(v_agg_trace_sero_states)
  
  m_M_agg_trace       <- array(0, dim = c((n_t + 1), n_agg_trace_states),
                               dimnames = list(0:n_t, v_agg_trace_states))
  m_M_agg_trace_death <- array(0, dim = c((n_t + 1), n_agg_trace_death_states),
                               dimnames = list(0:n_t, v_agg_trace_death_states))
  m_M_agg_trace_sero  <- array(0, dim = c((n_t + 1), n_agg_trace_sero_states),
                               dimnames = list(0:n_t, v_agg_trace_sero_states))
  
  for (i in 1:n_t){
    #m_M_agg_trace[i, "Alive"] <- sum(m_M_trace[i, ])
    m_M_agg_trace[i, "BUP"]   <- sum(m_M_trace[i, BUP])
    m_M_agg_trace[i, "BUPC"]  <- sum(m_M_trace[i, BUPC])
    m_M_agg_trace[i, "MET"]   <- sum(m_M_trace[i, MET])
    m_M_agg_trace[i, "METC"]  <- sum(m_M_trace[i, METC])
    m_M_agg_trace[i, "REL"]   <- sum(m_M_trace[i, REL])
    m_M_agg_trace[i, "ABS"]   <- sum(m_M_trace[i, ABS])
    m_M_agg_trace[i, "ODN"]   <- sum(m_M_trace[i, ODN])
    m_M_agg_trace[i, "ODF"]   <- sum(m_M_trace[i, ODF])
    m_M_agg_trace[i, "Death"] <- 1 - sum(m_M_trace[i, ])
  }
  # Aggregated state occupancy matrix
  write.csv(m_M_agg_trace,"outputs/trace/m_M_agg_trace.csv", row.names = TRUE)
  
  for (i in 1:n_t){
    m_M_agg_trace_death[i, "Total"] <- sum(m_M_trace_cumsum_death[i, ])
    m_M_agg_trace_death[i, "ODN"]   <- sum(m_M_trace_cumsum_death[i, ODN])
    m_M_agg_trace_death[i, "ODF"]   <- sum(m_M_trace_cumsum_death[i, ODF])
    m_M_agg_trace_death[i, "REL"]   <- sum(m_M_trace_cumsum_death[i, REL])
    m_M_agg_trace_death[i, "BUP"]   <- sum(m_M_trace_cumsum_death[i, BUP])
    m_M_agg_trace_death[i, "BUPC"]  <- sum(m_M_trace_cumsum_death[i, BUPC])
    m_M_agg_trace_death[i, "MET"]   <- sum(m_M_trace_cumsum_death[i, MET])
    m_M_agg_trace_death[i, "METC"]  <- sum(m_M_trace_cumsum_death[i, METC])
    m_M_agg_trace_death[i, "ABS"]   <- sum(m_M_trace_cumsum_death[i, ABS])
  }
  
  for (i in 1:n_t){
    m_M_agg_trace_sero[i, "NEG-Alive"] <- sum(m_M_trace[i, NEG])
    m_M_agg_trace_sero[i, "HIV-Alive"] <- sum(m_M_trace[i, HIV])
    m_M_agg_trace_sero[i, "HCV-Alive"] <- sum(m_M_trace[i, HCV])
    m_M_agg_trace_sero[i, "COI-Alive"] <- sum(m_M_trace[i, COI])
    m_M_agg_trace_sero[i, "NEG-Dead"]  <- sum(m_M_trace_cumsum_death[i, NEG])
    m_M_agg_trace_sero[i, "HIV-Dead"]  <- sum(m_M_trace_cumsum_death[i, HIV])
    m_M_agg_trace_sero[i, "HCV-Dead"]  <- sum(m_M_trace_cumsum_death[i, HCV])
    m_M_agg_trace_sero[i, "COI-Dead"]  <- sum(m_M_trace_cumsum_death[i, COI])
  }
  
  return(list(l_index_s = l_index_s,
              #a_TDP = a_TDP,
              m_M_trace = m_M_trace,
              m_M_agg_trace = m_M_agg_trace,
              m_M_agg_trace_death = m_M_agg_trace_death,
              m_M_agg_trace_sero = m_M_agg_trace_sero))
  }
 )
}

#' Check if transition array is valid
#'
#' \code{check_transition_probability} checks if individual transition probabilities are in range \[0, 1\].
#'
#' @param a_P A transition probability array.
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. 
#' Default = FALSE
#'
#' @return
#' This function stops if transition probability array is not valid and shows which entries are invalid
#' @import utils
#' @export
check_transition_probability <- function(a_P,
                                         err_stop = FALSE, 
                                         verbose = FALSE) {
  
  m_indices_notvalid <- arrayInd(which(a_P < 0 | a_P > 1), 
                                 dim(a_P))
  
  if(dim(m_indices_notvalid)[1] != 0){
    v_rows_notval   <- rownames(a_P)[m_indices_notvalid[, 1]]
    v_cols_notval   <- colnames(a_P)[m_indices_notvalid[, 2]]
    v_cycles_notval <- dimnames(a_P)[[3]][m_indices_notvalid[, 3]]
    
    df_notvalid <- data.frame(`Transition probabilities not valid:` = 
                                matrix(paste0(paste(v_rows_notval, v_cols_notval, sep = "->"),
                                              "; at cycle ",
                                              v_cycles_notval), ncol = 1), 
                              check.names = FALSE)
    
    if(err_stop) {
      stop("Not valid transition probabilities\n",
           paste(capture.output(df_notvalid), collapse = "\n"))
    }
    
    if(verbose){
      warning("Not valid transition probabilities\n",
              paste(capture.output(df_notvalid), collapse = "\n"))
    } 
  }
}

#' Check if the sum of transition probabilities from each state are equal to one. 
#'
#' \code{check_sum_of_transition_array} checks if each of the rows of the 
#' transition matrices sum to one. 
#' 
#' @param a_P Transition probability array.
#' @param n_states Number of health states.
#' @param n_t Number of time periods.
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE.
#' @return 
#' The transition probability array and the cohort trace matrix.
#' @import dplyr
#' @export
check_sum_of_transition_array <- function(a_P,
                                          n_states,
                                          n_t,  
                                          err_stop = FALSE, 
                                          verbose = FALSE) {
  
  valid <- (apply(a_P, 3, function(x) sum(rowSums(x))) == n_states)
  if (!isTRUE(all_equal(as.numeric(sum(valid)), as.numeric(n_t)))) {
    if(err_stop) {
      stop("This is not a valid transition Matrix")
    }
    
    if(verbose){
      warning("This is not a valid transition Matrix")
    } 
  }
}
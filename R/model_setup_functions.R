#' Markov model
#'
#' \code{markov_model} implements the main model functions to calculate Markov trace.
#'
#' @param l_params_all List with all parameters
#' @param err_stop Logical variable to stop model run if transition array is invalid, if TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE
#' @return 
#' a_TDP: Transition probability array
#' m_M_trace: Fully stratified markov cohort trace
#' m_M_agg_trace: Aggregated markov trace over base health states
#' m_M_agg_trace_death: State-specific mortality from each health state
#' m_M_agg_trace_sero: HIV seroconversions from each health state
#' @export
markov_model <- function(l_params_all, err_stop = FALSE, verbose = FALSE, checks = FALSE){
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
  # Maximum model periods
  n_t <- (n_age_max - n_age_init) * n_per # convert years
  
  df_flat <- expand.grid(l_dim_s) #combine all elements together into vector of health states
  df_flat <- rename(df_flat, BASE    = Var1, 
                             INJECT  = Var2, 
                             EP      = Var3, 
                             SERO    = Var4)

  # Create index of states to populate transition matrices
  # All treatment
  TX <- df_flat$BASE == "BUP" | df_flat$BASE == "BUPC" | df_flat$BASE == "MET" | df_flat$BASE == "METC"
  
  # All out-of-treatment (incl ABS)
  OOT <- df_flat$BASE == "REL" | df_flat$BASE == "ODN" | df_flat$BASE == "ODF" | df_flat$BASE == "ABS"
  
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
  #' @param first_month Logical parameter to switch between month 1 and month 2+ for parameter estimation
  #' @param fatal Logical parameter to switch between fatal/non-fatal overdose
  #' @param injection Logical parameter to adjust rate calculation for injection/non-injection use
  #' 
  #' @return 
  #' `p_OD` monthly probability of fatal or non-fatal overdose from a given health state
  #' @export
  p_OD <- function(rate = rate,
                   rate_fatal = n_fatal_OD,
                   rate_fent = n_fent_OD,
                   multiplier = multiplier,
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
    # Convert input monthly rates to monthly probabilities - multiply rates by first month multiplier before converting
    if (injection == TRUE && first_month == TRUE){
      p_base_OD <- 1 - exp(-(rate * n_INJ_OD_mult * multiplier))
      p_fent_OD <- 1 - exp(-(rate_fent * multiplier))
    }
    else if (injection == TRUE && first_month == FALSE){
      p_base_OD <- 1 - exp(-(rate * n_INJ_OD_mult))
      p_fent_OD <- 1 - exp(-(rate_fent))
    }
    else if (injection == FALSE && first_month == TRUE){
      p_base_OD <- 1 - exp(-(rate * multiplier))
      p_fent_OD <- 1 - exp(-(rate_fent * multiplier))
    }
    else if (injection == FALSE && first_month == FALSE){
      p_base_OD <- 1 - exp(-(rate))
      p_fent_OD <- 1 - exp(-(rate_fent))
    }

    # Naloxone effect on fatal overdose
    p_fatal_OD_NX <- p_fatal_OD * (1 - p_NX_rev)
    # Probability of fentanyl exposure (adjusted for injection/non-injection)
    if (injection){
    p_fent_exp <- p_fent_exp
    }  else{
      p_fent_exp <- p_fent_exp * p_ni_fent_reduction
    }
    
    # Calculate fatal and non-fatal overdose probabilities
    if (fatal){
      p_OD <- ((p_base_OD * (1 - p_fent_exp)) + (p_fent_OD * (p_fent_exp))) * p_fatal_OD_NX
    } else{
      p_OD <- ((p_base_OD * (1 - p_fent_exp)) + (p_fent_OD * (p_fent_exp))) * (1 - p_fatal_OD_NX)
    }
    return(p_OD)
  }

  # Module to calculate probability of overdose from states
  # Probability of overdose
  # Non-injection
  # Non-fatal
  p_BUP_ODN_NI  <- p_MET_ODN_NI <- p_OD(rate = n_TX_OD, first_month = FALSE, fatal = FALSE, injection = FALSE)
  p_BUPC_ODN_NI <- p_METC_ODN_NI <- p_OD(rate = n_TXC_OD, first_month = FALSE, fatal = FALSE, injection = FALSE)
  p_REL_ODN_NI  <- p_ODN_ODN_NI <- p_OD(rate = n_REL_OD, first_month = FALSE, fatal = FALSE, injection = FALSE)
  p_ABS_ODN_NI  <- p_OD(rate = n_ABS_OD, first_month = FALSE, fatal = FALSE, injection = FALSE)

  # Fatal
  p_BUP_ODF_NI  <- p_MET_ODF_NI <- p_OD(rate = n_TX_OD, first_month = FALSE, fatal = TRUE, injection = FALSE)
  p_BUPC_ODF_NI <- p_METC_ODF_NI <- p_OD(rate = n_TXC_OD, first_month = FALSE, fatal = TRUE, injection = FALSE)
  p_REL_ODF_NI  <- p_ODN_ODF_NI <- p_OD(rate = n_REL_OD, first_month = FALSE, fatal = TRUE, injection = FALSE)
  p_ABS_ODF_NI  <- p_OD(rate = n_ABS_OD, first_month = FALSE, fatal = TRUE, injection = FALSE)

  # Injection
  # Non-fatal
  p_BUP_ODN_INJ  <- p_MET_ODN_INJ <- p_OD(rate = n_TX_OD, first_month = FALSE, fatal = FALSE, injection = TRUE)
  p_BUPC_ODN_INJ <- p_METC_ODN_INJ <- p_OD(rate = n_TXC_OD, first_month = FALSE, fatal = FALSE, injection = TRUE)
  p_REL_ODN_INJ  <- p_ODN_ODN_INJ <- p_OD(rate = n_REL_OD, first_month = FALSE, fatal = FALSE, injection = TRUE)
  p_ABS_ODN_INJ  <- p_OD(rate = n_ABS_OD, first_month = FALSE, fatal = FALSE, injection = TRUE)

  # Fatal
  p_BUP_ODF_INJ  <- p_MET_ODF_INJ <- p_OD(rate = n_TX_OD, first_month = FALSE, fatal = TRUE, injection = TRUE)
  p_BUPC_ODF_INJ <- p_METC_ODF_INJ <- p_OD(rate = n_TXC_OD, first_month = FALSE, fatal = TRUE, injection = TRUE)
  p_REL_ODF_INJ  <- p_ODN_ODF_INJ <- p_OD(rate = n_REL_OD, first_month = FALSE, fatal = TRUE, injection = TRUE)
  p_ABS_ODF_INJ  <- p_OD(rate = n_ABS_OD, first_month = FALSE, fatal = TRUE, injection = TRUE)

  # Probability of overdose (first month multiplier)
  # No multiplier for OD -> OD as individuals have already spent at least 1 month in relapse
  # Non-injection
  # Non-fatal
  p_BUP_ODN_NI_4wk  <- p_MET_ODN_NI_4wk <- p_OD(rate = n_TX_OD, multiplier = n_TX_OD_mult, first_month = TRUE, fatal = FALSE, injection = FALSE)
  p_BUPC_ODN_NI_4wk <- p_METC_ODN_NI_4wk <- p_OD(rate = n_TXC_OD, multiplier = n_TXC_OD_mult, first_month = TRUE, fatal = FALSE, injection = FALSE)
  p_REL_ODN_NI_4wk  <- p_OD(rate = n_REL_OD, multiplier = n_REL_OD_mult, first_month = TRUE, fatal = FALSE, injection = FALSE)
  p_ABS_ODN_NI_4wk  <- p_OD(rate = n_ABS_OD, multiplier = n_ABS_OD_mult, first_month = TRUE, fatal = FALSE, injection = FALSE)
  
  # Fatal
  p_BUP_ODF_NI_4wk  <- p_MET_ODF_NI_4wk <- p_OD(rate = n_TX_OD, multiplier = n_TX_OD_mult, first_month = TRUE, fatal = TRUE, injection = FALSE)
  p_BUPC_ODF_NI_4wk <- p_METC_ODF_NI_4wk <- p_OD(rate = n_TXC_OD, multiplier = n_TXC_OD_mult, first_month = TRUE, fatal = TRUE, injection = FALSE)
  p_REL_ODF_NI_4wk  <- p_OD(rate = n_REL_OD, multiplier = n_REL_OD_mult, first_month = TRUE, fatal = TRUE, injection = FALSE)
  p_ABS_ODF_NI_4wk  <- p_OD(rate = n_ABS_OD, multiplier = n_ABS_OD_mult, first_month = TRUE, fatal = TRUE, injection = FALSE)
  
  # Injection
  # Non-fatal
  p_BUP_ODN_INJ_4wk  <- p_MET_ODN_INJ_4wk <- p_OD(rate = n_TX_OD, multiplier = n_TX_OD_mult, first_month = TRUE, fatal = FALSE, injection = TRUE)
  p_BUPC_ODN_INJ_4wk <- p_METC_ODN_INJ_4wk <- p_OD(rate = n_TXC_OD, multiplier = n_TXC_OD_mult, first_month = TRUE, fatal = FALSE, injection = TRUE)
  p_REL_ODN_INJ_4wk  <- p_OD(rate = n_REL_OD, multiplier = n_REL_OD_mult, first_month = TRUE, fatal = FALSE, injection = TRUE)
  p_ABS_ODN_INJ_4wk  <- p_OD(rate = n_ABS_OD, multiplier = n_ABS_OD_mult, first_month = TRUE, fatal = FALSE, injection = TRUE)
  
  # Fatal
  p_BUP_ODF_INJ_4wk  <- p_MET_ODF_INJ_4wk <- p_OD(rate = n_TX_OD, multiplier = n_TX_OD_mult, first_month = TRUE, fatal = TRUE, injection = TRUE)
  p_BUPC_ODF_INJ_4wk <- p_METC_ODF_INJ_4wk <- p_OD(rate = n_TXC_OD, multiplier = n_TXC_OD_mult, first_month = TRUE, fatal = TRUE, injection = TRUE)
  p_REL_ODF_INJ_4wk  <- p_OD(rate = n_REL_OD, multiplier = n_REL_OD_mult, first_month = TRUE, fatal = TRUE, injection = TRUE)
  p_ABS_ODF_INJ_4wk  <- p_OD(rate = n_ABS_OD, multiplier = n_ABS_OD_mult, first_month = TRUE, fatal = TRUE, injection = TRUE)

  #### Time-dependent survival probabilities ####
    # Empty 2-D matrix
    m_TDP <- array(0, dim = c(n_states, n_t),
                      dimnames = list(v_n_states, 1:n_t))

    # Probability of remaining in health state
    # All remain in fatal overdose, remain probability = 1
    for(i in 1:n_t){
      # Non-injection
        # Episode 1
        m_TDP[EP1 & BUP & NI, i]  <- as.vector(exp(p_frailty_BUP_NI_1 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI)))) # vector of remain probabilities
        m_TDP[EP1 & BUPC & NI, i] <- as.vector(exp(p_frailty_BUPC_NI_1 * p_weibull_scale_BUPC_NI * (((i - 1)^p_weibull_shape_BUPC_NI) - (i^p_weibull_shape_BUPC_NI))))
        m_TDP[EP1 & MET & NI, i]  <- as.vector(exp(p_frailty_MET_NI_1 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        m_TDP[EP1 & METC & NI, i] <- as.vector(exp(p_frailty_METC_NI_1 * p_weibull_scale_METC_NI * (((i - 1)^p_weibull_shape_METC_NI) - (i^p_weibull_shape_METC_NI))))
        m_TDP[EP1 & ABS & NI, i]  <- as.vector(exp(p_frailty_ABS_NI_1 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP1 & REL & NI, i]  <- as.vector(exp(p_frailty_REL_NI_1 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
        m_TDP[EP1 & ODN & NI, i]  <- p_ODN_ODN_NI #use same probability as overdose from relapse (no episode multipliers at this point)
        m_TDP[EP1 & ODF & NI, i]  <- 1 # all remain in ODF
        # Episode 2
        m_TDP[EP2 & BUP & NI, i]  <- as.vector(exp(p_frailty_BUP_NI_2 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI))))
        m_TDP[EP2 & BUPC & NI, i] <- as.vector(exp(p_frailty_BUPC_NI_2 * p_weibull_scale_BUPC_NI * (((i - 1)^p_weibull_shape_BUPC_NI) - (i^p_weibull_shape_BUPC_NI))))
        m_TDP[EP2 & MET & NI, i]  <- as.vector(exp(p_frailty_MET_NI_2 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        m_TDP[EP2 & METC & NI, i] <- as.vector(exp(p_frailty_METC_NI_2 * p_weibull_scale_METC_NI * (((i - 1)^p_weibull_shape_METC_NI) - (i^p_weibull_shape_METC_NI))))
        m_TDP[EP2 & ABS & NI, i]  <- as.vector(exp(p_frailty_ABS_NI_2 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP2 & REL & NI, i]  <- as.vector(exp(p_frailty_REL_NI_2 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
        m_TDP[EP2 & ODN & NI, i]  <- p_ODN_ODN_NI #use same probability as overdose from relapse (no episode multipliers at this point)
        m_TDP[EP2 & ODF & NI, i]  <- 1
        # Episode 3
        m_TDP[EP3 & BUP & NI, i] <- as.vector(exp(p_frailty_BUP_NI_3 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI))))
        m_TDP[EP3 & BUPC & NI, i] <- as.vector(exp(p_frailty_BUPC_NI_3 * p_weibull_scale_BUPC_NI * (((i - 1)^p_weibull_shape_BUPC_NI) - (i^p_weibull_shape_BUPC_NI))))
        m_TDP[EP3 & MET & NI, i] <- as.vector(exp(p_frailty_MET_NI_3 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        m_TDP[EP3 & METC & NI, i] <- as.vector(exp(p_frailty_METC_NI_3 * p_weibull_scale_METC_NI * (((i - 1)^p_weibull_shape_METC_NI) - (i^p_weibull_shape_METC_NI))))
        m_TDP[EP3 & ABS & NI, i] <- as.vector(exp(p_frailty_ABS_NI_3 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP3 & REL & NI, i] <- as.vector(exp(p_frailty_REL_NI_3 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
        m_TDP[EP3 & ODN & NI, i] <- p_ODN_ODN_NI #use same probability as overdose from relapse (no episode multipliers at this point)
        m_TDP[EP3 & ODF & NI, i] <- 1
        
    # Injection
        # Episode 1
        m_TDP[EP1 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_1 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ)))) # vector of remain probabilities 
        m_TDP[EP1 & BUPC & INJ, i] <- as.vector(exp(p_frailty_BUPC_INJ_1 * p_weibull_scale_BUPC_INJ * (((i - 1)^p_weibull_shape_BUPC_INJ) - (i^p_weibull_shape_BUPC_INJ))))
        m_TDP[EP1 & MET & INJ, i] <- as.vector(exp(p_frailty_MET_INJ_1 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        m_TDP[EP1 & METC & INJ, i] <- as.vector(exp(p_frailty_METC_INJ_1 * p_weibull_scale_METC_INJ * (((i - 1)^p_weibull_shape_METC_INJ) - (i^p_weibull_shape_METC_INJ))))
        m_TDP[EP1 & ABS & INJ, i] <- as.vector(exp(p_frailty_ABS_INJ_1 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP1 & REL & INJ, i] <- as.vector(exp(p_frailty_REL_INJ_1 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
        m_TDP[EP1 & ODN & INJ, i] <- p_ODN_ODN_INJ #use same probability as overdose from relapse (no episode multipliers at this point)
        m_TDP[EP1 & ODF & INJ, i] <- 1
        # Episode 2
        m_TDP[EP2 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_2 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ))))
        m_TDP[EP2 & BUPC & INJ, i] <- as.vector(exp(p_frailty_BUPC_INJ_2 * p_weibull_scale_BUPC_INJ * (((i - 1)^p_weibull_shape_BUPC_INJ) - (i^p_weibull_shape_BUPC_INJ))))
        m_TDP[EP2 & MET & INJ, i] <- as.vector(exp(p_frailty_MET_INJ_2 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        m_TDP[EP2 & METC & INJ, i] <- as.vector(exp(p_frailty_METC_INJ_2 * p_weibull_scale_METC_INJ * (((i - 1)^p_weibull_shape_METC_INJ) - (i^p_weibull_shape_METC_INJ))))
        m_TDP[EP2 & ABS & INJ, i] <- as.vector(exp(p_frailty_ABS_INJ_2 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP2 & REL & INJ, i] <- as.vector(exp(p_frailty_REL_INJ_2 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
        m_TDP[EP2 & ODN & INJ, i] <- p_ODN_ODN_INJ #use same probability as overdose from relapse (no episode multipliers at this point)
        m_TDP[EP2 & ODF & INJ, i] <- 1
        # Episode 3
        m_TDP[EP3 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_3 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ))))
        m_TDP[EP3 & BUPC & INJ, i] <- as.vector(exp(p_frailty_BUPC_INJ_3 * p_weibull_scale_BUPC_INJ * (((i - 1)^p_weibull_shape_BUPC_INJ) - (i^p_weibull_shape_BUPC_INJ))))
        m_TDP[EP3 & MET & INJ, i] <- as.vector(exp(p_frailty_MET_INJ_3 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        m_TDP[EP3 & METC & INJ, i] <- as.vector(exp(p_frailty_METC_INJ_3 * p_weibull_scale_METC_INJ * (((i - 1)^p_weibull_shape_METC_INJ) - (i^p_weibull_shape_METC_INJ))))
        m_TDP[EP3 & ABS & INJ, i] <- as.vector(exp(p_frailty_ABS_INJ_3 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP3 & REL & INJ, i] <- as.vector(exp(p_frailty_REL_INJ_3 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
        m_TDP[EP3 & ODN & INJ, i] <- p_ODN_ODN_INJ #use same probability as overdose from relapse (no episode multipliers at this point)
        m_TDP[EP3 & ODF & INJ, i] <- 1
  }

  # Probability of state-exit
  m_leave <- 1 - m_TDP
  
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

  # Alive probability in each period
  m_alive <- 1 - m_mort

  #### Unconditional transition probabilities ####
  # Empty 2-D unconditional transition matrix (from states, to states)
  # Create 2 matrices: 1 - First four weeks (higher OD prob)
  m_UP <- m_UP_4wk <- array(0, dim = c(n_states, n_states),
                            dimnames = list(v_n_states, v_n_states))
  # Populate unconditional transition matrix
  # Overdose probability populated first, accounting for higher probability of overdose transition in first 4 weeks of BUP/MET/REL
  # Non-Injection
  # From BUP
  # Overall
  m_UP[BUP & NI, BUPC & NI] <- p_BUP_BUPC_NI * (1 - p_BUP_ODN_NI - p_BUP_ODF_NI)
  m_UP[BUP & NI, MET & NI]  <- p_BUP_MET_NI * (1 - p_BUP_ODN_NI - p_BUP_ODF_NI)
  m_UP[BUP & NI, METC & NI] <- p_BUP_METC_NI * (1 - p_BUP_ODN_NI - p_BUP_ODF_NI)
  m_UP[BUP & NI, ABS & NI]  <- p_BUP_ABS_NI * (1 - p_BUP_ODN_NI - p_BUP_ODF_NI)
  m_UP[BUP & NI, REL & NI]  <- p_BUP_REL_NI * (1 - p_BUP_ODN_NI - p_BUP_ODF_NI)
  m_UP[BUP & NI, ODN & NI]  <- p_BUP_ODN_NI
  m_UP[BUP & NI, ODF & NI]  <- p_BUP_ODF_NI
  # First 4 weeks
  m_UP_4wk[BUP & NI, BUPC & NI] <- p_BUP_BUPC_NI * (1 - p_BUP_ODN_NI_4wk - p_BUP_ODF_NI_4wk)
  m_UP_4wk[BUP & NI, MET & NI]  <- p_BUP_MET_NI * (1 - p_BUP_ODN_NI_4wk - p_BUP_ODF_NI_4wk)
  m_UP_4wk[BUP & NI, METC & NI] <- p_BUP_METC_NI * (1 - p_BUP_ODN_NI_4wk - p_BUP_ODF_NI_4wk)
  m_UP_4wk[BUP & NI, ABS & NI]  <- p_BUP_ABS_NI * (1 - p_BUP_ODN_NI_4wk - p_BUP_ODF_NI_4wk)
  m_UP_4wk[BUP & NI, REL & NI]  <- p_BUP_REL_NI * (1 - p_BUP_ODN_NI_4wk - p_BUP_ODF_NI_4wk)
  m_UP_4wk[BUP & NI, ODN & NI]  <- p_BUP_ODN_NI_4wk
  m_UP_4wk[BUP & NI, ODF & NI]  <- p_BUP_ODF_NI_4wk
  
  # From BUPC
  # Overall
  m_UP[BUPC & NI, BUP & NI]  <- p_BUPC_BUP_NI * (1 - p_BUPC_ODN_NI - p_BUPC_ODF_NI)
  m_UP[BUPC & NI, MET & NI]  <- p_BUPC_MET_NI * (1 - p_BUPC_ODN_NI - p_BUPC_ODF_NI)
  m_UP[BUPC & NI, METC & NI] <- p_BUPC_METC_NI * (1 - p_BUPC_ODN_NI - p_BUPC_ODF_NI)
  m_UP[BUPC & NI, ABS & NI]  <- p_BUPC_ABS_NI * (1 - p_BUPC_ODN_NI - p_BUPC_ODF_NI)
  m_UP[BUPC & NI, REL & NI]  <- p_BUPC_REL_NI * (1 - p_BUPC_ODN_NI - p_BUPC_ODF_NI)
  m_UP[BUPC & NI, ODN & NI]  <- p_BUPC_ODN_NI
  m_UP[BUPC & NI, ODF & NI]  <- p_BUPC_ODF_NI
  # First 4 weeks
  m_UP_4wk[BUPC & NI, BUP & NI]  <- p_BUPC_BUP_NI * (1 - p_BUPC_ODN_NI_4wk - p_BUPC_ODF_NI_4wk)
  m_UP_4wk[BUPC & NI, MET & NI]  <- p_BUPC_MET_NI * (1 - p_BUPC_ODN_NI_4wk - p_BUPC_ODF_NI_4wk)
  m_UP_4wk[BUPC & NI, METC & NI] <- p_BUPC_METC_NI * (1 - p_BUPC_ODN_NI_4wk - p_BUPC_ODF_NI_4wk)
  m_UP_4wk[BUPC & NI, ABS & NI]  <- p_BUPC_ABS_NI * (1 - p_BUPC_ODN_NI_4wk - p_BUPC_ODF_NI_4wk)
  m_UP_4wk[BUPC & NI, REL & NI]  <- p_BUPC_REL_NI * (1 - p_BUPC_ODN_NI_4wk - p_BUPC_ODF_NI_4wk)
  m_UP_4wk[BUPC & NI, ODN & NI]  <- p_BUPC_ODN_NI_4wk
  m_UP_4wk[BUPC & NI, ODF & NI]  <- p_BUPC_ODF_NI_4wk
  
  # From MET
  # Overall
  m_UP[MET & NI, METC & NI] <- p_MET_METC_NI * (1 - p_MET_ODN_NI - p_MET_ODF_NI)
  m_UP[MET & NI, BUP & NI]  <- p_MET_BUP_NI * (1 - p_MET_ODN_NI - p_MET_ODF_NI)
  m_UP[MET & NI, BUPC & NI] <- p_MET_BUPC_NI * (1 - p_MET_ODN_NI - p_MET_ODF_NI)
  m_UP[MET & NI, ABS & NI]  <- p_MET_ABS_NI * (1 - p_MET_ODN_NI - p_MET_ODF_NI)
  m_UP[MET & NI, REL & NI]  <- p_MET_REL_NI * (1 - p_MET_ODN_NI - p_MET_ODF_NI)
  m_UP[MET & NI, ODN & NI]  <- p_MET_ODN_NI
  m_UP[MET & NI, ODF & NI]  <- p_MET_ODF_NI
  # First 4 weeks
  m_UP_4wk[MET & NI, METC & NI] <- p_MET_METC_NI * (1 - p_MET_ODN_NI_4wk - p_MET_ODF_NI_4wk)
  m_UP_4wk[MET & NI, BUP & NI]  <- p_MET_BUP_NI * (1 - p_MET_ODN_NI_4wk - p_MET_ODF_NI_4wk)
  m_UP_4wk[MET & NI, BUPC & NI] <- p_MET_BUPC_NI * (1 - p_MET_ODN_NI_4wk - p_MET_ODF_NI_4wk)
  m_UP_4wk[MET & NI, ABS & NI]  <- p_MET_ABS_NI * (1 - p_MET_ODN_NI_4wk - p_MET_ODF_NI_4wk)
  m_UP_4wk[MET & NI, REL & NI]  <- p_MET_REL_NI * (1 - p_MET_ODN_NI_4wk - p_MET_ODF_NI_4wk)
  m_UP_4wk[MET & NI, ODN & NI]  <- p_MET_ODN_NI_4wk
  m_UP_4wk[MET & NI, ODF & NI]  <- p_MET_ODF_NI_4wk
  
  # From METC
  # Overall
  m_UP[METC & NI, MET & NI]  <- p_METC_MET_NI * (1 - p_METC_ODN_NI - p_METC_ODF_NI)
  m_UP[METC & NI, BUP & NI]  <- p_METC_BUP_NI * (1 - p_METC_ODN_NI - p_METC_ODF_NI)
  m_UP[METC & NI, BUPC & NI] <- p_METC_BUPC_NI * (1 - p_METC_ODN_NI - p_METC_ODF_NI)
  m_UP[METC & NI, ABS & NI]  <- p_METC_ABS_NI * (1 - p_METC_ODN_NI - p_METC_ODF_NI)
  m_UP[METC & NI, REL & NI]  <- p_METC_REL_NI * (1 - p_METC_ODN_NI - p_METC_ODF_NI)
  m_UP[METC & NI, ODN & NI]  <- p_METC_ODN_NI
  m_UP[METC & NI, ODF & NI]  <- p_METC_ODF_NI
  # First 4 weeks
  m_UP_4wk[METC & NI, MET & NI]  <- p_METC_MET_NI * (1 - p_METC_ODN_NI_4wk - p_METC_ODF_NI_4wk)
  m_UP_4wk[METC & NI, BUP & NI]  <- p_METC_BUP_NI * (1 - p_METC_ODN_NI_4wk - p_METC_ODF_NI_4wk)
  m_UP_4wk[METC & NI, BUPC & NI] <- p_METC_BUPC_NI * (1 - p_METC_ODN_NI_4wk - p_METC_ODF_NI_4wk)
  m_UP_4wk[METC & NI, ABS & NI]  <- p_METC_ABS_NI * (1 - p_METC_ODN_NI_4wk - p_METC_ODF_NI_4wk)
  m_UP_4wk[METC & NI, REL & NI]  <- p_METC_REL_NI * (1 - p_METC_ODN_NI_4wk - p_METC_ODF_NI_4wk)
  m_UP_4wk[METC & NI, ODN & NI]  <- p_METC_ODN_NI_4wk
  m_UP_4wk[METC & NI, ODF & NI]  <- p_METC_ODF_NI_4wk
  
  # From ABS
  # Overall
  m_UP[ABS & NI, REL & NI]  <- p_ABS_REL_NI * (1 - p_ABS_ODN_NI - p_ABS_ODF_NI)
  m_UP[ABS & NI, MET & NI]  <- p_ABS_MET_NI * (1 - p_ABS_ODN_NI - p_ABS_ODF_NI)
  m_UP[ABS & NI, METC & NI] <- p_ABS_METC_NI * (1 - p_ABS_ODN_NI - p_ABS_ODF_NI)
  m_UP[ABS & NI, BUP & NI]  <- p_ABS_BUP_NI * (1 - p_ABS_ODN_NI - p_ABS_ODF_NI)
  m_UP[ABS & NI, BUPC & NI] <- p_ABS_BUPC_NI * (1 - p_ABS_ODN_NI - p_ABS_ODF_NI)
  m_UP[ABS & NI, ODN & NI]  <- p_ABS_ODN_NI
  m_UP[ABS & NI, ODF & NI]  <- p_ABS_ODF_NI
  # First 4 weeks
  m_UP_4wk[ABS & NI, REL & NI]  <- p_ABS_REL_NI * (1 - p_ABS_ODN_NI_4wk - p_ABS_ODF_NI_4wk)
  m_UP_4wk[ABS & NI, MET & NI]  <- p_ABS_MET_NI * (1 - p_ABS_ODN_NI_4wk - p_ABS_ODF_NI_4wk)
  m_UP_4wk[ABS & NI, METC & NI] <- p_ABS_METC_NI * (1 - p_ABS_ODN_NI_4wk - p_ABS_ODF_NI_4wk)
  m_UP_4wk[ABS & NI, BUP & NI]  <- p_ABS_BUP_NI * (1 - p_ABS_ODN_NI_4wk - p_ABS_ODF_NI_4wk)
  m_UP_4wk[ABS & NI, BUPC & NI] <- p_ABS_BUPC_NI * (1 - p_ABS_ODN_NI_4wk - p_ABS_ODF_NI_4wk)
  m_UP_4wk[ABS & NI, ODN & NI]  <- p_ABS_ODN_NI_4wk
  m_UP_4wk[ABS & NI, ODF & NI]  <- p_ABS_ODF_NI_4wk
  
  # From REL
  # Overall
  m_UP[REL & NI, MET & NI]  <- p_REL_MET_NI * (1 - p_REL_ODN_NI - p_REL_ODF_NI)
  m_UP[REL & NI, METC & NI] <- p_REL_METC_NI * (1 - p_REL_ODN_NI - p_REL_ODF_NI)
  m_UP[REL & NI, BUP & NI]  <- p_REL_BUP_NI * (1 - p_REL_ODN_NI - p_REL_ODF_NI)
  m_UP[REL & NI, BUPC & NI] <- p_REL_BUPC_NI * (1 - p_REL_ODN_NI - p_REL_ODF_NI)
  m_UP[REL & NI, ABS & NI]  <- p_REL_ABS_NI * (1 - p_REL_ODN_NI - p_REL_ODF_NI)
  m_UP[REL & NI, ODN & NI]  <- p_REL_ODN_NI
  m_UP[REL & NI, ODF & NI]  <- p_REL_ODF_NI
  # First 4 weeks
  m_UP_4wk[REL & NI, MET & NI]  <- p_REL_MET_NI * (1 - p_REL_ODN_NI_4wk - p_REL_ODF_NI_4wk)
  m_UP_4wk[REL & NI, METC & NI] <- p_REL_METC_NI * (1 - p_REL_ODN_NI_4wk - p_REL_ODF_NI_4wk)
  m_UP_4wk[REL & NI, BUP & NI]  <- p_REL_BUP_NI * (1 - p_REL_ODN_NI_4wk - p_REL_ODF_NI_4wk)
  m_UP_4wk[REL & NI, BUPC & NI] <- p_REL_BUPC_NI * (1 - p_REL_ODN_NI_4wk - p_REL_ODF_NI_4wk)
  m_UP_4wk[REL & NI, ABS & NI]  <- p_REL_ABS_NI * (1 - p_REL_ODN_NI_4wk - p_REL_ODF_NI_4wk)
  m_UP_4wk[REL & NI, ODN & NI]  <- p_REL_ODN_NI_4wk
  m_UP_4wk[REL & NI, ODF & NI]  <- p_REL_ODF_NI_4wk
  
  # From OD (first 4 weeks same)
  m_UP[ODN & NI, MET & NI]  <- m_UP_4wk[ODN & NI, MET & NI] <- p_ODN_MET_NI
  m_UP[ODN & NI, METC & NI] <- m_UP_4wk[ODN & NI, METC & NI] <- p_ODN_METC_NI
  m_UP[ODN & NI, BUP & NI]  <- m_UP_4wk[ODN & NI, BUP & NI] <- p_ODN_BUP_NI
  m_UP[ODN & NI, BUPC & NI] <- m_UP_4wk[ODN & NI, BUPC & NI] <- p_ODN_BUPC_NI
  m_UP[ODN & NI, ABS & NI]  <- m_UP_4wk[ODN & NI, ABS & NI] <- p_ODN_ABS_NI
  m_UP[ODN & NI, REL & NI]  <- m_UP_4wk[ODN & NI, REL & NI] <- p_ODN_REL_NI

  # Injection
  # From BUP
  # Overall
  m_UP[BUP & INJ, BUPC & INJ] <- p_BUP_BUPC_INJ * (1 - p_BUP_ODN_INJ - p_BUP_ODF_INJ)
  m_UP[BUP & INJ, MET & INJ]  <- p_BUP_MET_INJ * (1 - p_BUP_ODN_INJ - p_BUP_ODF_INJ)
  m_UP[BUP & INJ, METC & INJ] <- p_BUP_METC_INJ * (1 - p_BUP_ODN_INJ - p_BUP_ODF_INJ)
  m_UP[BUP & INJ, ABS & INJ]  <- p_BUP_ABS_INJ * (1 - p_BUP_ODN_INJ - p_BUP_ODF_INJ)
  m_UP[BUP & INJ, REL & INJ]  <- p_BUP_REL_INJ * (1 - p_BUP_ODN_INJ - p_BUP_ODF_INJ)
  m_UP[BUP & INJ, ODN & INJ]  <- p_BUP_ODN_INJ
  m_UP[BUP & INJ, ODF & INJ]  <- p_BUP_ODF_INJ
  # First 4 weeks
  m_UP_4wk[BUP & INJ, BUPC & INJ] <- p_BUP_BUPC_INJ * (1 - p_BUP_ODN_INJ_4wk - p_BUP_ODF_INJ_4wk)
  m_UP_4wk[BUP & INJ, MET & INJ]  <- p_BUP_MET_INJ * (1 - p_BUP_ODN_INJ_4wk - p_BUP_ODF_INJ_4wk)
  m_UP_4wk[BUP & INJ, METC & INJ] <- p_BUP_METC_INJ * (1 - p_BUP_ODN_INJ_4wk - p_BUP_ODF_INJ_4wk)
  m_UP_4wk[BUP & INJ, ABS & INJ]  <- p_BUP_ABS_INJ * (1 - p_BUP_ODN_INJ_4wk - p_BUP_ODF_INJ_4wk)
  m_UP_4wk[BUP & INJ, REL & INJ]  <- p_BUP_REL_INJ * (1 - p_BUP_ODN_INJ_4wk - p_BUP_ODF_INJ_4wk)
  m_UP_4wk[BUP & INJ, ODN & INJ]  <- p_BUP_ODN_INJ_4wk
  m_UP_4wk[BUP & INJ, ODF & INJ]  <- p_BUP_ODF_INJ_4wk
  
  # From BUPC
  # Overall
  m_UP[BUPC & INJ, BUP & INJ]  <- p_BUPC_BUP_INJ * (1 - p_BUPC_ODN_INJ - p_BUPC_ODF_INJ)
  m_UP[BUPC & INJ, MET & INJ]  <- p_BUPC_MET_INJ * (1 - p_BUPC_ODN_INJ - p_BUPC_ODF_INJ)
  m_UP[BUPC & INJ, METC & INJ] <- p_BUPC_METC_INJ * (1 - p_BUPC_ODN_INJ - p_BUPC_ODF_INJ)
  m_UP[BUPC & INJ, ABS & INJ]  <- p_BUPC_ABS_INJ * (1 - p_BUPC_ODN_INJ - p_BUPC_ODF_INJ)
  m_UP[BUPC & INJ, REL & INJ]  <- p_BUPC_REL_INJ * (1 - p_BUPC_ODN_INJ - p_BUPC_ODF_INJ)
  m_UP[BUPC & INJ, ODN & INJ]  <- p_BUPC_ODN_INJ
  m_UP[BUPC & INJ, ODF & INJ]  <- p_BUPC_ODF_INJ
  # First 4 weeks
  m_UP_4wk[BUPC & INJ, BUP & INJ]  <- p_BUPC_BUP_INJ * (1 - p_BUPC_ODN_INJ_4wk - p_BUPC_ODF_INJ_4wk)
  m_UP_4wk[BUPC & INJ, MET & INJ]  <- p_BUPC_MET_INJ * (1 - p_BUPC_ODN_INJ_4wk - p_BUPC_ODF_INJ_4wk)
  m_UP_4wk[BUPC & INJ, METC & INJ] <- p_BUPC_METC_INJ * (1 - p_BUPC_ODN_INJ_4wk - p_BUPC_ODF_INJ_4wk)
  m_UP_4wk[BUPC & INJ, ABS & INJ]  <- p_BUPC_ABS_INJ * (1 - p_BUPC_ODN_INJ_4wk - p_BUPC_ODF_INJ_4wk)
  m_UP_4wk[BUPC & INJ, REL & INJ]  <- p_BUPC_REL_INJ * (1 - p_BUPC_ODN_INJ_4wk - p_BUPC_ODF_INJ_4wk)
  m_UP_4wk[BUPC & INJ, ODN & INJ]  <- p_BUPC_ODN_INJ_4wk
  m_UP_4wk[BUPC & INJ, ODF & INJ]  <- p_BUPC_ODF_INJ_4wk
  
  # From MET
  # Overall
  m_UP[MET & INJ, METC & INJ] <- p_MET_METC_INJ * (1 - p_MET_ODN_INJ - p_MET_ODF_INJ)
  m_UP[MET & INJ, BUP & INJ]  <- p_MET_BUP_INJ * (1 - p_MET_ODN_INJ - p_MET_ODF_INJ)
  m_UP[MET & INJ, BUPC & INJ] <- p_MET_BUPC_INJ * (1 - p_MET_ODN_INJ - p_MET_ODF_INJ)
  m_UP[MET & INJ, ABS & INJ]  <- p_MET_ABS_INJ * (1 - p_MET_ODN_INJ - p_MET_ODF_INJ)
  m_UP[MET & INJ, REL & INJ]  <- p_MET_REL_INJ * (1 - p_MET_ODN_INJ - p_MET_ODF_INJ)
  m_UP[MET & INJ, ODN & INJ]  <- p_MET_ODN_INJ
  m_UP[MET & INJ, ODF & INJ]  <- p_MET_ODF_INJ
  # First 4 weeks
  m_UP_4wk[MET & INJ, METC & INJ] <- p_MET_METC_INJ * (1 - p_MET_ODN_INJ_4wk - p_MET_ODF_INJ_4wk)
  m_UP_4wk[MET & INJ, BUP & INJ]  <- p_MET_BUP_INJ * (1 - p_MET_ODN_INJ_4wk - p_MET_ODF_INJ_4wk)
  m_UP_4wk[MET & INJ, BUPC & INJ] <- p_MET_BUPC_INJ * (1 - p_MET_ODN_INJ_4wk - p_MET_ODF_INJ_4wk)
  m_UP_4wk[MET & INJ, ABS & INJ]  <- p_MET_ABS_INJ * (1 - p_MET_ODN_INJ_4wk - p_MET_ODF_INJ_4wk)
  m_UP_4wk[MET & INJ, REL & INJ]  <- p_MET_REL_INJ * (1 - p_MET_ODN_INJ_4wk - p_MET_ODF_INJ_4wk)
  m_UP_4wk[MET & INJ, ODN & INJ]  <- p_MET_ODN_INJ_4wk
  m_UP_4wk[MET & INJ, ODF & INJ]  <- p_MET_ODF_INJ_4wk
  
  # From METC
  # Overall
  m_UP[METC & INJ, MET & INJ]  <- p_METC_MET_INJ * (1 - p_METC_ODN_INJ - p_METC_ODF_INJ)
  m_UP[METC & INJ, BUP & INJ]  <- p_METC_BUP_INJ * (1 - p_METC_ODN_INJ - p_METC_ODF_INJ)
  m_UP[METC & INJ, BUPC & INJ] <- p_METC_BUPC_INJ * (1 - p_METC_ODN_INJ - p_METC_ODF_INJ)
  m_UP[METC & INJ, ABS & INJ]  <- p_METC_ABS_INJ * (1 - p_METC_ODN_INJ - p_METC_ODF_INJ)
  m_UP[METC & INJ, REL & INJ]  <- p_METC_REL_INJ * (1 - p_METC_ODN_INJ - p_METC_ODF_INJ)
  m_UP[METC & INJ, ODN & INJ]  <- p_METC_ODN_INJ
  m_UP[METC & INJ, ODF & INJ]  <- p_METC_ODF_INJ
  # First 4 weeks
  m_UP_4wk[METC & INJ, MET & INJ]  <- p_METC_MET_INJ * (1 - p_METC_ODN_INJ_4wk - p_METC_ODF_INJ_4wk)
  m_UP_4wk[METC & INJ, BUP & INJ]  <- p_METC_BUP_INJ * (1 - p_METC_ODN_INJ_4wk - p_METC_ODF_INJ_4wk)
  m_UP_4wk[METC & INJ, BUPC & INJ] <- p_METC_BUPC_INJ * (1 - p_METC_ODN_INJ_4wk - p_METC_ODF_INJ_4wk)
  m_UP_4wk[METC & INJ, ABS & INJ]  <- p_METC_ABS_INJ * (1 - p_METC_ODN_INJ_4wk - p_METC_ODF_INJ_4wk)
  m_UP_4wk[METC & INJ, REL & INJ]  <- p_METC_REL_INJ * (1 - p_METC_ODN_INJ_4wk - p_METC_ODF_INJ_4wk)
  m_UP_4wk[METC & INJ, ODN & INJ]  <- p_METC_ODN_INJ_4wk
  m_UP_4wk[METC & INJ, ODF & INJ]  <- p_METC_ODF_INJ_4wk
  
  # From ABS
  # Overall
  m_UP[ABS & INJ, REL & INJ]  <- p_ABS_REL_INJ * (1 - p_ABS_ODN_INJ - p_ABS_ODF_INJ)
  m_UP[ABS & INJ, MET & INJ]  <- p_ABS_MET_INJ * (1 - p_ABS_ODN_INJ - p_ABS_ODF_INJ)
  m_UP[ABS & INJ, METC & INJ] <- p_ABS_METC_INJ * (1 - p_ABS_ODN_INJ - p_ABS_ODF_INJ)
  m_UP[ABS & INJ, BUP & INJ]  <- p_ABS_BUP_INJ * (1 - p_ABS_ODN_INJ - p_ABS_ODF_INJ)
  m_UP[ABS & INJ, BUPC & INJ] <- p_ABS_BUPC_INJ * (1 - p_ABS_ODN_INJ - p_ABS_ODF_INJ)
  m_UP[ABS & INJ, ODN & INJ]  <- p_ABS_ODN_INJ
  m_UP[ABS & INJ, ODF & INJ]  <- p_ABS_ODF_INJ
  # First 4 weeks
  m_UP_4wk[ABS & INJ, REL & INJ]  <- p_ABS_REL_INJ * (1 - p_ABS_ODN_INJ_4wk - p_ABS_ODF_INJ_4wk)
  m_UP_4wk[ABS & INJ, MET & INJ]  <- p_ABS_MET_INJ * (1 - p_ABS_ODN_INJ_4wk - p_ABS_ODF_INJ_4wk)
  m_UP_4wk[ABS & INJ, METC & INJ] <- p_ABS_METC_INJ * (1 - p_ABS_ODN_INJ_4wk - p_ABS_ODF_INJ_4wk)
  m_UP_4wk[ABS & INJ, BUP & INJ]  <- p_ABS_BUP_INJ * (1 - p_ABS_ODN_INJ_4wk - p_ABS_ODF_INJ_4wk)
  m_UP_4wk[ABS & INJ, BUPC & INJ] <- p_ABS_BUPC_INJ * (1 - p_ABS_ODN_INJ_4wk - p_ABS_ODF_INJ_4wk)
  m_UP_4wk[ABS & INJ, ODN & INJ]  <- p_ABS_ODN_INJ_4wk
  m_UP_4wk[ABS & INJ, ODF & INJ]  <- p_ABS_ODF_INJ_4wk
  
  # From REL
  # Overall
  m_UP[REL & INJ, MET & INJ]  <- p_REL_MET_INJ * (1 - p_REL_ODN_INJ - p_REL_ODF_INJ)
  m_UP[REL & INJ, METC & INJ] <- p_REL_METC_INJ * (1 - p_REL_ODN_INJ - p_REL_ODF_INJ)
  m_UP[REL & INJ, BUP & INJ]  <- p_REL_BUP_INJ * (1 - p_REL_ODN_INJ - p_REL_ODF_INJ)
  m_UP[REL & INJ, BUPC & INJ] <- p_REL_BUPC_INJ * (1 - p_REL_ODN_INJ - p_REL_ODF_INJ)
  m_UP[REL & INJ, ABS & INJ]  <- p_REL_ABS_INJ * (1 - p_REL_ODN_INJ - p_REL_ODF_INJ)
  m_UP[REL & INJ, ODN & INJ]  <- p_REL_ODN_INJ
  m_UP[REL & INJ, ODF & INJ]  <- p_REL_ODF_INJ
  # First 4 weeks
  m_UP_4wk[REL & INJ, MET & INJ]  <- p_REL_MET_INJ * (1 - p_REL_ODN_INJ_4wk - p_REL_ODF_INJ_4wk)
  m_UP_4wk[REL & INJ, METC & INJ] <- p_REL_METC_INJ * (1 - p_REL_ODN_INJ_4wk - p_REL_ODF_INJ_4wk)
  m_UP_4wk[REL & INJ, BUP & INJ]  <- p_REL_BUP_INJ * (1 - p_REL_ODN_INJ_4wk - p_REL_ODF_INJ_4wk)
  m_UP_4wk[REL & INJ, BUPC & INJ] <- p_REL_BUPC_INJ * (1 - p_REL_ODN_INJ_4wk - p_REL_ODF_INJ_4wk)
  m_UP_4wk[REL & INJ, ABS & INJ]  <- p_REL_ABS_INJ * (1 - p_REL_ODN_INJ_4wk - p_REL_ODF_INJ_4wk)
  m_UP_4wk[REL & INJ, ODN & INJ]  <- p_REL_ODN_INJ_4wk
  m_UP_4wk[REL & INJ, ODF & INJ]  <- p_REL_ODF_INJ_4wk
  
  # From OD (first 4 wks same)
  m_UP[ODN & INJ, MET & INJ]  <- m_UP_4wk[ODN & INJ, MET & INJ] <- p_ODN_MET_INJ
  m_UP[ODN & INJ, METC & INJ] <- m_UP_4wk[ODN & INJ, METC & INJ] <- p_ODN_METC_INJ
  m_UP[ODN & INJ, BUP & INJ]  <- m_UP_4wk[ODN & INJ, BUP & INJ] <- p_ODN_BUP_INJ
  m_UP[ODN & INJ, BUPC & INJ] <- m_UP_4wk[ODN & INJ, BUPC & INJ] <- p_ODN_BUPC_INJ
  m_UP[ODN & INJ, ABS & INJ]  <- m_UP_4wk[ODN & INJ, ABS & INJ] <- p_ODN_ABS_INJ
  m_UP[ODN & INJ, REL & INJ]  <- m_UP_4wk[ODN & INJ, REL & INJ] <- p_ODN_REL_INJ

  #### Create full time-dependent transition array ####
  # Empty 3-D array
  a_TDP <- array(0, dim = c(n_states, n_states, n_t),
                   dimnames = list(v_n_states, v_n_states, 1:n_t))
  
  # Add transitions conditional on state-exit (m_leave = 1 - remain)
  # Modified transitions for first month
  for (i in 1:4){
    a_TDP[, , i] <- m_UP_4wk * m_leave[, i]
  }
  # All transitions 1+ months
  for (i in 5:n_t){
    a_TDP[, , i] <- m_UP * m_leave[, i]
  }

  # Add time-dependent remain probabilities
  for (i in 1:n_t){
    # Non-injection
    # Episode 1
    a_TDP[EP1 & BUP & NI, EP1 & BUP & NI, i]   <- m_TDP[EP1 & BUP & NI, i]
    a_TDP[EP1 & BUPC & NI, EP1 & BUPC & NI, i] <- m_TDP[EP1 & BUPC & NI, i]
    a_TDP[EP1 & MET & NI, EP1 & MET & NI, i]   <- m_TDP[EP1 & MET & NI, i]
    a_TDP[EP1 & METC & NI, EP1 & METC & NI, i] <- m_TDP[EP1 & METC & NI, i]
    a_TDP[EP1 & ABS & NI, EP1 & ABS & NI, i]   <- m_TDP[EP1 & ABS & NI, i]
    a_TDP[EP1 & REL & NI, EP1 & REL & NI, i]   <- m_TDP[EP1 & REL & NI, i]
    a_TDP[EP1 & ODN & NI, EP1 & ODN & NI, i]   <- m_TDP[EP1 & ODN & NI, i]
    a_TDP[EP1 & ODF & NI, EP1 & ODF & NI, i]   <- m_TDP[EP1 & ODF & NI, i]
    # Episode 2
    a_TDP[EP2 & BUP & NI, EP2 & BUP & NI, i]   <- m_TDP[EP2 & BUP & NI, i]
    a_TDP[EP2 & BUPC & NI, EP2 & BUPC & NI, i] <- m_TDP[EP2 & BUPC & NI, i]
    a_TDP[EP2 & MET & NI, EP2 & MET & NI, i]   <- m_TDP[EP2 & MET & NI, i]
    a_TDP[EP2 & METC & NI, EP2 & METC & NI, i] <- m_TDP[EP2 & METC & NI, i]
    a_TDP[EP2 & ABS & NI, EP2 & ABS & NI, i]   <- m_TDP[EP2 & ABS & NI, i]
    a_TDP[EP2 & REL & NI, EP2 & REL & NI, i]   <- m_TDP[EP2 & REL & NI, i]
    a_TDP[EP2 & ODN & NI, EP2 & ODN & NI, i]   <- m_TDP[EP2 & ODN & NI, i]
    a_TDP[EP2 & ODF & NI, EP2 & ODF & NI, i]   <- m_TDP[EP2 & ODF & NI, i]
    # Episode 3
    a_TDP[EP3 & BUP & NI, EP3 & BUP & NI, i]   <- m_TDP[EP3 & BUP & NI, i]
    a_TDP[EP3 & BUPC & NI, EP3 & BUPC & NI, i] <- m_TDP[EP3 & BUPC & NI, i]
    a_TDP[EP3 & MET & NI, EP3 & MET & NI, i]   <- m_TDP[EP3 & MET & NI, i]
    a_TDP[EP3 & METC & NI, EP3 & METC & NI, i] <- m_TDP[EP3 & METC & NI, i]
    a_TDP[EP3 & ABS & NI, EP3 & ABS & NI, i]   <- m_TDP[EP3 & ABS & NI, i]
    a_TDP[EP3 & REL & NI, EP3 & REL & NI, i]   <- m_TDP[EP3 & REL & NI, i]
    a_TDP[EP3 & ODN & NI, EP3 & ODN & NI, i]   <- m_TDP[EP3 & ODN & NI, i]
    a_TDP[EP3 & ODF & NI, EP3 & ODF & NI, i]   <- m_TDP[EP3 & ODF & NI, i]
    
    # Injection
    # Episode 1
    a_TDP[EP1 & BUP & INJ, EP1 & BUP & INJ, i]   <- m_TDP[EP1 & BUP & INJ, i]
    a_TDP[EP1 & BUPC & INJ, EP1 & BUPC & INJ, i] <- m_TDP[EP1 & BUPC & INJ, i]
    a_TDP[EP1 & MET & INJ, EP1 & MET & INJ, i]   <- m_TDP[EP1 & MET & INJ, i]
    a_TDP[EP1 & METC & INJ, EP1 & METC & INJ, i] <- m_TDP[EP1 & METC & INJ, i]
    a_TDP[EP1 & ABS & INJ, EP1 & ABS & INJ, i]   <- m_TDP[EP1 & ABS & INJ, i]
    a_TDP[EP1 & REL & INJ, EP1 & REL & INJ, i]   <- m_TDP[EP1 & REL & INJ, i]
    a_TDP[EP1 & ODN & INJ, EP1 & ODN & INJ, i]   <- m_TDP[EP1 & ODN & INJ, i]
    a_TDP[EP1 & ODF & INJ, EP1 & ODF & INJ, i]   <- m_TDP[EP1 & ODF & INJ, i]
    # Episode 2
    a_TDP[EP2 & BUP & INJ, EP2 & BUP & INJ, i]   <- m_TDP[EP2 & BUP & INJ, i]
    a_TDP[EP2 & BUPC & INJ, EP2 & BUPC & INJ, i] <- m_TDP[EP2 & BUPC & INJ, i]
    a_TDP[EP2 & MET & INJ, EP2 & MET & INJ, i]   <- m_TDP[EP2 & MET & INJ, i]
    a_TDP[EP2 & METC & INJ, EP2 & METC & INJ, i] <- m_TDP[EP2 & METC & INJ, i]
    a_TDP[EP2 & ABS & INJ, EP2 & ABS & INJ, i]   <- m_TDP[EP2 & ABS & INJ, i]
    a_TDP[EP2 & REL & INJ, EP2 & REL & INJ, i]   <- m_TDP[EP2 & REL & INJ, i]
    a_TDP[EP2 & ODN & INJ, EP2 & ODN & INJ, i]   <- m_TDP[EP2 & ODN & INJ, i]
    a_TDP[EP2 & ODF & INJ, EP2 & ODF & INJ, i]   <- m_TDP[EP2 & ODF & INJ, i]
    # Episode 3
    a_TDP[EP3 & BUP & INJ, EP3 & BUP & INJ, i]   <- m_TDP[EP3 & BUP & INJ, i]
    a_TDP[EP3 & BUPC & INJ, EP3 & BUPC & INJ, i] <- m_TDP[EP3 & BUPC & INJ, i]
    a_TDP[EP3 & MET & INJ, EP3 & MET & INJ, i]   <- m_TDP[EP3 & MET & INJ, i]
    a_TDP[EP3 & METC & INJ, EP3 & METC & INJ, i] <- m_TDP[EP3 & METC & INJ, i]
    a_TDP[EP3 & ABS & INJ, EP3 & ABS & INJ, i]   <- m_TDP[EP3 & ABS & INJ, i]
    a_TDP[EP3 & REL & INJ, EP3 & REL & INJ, i]   <- m_TDP[EP3 & REL & INJ, i]
    a_TDP[EP3 & ODN & INJ, EP3 & ODN & INJ, i]   <- m_TDP[EP3 & ODN & INJ, i]
    a_TDP[EP3 & ODF & INJ, EP3 & ODF & INJ, i]   <- m_TDP[EP3 & ODF & INJ, i]
  }
  
  initial <- a_TDP[, , 1]

  #### Seroconversion ####
  # Apply seroconversion probability to re-weight NEG -> POS for to-states each time period
  # Probabilities applied equally across POS/NEG initially, re-weight by sero prob
  # Non-injection
  # BUP
  # From NEG
  a_TDP[NEG & NI, BUP & NI & NEG, ]  <- a_TDP[NEG & NI, BUP & NI & NEG, ] * (1 - p_HIV_BUP_NI - p_HCV_BUP_NI)
  a_TDP[NEG & NI, BUP & NI & HIV, ]  <- a_TDP[NEG & NI, BUP & NI & HIV, ] * p_HIV_BUP_NI
  a_TDP[NEG & NI, BUP & NI & HCV, ]  <- a_TDP[NEG & NI, BUP & NI & HCV, ] * p_HCV_BUP_NI
  # From HIV
  a_TDP[HIV & NI, BUP & NI & HIV, ]  <- a_TDP[HIV & NI, BUP & NI & HIV, ] * (1 - p_HIV_HCV_BUP_NI)
  a_TDP[HIV & NI, BUP & NI & COI, ]  <- a_TDP[HIV & NI, BUP & NI & COI, ] * p_HIV_HCV_BUP_NI # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & NI, BUP & NI & HCV, ]  <- a_TDP[HCV & NI, BUP & NI & HCV, ] * (1 - p_HCV_HIV_BUP_NI)
  a_TDP[HCV & NI, BUP & NI & COI, ]  <- a_TDP[HCV & NI, BUP & NI & COI, ] * p_HCV_HIV_BUP_NI # Probability of HIV conditional on HCV
  
  # BUPC
  # From NEG
  a_TDP[NEG & NI, BUPC & NI & NEG, ]  <- a_TDP[NEG & NI, BUPC & NI & NEG, ] * (1 - p_HIV_BUPC_NI - p_HCV_BUPC_NI)
  a_TDP[NEG & NI, BUPC & NI & HIV, ]  <- a_TDP[NEG & NI, BUPC & NI & HIV, ] * p_HIV_BUPC_NI
  a_TDP[NEG & NI, BUPC & NI & HCV, ]  <- a_TDP[NEG & NI, BUPC & NI & HCV, ] * p_HCV_BUPC_NI
  # From HIV
  a_TDP[HIV & NI, BUPC & NI & HIV, ]  <- a_TDP[HIV & NI, BUPC & NI & HIV, ] * (1 - p_HIV_HCV_BUPC_NI)
  a_TDP[HIV & NI, BUPC & NI & COI, ]  <- a_TDP[HIV & NI, BUPC & NI & COI, ] * p_HIV_HCV_BUPC_NI # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & NI, BUPC & NI & HCV, ]  <- a_TDP[HCV & NI, BUPC & NI & HCV, ] * (1 - p_HCV_HIV_BUPC_NI)
  a_TDP[HCV & NI, BUPC & NI & COI, ]  <- a_TDP[HCV & NI, BUPC & NI & COI, ] * p_HCV_HIV_BUPC_NI # Probability of HIV conditional on HCV
  
  # MET
  # From NEG
  a_TDP[NEG & NI, MET & NI & NEG, ]  <- a_TDP[NEG & NI, MET & NI & NEG, ] * (1 - p_HIV_MET_NI - p_HCV_MET_NI)
  a_TDP[NEG & NI, MET & NI & HIV, ]  <- a_TDP[NEG & NI, MET & NI & HIV, ] * p_HIV_MET_NI
  a_TDP[NEG & NI, MET & NI & HCV, ]  <- a_TDP[NEG & NI, MET & NI & HCV, ] * p_HCV_MET_NI
  # From HIV
  a_TDP[HIV & NI, MET & NI & HIV, ]  <- a_TDP[HIV & NI, MET & NI & HIV, ] * (1 - p_HIV_HCV_MET_NI)
  a_TDP[HIV & NI, MET & NI & COI, ]  <- a_TDP[HIV & NI, MET & NI & COI, ] * p_HIV_HCV_MET_NI # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & NI, MET & NI & HCV, ]  <- a_TDP[HCV & NI, MET & NI & HCV, ] * (1 - p_HCV_HIV_MET_NI)
  a_TDP[HCV & NI, MET & NI & COI, ]  <- a_TDP[HCV & NI, MET & NI & COI, ] * p_HCV_HIV_MET_NI # Probability of HIV conditional on HCV
  
  # METC
  # From NEG
  a_TDP[NEG & NI, METC & NI & NEG, ]  <- a_TDP[NEG & NI, METC & NI & NEG, ] * (1 - p_HIV_METC_NI - p_HCV_METC_NI)
  a_TDP[NEG & NI, METC & NI & HIV, ]  <- a_TDP[NEG & NI, METC & NI & HIV, ] * p_HIV_METC_NI
  a_TDP[NEG & NI, METC & NI & HCV, ]  <- a_TDP[NEG & NI, METC & NI & HCV, ] * p_HCV_METC_NI
  # From HIV
  a_TDP[HIV & NI, METC & NI & HIV, ]  <- a_TDP[HIV & NI, METC & NI & HIV, ] * (1 - p_HIV_HCV_METC_NI)
  a_TDP[HIV & NI, METC & NI & COI, ]  <- a_TDP[HIV & NI, METC & NI & COI, ] * p_HIV_HCV_METC_NI # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & NI, METC & NI & HCV, ]  <- a_TDP[HCV & NI, METC & NI & HCV, ] * (1 - p_HCV_HIV_METC_NI)
  a_TDP[HCV & NI, METC & NI & COI, ]  <- a_TDP[HCV & NI, METC & NI & COI, ] * p_HCV_HIV_METC_NI # Probability of HIV conditional on HCV
  
  # REL
  # From NEG
  a_TDP[NEG & NI, REL & NI & NEG, ]  <- a_TDP[NEG & NI, REL & NI & NEG, ] * (1 - p_HIV_REL_NI - p_HCV_REL_NI)
  a_TDP[NEG & NI, REL & NI & HIV, ]  <- a_TDP[NEG & NI, REL & NI & HIV, ] * p_HIV_REL_NI
  a_TDP[NEG & NI, REL & NI & HCV, ]  <- a_TDP[NEG & NI, REL & NI & HCV, ] * p_HCV_REL_NI
  # From HIV
  a_TDP[HIV & NI, REL & NI & HIV, ]  <- a_TDP[HIV & NI, REL & NI & HIV, ] * (1 - p_HIV_HCV_REL_NI)
  a_TDP[HIV & NI, REL & NI & COI, ]  <- a_TDP[HIV & NI, REL & NI & COI, ] * p_HIV_HCV_REL_NI # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & NI, REL & NI & HCV, ]  <- a_TDP[HCV & NI, REL & NI & HCV, ] * (1 - p_HCV_HIV_REL_NI)
  a_TDP[HCV & NI, REL & NI & COI, ]  <- a_TDP[HCV & NI, REL & NI & COI, ] * p_HCV_HIV_REL_NI # Probability of HIV conditional on HCV
  
  # ODN
  # From NEG
  a_TDP[NEG & NI, ODN & NI & NEG, ]  <- a_TDP[NEG & NI, ODN & NI & NEG, ] * (1 - p_HIV_ODN_NI - p_HCV_ODN_NI)
  a_TDP[NEG & NI, ODN & NI & HIV, ]  <- a_TDP[NEG & NI, ODN & NI & HIV, ] * p_HIV_ODN_NI
  a_TDP[NEG & NI, ODN & NI & HCV, ]  <- a_TDP[NEG & NI, ODN & NI & HCV, ] * p_HCV_ODN_NI
  # From HIV
  a_TDP[HIV & NI, ODN & NI & HIV, ]  <- a_TDP[HIV & NI, ODN & NI & HIV, ] * (1 - p_HIV_HCV_ODN_NI)
  a_TDP[HIV & NI, ODN & NI & COI, ]  <- a_TDP[HIV & NI, ODN & NI & COI, ] * p_HIV_HCV_ODN_NI # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & NI, ODN & NI & HCV, ]  <- a_TDP[HCV & NI, ODN & NI & HCV, ] * (1 - p_HCV_HIV_ODN_NI)
  a_TDP[HCV & NI, ODN & NI & COI, ]  <- a_TDP[HCV & NI, ODN & NI & COI, ] * p_HCV_HIV_ODN_NI # Probability of HIV conditional on HCV

  # ODF
  # From NEG
  a_TDP[NEG & NI, ODF & NI & NEG, ]  <- a_TDP[NEG & NI, ODF & NI & NEG, ] * 1
  a_TDP[NEG & NI, ODF & NI & HIV, ]  <- a_TDP[NEG & NI, ODF & NI & HIV, ] * 0
  a_TDP[NEG & NI, ODF & NI & HCV, ]  <- a_TDP[NEG & NI, ODF & NI & HCV, ] * 0
  # From HIV
  a_TDP[HIV & NI, ODF & NI & HIV, ]  <- a_TDP[HIV & NI, ODF & NI & HIV, ] * 1
  a_TDP[HIV & NI, ODF & NI & COI, ]  <- a_TDP[HIV & NI, ODF & NI & COI, ] * 0
  # From HCV
  a_TDP[HCV & NI, ODF & NI & HCV, ]  <- a_TDP[HCV & NI, ODF & NI & HCV, ] * 1
  a_TDP[HCV & NI, ODF & NI & COI, ]  <- a_TDP[HCV & NI, ODF & NI & COI, ] * 0
  
  
  # ABS
  # From NEG
  a_TDP[NEG & NI, ABS & NI & NEG, ]  <- a_TDP[NEG & NI, ABS & NI & NEG, ] * (1 - p_HIV_ABS_NI - p_HCV_ABS_NI)
  a_TDP[NEG & NI, ABS & NI & HIV, ]  <- a_TDP[NEG & NI, ABS & NI & HIV, ] * p_HIV_ABS_NI
  a_TDP[NEG & NI, ABS & NI & HCV, ]  <- a_TDP[NEG & NI, ABS & NI & HCV, ] * p_HCV_ABS_NI
  # From HIV
  a_TDP[HIV & NI, ABS & NI & HIV, ]  <- a_TDP[HIV & NI, ABS & NI & HIV, ] * (1 - p_HIV_HCV_ABS_NI)
  a_TDP[HIV & NI, ABS & NI & COI, ]  <- a_TDP[HIV & NI, ABS & NI & COI, ] * p_HIV_HCV_ABS_NI # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & NI, ABS & NI & HCV, ]  <- a_TDP[HCV & NI, ABS & NI & HCV, ] * (1 - p_HCV_HIV_ABS_NI)
  a_TDP[HCV & NI, ABS & NI & COI, ]  <- a_TDP[HCV & NI, ABS & NI & COI, ] * p_HCV_HIV_ABS_NI # Probability of HIV conditional on HCV

  
  # Injection
  # BUP
  # From NEG
  a_TDP[NEG & INJ, BUP & INJ & NEG, ]  <- a_TDP[NEG & INJ, BUP & INJ & NEG, ] * (1 - p_HIV_BUP_INJ - p_HCV_BUP_INJ)
  a_TDP[NEG & INJ, BUP & INJ & HIV, ]  <- a_TDP[NEG & INJ, BUP & INJ & HIV, ] * p_HIV_BUP_INJ
  a_TDP[NEG & INJ, BUP & INJ & HCV, ]  <- a_TDP[NEG & INJ, BUP & INJ & HCV, ] * p_HCV_BUP_INJ
  # From HIV
  a_TDP[HIV & INJ, BUP & INJ & HIV, ]  <- a_TDP[HIV & INJ, BUP & INJ & HIV, ] * (1 - p_HIV_HCV_BUP_INJ)
  a_TDP[HIV & INJ, BUP & INJ & COI, ]  <- a_TDP[HIV & INJ, BUP & INJ & COI, ] * p_HIV_HCV_BUP_INJ # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & INJ, BUP & INJ & HCV, ]  <- a_TDP[HCV & INJ, BUP & INJ & HCV, ] * (1 - p_HCV_HIV_BUP_INJ)
  a_TDP[HCV & INJ, BUP & INJ & COI, ]  <- a_TDP[HCV & INJ, BUP & INJ & COI, ] * p_HCV_HIV_BUP_INJ # Probability of HIV conditional on HCV
  
  # BUPC
  # From NEG
  a_TDP[NEG & INJ, BUPC & INJ & NEG, ]  <- a_TDP[NEG & INJ, BUPC & INJ & NEG, ] * (1 - p_HIV_BUPC_INJ - p_HCV_BUPC_INJ)
  a_TDP[NEG & INJ, BUPC & INJ & HIV, ]  <- a_TDP[NEG & INJ, BUPC & INJ & HIV, ] * p_HIV_BUPC_INJ
  a_TDP[NEG & INJ, BUPC & INJ & HCV, ]  <- a_TDP[NEG & INJ, BUPC & INJ & HCV, ] * p_HCV_BUPC_INJ
  # From HIV
  a_TDP[HIV & INJ, BUPC & INJ & HIV, ]  <- a_TDP[HIV & INJ, BUPC & INJ & HIV, ] * (1 - p_HIV_HCV_BUPC_INJ)
  a_TDP[HIV & INJ, BUPC & INJ & COI, ]  <- a_TDP[HIV & INJ, BUPC & INJ & COI, ] * p_HIV_HCV_BUPC_INJ # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & INJ, BUPC & INJ & HCV, ]  <- a_TDP[HCV & INJ, BUPC & INJ & HCV, ] * (1 - p_HCV_HIV_BUPC_INJ)
  a_TDP[HCV & INJ, BUPC & INJ & COI, ]  <- a_TDP[HCV & INJ, BUPC & INJ & COI, ] * p_HCV_HIV_BUPC_INJ # Probability of HIV conditional on HCV
  
  # MET
  # From NEG
  a_TDP[NEG & INJ, MET & INJ & NEG, ]  <- a_TDP[NEG & INJ, MET & INJ & NEG, ] * (1 - p_HIV_MET_INJ - p_HCV_MET_INJ)
  a_TDP[NEG & INJ, MET & INJ & HIV, ]  <- a_TDP[NEG & INJ, MET & INJ & HIV, ] * p_HIV_MET_INJ
  a_TDP[NEG & INJ, MET & INJ & HCV, ]  <- a_TDP[NEG & INJ, MET & INJ & HCV, ] * p_HCV_MET_INJ
  # From HIV
  a_TDP[HIV & INJ, MET & INJ & HIV, ]  <- a_TDP[HIV & INJ, MET & INJ & HIV, ] * (1 - p_HIV_HCV_MET_INJ)
  a_TDP[HIV & INJ, MET & INJ & COI, ]  <- a_TDP[HIV & INJ, MET & INJ & COI, ] * p_HIV_HCV_MET_INJ # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & INJ, MET & INJ & HCV, ]  <- a_TDP[HCV & INJ, MET & INJ & HCV, ] * (1 - p_HCV_HIV_MET_INJ)
  a_TDP[HCV & INJ, MET & INJ & COI, ]  <- a_TDP[HCV & INJ, MET & INJ & COI, ] * p_HCV_HIV_MET_INJ # Probability of HIV conditional on HCV
  
  # METC
  # From NEG
  a_TDP[NEG & INJ, METC & INJ & NEG, ]  <- a_TDP[NEG & INJ, METC & INJ & NEG, ] * (1 - p_HIV_METC_INJ - p_HCV_METC_INJ)
  a_TDP[NEG & INJ, METC & INJ & HIV, ]  <- a_TDP[NEG & INJ, METC & INJ & HIV, ] * p_HIV_METC_INJ
  a_TDP[NEG & INJ, METC & INJ & HCV, ]  <- a_TDP[NEG & INJ, METC & INJ & HCV, ] * p_HCV_METC_INJ
  # From HIV
  a_TDP[HIV & INJ, METC & INJ & HIV, ]  <- a_TDP[HIV & INJ, METC & INJ & HIV, ] * (1 - p_HIV_HCV_METC_INJ)
  a_TDP[HIV & INJ, METC & INJ & COI, ]  <- a_TDP[HIV & INJ, METC & INJ & COI, ] * p_HIV_HCV_METC_INJ # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & INJ, METC & INJ & HCV, ]  <- a_TDP[HCV & INJ, METC & INJ & HCV, ] * (1 - p_HCV_HIV_METC_INJ)
  a_TDP[HCV & INJ, METC & INJ & COI, ]  <- a_TDP[HCV & INJ, METC & INJ & COI, ] * p_HCV_HIV_METC_INJ # Probability of HIV conditional on HCV
  
  # REL
  # From NEG
  a_TDP[NEG & INJ, REL & INJ & NEG, ]  <- a_TDP[NEG & INJ, REL & INJ & NEG, ] * (1 - p_HIV_REL_INJ - p_HCV_REL_INJ)
  a_TDP[NEG & INJ, REL & INJ & HIV, ]  <- a_TDP[NEG & INJ, REL & INJ & HIV, ] * p_HIV_REL_INJ
  a_TDP[NEG & INJ, REL & INJ & HCV, ]  <- a_TDP[NEG & INJ, REL & INJ & HCV, ] * p_HCV_REL_INJ
  # From HIV
  a_TDP[HIV & INJ, REL & INJ & HIV, ]  <- a_TDP[HIV & INJ, REL & INJ & HIV, ] * (1 - p_HIV_HCV_REL_INJ)
  a_TDP[HIV & INJ, REL & INJ & COI, ]  <- a_TDP[HIV & INJ, REL & INJ & COI, ] * p_HIV_HCV_REL_INJ # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & INJ, REL & INJ & HCV, ]  <- a_TDP[HCV & INJ, REL & INJ & HCV, ] * (1 - p_HCV_HIV_REL_INJ)
  a_TDP[HCV & INJ, REL & INJ & COI, ]  <- a_TDP[HCV & INJ, REL & INJ & COI, ] * p_HCV_HIV_REL_INJ # Probability of HIV conditional on HCV
  
  # ODN
  # From NEG
  a_TDP[NEG & INJ, ODN & INJ & NEG, ]  <- a_TDP[NEG & INJ, ODN & INJ & NEG, ] * (1 - p_HIV_ODN_INJ - p_HCV_ODN_INJ)
  a_TDP[NEG & INJ, ODN & INJ & HIV, ]  <- a_TDP[NEG & INJ, ODN & INJ & HIV, ] * p_HIV_ODN_INJ
  a_TDP[NEG & INJ, ODN & INJ & HCV, ]  <- a_TDP[NEG & INJ, ODN & INJ & HCV, ] * p_HCV_ODN_INJ
  # From HIV
  a_TDP[HIV & INJ, ODN & INJ & HIV, ]  <- a_TDP[HIV & INJ, ODN & INJ & HIV, ] * (1 - p_HIV_HCV_ODN_INJ)
  a_TDP[HIV & INJ, ODN & INJ & COI, ]  <- a_TDP[HIV & INJ, ODN & INJ & COI, ] * p_HIV_HCV_ODN_INJ # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & INJ, ODN & INJ & HCV, ]  <- a_TDP[HCV & INJ, ODN & INJ & HCV, ] * (1 - p_HCV_HIV_ODN_INJ)
  a_TDP[HCV & INJ, ODN & INJ & COI, ]  <- a_TDP[HCV & INJ, ODN & INJ & COI, ] * p_HCV_HIV_ODN_INJ # Probability of HIV conditional on HCV
  
  # ODF
  # From NEG
  a_TDP[NEG & INJ, ODF & INJ & NEG, ]  <- a_TDP[NEG & INJ, ODF & INJ & NEG, ] * 1
  a_TDP[NEG & INJ, ODF & INJ & HIV, ]  <- a_TDP[NEG & INJ, ODF & INJ & HIV, ] * 0
  a_TDP[NEG & INJ, ODF & INJ & HCV, ]  <- a_TDP[NEG & INJ, ODF & INJ & HCV, ] * 0
  # From HIV
  a_TDP[HIV & INJ, ODF & INJ & HIV, ]  <- a_TDP[HIV & INJ, ODF & INJ & HIV, ] * 1
  a_TDP[HIV & INJ, ODF & INJ & COI, ]  <- a_TDP[HIV & INJ, ODF & INJ & COI, ] * 0 # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & INJ, ODF & INJ & HCV, ]  <- a_TDP[HCV & INJ, ODF & INJ & HCV, ] * 1
  a_TDP[HCV & INJ, ODF & INJ & COI, ]  <- a_TDP[HCV & INJ, ODF & INJ & COI, ] * 0 # Probability of HIV conditional on HCV
  
  # ABS
  # From NEG
  a_TDP[NEG & INJ, ABS & INJ & NEG, ]  <- a_TDP[NEG & INJ, ABS & INJ & NEG, ] * (1 - p_HIV_ABS_INJ - p_HCV_ABS_INJ)
  a_TDP[NEG & INJ, ABS & INJ & HIV, ]  <- a_TDP[NEG & INJ, ABS & INJ & HIV, ] * p_HIV_ABS_INJ
  a_TDP[NEG & INJ, ABS & INJ & HCV, ]  <- a_TDP[NEG & INJ, ABS & INJ & HCV, ] * p_HCV_ABS_INJ
  # From HIV
  a_TDP[HIV & INJ, ABS & INJ & HIV, ]  <- a_TDP[HIV & INJ, ABS & INJ & HIV, ] * (1 - p_HIV_HCV_ABS_INJ)
  a_TDP[HIV & INJ, ABS & INJ & COI, ]  <- a_TDP[HIV & INJ, ABS & INJ & COI, ] * p_HIV_HCV_ABS_INJ # Probability of HCV conditional on HIV
  # From HCV
  a_TDP[HCV & INJ, ABS & INJ & HCV, ]  <- a_TDP[HCV & INJ, ABS & INJ & HCV, ] * (1 - p_HCV_HIV_ABS_INJ)
  a_TDP[HCV & INJ, ABS & INJ & COI, ]  <- a_TDP[HCV & INJ, ABS & INJ & COI, ] * p_HCV_HIV_ABS_INJ # Probability of HIV conditional on HCV

  # Disallowed transitions (ensure that impossible transitions are set to zero)
  # Episode rules
  a_TDP[EP1, EP3, ] = 0
  a_TDP[EP2, EP1, ] = 0
  a_TDP[EP3, EP1, ] = 0
  a_TDP[EP3, EP2, ] = 0
  # Seroconversions
  a_TDP[HIV, NEG, ] = 0
  a_TDP[HCV, NEG, ] = 0 # disallowing potential transitions from COI to HIV-only (i.e. HCV cure), calculated within overall HCV infection rate
  a_TDP[COI, NEG, ] = 0
  a_TDP[HIV, HCV, ] = 0
  a_TDP[HCV, HIV, ] = 0
  a_TDP[COI, HIV, ] = 0 # disallowing potential transitions from COI to HIV-only (i.e. HCV cure), calculated within overall HCV infection rate
  a_TDP[COI, HCV, ] = 0
  a_TDP[COI, NEG, ] = 0
  a_TDP[NEG, COI, ] = 0
  # Abstinence directly to treatment
  #a_TDP[ABS, TX, ]  = 0
  # Conditional transitions
  # Next episode with out-of-treatment(OOT) EPi -> treatment(TX) EP(i+1)
  a_TDP[TX & EP1, TX & EP2, ] = 0
  a_TDP[TX & EP1, TX & EP3, ] = 0
  a_TDP[TX & EP2, TX & EP3, ] = 0
  a_TDP[OOT & EP1, OOT & EP2, ] = 0
  a_TDP[OOT & EP1, OOT & EP3, ] = 0
  a_TDP[OOT & EP2, OOT & EP3, ] = 0
  a_TDP[TX & EP1, OOT & EP2, ] = 0
  a_TDP[TX & EP2, OOT & EP3, ] = 0
  a_TDP[OOT & EP1, TX & EP1, ] = 0
  a_TDP[OOT & EP2, TX & EP2, ] = 0

  k_hiv <- a_TDP[, , 50]
  k_hiv2<- a_TDP[, , 710]

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
  
  # Episode 2
  #v_s_init[BUP & EP2] <- v_init_dist["pe", "BUP"] # Empirically observed proportions from base states
  #v_s_init[MET & EP2] <- v_init_dist["pe", "MET"]
  #v_s_init[REL & EP2] <- v_init_dist["pe", "REL"]
  #v_s_init[ODN & EP2] <- v_init_dist["pe", "ODN"]
  #v_s_init[ODF & EP2] <- v_init_dist["pe", "ODF"]
  #v_s_init[ABS & EP2] <- v_init_dist["pe", "ABS"]
  
  # Episode 3+
  #v_s_init[BUP & EP3] <- v_init_dist["pe", "BUP"] # Empirically observed proportions from base states
  #v_s_init[MET & EP3] <- v_init_dist["pe", "MET"]
  #v_s_init[REL & EP3] <- v_init_dist["pe", "REL"]
  #v_s_init[ODN & EP3] <- v_init_dist["pe", "ODN"]
  #v_s_init[ODF & EP3] <- v_init_dist["pe", "ODF"]
  #v_s_init[ABS & EP3] <- v_init_dist["pe", "ABS"]
  
  # Distribute by injection/non-injection
  v_s_init[NI]  <- v_s_init[NI] * (1 - n_INJ)
  v_s_init[INJ] <- v_s_init[INJ] * n_INJ
  
  # Distribute HIV/HCV/COI
  v_s_init[NEG] <- v_s_init[NEG] * (1 - n_HIV - n_HCV - n_COI)
  v_s_init[HIV] <- v_s_init[HIV] * n_HIV
  v_s_init[HCV] <- v_s_init[HCV] * n_HCV
  v_s_init[COI] <- v_s_init[COI] * n_COI

  # Create Markov Trace
    # Initialize population
    a_M_trace <- array(0, dim = c((n_t + 1), n_states, (n_t + 1)),
                       dimnames = list(0:n_t, v_n_states, 0:n_t))
    a_M_trace[1, , 1] <- v_s_init

    # All model time periods
      for(i in 2:(n_t)){
        # Time spent in given health state
        for(j in 1:(i - 1)){
          #state-time-dependent transition probability (j) * age (model-time)-specific mortality (i)
          m_sojourn <- a_TDP[, , j] * m_alive[, i - 1]
          
          v_current_state <- as.vector(a_M_trace[i - 1, , j])
          
          v_same_state <- as.vector(v_current_state * diag(m_sojourn))
          
          a_M_trace[i, ,j + 1] <- v_same_state
          
          diag(m_sojourn) <- 0
          
          v_new_state <- as.vector(v_current_state %*% m_sojourn)
          
          a_M_trace[i, ,1] <- v_new_state + a_M_trace[i, ,1]
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
  write.csv(m_M_agg_trace,"checks/m_M_agg_trace.csv", row.names = TRUE)
  
  # Model diagnostic outputs
  # Output csv checks if true
  if(checks){
    # Time dependent state-exit probabilities (from weibull estimates)
    write.csv(m_TDP,"C:/Users/Benjamin/Desktop/m_TDP.csv", row.names = TRUE)
    # Mortality matrix
    write.csv(m_mort,"checks/m_mort.csv", row.names = TRUE)
    # Unconditional state transitions
    write.csv(m_UP,"checks/m_UP.csv", row.names = TRUE)
    # First month
    write.csv(m_UP_4wk,"checks/m_UP_4wk.csv", row.names = TRUE)
    # First period full array
    write.csv(initial,"checks/initial_full.csv", row.names = TRUE)
    # Array at time = 50 months
    write.csv(k_hiv,"checks/K2.csv", row.names = TRUE)
    # Array at time = 710 months
    write.csv(k_hiv2,"checks/K3.csv", row.names = TRUE)
    # Initial state occupancy at model initiation
    write.csv(v_s_init,"checks/v_s_init.csv", row.names = TRUE)
  } else{}
  
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
              a_TDP = a_TDP,
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
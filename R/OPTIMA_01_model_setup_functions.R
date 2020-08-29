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
markov_model <- function(l_params_all, err_stop = FALSE, verbose = FALSE){
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
  BASE <- l_dim_s[[1]] <- c("MET", "BUP", "ABS", "REL", "ODN", "ODF")
  # Injection/non-injection stratification
  INJECT <- l_dim_s[[2]] <- c("NI", "INJ")
  # Episodes (1-3)
  EP <-  l_dim_s[[3]] <- c("1", "2", "3")
  # HIV status
  HIV <- l_dim_s[[4]] <- c("HIV", "HCV", "COI", "NEG")
  n_t <- (n_age_max - n_age_init) * 52 # convert years into weeks
  
  df_flat <- expand.grid(l_dim_s) #combine all elements together into vector of health states
  df_flat <- rename(df_flat, BASE    = Var1, 
                             INJECT  = Var2, 
                             EP      = Var3, 
                             HIV     = Var4)

  # Create index of states to populate transition matrices
  # All treatment
  TX <- df_flat$BASE == "BUP" | df_flat$BASE == "MET"
  
  # All out-of-treatment (incl ABS)
  OOT <- df_flat$BASE == "REL" | df_flat$BASE == "ODN" | df_flat$BASE == "ODF" | df_flat$BASE == "ABS"
  
  # Buprenorphine
  #all_BUP <- df_flat$BASE == "BUP" | df_flat$BASE == "BUP1"
  BUP     <- df_flat$BASE == "BUP"
  #BUP1    <- df_flat$BASE == "BUP1"
  
  # Methadone
  #all_MET <- df_flat$BASE == "MET" | df_flat$BASE == "MET1"
  MET     <- df_flat$BASE == "MET"
  #MET1    <- df_flat$BASE == "MET1"
  
  # Relapse
  #all_REL <- df_flat$BASE == "REL" | df_flat$BASE == "REL1"
  REL <- df_flat$BASE == "REL"
  #REL1 <- df_flat$BASE == "REL1"

  # Overdose
  all_OD <- df_flat$BASE == "ODN" | df_flat$BASE == "ODF"
  ODN <- df_flat$BASE == "ODN"
  ODF <- df_flat$BASE == "ODF"
  
  # Abstinence
  ABS <- df_flat$BASE == "ABS"
  
  # Serostatus
  NEG <- df_flat$HIV == "NEG"
  HIV <- df_flat$HIV == "HIV"
  HCV <- df_flat$HIV == "HCV"
  COI <- df_flat$HIV == "COI"
  
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
                     BUP = BUP, 
                     MET = MET, 
                     REL = REL, 
                     all_OD = all_OD, ODN = ODN, ODF = ODF, 
                     ABS = ABS, 
                     NEG = NEG, HIV = HIV, HCV = HCV, COI = COI, 
                     INJ = INJ, NI = NI, 
                     EP1 = EP1, EP2 = EP2, EP3 = EP3)
  
  #### Overdose probability ####
  #' Probability of non-fatal and fatal overdose
  #'
  #' \code{p_OD} is used to calculate weekly overdose probabilities from health states. This function also requires additional overdose/fentanyl/naloxone parameters included in `l_params_all`
  #'
  #' @param p_from_state_OD Baseline overdose probability from each health state
  #' @param fatal Logical parameter to switch between fatal/non-fatal overdose
  #' 
  #' @return 
  #' `p_OD` weekly probability of fatal or non-fatal overdose from a given health state
  #' @export
  p_OD <- function(p_from_state_OD = p_from_state_OD,
                   fatal = FALSE){
    
    # Probability of successful naloxone use
    p_NX_rev <- (p_witness * p_NX_used * p_NX_success)
    
    # Probability of mortality from overdose accounting for baseline overdose fatality and effectiveness of naloxone
    # Subsets overdose into fatal and non-fatal, conditional on different parameters
    p_fatal_OD_NX <- p_fatal_OD * (1 - p_NX_rev)
    
    if (fatal = FALSE){
      p_OD <- ((from_state * (1 - p_fent_exp)) + (p_fent_OD * (p_fent_exp))) * (1 - p_fatal_OD_NX)
    } else if (fatal = TRUE){
        p_OD <- ((from_state * (1 - p_fent_exp)) + (p_fent_OD * (p_fent_exp))) * p_fatal_OD_NX
    }
    return(p_OD)
  }

  # Module to calculate probability of overdose from states
  # Probability of overdose
  # Non-injection
  p_BUP_ODN_NI  <- p_OD(p_from_state_OD = p_BUP_OD_NI, fatal = FALSE)
  p_MET_ODN_NI  <- p_OD(p_from_state_OD = p_MET_OD_NI, fatal = FALSE)
  p_REL_ODN_NI  <- p_OD(p_from_state_OD = p_REL_OD_NI, fatal = FALSE)
  p_ABS_ODN_NI  <- p_OD(p_from_state_OD = p_ABS_OD_NI, fatal = FALSE)
  p_BUP_ODF_NI  <- p_OD(p_from_state_OD = p_BUP_OD_NI, fatal = TRUE)
  p_MET_ODF_NI  <- p_OD(p_from_state_OD = p_MET_OD_NI, fatal = TRUE)
  p_REL_ODF_NI  <- p_OD(p_from_state_OD = p_REL_OD_NI, fatal = TRUE)
  p_ABS_ODF_NI  <- p_OD(p_from_state_OD = p_ABS_OD_NI, fatal = TRUE)
  
  # Injection
  p_BUP_ODN_INJ <- p_OD(p_from_state_OD = p_BUP_OD_INJ, fatal = FALSE)
  p_MET_ODN_INJ <- p_OD(p_from_state_OD = p_MET_OD_INJ, fatal = FALSE)
  p_REL_ODN_INJ <- p_OD(p_from_state_OD = p_REL_OD_INJ, fatal = FALSE)
  p_ABS_ODN_INJ <- p_OD(p_from_state_OD = p_ABS_OD_INJ, fatal = FALSE)
  p_BUP_ODF_INJ <- p_OD(p_from_state_OD = p_BUP_OD_INJ, fatal = TRUE)
  p_MET_ODF_INJ <- p_OD(p_from_state_OD = p_MET_OD_INJ, fatal = TRUE)
  p_REL_ODF_INJ <- p_OD(p_from_state_OD = p_REL_OD_INJ, fatal = TRUE)
  p_ABS_ODF_INJ <- p_OD(p_from_state_OD = p_ABS_OD_INJ, fatal = TRUE)
  
  #### Time-dependent survival probabilities ####
    # Empty 2-D matrix
    m_TDP <- array(0, dim = c(n_states, n_t),
                      dimnames = list(v_n_states, 1:n_t))

    # Probability of remaining in health state
    # All transition out of non-fatal overdose after 1 week, probs already set to zero
    # All remain in fatal overdose, absorbing state
    for(i in 1:n_t){
      # Non-injection
        # Episode 1
        m_TDP[EP1 & BUP & NI, i] <- as.vector(exp(p_frailty_BUP_NI_1 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI)))) # vector of remain probabilities
        m_TDP[EP1 & MET & NI, i] <- as.vector(exp(p_frailty_MET_NI_1 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        m_TDP[EP1 & ABS & NI, i] <- as.vector(exp(p_frailty_ABS_NI_1 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP1 & REL & NI, i] <- as.vector(exp(p_frailty_REL_NI_1 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
        m_TDP[EP1 & ODF & NI, i] <- rep(1, n_t)  
        # Episode 2
        m_TDP[EP2 & BUP & NI, i] <- as.vector(exp(p_frailty_BUP_NI_2 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI))))
        m_TDP[EP2 & MET & NI, i] <- as.vector(exp(p_frailty_MET_NI_2 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        m_TDP[EP2 & ABS & NI, i] <- as.vector(exp(p_frailty_ABS_NI_2 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP2 & REL & NI, i] <- as.vector(exp(p_frailty_REL_NI_2 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
        m_TDP[EP2 & ODF & NI, i] <- rep(1, n_t)
        # Episode 3
        m_TDP[EP3 & BUP & NI, i] <- as.vector(exp(p_frailty_BUP_NI_3 * p_weibull_scale_BUP_NI * (((i - 1)^p_weibull_shape_BUP_NI) - (i^p_weibull_shape_BUP_NI))))
        m_TDP[EP3 & MET & NI, i] <- as.vector(exp(p_frailty_MET_NI_3 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        m_TDP[EP3 & ABS & NI, i] <- as.vector(exp(p_frailty_ABS_NI_3 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP3 & REL & NI, i] <- as.vector(exp(p_frailty_REL_NI_3 * p_weibull_scale_REL_NI * (((i - 1)^p_weibull_shape_REL_NI) - (i^p_weibull_shape_REL_NI))))
        m_TDP[EP3 & ODF & NI, i] <- rep(1, n_t)
        
    # Injection
        # Episode 1
        m_TDP[EP1 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_1 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ)))) # vector of remain probabilities 
        m_TDP[EP1 & MET & INJ, i] <- as.vector(exp(p_frailty_MET_INJ_1 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        m_TDP[EP1 & ABS & INJ, i] <- as.vector(exp(p_frailty_ABS_INJ_1 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP1 & REL & INJ, i] <- as.vector(exp(p_frailty_REL_INJ_1 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
        m_TDP[EP1 & ODF & INJ, i] <- rep(1, n_t)
        # Episode 2
        m_TDP[EP2 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_2 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ))))
        m_TDP[EP2 & MET & INJ, i] <- as.vector(exp(p_frailty_MET_INJ_2 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        m_TDP[EP2 & ABS & INJ, i] <- as.vector(exp(p_frailty_ABS_INJ_2 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP2 & REL & INJ, i] <- as.vector(exp(p_frailty_REL_INJ_2 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
        m_TDP[EP2 & ODF & INJ, i] <- rep(1, n_t)
        # Episode 3
        m_TDP[EP3 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_3 * p_weibull_scale_BUP_INJ * (((i - 1)^p_weibull_shape_BUP_INJ) - (i^p_weibull_shape_BUP_INJ))))
        m_TDP[EP3 & MET & INJ, i] <- as.vector(exp(p_frailty_MET_INJ_3 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        m_TDP[EP3 & ABS & INJ, i] <- as.vector(exp(p_frailty_ABS_INJ_3 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP3 & REL & INJ, i] <- as.vector(exp(p_frailty_REL_INJ_3 * p_weibull_scale_REL_INJ * (((i - 1)^p_weibull_shape_REL_INJ) - (i^p_weibull_shape_REL_INJ))))
        m_TDP[EP3 & ODF & INJ, i] <- rep(1, n_t)
  }

  # Probability of state-exit
  m_leave <- 1 - m_TDP
  
  #### Mortality ####
  #' Weekly mortality estimates
  #'
  #' \code{v_mort} is used to populate mortality probability vectors.
  #'
  #' @param hr
  #' @return 
  #' Mortality vectors for each age applied to 52 weeks, includes state-specific hr.
  #' Overdose deaths tracked separately as "ODF"
  #' @export
  v_mort <- function(hr = hr){
    v_mort <- rep((1 - exp(-v_r_mort_by_age[n_age_init:(n_age_max - 1), ] * (1/52) * hr)), each = 52)
    return(v_mort)
  }
  # Non-injection
  v_mort_BUP_NI     <- v_mort(hr = hr_BUP_NI)
  v_mort_MET_NI     <- v_mort(hr = hr_MET_NI)
  v_mort_REL_NI     <- v_mort(hr = hr_REL_NI)
  v_mort_ODN_NI     <- v_mort(hr = hr_ODN_NI) # Mortality equal for REL/ODN (fatal overdoses counted separately)
  v_mort_ODF_NI     <- rep(0, n_t) # mortality transition = 0 as death already tracked in ODF
  v_mort_ABS_NEG_NI <- v_mort(hr = hr_ABS_NI)
  v_mort_ABS_POS_NI <- v_mort(hr = hr_HIV_NI)
  
  # Injection
  v_mort_BUP_INJ     <- v_mort(hr = hr_BUP_INJ)
  v_mort_MET_INJ     <- v_mort(hr = hr_MET_INJ)
  v_mort_REL_INJ     <- v_mort(hr = hr_REL_INJ)
  v_mort_ODN_INJ     <- v_mort(hr = hr_ODN_INJ) # Mortality equal for REL/ODN (fatal overdoses counted separately)
  v_mort_ODF_INJ     <- rep(0, n_t) # mortality transition = 0 as death already tracked in ODF
  v_mort_ABS_NEG_INJ <- v_mort(hr = hr_ABS_INJ)
  v_mort_ABS_POS_INJ <- v_mort(hr = hr_HIV_INJ)

  # Create empty mortality matrix
  m_mort <- array(0, dim = c(n_states, n_t),
                  dimnames = list(v_n_states, 1:n_t))
  # Populate mortality matrix (weekly death probability from each state)
  for (i in 1:n_t){
    # Non-injection
    m_mort[BUP & NI, i]       <- v_mort_BUP_NI[i]
    m_mort[MET & NI, i]       <- v_mort_MET_NI[i]
    m_mort[REL & NI, i]       <- v_mort_REL_NI[i]
    m_mort[ODN & NI, i]       <- v_mort_ODN_NI[i] # using background excess mortality for relapse in non-fatal overdose
    m_mort[ODF & NI, i]       <- v_mort_ODF_NI[i] # already counted as fatal overdoses, e.g. transition to death = 1
    m_mort[ABS & NI & NEG, i] <- v_mort_ABS_NEG_NI[i]
    m_mort[ABS & NI & POS, i] <- v_mort_ABS_POS_NI[i]
    # Injection
    m_mort[BUP & INJ, i]       <- v_mort_BUP_INJ[i]
    m_mort[MET & INJ, i]       <- v_mort_MET_INJ[i]
    m_mort[REL & INJ, i]       <- v_mort_REL_INJ[i]
    m_mort[ODN & INJ, i]       <- v_mort_ODN_INJ[i] # using background excess mortality for relapse in non-fatal overdose
    m_mort[ODF & INJ, i]       <- v_mort_ODF_INJ[i] # already counted as fatal overdoses, e.g. transition to death = 1
    m_mort[ABS & INJ & NEG, i] <- v_mort_ABS_NEG_INJ[i]
    m_mort[ABS & INJ & POS, i] <- v_mort_ABS_POS_INJ[i]
  }

  # Alive probability in each period
  m_alive <- 1 - m_mort

  #### Unconditional transition probabilities ####
  # Empty 2-D unconditional transition matrix (from states, to states)
  m_UP <- m_UP_4wk <- array(0, dim = c(n_states, n_states),
                            dimnames = list(v_n_states, v_n_states))
  # Populate unconditional transition matrix
  # Overdose probability populated first, accounting for higher probability of overdose transition in first 4 weeks of BUP/MET/REL
  # Non-Injection
  # From BUP
  # Overall
  m_UP[BUP & NI, MET & NI] <- p_BUP_MET_NI * (1 - p_BUP_ODN_NI - p_BUP_ODF_NI)
  m_UP[BUP & NI, ABS & NI] <- p_BUP_ABS_NI * (1 - p_BUP_ODN_NI - p_BUP_ODF_NI)
  m_UP[BUP & NI, REL & NI] <- p_BUP_REL_NI * (1 - p_BUP_ODN_NI - p_BUP_ODF_NI)
  m_UP[BUP & NI, ODN & NI] <- p_BUP_ODN_NI
  m_UP[BUP & NI, ODF & NI] <- p_BUP_ODF_NI
  # First 4 weeks
  m_UP_4wk[BUP & NI, MET & NI] <- p_BUP_MET_NI * (1 - (p_BUP_ODN_NI * p_BUP_OD_mult) - (p_BUP_ODF_NI * p_BUP_OD_mult))
  m_UP_4wk[BUP & NI, ABS & NI] <- p_BUP_ABS_NI * (1 - (p_BUP_ODN_NI * p_BUP_OD_mult) - (p_BUP_ODF_NI * p_BUP_OD_mult))
  m_UP_4wk[BUP & NI, REL & NI] <- p_BUP_REL_NI * (1 - (p_BUP_ODN_NI * p_BUP_OD_mult) - (p_BUP_ODF_NI * p_BUP_OD_mult))
  m_UP_4wk[BUP & NI, ODN & NI] <- p_BUP_ODN_NI * p_BUP_OD_mult # Consider estimating probability directly to avoid multiplier issue
  m_UP_4wk[BUP & NI, ODF & NI] <- p_BUP_ODF_NI * p_BUP_OD_mult # Consider estimating probability directly to avoid multiplier issue
  
  # From MET
  # Overall
  m_UP[MET & NI, BUP & NI] <- p_MET_BUP_NI * (1 - p_MET_ODN_NI - p_MET_ODF_NI)
  m_UP[MET & NI, ABS & NI] <- p_MET_ABS_NI * (1 - p_MET_ODN_NI - p_MET_ODF_NI)
  m_UP[MET & NI, REL & NI] <- p_MET_REL_NI * (1 - p_MET_ODN_NI - p_MET_ODF_NI)
  m_UP[MET & NI, ODN & NI] <- p_MET_ODN_NI
  m_UP[MET & NI, ODF & NI] <- p_MET_ODF_NI
  # First 4 weeks
  m_UP_4wk[MET & NI, BUP & NI] <- p_MET_BUP_NI * (1 - (p_MET_ODN_NI * p_MET_OD_mult) - (p_MET_ODF_NI * p_MET_OD_mult))
  m_UP_4wk[MET & NI, ABS & NI] <- p_MET_ABS_NI * (1 - (p_MET_ODN_NI * p_MET_OD_mult) - (p_MET_ODF_NI * p_MET_OD_mult))
  m_UP_4wk[MET & NI, REL & NI] <- p_MET_REL_NI * (1 - (p_MET_ODN_NI * p_MET_OD_mult) - (p_MET_ODF_NI * p_MET_OD_mult))
  m_UP_4wk[MET & NI, ODN & NI] <- p_MET_ODN_NI * p_MET_OD_mult # Consider estimating probability directly to avoid multiplier issue
  m_UP_4wk[MET & NI, ODF & NI] <- p_MET_ODF_NI * p_MET_OD_mult # Consider estimating probability directly to avoid multiplier issue
  
  # From ABS
  m_UP[ABS & NI, REL & NI] <- p_ABS_REL_NI * (1 - p_ABS_OD_NI)
  m_UP[ABS & NI, ODN & NI] <- p_ABS_OD_NI * (1 - p_fatal_OD_NI)
  m_UP[ABS & NI, ODF & NI] <- p_ABS_OD_NI * p_fatal_OD_NI
  
  # From REL
  # Overall
  m_UP[REL & NI, MET & NI] <- p_REL_MET_NI * (1 - p_REL_ODN_NI - p_REL_ODF_NI)
  m_UP[REL & NI, BUP & NI] <- p_REL_BUP_NI * (1 - p_REL_ODN_NI - p_REL_ODF_NI)
  m_UP[REL & NI, ABS & NI] <- p_REL_ABS_NI * (1 - p_REL_ODN_NI - p_REL_ODF_NI)
  m_UP[REL & NI, ODN & NI] <- p_REL_ODN_NI
  m_UP[REL & NI, ODF & NI] <- p_REL_ODF_NI
  # First 4 weeks
  m_UP_4wk[REL & NI, MET & NI] <- p_REL_MET_NI * (1 - (p_REL_ODN_NI * p_REL_OD_mult) - (p_REL_ODF_NI * p_REL_OD_mult))
  m_UP_4wk[REL & NI, BUP & NI] <- p_REL_BUP_NI * (1 - (p_REL_ODN_NI * p_REL_OD_mult) - (p_REL_ODF_NI * p_REL_OD_mult))
  m_UP_4wk[REL & NI, ABS & NI] <- p_REL_ABS_NI * (1 - (p_REL_ODN_NI * p_REL_OD_mult) - (p_REL_ODF_NI * p_REL_OD_mult))
  m_UP_4wk[REL & NI, ODN & NI] <- p_REL_ODN_NI * p_REL_OD_mult # Consider estimating probability directly to avoid multiplier issue
  m_UP_4wk[REL & NI, ODF & NI] <- p_REL_ODF_NI * p_REL_OD_mult # Consider estimating probability directly to avoid multiplier issue
  
  # From OD
  m_UP[ODN & NI, MET & NI] <- p_ODN_MET_NI
  m_UP[ODN & NI, BUP & NI] <- p_ODN_BUP_NI
  m_UP[ODN & NI, ABS & NI] <- p_ODN_ABS_NI
  m_UP[ODN & NI, REL & NI] <- p_ODN_REL_NI

  # Injection
  # From BUP
  # Overall
  m_UP[BUP & INJ, MET & INJ] <- p_BUP_MET_INJ * (1 - p_BUP_ODN_INJ - p_BUP_ODF_INJ)
  m_UP[BUP & INJ, ABS & INJ] <- p_BUP_ABS_INJ * (1 - p_BUP_ODN_INJ - p_BUP_ODF_INJ)
  m_UP[BUP & INJ, REL & INJ] <- p_BUP_REL_INJ * (1 - p_BUP_ODN_INJ - p_BUP_ODF_INJ)
  m_UP[BUP & INJ, ODN & INJ] <- p_BUP_ODN_INJ
  m_UP[BUP & INJ, ODF & INJ] <- p_BUP_ODF_INJ
  # First 4 weeks
  m_UP_4wk[BUP & INJ, MET & INJ] <- p_BUP_MET_INJ * (1 - (p_BUP_ODN_INJ * p_BUP_OD_mult) - (p_BUP_ODF_INJ * p_BUP_OD_mult))
  m_UP_4wk[BUP & INJ, ABS & INJ] <- p_BUP_ABS_INJ * (1 - (p_BUP_ODN_INJ * p_BUP_OD_mult) - (p_BUP_ODF_INJ * p_BUP_OD_mult))
  m_UP_4wk[BUP & INJ, REL & INJ] <- p_BUP_REL_INJ * (1 - (p_BUP_ODN_INJ * p_BUP_OD_mult) - (p_BUP_ODF_INJ * p_BUP_OD_mult))
  m_UP_4wk[BUP & INJ, ODN & INJ] <- p_BUP_ODN_INJ * p_BUP_OD_mult
  m_UP_4wk[BUP & INJ, ODF & INJ] <- p_BUP_ODF_INJ * p_BUP_OD_mult
  
  # From MET
  # Overall
  m_UP[MET & INJ, BUP & INJ] <- p_MET_BUP_INJ * (1 - p_MET_ODN_INJ - p_MET_ODF_INJ)
  m_UP[MET & INJ, ABS & INJ] <- p_MET_ABS_INJ * (1 - p_MET_ODN_INJ - p_MET_ODF_INJ)
  m_UP[MET & INJ, REL & INJ] <- p_MET_REL_INJ * (1 - p_MET_ODN_INJ - p_MET_ODF_INJ)
  m_UP[MET & INJ, ODN & INJ] <- p_MET_ODN_INJ
  m_UP[MET & INJ, ODF & INJ] <- p_MET_ODF_INJ
  # First 4 weeks
  m_UP_4wk[MET & INJ, BUP & INJ] <- p_MET_BUP_INJ * (1 - (p_MET_ODN_INJ * p_MET_OD_mult) - (p_MET_ODF_INJ * p_MET_OD_mult))
  m_UP_4wk[MET & INJ, ABS & INJ] <- p_MET_ABS_INJ * (1 - (p_MET_ODN_INJ * p_MET_OD_mult) - (p_MET_ODF_INJ * p_MET_OD_mult))
  m_UP_4wk[MET & INJ, REL & INJ] <- p_MET_REL_INJ * (1 - (p_MET_ODN_INJ * p_MET_OD_mult) - (p_MET_ODF_INJ * p_MET_OD_mult))
  m_UP_4wk[MET & INJ, ODN & INJ] <- p_MET_ODN_INJ * p_MET_OD_mult # Consider estimating probability directly to avoid multiplier issue
  m_UP_4wk[MET & INJ, ODF & INJ] <- p_MET_ODF_INJ * p_MET_OD_mult # Consider estimating probability directly to avoid multiplier issue
  
  # From ABS
  m_UP[ABS & INJ, REL & INJ] <- p_ABS_REL_INJ * (1 - p_ABS_OD_INJ)
  m_UP[ABS & INJ, ODN & INJ] <- p_ABS_OD_INJ * (1 - p_fatal_OD_INJ)
  m_UP[ABS & INJ, ODF & INJ] <- p_ABS_OD_INJ * p_fatal_OD_INJ
  
  # From REL
  # Overall
  m_UP[REL & INJ, MET & INJ] <- p_REL_MET_INJ * (1 - p_REL_ODN_INJ - p_REL_ODF_INJ)
  m_UP[REL & INJ, BUP & INJ] <- p_REL_BUP_INJ * (1 - p_REL_ODN_INJ - p_REL_ODF_INJ)
  m_UP[REL & INJ, ABS & INJ] <- p_REL_ABS_INJ * (1 - p_REL_ODN_INJ - p_REL_ODF_INJ)
  m_UP[REL & INJ, ODN & INJ] <- p_REL_ODN_INJ
  m_UP[REL & INJ, ODF & INJ] <- p_REL_ODF_INJ
  # First 4 weeks
  m_UP_4wk[REL & INJ, MET & INJ] <- p_REL_MET_INJ * (1 - (p_REL_ODN_INJ * p_REL_OD_mult) - (p_REL_ODF_INJ * p_REL_OD_mult))
  m_UP_4wk[REL & INJ, BUP & INJ] <- p_REL_BUP_INJ * (1 - (p_REL_ODN_INJ * p_REL_OD_mult) - (p_REL_ODF_INJ * p_REL_OD_mult))
  m_UP_4wk[REL & INJ, ABS & INJ] <- p_REL_ABS_INJ * (1 - (p_REL_ODN_INJ * p_REL_OD_mult) - (p_REL_ODF_INJ * p_REL_OD_mult))
  m_UP_4wk[REL & INJ, ODN & INJ] <- p_REL_ODN_INJ * p_REL_OD_mult # Consider estimating probability directly to avoid multiplier issue
  m_UP_4wk[REL & INJ, ODF & INJ] <- p_REL_ODF_INJ * p_REL_OD_mult # Consider estimating probability directly to avoid multiplier issue
  
  # From OD
  m_UP[ODN & INJ, MET & INJ] <- p_ODN_MET_INJ
  m_UP[ODN & INJ, BUP & INJ] <- p_ODN_BUP_INJ
  m_UP[ODN & INJ, ABS & INJ] <- p_ODN_ABS_INJ
  m_UP[ODN & INJ, REL & INJ] <- p_ODN_REL_INJ

  # Checks
  #write.csv(m_UP,"C:/Users/Benjamin/Desktop/m_UP.csv", row.names = TRUE)
  
  #### Create full time-dependent transition array ####
  # Empty 3-D array
  a_TDP <- array(0, dim = c(n_states, n_states, n_t),
                   dimnames = list(v_n_states, v_n_states, 1:n_t))
  
  # Add transitions conditional on state-exit (m_leave = 1 - remain)
  # Modified transitions for first 4-weeks
  for (i in 1:4){
    a_TDP[, , i] <- m_UP_4wk * m_leave[, i]
  }
  # All transitions 5+ weeks
  for (i in 5:n_t){
    a_TDP[, , i] <- m_UP * m_leave[, i]
  }

  # Add time-dependent remain probabilities
  for (i in 1:n_t){
    for (j in 1:n_states){
   a_TDP[j, j, i] <- m_TDP[j, i]
    } 
  }

  # Add NEG -> POS remain probabilities
  # To-do: See if there is a better way to do this
  # AUGUST 28, 2020 **UPDATE SEROCONVERSION SECTION TO INCLUDE HCV AND COI**
  for (i in 1:n_t){
    # Non-injection
    # BUP
    a_TDP[BUP & NI & EP1 & NEG, BUP & NI & EP1 & HIV, i] <- m_TDP[BUP & NI & EP1 & NEG, i]
    a_TDP[BUP & NI & EP2 & NEG, BUP & NI & EP2 & HIV, i] <- m_TDP[BUP & NI & EP2 & NEG, i]
    a_TDP[BUP & NI & EP3 & NEG, BUP & NI & EP3 & HIV, i] <- m_TDP[BUP & NI & EP3 & NEG, i]
    
    a_TDP[BUP & NI & EP1 & NEG, BUP & NI & EP1 & HCV, i] <- m_TDP[BUP & NI & EP1 & NEG, i]
    a_TDP[BUP & NI & EP2 & NEG, BUP & NI & EP2 & HCV, i] <- m_TDP[BUP & NI & EP2 & NEG, i]
    a_TDP[BUP & NI & EP3 & NEG, BUP & NI & EP3 & HCV, i] <- m_TDP[BUP & NI & EP3 & NEG, i]
    
    a_TDP[BUP & NI & EP1 & NEG, BUP & NI & EP1 & COI, i] <- m_TDP[BUP & NI & EP1 & NEG, i]
    a_TDP[BUP & NI & EP2 & NEG, BUP & NI & EP2 & COI, i] <- m_TDP[BUP & NI & EP2 & NEG, i]
    a_TDP[BUP & NI & EP3 & NEG, BUP & NI & EP3 & COI, i] <- m_TDP[BUP & NI & EP3 & NEG, i]
    
    a_TDP[BUP & NI & EP1 & HIV, BUP & NI & EP1 & COI, i] <- m_TDP[BUP & NI & EP1 & NEG, i]
    a_TDP[BUP & NI & EP2 & HIV, BUP & NI & EP2 & COI, i] <- m_TDP[BUP & NI & EP2 & NEG, i]
    a_TDP[BUP & NI & EP3 & HIV, BUP & NI & EP3 & COI, i] <- m_TDP[BUP & NI & EP3 & NEG, i]
    
    a_TDP[BUP & NI & EP1 & HCV, BUP & NI & EP1 & COI, i] <- m_TDP[BUP & NI & EP1 & NEG, i]
    a_TDP[BUP & NI & EP2 & HCV, BUP & NI & EP2 & COI, i] <- m_TDP[BUP & NI & EP2 & NEG, i]
    a_TDP[BUP & NI & EP3 & HCV, BUP & NI & EP3 & COI, i] <- m_TDP[BUP & NI & EP3 & NEG, i]
    
    # MET
    a_TDP[MET & NI & EP1 & NEG, MET & NI & EP1 & HIV, i] <- m_TDP[MET & NI & EP1 & NEG, i]
    a_TDP[MET & NI & EP2 & NEG, MET & NI & EP2 & HIV, i] <- m_TDP[MET & NI & EP2 & NEG, i]
    a_TDP[MET & NI & EP3 & NEG, MET & NI & EP3 & HIV, i] <- m_TDP[MET & NI & EP3 & NEG, i]
    
    a_TDP[MET & NI & EP1 & NEG, MET & NI & EP1 & HCV, i] <- m_TDP[MET & NI & EP1 & NEG, i]
    a_TDP[MET & NI & EP2 & NEG, MET & NI & EP2 & HCV, i] <- m_TDP[MET & NI & EP2 & NEG, i]
    a_TDP[MET & NI & EP3 & NEG, MET & NI & EP3 & HCV, i] <- m_TDP[MET & NI & EP3 & NEG, i]
    
    a_TDP[MET & NI & EP1 & NEG, MET & NI & EP1 & COI, i] <- m_TDP[MET & NI & EP1 & NEG, i]
    a_TDP[MET & NI & EP2 & NEG, MET & NI & EP2 & COI, i] <- m_TDP[MET & NI & EP2 & NEG, i]
    a_TDP[MET & NI & EP3 & NEG, MET & NI & EP3 & COI, i] <- m_TDP[MET & NI & EP3 & NEG, i]
    
    a_TDP[MET & NI & EP1 & HIV, MET & NI & EP1 & COI, i] <- m_TDP[MET & NI & EP1 & NEG, i]
    a_TDP[MET & NI & EP2 & HIV, MET & NI & EP2 & COI, i] <- m_TDP[MET & NI & EP2 & NEG, i]
    a_TDP[MET & NI & EP3 & HIV, MET & NI & EP3 & COI, i] <- m_TDP[MET & NI & EP3 & NEG, i]
    
    a_TDP[MET & NI & EP1 & HCV, MET & NI & EP1 & COI, i] <- m_TDP[MET & NI & EP1 & NEG, i]
    a_TDP[MET & NI & EP2 & HCV, MET & NI & EP2 & COI, i] <- m_TDP[MET & NI & EP2 & NEG, i]
    a_TDP[MET & NI & EP3 & HCV, MET & NI & EP3 & COI, i] <- m_TDP[MET & NI & EP3 & NEG, i]
    
    # REL
    a_TDP[REL & NI & EP1 & NEG, REL & NI & EP1 & HIV, i] <- m_TDP[REL & NI & EP1 & NEG, i]
    a_TDP[REL & NI & EP2 & NEG, REL & NI & EP2 & HIV, i] <- m_TDP[REL & NI & EP2 & NEG, i]
    a_TDP[REL & NI & EP3 & NEG, REL & NI & EP3 & HIV, i] <- m_TDP[REL & NI & EP3 & NEG, i]
    
    a_TDP[REL & NI & EP1 & NEG, REL & NI & EP1 & HCV, i] <- m_TDP[REL & NI & EP1 & NEG, i]
    a_TDP[REL & NI & EP2 & NEG, REL & NI & EP2 & HCV, i] <- m_TDP[REL & NI & EP2 & NEG, i]
    a_TDP[REL & NI & EP3 & NEG, REL & NI & EP3 & HCV, i] <- m_TDP[REL & NI & EP3 & NEG, i]
    
    a_TDP[REL & NI & EP1 & NEG, REL & NI & EP1 & COI, i] <- m_TDP[REL & NI & EP1 & NEG, i]
    a_TDP[REL & NI & EP2 & NEG, REL & NI & EP2 & COI, i] <- m_TDP[REL & NI & EP2 & NEG, i]
    a_TDP[REL & NI & EP3 & NEG, REL & NI & EP3 & COI, i] <- m_TDP[REL & NI & EP3 & NEG, i]
    
    a_TDP[REL & NI & EP1 & HIV, REL & NI & EP1 & COI, i] <- m_TDP[REL & NI & EP1 & NEG, i]
    a_TDP[REL & NI & EP2 & HIV, REL & NI & EP2 & COI, i] <- m_TDP[REL & NI & EP2 & NEG, i]
    a_TDP[REL & NI & EP3 & HIV, REL & NI & EP3 & COI, i] <- m_TDP[REL & NI & EP3 & NEG, i]
    
    a_TDP[REL & NI & EP1 & HCV, REL & NI & EP1 & COI, i] <- m_TDP[REL & NI & EP1 & NEG, i]
    a_TDP[REL & NI & EP2 & HCV, REL & NI & EP2 & COI, i] <- m_TDP[REL & NI & EP2 & NEG, i]
    a_TDP[REL & NI & EP3 & HCV, REL & NI & EP3 & COI, i] <- m_TDP[REL & NI & EP3 & NEG, i]
    
    # OD
    #a_TDP[OD & NI & EP1 & NEG, OD & NI & EP1 & POS, i] <- m_TDP[OD & NI & EP1 & NEG, i]
    #a_TDP[OD & NI & EP2 & NEG, OD & NI & EP2 & POS, i] <- m_TDP[OD & NI & EP2 & NEG, i]
    #a_TDP[OD & NI & EP3 & NEG, OD & NI & EP3 & POS, i] <- m_TDP[OD & NI & EP3 & NEG, i]
    
    # ABS
    a_TDP[ABS & NI & EP1 & NEG, ABS & NI & EP1 & HIV, i] <- m_TDP[ABS & NI & EP1 & NEG, i]
    a_TDP[ABS & NI & EP2 & NEG, ABS & NI & EP2 & HIV, i] <- m_TDP[ABS & NI & EP2 & NEG, i]
    a_TDP[ABS & NI & EP3 & NEG, ABS & NI & EP3 & HIV, i] <- m_TDP[ABS & NI & EP3 & NEG, i]
    
    a_TDP[ABS & NI & EP1 & NEG, ABS & NI & EP1 & HCV, i] <- m_TDP[ABS & NI & EP1 & NEG, i]
    a_TDP[ABS & NI & EP2 & NEG, ABS & NI & EP2 & HCV, i] <- m_TDP[ABS & NI & EP2 & NEG, i]
    a_TDP[ABS & NI & EP3 & NEG, ABS & NI & EP3 & HCV, i] <- m_TDP[ABS & NI & EP3 & NEG, i]
    
    a_TDP[ABS & NI & EP1 & NEG, ABS & NI & EP1 & COI, i] <- m_TDP[ABS & NI & EP1 & NEG, i]
    a_TDP[ABS & NI & EP2 & NEG, ABS & NI & EP2 & COI, i] <- m_TDP[ABS & NI & EP2 & NEG, i]
    a_TDP[ABS & NI & EP3 & NEG, ABS & NI & EP3 & COI, i] <- m_TDP[ABS & NI & EP3 & NEG, i]
    
    a_TDP[ABS & NI & EP1 & HIV, ABS & NI & EP1 & COI, i] <- m_TDP[ABS & NI & EP1 & NEG, i]
    a_TDP[ABS & NI & EP2 & HIV, ABS & NI & EP2 & COI, i] <- m_TDP[ABS & NI & EP2 & NEG, i]
    a_TDP[ABS & NI & EP3 & HIV, ABS & NI & EP3 & COI, i] <- m_TDP[ABS & NI & EP3 & NEG, i]
    
    a_TDP[ABS & NI & EP1 & HCV, ABS & NI & EP1 & COI, i] <- m_TDP[ABS & NI & EP1 & NEG, i]
    a_TDP[ABS & NI & EP2 & HCV, ABS & NI & EP2 & COI, i] <- m_TDP[ABS & NI & EP2 & NEG, i]
    a_TDP[ABS & NI & EP3 & HCV, ABS & NI & EP3 & COI, i] <- m_TDP[ABS & NI & EP3 & NEG, i]

    # Injection
    # BUP
    a_TDP[BUP & INJ & EP1 & NEG, BUP & INJ & EP1 & HIV, i] <- m_TDP[BUP & INJ & EP1 & NEG, i]
    a_TDP[BUP & INJ & EP2 & NEG, BUP & INJ & EP2 & HIV, i] <- m_TDP[BUP & INJ & EP2 & NEG, i]
    a_TDP[BUP & INJ & EP3 & NEG, BUP & INJ & EP3 & HIV, i] <- m_TDP[BUP & INJ & EP3 & NEG, i]
    
    a_TDP[BUP & INJ & EP1 & NEG, BUP & INJ & EP1 & HCV, i] <- m_TDP[BUP & INJ & EP1 & NEG, i]
    a_TDP[BUP & INJ & EP2 & NEG, BUP & INJ & EP2 & HCV, i] <- m_TDP[BUP & INJ & EP2 & NEG, i]
    a_TDP[BUP & INJ & EP3 & NEG, BUP & INJ & EP3 & HCV, i] <- m_TDP[BUP & INJ & EP3 & NEG, i]
    
    a_TDP[BUP & INJ & EP1 & NEG, BUP & INJ & EP1 & COI, i] <- m_TDP[BUP & INJ & EP1 & NEG, i]
    a_TDP[BUP & INJ & EP2 & NEG, BUP & INJ & EP2 & COI, i] <- m_TDP[BUP & INJ & EP2 & NEG, i]
    a_TDP[BUP & INJ & EP3 & NEG, BUP & INJ & EP3 & COI, i] <- m_TDP[BUP & INJ & EP3 & NEG, i]
    
    a_TDP[BUP & INJ & EP1 & HIV, BUP & INJ & EP1 & COI, i] <- m_TDP[BUP & INJ & EP1 & NEG, i]
    a_TDP[BUP & INJ & EP2 & HIV, BUP & INJ & EP2 & COI, i] <- m_TDP[BUP & INJ & EP2 & NEG, i]
    a_TDP[BUP & INJ & EP3 & HIV, BUP & INJ & EP3 & COI, i] <- m_TDP[BUP & INJ & EP3 & NEG, i]
    
    a_TDP[BUP & INJ & EP1 & HCV, BUP & INJ & EP1 & COI, i] <- m_TDP[BUP & INJ & EP1 & NEG, i]
    a_TDP[BUP & INJ & EP2 & HCV, BUP & INJ & EP2 & COI, i] <- m_TDP[BUP & INJ & EP2 & NEG, i]
    a_TDP[BUP & INJ & EP3 & HCV, BUP & INJ & EP3 & COI, i] <- m_TDP[BUP & INJ & EP3 & NEG, i]
    
    # MET
    a_TDP[MET & INJ & EP1 & NEG, MET & INJ & EP1 & HIV, i] <- m_TDP[MET & INJ & EP1 & NEG, i]
    a_TDP[MET & INJ & EP2 & NEG, MET & INJ & EP2 & HIV, i] <- m_TDP[MET & INJ & EP2 & NEG, i]
    a_TDP[MET & INJ & EP3 & NEG, MET & INJ & EP3 & HIV, i] <- m_TDP[MET & INJ & EP3 & NEG, i]
    
    a_TDP[MET & INJ & EP1 & NEG, MET & INJ & EP1 & HCV, i] <- m_TDP[MET & INJ & EP1 & NEG, i]
    a_TDP[MET & INJ & EP2 & NEG, MET & INJ & EP2 & HCV, i] <- m_TDP[MET & INJ & EP2 & NEG, i]
    a_TDP[MET & INJ & EP3 & NEG, MET & INJ & EP3 & HCV, i] <- m_TDP[MET & INJ & EP3 & NEG, i]
    
    a_TDP[MET & INJ & EP1 & NEG, MET & INJ & EP1 & COI, i] <- m_TDP[MET & INJ & EP1 & NEG, i]
    a_TDP[MET & INJ & EP2 & NEG, MET & INJ & EP2 & COI, i] <- m_TDP[MET & INJ & EP2 & NEG, i]
    a_TDP[MET & INJ & EP3 & NEG, MET & INJ & EP3 & COI, i] <- m_TDP[MET & INJ & EP3 & NEG, i]
    
    a_TDP[MET & INJ & EP1 & HIV, MET & INJ & EP1 & COI, i] <- m_TDP[MET & INJ & EP1 & NEG, i]
    a_TDP[MET & INJ & EP2 & HIV, MET & INJ & EP2 & COI, i] <- m_TDP[MET & INJ & EP2 & NEG, i]
    a_TDP[MET & INJ & EP3 & HIV, MET & INJ & EP3 & COI, i] <- m_TDP[MET & INJ & EP3 & NEG, i]
    
    a_TDP[MET & INJ & EP1 & HCV, MET & INJ & EP1 & COI, i] <- m_TDP[MET & INJ & EP1 & NEG, i]
    a_TDP[MET & INJ & EP2 & HCV, MET & INJ & EP2 & COI, i] <- m_TDP[MET & INJ & EP2 & NEG, i]
    a_TDP[MET & INJ & EP3 & HCV, MET & INJ & EP3 & COI, i] <- m_TDP[MET & INJ & EP3 & NEG, i]
    
    # REL
    a_TDP[REL & INJ & EP1 & NEG, REL & INJ & EP1 & HIV, i] <- m_TDP[REL & INJ & EP1 & NEG, i]
    a_TDP[REL & INJ & EP2 & NEG, REL & INJ & EP2 & HIV, i] <- m_TDP[REL & INJ & EP2 & NEG, i]
    a_TDP[REL & INJ & EP3 & NEG, REL & INJ & EP3 & HIV, i] <- m_TDP[REL & INJ & EP3 & NEG, i]
    
    a_TDP[REL & INJ & EP1 & NEG, REL & INJ & EP1 & HCV, i] <- m_TDP[REL & INJ & EP1 & NEG, i]
    a_TDP[REL & INJ & EP2 & NEG, REL & INJ & EP2 & HCV, i] <- m_TDP[REL & INJ & EP2 & NEG, i]
    a_TDP[REL & INJ & EP3 & NEG, REL & INJ & EP3 & HCV, i] <- m_TDP[REL & INJ & EP3 & NEG, i]
    
    a_TDP[REL & INJ & EP1 & NEG, REL & INJ & EP1 & COI, i] <- m_TDP[REL & INJ & EP1 & NEG, i]
    a_TDP[REL & INJ & EP2 & NEG, REL & INJ & EP2 & COI, i] <- m_TDP[REL & INJ & EP2 & NEG, i]
    a_TDP[REL & INJ & EP3 & NEG, REL & INJ & EP3 & COI, i] <- m_TDP[REL & INJ & EP3 & NEG, i]
    
    a_TDP[REL & INJ & EP1 & HIV, REL & INJ & EP1 & COI, i] <- m_TDP[REL & INJ & EP1 & NEG, i]
    a_TDP[REL & INJ & EP2 & HIV, REL & INJ & EP2 & COI, i] <- m_TDP[REL & INJ & EP2 & NEG, i]
    a_TDP[REL & INJ & EP3 & HIV, REL & INJ & EP3 & COI, i] <- m_TDP[REL & INJ & EP3 & NEG, i]
    
    a_TDP[REL & INJ & EP1 & HCV, REL & INJ & EP1 & COI, i] <- m_TDP[REL & INJ & EP1 & NEG, i]
    a_TDP[REL & INJ & EP2 & HCV, REL & INJ & EP2 & COI, i] <- m_TDP[REL & INJ & EP2 & NEG, i]
    a_TDP[REL & INJ & EP3 & HCV, REL & INJ & EP3 & COI, i] <- m_TDP[REL & INJ & EP3 & NEG, i]
    
    # OD
    #a_TDP[OD & INJ & EP1 & NEG, OD & INJ & EP1 & POS, i] <- m_TDP[OD & INJ & EP1 & NEG, i]
    #a_TDP[OD & INJ & EP2 & NEG, OD & INJ & EP2 & POS, i] <- m_TDP[OD & INJ & EP2 & NEG, i]
    #a_TDP[OD & INJ & EP3 & NEG, OD & INJ & EP3 & POS, i] <- m_TDP[OD & INJ & EP3 & NEG, i]
    
    # ABS
    a_TDP[ABS & INJ & EP1 & NEG, ABS & INJ & EP1 & HIV, i] <- m_TDP[ABS & INJ & EP1 & NEG, i]
    a_TDP[ABS & INJ & EP2 & NEG, ABS & INJ & EP2 & HIV, i] <- m_TDP[ABS & INJ & EP2 & NEG, i]
    a_TDP[ABS & INJ & EP3 & NEG, ABS & INJ & EP3 & HIV, i] <- m_TDP[ABS & INJ & EP3 & NEG, i]
    
    a_TDP[ABS & INJ & EP1 & NEG, ABS & INJ & EP1 & HCV, i] <- m_TDP[ABS & INJ & EP1 & NEG, i]
    a_TDP[ABS & INJ & EP2 & NEG, ABS & INJ & EP2 & HCV, i] <- m_TDP[ABS & INJ & EP2 & NEG, i]
    a_TDP[ABS & INJ & EP3 & NEG, ABS & INJ & EP3 & HCV, i] <- m_TDP[ABS & INJ & EP3 & NEG, i]
    
    a_TDP[ABS & INJ & EP1 & NEG, ABS & INJ & EP1 & COI, i] <- m_TDP[ABS & INJ & EP1 & NEG, i]
    a_TDP[ABS & INJ & EP2 & NEG, ABS & INJ & EP2 & COI, i] <- m_TDP[ABS & INJ & EP2 & NEG, i]
    a_TDP[ABS & INJ & EP3 & NEG, ABS & INJ & EP3 & COI, i] <- m_TDP[ABS & INJ & EP3 & NEG, i]
    
    a_TDP[ABS & INJ & EP1 & HIV, ABS & INJ & EP1 & COI, i] <- m_TDP[ABS & INJ & EP1 & NEG, i]
    a_TDP[ABS & INJ & EP2 & HIV, ABS & INJ & EP2 & COI, i] <- m_TDP[ABS & INJ & EP2 & NEG, i]
    a_TDP[ABS & INJ & EP3 & HIV, ABS & INJ & EP3 & COI, i] <- m_TDP[ABS & INJ & EP3 & NEG, i]
    
    a_TDP[ABS & INJ & EP1 & HCV, ABS & INJ & EP1 & COI, i] <- m_TDP[ABS & INJ & EP1 & NEG, i]
    a_TDP[ABS & INJ & EP2 & HCV, ABS & INJ & EP2 & COI, i] <- m_TDP[ABS & INJ & EP2 & NEG, i]
    a_TDP[ABS & INJ & EP3 & HCV, ABS & INJ & EP3 & COI, i] <- m_TDP[ABS & INJ & EP3 & NEG, i]
  }

  #### Seroconversion ####
  # Apply seroconversion probability to re-weight NEG -> POS for to-states each time period
  # Probabilities applied equally across POS/NEG initially, re-weight by sero prob
  # Currently applies to HIV, need to expand state-space for HCV
  
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

  # Episode rules
  # Disallowed transitions
  a_TDP[EP1, EP3, ] = 0
  a_TDP[EP2, EP1, ] = 0
  a_TDP[EP3, EP1, ] = 0
  a_TDP[EP3, EP2, ] = 0
  a_TDP[HIV, NEG, ] = 0
  a_TDP[HCV, NEG, ] = 0 # disallowing potential transitions from COI to HIV-only (i.e. HCV cure), calculated within overall HCV infection rate
  a_TDP[COI, NEG, ] = 0
  a_TDP[HIV, HCV, ] = 0
  a_TDP[HCV, HIV, ] = 0
  a_TDP[COI, HIV, ] = 0 # disallowing potential transitions from COI to HIV-only (i.e. HCV cure), calculated within overall HCV infection rate
  a_TDP[COI, HCV, ] = 0
  a_TDP[COI, NEG, ] = 0
  a_TDP[ABS, TX, ]  = 0
  #a_TDP[BUP1, BUP1, ]  = 0
  #a_TDP[MET1, MET1, ]  = 0
  #a_TDP[REL1, REL1, ]  = 0
  # Conditional transitions
  # Next episode with out-of-treatment(OOT) EPi -> treatment(TX) EP(i+1)
  a_TDP[TX & EP1, OOT & EP2, ] = 0
  a_TDP[TX & EP2, OOT & EP3, ] = 0
  a_TDP[OOT & EP1, TX & EP1, ] = 0
  a_TDP[OOT & EP2, TX & EP2, ] = 0

  #### Check transition array ####
  check_transition_probability(a_P = a_TDP, err_stop = err_stop, verbose = verbose) # check all probs [0, 1]
  check_sum_of_transition_array(a_P = a_TDP, n_states = n_states, n_t = n_t, err_stop = err_stop, verbose = verbose) # check prob sums = 1

  #### Run Markov model ####
  # Create empty initial state vectors
  v_s_init <- rep(0, n_states)
  names(v_s_init) <- v_n_states

  #### Set initial state vector ####
  # Baseline
  # Populate first episode in base states
  #v_s_init[BUP1 & EP1] <- v_init_dist["pe", "BUP1"] # Empirically observed proportions from base states
  v_s_init[BUP & EP1]  <- v_init_dist["pe", "BUP"]
  #v_s_init[MET1 & EP1] <- v_init_dist["pe", "MET1"]
  v_s_init[MET & EP1]  <- v_init_dist["pe", "MET"]
  #v_s_init[REL1 & EP1] <- v_init_dist["pe", "REL1"]
  v_s_init[REL & EP1]  <- v_init_dist["pe", "REL"]
  v_s_init[ODN & EP1]   <- v_init_dist["pe", "ODN"]
  v_s_init[ODF & EP1]   <- v_init_dist["pe", "ODF"]
  v_s_init[ABS & EP1]  <- v_init_dist["pe", "ABS"]
  
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
    m_M_trace_death[i, ] <- m_M_trace[i - 1, ] * m_mort[, i - 1] # State-specific deaths at each time point as function of alive in t-1
  }
  m_M_trace_cumsum_death <- apply(m_M_trace_death, 2, cumsum) # Cumulative deaths at each time point (use m_M_trace_death for individual period deaths)

  #### Create aggregated trace matrices ####
  v_agg_trace_states <- c("Alive", "Death", "ODN", "ODF", "REL", "BUP", "MET", "ABS") # states to aggregate
  v_agg_trace_death_states <- c("Total", "ODN", "ODF", "REL", "BUP", "MET", "ABS") # states to aggregate
  v_agg_trace_sero_states <- c("HIV - Alive", "HIV - Dead") # states to aggregate
  
  n_agg_trace_states <- length(v_agg_trace_states)
  n_agg_trace_death_states <- length(v_agg_trace_death_states)
  n_agg_trace_sero_states <- length(v_agg_trace_sero_states)
  
  m_M_agg_trace <- array(0, dim = c((n_t + 1), n_agg_trace_states),
                         dimnames = list(0:n_t, v_agg_trace_states))
  m_M_agg_trace_death <- array(0, dim = c((n_t + 1), n_agg_trace_death_states),
                               dimnames = list(0:n_t, v_agg_trace_death_states))
  m_M_agg_trace_sero  <- array(0, dim = c((n_t + 1), n_agg_trace_sero_states),
                               dimnames = list(0:n_t, v_agg_trace_sero_states))
  
  for (i in 1:n_t){
    m_M_agg_trace[i, "Alive"] <- sum(m_M_trace[i, ])
    #m_M_agg_trace[i, "BUP1"]  <- sum(m_M_trace[i, BUP1])
    m_M_agg_trace[i, "BUP"]   <- sum(m_M_trace[i, BUP])
    #m_M_agg_trace[i, "MET1"]  <- sum(m_M_trace[i, MET1])
    m_M_agg_trace[i, "MET"]   <- sum(m_M_trace[i, MET])
    #m_M_agg_trace[i, "REL1"]  <- sum(m_M_trace[i, REL1])
    m_M_agg_trace[i, "REL"]   <- sum(m_M_trace[i, REL])
    m_M_agg_trace[i, "ABS"]   <- sum(m_M_trace[i, ABS])
    m_M_agg_trace[i, "ODN"]    <- sum(m_M_trace[i, ODN])
    m_M_agg_trace[i, "ODF"]    <- sum(m_M_trace[i, ODF])
    m_M_agg_trace[i, "Death"] <- 1 - sum(m_M_trace[i, ])
  }
  
  for (i in 1:n_t){
    m_M_agg_trace_death[i, "Total"] <- sum(m_M_trace_cumsum_death[i, ])
    m_M_agg_trace_death[i, "ODN"]    <- sum(m_M_trace_cumsum_death[i, ODN])
    m_M_agg_trace_death[i, "ODF"]    <- sum(m_M_trace_cumsum_death[i, ODF])
    #m_M_agg_trace_death[i, "REL1"]  <- sum(m_M_trace_cumsum_death[i, REL1])
    m_M_agg_trace_death[i, "REL"]   <- sum(m_M_trace_cumsum_death[i, REL])
    #m_M_agg_trace_death[i, "BUP1"]  <- sum(m_M_trace_cumsum_death[i, BUP1])
    m_M_agg_trace_death[i, "BUP"]   <- sum(m_M_trace_cumsum_death[i, BUP])
    #m_M_agg_trace_death[i, "MET1"]  <- sum(m_M_trace_cumsum_death[i, MET1])
    m_M_agg_trace_death[i, "MET"]   <- sum(m_M_trace_cumsum_death[i, MET])
    m_M_agg_trace_death[i, "ABS"]   <- sum(m_M_trace_cumsum_death[i, ABS])
  }
  
  for (i in 1:n_t){
    m_M_agg_trace_sero[i, "HIV - Alive"] <- sum(m_M_trace[i, POS])
    m_M_agg_trace_sero[i, "HIV - Dead"]  <- sum(m_M_trace_cumsum_death[i, POS])
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
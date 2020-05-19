#' Markov model
#'
#' \code{markov_model} implements the main model function.
#'
#' @param l_params_all List with all parameters
#' @param err_stop Logical variable to stop model run if transition array is invalid, if TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. Default = FALSE
#' @return 
#' a_TDP: Transition probability array
#' m_M_trace: Full markov cohort trace
#' m_M_agg_trace: Aggregated trace over base health states
#' @export
markov_model_CA <- function(l_params_all, err_stop = FALSE, verbose = FALSE){
  ### Definition:
  ##   Markov model implementation function
  ### Arguments:  
  ##   l_params_all: List with all parameters
  ##   verbose: Logical variable to indicate print out of messages
  ### Returns:
  ##   a_TDP: Transition probability array.
  ##   m_M_trace: Matrix cohort trace.
  ##   m_M_agg_trace: Aggregated trace over base health states.
  ##
  with(as.list(l_params_all), {

  #### Set up model states ####
  n_t <- (n_age_max - n_age_init) * 12 # modeling time horizon in months
  l_dim_s  <- list() # list of base states
  
  # Base health states
  BASE <- l_dim_s[[1]] <- c("MET", "ABS", "REL1", "REL", "OD")
  # Injection/non-injection stratification
  INJECT <- l_dim_s[[2]] <- c("NI", "INJ")
  # Episodes (1-3)
  EP <-  l_dim_s[[3]] <- c("1", "2", "3")
  # HIV status
  HIV <- l_dim_s[[4]] <- c("POS", "NEG")
  n_t <- (n_age_max - n_age_init) * 12
  
  df_flat <- expand.grid(l_dim_s) #combine all elements together into vector of health states
  df_flat <- rename(df_flat, BASE    = Var1, 
                             INJECT  = Var2, 
                             EP      = Var3, 
                             HIV     = Var4)

  # Create index of states to populate transition matrices
  # All treatment
  TX <- df_flat$BASE == "MET"
  
  # All out-of-treatment (incl ABS)
  OOT <- df_flat$BASE == "REL1" | df_flat$BASE == "REL" | df_flat$BASE == "OD" | df_flat$BASE == "ABS"

  # Methadone
  all_MET <- df_flat$BASE == "MET"
  MET     <- df_flat$BASE == "MET"

  
  # Relapse
  all_REL <- df_flat$BASE == "REL" | df_flat$BASE == "REL1"
  REL <- df_flat$BASE == "REL"
  REL1 <- df_flat$BASE == "REL1"

  # Overdose
  OD <- df_flat$BASE == "OD"
  
  # Abstinence
  ABS <- df_flat$BASE == "ABS"
  
  # HIV status
  NEG <- df_flat$HIV == "NEG"
  POS <- df_flat$HIV == "POS"
  
  # Injection
  INJ <- df_flat$INJECT == "INJ"
  NI <- df_flat$INJECT == "NI"
  
  # Episodes
  EP1 <- df_flat$EP == "1"
  EP2 <- df_flat$EP == "2"
  EP3 <- df_flat$EP == "3"

  df_n <- unite(df_flat, newCol) # combine columns into one data frame of all health states (8 states * 2 inj * 2 HIV * 3 Episodes)
  v_n_states <- df_n[,1] # convert df into vector
  n_states <- length(v_n_states) # total number of health states
  l_index_s  <- list(TX = TX, OOT = OOT, 
                     all_MET = all_MET, MET = MET, 
                     all_REL = all_REL, REL = REL, REL1 = REL1, 
                     OD = OD, ABS = ABS, 
                     NEG = NEG, POS = POS, 
                     INJ = INJ, NI = NI, 
                     EP1 = EP1, EP2 = EP2, EP3 = EP3)
  
  #### Time-dependent survival probabilities ####
    # Empty 2-D matrix
    m_TDP <- array(0, dim = c(n_states, n_t),
                      dimnames = list(v_n_states, 1:n_t))

    # Probability of remaining in health state
    # For REL1: p_remain = 0
    for(i in 1:n_t){
      t <- i+1 # Shift REL ahead 1 month to adjust for REL1
      # Non-injection
        #m_TDP[REL1 & NI, i] <- 0
        # Episode 1
        m_TDP[EP1 & MET & NI, i] <- as.vector(exp(p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        m_TDP[EP1 & ABS & NI, i] <- as.vector(exp(p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP1 & REL & NI, i] <- as.vector(exp(p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
        m_TDP[EP1 & OD & NI, i]  <- as.vector(exp(p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))  
        # Episode 2
        m_TDP[EP2 & MET & NI, i] <- as.vector(exp(p_frailty_MET_NI_2 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        m_TDP[EP2 & ABS & NI, i] <- as.vector(exp(p_frailty_ABS_NI_2 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP2 & REL & NI, i] <- as.vector(exp(p_frailty_REL_NI_2 * p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
        m_TDP[EP2 & OD & NI, i]  <- as.vector(exp(p_frailty_OD_NI_2  * p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))
        # Episode 3
        m_TDP[EP3 & MET & NI, i] <- as.vector(exp(p_frailty_MET_NI_3 * p_weibull_scale_MET_NI * (((i - 1)^p_weibull_shape_MET_NI) - (i^p_weibull_shape_MET_NI))))
        m_TDP[EP3 & ABS & NI, i] <- as.vector(exp(p_frailty_ABS_NI_3 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP3 & REL & NI, i] <- as.vector(exp(p_frailty_REL_NI_3 * p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
        m_TDP[EP3 & OD & NI, i]  <- as.vector(exp(p_frailty_OD_NI_3  * p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))
        
    # Injection
        #m_TDP[REL1 & INJ, i] <- 0
        # Episode 1
        m_TDP[EP1 & MET & INJ, i] <- as.vector(exp(p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        m_TDP[EP1 & ABS & INJ, i] <- as.vector(exp(p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP1 & REL & INJ, i] <- as.vector(exp(p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
        m_TDP[EP1 & OD & INJ, i]  <- as.vector(exp(p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ)))) 
        # Episode 2
        m_TDP[EP2 & MET & INJ, i] <- as.vector(exp(p_frailty_MET_INJ_2 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        m_TDP[EP2 & ABS & INJ, i] <- as.vector(exp(p_frailty_ABS_INJ_2 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP2 & REL & INJ, i] <- as.vector(exp(p_frailty_REL_INJ_2 * p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
        m_TDP[EP2 & OD & INJ, i]  <- as.vector(exp(p_frailty_OD_INJ_2  * p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ)))) 
        # Episode 3
        m_TDP[EP3 & MET & INJ, i] <- as.vector(exp(p_frailty_MET_INJ_3 * p_weibull_scale_MET_INJ * (((i - 1)^p_weibull_shape_MET_INJ) - (i^p_weibull_shape_MET_INJ))))
        m_TDP[EP3 & ABS & INJ, i] <- as.vector(exp(p_frailty_ABS_INJ_3 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP3 & REL & INJ, i] <- as.vector(exp(p_frailty_REL_INJ_3 * p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
        m_TDP[EP3 & OD & INJ, i]  <- as.vector(exp(p_frailty_OD_INJ_3  * p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ)))) 
  }

  # Probability of state-exit
  m_leave <- 1 - m_TDP
  
  write.csv(m_TDP,"C:/Users/Benjamin/Desktop/m_TDP.csv", row.names = TRUE)
  write.csv(m_leave,"C:/Users/Benjamin/Desktop/m_leave.csv", row.names = TRUE)

  #### Mortality ####
  # Monthly mortality for each age applied to 12 months, includes state-specific hr
  v_mort_NI <- function(hr = hr){
    #v_mort <- rep((1 - exp(-v_r_mort_by_age[n_age_init:(n_age_max - 1), ] * (1/12) * hr)), each = 12) # all three functions produce identical results
    #v_mort <- rep(1 - (1 - (1 - exp(-v_r_mort_by_age[n_age_init:(n_age_max - 1), ] * hr)))^(1/12), each = 12)
    v_mort_NI <- rep((1 - (1 - (1 - exp(-v_r_mort_by_age_NI[n_age_init:(n_age_max - 1), ])))^(1/12)) * hr, each = 12) # test
    return(v_mort_NI)
  }
  v_mort_INJ <- function(hr = hr){
    v_mort_INJ <- rep((1 - (1 - (1 - exp(-v_r_mort_by_age_INJ[n_age_init:(n_age_max - 1), ])))^(1/12)) * hr, each = 12) # test
    return(v_mort_INJ)
  }
  # Non-injection
  v_mort_MET_NI       <- v_mort_NI(hr = hr_MET_NI)
  v_mort_REL1_NI      <- v_mort_NI(hr = hr_REL1_NI)
  v_mort_REL_NI       <- v_mort_NI(hr = hr_REL_NI)
  v_mort_OD_NEG_NI    <- v_mort_NI(hr = hr_OD_NI)
  v_mort_OD_POS_NI    <- v_mort_NI(hr = hr_HIV_OD_NI) # REMOVE THIS AFTER CA REP
  v_mort_ABS_NEG_NI   <- v_mort_NI(hr = hr_ABS_NI)
  v_mort_ABS_POS_NI   <- v_mort_NI(hr = hr_HIV_NI)
  # Injection
  v_mort_MET_INJ      <- v_mort_INJ(hr = hr_MET_INJ)
  v_mort_REL1_INJ     <- v_mort_INJ(hr = hr_REL1_INJ)
  v_mort_REL_INJ      <- v_mort_INJ(hr = hr_REL_INJ)
  v_mort_OD_NEG_INJ   <- v_mort_INJ(hr = hr_OD_INJ)
  v_mort_OD_POS_INJ   <- v_mort_INJ(hr = hr_HIV_OD_INJ) # REMOVE THIS AFTER CA REP
  v_mort_ABS_NEG_INJ  <- v_mort_INJ(hr = hr_ABS_INJ)
  v_mort_ABS_POS_INJ  <- v_mort_INJ(hr = hr_HIV_INJ)

  # Create empty mortality matrix
  m_mort <- array(0, dim = c(n_states, n_t),
                  dimnames = list(v_n_states, 1:n_t))
  # Populate mortality matrix (monthly death probability from each state)
  for (i in 1:n_t){
    # Non-injection
    m_mort[MET & NI, i]         <- v_mort_MET_NI[i]
    m_mort[REL1 & NI, i]        <- v_mort_REL1_NI[i]
    m_mort[REL & NI, i]         <- v_mort_REL_NI[i]
    m_mort[OD & NI & NEG, i]    <- v_mort_OD_NEG_NI[i]
    m_mort[OD & NI & POS, i]    <- v_mort_OD_POS_NI[i]
    m_mort[ABS & NI & NEG, i]   <- v_mort_ABS_NEG_NI[i]
    m_mort[ABS & NI & POS, i]   <- v_mort_ABS_POS_NI[i]
    # Injection
    m_mort[MET & INJ, i]        <- v_mort_MET_INJ[i]
    m_mort[REL1 & INJ, i]       <- v_mort_REL1_INJ[i]
    m_mort[REL & INJ, i]        <- v_mort_REL_INJ[i]
    m_mort[OD & INJ & NEG, i]   <- v_mort_OD_NEG_INJ[i]
    m_mort[OD & INJ & POS, i]   <- v_mort_OD_POS_INJ[i]
    m_mort[ABS & INJ & NEG, i]  <- v_mort_ABS_NEG_INJ[i]
    m_mort[ABS & INJ & POS, i]  <- v_mort_ABS_POS_INJ[i]
  }
  
  # Checks
  write.csv(m_mort,"C:/Users/Benjamin/Desktop/m_mort.csv", row.names = TRUE)

  # Alive probability in each period
  m_alive <- 1 - m_mort

  #### Unconditional transition probabilities ####
  # Empty 2-D unconditional transition matrix (from states, to states)
  m_UP <- array(0, dim = c(n_states, n_states),
                dimnames = list(v_n_states, v_n_states))
  # Populate unconditional transition matrix
  # Non-Injection
  # From MET
  m_UP[MET & NI, ABS & NI] <- p_MET_ABS_NI
  m_UP[MET & NI, REL1 & NI] <- p_MET_REL1_NI
  m_UP[MET & NI, OD & NI] <- p_MET_OD_NI
  
  # From ABS
  m_UP[ABS & NI, REL1 & NI] <- p_ABS_REL1_NI
  m_UP[ABS & NI, OD & NI] <- p_ABS_OD_NI
  
  # From REL1
  m_UP[REL1 & NI, REL & NI] <- p_REL1_REL_NI
  m_UP[REL1 & NI, MET & NI] <- p_REL1_MET_NI
  m_UP[REL1 & NI, ABS & NI] <- p_REL1_ABS_NI
  m_UP[REL1 & NI, OD & NI] <- p_REL1_OD_NI
  
  # From REL
  m_UP[REL & NI, MET & NI] <- p_REL_MET_NI
  m_UP[REL & NI, ABS & NI] <- p_REL_ABS_NI
  m_UP[REL & NI, OD & NI] <- p_REL_OD_NI
  # From OD
  m_UP[OD & NI, MET & NI] <- p_OD_MET_NI
  m_UP[OD & NI, ABS & NI] <- p_OD_ABS_NI
  m_UP[OD & NI, REL1 & NI] <- p_OD_REL1_NI

  # Injection
  # From MET
  m_UP[MET & INJ, ABS & INJ] <- p_MET_ABS_INJ
  m_UP[MET & INJ, REL1 & INJ] <- p_MET_REL1_INJ
  m_UP[MET & INJ, OD & INJ] <- p_MET_OD_INJ
  
  # From ABS
  m_UP[ABS & INJ, REL1 & INJ] <- p_ABS_REL1_INJ
  m_UP[ABS & INJ, OD & INJ] <- p_ABS_OD_INJ
  
  # From REL1
  m_UP[REL1 & INJ, REL & INJ] <- p_REL1_REL_INJ
  m_UP[REL1 & INJ, MET & INJ] <- p_REL1_MET_INJ
  m_UP[REL1 & INJ, ABS & INJ] <- p_REL1_ABS_INJ
  m_UP[REL1 & INJ, OD & INJ] <- p_REL1_OD_INJ

  # From REL
  m_UP[REL & INJ, MET & INJ] <- p_REL_MET_INJ
  m_UP[REL & INJ, ABS & INJ] <- p_REL_ABS_INJ
  m_UP[REL & INJ, OD & INJ] <- p_REL_OD_INJ
  
  # From OD
  m_UP[OD & INJ, MET & INJ] <- p_OD_MET_INJ
  m_UP[OD & INJ, ABS & INJ] <- p_OD_ABS_INJ
  m_UP[OD & INJ, REL1 & INJ] <- p_OD_REL1_INJ

  # Apply transition rules
  # Episode rules
  # Disallowed transitions
  m_UP[EP1, EP3] = 0
  m_UP[EP2, EP1] = 0
  m_UP[EP3, EP1] = 0
  m_UP[EP3, EP2] = 0
  m_UP[POS, NEG] = 0
  m_UP[ABS, TX]  = 0
  m_UP[REL1, REL1]  = 0
  
  # Conditional transitions
  # Maintain cycles (initiate episode i+1 with OOT -> TX)
  m_UP[TX & EP1, TX & EP2]   = 0
  m_UP[TX & EP2, TX & EP3]   = 0
  m_UP[TX & EP1, OOT & EP2]  = 0
  m_UP[TX & EP2, OOT & EP3]  = 0
  m_UP[OOT & EP1, OOT & EP2] = 0
  m_UP[OOT & EP2, OOT & EP3] = 0
  m_UP[OOT & EP1, TX & EP1]  = 0
  m_UP[OOT & EP2, TX & EP2]  = 0

  # Checks
  write.csv(m_UP,"C:/Users/Benjamin/Desktop/m_UP.csv", row.names = TRUE)
  
  #### Create full time-dependent transition array ####
  # Empty 3-D array
  a_TDP <- array(0, dim = c(n_states, n_states, n_t),
                   dimnames = list(v_n_states, v_n_states, 1:n_t))
  
  # Add transitions conditional on state-exit (m_leave = 1 - remain)
  for (i in 1:n_t){
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
  for (i in 1:n_t){
    # Non-injection
    a_TDP[MET & NI & EP1 & NEG, MET & NI & EP1 & POS, i] <- m_TDP[MET & NI & EP1 & NEG, i]
    a_TDP[MET & NI & EP2 & NEG, MET & NI & EP2 & POS, i] <- m_TDP[MET & NI & EP2 & NEG, i]
    a_TDP[MET & NI & EP3 & NEG, MET & NI & EP3 & POS, i] <- m_TDP[MET & NI & EP3 & NEG, i]
    
    a_TDP[REL & NI & EP1 & NEG, REL & NI & EP1 & POS, i] <- m_TDP[REL & NI & EP1 & NEG, i]
    a_TDP[REL & NI & EP2 & NEG, REL & NI & EP2 & POS, i] <- m_TDP[REL & NI & EP2 & NEG, i]
    a_TDP[REL & NI & EP3 & NEG, REL & NI & EP3 & POS, i] <- m_TDP[REL & NI & EP3 & NEG, i]
    
    a_TDP[OD & NI & EP1 & NEG, OD & NI & EP1 & POS, i] <- m_TDP[OD & NI & EP1 & NEG, i]
    a_TDP[OD & NI & EP2 & NEG, OD & NI & EP2 & POS, i] <- m_TDP[OD & NI & EP2 & NEG, i]
    a_TDP[OD & NI & EP3 & NEG, OD & NI & EP3 & POS, i] <- m_TDP[OD & NI & EP3 & NEG, i]
    
    a_TDP[ABS & NI & EP1 & NEG, ABS & NI & EP1 & POS, i] <- m_TDP[ABS & NI & EP1 & NEG, i]
    a_TDP[ABS & NI & EP2 & NEG, ABS & NI & EP2 & POS, i] <- m_TDP[ABS & NI & EP2 & NEG, i]
    a_TDP[ABS & NI & EP3 & NEG, ABS & NI & EP3 & POS, i] <- m_TDP[ABS & NI & EP3 & NEG, i]

    # Injection

    a_TDP[MET & INJ & EP1 & NEG, MET & INJ & EP1 & POS, i] <- m_TDP[MET & INJ & EP1 & NEG, i]
    a_TDP[MET & INJ & EP2 & NEG, MET & INJ & EP2 & POS, i] <- m_TDP[MET & INJ & EP2 & NEG, i]
    a_TDP[MET & INJ & EP3 & NEG, MET & INJ & EP3 & POS, i] <- m_TDP[MET & INJ & EP3 & NEG, i]
    
    a_TDP[REL & INJ & EP1 & NEG, REL & INJ & EP1 & POS, i] <- m_TDP[REL & INJ & EP1 & NEG, i]
    a_TDP[REL & INJ & EP2 & NEG, REL & INJ & EP2 & POS, i] <- m_TDP[REL & INJ & EP2 & NEG, i]
    a_TDP[REL & INJ & EP3 & NEG, REL & INJ & EP3 & POS, i] <- m_TDP[REL & INJ & EP3 & NEG, i]
    
    a_TDP[OD & INJ & EP1 & NEG, OD & INJ & EP1 & POS, i] <- m_TDP[OD & INJ & EP1 & NEG, i]
    a_TDP[OD & INJ & EP2 & NEG, OD & INJ & EP2 & POS, i] <- m_TDP[OD & INJ & EP2 & NEG, i]
    a_TDP[OD & INJ & EP3 & NEG, OD & INJ & EP3 & POS, i] <- m_TDP[OD & INJ & EP3 & NEG, i]
    
    a_TDP[ABS & INJ & EP1 & NEG, ABS & INJ & EP1 & POS, i] <- m_TDP[ABS & INJ & EP1 & NEG, i]
    a_TDP[ABS & INJ & EP2 & NEG, ABS & INJ & EP2 & POS, i] <- m_TDP[ABS & INJ & EP2 & NEG, i]
    a_TDP[ABS & INJ & EP3 & NEG, ABS & INJ & EP3 & POS, i] <- m_TDP[ABS & INJ & EP3 & NEG, i]
  }
  write.csv(a_TDP[, ,50],"C:/Users/Benjamin/Desktop/a_TDP.csv", row.names = TRUE)
  #### Seroconversion ####
  # Apply seroconversion probability to re-weight NEG -> POS for to-states each time period
  # Probabilities applied equally across POS/NEG initially, re-weight by sero prob
  # Non-injection

  a_TDP[NEG & NI, MET & NI & NEG, ]  <- a_TDP[NEG & NI, MET & NI & NEG, ] * (1 - p_sero_MET_NI)
  a_TDP[NEG & NI, MET & NI & POS, ]  <- a_TDP[NEG & NI, MET & NI & POS, ] * p_sero_MET_NI
  a_TDP[NEG & NI, REL1 & NI & NEG, ] <- a_TDP[NEG & NI, REL1 & NI & NEG, ] * (1 - p_sero_REL1_NI)
  a_TDP[NEG & NI, REL1 & NI & POS, ] <- a_TDP[NEG & NI, REL1 & NI & POS, ] * p_sero_REL1_NI
  a_TDP[NEG & NI, REL & NI & NEG, ]  <- a_TDP[NEG & NI, REL & NI & NEG, ] * (1 - p_sero_REL_NI)
  a_TDP[NEG & NI, REL & NI & POS, ]  <- a_TDP[NEG & NI, REL & NI & POS, ] * p_sero_REL_NI
  a_TDP[NEG & NI, OD & NI & NEG, ]   <- a_TDP[NEG & NI, OD & NI & NEG, ] * (1 - p_sero_OD_NI)
  a_TDP[NEG & NI, OD & NI & POS, ]   <- a_TDP[NEG & NI, OD & NI & POS, ] * p_sero_OD_NI
  a_TDP[NEG & NI, ABS & NI & NEG, ]  <- a_TDP[NEG & NI, ABS & NI & NEG, ] * (1 - p_sero_ABS_NI)
  a_TDP[NEG & NI, ABS & NI & POS, ]  <- a_TDP[NEG & NI, ABS & NI & POS, ] * p_sero_ABS_NI

  # Injection

  a_TDP[NEG & INJ, MET & INJ & NEG, ]  <- a_TDP[NEG & INJ, MET & INJ & NEG, ] * (1 - p_sero_MET_INJ)
  a_TDP[NEG & INJ, MET & INJ & POS, ]  <- a_TDP[NEG & INJ, MET & INJ & POS, ] * p_sero_MET_INJ
  a_TDP[NEG & INJ, REL1 & INJ & NEG, ] <- a_TDP[NEG & INJ, REL1 & INJ & NEG, ] * (1 - p_sero_REL1_INJ)
  a_TDP[NEG & INJ, REL1 & INJ & POS, ] <- a_TDP[NEG & INJ, REL1 & INJ & POS, ] * p_sero_REL1_INJ
  a_TDP[NEG & INJ, REL & INJ & NEG, ]  <- a_TDP[NEG & INJ, REL & INJ & NEG, ] * (1 - p_sero_REL_INJ)
  a_TDP[NEG & INJ, REL & INJ & POS, ]  <- a_TDP[NEG & INJ, REL & INJ & POS, ] * p_sero_REL_INJ
  a_TDP[NEG & INJ, OD & INJ & NEG, ]   <- a_TDP[NEG & INJ, OD & INJ & NEG, ] * (1 - p_sero_OD_INJ)
  a_TDP[NEG & INJ, OD & INJ & POS, ]   <- a_TDP[NEG & INJ, OD & INJ & POS, ] * p_sero_OD_INJ
  a_TDP[NEG & INJ, ABS & INJ & NEG, ]  <- a_TDP[NEG & INJ, ABS & INJ & NEG, ] * (1 - p_sero_ABS_INJ)
  a_TDP[NEG & INJ, ABS & INJ & POS, ]  <- a_TDP[NEG & INJ, ABS & INJ & POS, ] * p_sero_ABS_INJ

  # Episode rules
  # Disallowed transitions
  a_TDP[EP1, EP3, ] = 0
  a_TDP[EP2, EP1, ] = 0
  a_TDP[EP3, EP1, ] = 0
  a_TDP[EP3, EP2, ] = 0
  a_TDP[POS, NEG, ] = 0
  a_TDP[ABS, TX, ]  = 0
  a_TDP[REL1, REL1, ]  = 0
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
  v_s_init[MET & EP1]  <- v_init_dist["pe", "MET"]
  v_s_init[REL1 & EP1] <- v_init_dist["pe", "REL1"]
  v_s_init[REL & EP1]  <- v_init_dist["pe", "REL"]
  v_s_init[OD & EP1]   <- v_init_dist["pe", "OD"]
  v_s_init[ABS & EP1]  <- v_init_dist["pe", "ABS"]
  
  # Distribute by injection/non-injection
  v_s_init[NI]  <- v_s_init[NI] * (1 - p_INJ)
  v_s_init[INJ] <- v_s_init[INJ] * p_INJ
  
  # Distribute HIV+/-
  v_s_init[NEG] <- v_s_init[NEG] * (1 - p_HIV_POS)
  v_s_init[POS] <- v_s_init[POS] * p_HIV_POS
  
  write.csv(v_s_init[],"C:/Users/Benjamin/Desktop/v_s_init.csv", row.names = TRUE)

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
        
        v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # diag returns probability of remaining
        
        a_M_trace[i, ,j + 1] <- v_same_state # add remain to trace
        
        diag(m_sojourn) <- 0
        
        v_new_state <- as.vector(v_current_state %*% m_sojourn) # individuals transitioning to new state
        
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
  v_agg_trace_states <- c("Alive", "Death", "OD", "REL1", "REL", "MET", "ABS") # states to aggregate
  v_agg_trace_death_states <- c("Total", "OD", "REL1", "REL", "MET", "ABS") # states to aggregate
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
    m_M_agg_trace[i, "MET"]   <- sum(m_M_trace[i, MET])
    m_M_agg_trace[i, "REL1"]  <- sum(m_M_trace[i, REL1])
    m_M_agg_trace[i, "REL"]   <- sum(m_M_trace[i, REL])
    m_M_agg_trace[i, "ABS"]   <- sum(m_M_trace[i, ABS])
    m_M_agg_trace[i, "OD"]    <- sum(m_M_trace[i, OD])
    m_M_agg_trace[i, "Death"] <- 1 - sum(m_M_trace[i, ])
  }
  
  for (i in 1:n_t){
    m_M_agg_trace_death[i, "Total"] <- sum(m_M_trace_cumsum_death[i, ])
    m_M_agg_trace_death[i, "OD"]    <- sum(m_M_trace_cumsum_death[i, OD])
    m_M_agg_trace_death[i, "REL1"]  <- sum(m_M_trace_cumsum_death[i, REL1])
    m_M_agg_trace_death[i, "REL"]   <- sum(m_M_trace_cumsum_death[i, REL])
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
#' \code{check_transition_probability} checks if individual transition probabilities are in \[0, 1\].
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

#' Check if the sum of transition probabilities equal to one. 
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
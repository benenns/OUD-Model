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
markov_model <- function(l_params_all, err_stop = FALSE, verbose = FALSE){
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
  BASE <- l_dim_s[[1]] <- c("MET1", "MET", "BUP1", "BUP", "ABS", "REL1", "REL", "OD")
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
  TX <- df_flat$BASE == "BUP" | df_flat$BASE == "BUP1" | df_flat$BASE == "MET" | df_flat$BASE == "MET1"
  
  # All out-of-treatment (incl ABS)
  OOT <- df_flat$BASE == "REL1" | df_flat$BASE == "REL" | df_flat$BASE == "OD" | df_flat$BASE == "ABS"
  
  # Buprenorphine
  all_BUP <- df_flat$BASE == "BUP" | df_flat$BASE == "BUP1"
  BUP     <- df_flat$BASE == "BUP"
  BUP1    <- df_flat$BASE == "BUP1"
  
  # Methadone
  all_MET <- df_flat$BASE == "MET" | df_flat$BASE == "MET1"
  MET     <- df_flat$BASE == "MET"
  MET1    <- df_flat$BASE == "MET1"
  
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
  
  #### Time-dependent survival probabilities ####
    # Empty 2-D matrix
    m_TDP <- array(0, dim = c(n_states, n_t),
                      dimnames = list(v_n_states, 1:n_t))

    # Probability of remaining in health state
    # For BUP1, MET1 and REL1: p_remain = 0
    for(i in 1:n_t){
      t <- i+1 # Shift BUP, MET, REL ahead 1 month to adjust for BUP1, MET1, REL1
      # Non-injection
        # Episode 1
        m_TDP[EP1 & BUP & NI, i] <- as.vector(exp(p_frailty_BUP_NI_1 * p_weibull_scale_BUP_NI * (((t - 1)^p_weibull_shape_BUP_NI) - (t^p_weibull_shape_BUP_NI)))) # vector of remain probabilities
        m_TDP[EP1 & MET & NI, i] <- as.vector(exp(p_frailty_MET_NI_1 * p_weibull_scale_MET_NI * (((t - 1)^p_weibull_shape_MET_NI) - (t^p_weibull_shape_MET_NI))))
        m_TDP[EP1 & ABS & NI, i] <- as.vector(exp(p_frailty_ABS_NI_1 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP1 & REL & NI, i] <- as.vector(exp(p_frailty_REL_NI_1 * p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
        m_TDP[EP1 & OD & NI, i]  <- as.vector(exp(p_frailty_OD_NI_1  * p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))  
        # Episode 2
        m_TDP[EP2 & BUP & NI, i] <- as.vector(exp(p_frailty_BUP_NI_2 * p_weibull_scale_BUP_NI * (((t - 1)^p_weibull_shape_BUP_NI) - (t^p_weibull_shape_BUP_NI))))
        m_TDP[EP2 & MET & NI, i] <- as.vector(exp(p_frailty_MET_NI_2 * p_weibull_scale_MET_NI * (((t - 1)^p_weibull_shape_MET_NI) - (t^p_weibull_shape_MET_NI))))
        m_TDP[EP2 & ABS & NI, i] <- as.vector(exp(p_frailty_ABS_NI_2 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP2 & REL & NI, i] <- as.vector(exp(p_frailty_REL_NI_2 * p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
        m_TDP[EP2 & OD & NI, i]  <- as.vector(exp(p_frailty_OD_NI_2  * p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))
        # Episode 3
        m_TDP[EP3 & BUP & NI, i] <- as.vector(exp(p_frailty_BUP_NI_3 * p_weibull_scale_BUP_NI * (((t - 1)^p_weibull_shape_BUP_NI) - (t^p_weibull_shape_BUP_NI))))
        m_TDP[EP3 & MET & NI, i] <- as.vector(exp(p_frailty_MET_NI_3 * p_weibull_scale_MET_NI * (((t - 1)^p_weibull_shape_MET_NI) - (t^p_weibull_shape_MET_NI))))
        m_TDP[EP3 & ABS & NI, i] <- as.vector(exp(p_frailty_ABS_NI_3 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
        m_TDP[EP3 & REL & NI, i] <- as.vector(exp(p_frailty_REL_NI_3 * p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
        m_TDP[EP3 & OD & NI, i]  <- as.vector(exp(p_frailty_OD_NI_3  * p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))
        
    # Injection
        # Episode 1
        m_TDP[EP1 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_1 * p_weibull_scale_BUP_INJ * (((t - 1)^p_weibull_shape_BUP_INJ) - (t^p_weibull_shape_BUP_INJ)))) # vector of remain probabilities 
        m_TDP[EP1 & MET & INJ, i] <- as.vector(exp(p_frailty_MET_INJ_1 * p_weibull_scale_MET_INJ * (((t - 1)^p_weibull_shape_MET_INJ) - (t^p_weibull_shape_MET_INJ))))
        m_TDP[EP1 & ABS & INJ, i] <- as.vector(exp(p_frailty_ABS_INJ_1 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP1 & REL & INJ, i] <- as.vector(exp(p_frailty_REL_INJ_1 * p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
        m_TDP[EP1 & OD & INJ, i]  <- as.vector(exp(p_frailty_OD_INJ_1  * p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ)))) 
        # Episode 2
        m_TDP[EP2 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_2 * p_weibull_scale_BUP_INJ * (((t - 1)^p_weibull_shape_BUP_INJ) - (t^p_weibull_shape_BUP_INJ))))
        m_TDP[EP2 & MET & INJ, i] <- as.vector(exp(p_frailty_MET_INJ_2 * p_weibull_scale_MET_INJ * (((t - 1)^p_weibull_shape_MET_INJ) - (t^p_weibull_shape_MET_INJ))))
        m_TDP[EP2 & ABS & INJ, i] <- as.vector(exp(p_frailty_ABS_INJ_2 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP2 & REL & INJ, i] <- as.vector(exp(p_frailty_REL_INJ_2 * p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
        m_TDP[EP2 & OD & INJ, i]  <- as.vector(exp(p_frailty_OD_INJ_2  * p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ)))) 
        # Episode 3
        m_TDP[EP3 & BUP & INJ, i] <- as.vector(exp(p_frailty_BUP_INJ_3 * p_weibull_scale_BUP_INJ * (((t - 1)^p_weibull_shape_BUP_INJ) - (t^p_weibull_shape_BUP_INJ))))
        m_TDP[EP3 & MET & INJ, i] <- as.vector(exp(p_frailty_MET_INJ_3 * p_weibull_scale_MET_INJ * (((t - 1)^p_weibull_shape_MET_INJ) - (t^p_weibull_shape_MET_INJ))))
        m_TDP[EP3 & ABS & INJ, i] <- as.vector(exp(p_frailty_ABS_INJ_3 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
        m_TDP[EP3 & REL & INJ, i] <- as.vector(exp(p_frailty_REL_INJ_3 * p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
        m_TDP[EP3 & OD & INJ, i]  <- as.vector(exp(p_frailty_OD_INJ_3  * p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ)))) 
  }

  # Probability of state-exit
  m_leave <- 1 - m_TDP

  # Mortality
  # Monthly mortality for each age applied to 12 months, includes state-specific hr
  v_mort <- function(hr = hr){
    v_mort <- rep((1 - exp(-v_r_mort_by_age[n_age_init:(n_age_max - 1), ] * (1/12) * hr)), each = 12)
    return(v_mort)
  }
  # Non-injection
  v_mort_BUP1_NI      <- v_mort(hr = hr_BUP1_NI)
  v_mort_BUP_NI       <- v_mort(hr = hr_BUP_NI)
  v_mort_MET1_NI      <- v_mort(hr = hr_MET1_NI)
  v_mort_MET_NI       <- v_mort(hr = hr_MET_NI)
  v_mort_REL1_NI      <- v_mort(hr = hr_REL1_NI)
  v_mort_REL_NI       <- v_mort(hr = hr_REL_NI)
  v_mort_OD_NI        <- v_mort(hr = hr_OD_NI)
  v_mort_ABS_NEG_NI   <- v_mort(hr = hr_ABS_NI)
  v_mort_ABS_POS_NI   <- v_mort(hr = hr_HIV_NI)
  # Injection
  v_mort_BUP1_INJ     <- v_mort(hr = hr_BUP1_INJ)
  v_mort_BUP_INJ      <- v_mort(hr = hr_BUP_INJ)
  v_mort_MET1_INJ     <- v_mort(hr = hr_MET1_INJ)
  v_mort_MET_INJ      <- v_mort(hr = hr_MET_INJ)
  v_mort_REL1_INJ     <- v_mort(hr = hr_REL1_INJ)
  v_mort_REL_INJ      <- v_mort(hr = hr_REL_INJ)
  v_mort_OD_INJ       <- v_mort(hr = hr_OD_INJ)
  v_mort_ABS_NEG_INJ  <- v_mort(hr = hr_ABS_INJ)
  v_mort_ABS_POS_INJ  <- v_mort(hr = hr_HIV_INJ)

  # Create empty mortality matrix
  m_mort <- array(0, dim = c(n_states, n_t),
                  dimnames = list(v_n_states, 1:n_t))
  # Populate mortality matrix (monthly death probability from each state)
  for (i in 1:n_t){
    # Non-injection
    m_mort[BUP1 & NI, i]        <- v_mort_BUP1_NI[i]
    m_mort[BUP & NI, i]         <- v_mort_BUP_NI[i]
    m_mort[MET1 & NI, i]        <- v_mort_MET1_NI[i]
    m_mort[MET & NI, i]         <- v_mort_MET_NI[i]
    m_mort[REL1 & NI, i]        <- v_mort_REL1_NI[i]
    m_mort[REL & NI, i]         <- v_mort_REL_NI[i]
    m_mort[OD & NI, i]          <- v_mort_OD_NI[i]
    m_mort[ABS & NI & NEG, i]   <- v_mort_ABS_NEG_NI[i]
    m_mort[ABS & NI & POS, i]   <- v_mort_ABS_POS_NI[i]
    # Injection
    m_mort[BUP1 & INJ, i]       <- v_mort_BUP1_INJ[i]
    m_mort[BUP & INJ, i]        <- v_mort_BUP_INJ[i]
    m_mort[MET1 & INJ, i]       <- v_mort_MET1_INJ[i]
    m_mort[MET & INJ, i]        <- v_mort_MET_INJ[i]
    m_mort[REL1 & INJ, i]       <- v_mort_REL1_INJ[i]
    m_mort[REL & INJ, i]        <- v_mort_REL_INJ[i]
    m_mort[OD & INJ, i]         <- v_mort_OD_INJ[i]
    m_mort[ABS & INJ & NEG, i]  <- v_mort_ABS_NEG_INJ[i]
    m_mort[ABS & INJ & POS, i]  <- v_mort_ABS_POS_INJ[i]
  }

  # Alive probability in each period
  m_alive <- 1 - m_mort

  # Alive probability for every state and time period
  #v_alive <- function(from_state, n_mort_per){
  #  v_alive_prob <- m_alive[from_state, n_mort_per]
  #return(m_alive)
  #}
  # Count deaths
  #v_death <- function(from_state, n_mort_per){
  #  v_death_prob <- m_mort[from_state, n_mort_per]
  #  return(m_mort)
  #}

  #### Unconditional transition probabilities ####
  # Empty 2-D unconditional transition matrix (from states, to states)
  m_UP <- array(0, dim = c(n_states, n_states),
                dimnames = list(v_n_states, v_n_states))
  # Populate unconditional transition matrix
  # Transitions from BUP1, MET1, REL1 adjusted to account for "remain" (i.e.  BUP1 -> BUP; MET1 -> MET; REL1 -> REL)
  # which change by episode # and are equal to 1st month remain for BUP, MET and REL
  p_BUP1_BUP_NI_1  <- exp(p_frailty_BUP_NI_1 * p_weibull_scale_BUP_NI * (((0)^p_weibull_shape_BUP_NI) - (1^p_weibull_shape_BUP_NI))) # Transition from BUP1 -> BUP = Month-1(BUP) -> Month-2(BUP)
  p_BUP1_BUP_NI_2  <- exp(p_frailty_BUP_NI_2 * p_weibull_scale_BUP_NI * (((0)^p_weibull_shape_BUP_NI) - (1^p_weibull_shape_BUP_NI)))
  p_BUP1_BUP_NI_3  <- exp(p_frailty_BUP_NI_3 * p_weibull_scale_BUP_NI * (((0)^p_weibull_shape_BUP_NI) - (1^p_weibull_shape_BUP_NI)))
  p_MET1_MET_NI_1  <- exp(p_frailty_MET_NI_1 * p_weibull_scale_MET_NI * (((0)^p_weibull_shape_MET_NI) - (1^p_weibull_shape_MET_NI))) # Transition from BUP1 -> BUP = Month-1(MET) -> Month-2(MET)
  p_MET1_MET_NI_2  <- exp(p_frailty_MET_NI_2 * p_weibull_scale_MET_NI * (((0)^p_weibull_shape_MET_NI) - (1^p_weibull_shape_MET_NI)))
  p_MET1_MET_NI_3  <- exp(p_frailty_MET_NI_3 * p_weibull_scale_MET_NI * (((0)^p_weibull_shape_MET_NI) - (1^p_weibull_shape_MET_NI)))
  p_REL1_REL_NI_1  <- exp(p_frailty_REL_NI_1 * p_weibull_scale_REL_NI * (((0)^p_weibull_shape_REL_NI) - (1^p_weibull_shape_REL_NI))) # Transition from BUP1 -> BUP = Month-1(REL) -> Month-2(REL)
  p_REL1_REL_NI_2  <- exp(p_frailty_REL_NI_2 * p_weibull_scale_REL_NI * (((0)^p_weibull_shape_REL_NI) - (1^p_weibull_shape_REL_NI)))
  p_REL1_REL_NI_3  <- exp(p_frailty_REL_NI_3 * p_weibull_scale_REL_NI * (((0)^p_weibull_shape_REL_NI) - (1^p_weibull_shape_REL_NI)))

  p_BUP1_BUP_INJ_1  <- exp(p_frailty_BUP_INJ_1 * p_weibull_scale_BUP_INJ * (((0)^p_weibull_shape_BUP_INJ) - (1^p_weibull_shape_BUP_INJ)))
  p_BUP1_BUP_INJ_2  <- exp(p_frailty_BUP_INJ_2 * p_weibull_scale_BUP_INJ * (((0)^p_weibull_shape_BUP_INJ) - (1^p_weibull_shape_BUP_INJ)))
  p_BUP1_BUP_INJ_3  <- exp(p_frailty_BUP_INJ_3 * p_weibull_scale_BUP_INJ * (((0)^p_weibull_shape_BUP_INJ) - (1^p_weibull_shape_BUP_INJ)))
  p_MET1_MET_INJ_1  <- exp(p_frailty_MET_INJ_1 * p_weibull_scale_MET_INJ * (((0)^p_weibull_shape_MET_INJ) - (1^p_weibull_shape_MET_INJ)))
  p_MET1_MET_INJ_2  <- exp(p_frailty_MET_INJ_2 * p_weibull_scale_MET_INJ * (((0)^p_weibull_shape_MET_INJ) - (1^p_weibull_shape_MET_INJ)))
  p_MET1_MET_INJ_3  <- exp(p_frailty_MET_INJ_3 * p_weibull_scale_MET_INJ * (((0)^p_weibull_shape_MET_INJ) - (1^p_weibull_shape_MET_INJ)))
  p_REL1_REL_INJ_1  <- exp(p_frailty_REL_INJ_1 * p_weibull_scale_REL_INJ * (((0)^p_weibull_shape_REL_INJ) - (1^p_weibull_shape_REL_INJ)))
  p_REL1_REL_INJ_2  <- exp(p_frailty_REL_INJ_2 * p_weibull_scale_REL_INJ * (((0)^p_weibull_shape_REL_INJ) - (1^p_weibull_shape_REL_INJ)))
  p_REL1_REL_INJ_3  <- exp(p_frailty_REL_INJ_3 * p_weibull_scale_REL_INJ * (((0)^p_weibull_shape_REL_INJ) - (1^p_weibull_shape_REL_INJ)))

  # Non-Injection
  # From BUP1
  m_UP[BUP1 & NI & EP1, BUP & NI & EP1] <- p_BUP1_BUP_NI_1
  m_UP[BUP1 & NI & EP2, BUP & NI & EP2] <- p_BUP1_BUP_NI_2
  m_UP[BUP1 & NI & EP3, BUP & NI & EP3] <- p_BUP1_BUP_NI_3
  
  m_UP[BUP1 & NI & EP1, MET1 & NI & EP1] <- p_BUP1_MET1_NI * (1 - p_BUP1_BUP_NI_1)
  m_UP[BUP1 & NI & EP2, MET1 & NI & EP2] <- p_BUP1_MET1_NI * (1 - p_BUP1_BUP_NI_2)
  m_UP[BUP1 & NI & EP3, MET1 & NI & EP3] <- p_BUP1_MET1_NI * (1 - p_BUP1_BUP_NI_3)
  
  m_UP[BUP1 & NI & EP1, ABS & NI & EP1] <- p_BUP1_ABS_NI * (1 - p_BUP1_BUP_NI_1)
  m_UP[BUP1 & NI & EP2, ABS & NI & EP2] <- p_BUP1_ABS_NI * (1 - p_BUP1_BUP_NI_2)
  m_UP[BUP1 & NI & EP3, ABS & NI & EP3] <- p_BUP1_ABS_NI * (1 - p_BUP1_BUP_NI_3)
  
  m_UP[BUP1 & NI & EP1, REL1 & NI & EP1] <- p_BUP1_REL1_NI * (1 - p_BUP1_BUP_NI_1)
  m_UP[BUP1 & NI & EP2, REL1 & NI & EP2] <- p_BUP1_REL1_NI * (1 - p_BUP1_BUP_NI_2)
  m_UP[BUP1 & NI & EP3, REL1 & NI & EP3] <- p_BUP1_REL1_NI * (1 - p_BUP1_BUP_NI_3)
  
  m_UP[BUP1 & NI & EP1, OD & NI & EP1] <- p_BUP1_OD_NI * (1 - p_BUP1_BUP_NI_1)
  m_UP[BUP1 & NI & EP2, OD & NI & EP2] <- p_BUP1_OD_NI * (1 - p_BUP1_BUP_NI_2)
  m_UP[BUP1 & NI & EP3, OD & NI & EP3] <- p_BUP1_OD_NI * (1 - p_BUP1_BUP_NI_3)
  
  # From BUP
  m_UP[BUP & NI, MET1 & NI] <- p_BUP_MET1_NI
  m_UP[BUP & NI, ABS & NI] <- p_BUP_ABS_NI
  m_UP[BUP & NI, REL1 & NI] <- p_BUP_REL1_NI
  m_UP[BUP & NI, OD & NI] <- p_BUP_OD_NI
  
  # From MET1
  m_UP[MET1 & NI & EP1, MET & NI & EP1] <- p_MET1_MET_NI_1
  m_UP[MET1 & NI & EP2, MET & NI & EP2] <- p_MET1_MET_NI_2
  m_UP[MET1 & NI & EP3, MET & NI & EP3] <- p_MET1_MET_NI_3
  
  m_UP[MET1 & NI & EP1, BUP1 & NI & EP1] <- p_MET1_BUP1_NI * (1 - p_MET1_MET_NI_1)
  m_UP[MET1 & NI & EP2, BUP1 & NI & EP2] <- p_MET1_BUP1_NI * (1 - p_MET1_MET_NI_2)
  m_UP[MET1 & NI & EP3, BUP1 & NI & EP3] <- p_MET1_BUP1_NI * (1 - p_MET1_MET_NI_3)
  
  m_UP[MET1 & NI & EP1, ABS & NI & EP1] <- p_MET1_ABS_NI * (1 - p_MET1_MET_NI_1)
  m_UP[MET1 & NI & EP2, ABS & NI & EP2] <- p_MET1_ABS_NI * (1 - p_MET1_MET_NI_2)
  m_UP[MET1 & NI & EP3, ABS & NI & EP3] <- p_MET1_ABS_NI * (1 - p_MET1_MET_NI_3)
  
  m_UP[MET1 & NI & EP1, REL1 & NI & EP1] <- p_MET1_REL1_NI * (1 - p_MET1_MET_NI_1)
  m_UP[MET1 & NI & EP2, REL1 & NI & EP2] <- p_MET1_REL1_NI * (1 - p_MET1_MET_NI_2)
  m_UP[MET1 & NI & EP3, REL1 & NI & EP3] <- p_MET1_REL1_NI * (1 - p_MET1_MET_NI_3)
  
  m_UP[MET1 & NI & EP1, OD & NI & EP1] <- p_MET1_OD_NI * (1 - p_MET1_MET_NI_1)
  m_UP[MET1 & NI & EP2, OD & NI & EP2] <- p_MET1_OD_NI * (1 - p_MET1_MET_NI_2)
  m_UP[MET1 & NI & EP3, OD & NI & EP3] <- p_MET1_OD_NI * (1 - p_MET1_MET_NI_3)
  
  # From MET
  m_UP[MET & NI, BUP1 & NI] <- p_MET_BUP1_NI
  m_UP[MET & NI, ABS & NI] <- p_MET_ABS_NI
  m_UP[MET & NI, REL1 & NI] <- p_MET_REL1_NI
  m_UP[MET & NI, OD & NI] <- p_MET_OD_NI
  
  # From ABS
  m_UP[ABS & NI, REL1 & NI] <- p_ABS_REL1_NI
  m_UP[ABS & NI, OD & NI] <- p_ABS_OD_NI
  
  # From REL1
  m_UP[REL1 & NI & EP1, REL & NI & EP1] <- p_REL1_REL_NI_1
  m_UP[REL1 & NI & EP2, REL & NI & EP2] <- p_REL1_REL_NI_2
  m_UP[REL1 & NI & EP3, REL & NI & EP3] <- p_REL1_REL_NI_3
  
  # Transitions from EP(i) OOT -> EP(i+1) TX
  m_UP[REL1 & NI & EP1, BUP1 & NI & EP2] <- p_REL1_BUP1_NI * (1 - p_REL1_REL_NI_1)
  m_UP[REL1 & NI & EP2, BUP1 & NI & EP3] <- p_REL1_BUP1_NI * (1 - p_REL1_REL_NI_2)
  m_UP[REL1 & NI & EP3, BUP1 & NI & EP3] <- p_REL1_BUP1_NI * (1 - p_REL1_REL_NI_3)
  m_UP[REL1 & NI & EP1, MET1 & NI & EP2] <- p_REL1_MET1_NI * (1 - p_REL1_REL_NI_1)
  m_UP[REL1 & NI & EP2, MET1 & NI & EP3] <- p_REL1_MET1_NI * (1 - p_REL1_REL_NI_2)
  m_UP[REL1 & NI & EP3, MET1 & NI & EP3] <- p_REL1_MET1_NI * (1 - p_REL1_REL_NI_3)
  
  m_UP[REL1 & NI & EP1, ABS & NI & EP1] <- p_REL1_ABS_NI * (1 - p_REL1_REL_NI_1)
  m_UP[REL1 & NI & EP2, ABS & NI & EP2] <- p_REL1_ABS_NI * (1 - p_REL1_REL_NI_2)
  m_UP[REL1 & NI & EP3, ABS & NI & EP3] <- p_REL1_ABS_NI * (1 - p_REL1_REL_NI_3)
  m_UP[REL1 & NI & EP1, OD & NI & EP1] <- p_REL1_OD_NI * (1 - p_REL1_REL_NI_1)
  m_UP[REL1 & NI & EP2, OD & NI & EP2] <- p_REL1_OD_NI * (1 - p_REL1_REL_NI_2)
  m_UP[REL1 & NI & EP3, OD & NI & EP3] <- p_REL1_OD_NI * (1 - p_REL1_REL_NI_3)
  
  # From REL
  m_UP[REL & NI, MET1 & NI] <- p_REL_MET1_NI
  m_UP[REL & NI, BUP1 & NI] <- p_REL_BUP1_NI
  m_UP[REL & NI, ABS & NI] <- p_REL_ABS_NI
  m_UP[REL & NI, OD & NI] <- p_REL_OD_NI
  # From OD
  m_UP[OD & NI, MET1 & NI] <- p_OD_MET1_NI
  m_UP[OD & NI, BUP1 & NI] <- p_OD_BUP1_NI
  m_UP[OD & NI, ABS & NI] <- p_OD_ABS_NI
  m_UP[OD & NI, REL1 & NI] <- p_OD_REL1_NI

  # Injection
  # From BUP1
  m_UP[BUP1 & INJ & EP1, BUP & INJ & EP1] <- p_BUP1_BUP_INJ_1
  m_UP[BUP1 & INJ & EP2, BUP & INJ & EP2] <- p_BUP1_BUP_INJ_2
  m_UP[BUP1 & INJ & EP3, BUP & INJ & EP3] <- p_BUP1_BUP_INJ_3
  
  m_UP[BUP1 & INJ & EP1, MET1 & INJ & EP1] <- p_BUP1_MET1_INJ * (1 - p_BUP1_BUP_INJ_1)
  m_UP[BUP1 & INJ & EP2, MET1 & INJ & EP2] <- p_BUP1_MET1_INJ * (1 - p_BUP1_BUP_INJ_2)
  m_UP[BUP1 & INJ & EP3, MET1 & INJ & EP3] <- p_BUP1_MET1_INJ * (1 - p_BUP1_BUP_INJ_3)
  
  m_UP[BUP1 & INJ & EP1, ABS & INJ & EP1] <- p_BUP1_ABS_INJ * (1 - p_BUP1_BUP_INJ_1)
  m_UP[BUP1 & INJ & EP2, ABS & INJ & EP2] <- p_BUP1_ABS_INJ * (1 - p_BUP1_BUP_INJ_2)
  m_UP[BUP1 & INJ & EP3, ABS & INJ & EP3] <- p_BUP1_ABS_INJ * (1 - p_BUP1_BUP_INJ_3)
  
  m_UP[BUP1 & INJ & EP1, REL1 & INJ & EP1] <- p_BUP1_REL1_INJ * (1 - p_BUP1_BUP_INJ_1)
  m_UP[BUP1 & INJ & EP2, REL1 & INJ & EP2] <- p_BUP1_REL1_INJ * (1 - p_BUP1_BUP_INJ_2)
  m_UP[BUP1 & INJ & EP3, REL1 & INJ & EP3] <- p_BUP1_REL1_INJ * (1 - p_BUP1_BUP_INJ_3)
  
  m_UP[BUP1 & INJ & EP1, OD & INJ & EP1] <- p_BUP1_OD_INJ * (1 - p_BUP1_BUP_INJ_1)
  m_UP[BUP1 & INJ & EP2, OD & INJ & EP2] <- p_BUP1_OD_INJ * (1 - p_BUP1_BUP_INJ_2)
  m_UP[BUP1 & INJ & EP3, OD & INJ & EP3] <- p_BUP1_OD_INJ * (1 - p_BUP1_BUP_INJ_3)
  
  # From BUP
  m_UP[BUP & INJ, MET1 & INJ] <- p_BUP_MET1_INJ
  m_UP[BUP & INJ, ABS & INJ] <- p_BUP_ABS_INJ
  m_UP[BUP & INJ, REL1 & INJ] <- p_BUP_REL1_INJ
  m_UP[BUP & INJ, OD & INJ] <- p_BUP_OD_INJ
  
  # From MET1
  m_UP[MET1 & INJ & EP1, MET & INJ & EP1] <- p_MET1_MET_INJ_1
  m_UP[MET1 & INJ & EP2, MET & INJ & EP2] <- p_MET1_MET_INJ_2
  m_UP[MET1 & INJ & EP3, MET & INJ & EP3] <- p_MET1_MET_INJ_3
  
  m_UP[MET1 & INJ & EP1, BUP1 & INJ & EP1] <- p_MET1_BUP1_INJ * (1 - p_MET1_MET_INJ_1)
  m_UP[MET1 & INJ & EP2, BUP1 & INJ & EP2] <- p_MET1_BUP1_INJ * (1 - p_MET1_MET_INJ_2)
  m_UP[MET1 & INJ & EP3, BUP1 & INJ & EP3] <- p_MET1_BUP1_INJ * (1 - p_MET1_MET_INJ_3)
  
  m_UP[MET1 & INJ & EP1, ABS & INJ & EP1] <- p_MET1_ABS_INJ * (1 - p_MET1_MET_INJ_1)
  m_UP[MET1 & INJ & EP2, ABS & INJ & EP2] <- p_MET1_ABS_INJ * (1 - p_MET1_MET_INJ_2)
  m_UP[MET1 & INJ & EP3, ABS & INJ & EP3] <- p_MET1_ABS_INJ * (1 - p_MET1_MET_INJ_3)
  
  m_UP[MET1 & INJ & EP1, REL1 & INJ & EP1] <- p_MET1_REL1_INJ * (1 - p_MET1_MET_INJ_1)
  m_UP[MET1 & INJ & EP2, REL1 & INJ & EP2] <- p_MET1_REL1_INJ * (1 - p_MET1_MET_INJ_2)
  m_UP[MET1 & INJ & EP3, REL1 & INJ & EP3] <- p_MET1_REL1_INJ * (1 - p_MET1_MET_INJ_3)
  
  m_UP[MET1 & INJ & EP1, OD & INJ & EP1] <- p_MET1_OD_INJ * (1 - p_MET1_MET_INJ_1)
  m_UP[MET1 & INJ & EP2, OD & INJ & EP2] <- p_MET1_OD_INJ * (1 - p_MET1_MET_INJ_2)
  m_UP[MET1 & INJ & EP3, OD & INJ & EP3] <- p_MET1_OD_INJ * (1 - p_MET1_MET_INJ_3)
  
  # From MET
  m_UP[MET & INJ, BUP1 & INJ] <- p_MET_BUP1_INJ
  m_UP[MET & INJ, ABS & INJ] <- p_MET_ABS_INJ
  m_UP[MET & INJ, REL1 & INJ] <- p_MET_REL1_INJ
  m_UP[MET & INJ, OD & INJ] <- p_MET_OD_INJ
  
  # From ABS
  m_UP[ABS & INJ, REL1 & INJ] <- p_ABS_REL1_INJ
  m_UP[ABS & INJ, OD & INJ] <- p_ABS_OD_INJ
  
  # From REL1
  m_UP[REL1 & INJ & EP1, REL & INJ & EP1] <- p_REL1_REL_INJ_1
  m_UP[REL1 & INJ & EP2, REL & INJ & EP2] <- p_REL1_REL_INJ_2
  m_UP[REL1 & INJ & EP3, REL & INJ & EP3] <- p_REL1_REL_INJ_3
  
  # Transitions from EP(i) OOT -> EP(i+1) TX
  m_UP[REL1 & INJ & EP1, BUP1 & INJ & EP2] <- p_REL1_BUP1_INJ * (1 - p_REL1_REL_INJ_1)
  m_UP[REL1 & INJ & EP2, BUP1 & INJ & EP3] <- p_REL1_BUP1_INJ * (1 - p_REL1_REL_INJ_2)
  m_UP[REL1 & INJ & EP3, BUP1 & INJ & EP3] <- p_REL1_BUP1_INJ * (1 - p_REL1_REL_INJ_3)
  m_UP[REL1 & INJ & EP1, MET1 & INJ & EP2] <- p_REL1_MET1_INJ * (1 - p_REL1_REL_INJ_1)
  m_UP[REL1 & INJ & EP2, MET1 & INJ & EP3] <- p_REL1_MET1_INJ * (1 - p_REL1_REL_INJ_2)
  m_UP[REL1 & INJ & EP3, MET1 & INJ & EP3] <- p_REL1_MET1_INJ * (1 - p_REL1_REL_INJ_3)
  
  m_UP[REL1 & INJ & EP1, ABS & INJ & EP1] <- p_REL1_ABS_INJ * (1 - p_REL1_REL_INJ_1)
  m_UP[REL1 & INJ & EP2, ABS & INJ & EP2] <- p_REL1_ABS_INJ * (1 - p_REL1_REL_INJ_2)
  m_UP[REL1 & INJ & EP3, ABS & INJ & EP3] <- p_REL1_ABS_INJ * (1 - p_REL1_REL_INJ_3)
  m_UP[REL1 & INJ & EP1, OD & INJ & EP1] <- p_REL1_OD_INJ * (1 - p_REL1_REL_INJ_1)
  m_UP[REL1 & INJ & EP2, OD & INJ & EP2] <- p_REL1_OD_INJ * (1 - p_REL1_REL_INJ_2)
  m_UP[REL1 & INJ & EP3, OD & INJ & EP3] <- p_REL1_OD_INJ * (1 - p_REL1_REL_INJ_3)
  
  # From REL
  m_UP[REL & INJ, MET1 & INJ] <- p_REL_MET1_INJ
  m_UP[REL & INJ, BUP1 & INJ] <- p_REL_BUP1_INJ
  m_UP[REL & INJ, ABS & INJ] <- p_REL_ABS_INJ
  m_UP[REL & INJ, OD & INJ] <- p_REL_OD_INJ
  
  # From OD
  m_UP[OD & INJ, MET1 & INJ] <- p_OD_MET1_INJ
  m_UP[OD & INJ, BUP1 & INJ] <- p_OD_BUP1_INJ
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
  m_UP[BUP1, BUP1]  = 0
  m_UP[MET1, MET1]  = 0
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
  #write.csv(m_UP,"C:/Users/Benjamin/Desktop/m_UP.csv", row.names = TRUE)
  
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
    a_TDP[BUP & NI & EP1 & NEG, BUP & NI & EP1 & POS, i] <- m_TDP[BUP & NI & EP1 & NEG, i]
    a_TDP[BUP & NI & EP2 & NEG, BUP & NI & EP2 & POS, i] <- m_TDP[BUP & NI & EP2 & NEG, i]
    a_TDP[BUP & NI & EP3 & NEG, BUP & NI & EP3 & POS, i] <- m_TDP[BUP & NI & EP3 & NEG, i]
    
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
    a_TDP[BUP & INJ & EP1 & NEG, BUP & INJ & EP1 & POS, i] <- m_TDP[BUP & INJ & EP1 & NEG, i]
    a_TDP[BUP & INJ & EP2 & NEG, BUP & INJ & EP2 & POS, i] <- m_TDP[BUP & INJ & EP2 & NEG, i]
    a_TDP[BUP & INJ & EP3 & NEG, BUP & INJ & EP3 & POS, i] <- m_TDP[BUP & INJ & EP3 & NEG, i]
    
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

  #### Seroconversion ####
  # Apply seroconversion probability to re-weight NEG -> POS for to-states each time period
  # Probabilities applied equally across POS/NEG initially, re-weight by sero prob
  # Non-injection
  a_TDP[NEG & NI, BUP1 & NI & NEG, ] <- a_TDP[NEG & NI, BUP1 & NI & NEG, ] * (1 - p_sero_BUP1_NI)
  a_TDP[NEG & NI, BUP1 & NI & POS, ] <- a_TDP[NEG & NI, BUP1 & NI & POS, ] * p_sero_BUP1_NI
  a_TDP[NEG & NI, BUP & NI & NEG, ]  <- a_TDP[NEG & NI, BUP & NI & NEG, ] * (1 - p_sero_BUP_NI)
  a_TDP[NEG & NI, BUP & NI & POS, ]  <- a_TDP[NEG & NI, BUP & NI & POS, ] * p_sero_BUP_NI
  a_TDP[NEG & NI, MET1 & NI & NEG, ] <- a_TDP[NEG & NI, MET1 & NI & NEG, ] * (1 - p_sero_MET1_NI)
  a_TDP[NEG & NI, MET1 & NI & POS, ] <- a_TDP[NEG & NI, MET1 & NI & POS, ] * p_sero_MET1_NI
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
  a_TDP[NEG & INJ, BUP1 & INJ & NEG, ] <- a_TDP[NEG & INJ, BUP1 & INJ & NEG, ] * (1 - p_sero_BUP1_INJ)
  a_TDP[NEG & INJ, BUP1 & INJ & POS, ] <- a_TDP[NEG & INJ, BUP1 & INJ & POS, ] * p_sero_BUP1_INJ
  a_TDP[NEG & INJ, BUP & INJ & NEG, ]  <- a_TDP[NEG & INJ, BUP & INJ & NEG, ] * (1 - p_sero_BUP_INJ)
  a_TDP[NEG & INJ, BUP & INJ & POS, ]  <- a_TDP[NEG & INJ, BUP & INJ & POS, ] * p_sero_BUP_INJ
  a_TDP[NEG & INJ, MET1 & INJ & NEG, ] <- a_TDP[NEG & INJ, MET1 & INJ & NEG, ] * (1 - p_sero_MET1_INJ)
  a_TDP[NEG & INJ, MET1 & INJ & POS, ] <- a_TDP[NEG & INJ, MET1 & INJ & POS, ] * p_sero_MET1_INJ
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
  a_TDP[BUP1, BUP1, ]  = 0
  a_TDP[MET1, MET1, ]  = 0
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
  v_s_init_BL <- v_s_init_BUP <- v_s_init_MET <- rep(0, n_states)
  names(v_s_init_BL) <- names(v_s_init_BUP) <- names(v_s_init_MET) <- v_n_states

  # Set initial state vector
  # Baseline
  v_s_init_BL[BUP1 & NI & EP1 & NEG]  <- 0.5 # Empirically observed proportions from base state
  v_s_init_BL[MET1 & NI & EP1 & NEG]  <- 0.5
  
  # BUP
  v_s_init_BUP[BUP] = 1/sum(BUP) # Set all BUP states equal, and sum to 1
  # MET
  v_s_init_MET[MET] = 1/sum(MET) # Set all MET states equal, and sum to 1

  # Create Markov Trace
    # Initialize population
    a_M_trace <- array(0, dim = c((n_t + 1), n_states, (n_t + 1)),
                       dimnames = list(0:n_t, v_n_states, 0:n_t))
    a_M_trace[1, , 1] <- v_s_init_BL

    # All model time periods
      for(i in 2:(n_t)){ #added + 1
      # Time spent in given health state
      for(j in 1:(i - 1)){
        #state-time-dependent transition probability (j) * age (model-time)-specific mortality (i)
        m_sojourn <- a_TDP[, , j] * m_alive[, i]
        
        v_current_state <- as.vector(a_M_trace[i - 1, , j])
        
        v_same_state <- as.vector(v_current_state * diag(m_sojourn))
        
        a_M_trace[i, ,j + 1] <- v_same_state
        
        diag(m_sojourn) <- 0
        
        v_new_state <- as.vector(v_current_state %*% m_sojourn)
        
        a_M_trace[i, ,1] <- v_new_state + a_M_trace[i, ,1]
        }
  }
  m_M_trace <- array(0, dim = c((n_t + 1), n_states),
                     dimnames = list(0:n_t, v_n_states))
  for (i in 1:n_t){
    m_M_trace[i, ] <- rowSums(a_M_trace[i, ,])
  }

  #### Create aggregated trace matrix ####
  v_agg_trace_states <- c("Alive", "Death", "OD", "REL1", "REL", "BUP1", "BUP", "MET1", "MET", "ABS") # states to aggregate
  n_agg_trace_states <- length(v_agg_trace_states)
  m_M_agg_trace <- array(0, dim = c((n_t + 1), n_agg_trace_states),
                         dimnames = list(0:n_t, v_agg_trace_states))
  
  for (i in 1:n_t){
    m_M_agg_trace[i, "Alive"] <- sum(m_M_trace[i, ])
    m_M_agg_trace[i, "BUP1"]  <- sum(m_M_trace[i, BUP1])
    m_M_agg_trace[i, "BUP"]   <- sum(m_M_trace[i, BUP])
    m_M_agg_trace[i, "MET1"]  <- sum(m_M_trace[i, MET1])
    m_M_agg_trace[i, "MET"]   <- sum(m_M_trace[i, MET])
    m_M_agg_trace[i, "REL1"]  <- sum(m_M_trace[i, REL1])
    m_M_agg_trace[i, "REL"]   <- sum(m_M_trace[i, REL])
    m_M_agg_trace[i, "ABS"]   <- sum(m_M_trace[i, ABS])
    m_M_agg_trace[i, "OD"]    <- sum(m_M_trace[i, OD])
    m_M_agg_trace[i, "Death"] <- 1 - sum(m_M_trace[i, ])
  }
  #write.csv(m_M_agg_trace,"C:/Users/Benjamin/Desktop/trace.csv", row.names = TRUE)
  return(list(a_TDP = a_TDP,
              m_M_trace = m_M_trace,
              m_M_agg_trace = m_M_agg_trace))
  }
 )
}

#' Check if transition array is valid
#'
#' \code{check_transition_probability} checks if transition probabilities are in \[0, 1\].
#'
#' @param a_P A transition probability array.
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. 
#' Default = FALSE
#'
#' @return
#' This function stops if transition probability array is not valid and shows 
#' what are the entries that are not valid
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
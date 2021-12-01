#  BUP1/BUP/MET1/MET episodes
# Episode 1
EP1_BUP <- df_flat$BASE == "BUP" & df_flat$EP == "1"
EP1_MET <- df_flat$BASE == "MET" & df_flat$EP == "1"
EP1_TX <- df_flat$BASE == "BUP1" & df_flat$EP == "1" | df_flat$BASE == "BUP" & df_flat$EP == "1" | df_flat$BASE == "MET1" & df_flat$EP == "1" | df_flat$BASE == "MET" & df_flat$EP == "1"
# Episode 2
EP2_BUP <- df_flat$BASE == "BUP" & df_flat$EP == "2"
EP2_MET <- df_flat$BASE == "MET" & df_flat$EP == "2"
EP2_TX <- df_flat$BASE == "BUP1" & df_flat$EP == "2" | df_flat$BASE == "BUP" & df_flat$EP == "2" | df_flat$BASE == "MET1" & df_flat$EP == "2" | df_flat$BASE == "MET" & df_flat$EP == "2"
# Episode 3
EP3_BUP <- df_flat$BASE == "BUP" & df_flat$EP == "3"
EP3_TX <- df_flat$BASE == "BUP1" & df_flat$EP == "3" | df_flat$BASE == "BUP" & df_flat$EP == "3" | df_flat$BASE == "MET1" & df_flat$EP == "3" | df_flat$BASE == "MET" & df_flat$EP == "3"

#  REL1/REL episodes
# Episode 1
EP1_REL <- df_flat$BASE == "REL1" & df_flat$EP == "1" | df_flat$BASE == "REL" & df_flat$EP == "1" 
# Episode 2
EP2_REL <- df_flat$BASE == "REL1" & df_flat$EP == "2" | df_flat$BASE == "REL" & df_flat$EP == "2" 
# Episode 3
EP3_REL <- df_flat$BASE == "REL1" & df_flat$EP == "3" | df_flat$BASE == "REL" & df_flat$EP == "3" 

#  OD episodes
# Episode 1
EP1_OD <- df_flat$BASE == "OD" & df_flat$EP == "1"
# Episode 2
EP2_OD <- df_flat$BASE == "OD" & df_flat$EP == "2"
# Episode 3
EP3_OD <- df_flat$BASE == "OD" & df_flat$EP == "3"




# Time-dependent survival probabilities
for(i in 1:n_t){
  #if(length(grep(("AR1"), From.baseline, fixed = FALSE))==0){
  # Time-independent
  #m_TD_remain[ i, BUP1] <- 0
  #m_TD_remain[ i, MET1] <- 0
  #m_TD_remain[ i, REL1] <- 0
  m_TD_remain [ i, D] <- 1
  
  # Episode 1
  # Shift ahead 1 period to account for BUP1; MET1; REL1
  m_TD_remain[i, EP1 & BUP] <- exp(v_frailty["BUP_1"]*v_wb_scale["BUP_1"]*((((i+1)-1)^v_wb_shape["BUP_1"])-((i+1)^v_wb_shape["BUP_1"]))) # (survival curve at time i)/(survival curve at time i-1) 
  m_TD_remain[i, EP1 & MET] <- exp(v_frailty["MET_1"]*v_wb_scale["MET_1"]*((((i+1)-1)^v_wb_shape["MET_1"])-((i+1)^v_wb_shape["MET_1"])))
  m_TD_remain[i, EP1 & REL] <- exp(v_frailty["REL_1"]*v_wb_scale["REL_1"]*((((i+1)-1)^v_wb_shape["REL_1"])-((i+1)^v_wb_shape["REL_1"])))
  
  m_TD_remain[i, EP1 & OD]    <- exp(v_frailty["OD_1"]*v_wb_scale["OD_1"]*((((i-1)^v_wb_shape["OD_1"])-(i^v_wb_shape["OD_1"]))))
  m_TD_remain[i, EP1 & ABS]   <- exp(v_frailty["ABS_1"]*v_wb_scale["ABS_1"]*((((i-1)^v_wb_shape["ABS_1"])-(i^v_wb_shape["ABS_1"]))))
  
  # Episode 2
  # Shift ahead 1 period to account for BUP1; MET1; REL1
  m_TD_remain[i, EP2 & BUP] <- exp(v_frailty["BUP_2"]*v_wb_scale["BUP_2"]*((((i+1)-1)^v_wb_shape["BUP_2"])-((i+1)^v_wb_shape["BUP_2"]))) # (survival curve at time i)/(survival curve at time i-1) 
  m_TD_remain[i, EP2 & MET] <- exp(v_frailty["MET_2"]*v_wb_scale["MET_2"]*((((i+1)-1)^v_wb_shape["MET_2"])-((i+1)^v_wb_shape["MET_2"])))
  m_TD_remain[i, EP2 & REL] <- exp(v_frailty["REL_2"]*v_wb_scale["REL_2"]*((((i+1)-1)^v_wb_shape["REL_2"])-((i+1)^v_wb_shape["REL_2"])))
  
  m_TD_remain[i, EP2 & OD]    <- exp(v_frailty["OD_2"]*v_wb_scale["OD_2"]*((((i-1)^v_wb_shape["OD_2"])-(i^v_wb_shape["OD_2"]))))
  m_TD_remain[i, EP2 & ABS]   <- exp(v_frailty["ABS_2"]*v_wb_scale["ABS_2"]*((((i-1)^v_wb_shape["ABS_2"])-(i^v_wb_shape["ABS_2"]))))
  
  # Episode 3
  # Shift ahead 1 period to account for BUP1; MET1; REL1
  m_TD_remain[i, EP3 & BUP] <- exp(v_frailty["BUP_3"]*v_wb_scale["BUP_3"]*((((i+1)-1)^v_wb_shape["BUP_3"])-((i+1)^v_wb_shape["BUP_3"]))) # (survival curve at time i)/(survival curve at time i-1) 
  m_TD_remain[i, EP3 & MET] <- exp(v_frailty["MET_3"]*v_wb_scale["MET_3"]*((((i+1)-1)^v_wb_shape["MET_3"])-((i+1)^v_wb_shape["MET_3"])))
  m_TD_remain[i, EP3 & REL] <- exp(v_frailty["REL_3"]*v_wb_scale["REL_3"]*((((i+1)-1)^v_wb_shape["REL_3"])-((i+1)^v_wb_shape["REL_3"])))
  
  m_TD_remain[i, EP3 & OD]    <- exp(v_frailty["OD_3"]*v_wb_scale["OD_3"]*((((i-1)^v_wb_shape["OD_3"])-(i^v_wb_shape["OD_3"]))))
  m_TD_remain[i, EP3 & ABS]   <- exp(v_frailty["ABS_3"]*v_wb_scale["ABS_3"]*((((i-1)^v_wb_shape["ABS_3"])-(i^v_wb_shape["ABS_3"]))))
}   
m_TD_remain      




#######################################################################################################################################
# List of impossible transitions
# Any "INJ" -> "NI" - Not allowing switch from injection to non-injection
# Any "NI" -> "INJ" - Not allowing switch from injection to non-injection (for now)
# Any "POS" -> Any "NEG" - No transition from HIV+ to HIV-
# Any "REL1", "REL", "OD" -> Any "ABS" - Can only access ABS from any treatment staten (modify in SA)
# Any "D" -> Any state - No zombies
# Any state other than "MET1" -> "MET"
# Any state other than "BUP1" -> "BUP"
# Any state other than "REL1" -> "REL"
# Same state "1" -> Same state "2/3" - Episode can only increase when transitioning from another state (think about whether to have special rule for BUP and MET to keep episode consistent)
# Same state "2" -> Same state "1/3"
# Same state "3" -> Same state "1/2"
#######################################################################################################################################

################################################################
### Time-dependent probability of remaining in health states ###
################################################################
m_TD_remain <- matrix(0, 
                      nrow = (n_t + 1), ncol = n_states, 
                      dimnames = list(0:n_t, v_n))
# Time-dependent survival probabilities
for(i in 1:n_t){
  #if(length(grep(("AR1"), From.baseline, fixed = FALSE))==0){
  # Time-independent
  #m_TD_remain[ i, BUP1] <- 0
  #m_TD_remain[ i, MET1] <- 0
  #m_TD_remain[ i, REL1] <- 0
  m_TD_remain [ i, D] <- 1
  
  # Episode 1
  # Shift ahead 1 period to account for BUP1; MET1; REL1
  m_TD_remain[i, EP1 & BUP] <- exp(v_frailty["BUP_1"]*v_wb_scale["BUP_1"]*((((i+1)-1)^v_wb_shape["BUP_1"])-((i+1)^v_wb_shape["BUP_1"]))) # (survival curve at time i)/(survival curve at time i-1) 
  m_TD_remain[i, EP1 & MET] <- exp(v_frailty["MET_1"]*v_wb_scale["MET_1"]*((((i+1)-1)^v_wb_shape["MET_1"])-((i+1)^v_wb_shape["MET_1"])))
  m_TD_remain[i, EP1 & REL] <- exp(v_frailty["REL_1"]*v_wb_scale["REL_1"]*((((i+1)-1)^v_wb_shape["REL_1"])-((i+1)^v_wb_shape["REL_1"])))
  
  m_TD_remain[i, EP1 & OD]    <- exp(v_frailty["OD_1"]*v_wb_scale["OD_1"]*((((i-1)^v_wb_shape["OD_1"])-(i^v_wb_shape["OD_1"]))))
  m_TD_remain[i, EP1 & ABS]   <- exp(v_frailty["ABS_1"]*v_wb_scale["ABS_1"]*((((i-1)^v_wb_shape["ABS_1"])-(i^v_wb_shape["ABS_1"]))))
  
  # Episode 2
  # Shift ahead 1 period to account for BUP1; MET1; REL1
  m_TD_remain[i, EP2 & BUP] <- exp(v_frailty["BUP_2"]*v_wb_scale["BUP_2"]*((((i+1)-1)^v_wb_shape["BUP_2"])-((i+1)^v_wb_shape["BUP_2"]))) # (survival curve at time i)/(survival curve at time i-1) 
  m_TD_remain[i, EP2 & MET] <- exp(v_frailty["MET_2"]*v_wb_scale["MET_2"]*((((i+1)-1)^v_wb_shape["MET_2"])-((i+1)^v_wb_shape["MET_2"])))
  m_TD_remain[i, EP2 & REL] <- exp(v_frailty["REL_2"]*v_wb_scale["REL_2"]*((((i+1)-1)^v_wb_shape["REL_2"])-((i+1)^v_wb_shape["REL_2"])))
  
  m_TD_remain[i, EP2 & OD]    <- exp(v_frailty["OD_2"]*v_wb_scale["OD_2"]*((((i-1)^v_wb_shape["OD_2"])-(i^v_wb_shape["OD_2"]))))
  m_TD_remain[i, EP2 & ABS]   <- exp(v_frailty["ABS_2"]*v_wb_scale["ABS_2"]*((((i-1)^v_wb_shape["ABS_2"])-(i^v_wb_shape["ABS_2"]))))
  
  # Episode 3
  # Shift ahead 1 period to account for BUP1; MET1; REL1
  m_TD_remain[i, EP3 & BUP] <- exp(v_frailty["BUP_3"]*v_wb_scale["BUP_3"]*((((i+1)-1)^v_wb_shape["BUP_3"])-((i+1)^v_wb_shape["BUP_3"]))) # (survival curve at time i)/(survival curve at time i-1) 
  m_TD_remain[i, EP3 & MET] <- exp(v_frailty["MET_3"]*v_wb_scale["MET_3"]*((((i+1)-1)^v_wb_shape["MET_3"])-((i+1)^v_wb_shape["MET_3"])))
  m_TD_remain[i, EP3 & REL] <- exp(v_frailty["REL_3"]*v_wb_scale["REL_3"]*((((i+1)-1)^v_wb_shape["REL_3"])-((i+1)^v_wb_shape["REL_3"])))
  
  m_TD_remain[i, EP3 & OD]    <- exp(v_frailty["OD_3"]*v_wb_scale["OD_3"]*((((i-1)^v_wb_shape["OD_3"])-(i^v_wb_shape["OD_3"]))))
  m_TD_remain[i, EP3 & ABS]   <- exp(v_frailty["ABS_3"]*v_wb_scale["ABS_3"]*((((i-1)^v_wb_shape["ABS_3"])-(i^v_wb_shape["ABS_3"]))))
}   
m_TD_remain      


v_ %>% left_join(v_state_rules)







# HIV Negative
# ABS (same as baseline)
#v_p_ABS_D_age_NEG_yr   <- (1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_age_max - 1)])) # Yearly probability
v_p_ABS_D_age_NEG   <- (1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_age_max - 1)] * (1/12)))
test <- rep(v_p_ABS_D_age_NEG, each = 12)
test
test1 <- rep(n_age_init:(n_age_max - 1), each = 12)
test1
# BUP
v_p_BUP1_D_age_NEG  <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_age_max - 1)] * (1/12) * hr_BUP1)
v_p_BUP_D_age_NEG   <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_age_max - 1)] * (1/12) * hr_BUP)
# MET
v_p_MET1_D_age_NEG  <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_age_max - 1)] * (1/12) * hr_MET1)    
v_p_MET_D_age_NEG   <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_age_max - 1)] * (1/12) * hr_MET)
#REL
v_p_REL1_D_age_NEG  <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_age_max - 1)] * (1/12) * hr_REL1)    
v_p_REL_D_age_NEG   <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_age_max - 1)] * (1/12) * hr_REL)
#OD
v_p_OD_D_age_NEG    <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_age_max - 1)] * (1/12) * hr_OD)






# General transition rules
# Fill transition array first, then adjust for impossible transitions
# Death rules
a_P[D, , ] = 0 # From death  = 0; remain in death across all strata and all time points
# HIV rules
a_P[POS, NEG, ] = 0 # No transition from HIV+ to HIV-
# Injection rules
a_P[INJ, NI] = 0
a_P[NI, INJ] = 0

# Episode rules
# Disallowed transitions
a_P[EP1, EP3, ] = 0
a_P[EP2, EP1, ] = 0
a_P[EP3, EP1, ] = 0
a_P[EP3, EP2, ] = 0

# Conditional transitions
# Maintain cycles (initiate cycle 2 with OOT -> TX)
a_P[TX & EP1, OOT & EP2, ] = 0
a_P[TX & EP2, OOT & EP3, ] = 0
a_P[OOT & EP1, TX & EP1, ] = 0
a_P[OOT & EP2, TX & EP2, ] = 0

# Remain
a_P[BUP, BUP, ] <- a_P[BUP1, BUP, ] <- p_BUP_BUP

#Examples
a_P[NI, INJ, ] = 0.5 # set any cells with NI -> INJ to 0.5
a_P[NI & POS, INJ & POS, ] = 0.25 #


# Create empty mortality matrix (from_states X n_periods)
# CAN PROBABLY DROP
m_mort <- array(0, dim = c(n_states_from, n_t),
                dimnames = list(v_n_from, 0:(n_t - 1)))

for (i in 1:n_t){
  m_mort[BUP1, i] <- v_mort_BUP1_NEG[i]
  m_mort[BUP, i]  <- v_mort_BUP_NEG[i]
  m_mort[MET1, i] <- v_mort_MET1_NEG[i]
  m_mort[MET, i]  <- v_mort_MET_NEG[i]
  m_mort[REL1, i] <- v_mort_REL1_NEG[i]
  m_mort[REL, i]  <- v_mort_REL_NEG[i]
  m_mort[OD, i]   <- v_mort_OD_NEG[i]
  m_mort[ABS, i]  <- v_mort_ABS_NEG[i]
}


# Add death probabilities
for(i in 1:n_t){
  # Non-injection
  a_P[BUP1 & NI, "D", i] <- v_mort_BUP1_NI_NEG[i]
  a_P[BUP & NI, "D", i]  <- v_mort_BUP_NI_NEG[i] 
  a_P[MET1 & NI, "D", i] <- v_mort_MET1_NI_NEG[i]
  a_P[MET & NI, "D", i]  <- v_mort_MET_NI_NEG[i] 
  a_P[REL1 & NI, "D", i] <- v_mort_REL1_NI_NEG[i]
  a_P[REL & NI, "D", i]  <- v_mort_REL_NI_NEG[i] 
  a_P[OD & NI, "D", i]   <- v_mort_OD_NI_NEG[i]  
  a_P[ABS & NI, "D", i]  <- v_mort_ABS_NI_NEG[i] 
  
  #Injection
  a_P[BUP1 & INJ, "D", i] <- v_mort_BUP1_INJ_NEG[i]
  a_P[BUP & INJ, "D", i]  <- v_mort_BUP_INJ_NEG[i] 
  a_P[MET1 & INJ, "D", i] <- v_mort_MET1_INJ_NEG[i]
  a_P[MET & INJ, "D", i]  <- v_mort_MET_INJ_NEG[i] 
  a_P[REL1 & INJ, "D", i] <- v_mort_REL1_INJ_NEG[i]
  a_P[REL & INJ, "D", i]  <- v_mort_REL_INJ_NEG[i] 
  a_P[OD & INJ, "D", i]   <- v_mort_OD_INJ_NEG[i]  
  a_P[ABS & INJ, "D", i]  <- v_mort_ABS_INJ_NEG[i] 
}



# Add remain probabilities
for(i in 1:n_t){
  
}

#############################################
### Time-dependent survival probabilities ###
#############################################
# Empty 2-D matrix
m_TDP <- array(0, dim = c(n_states, n_t),
               dimnames = list(v_n_states, 0:(n_t - 1)))

# Create dummy data
m_frailty <- array(0.4, dim = c(n_EP, n_BASE, n_INJECT),
                   dimnames = list(EP, BASE, INJECT))
m_weibull_scale <- array(0.7, dim = c(n_BASE, n_INJECT),
                         dimnames = list(BASE, INJECT))
m_weibull_shape <- array(0.7, dim = c(n_BASE, n_INJECT),
                         dimnames = list(BASE, INJECT))


# Probability of remaining in given health state
for(i in 1:n_t){
  t <- i+1 # Shift BUP, MET, REL ahead 1 month to adjust for BUP1, MET1, REL1
  # Non-injection
  # Episode 1
  m_TDP[EP1 & BUP & NI, i] <- v_TDP_BUP_NI_1[i] # vector of remain probabilities 
  m_TDP[EP1 & MET & NI, i] <- v_TDP_MET_NI_1[i]
  m_TDP[EP1 & ABS & NI, i] <- v_TDP_ABS_NI_1[i]
  m_TDP[EP1 & REL & NI, i] <- v_TDP_REL_NI_1[i]
  m_TDP[EP1 & OD & NI, i]  <- v_TDP_OD_NI_1[i]  
  # Episode 2
  v_TDP_BUP_NI_2[i] <- m_TDP[EP2 & BUP & NI, i] <- exp(m_frailty["2", "BUP", "NI"] * m_weibull_scale["BUP", "NI"] * (((t-1)^m_weibull_shape["BUP", "NI"]) - (t^m_weibull_shape["BUP", "NI"]))) # (survival curve at time i)/(survival curve at time i-1) 
  v_TDP_MET_NI_2[i] <- m_TDP[EP2 & MET & NI, i] <- exp(m_frailty["2", "MET", "NI"] * m_weibull_scale["MET", "NI"] * (((t-1)^m_weibull_shape["MET", "NI"]) - (t^m_weibull_shape["MET", "NI"])))
  v_TDP_ABS_NI_2[i] <- m_TDP[EP2 & ABS & NI, i] <- exp(m_frailty["2", "ABS", "NI"] * m_weibull_scale["ABS", "NI"] * (((i-1)^m_weibull_shape["ABS", "NI"]) - (i^m_weibull_shape["ABS", "NI"])))
  v_TDP_REL_NI_2[i] <- m_TDP[EP2 & REL & NI, i] <- exp(m_frailty["2", "REL", "NI"] * m_weibull_scale["REL", "NI"] * (((t-1)^m_weibull_shape["REL", "NI"]) - (t^m_weibull_shape["REL", "NI"])))
  v_TDP_OD_NI_2[i]  <- m_TDP[EP2 & OD & NI, i]  <- exp(m_frailty["2", "OD" , "NI"] * m_weibull_scale["OD", "NI"] * (((i-1)^m_weibull_shape["OD", "NI"]) - (i^m_weibull_shape["OD", "NI"])))
  # Episode 3
  v_TDP_BUP_NI_3[i] <- m_TDP[EP3 & BUP & NI, i] <- exp(m_frailty["3", "BUP", "NI"] * m_weibull_scale["BUP", "NI"] * (((t-1)^m_weibull_shape["BUP", "NI"]) - (t^m_weibull_shape["BUP", "NI"]))) # (survival curve at time i)/(survival curve at time i-1) 
  v_TDP_MET_NI_3[i] <- m_TDP[EP3 & MET & NI, i] <- exp(m_frailty["3", "MET", "NI"] * m_weibull_scale["MET", "NI"] * (((t-1)^m_weibull_shape["MET", "NI"]) - (t^m_weibull_shape["MET", "NI"])))
  v_TDP_ABS_NI_3[i] <- m_TDP[EP3 & ABS & NI, i] <- exp(m_frailty["3", "ABS", "NI"] * m_weibull_scale["ABS", "NI"] * (((i-1)^m_weibull_shape["ABS", "NI"]) - (i^m_weibull_shape["ABS", "NI"])))
  v_TDP_REL_NI_3[i] <- m_TDP[EP3 & REL & NI, i] <- exp(m_frailty["3", "REL", "NI"] * m_weibull_scale["REL", "NI"] * (((t-1)^m_weibull_shape["REL", "NI"]) - (t^m_weibull_shape["REL", "NI"])))
  v_TDP_OD_NI_3[i]  <- m_TDP[EP3 & OD & NI, i]  <- exp(m_frailty["3", "OD" , "NI"] * m_weibull_scale["OD", "NI"] * (((i-1)^m_weibull_shape["OD", "NI"]) - (i^m_weibull_shape["OD", "NI"])))
  
  # Injection
  # Episode 1
  v_TDP_BUP_INJ_1[i] <- m_TDP[EP1 & BUP & INJ, i] <- exp(m_frailty["1", "BUP", "INJ"] * m_weibull_scale["BUP", "INJ"] * (((t-1)^m_weibull_shape["BUP", "INJ"]) - (t^m_weibull_shape["BUP", "INJ"]))) # (survival curve at time i)/(survival curve at time i-1) 
  v_TDP_MET_INJ_1[i] <- m_TDP[EP1 & MET & INJ, i] <- exp(m_frailty["1", "MET", "INJ"] * m_weibull_scale["MET", "INJ"] * (((t-1)^m_weibull_shape["MET", "INJ"]) - (t^m_weibull_shape["MET", "INJ"])))
  v_TDP_ABS_INJ_1[i] <- m_TDP[EP1 & ABS & INJ, i] <- exp(m_frailty["1", "ABS", "INJ"] * m_weibull_scale["ABS", "INJ"] * (((i-1)^m_weibull_shape["ABS", "INJ"]) - (i^m_weibull_shape["ABS", "INJ"])))
  v_TDP_REL_INJ_1[i] <- m_TDP[EP1 & REL & INJ, i] <- exp(m_frailty["1", "REL", "INJ"] * m_weibull_scale["REL", "INJ"] * (((t-1)^m_weibull_shape["REL", "INJ"]) - (t^m_weibull_shape["REL", "INJ"])))
  v_TDP_OD_INJ_1[i]  <- m_TDP[EP1 & OD & INJ, i]  <- exp(m_frailty["1", "OD" , "INJ"] * m_weibull_scale["OD", "INJ"] * (((i-1)^m_weibull_shape["OD", "INJ"]) - (i^m_weibull_shape["OD", "INJ"])))
  # Episode 2
  v_TDP_BUP_INJ_2[i] <- m_TDP[EP2 & BUP & INJ, i] <- exp(m_frailty["2", "BUP", "INJ"] * m_weibull_scale["BUP", "INJ"] * (((t-1)^m_weibull_shape["BUP", "INJ"]) - (t^m_weibull_shape["BUP", "INJ"]))) # (survival curve at time i)/(survival curve at time i-1) 
  v_TDP_MET_INJ_2[i] <- m_TDP[EP2 & MET & INJ, i] <- exp(m_frailty["2", "MET", "INJ"] * m_weibull_scale["MET", "INJ"] * (((t-1)^m_weibull_shape["MET", "INJ"]) - (t^m_weibull_shape["MET", "INJ"])))
  v_TDP_ABS_INJ_2[i] <- m_TDP[EP2 & ABS & INJ, i] <- exp(m_frailty["2", "ABS", "INJ"] * m_weibull_scale["ABS", "INJ"] * (((i-1)^m_weibull_shape["ABS", "INJ"]) - (i^m_weibull_shape["ABS", "INJ"])))
  v_TDP_REL_INJ_2[i] <- m_TDP[EP2 & REL & INJ, i] <- exp(m_frailty["2", "REL", "INJ"] * m_weibull_scale["REL", "INJ"] * (((t-1)^m_weibull_shape["REL", "INJ"]) - (t^m_weibull_shape["REL", "INJ"])))
  v_TDP_OD_INJ_2[i]  <- m_TDP[EP2 & OD & INJ, i]  <- exp(m_frailty["2", "OD" , "INJ"] * m_weibull_scale["OD", "INJ"] * (((i-1)^m_weibull_shape["OD", "INJ"]) - (i^m_weibull_shape["OD", "INJ"])))
  # Episode 3
  v_TDP_BUP_INJ_3[i] <- m_TDP[EP3 & BUP & INJ, i] <- exp(m_frailty["3", "BUP", "INJ"] * m_weibull_scale["BUP", "INJ"] * (((t-1)^m_weibull_shape["BUP", "INJ"]) - (t^m_weibull_shape["BUP", "INJ"]))) # (survival curve at time i)/(survival curve at time i-1) 
  v_TDP_MET_INJ_3[i] <- m_TDP[EP3 & MET & INJ, i] <- exp(m_frailty["3", "MET", "INJ"] * m_weibull_scale["MET", "INJ"] * (((t-1)^m_weibull_shape["MET", "INJ"]) - (t^m_weibull_shape["MET", "INJ"])))
  v_TDP_ABS_INJ_3[i] <- m_TDP[EP3 & ABS & INJ, i] <- exp(m_frailty["3", "ABS", "INJ"] * m_weibull_scale["ABS", "INJ"] * (((i-1)^m_weibull_shape["ABS", "INJ"]) - (i^m_weibull_shape["ABS", "INJ"])))
  v_TDP_REL_INJ_3[i] <- m_TDP[EP3 & REL & INJ, i] <- exp(m_frailty["3", "REL", "INJ"] * m_weibull_scale["REL", "INJ"] * (((t-1)^m_weibull_shape["REL", "INJ"]) - (t^m_weibull_shape["REL", "INJ"])))
  v_TDP_OD_INJ_3[i]  <- m_TDP[EP3 & OD & INJ, i]  <- exp(m_frailty["3", "OD" , "INJ"] * m_weibull_scale["OD", "INJ"] * (((i-1)^m_weibull_shape["OD", "INJ"]) - (i^m_weibull_shape["OD", "INJ"])))
}

# Create dummy data
m_frailty <- array(0.4, dim = c(n_EP, n_BASE, n_INJECT),
                   dimnames = list(EP, BASE, INJECT))
m_weibull_scale <- array(0.7, dim = c(n_BASE, n_INJECT),
                         dimnames = list(BASE, INJECT))
m_weibull_shape <- array(0.7, dim = c(n_BASE, n_INJECT),
                         dimnames = list(BASE, INJECT))



###############################
### Age-dependent mortality ###
###############################
# Baseline mortality by age
lt_can_2018 <- read.csv("data/01_all_cause_mortality.csv")
#v_mortality_age <- lt_can_2018[which(lt_can_2018$Age >= n_age_init & lt_can_2018$Age < n_age_max), ] 

v_r_mort_by_age <- lt_can_2018 %>%
  #subset(Age >= n_age_init & Age < n_age_max) %>%
  #filter(Age >= n_age_init & Age < n_age_max) %>%
  select(Total) %>%
  as.matrix()

# Hazard ratios for death probability
hr_s <- read.csv("data/01_death_hr.csv", header = TRUE)
#hr_s
#hr_BUP1 <- hr_s["BUP1", "HR"] figure out why this doesn't work
hr_BUP1_NI <- hr_s[1, 2]
hr_BUP_NI  <- hr_s[2, 2]
hr_MET1_NI <- hr_s[3, 2]
hr_MET_NI  <- hr_s[4, 2]
hr_REL1_NI <- hr_s[5, 2]
hr_REL_NI  <- hr_s[6, 2]
hr_OD_NI   <- hr_s[7, 2]
hr_ABS_NI  <- hr_s[8, 2]

hr_BUP1_INJ <- hr_s[9, 2]
hr_BUP_INJ  <- hr_s[10, 2]
hr_MET1_INJ <- hr_s[11, 2]
hr_MET_INJ  <- hr_s[12, 2]
hr_REL1_INJ <- hr_s[13, 2]
hr_REL_INJ  <- hr_s[14, 2]
hr_OD_INJ   <- hr_s[15, 2]
hr_ABS_INJ  <- hr_s[16, 2]

hr_ABS_HIV  <- hr_s[17, 2]

########################## 
### HIV Seroconversion ###
##########################
p_sero <- read.csv("data/01_hiv_sero.csv", header = TRUE)
p_sero

# Only applied to injection
p_sero_BUP1_NI <- p_sero[1, 2]
p_sero_BUP_NI  <- p_sero[2, 2]
p_sero_MET1_NI <- p_sero[3, 2]
p_sero_MET_NI  <- p_sero[4, 2]
p_sero_REL1_NI <- p_sero[5, 2]
p_sero_REL_NI  <- p_sero[6, 2]
p_sero_OD_NI   <- p_sero[7, 2]
p_sero_ABS_NI  <- p_sero[8, 2]

########################
### Model Parameters ###
########################

###############################
### Initial characteristics ###
###############################

#v_initial_BUP <- read.csv("data/01_initial_BUP.csv")
#v_initial_MET <- read.csv("data/01_initial_MET.csv")

n_age_init <- 35 # age at baseline
n_age_max <- 95 # maximum age of follow up

n_t <- (n_age_max - n_age_init) * 12 # modeling time horizon in months
p_discount <- 0.03

p_male_BUP <- 0.35
p_male_MET <- 0.40

p_HIV_POS <- 0.05 # % of HIV-positive individuals
p_HIV_ART <- 0.75 # % of HIV-positive on-ART

#######################################################  
### Transitions conditional on leaving and survival ###
#######################################################
# a_P[i, j, k] = Transition probibility array [from, to, time(month)]
# Transition array already populated with zeros, only need to define possible transitions
# Only mortality probability and remain probability changing over time
# Other transitions re-weighted by probability of leaving

# Transition probabilities among all base states (divide by n(strata) if not strata-specific)
#n_strata <-  # number of non-base state strata




### Function to check transition array
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

a <- as.vector()
b <- as.vector()

for (i in 1:n_t){
  a[i] <- sum(a_TDP[BUP1 & NI & EP1,,i])
  b[i] <- sum(a_TDP[BUP1 & INJ & EP1,,i])
}

check_transition_probability(a_TDP, err_stop = FALSE, verbose = TRUE)
check_sum_of_transition_array(a_TDP, n_states, n_t, err_stop = FALSE, verbose = TRUE)

# Set first row of m.M with the initial state vector
# Baseline
m_M_BL[1, ]  <- v_s_init_BL
# BUP
m_M_BUP[1, ] <- v_s_init_BUP
# MET
m_M_MET[1, ] <- v_s_init_MET


# Iterate over all time periods
for(i in 1:n_t){
  # BUP
  m_M_BUP[i + 1, ] <- m_M_BUP[i, ] %*% a_P_BUP[, , i]
  # MET
  m_M_MET[i + 1, ] <- m_M_MET[i, ] %*% a_P_MET[, , i]
  # BL
  m_M_BL[i + 1, ]  <- m_M_BL[i, ] %*% a_P_BL[, , i]
  
  # Create empty initial state matrices
  # BUP/MET
  #m_M_BL <- m_M_BUP <- m_M_MET <- matrix(0, 
  #                                       nrow = (n_t + 1), ncol = n_states, 
  #                                       dimnames = list(0:n_t, v_n))
  
  
  # Iterate over all time periods
  for(t in 1:n_t){
    # BUP
    m_M_BUP[t + 1, ] <- m_M_BUP[t, ] %*% a_P_BUP[, , t]
    # MET
    m_M_MET[t + 1, ] <- m_M_MET[t, ] %*% a_P_MET[, , t]
  }
  m_M_BUP
  m_M_MET
  
  
  
  ############################################
  #### Check if transition array is valid ####
  ############################################
  check_sum_of_transition_array <- function(a_P,
                                            n_states_to,
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
  
  check_transition_probability(a_P, err_stop = err_stop, verbose = verbose)
  check_sum_of_transition_array(a_P, n_states_to, n_t, err_stop = err_stop, verbose = verbose)
  
  
  cbind(TOT = colSums(m_M_trace[i, ]), # total across all states at each time point
        
        BUP = colSums(m_M_trace[i, BUP]), # totals in base states (no strat)
        MET = colSums(m_M_trace[i, MET]),
        REL = colSums(m_M_trace[i, REL]),
        OD  = colSums(m_M_trace[i, OD]),
        ABS = colSums(m_M_trace[i, ABS]),
        
        HIV = colSums(m_M_trace[, POS]), # total with HIV
        
        INJ = colSums(m_M_trace[, INJ]))
  
  v_deaths <- array(0, dim = c(n_t, 1))
  for (i in 1:n_t){
    deaths[i,] <- as.vector(1 - rowSums(m_M_trace[i,]))
  }
  deaths
  
  
  
  v_agg_trace_episodes <- c("EP1", "EP2", "EP3")
  n_agg_trace_episodes <- length(v_agg_trace_episodes)
  m_M_agg_trace_ep <- array(0, dim = c(n_t + 1, n_agg_trace_episodes),
                            dimnames = list(0:n_t, v_agg_trace_episodes))
  
  for (i in 1:n_t){
    m_M_agg_trace_ep[i, "EP1"] <- sum(m_M_trace[i, EP1])  
    m_M_agg_trace_ep[i, "EP2"] <- sum(m_M_trace[i, EP2])
    m_M_agg_trace_ep[i, "EP3"] <- sum(m_M_trace[i, EP3])
  }
  
  
  
  # Non-injection
  a_TDP[BUP1 & NI & NEG, BUP1 & NI & NEG, ]  <- a_TDP[BUP1 & NI & NEG, BUP1 & NI & NEG, ] * (1 - p_sero_BUP1_NI)
  a_TDP[BUP1 & NI & NEG, BUP1 & NI & POS, ]  <- a_TDP[BUP1 & NI & NEG, BUP1 & NI & POS, ] * p_sero_BUP1_NI
  a_TDP[BUP & NI & NEG , BUP & NI & NEG, ]   <- a_TDP[BUP & NI & NEG , BUP & NI & NEG, ] * (1 - p_sero_BUP_NI)
  a_TDP[BUP & NI & NEG , BUP & NI & POS, ]   <- a_TDP[BUP & NI & NEG , BUP & NI & POS, ] * p_sero_BUP_NI
  a_TDP[MET1 & NI & NEG, MET1 & NI & NEG, ]  <- a_TDP[MET1 & NI & NEG, MET1 & NI & NEG, ] * (1 - p_sero_MET1_NI)
  a_TDP[MET1 & NI & NEG, MET1 & NI & POS, ]  <- a_TDP[MET1 & NI & NEG, MET1 & NI & POS, ] * p_sero_MET1_NI
  a_TDP[MET & NI & NEG , MET & NI & NEG, ]   <- a_TDP[MET & NI & NEG , MET & NI & NEG, ] * (1 - p_sero_MET_NI)
  a_TDP[MET & NI & NEG , MET & NI & POS, ]   <- a_TDP[MET & NI & NEG , MET & NI & POS, ] * p_sero_MET_NI
  a_TDP[REL1 & NI & NEG, REL1 & NI & NEG, ]  <- a_TDP[REL1 & NI & NEG, REL1 & NI & NEG, ] * (1 - p_sero_REL1_NI)
  a_TDP[REL1 & NI & NEG, REL1 & NI & POS, ]  <- a_TDP[REL1 & NI & NEG, REL1 & NI & POS, ] * p_sero_REL1_NI
  a_TDP[REL & NI & NEG , REL & NI & NEG, ]   <- a_TDP[REL & NI & NEG , REL & NI & NEG, ] * (1 - p_sero_REL_NI)
  a_TDP[REL & NI & NEG , REL & NI & POS, ]   <- a_TDP[REL & NI & NEG , REL & NI & POS, ] * p_sero_REL_NI
  a_TDP[OD & NI & NEG  , OD & NI & NEG, ]    <- a_TDP[OD & NI & NEG  , OD & NI & NEG, ] * (1 - p_sero_OD_NI)
  a_TDP[OD & NI & NEG  , OD & NI & POS, ]    <- a_TDP[OD & NI & NEG  , OD & NI & POS, ] * p_sero_OD_NI
  a_TDP[ABS & NI & NEG , ABS & NI & NEG, ]   <- a_TDP[ABS & NI & NEG , ABS & NI & NEG, ] * (1 - p_sero_ABS_NI)
  a_TDP[ABS & NI & NEG , ABS & NI & POS, ]   <- a_TDP[ABS & NI & NEG , ABS & NI & POS, ] * p_sero_ABS_NI
  
  # Injection
  a_TDP[BUP1 & INJ & NEG, BUP1 & INJ & NEG, ]  <- a_TDP[BUP1 & INJ & NEG, BUP1 & INJ & NEG, ] * (1 - p_sero_BUP1_INJ)
  a_TDP[BUP1 & INJ & NEG, BUP1 & INJ & POS, ]  <- a_TDP[BUP1 & INJ & NEG, BUP1 & INJ & POS, ] * p_sero_BUP1_INJ
  a_TDP[BUP & INJ & NEG , BUP & INJ & NEG, ]   <- a_TDP[BUP & INJ & NEG , BUP & INJ & NEG, ] * (1 - p_sero_BUP_INJ)
  a_TDP[BUP & INJ & NEG , BUP & INJ & POS, ]   <- a_TDP[BUP & INJ & NEG , BUP & INJ & POS, ] * p_sero_BUP_INJ
  a_TDP[MET1 & INJ & NEG, MET1 & INJ & NEG, ]  <- a_TDP[MET1 & INJ & NEG, MET1 & INJ & NEG, ] * (1 - p_sero_MET1_INJ)
  a_TDP[MET1 & INJ & NEG, MET1 & INJ & POS, ]  <- a_TDP[MET1 & INJ & NEG, MET1 & INJ & POS, ] * p_sero_MET1_INJ
  a_TDP[MET & INJ & NEG , MET & INJ & NEG, ]   <- a_TDP[MET & INJ & NEG , MET & INJ & NEG, ] * (1 - p_sero_MET_INJ)
  a_TDP[MET & INJ & NEG , MET & INJ & POS, ]   <- a_TDP[MET & INJ & NEG , MET & INJ & POS, ] * p_sero_MET_INJ
  a_TDP[REL1 & INJ & NEG, REL1 & INJ & NEG, ]  <- a_TDP[REL1 & INJ & NEG, REL1 & INJ & NEG, ] * (1 - p_sero_REL1_INJ)
  a_TDP[REL1 & INJ & NEG, REL1 & INJ & POS, ]  <- a_TDP[REL1 & INJ & NEG, REL1 & INJ & POS, ] * p_sero_REL1_INJ
  a_TDP[REL & INJ & NEG , REL & INJ & NEG, ]   <- a_TDP[REL & INJ & NEG , REL & INJ & NEG, ] * (1 - p_sero_REL_INJ)
  a_TDP[REL & INJ & NEG , REL & INJ & POS, ]   <- a_TDP[REL & INJ & NEG , REL & INJ & POS, ] * p_sero_REL_INJ
  a_TDP[OD & INJ & NEG  , OD & INJ & NEG, ]    <- a_TDP[OD & INJ & NEG  , OD & INJ & NEG, ] * (1 - p_sero_OD_INJ)
  a_TDP[OD & INJ & NEG  , OD & INJ & POS, ]    <- a_TDP[OD & INJ & NEG  , OD & INJ & POS, ] * p_sero_OD_INJ
  a_TDP[ABS & INJ & NEG , ABS & INJ & NEG, ]   <- a_TDP[ABS & INJ & NEG , ABS & INJ & NEG, ] * (1 - p_sero_ABS_INJ)
  a_TDP[ABS & INJ & NEG , ABS & INJ & POS, ]   <- a_TDP[ABS & INJ & NEG , ABS & INJ & POS, ] * p_sero_ABS_INJ
  
  
  # Require library 'optim'
  # Create wrapper function around model
  fr <- function(x, data){   
    x1 <- x[1] # Free parameters to calibrate
    x2 <- x[2]
    outcome <- run_model(x...) # code to run model and return outputs required by calibration routine
    y <- (w_o * (outcome - obs_outcome)^2) # weighted GOF function
    return(y)
  }
  
  load_mort_data <- function(file = NULL){
    # Load mortality data from file
    if(!is.null(file)) {
      df_r_mort_by_age <- read.csv(file = file)}
    else{
      df_r_mort_by_age <- all_cause_mortality
    }
    # Vector with mortality rates
    v_r_mort_by_age  <- as.matrix(dplyr::select(df_r_mort_by_age, .data$Total))
    
    return(v_r_mort_by_age)
  }
  
  # Episode 2
  v_TDP_BUP_NI_2 <- vector()
  v_TDP_MET_NI_2 <- vector()
  v_TDP_ABS_NI_2 <- vector()
  v_TDP_REL_NI_2 <- vector()
  v_TDP_OD_NI_2  <- vector()
  # Episode 3
  v_TDP_BUP_NI_3 <- vector()
  v_TDP_MET_NI_3 <- vector()
  v_TDP_ABS_NI_3 <- vector()
  v_TDP_REL_NI_3 <- vector()
  v_TDP_OD_NI_3  <- vector()
  
  # Injection
  # Episode 1
  v_TDP_BUP_INJ_1 <- vector()
  v_TDP_MET_INJ_1 <- vector()
  v_TDP_ABS_INJ_1 <- vector()
  v_TDP_REL_INJ_1 <- vector()
  v_TDP_OD_INJ_1  <- vector()
  # Episode 2
  v_TDP_BUP_INJ_2 <- vector()
  v_TDP_MET_INJ_2 <- vector()
  v_TDP_ABS_INJ_2 <- vector()
  v_TDP_REL_INJ_2 <- vector()
  v_TDP_OD_INJ_2  <- vector()
  # Episode 3
  v_TDP_BUP_INJ_3 <- vector()
  v_TDP_MET_INJ_3 <- vector()
  v_TDP_ABS_INJ_3 <- vector()
  v_TDP_REL_INJ_3 <- vector()
  v_TDP_OD_INJ_3  <- vector()
  
  load_all_params <- function(file.init = NULL,
                              file.mort = NULL,
                              file.weibull = NULL,
                              file.unconditional = NULL,
                              file.sero = NULL
  ){ # User defined
    #### Load initial set of initial parameters from .csv file ####
    if(!is.null(file.init)) {
      df_params_init <- read.csv(file = file.init)
    } else{
      df_params_init <- df_params_init
    }
    
    #### All-cause age-specific mortality from .csv file ####
    v_r_mort_by_age <- load_mort_data(file = file.mort)
    
    l_params_all <- with(as.list(df_params_init), {
      #### General setup ####
      v_names_str <- c("No Treatment", "Treatment")  # CEA strategies
      n_str       <- length(v_names_str) # Number of strategies
      v_age_names <- n_age_init:(n_age_init + n_t - 1) # vector with age names
      v_n <- c("H", "S1", "S2", "D")  # vector with the 4 health states of the model:
      # Healthy (H), Sick (S1), Sicker (S2), Dead (D)
      n_states <- length(v_n)         # number of health states 
      v_s_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
      #### Create list with all parameters ####
      
      
      
      l_params_all <- list(
        
        v_names_str = v_names_str,
        n_str       = n_str      ,
        n_age_init  = n_age_init, 
        n_t         = n_t       , 
        v_age_names = v_age_names,
        v_n = v_n,
        n_states = n_states,
        v_s_init = c(H = 1, S1 = 0, S2 = 0, D = 0),
        v_r_mort_by_age = v_r_mort_by_age
        
      )
      
      
      
      
      
      
      
      return(l_params_all)
    }
    )
    
    l_params_all <- c(l_params_all, 
                      df_params_init) # Add initial set of parameters
  }
  
  
  
  #### START OF WORKING CODE ####
  n_age_init <- 35 # age at baseline
  n_age_max <- 95 # maximum age of follow up
  
  n_t <- (n_age_max - n_age_init) * 12 # modeling time horizon in months
  p_discount <- 0.03
  
  p_male_BUP <- 0.35
  p_male_MET <- 0.40
  
  p_HIV_POS <- 0.05 # % of HIV-positive individuals
  p_HIV_ART <- 0.75 # % of HIV-positive on-ART
  
  
  ### Time-dependent survival probabilities ###
  # Initialize empty vectors
  v_TDP_BUP_NI_1 = v_TDP_MET_NI_1 = v_TDP_ABS_NI_1 = v_TDP_REL_NI_1 = v_TDP_OD_NI_1 = vector(),
  v_TDP_BUP_NI_2 = v_TDP_MET_NI_2 = v_TDP_ABS_NI_2 = v_TDP_REL_NI_2 = v_TDP_OD_NI_2 = vector(),
  v_TDP_BUP_NI_3 = v_TDP_MET_NI_3 = v_TDP_ABS_NI_3 = v_TDP_REL_NI_3 = v_TDP_OD_NI_3 = vector(),
  v_TDP_BUP_INJ_1 = v_TDP_MET_INJ_1 = v_TDP_ABS_INJ_1 = v_TDP_REL_INJ_1 = v_TDP_OD_INJ_1 = vector(),
  v_TDP_BUP_INJ_2 = v_TDP_MET_INJ_2 = v_TDP_ABS_INJ_2 = v_TDP_REL_INJ_2 = v_TDP_OD_INJ_2 = vector(),
  v_TDP_BUP_INJ_3 = v_TDP_MET_INJ_3 = v_TDP_ABS_INJ_3 = v_TDP_REL_INJ_3 = v_TDP_OD_INJ_3 = vector(),
  # Probability of remaining in given health state
  for(i in 1:n_t){
    t <- i+1 # Shift BUP, MET, REL ahead 1 month to adjust for BUP1, MET1, REL1
    # Non-injection
    # Episode 1
    v_TDP_BUP_NI_1[i] = as.vector(exp(p_frailty_BUP_NI_1 * p_weibull_scale_BUP_NI * (((t - 1)^p_weibull_shape_BUP_NI) - (t^p_weibull_shape_BUP_NI)))) # (survival curve at time i)/(survival curve at time i-1) 
    v_TDP_MET_NI_1[i] = as.vector(exp(p_frailty_MET_NI_1 * p_weibull_scale_MET_NI * (((t - 1)^p_weibull_shape_MET_NI) - (t^p_weibull_shape_MET_NI))))
    v_TDP_ABS_NI_1[i] = as.vector(exp(p_frailty_ABS_NI_1 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
    v_TDP_REL_NI_1[i] = as.vector(exp(p_frailty_REL_NI_1 * p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
    v_TDP_OD_NI_1[i]  = as.vector(exp(p_frailty_OD_NI_1  * p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))
    # Episode 2
    v_TDP_BUP_NI_2[i] = as.vector(exp(p_frailty_BUP_NI_2 * p_weibull_scale_BUP_NI * (((t - 1)^p_weibull_shape_BUP_NI) - (t^p_weibull_shape_BUP_NI)))) # (survival curve at time i)/(survival curve at time i-1) 
    v_TDP_MET_NI_2[i] = as.vector(exp(p_frailty_MET_NI_2 * p_weibull_scale_MET_NI * (((t - 1)^p_weibull_shape_MET_NI) - (t^p_weibull_shape_MET_NI))))
    v_TDP_ABS_NI_2[i] = as.vector(exp(p_frailty_ABS_NI_2 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
    v_TDP_REL_NI_2[i] = as.vector(exp(p_frailty_REL_NI_2 * p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
    v_TDP_OD_NI_2[i]  = as.vector(exp(p_frailty_OD_NI_2  * p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))
    # Episode 3
    v_TDP_BUP_NI_3[i] = as.vector(exp(p_frailty_BUP_NI_3 * p_weibull_scale_BUP_NI * (((t - 1)^p_weibull_shape_BUP_NI) - (t^p_weibull_shape_BUP_NI)))) # (survival curve at time i)/(survival curve at time i-1) 
    v_TDP_MET_NI_3[i] = as.vector(exp(p_frailty_MET_NI_3 * p_weibull_scale_MET_NI * (((t - 1)^p_weibull_shape_MET_NI) - (t^p_weibull_shape_MET_NI))))
    v_TDP_ABS_NI_3[i] = as.vector(exp(p_frailty_ABS_NI_3 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
    v_TDP_REL_NI_3[i] = as.vector(exp(p_frailty_REL_NI_3 * p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
    v_TDP_OD_NI_3[i]  = as.vector(exp(p_frailty_OD_NI_3  * p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))
    
    # Injection
    # Episode 1
    v_TDP_BUP_INJ_1[i] = as.vector(exp(p_frailty_BUP_INJ_1 * p_weibull_scale_BUP_INJ * (((t - 1)^p_weibull_shape_BUP_INJ) - (t^p_weibull_shape_BUP_INJ)))) # (survival curve at time i)/(survival curve at time i-1) 
    v_TDP_MET_INJ_1[i] = as.vector(exp(p_frailty_MET_INJ_1 * p_weibull_scale_MET_INJ * (((t - 1)^p_weibull_shape_MET_INJ) - (t^p_weibull_shape_MET_INJ))))
    v_TDP_ABS_INJ_1[i] = as.vector(exp(p_frailty_ABS_INJ_1 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
    v_TDP_REL_INJ_1[i] = as.vector(exp(p_frailty_REL_INJ_1 * p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
    v_TDP_OD_INJ_1[i]  = as.vector(exp(p_frailty_OD_INJ_1  * p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ))))
    # Episode 2
    v_TDP_BUP_INJ_2[i] = as.vector(exp(p_frailty_BUP_INJ_2 * p_weibull_scale_BUP_INJ * (((t - 1)^p_weibull_shape_BUP_INJ) - (t^p_weibull_shape_BUP_INJ)))) # (survival curve at time i)/(survival curve at time i-1) 
    v_TDP_MET_INJ_2[i] = as.vector(exp(p_frailty_MET_INJ_2 * p_weibull_scale_MET_INJ * (((t - 1)^p_weibull_shape_MET_INJ) - (t^p_weibull_shape_MET_INJ))))
    v_TDP_ABS_INJ_2[i] = as.vector(exp(p_frailty_ABS_INJ_2 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
    v_TDP_REL_INJ_2[i] = as.vector(exp(p_frailty_REL_INJ_2 * p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
    v_TDP_OD_INJ_2[i]  = as.vector(exp(p_frailty_OD_INJ_2  * p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ))))
    # Episode 3
    v_TDP_BUP_INJ_3[i] = as.vector(exp(p_frailty_BUP_INJ_3 * p_weibull_scale_BUP_INJ * (((t - 1)^p_weibull_shape_BUP_INJ) - (t^p_weibull_shape_BUP_INJ)))) # (survival curve at time i)/(survival curve at time i-1) 
    v_TDP_MET_INJ_3[i] = as.vector(exp(p_frailty_MET_INJ_3 * p_weibull_scale_MET_INJ * (((t - 1)^p_weibull_shape_MET_INJ) - (t^p_weibull_shape_MET_INJ))))
    v_TDP_ABS_INJ_3[i] = as.vector(exp(p_frailty_ABS_INJ_3 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
    v_TDP_REL_INJ_3[i] = as.vector(exp(p_frailty_REL_INJ_3 * p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
    v_TDP_OD_INJ_3[i]  = as.vector(exp(p_frailty_OD_INJ_3  * p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ))))
  }
  
  # Aggregated trace matrix
  #id=factor(c("Second", "First", "Third"), levels=c("First", "Second", "Third"))
  v_agg_trace_states <- c("Alive", "Death", "OD", "REL1", "REL", "BUP1", "BUP", "MET1", "MET", "ABS")
  n_agg_trace_states <- length(v_agg_trace_states)
  m_M_agg_trace <- array(0, dim = c((n_t + 1), n_agg_trace_states),
                         dimnames = list(0:n_t, v_agg_trace_states))
  
  for (i in 1:n_t){
    m_M_agg_trace[i, "Alive"] <- sum(m_M_trace[i, ])
    m_M_agg_trace[i, "BUP1"]   <- sum(m_M_trace[i, BUP1])
    m_M_agg_trace[i, "BUP"]   <- sum(m_M_trace[i, BUP])
    m_M_agg_trace[i, "MET1"]   <- sum(m_M_trace[i, MET1])
    m_M_agg_trace[i, "MET"]   <- sum(m_M_trace[i, MET])
    m_M_agg_trace[i, "REL1"]   <- sum(m_M_trace[i, REL1])
    m_M_agg_trace[i, "REL"]   <- sum(m_M_trace[i, REL])
    m_M_agg_trace[i, "ABS"]   <- sum(m_M_trace[i, ABS])
    m_M_agg_trace[i, "OD"]    <- sum(m_M_trace[i, OD])
    #m_M_agg_trace[i, "HIV"]  <- sum(m_M_trace[i, POS]) # Need cumulative sum for HIV
    m_M_agg_trace[i, "Death"] <- 1 - sum(m_M_trace[i, ])
  }
  
  #### Create plots ####
  # Run model
  l_out_markov <- markov_model(l_params_all = l_params_all)
  
  # Prepare data
  df_M_agg_trace <- as.data.frame(l_out_markov$m_M_agg_trace)
  df_M_agg_trace$month <- as.numeric(rownames(df_M_agg_trace))
  df_M_agg_trace_plot <- df_M_agg_trace %>% gather(state, proportion, "Death", "OD", "REL1", "REL", "BUP1", "BUP", "MET1", "MET", "ABS") # health states to plot
  df_M_agg_state_time <- df_M_agg_trace %>% gather(state, proportion, "OD", "REL1", "REL", "BUP1", "BUP", "MET1", "MET", "ABS") # alive health states to plot
  df_M_agg_state_time <- df_M_agg_state_time %>% 
    group_by(state) %>% 
    summarise_each(funs(sum), proportion) %>%
    mutate(percentage = round((proportion / sum(proportion)) * 100,1))
  
  # Preserve order for plotting
  state_order_trace <- factor(df_M_agg_trace_plot$state, levels = c("Death", "OD", "REL1", "REL", "BUP1", "BUP", "MET1", "MET", "ABS"))
  state_order_time  <- factor(df_M_agg_state_time$state, levels = c("OD", "REL1", "REL", "BUP1", "BUP", "MET1", "MET", "ABS"))
  state_colours_trace <- c("#d9d9d9", "#d53e4f", "#f46d43", "#fdae61", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd")
  state_colours_time <- c("#d53e4f", "#f46d43", "#fdae61", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd")
  
  ### Markov trace plots ###
  main_states_trace_plot <- ggplot(df_M_agg_trace_plot, aes(x = month, y = proportion, fill = state_order_trace)) + 
    theme_bw() +
    geom_area() +
    scale_fill_manual(values = state_colours_trace)
  pdf("Plots/Markov Trace/trace_states.pdf")
  main_states_trace_plot
  dev.off()
  
  ### Time spent in health states ###
  main_states_time <- ggplot(df_M_agg_state_time, aes(x = state_order_time, y = proportion, fill = state_order_time)) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab("Health State") + ylab("Time") +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = state_colours_time) +
    geom_text(aes(label = percentage), hjust = -0.25, size = 3.5) +
    coord_flip(ylim = c(0, 720))
  pdf(file = "Plots/Markov Trace/time_states.pdf", width = 8, height = 4)
  main_states_time
  dev.off()
  
  #######################################
  #### State and Transition Outcomes ####
  #######################################
  ### Cost-effectiveness analysis ###
  ### Vector of total costs for both strategies
  v_ted_cost <- c(n_totcost_UC, n_totcost_Tr)
  ### Vector of effectiveness for both strategies
  v_ted_qaly <- c(n_totqaly_UC, n_totqaly_Tr)
  
  ### Calculate incremental cost-effectiveness ratios (ICERs)
  df_cea <- calculate_icers(cost = v_ted_cost,
                            effect = v_ted_qaly,
                            strategies = v_names_str)
  df_cea
  ### Create CEA table
  table_cea <- df_cea
  ## Format column names
  colnames(table_cea)[2:6] <- c("Costs ($)", "QALYs",
                                "Incremental Costs ($)", "Incremental QALYs",
                                "ICER ($/QALY)") # name the columns
  ## Format rows
  table_cea$`Costs ($)` <- comma(round(table_cea$`Costs ($)`, 0))
  table_cea$`Incremental Costs ($)` <- comma(round(table_cea$`Incremental Costs ($)`, 0))
  table_cea$QALYs <- round(table_cea$QALYs, 3)
  table_cea$`Incremental QALYs` <- round(table_cea$`Incremental QALYs`, 3)
  table_cea$`ICER ($/QALY)` <- comma(round(table_cea$`ICER ($/QALY)`, 0))
  table_cea
  ### CEA frontier
  plot(df_cea) +
    expand_limits(x = 20.8)
  
  
  #k_hiv <- a_TDP[, , 50]
  #k_hiv2<- a_TDP[, , 710]
  
  #write.csv(k_hiv,"C:/Users/Benjamin/Desktop/K2.csv", row.names = TRUE)
  #write.csv(k_hiv2,"C:/Users/Benjamin/Desktop/K3.csv", row.names = TRUE)
  
  
  c_BUP_NI_crime = df_crime_costs[, "BUP_NI"],
  c_MET_NI_crime = df_crime_costs[, "MET_NI"],
  c_REL_NI_crime = df_crime_costs[, "REL_NI"],
  c_OD_NI_crime = df_crime_costs[, "OD_NI"],
  c_ABS_NI_crime = df_crime_costs[, "ABS_NI"],
  
  c_BUP_INJ_crime = df_crime_costs[, "BUP_INJ"],
  c_MET_INJ_crime = df_crime_costs[, "MET_INJ"],
  c_REL_INJ_crime = df_crime_costs[, "REL_INJ"],
  c_OD_INJ_crime = df_crime_costs[, "OD_INJ"],
  c_ABS_INJ_crime = df_crime_costs[, "ABS_INJ"],
  
  
  gh <- outcomes(l_params_all = l_params_all, l_out_markov = l_out_markov)
  ty <- gh$m_TOTAL_costs_states
  a <- gh$m_TX_costs
  b <- gh$m_HRU_costs
  c <- gh$m_HIV_costs
  d <- gh$m_crime_costs
  write.csv(ty,"C:/Users/Benjamin/Desktop/total_costs.csv", row.names = TRUE)
  write.csv(l_out_markov$m_M_trace,"C:/Users/Benjamin/Desktop/trace.csv", row.names = TRUE)
  write.csv(a,"C:/Users/Benjamin/Desktop/tx.csv", row.names = TRUE)
  write.csv(b,"C:/Users/Benjamin/Desktop/hru.csv", row.names = TRUE)
  write.csv(c,"C:/Users/Benjamin/Desktop/hiv.csv", row.names = TRUE)
  write.csv(d,"C:/Users/Benjamin/Desktop/crime.csv", row.names = TRUE)
  
  
  
  # Probability of successful naloxone use
  p_NX_rev <- (p_witness * p_NX_used * p_NX_success)
  
  # Probability of mortality from overdose accounting for baseline overdose fatality and effectiveness of naloxone
  # Subsets overdose into fatal and non-fatal, conditional on different parameters
  p_fatal_OD_NX <- p_fatal_OD * (1 - p_NX_rev)
  
  #' @param p_witness Probability of witnessed/attended overdose
  #' @param p_NX_used Probability of naloxone use, given witnessed overdose
  #' @param p_NX_success Probability of naloxone effectiveness, given use
  #' @param p_fatal_OD Probability of fatal overdose, conditional on overdose 
  #' 

  # Add NEG -> POS remain probabilities
  # To-do: See if there is a better way to do this
  # AUGUST 28, 2020 **UPDATE SEROCONVERSION SECTION TO INCLUDE HCV AND COI**
  for (i in 1:n_t){
    # Non-injection
    # BUP
    
    # FIX THIS SECTION FOR ALL PARTS!! MAKE SURE ALL SEROSTATUS ARE POPULATED
    #a_TDP[BUP & NI & EP1 & NEG, BUP & NI & EP1 & NEG, i] <- m_TDP[BUP & NI & EP1 & NEG, i]
    #a_TDP[BUP & NI & EP2 & NEG, BUP & NI & EP2 & NEG, i] <- m_TDP[BUP & NI & EP2 & NEG, i]
    #a_TDP[BUP & NI & EP3 & NEG, BUP & NI & EP3 & NEG, i] <- m_TDP[BUP & NI & EP3 & NEG, i]
    
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
  
  
  #for (j in 1:n_states){
  # a_TDP[j, j, i] <- m_TDP[j, i] # same-state to same-state
  #a_TDP[j, j + n_BASE_INJECT_EP, i] <- m_TDP[j, i] # 
  #a_TDP[j, j + (2 * n_BASE_INJECT_EP), i] <- m_TDP[j, i]
  #a_TDP[j, j + (3 * n_BASE_INJECT_EP), i] <- m_TDP[j, i]
  #} 
  
  
  
  #PLOTTING CALIBRATION PARAMETER ESTIMATES
  # Base overdose rates
  df_calib_post_plot_base  <- gather(df_calib_prior_post, parameter, draw, factor_key = TRUE) %>% filter(parameter == "n_TX_OD_post" | parameter == "n_TX_OD_prior" | parameter == "n_TXC_OD_post" | parameter == "n_TXC_OD_prior" | parameter == "n_REL_OD_post" | parameter == "n_REL_OD_prior")
  base_order <- factor(df_calib_post_plot_base$parameter, levels = c("n_TX_OD_post", "n_TX_OD_prior", "n_TXC_OD_post", "n_TXC_OD_prior", "n_REL_OD_post", "n_REL_OD_prior"))
  
  # Overdose rate multipliers
  df_calib_post_plot_mult  <- gather(df_calib_prior_post, parameter, draw, factor_key = TRUE) %>% filter(parameter == "n_TX_OD_mult_post" | parameter == "n_TX_OD_mult_prior" | parameter == "n_TXC_OD_mult_prior" | parameter == "n_TXC_OD_mult_prior" | parameter == "n_REL_OD_mult_post" | parameter == "n_REL_OD_mult_prior" | parameter == "n_INJ_OD_mult_post" | parameter == "n_INJ_OD_mult_prior")
  
  # Fatal overdose rate
  df_calib_post_plot_fatal <- gather(df_calib_prior_post, parameter, draw, factor_key = TRUE) %>% filter(parameter == "n_fatal_OD_post" | parameter == "n_fatal_OD_prior")
  
  # Base overdose rates
  cali_base_od_prior_post <- ggplot(df_calib_post_plot_base, 
                                    aes(x = draw, y = parameter)) +
    #geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
    geom_density_ridges(aes(fill = parameter)) +
    #scale_fill_viridis_c(name = "Monthly rate", option = "C") +
    scale_fill_manual(values = c("#2ca25f", "#e5f5f9", "#3182bd", "#deebf7", "#e6550d", "#fee6ce")) +
    labs(title = 'Base Overdose Rates')
  
  # Overdose rate multipliers
  cali_mult_od_prior_post <- ggplot(df_calib_post_plot_mult, 
                                    aes(x = draw, y = parameter)) +
    #geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
    geom_density_ridges(aes(fill = parameter)) +
    #scale_fill_viridis_c(name = "Monthly rate", option = "C") +
    scale_fill_manual(values = c("#2ca25f", "#e5f5f9", "#3182bd", "#deebf7", "#e6550d", "#fee6ce")) +
    labs(title = 'Base Overdose Rates')
  
  # Fatal overdose rates
  cali_fatal_od_prior_post <- ggplot(df_calib_post_plot_fatal, 
                                     aes(x = draw, y = parameter)) +
    #geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
    geom_density_ridges(aes(fill = parameter)) +
    #scale_fill_viridis_c(name = "Monthly rate", option = "C") +
    scale_fill_manual(values = c("#2ca25f", "#e5f5f9", "#3182bd", "#deebf7", "#e6550d", "#fee6ce")) +
    labs(title = 'Base Overdose Rates')
  
  # Output plots
  pdf("Plots/Calibration/cali_base_od_prior_post.pdf", width = 8, height = 6)
  cali_base_od_prior_post
  dev.off()
  
  df_calib_post  <- as.data.frame(m_calib_post) %>% mutate(ID = row_number()) %>% rename(n_TX_OD_post  = n_TX_OD,
                                                                                         n_TXC_OD_post = n_TXC_OD,
                                                                                         n_REL_OD_post = n_REL_OD,
                                                                                         n_TX_OD_mult_post = n_TX_OD_mult,
                                                                                         n_TXC_OD_mult_post = n_TXC_OD_mult,
                                                                                         n_REL_OD_mult_post = n_REL_OD_mult,
                                                                                         n_INJ_OD_mult_post = n_INJ_OD_mult,
                                                                                         n_fatal_OD_post = n_fatal_od) #make sure this label is changed in next iteration ("OD")
  df_calib_prior <- as.data.frame(m_calib_prior) %>% mutate(ID = row_number()) %>% rename(n_TX_OD_prior  = n_TX_OD,
                                                                                          n_TXC_OD_prior = n_TXC_OD,
                                                                                          n_REL_OD_prior = n_REL_OD,
                                                                                          n_TX_OD_mult_prior = n_TX_OD_mult,
                                                                                          n_TXC_OD_mult_prior = n_TXC_OD_mult,
                                                                                          n_REL_OD_mult_prior = n_REL_OD_mult,
                                                                                          n_INJ_OD_mult_prior = n_INJ_OD_mult,
                                                                                          n_fatal_OD_prior = n_fatal_OD)
  #df_calib_prior_post <- merge(df_calib_prior, df_calib_post, by.x = "ID", by.y = "ID") #merge prior and posterior estimates
  
  q0 = min(value),
  q01 = quantile(value, probs = .01), q02 = quantile(value, probs = .02), q03 = quantile(value, probs = .03), q04 = quantile(value, probs = .04), q05 = quantile(value, probs = .05),
  q06 = quantile(value, probs = .06), q07 = quantile(value, probs = .07), q08 = quantile(value, probs = .08), q09 = quantile(value, probs = .09), q10 = quantile(value, probs = .1),
  q11 = quantile(value, probs = .11), q12 = quantile(value, probs = .12), q13 = quantile(value, probs = .13), q14 = quantile(value, probs = .14), q15 = quantile(value, probs = .15),
  q15 = quantile(value, probs = .15), q16 = quantile(value, probs = .16), q17 = quantile(value, probs = .17), q18 = quantile(value, probs = .18), q19 = quantile(value, probs = .19),
  q20 = quantile(value, probs = .20), q21 = quantile(value, probs = .21), q22 = quantile(value, probs = .22), q23 = quantile(value, probs = .23), q24 = quantile(value, probs = .24),
  q25 = quantile(value, probs = .25), q26 = quantile(value, probs = .26), q27 = quantile(value, probs = .27), q28 = quantile(value, probs = .28), q29 = quantile(value, probs = .29),
  q30 = quantile(value, probs = .30), q31 = quantile(value, probs = .31), q32 = quantile(value, probs = .32), q33 = quantile(value, probs = .33), q34 = quantile(value, probs = .34),
  q35 = quantile(value, probs = .35), q36 = quantile(value, probs = .36), q37 = quantile(value, probs = .37), q38 = quantile(value, probs = .38), q39 = quantile(value, probs = .39),
  q40 = quantile(value, probs = .40), q41 = quantile(value, probs = .41), q42 = quantile(value, probs = .42), q43 = quantile(value, probs = .43), q44 = quantile(value, probs = .44),
  q45 = quantile(value, probs = .45), q46 = quantile(value, probs = .46), q47 = quantile(value, probs = .47), q48 = quantile(value, probs = .48), q49 = quantile(value, probs = .49),
  q50 = quantile(value, probs = .50), q51 = quantile(value, probs = .51), q52 = quantile(value, probs = .52), q53 = quantile(value, probs = .53), q54 = quantile(value, probs = .54),
  q55 = quantile(value, probs = .55), q56 = quantile(value, probs = .56), q57 = quantile(value, probs = .57), q58 = quantile(value, probs = .58), q59 = quantile(value, probs = .59),
  q60 = quantile(value, probs = .60), q61 = quantile(value, probs = .61), q62 = quantile(value, probs = .62), q63 = quantile(value, probs = .63), q64 = quantile(value, probs = .64),
  q65 = quantile(value, probs = .65),
  q70 = quantile(value, probs = .7),
  q75 = quantile(value, probs = .75),
  q80 = quantile(value, probs = .8),
  q85 = quantile(value, probs = .85),
  q90 = quantile(value, probs = .9),
  q95 = quantile(value, probs = .95),
  q100 = max(value))
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

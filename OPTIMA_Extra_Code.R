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

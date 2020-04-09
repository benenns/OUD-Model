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

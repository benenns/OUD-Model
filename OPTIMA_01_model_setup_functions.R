###############################################################################
## OUD Cohort Model (OPTIMA Trial Adaptation)
## Deterministic - Model Functions from: "OPTIMA_00_input_parameter_values.R"
############################## Initial setup ##################################
#rm(list = ls())  # remove any variables in R's memory 
library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(dampack)  # for CEA and calculate ICERs
library(tidyverse)
###############################################################################
####################################
# Set up base to/from model states #
####################################
l_dim_s  <- list() # list of base states
#l_dim_hiv   <- list() # list of base states
#l_dim_death <- list() # list of base states

# Base health states
BASE <- l_dim_s[[1]] <- c("MET1", "MET", "BUP1", "BUP", "ABS", "REL1", "REL", "OD")
n_BASE <- length(BASE)

# Injection/non-injection stratification
INJECT <- l_dim_s[[2]] <- c("NI", "INJ")
n_INJECT <- length(INJECT)

# Episodes (1-3)
EP <-  l_dim_s[[3]] <- c("1", "2", "3")
n_EP <- length(EP)

# HIV status
HIV <- l_dim_s[[4]] <- c("POS", "NEG")

#n_dim_s <- length(l_dim_s)
#n_dim_s
n_t <- (n_age_max - n_age_init) * 12

df_flat <- expand.grid(l_dim_s) #combine all elements together into vector of health states
df_flat <- rename(df_flat, BASE    = Var1, 
                           INJECT  = Var2, 
                           HIV     = Var3, 
                           EP      = Var4)

################ Create index of states to populate transition matrices ################
# All treatment
TX <- df_flat$BASE == "BUP" | df_flat$BASE == "BUP1" | df_flat$BASE == "MET" | df_flat$BASE == "MET1"
# All out-of-treatment (excl ABS)
OOT <- df_flat$BASE == "REL1" | df_flat$BASE == "REL" | df_flat$BASE == "OD"

# All BUP
all_BUP <- df_flat$BASE == "BUP" | df_flat$BASE == "BUP1"
BUP     <- df_flat$BASE == "BUP"
BUP1    <- df_flat$BASE == "BUP1"
# All MET
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

#######
df_n <- unite(df_flat, newCol) # combine columns into one data frame of all health states (1 X [9 states * 2 inj * 2 HIV * 3 Episodes])
v_n_states <- df_n[,1] # convert df into vector
n_states <- length(v_n_states) # total number of health states

#############################################
### Time-dependent survival probabilities ###
#############################################
  # Empty 2-D matrix
  m_TDP <- array(0, dim = c(n_states, n_t),
                    dimnames = list(v_n_states, 0:(n_t - 1)))

  # Probability of remaining in given health state
  # Only populate for allowed transitions
  for(i in 1:n_t){
    #t <- i+1 # Shift BUP, MET, REL ahead 1 month to adjust for BUP1, MET1, REL1
    # Non-injection
      # Episode 1
        m_TDP[EP1 & BUP & NI, i] <- v_TDP_BUP_NI_1[i] # vector of remain probabilities
        m_TDP[EP1 & MET & NI, i] <- v_TDP_MET_NI_1[i]
        m_TDP[EP1 & ABS & NI, i] <- v_TDP_ABS_NI_1[i]
        m_TDP[EP1 & REL & NI, i] <- v_TDP_REL_NI_1[i]
        m_TDP[EP1 & OD & NI, i]  <- v_TDP_OD_NI_1[i]  
      # Episode 2
        m_TDP[EP2 & BUP & NI, i] <- v_TDP_BUP_NI_2[i]
        m_TDP[EP2 & MET & NI, i] <- v_TDP_MET_NI_2[i]
        m_TDP[EP2 & ABS & NI, i] <- v_TDP_ABS_NI_2[i]
        m_TDP[EP2 & REL & NI, i] <- v_TDP_REL_NI_2[i]
        m_TDP[EP2 & OD & NI, i]  <- v_TDP_OD_NI_2[i]
      # Episode 3
        m_TDP[EP3 & BUP & NI, i] <- v_TDP_BUP_NI_3[i]
        m_TDP[EP3 & MET & NI, i] <- v_TDP_MET_NI_3[i]
        m_TDP[EP3 & ABS & NI, i] <- v_TDP_ABS_NI_3[i]
        m_TDP[EP3 & REL & NI, i] <- v_TDP_REL_NI_3[i]
        m_TDP[EP3 & OD & NI, i]  <- v_TDP_OD_NI_3[i]
        
    # Injection
      # Episode 1
        m_TDP[EP1 & BUP & INJ, i] <- v_TDP_BUP_INJ_1[i] # vector of remain probabilities 
        m_TDP[EP1 & MET & INJ, i] <- v_TDP_MET_INJ_1[i]
        m_TDP[EP1 & ABS & INJ, i] <- v_TDP_ABS_INJ_1[i]
        m_TDP[EP1 & REL & INJ, i] <- v_TDP_REL_INJ_1[i]
        m_TDP[EP1 & OD & INJ, i]  <- v_TDP_OD_INJ_1[i]  
      # Episode 2
        m_TDP[EP2 & BUP & INJ, i] <- v_TDP_BUP_INJ_2[i]
        m_TDP[EP2 & MET & INJ, i] <- v_TDP_MET_INJ_2[i]
        m_TDP[EP2 & ABS & INJ, i] <- v_TDP_ABS_INJ_2[i]
        m_TDP[EP2 & REL & INJ, i] <- v_TDP_REL_INJ_2[i]
        m_TDP[EP2 & OD & INJ, i]  <- v_TDP_OD_INJ_2[i]  
      # Episode 3
        m_TDP[EP3 & BUP & INJ, i] <- v_TDP_BUP_INJ_3[i]
        m_TDP[EP3 & MET & INJ, i] <- v_TDP_MET_INJ_3[i]
        m_TDP[EP3 & ABS & INJ, i] <- v_TDP_ABS_INJ_3[i]
        m_TDP[EP3 & REL & INJ, i] <- v_TDP_REL_INJ_3[i]
        m_TDP[EP3 & OD & INJ, i]  <- v_TDP_OD_INJ_3[i]  
  }
m_TDP

# Probability of from-state-exit
m_leave <- 1 - m_TDP

#################
### Mortality ###
#################
# Monthly mortality vector for all model periods by state
# Monthly mortality for each age applied to 12 periods
v_mort <- function(hr, hiv){
  v_mort <- rep((1 - exp(-v_r_mort_by_age[n_age_init:(n_age_max - 1), ] * (1/12) * hr)), each = 12)
  return(v_mort)
}

v_mort_BUP1_NI <- v_mort(hr_BUP1_NI)
v_mort_BUP_NI  <- v_mort(hr_BUP_NI)
v_mort_MET1_NI <- v_mort(hr_MET1_NI)
v_mort_MET_NI  <- v_mort(hr_MET_NI)
v_mort_REL1_NI <- v_mort(hr_REL1_NI)
v_mort_REL_NI  <- v_mort(hr_REL_NI)
v_mort_OD_NI   <- v_mort(hr_OD_NI)
v_mort_ABS_NEG_NI  <- v_mort(hr_ABS_NI)
v_mort_ABS_POS_NI  <- v_mort(hr_ABS_NI_HIV)

v_mort_BUP1_INJ <- v_mort(hr_BUP1_INJ)
v_mort_BUP_INJ  <- v_mort(hr_BUP_INJ)
v_mort_MET1_INJ <- v_mort(hr_MET1_INJ)
v_mort_MET_INJ  <- v_mort(hr_MET_INJ)
v_mort_REL1_INJ <- v_mort(hr_REL1_INJ)
v_mort_REL_INJ  <- v_mort(hr_REL_INJ)
v_mort_OD_INJ   <- v_mort(hr_OD_INJ)
v_mort_ABS_NEG_INJ  <- v_mort(hr_ABS_INJ)
v_mort_ABS_POS_INJ  <- v_mort(hr_ABS_INJ_HIV)

# Create empty mortality matrix (from_states X n_periods)
m_mort <- array(0, dim = c(n_states, n_t),
                dimnames = list(v_n_states, 0:(n_t - 1)))
# Populate mortality matrix (monthly death probability from each state)
for (i in 1:n_t){
  m_mort[BUP1 & NI, i] <- v_mort_BUP1_NEG_NI[i]
  m_mort[BUP & NI, i]  <- v_mort_BUP_NEG_NI[i]
  m_mort[MET1 & NI, i] <- v_mort_MET1_NEG_NI[i]
  m_mort[MET & NI, i]  <- v_mort_MET_NEG_NI[i]
  m_mort[REL1 & NI, i] <- v_mort_REL1_NEG_NI[i]
  m_mort[REL & NI, i]  <- v_mort_REL_NEG_NI[i]
  m_mort[OD & NI, i]   <- v_mort_OD_NEG_NI[i]
  m_mort[ABS & NI, i]  <- v_mort_ABS_NEG_NI[i]
  
  m_mort[BUP1 & INJ, i] <- v_mort_BUP1_NEG_INJ[i]
  m_mort[BUP & INJ, i]  <- v_mort_BUP_NEG_INJ[i]
  m_mort[MET1 & INJ, i] <- v_mort_MET1_NEG_INJ[i]
  m_mort[MET & INJ, i]  <- v_mort_MET_NEG_INJ[i]
  m_mort[REL1 & INJ, i] <- v_mort_REL1_NEG_INJ[i]
  m_mort[REL & INJ, i]  <- v_mort_REL_NEG_INJ[i]
  m_mort[OD & INJ, i]   <- v_mort_OD_NEG_INJ[i]
  m_mort[ABS & INJ, i]  <- v_mort_ABS_NEG_INJ[i]
}

# Alive probability in each period
m_alive <- 1 - m_mort

# Alive probability for every state and time period
# Alive probs in from states
v_alive <- function(from_state, n_mort_per){
  v_alive_prob <- m_alive[from_state, n_mort_per]
return(m_alive)
}

# Count deaths
v_death <- function(from_state, n_mort_per){
  v_death_prob <- m_mort[from_state, n_mort_per]
  return(m_mort)
}


##############################################
### Unconditional transition probabilities ###
##############################################
# Base state transition probabilities (from -> to)

##### Empty 2-D unconditional transition matrix (from states, to states) #######
m_UP <- array(0, dim = c(n_states, n_states),
              dimnames = list(v_n_states, v_n_states))
# Populate unconditional transition probabilities
######## Non-Injection #########
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
m_UP[MET1 & NI & EP1, MET1 & NI & EP1] <- p_MET1_MET1_NI * (1 - p_MET1_MET_NI_1)
m_UP[MET1 & NI & EP2, MET1 & NI & EP2] <- p_MET1_MET1_NI * (1 - p_MET1_MET_NI_2)
m_UP[MET1 & NI & EP3, MET1 & NI & EP3] <- p_MET1_MET1_NI * (1 - p_MET1_MET_NI_3)
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
m_UP[REL1 & NI & EP1, MET1 & NI & EP1] <- p_REL1_MET1_NI * (1 - p_REL1_REL_NI_1)
m_UP[REL1 & NI & EP2, MET1 & NI & EP2] <- p_REL1_MET1_NI * (1 - p_REL1_REL_NI_2)
m_UP[REL1 & NI & EP3, MET1 & NI & EP3] <- p_REL1_MET1_NI * (1 - p_REL1_REL_NI_3)
m_UP[REL1 & NI & EP1, ABS & NI & EP1] <- p_REL1_ABS_NI * (1 - p_REL1_REL_NI_1)
m_UP[REL1 & NI & EP2, ABS & NI & EP2] <- p_REL1_ABS_NI * (1 - p_REL1_REL_NI_2)
m_UP[REL1 & NI & EP3, ABS & NI & EP3] <- p_REL1_ABS_NI * (1 - p_REL1_REL_NI_3)
m_UP[REL1 & NI & EP1, REL1 & NI & EP1] <- p_REL1_REL1_NI * (1 - p_REL1_REL_NI_1)
m_UP[REL1 & NI & EP2, REL1 & NI & EP2] <- p_REL1_REL1_NI * (1 - p_REL1_REL_NI_2)
m_UP[REL1 & NI & EP3, REL1 & NI & EP3] <- p_REL1_REL1_NI * (1 - p_REL1_REL_NI_3)
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

######## Injection ##########
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
m_UP[MET1 & INJ & EP1, MET1 & INJ & EP1] <- p_MET1_MET1_INJ * (1 - p_MET1_MET_INJ_1)
m_UP[MET1 & INJ & EP2, MET1 & INJ & EP2] <- p_MET1_MET1_INJ * (1 - p_MET1_MET_INJ_2)
m_UP[MET1 & INJ & EP3, MET1 & INJ & EP3] <- p_MET1_MET1_INJ * (1 - p_MET1_MET_INJ_3)
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
m_UP[REL1 & INJ & EP1, MET1 & INJ & EP1] <- p_REL1_MET1_INJ * (1 - p_REL1_REL_INJ_1)
m_UP[REL1 & INJ & EP2, MET1 & INJ & EP2] <- p_REL1_MET1_INJ * (1 - p_REL1_REL_INJ_2)
m_UP[REL1 & INJ & EP3, MET1 & INJ & EP3] <- p_REL1_MET1_INJ * (1 - p_REL1_REL_INJ_3)
m_UP[REL1 & INJ & EP1, ABS & INJ & EP1] <- p_REL1_ABS_INJ * (1 - p_REL1_REL_INJ_1)
m_UP[REL1 & INJ & EP2, ABS & INJ & EP2] <- p_REL1_ABS_INJ * (1 - p_REL1_REL_INJ_2)
m_UP[REL1 & INJ & EP3, ABS & INJ & EP3] <- p_REL1_ABS_INJ * (1 - p_REL1_REL_INJ_3)
m_UP[REL1 & INJ & EP1, REL1 & INJ & EP1] <- p_REL1_REL1_INJ * (1 - p_REL1_REL_INJ_1)
m_UP[REL1 & INJ & EP2, REL1 & INJ & EP2] <- p_REL1_REL1_INJ * (1 - p_REL1_REL_INJ_2)
m_UP[REL1 & INJ & EP3, REL1 & INJ & EP3] <- p_REL1_REL1_INJ * (1 - p_REL1_REL_INJ_3)
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

### Apply transition rules #####
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
# Maintain cycles (initiate cycle i+1 with OOT -> TX)
m_UP[TX & EP1, OOT & EP2] = 0
m_UP[TX & EP2, OOT & EP3] = 0
m_UP[OOT & EP1, TX & EP1] = 0
m_UP[OOT & EP2, TX & EP2] = 0

##### Empty 3-D matrix (from states, to states, time periods)
a_TDP <- array(0, dim = c(n_states, n_states, n_t),
                 dimnames = list(v_n_states, v_n_states, 0:(n_t - 1)))

# Add transitions conditional on exit
for (i in 1:n_t){
 a_TDP[, , i] <- m_UP * m_leave[, i]
}

# Add time-dependent remain probabilities
for (i in 1:n_t){
  for (j in 1:n_states){
 a_TDP[j, j, i] <- m_TDP[j, i]
  } 
}
a_TDP



### Apply transition rules #####
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
# Maintain cycles (initiate cycle i+1 with OOT -> TX)
a_TDP[TX & EP1, OOT & EP2, ] = 0
a_TDP[TX & EP2, OOT & EP3, ] = 0
a_TDP[OOT & EP1, TX & EP1, ] = 0
a_TDP[OOT & EP2, TX & EP2, ] = 0

#######################################################  
### Transitions conditional on leaving and survival ###
#######################################################
# a_P[i, j, k] = Transition probibility array [from, to, time(month)]
# Transition array already populated with zeros, only need to define possible transitions
# Only mortality probability and remain probability changing over time
# Other transitions re-weighted by probability of leaving

# Transition probabilities among all base states (divide by n(strata) if not strata-specific)
#n_strata <-  # number of non-base state strata

##########################
#### Run Markov model ####
##########################
# Create empty initial state vectors
# BUP/MET
v_s_init_BL <- v_s_init_BUP <- v_s_init_MET <- rep(0, n_states) # initialize first trace vector of zeros
names(v_s_init_BL) <- names(v_s_init_BUP) <- names(v_s_init_MET) <- v_n_states

# Create empty initial state matrices
# BUP/MET
#m_M_BL <- m_M_BUP <- m_M_MET <- matrix(0, 
#                                       nrow = (n_t + 1), ncol = n_states, 
#                                       dimnames = list(0:n_t, v_n))

################################
### Set initial state vector ###
################################
# Baseline
v_s_init_BL[BUP]  = 0.5/sum(BUP) # Empirically observed proportions from base state
v_s_init_BL[MET]  = 0.5/sum(MET)

# BUP
v_s_init_BUP[BUP] = 1/sum(BUP) # Set all BUP states equal, and sum to 1
# MET
v_s_init_MET[MET] = 1/sum(MET) # Set all MET states equal, and sum to 1
sum(v_s_init_BUP)


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
}

#################################
### Hawkins sojourn fucnction ###
#################################
# Create transition probability array of zeros (from-states * to-states * periods (TDP) * periods (death))
#a_P <- array(0, dim = c(n_states_from, n_states_to, n_t_per, n_mort_per),
#             dimnames = list(v_n_from, v_n_to, 0:(n_t - 1), 0:(n_mort - 1)))

  # Initialize population
    #MarkovTrace[1,,1] <- c(as.vector(sapply(Init.dist, function(x) c(x, rep(0, (NumberOfEpisodes - 1))))), 0, 0)
    a_M_trace <- array(0, dim = c((n_t + 1), n_states, (n_t + 1)),
                       dimnames = list(0:n_t, v_n_states, 0:n_t))
    a_M_trace[1, , 1] <- v_s_init_BL
    a_M_trace

    # All model time periods
      for(i in 2:(n_t)){
      # Time spent in given health state
      for(j in 1:(i - 1)){
        #state-time-dependent transition probability (j) * age (model-time)-specific mortality (i)
        m_sojourn <- a_P_tdp[, , j] * m_alive[, i]
        
        v_current_state <- as.vector(a_M_trace[i - 1, , j])
        v_same_state <- as.vector(v_current_state * diag(m_sojourn))
        a_M_trace[i, ,j + 1] <- v_same_state 
    
        diag(m_sojourn) <- 0
        
        v_new_state <- as.vector(v_current_state %*% m_sojourn)
        a_M_trace[i,,1] <- v_new_state + a_M_trace[i,,1]
      }
    #out.trace <- list(Trace = MarkovTrace, 
     #                 HIVseroconversion = sero.out)
}
a_M_trace

m_M_trace <- array(0, dim = c(n_t + 1, n_states),
                   dimnames = list(0:n_t, v_n_states))
for (i in 1:n_t){
  m_M_trace[i, ] <- rowSums(a_M_trace[i, ,])
}

v_deaths <- array(0, dim = c(n_t, 1))
for (i in 1:n_t){
  deaths[i,] <- as.vector(1 - rowSums(m_M_trace[i,]))
}
deaths


########################################
### Create aggregated trace matrices ###
########################################
# Iterate over all time periods
for(t in 1:n_t){
  # BUP
  m_M_BUP[t + 1, ] <- m_M_BUP[t, ] %*% a_P_BUP[, , t]
  # MET
  m_M_MET[t + 1, ] <- m_M_MET[t, ] %*% a_P_MET[, , t]
}
m_M_BUP
m_M_MET

# Aggregated trace matrix
function(m_M){
m_M_agg_trace <- cbind(TOT = rowSums(m_M[, ]), # total across all states at each time point
  
                       BUP = rowSums(m_M[, BUP]), # totals in base states (no strat)
                       MET = rowSums(m_M[, MET]),
                       REL = rowSums(m_M[, REL]),
                       OD  = rowSums(m_M[, OD]),
                       ABS = rowSums(m_M[, ABS]),
                       
                       HIV = rowSums(m_M[, HIV]), # total with HIV
                       
                       INJ = rowSums(m_M[, INJ])
                       )
return(m_M_agg_trace)
}

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

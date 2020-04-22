##################################### Initial setup ###########################
#rm(list = ls())  # remove any variables in R's memory 
library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(dampack)  # for CEA and calculate ICERs
library(tidyverse)
###############################################################################

########################
### Model Parameters ###
########################
# (1) Populate list in "parameter space" of input parameters
#     - Use list
# Re-code as read-in from excel sheet

#l_params_all <- list() # list of model inputs

###############################
### Initial characteristics ###
###############################

v_initial_BUP <- read.csv("data/01_initial_BUP.csv")
v_initial_MET <- read.csv("data/01_initial_MET.csv")

n_age_init <- 35 # age at baseline
n_age_max <- 95 # maximum age of follow up

n_t <- (n_age_max - n_age_init) * 12 # modeling time horizon in months
p_discount <- 0.03

p_male_BUP <- 0.35
p_male_MET <- 0.40

p_HIV_POS <- 0.05 # % of HIV-positive individuals
p_HIV_ART <- 0.75 # % of HIV-positive on-ART

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
#HIV <- l_dim_s[[4]] <- c("POS", "NEG")

#n_dim_s <- length(l_dim_s)
#n_dim_s
n_t <- 720

df_flat <- expand.grid(l_dim_s) #combine all elements together into vector of health states
df_flat <- rename(df_flat, BASE    = Var1, 
                           INJECT  = Var2, 
                           #HIV    = Var3, 
                           EP      = Var3)

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
#v_n_from
#v_n_to   <- append(v_n_from, "D") # convert df into vector

n_states <- length(v_n_states) # total number of health states
#n_states_to   <- length(v_n_to)   # total number of health states incl death


# Create transition probability array of zeroes (108 from-states * 108 to-states * 720 months)
#a_P <- array(0, dim = c(n_states_from, n_states_to, n_t),
#             dimnames = list(v_n_from, v_n_to, 0:(n_t - 1)))


#############################################
### Time-dependent survival probabilities ###
#############################################
  # Empty 2-D matrix
  m_TDP <- array(0, dim = c(n_t, n_states_from),
                    dimnames = list(0:(n_t - 1), v_n_from))

  v_TDP_BUP_NI_1 <- c()
  
  # Create dummy data
  m_frailty <- array(0.4, dim = c(n_EP, n_BASE, n_INJECT),
                     dimnames = list(EP, BASE, INJECT))
  m_weibull_scale <- array(0.7, dim = c(n_BASE, n_INJECT),
                           dimnames = list(BASE, INJECT))
  m_weibull_shape <- array(0.7, dim = c(n_BASE, n_INJECT),
                           dimnames = list(BASE, INJECT))
  

  # Time-dependent survival probabilities
  for(i in 1:n_t){
    t <- i+1 # Shift BUP, MET, REL ahead 1 month to adjust for BUP1, MET1, REL1
    # Non-injection
      # Episode 1
        v_TDP_BUP_NI_1[i] <- m_TDP[i, EP1 & BUP & NI] <- exp(m_frailty["1", "BUP", "NI"] * m_weibull_scale["BUP", "NI"] * (((t-1)^m_weibull_shape["BUP", "NI"]) - (t^m_weibull_shape["BUP", "NI"]))) # (survival curve at time i)/(survival curve at time i-1) 
        v_TDP_MET_NI_1[i] <- m_TDP[i, EP1 & MET & NI] <- exp(m_frailty["1", "MET", "NI"] * m_weibull_scale["MET", "NI"] * (((t-1)^m_weibull_shape["MET", "NI"]) - (t^m_weibull_shape["MET", "NI"])))
        v_TDP_ABS_NI_1[i] <- m_TDP[i, EP1 & ABS & NI] <- exp(m_frailty["1", "ABS", "NI"] * m_weibull_scale["ABS", "NI"] * (((i-1)^m_weibull_shape["ABS", "NI"]) - (i^m_weibull_shape["ABS", "NI"])))
        v_TDP_REL_NI_1[i] <- m_TDP[i, EP1 & REL & NI] <- exp(m_frailty["1", "REL", "NI"] * m_weibull_scale["REL", "NI"] * (((t-1)^m_weibull_shape["REL", "NI"]) - (t^m_weibull_shape["REL", "NI"])))
        v_TDP_OD_NI_1[i]  <- m_TDP[i, EP1 & OD & NI]  <- exp(m_frailty["1", "OD" , "NI"] * m_weibull_scale["OD", "NI"] * (((i-1)^m_weibull_shape["OD", "NI"]) - (i^m_weibull_shape["OD", "NI"])))
      # Episode 2
        v_TDP_BUP_NI_2[i] <- m_TDP[i, EP2 & BUP & NI] <- exp(m_frailty["2", "BUP", "NI"] * m_weibull_scale["BUP", "NI"] * (((t-1)^m_weibull_shape["BUP", "NI"]) - (t^m_weibull_shape["BUP", "NI"]))) # (survival curve at time i)/(survival curve at time i-1) 
        v_TDP_MET_NI_2[i] <- m_TDP[i, EP2 & MET & NI] <- exp(m_frailty["2", "MET", "NI"] * m_weibull_scale["MET", "NI"] * (((t-1)^m_weibull_shape["MET", "NI"]) - (t^m_weibull_shape["MET", "NI"])))
        v_TDP_ABS_NI_2[i] <- m_TDP[i, EP2 & ABS & NI] <- exp(m_frailty["2", "ABS", "NI"] * m_weibull_scale["ABS", "NI"] * (((i-1)^m_weibull_shape["ABS", "NI"]) - (i^m_weibull_shape["ABS", "NI"])))
        v_TDP_REL_NI_2[i] <- m_TDP[i, EP2 & REL & NI] <- exp(m_frailty["2", "REL", "NI"] * m_weibull_scale["REL", "NI"] * (((t-1)^m_weibull_shape["REL", "NI"]) - (t^m_weibull_shape["REL", "NI"])))
        v_TDP_OD_NI_2[i]  <- m_TDP[i, EP2 & OD & NI]  <- exp(m_frailty["2", "OD" , "NI"] * m_weibull_scale["OD", "NI"] * (((i-1)^m_weibull_shape["OD", "NI"]) - (i^m_weibull_shape["OD", "NI"])))
      # Episode 3
        v_TDP_BUP_NI_3[i] <- m_TDP[i, EP3 & BUP & NI] <- exp(m_frailty["3", "BUP", "NI"] * m_weibull_scale["BUP", "NI"] * (((t-1)^m_weibull_shape["BUP", "NI"]) - (t^m_weibull_shape["BUP", "NI"]))) # (survival curve at time i)/(survival curve at time i-1) 
        v_TDP_MET_NI_3[i] <- m_TDP[i, EP3 & MET & NI] <- exp(m_frailty["3", "MET", "NI"] * m_weibull_scale["MET", "NI"] * (((t-1)^m_weibull_shape["MET", "NI"]) - (t^m_weibull_shape["MET", "NI"])))
        v_TDP_ABS_NI_3[i] <- m_TDP[i, EP3 & ABS & NI] <- exp(m_frailty["3", "ABS", "NI"] * m_weibull_scale["ABS", "NI"] * (((i-1)^m_weibull_shape["ABS", "NI"]) - (i^m_weibull_shape["ABS", "NI"])))
        v_TDP_REL_NI_3[i] <- m_TDP[i, EP3 & REL & NI] <- exp(m_frailty["3", "REL", "NI"] * m_weibull_scale["REL", "NI"] * (((t-1)^m_weibull_shape["REL", "NI"]) - (t^m_weibull_shape["REL", "NI"])))
        v_TDP_OD_NI_3[i]  <- m_TDP[i, EP3 & OD & NI]  <- exp(m_frailty["3", "OD" , "NI"] * m_weibull_scale["OD", "NI"] * (((i-1)^m_weibull_shape["OD", "NI"]) - (i^m_weibull_shape["OD", "NI"])))
        
    # Injection
      # Episode 1
        v_TDP_BUP_INJ_1[i] <- m_TDP[i, EP1 & BUP & INJ] <- exp(m_frailty["1", "BUP", "INJ"] * m_weibull_scale["BUP", "INJ"] * (((t-1)^m_weibull_shape["BUP", "INJ"]) - (t^m_weibull_shape["BUP", "INJ"]))) # (survival curve at time i)/(survival curve at time i-1) 
        v_TDP_MET_INJ_1[i] <- m_TDP[i, EP1 & MET & INJ] <- exp(m_frailty["1", "MET", "INJ"] * m_weibull_scale["MET", "INJ"] * (((t-1)^m_weibull_shape["MET", "INJ"]) - (t^m_weibull_shape["MET", "INJ"])))
        v_TDP_ABS_INJ_1[i] <- m_TDP[i, EP1 & ABS & INJ] <- exp(m_frailty["1", "ABS", "INJ"] * m_weibull_scale["ABS", "INJ"] * (((i-1)^m_weibull_shape["ABS", "INJ"]) - (i^m_weibull_shape["ABS", "INJ"])))
        v_TDP_REL_INJ_1[i] <- m_TDP[i, EP1 & REL & INJ] <- exp(m_frailty["1", "REL", "INJ"] * m_weibull_scale["REL", "INJ"] * (((t-1)^m_weibull_shape["REL", "INJ"]) - (t^m_weibull_shape["REL", "INJ"])))
        v_TDP_OD_INJ_1[i]  <- m_TDP[i, EP1 & OD & INJ]  <- exp(m_frailty["1", "OD" , "INJ"] * m_weibull_scale["OD", "INJ"] * (((i-1)^m_weibull_shape["OD", "INJ"]) - (i^m_weibull_shape["OD", "INJ"])))
      # Episode 2
        v_TDP_BUP_INJ_2[i] <- m_TDP[i, EP2 & BUP & INJ] <- exp(m_frailty["2", "BUP", "INJ"] * m_weibull_scale["BUP", "INJ"] * (((t-1)^m_weibull_shape["BUP", "INJ"]) - (t^m_weibull_shape["BUP", "INJ"]))) # (survival curve at time i)/(survival curve at time i-1) 
        v_TDP_MET_INJ_2[i] <- m_TDP[i, EP2 & MET & INJ] <- exp(m_frailty["2", "MET", "INJ"] * m_weibull_scale["MET", "INJ"] * (((t-1)^m_weibull_shape["MET", "INJ"]) - (t^m_weibull_shape["MET", "INJ"])))
        v_TDP_ABS_INJ_2[i] <- m_TDP[i, EP2 & ABS & INJ] <- exp(m_frailty["2", "ABS", "INJ"] * m_weibull_scale["ABS", "INJ"] * (((i-1)^m_weibull_shape["ABS", "INJ"]) - (i^m_weibull_shape["ABS", "INJ"])))
        v_TDP_REL_INJ_2[i] <- m_TDP[i, EP2 & REL & INJ] <- exp(m_frailty["2", "REL", "INJ"] * m_weibull_scale["REL", "INJ"] * (((t-1)^m_weibull_shape["REL", "INJ"]) - (t^m_weibull_shape["REL", "INJ"])))
        v_TDP_OD_INJ_2[i]  <- m_TDP[i, EP2 & OD & INJ]  <- exp(m_frailty["2", "OD" , "INJ"] * m_weibull_scale["OD", "INJ"] * (((i-1)^m_weibull_shape["OD", "INJ"]) - (i^m_weibull_shape["OD", "INJ"])))
      # Episode 3
        v_TDP_BUP_INJ_3[i] <- m_TDP[i, EP3 & BUP & INJ] <- exp(m_frailty["3", "BUP", "INJ"] * m_weibull_scale["BUP", "INJ"] * (((t-1)^m_weibull_shape["BUP", "INJ"]) - (t^m_weibull_shape["BUP", "INJ"]))) # (survival curve at time i)/(survival curve at time i-1) 
        v_TDP_MET_INJ_3[i] <- m_TDP[i, EP3 & MET & INJ] <- exp(m_frailty["3", "MET", "INJ"] * m_weibull_scale["MET", "INJ"] * (((t-1)^m_weibull_shape["MET", "INJ"]) - (t^m_weibull_shape["MET", "INJ"])))
        v_TDP_ABS_INJ_3[i] <- m_TDP[i, EP3 & ABS & INJ] <- exp(m_frailty["3", "ABS", "INJ"] * m_weibull_scale["ABS", "INJ"] * (((i-1)^m_weibull_shape["ABS", "INJ"]) - (i^m_weibull_shape["ABS", "INJ"])))
        v_TDP_REL_INJ_3[i] <- m_TDP[i, EP3 & REL & INJ] <- exp(m_frailty["3", "REL", "INJ"] * m_weibull_scale["REL", "INJ"] * (((t-1)^m_weibull_shape["REL", "INJ"]) - (t^m_weibull_shape["REL", "INJ"])))
        v_TDP_OD_INJ_3[i]  <- m_TDP[i, EP3 & OD & INJ]  <- exp(m_frailty["3", "OD" , "INJ"] * m_weibull_scale["OD", "INJ"] * (((i-1)^m_weibull_shape["OD", "INJ"]) - (i^m_weibull_shape["OD", "INJ"])))
  }


########################## 
### HIV Seroconversion ###
##########################  
# Only applied to injection
p_BUP1_HIV <- p_hiv[9, 2]
p_BUP_HIV  <- p_hiv[10, 2]
p_MET1_HIV <- p_hiv[11, 2]
p_MET_HIV  <- p_hiv[12, 2]
p_REL1_HIV <- p_hiv[13, 2]
p_REL_HIV  <- p_hiv[14, 2]
p_OD_HIV   <- p_hiv[15, 2]
p_ABS_HIV  <- p_hiv[16, 2]
  
  
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

# Monthly mortality for all model periods by state
# Monthly mortality for each age applied to 12 periods
v_mort <- function(hr, hiv){
  v_mort <- rep((1 - exp(-v_r_mort_by_age[n_age_init:(n_age_max - 1), ] * (1/12) * hr)), each = 12)
  return(v_mort)
}

v_mort_BUP1_NEG_NI <- v_mort(hr_BUP1_NI)
v_mort_BUP_NEG_NI  <- v_mort(hr_BUP_NI)
v_mort_MET1_NEG_NI <- v_mort(hr_MET1_NI)
v_mort_MET_NEG_NI  <- v_mort(hr_MET_NI)
v_mort_REL1_NEG_NI <- v_mort(hr_REL1_NI)
v_mort_REL_NEG_NI  <- v_mort(hr_REL_NI)
v_mort_OD_NEG_NI   <- v_mort(hr_OD_NI)
v_mort_ABS_NEG_NI  <- v_mort(hr_ABS_NI)

v_mort_BUP1_NEG_INJ <- v_mort(hr_BUP1_INJ)
v_mort_BUP_NEG_INJ  <- v_mort(hr_BUP_INJ)
v_mort_MET1_NEG_INJ <- v_mort(hr_MET1_INJ)
v_mort_MET_NEG_INJ  <- v_mort(hr_MET_INJ)
v_mort_REL1_NEG_INJ <- v_mort(hr_REL1_INJ)
v_mort_REL_NEG_INJ  <- v_mort(hr_REL_INJ)
v_mort_OD_NEG_INJ   <- v_mort(hr_OD_INJ)
v_mort_ABS_NEG_INJ  <- v_mort(hr_ABS_INJ)

# HIV Positive
# HIV-specific mortality applied only to abstinence

# Populate mortality matrix (monthly death probability from each state)
# Create empty mortality matrix (from_states X n_periods)
m_mort <- array(0, dim = c(n_states, n_t),
                dimnames = list(v_n, 0:(n_t - 1)))

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

# Alive probability for every state and time period
# Multiply by probs in every state
m_alive <- 1 - m_mort



##############################################
### Unconditional transition probabilities ###
##############################################
# Base state transition probabilities (from -> to)
# Used to re-weight probability of leaving after subtraction of remain + death

######## Non-Injection #########
# From BUP1
p_BUP1_BUP_NI_1  <- exp(m_frailty["1", "BUP", "NI"] * m_weibull_scale["BUP", "NI"] * (((0)^m_weibull_shape["BUP", "NI"]) - (1^m_weibull_shape["BUP", "NI"]))) # Transition from BUP1 -> BUP = Month-1(BUP) -> Month-2(BUP)
p_BUP1_BUP_NI_2  <- exp(m_frailty["2", "BUP", "NI"] * m_weibull_scale["BUP", "NI"] * (((0)^m_weibull_shape["BUP", "NI"]) - (1^m_weibull_shape["BUP", "NI"])))
p_BUP1_BUP_NI_3  <- exp(m_frailty["3", "BUP", "NI"] * m_weibull_scale["BUP", "NI"] * (((0)^m_weibull_shape["BUP", "NI"]) - (1^m_weibull_shape["BUP", "NI"])))
p_BUP1_MET1_NI <- 0.15
p_BUP1_ABS_NI  <- 0.23
p_BUP1_REL1_NI <- 0.18
p_BUP1_OD_NI   <- (1 - 0.15 - 0.23 - 0.18)

# From BUP
p_BUP_MET1_NI <- 0.05
p_BUP_ABS_NI  <- 0.23
p_BUP_REL1_NI <- 0.23
p_BUP_OD_NI   <- (1 - 0.05 - 0.23 - 0.23)

# From MET1
p_MET1_MET_NI_1  <- exp(m_frailty["1", "MET", "NI"] * m_weibull_scale["MET", "NI"] * (((0)^m_weibull_shape["MET", "NI"]) - (1^m_weibull_shape["MET", "NI"]))) # Transition from MET1 -> MET = Month-1(MET) -> Month-2(MET)
p_MET1_MET_NI_2  <- exp(m_frailty["2", "MET", "NI"] * m_weibull_scale["MET", "NI"] * (((0)^m_weibull_shape["MET", "NI"]) - (1^m_weibull_shape["MET", "NI"])))
p_MET1_MET_NI_3  <- exp(m_frailty["3", "MET", "NI"] * m_weibull_scale["MET", "NI"] * (((0)^m_weibull_shape["MET", "NI"]) - (1^m_weibull_shape["MET", "NI"])))
p_MET1_BUP1_NI <- 0.15
p_MET1_ABS_NI  <- 0.23
p_MET1_REL1_NI <- 0.18
p_MET1_OD_NI   <- (1 - 0.15 - 0.23 - 0.18) 
    
# From MET
p_MET_BUP1_NI <-
p_MET_ABS_NI  <-
p_MET_REL1_NI <-
p_MET_OD_NI   <-
  
# From ABS
p_ABS_MET1_NI <-
p_ABS_BUP1_NI <-
p_ABS_REL1_NI <-
p_ABS_OD_NI   <-
  
# From REL1
p_REL1_MET1_NI <-
p_REL1_BUP1_NI <-
p_REL1_ABS_NI  <-
p_REL1_OD_NI   <-
  
# From REL
p_REL_MET1_NI <-
p_REL_BUP1_NI <-
p_REL_ABS_NI  <-
p_REL_OD_NI   <-
  
# From OD
p_OD_MET1_NI <-
p_OD_BUP1_NI <-
p_OD_ABS_NI  <-
p_OD_OD_NI   <-

  
######## Injection ##########
# From BUP1
p_BUP1_MET1_INJ <-
p_BUP1_ABS_INJ  <-
p_BUP1_REL1_INJ <-
p_BUP1_OD_INJ   <-
  
# From BUP
p_BUP_MET1_INJ <-
p_BUP_ABS_INJ  <-
p_BUP_REL1_INJ <-
p_BUP_OD_INJ   <-
  
# From MET
p_MET_BUP1_INJ <-
p_MET_ABS_INJ  <-
p_MET_REL1_INJ <-
p_MET_OD_INJ   <-  
  
# From MET
p_MET_BUP1_INJ <-
p_MET_ABS_INJ  <-
p_MET_REL1_INJ <-
p_MET_OD_INJ   <-
  
# From ABS
p_ABS_MET1_INJ <-
p_ABS_BUP1_INJ <-
p_ABS_REL1_INJ <-
p_ABS_OD_INJ   <-
  
# From REL
p_REL_MET1_INJ <-
p_REL_BUP1_INJ <-
p_REL_ABS_INJ  <-
p_REL_OD_INJ   <-  
  
  
  
##### Empty 2-D matrix (from states X to states) #######
  m_UP <- array(0, dim = c(n_states, n_states),
                dimnames = list(v_n_states, v_n_states))






### Apply transition rules #####
# Episode rules
# Disallowed transitions
m_UP[EP1, EP3] = 0
m_UP[EP2, EP1] = 0
m_UP[EP3, EP1] = 0
m_UP[EP3, EP2] = 0

# Conditional transitions
# Maintain cycles (initiate cycle 2 with OOT -> TX)
m_UP[TX & EP1, OOT & EP2] = 0
m_UP[TX & EP2, OOT & EP3] = 0
m_UP[OOT & EP1, TX & EP1] = 0
m_UP[OOT & EP2, TX & EP2] = 0


  m_UP[, "D"] <- 69

#######################################################  
### Transitions conditional on leaving and survival ###
#######################################################
# a_P[i, j, k] = Transition probibility array [from, to, time(month)]
# Transition array already populated with zeros, only need to define possible transitions
# Only mortality probability and remain probability changing over time
# Other transitions re-weighted by probability of leaving

# Transition probabilities among all base states (divide by n(strata) if not strata-specific)
#n_strata <-  # number of non-base state strata
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
    
    
    
  death.monthly <- Alive.Probs.input[DeathPer, , drop = FALSE]
  Alive.Probs <- do.call(rbind, lapply(From.states, function(x) (1 - death.monthly[,x]) * Probs[x,]))
  rownames(Alive.Probs) <- From.states




#Sample code
#a <- c(1,2,3)
#a

#b <- rep(a, length.out = 108)
#b

##########################
#### Run Markov model ####
##########################
# Create empty initial state vectors
# BUP/MET
v_s_init_BL <- v_s_init_BUP <- v_s_init_MET <- rep(0, n_states_to) # initialize first trace vector of zeros
names(v_s_init_BL) <- names(v_s_init_BUP) <- names(v_s_init_MET) <- v_n_to

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
a_P <- array(0, dim = c(n_states_from, n_states_to, n_t, n_t),
             dimnames = list(v_n_from, v_n_to, 0:(n_t - 1), 0:(n_t - 1)))

  # Initialize population
    #MarkovTrace[1,,1] <- c(as.vector(sapply(Init.dist, function(x) c(x, rep(0, (NumberOfEpisodes - 1))))), 0, 0)
    a_M_trace <- array(0, dim = c((n_t + 1), n_states, (n_t + 1)),
                       dimnames = list(0:n_t, v_n_to, 0:n_t))
    a_M_trace[1, , 1] <- v_s_init_BL
    a_M_trace

    # All model time periods
    for (n_mort_per in 1:nrow(death)){
      for(i in 2:(n_t)){
        n_per <- i #Period counter
        
        n_mort_per <- min(which(Death[,1] == (n_age_init + ((i - 2) / 12))))
      # track mortality separately and re-calculate alive prob at each cycle
      # need to grab p_alive and p_dead for each overall period i from m_mort
        
        
      # Time spent in given health state
      for(j in 1:(i - 1)){
        m_sojourn <- a_P[, , j, i] #j tracks time in current state, i tracks overall periods for mortality
        
        
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
    }  


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

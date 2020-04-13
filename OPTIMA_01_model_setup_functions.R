##################################### Initial setup ###########################
rm(list = ls())  # remove any variables in R's memory 
library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(dampack)  # for CEA and calculate ICERs
library(tidyverse)
###############################################################################

#######################
# Set up model states #
#######################
l_dim_s <- list() # list of all states

# Base health states
BASE <- l_dim_s[[1]] <- c("MET1", "MET", "BUP1", "BUP", "ABS", "REL1", "REL", "OD", "D")
n_BASE <- length(BASE)
# Injection/non-injection stratification
INJ <- l_dim_s[[2]] <- c("NI", "INJ")
# HIV status
HIV <- l_dim_s[[3]] <- c("POS", "NEG")
# Episode number (1-3)
EP <-  l_dim_s[[4]] <- c("1", "2", "3")
n_EP <- length(EP)



n_dim_s <- length(l_dim_s)
#n_dim_s

df_flat <- expand.grid(l_dim_s) #combine all elements together into vector of health states
df_flat <- rename(df_flat, BASE = Var1, 
                           INJ  = Var2, 
                           HIV  = Var3, 
                           EP   = Var4)

df_n <- unite(df_flat, newCol) # combine columns into one data frame of all health states (1 X [9 states * 2 inj * 2 HIV * 3 Episodes])

v_n <- df_n[,1] # convert df into vector

n_states <- length(v_n) # total number of health states
#n_states


# Create transition probability array of zeroes (108 from-states * 108 to-states * 720 months)
a_P <- array(0, dim = c(n_states, n_states, n_t),
             dimnames = list(v_n, v_n, 0:(n_t - 1)))

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
# Create index of states to apply above transition rules
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
INJ <- df_flat$INJ == "INJ"
NI <- df_flat$INJ == "NI"

# Episodes
EP1 <- df_flat$EP == "1"
EP2 <- df_flat$EP == "2"
EP3 <- df_flat$EP == "3"

# All death
D <- df_flat$BASE == "D"



########################
### Model Parameters ###
########################
# (1) Populate list in "parameter space" of input parameters
#     - Use list
# Re-code as read-in from excel sheet

l_params_all <- list() # list of model inputs

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

#################
### Mortality ###
#################
# Baseline mortality by age
lt_can_2018 <- read.csv("data/01_all_cause_mortality.csv")
v_r_mort_by_age <- lt_can_2018 %>% 
  filter(Age >= n_age_init & Age < n_age_max) %>%
  select(Total) %>%
  as.matrix()

# Hazard ratios for death probability
hr_s <- read.csv("data/01_death_hr.csv", header = TRUE)
#hr_s
#hr_BUP1 <- hr_s["BUP1", "HR"] figure out why this doesn't work
hr_BUP1 <- hr_s[1, 2]
hr_BUP  <- hr_s[2, 2]
hr_MET1 <- hr_s[3, 2]
hr_MET  <- hr_s[4, 2]
hr_REL1 <- hr_s[5, 2]
hr_REL  <- hr_s[6, 2]
hr_OD   <- hr_s[7, 2]
hr_ABS  <- hr_s[8, 2]

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

# HIV Positive



##################################
### Probabilities of remaining ###
##################################
# Weibull parameters
#l <- 0.08 # scale
#g <- 1.1  # shape
# Weibull function
#p_ <- l * g * (1:n_t)^{g-1}

# Weibull shape
v_wb_shape   <- rep(0.69, n_BASE)
names(v_wb_shape) <- BASE
v_wb_shape

# Weibull scale
v_wb_scale   <- rep(0.69, n_BASE)
names(v_wb_scale) <- BASE
#v_wb_scale

# Frailty for episodes 2-3 (all frailties for ep1 are 1)
v_frailty <- rep(1, (n_BASE * n_EP))
names(v_frailty) <- c("MET1_1", "MET_1", "BUP1_1", "BUP_1", "ABS_1", "REL1_1", "REL_1", "OD_1", "D_1",
                      "MET1_2", "MET_2", "BUP1_2", "BUP_2", "ABS_2", "REL1_2", "REL_2", "OD_2", "D_2",
                      "MET1_3", "MET_3", "BUP1_3", "BUP_3", "ABS_3", "REL1_3", "REL_3", "OD_3", "D_3")
a <- v_frailty[EP1 & BUP]
b <- v_wb_scale[BUP]
c <- v_wb_shape["BUP"]
a
b
c
  
#######################################################  
### Transitions conditional on leaving and survival ###
#######################################################
# trans_prob = p * (1 - remain - death)

# (2) Translate parameters into model "state space"
# Remain in health states
# Probability of remaining in BUP, MET and REL == Probability of transitioning from "BUP1 -> BUP" or "BUP -> BUP"
v_BUP_remain <-
v_MET_remain <-



# Mortality in OD
a_P[OD, D]

for(i in 1:n_t){
a[i, ] <- exp(v_frailty["BUP_1"]*v_wb_scale["BUP"]*((((i+1)-1)^v_wb_shape["BUP"])-((i+1)^v_wb_shape["BUP"])))
#a <- v_frailty
}
a

############################################
### Fill in transition probability array ###
############################################
# a_P[i, j, k] = Transition probibility array [from, to, time(month)]
# Transition array already populated with zeros, only need to define possible transitions
# Only mortality probability and remain probability changing over time

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


# Transition probabilities among all base states (divide by n(strata) if not strata-specific)
n_strata <-  # number of non-base state strata


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


############################################
#### Check if transition array is valid ####
############################################
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

check_transition_probability(a_P, err_stop = err_stop, verbose = verbose)
check_sum_of_transition_array(a_P, n_states, n_t, err_stop = err_stop, verbose = verbose)

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
v_s_init_BL <- v_s_init_BUP <- v_s_init_MET <- rep(0, n_states) # initialize first trace vector of zeros
names(v_s_init_BL) <- names(v_s_init_BUP) <- names(v_s_init_MET) <- v_n

# Create empty initial state matrices
# BUP/MET
m_M_BL <- m_M_BUP <- m_M_MET <- matrix(0, 
                                       nrow = (n_t + 1), ncol = n_states, 
                                       dimnames = list(0:n_t, v_n))

################################
### Set initial state vector ###
################################
# Baseline
v_s_init_BL[BUP]  = 0.5/sum(BUP)
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




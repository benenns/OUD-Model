library(tidyverse)

#######################
# Set up model states #
#######################
l_dim_s <- list()

# Base health states
BASE <- l_dim_s[[1]] <- c("MET1", "MET", "BUP1", "BUP", "ABS", "REL1", "REL", "OD", "D")
# Injection/non-injection stratification
INJ <- l_dim_s[[2]] <- c("NI", "INJ")
# HIV status
HIV <- l_dim_s[[3]] <- c("POS", "NEG")
# Episode number (1-3)
EP <-  l_dim_s[[4]] <- c("1", "2", "3")

#l_dim_s

n_dim_s <- length(l_dim_s)
n_dim_s

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
# All BUP
BUP <- df_flat$BASE == "BUP" | df_flat$BASE == "BUP1"

# All MET
MET <- df_flat$BASE == "MET" | df_flat$BASE == "MET1"

# All HIV status
NEG <- df_flat$HIV == "NEG"
POS <- df_flat$HIV == "POS"

# All injection
INJ <- df_flat$INJ == "INJ"
NI <- df_flat$INJ == "NI"

# All death
D <- df_flat$BASE == "D"


#  BUP1/BUP/MET1/MET episodes
# Episode 1
EP1_TX <- df_flat$BASE == "BUP1" & df_flat$EP == "1" | df_flat$BASE == "BUP" & df_flat$EP == "1" | df_flat$BASE == "MET1" & df_flat$EP == "1" | df_flat$BASE == "MET" & df_flat$EP == "1"
# Episode 2
EP2_TX <- df_flat$BASE == "BUP1" & df_flat$EP == "2" | df_flat$BASE == "BUP" & df_flat$EP == "2" | df_flat$BASE == "MET1" & df_flat$EP == "2" | df_flat$BASE == "MET" & df_flat$EP == "2"
# Episode 3
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

########################
### Model Parameters ###
########################
# (1) Populate list in "parameter space" of input parameters
#     - Use list
# Re-code as read-in from excel sheet
###############################
### Initial characteristics ###
###############################
n_age_init <- 35 # age at baseline
n_age_max <- 95 # maximum age of follow up

n_t <- (n_age_max - n_age_init)*12 # modeling time horizon in months
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
  # filter(Age >= age & Age <= n_age_max) %>%
  select(Total) %>%
  as.matrix()

# HIV Negative
# ABS (same as baseline)
p_ABS_D_age_NEG   <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)]) # YEARLY PROBABILITY

# BUP
p_BUP1_D_age_NEG  <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_BUP)        
p_BUP_D_age_NEG   <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_BUP1)

# MET
p_MET1_D_age_NEG  <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_MET)        
p_MET_D_age_NEG   <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_MET1)



# HIV Positive
# BUP
p_BUP1_D_age_NEG  <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_BUP)        
p_BUP_D_age_NEG   <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_BUP1)

# MET
p_MET1_D_age_NEG  <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_MET)        
p_MET_D_age_NEG   <- 1 - exp(-v_r_mort_by_age[(n_age_init + 1) + 0:(n_t - 1)] * hr_MET1)


##################################
### Probabilities of remaining ###
##################################
# Weibull parameters
l <- 0.08 # scale
g <- 1.1  # shape
# Weibull function
p_ <- l * g * (1:n_t)^{g-1}

# Weibull shape
n_w_shape_BUP <- 
n_w_shape_MET <- 
n_w_shape_ABS <-
n_w_shape_REL <-  
n_w_shape_OD  <-
n_w_shape_D   <- # Stay in death

# Weibull scale
n_w_scale_BUP <- 
n_w_scale_MET <- 
n_w_scale_ABS <-
n_w_scale_REL <-  
n_w_scale_OD  <-
n_w_scale_D   <- # Stay in death


# Frailty for episodes 2-3 (all frailties for ep1 are 1)
n_frailty_BUP_2 <-  
n_frailty_BUP_3 <- 
  
n_frailty_MET_2 <- 
n_frailty_MET_3 <-
  
n_frailty_ABS_2 <-
n_frailty_ABS_3 <-
  
n_frailty_REL_2 <-
n_frailty_REL_3 <-
  
n_frailty_OD_2  <-
n_frailty_OD_3  <-  
  
#######################################################  
### Transitions conditional on leaving and survival ###
#######################################################
# trans_prob = p * (1 - remain - death)


  

# (2) Translate parameters into model "state space"
# Remain in health states
# Probability of remaining in BUP, MET and REL == Probability of transitioning from "BUP1 -> BUP" or "BUP -> BUP"
p_BUP1_BUP <-
p_BUP1_BUP <-



# Mortality in OD
a_P[OD1, D]




############################################
### Fill in transition probability array ###
############################################

v_ %>% left_join(v_state_rules)


# Transition probabilities among all base states (divide by n(strata) if not strata-specific)
n_strata <-  # number of non-base state strata


# General transition rules
a_P[D, D, ] = 1 # remain in death across all strata and all time points

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
a <- c(1,2,3)
a

b <- rep(a, length.out = 108)
b

##########################
#### Run Markov model ####
##########################
# Create initial state vectors
# BUP/MET
v_s_init_BUP <- v_s_init_MET <- rep(0, n_states) # initialize first trace vector of zeros
names(v_s_init_BUP) <- names(v_s_init_BUP) <- v_n
# MET
#v_s_init_MET <- rep(0, n_states) # initialize first trace vector of zeros
#names(v_s_init_MET) <- v_n

# Create initial state matrices
# BUP
m_M_BUP <- matrix(0, 
                  nrow = (n_t + 1), ncol = n_states, 
                  dimnames = list(0:n_t, v_n))
# MET
m_M_MET <- matrix(0, 
                  nrow = (n_t + 1), ncol = n_states, 
                  dimnames = list(0:n_t, v_n))

################################
### Set initial state values ###
################################
# BUP
v_s_init_BUP[BUP] = 1/sum(BUP) # Set all BUP states equal, and sum to 1
# MET
v_s_init_MET[MET] = 1/24 # Set all MET states equal, and sum to 1
sum(v_s_init_BUP)


# Set first row of m.M with the initial state vector
# BUP
m_M_BUP[1, ] <- v_s_init_BUP
# MET
m_M_MET[1, ] <- v_s_init_MET


######################################################
### Create trace matrices for all states over time ###
######################################################
# Iterate over all time periods
for(t in 1:n_t){
  # BUP
  m_M_BUP[t + 1, ] <- m_M_BUP[t, ] %*% a_P_BUP[, , t]
  # MET
  m_M_MET[t + 1, ] <- m_M_MET[t, ] %*% a_P_MET[, , t]
}
m_M_BUP
m_M_MET







######################################
#### State and Transition Rewards ####
######################################


########################### Cost-effectiveness analysis #######################
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




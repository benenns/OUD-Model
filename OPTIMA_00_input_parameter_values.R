##################################################################################
## OUD Cohort Model
## Deterministic - Set parameter values for: "OPTIMA_01_model_setup_functions.R"
##################################################################################
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

#########################
### Survival analysis ###
#########################
# Frailty terms
p_frailty_BUP_NI_1 <- 1
p_frailty_BUP_NI_2 <- 0.725
p_frailty_BUP_NI_3 <- 0.724
  
p_frailty_MET_NI_1 <- 1
p_frailty_MET_NI_2 <- 0.725
p_frailty_MET_NI_3 <- 0.724

p_frailty_ABS_NI_1 <- 1 
p_frailty_ABS_NI_2 <- 1.09
p_frailty_ABS_NI_3 <- 1.08 
  
p_frailty_REL_NI_1 <- 1 
p_frailty_REL_NI_2 <- 1.22 
p_frailty_REL_NI_3 <- 1.35
  
p_frailty_OD_NI_1 <- 1
p_frailty_OD_NI_2 <- 0.930
p_frailty_OD_NI_3 <- 0.910  

p_frailty_BUP_INJ_1 <- 1
p_frailty_BUP_INJ_2 <- 0.961
p_frailty_BUP_INJ_3 <- 0.959

p_frailty_MET_INJ_1 <- 1
p_frailty_MET_INJ_2 <- 0.961
p_frailty_MET_INJ_3 <- 0.959

p_frailty_ABS_INJ_1 <- 1 
p_frailty_ABS_INJ_2 <- 1.09 
p_frailty_ABS_INJ_3 <- 1.08 

p_frailty_REL_INJ_1 <- 1 
p_frailty_REL_INJ_2 <- 1.22 
p_frailty_REL_INJ_3 <- 1.35

p_frailty_OD_INJ_1 <- 1
p_frailty_OD_INJ_2 <- 0.930
p_frailty_OD_INJ_3 <- 0.910   
  
# Weibull scale  
p_weibull_scale_BUP_NI <- 0.153
p_weibull_scale_MET_NI <- 0.153
p_weibull_scale_ABS_NI <- 0.061
p_weibull_scale_REL_NI <- 0.110
p_weibull_scale_OD_NI  <- 0.750
  
p_weibull_scale_BUP_INJ <- 0.184
p_weibull_scale_MET_INJ <- 0.184
p_weibull_scale_ABS_INJ <- 0.089
p_weibull_scale_REL_INJ <- 0.091
p_weibull_scale_OD_INJ  <- 0.750

# Weibull shape
p_weibull_shape_BUP_NI <- 0.613
p_weibull_shape_MET_NI <- 0.613
p_weibull_shape_ABS_NI <- 0.800
p_weibull_shape_REL_NI <- 0.672
p_weibull_shape_OD_NI  <- 0.793
  
p_weibull_shape_BUP_INJ <- 0.623
p_weibull_shape_MET_INJ <- 0.623
p_weibull_shape_ABS_INJ <- 0.797
p_weibull_shape_REL_INJ <- 0.672
p_weibull_shape_OD_INJ  <- 0.793

#############################################
### Time-dependent survival probabilities ###
#############################################
# Initialize vectors
  v_TDP_BUP_NI_1 <- vector()
  v_TDP_MET_NI_1 <- vector()
  v_TDP_ABS_NI_1 <- vector()
  v_TDP_REL_NI_1 <- vector()
  v_TDP_OD_NI_1  <- vector()
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

# Probability of remaining in given health state
for(i in 1:n_t){
  t <- i+1 # Shift BUP, MET, REL ahead 1 month to adjust for BUP1, MET1, REL1
  # Non-injection
  # Episode 1
  v_TDP_BUP_NI_1[i] <- as.vector(exp(p_frailty_BUP_NI_1 * p_weibull_scale_BUP_NI * (((t - 1)^p_weibull_shape_BUP_NI) - (t^p_weibull_shape_BUP_NI)))) # (survival curve at time i)/(survival curve at time i-1) 
  v_TDP_MET_NI_1[i] <- as.vector(exp(p_frailty_MET_NI_1 * p_weibull_scale_MET_NI * (((t - 1)^p_weibull_shape_MET_NI) - (t^p_weibull_shape_MET_NI))))
  v_TDP_ABS_NI_1[i] <- as.vector(exp(p_frailty_ABS_NI_1 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
  v_TDP_REL_NI_1[i] <- as.vector(exp(p_frailty_REL_NI_1 * p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
  v_TDP_OD_NI_1[i]  <- as.vector(exp(p_frailty_OD_NI_1  * p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))
  # Episode 2
  v_TDP_BUP_NI_2[i] <- as.vector(exp(p_frailty_BUP_NI_2 * p_weibull_scale_BUP_NI * (((t - 1)^p_weibull_shape_BUP_NI) - (t^p_weibull_shape_BUP_NI)))) # (survival curve at time i)/(survival curve at time i-1) 
  v_TDP_MET_NI_2[i] <- as.vector(exp(p_frailty_MET_NI_2 * p_weibull_scale_MET_NI * (((t - 1)^p_weibull_shape_MET_NI) - (t^p_weibull_shape_MET_NI))))
  v_TDP_ABS_NI_2[i] <- as.vector(exp(p_frailty_ABS_NI_2 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
  v_TDP_REL_NI_2[i] <- as.vector(exp(p_frailty_REL_NI_2 * p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
  v_TDP_OD_NI_2[i]  <- as.vector(exp(p_frailty_OD_NI_2  * p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))
  # Episode 3
  v_TDP_BUP_NI_3[i] <- as.vector(exp(p_frailty_BUP_NI_3 * p_weibull_scale_BUP_NI * (((t - 1)^p_weibull_shape_BUP_NI) - (t^p_weibull_shape_BUP_NI)))) # (survival curve at time i)/(survival curve at time i-1) 
  v_TDP_MET_NI_3[i] <- as.vector(exp(p_frailty_MET_NI_3 * p_weibull_scale_MET_NI * (((t - 1)^p_weibull_shape_MET_NI) - (t^p_weibull_shape_MET_NI))))
  v_TDP_ABS_NI_3[i] <- as.vector(exp(p_frailty_ABS_NI_3 * p_weibull_scale_ABS_NI * (((i - 1)^p_weibull_shape_ABS_NI) - (i^p_weibull_shape_ABS_NI))))
  v_TDP_REL_NI_3[i] <- as.vector(exp(p_frailty_REL_NI_3 * p_weibull_scale_REL_NI * (((t - 1)^p_weibull_shape_REL_NI) - (t^p_weibull_shape_REL_NI))))
  v_TDP_OD_NI_3[i]  <- as.vector(exp(p_frailty_OD_NI_3  * p_weibull_scale_OD_NI  * (((i - 1)^p_weibull_shape_OD_NI) - (i^p_weibull_shape_OD_NI))))
  
  # Injection
  # Episode 1
  v_TDP_BUP_INJ_1[i] <- as.vector(exp(p_frailty_BUP_INJ_1 * p_weibull_scale_BUP_INJ * (((t - 1)^p_weibull_shape_BUP_INJ) - (t^p_weibull_shape_BUP_INJ)))) # (survival curve at time i)/(survival curve at time i-1) 
  v_TDP_MET_INJ_1[i] <- as.vector(exp(p_frailty_MET_INJ_1 * p_weibull_scale_MET_INJ * (((t - 1)^p_weibull_shape_MET_INJ) - (t^p_weibull_shape_MET_INJ))))
  v_TDP_ABS_INJ_1[i] <- as.vector(exp(p_frailty_ABS_INJ_1 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
  v_TDP_REL_INJ_1[i] <- as.vector(exp(p_frailty_REL_INJ_1 * p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
  v_TDP_OD_INJ_1[i]  <- as.vector(exp(p_frailty_OD_INJ_1  * p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ))))
  # Episode 2
  v_TDP_BUP_INJ_2[i] <- as.vector(exp(p_frailty_BUP_INJ_2 * p_weibull_scale_BUP_INJ * (((t - 1)^p_weibull_shape_BUP_INJ) - (t^p_weibull_shape_BUP_INJ)))) # (survival curve at time i)/(survival curve at time i-1) 
  v_TDP_MET_INJ_2[i] <- as.vector(exp(p_frailty_MET_INJ_2 * p_weibull_scale_MET_INJ * (((t - 1)^p_weibull_shape_MET_INJ) - (t^p_weibull_shape_MET_INJ))))
  v_TDP_ABS_INJ_2[i] <- as.vector(exp(p_frailty_ABS_INJ_2 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
  v_TDP_REL_INJ_2[i] <- as.vector(exp(p_frailty_REL_INJ_2 * p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
  v_TDP_OD_INJ_2[i]  <- as.vector(exp(p_frailty_OD_INJ_2  * p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ))))
  # Episode 3
  v_TDP_BUP_INJ_3[i] <- as.vector(exp(p_frailty_BUP_INJ_3 * p_weibull_scale_BUP_INJ * (((t - 1)^p_weibull_shape_BUP_INJ) - (t^p_weibull_shape_BUP_INJ)))) # (survival curve at time i)/(survival curve at time i-1) 
  v_TDP_MET_INJ_3[i] <- as.vector(exp(p_frailty_MET_INJ_3 * p_weibull_scale_MET_INJ * (((t - 1)^p_weibull_shape_MET_INJ) - (t^p_weibull_shape_MET_INJ))))
  v_TDP_ABS_INJ_3[i] <- as.vector(exp(p_frailty_ABS_INJ_3 * p_weibull_scale_ABS_INJ * (((i - 1)^p_weibull_shape_ABS_INJ) - (i^p_weibull_shape_ABS_INJ))))
  v_TDP_REL_INJ_3[i] <- as.vector(exp(p_frailty_REL_INJ_3 * p_weibull_scale_REL_INJ * (((t - 1)^p_weibull_shape_REL_INJ) - (t^p_weibull_shape_REL_INJ))))
  v_TDP_OD_INJ_3[i]  <- as.vector(exp(p_frailty_OD_INJ_3  * p_weibull_scale_OD_INJ  * (((i - 1)^p_weibull_shape_OD_INJ) - (i^p_weibull_shape_OD_INJ))))
}

##############################################
### Unconditional transition probabilities ###
##############################################
######## Non-Injection #########
# From BUP1
p_BUP1_BUP_NI_1  <- exp(p_frailty_BUP_NI_1 * p_weibull_scale_BUP_NI * (((0)^p_weibull_shape_BUP_NI) - (1^p_weibull_shape_BUP_NI))) # Transition from BUP1 -> BUP = Month-1(BUP) -> Month-2(BUP)
p_BUP1_BUP_NI_2  <- exp(p_frailty_BUP_NI_2 * p_weibull_scale_BUP_NI * (((0)^p_weibull_shape_BUP_NI) - (1^p_weibull_shape_BUP_NI)))
p_BUP1_BUP_NI_3  <- exp(p_frailty_BUP_NI_3 * p_weibull_scale_BUP_NI * (((0)^p_weibull_shape_BUP_NI) - (1^p_weibull_shape_BUP_NI)))
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
p_MET1_MET_NI_1  <- exp(p_frailty_MET_NI_1 * p_weibull_scale_MET_NI * (((0)^p_weibull_shape_MET_NI) - (1^p_weibull_shape_MET_NI)))
p_MET1_MET_NI_2  <- exp(p_frailty_MET_NI_2 * p_weibull_scale_MET_NI * (((0)^p_weibull_shape_MET_NI) - (1^p_weibull_shape_MET_NI)))
p_MET1_MET_NI_3  <- exp(p_frailty_MET_NI_3 * p_weibull_scale_MET_NI * (((0)^p_weibull_shape_MET_NI) - (1^p_weibull_shape_MET_NI)))
p_MET1_BUP1_NI <- 0.15
p_MET1_ABS_NI  <- 0.23
p_MET1_REL1_NI <- 0.18
p_MET1_OD_NI   <- (1 - 0.15 - 0.23 - 0.18) 

# From MET
p_MET_BUP1_NI <- 0.04
p_MET_ABS_NI  <- 0.32
p_MET_REL1_NI <- 0.18
p_MET_OD_NI   <- (1 - 0.04 - 0.32 - 0.18)

# From ABS
p_ABS_REL1_NI <- 0.5
p_ABS_OD_NI   <- 0.5

# From REL1
p_REL1_REL_NI_1  <- exp(p_frailty_REL_NI_1 * p_weibull_scale_REL_NI * (((0)^p_weibull_shape_REL_NI) - (1^p_weibull_shape_REL_NI)))
p_REL1_REL_NI_2  <- exp(p_frailty_REL_NI_2 * p_weibull_scale_REL_NI * (((0)^p_weibull_shape_REL_NI) - (1^p_weibull_shape_REL_NI)))
p_REL1_REL_NI_3  <- exp(p_frailty_REL_NI_3 * p_weibull_scale_REL_NI * (((0)^p_weibull_shape_REL_NI) - (1^p_weibull_shape_REL_NI)))
p_REL1_MET1_NI <- 0.04
p_REL1_BUP1_NI <- 0.32
p_REL1_ABS_NI  <- 0.18
p_REL1_OD_NI   <- (1 - 0.04 - 0.32 - 0.18)

# From REL
p_REL_MET1_NI <- 0.04
p_REL_BUP1_NI <- 0.32
p_REL_ABS_NI  <- 0.18
p_REL_OD_NI   <- (1 - 0.04 - 0.32 - 0.18)

# From OD
p_OD_MET1_NI <- 0.04
p_OD_BUP1_NI <- 0.32
p_OD_ABS_NI  <- 0.18
p_OD_REL1_NI   <- (1 - 0.04 - 0.32 - 0.18)

######## Injection ##########
# From BUP1
p_BUP1_BUP_INJ_1  <- exp(p_frailty_BUP_INJ_1 * p_weibull_scale_BUP_INJ * (((0)^p_weibull_shape_BUP_INJ) - (1^p_weibull_shape_BUP_INJ))) # Transition from BUP1 -> BUP = Month-1(BUP) -> Month-2(BUP)
p_BUP1_BUP_INJ_2  <- exp(p_frailty_BUP_INJ_2 * p_weibull_scale_BUP_INJ * (((0)^p_weibull_shape_BUP_INJ) - (1^p_weibull_shape_BUP_INJ)))
p_BUP1_BUP_INJ_3  <- exp(p_frailty_BUP_INJ_3 * p_weibull_scale_BUP_INJ * (((0)^p_weibull_shape_BUP_INJ) - (1^p_weibull_shape_BUP_INJ)))
p_BUP1_MET1_INJ <- 0.15
p_BUP1_ABS_INJ  <- 0.23
p_BUP1_REL1_INJ <- 0.18
p_BUP1_OD_INJ   <- (1 - 0.15 - 0.23 - 0.18)

# From BUP
p_BUP_MET1_INJ <- 0.05
p_BUP_ABS_INJ  <- 0.23
p_BUP_REL1_INJ <- 0.23
p_BUP_OD_INJ   <- (1 - 0.05 - 0.23 - 0.23)

# From MET1
p_MET1_MET_INJ_1  <- exp(p_frailty_MET_INJ_1 * p_weibull_scale_MET_INJ * (((0)^p_weibull_shape_MET_INJ) - (1^p_weibull_shape_MET_INJ)))
p_MET1_MET_INJ_2  <- exp(p_frailty_MET_INJ_2 * p_weibull_scale_MET_INJ * (((0)^p_weibull_shape_MET_INJ) - (1^p_weibull_shape_MET_INJ)))
p_MET1_MET_INJ_3  <- exp(p_frailty_MET_INJ_3 * p_weibull_scale_MET_INJ * (((0)^p_weibull_shape_MET_INJ) - (1^p_weibull_shape_MET_INJ)))
p_MET1_BUP1_INJ <- 0.15
p_MET1_ABS_INJ  <- 0.23
p_MET1_REL1_INJ <- 0.18
p_MET1_OD_INJ   <- (1 - 0.15 - 0.23 - 0.18) 

# From MET
p_MET_BUP1_INJ <- 0.04
p_MET_ABS_INJ  <- 0.32
p_MET_REL1_INJ <- 0.18
p_MET_OD_INJ   <- (1 - 0.04 - 0.32 - 0.18)

# From ABS
p_ABS_REL1_INJ <- 0.5
p_ABS_OD_INJ   <- 0.5

# From REL1
p_REL1_REL_INJ_1  <- exp(p_frailty_REL_INJ_1 * p_weibull_scale_REL_INJ * (((0)^p_weibull_shape_REL_INJ) - (1^p_weibull_shape_REL_INJ)))
p_REL1_REL_INJ_2  <- exp(p_frailty_REL_INJ_2 * p_weibull_scale_REL_INJ * (((0)^p_weibull_shape_REL_INJ) - (1^p_weibull_shape_REL_INJ)))
p_REL1_REL_INJ_3  <- exp(p_frailty_REL_INJ_3 * p_weibull_scale_REL_INJ * (((0)^p_weibull_shape_REL_INJ) - (1^p_weibull_shape_REL_INJ)))
p_REL1_MET1_INJ <- 0.04
p_REL1_BUP1_INJ <- 0.32
p_REL1_ABS_INJ  <- 0.18
p_REL1_OD_INJ   <- (1 - 0.04 - 0.32 - 0.18)

# From REL
p_REL_MET1_INJ <- 0.04
p_REL_BUP1_INJ <- 0.32
p_REL_ABS_INJ  <- 0.18
p_REL_OD_INJ   <- (1 - 0.04 - 0.32 - 0.18)

# From OD
p_OD_MET1_INJ <- 0.04
p_OD_BUP1_INJ <- 0.32
p_OD_ABS_INJ  <- 0.18
p_OD_REL1_INJ   <- (1 - 0.04 - 0.32 - 0.18)


########################## 
### HIV Seroconversion ###
##########################
# Non-injection
p_sero_BUP1_NI <- 0.000070
p_sero_BUP_NI  <- 0.000070
p_sero_MET1_NI <- 0.000070
p_sero_MET_NI  <- 0.000070
p_sero_REL1_NI <- 0.000926
p_sero_REL_NI  <- 0.000926
p_sero_OD_NI   <- 0.000926
p_sero_ABS_NI  <- 0.000018

# Injection
p_sero_BUP1_INJ <- 0.000234
p_sero_BUP_INJ  <- 0.000234
p_sero_MET1_INJ <- 0.000234
p_sero_MET_INJ  <- 0.000234
p_sero_REL1_INJ <- 0.00309
p_sero_REL_INJ  <- 0.00309
p_sero_OD_INJ   <- 0.00309
p_sero_ABS_INJ  <- 0.000058

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
#hr_s <- read.csv("data/01_death_hr.csv", header = TRUE)
#hr_s
# HIV-negative
hr_BUP1_NI <- 5
hr_BUP_NI  <- 1.332
hr_MET1_NI <- 8
hr_MET_NI  <- 1.332
hr_REL1_NI <- 14.615
hr_REL_NI  <- 3.922
hr_OD_NI   <- 40
hr_ABS_NI  <- 1
hr_ABS_NI_HIV <- 1.25

hr_BUP1_INJ <- 5
hr_BUP_INJ  <- 2.432
hr_MET1_INJ <- 8
hr_MET_INJ  <- 2.432
hr_REL1_INJ <- 26.689
hr_REL_INJ  <- 7.162
hr_OD_INJ   <- 50
hr_ABS_INJ  <- 1
hr_ABS_INJ_HIV  <- 1.47

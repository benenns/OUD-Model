#############
#TEST MODULE#
#############
dat <- as.matrix(read.csv("C:/Users/Ben/Desktop/Book1.csv", header = FALSE))
death <- as.matrix(read.csv("C:/Users/Ben/Desktop/Book2.csv", header = FALSE))
overdose <- as.matrix(read.csv("C:/Users/Ben/Desktop/Book3.csv", header = FALSE))

a_test <- array(0, dim = c(5, 5, 10),
                dimnames = list())
for (i in 1:10){
  a_test[, , i] <- dat[, ]
}

a_test_trace <- a_test_trace_d <- array(0, dim = c(11, 5, 11))
m_death <- array(0, dim = c(5, 10))
m_death <- death
m_overdose <- array(0, dim = c(5, 10))
m_overdose <- overdose

m_alive <- 1 - m_death
m_non_od <- 1 - m_overdose

v_init <- c(1, 0, 0, 0, 0)
a_test_trace[1, , 1] <- v_init

# All model time periods
for(i in 2:10){
  for(j in 1:(i - 1)){
    m_sojourn <- a_test[, , j] * m_alive[, i-1] * m_OD_time_adj[, i-1] # state-time transition matrix at time j, re-weight for model-time varying overdose at each time point
    
    v_current_state <- as.vector(a_test_trace[i - 1, , j]) # all in current state
    
    v_same_state <- as.vector(v_current_state * diag(m_sojourn)) # individuals remaining
    
    a_test_trace[i, ,j + 1] <- v_same_state # apply to next period
    
    diag(m_sojourn) <- 0 # reset remain to 0
    
    v_new_state <- as.vector(v_current_state %*% m_sojourn) # populate new states (excluding remaining)
    
    a_test_trace[i,,1] <- v_new_state + a_test_trace[i,,1] # add new state %'s to matrix
  }
}

m_M_trace <- array(0, dim = c(11, 5))
for (i in 1:10){
  m_M_trace[i, ] <- rowSums(a_test_trace[i, ,])
}
write.csv(m_M_trace,"C:/Users/Ben/Desktop/test_trace.csv", row.names = TRUE)

m_M_trace_d <- array(0, dim = c(11, 5))
for (i in 2:11){
m_M_trace_d[i, ] <- m_M_trace[i-1, ] * m_test_death[, i-1]
}
write.csv(m_M_trace_d,"C:/Users/Benjamin/Desktop/test_trace_d.csv", row.names = TRUE)

gd <- apply(m_M_trace_d, 2, cumsum)
write.csv(gd,"C:/Users/Benjamin/Desktop/test_trace_d.csv", row.names = TRUE)


a <- m_M_trace[1,]
b <- m_test_death[, 1]
a
b
c <- a * b
c



df_init_dist <- read.csv("C:/Users/Benjamin/Documents/GitHub/OUD-Model/Data/init_dist.csv", row.names = 1, header = TRUE) # Initial parameter values
v_init_dist = as.vector(df_init_dist["pe", ])


v_s_init <- rep(0, dim = c(n_states), 
                dimnames = list(v_n_states))

v_s_init <- rep(0, n_states)
names(v_s_init) <- v_n_states

# Set initial state vector
p_INJ <- 0.5
p_HIV_POS <- 0.05
n_t <- 720
# Baseline
v_s_init[BUP1 & EP1] <- v_init_dist["pe", "BUP1"] # Empirically observed proportions from base states
v_s_init[BUP & EP1]  <- v_init_dist["pe", "BUP"]
v_s_init[MET1 & EP1] <- v_init_dist["pe", "MET1"]
v_s_init[MET & EP1]  <- v_init_dist["pe", "MET"]
v_s_init[REL1 & EP1] <- v_init_dist["pe", "REL1"]
v_s_init[REL & EP1]  <- v_init_dist["pe", "REL"]
v_s_init[OD & EP1]   <- v_init_dist["pe", "OD"]
v_s_init[ABS & EP1]  <- v_init_dist["pe", "ABS"]

v_s_init[NI]  <- v_s_init[NI] * (1 - p_INJ)
v_s_init[INJ] <- v_s_init[INJ] * p_INJ

v_s_init[NEG] <- v_s_init[NEG] * (1 - p_HIV_POS)
v_s_init[POS] <- v_s_init[POS] * p_HIV_POS

g_unit <- sum(v_s_init)
write.csv(v_s_init,"C:/Users/Benjamin/Desktop/v_s_init.csv", row.names = TRUE)

# PSA Test Code
df_UP <- read.csv("C:/Users/Ben/Documents/GitHub/OUD-Model/Data/unconditional.csv", row.names = 1, header = TRUE) # Initial parameter values
m_dirichlet_UP = df_UP * 272 # weight matrix by sample size
v_dirichlet_UP_BUP = m_dirichlet_UP["BUP_NI",]
m_BUP_state = rdirichlet(50, c(v_dirichlet_UP_BUP["BUP_NI", "BUPC_NI"], v_dirichlet_UP_BUP["BUP_NI", "MET_NI"], v_dirichlet_UP_BUP["BUP_NI", "METC_NI"], v_dirichlet_UP_BUP["BUP_NI", "ABS_NI"], v_dirichlet_UP_BUP["BUP_NI", "REL_NI"]))
p_BUP_BUPC_NI = m_BUP_state[,1]

library(rBeta2009)

v_b <- matrix(0, nrow = 1, ncol = 5)
g <- ncol(v_b)
v_c <- matrix(0, nrow = 1, ncol = ncol(v_b))

n_sim <- 50

m_mort <- array(0, dim = c(n_sim, n_t),
                dimnames = list(v_n_states, 1:n_t))

m_dir <- as.matrix(rdirichlet(n_sim, c(1,2,3,4,5)))

v_1 <- matrix(1, nrow = n_sim, ncol = 1)

p_OD_MET_INJ  = m_dir[,1]
n_fent_OD = rgamma(n_sim, shape = 2, scale = 45)
u_BUP_NI_NEG  = truncnorm::rtruncnorm(n_sim, mean = 0.8, sd = 0.1, b = 1)


# Function to generate lognormal parameter
location <- function(m = m, s = s){
  log(m^2 / sqrt(s^2 + m^2))
}
shape <- function(m = m, s = s){
  shape <- sqrt(log(1 + (s^2 / m^2)))
}
hr_BUP_NI  = rlnorm(n_sim, location(m = 3, s = 0.5), shape(m = 3, s = 0.5))


df_2 <- data.frame(p_OD_MET_INJ, n_fent_OD, hr_BUP_NI, u_BUP_NI_NEG)

n_sim = 1000 
seed = 3730687 
n_samp = 250
file.death_hr = "data/death_hr.csv"
file.frailty = "data/frailty.csv"
file.weibull_scale = "data/weibull_scale.csv"
file.weibull_shape = "data/weibull_shape.csv"
file.unconditional = "data/unconditional.csv"
file.overdose = "data/overdose.csv"
file.hiv = "data/hiv_sero.csv"
file.hcv = "data/hcv_sero.csv"
file.costs = "data/costs.csv"
file.crime_costs = "data/crime_costs.csv"
file.qalys = "data/qalys.csv"
file.imis_output = "outputs/imis_output.RData"

t = c(3, 2, 1)
g = t

df_3 <- data.frame(
  t <- c(3, 2, 1),
  g <- t
)

df_4 <- data.frame(
  vec1 <- c(3, 2, 1),
  vec2 <- vec1
)

if (!missing(seed)) 
  set.seed(seed)
c_BUP_NI_crime = rgamma(n_sim, shape = df_crime_costs["shape", "BUP"], scale = df_crime_costs["scale", "BUP"])
set.seed(seed)
c_BUP_INJ_crime = rgamma(n_sim, shape = df_crime_costs["shape", "BUP"], scale = df_crime_costs["scale", "BUP"])


df_5 <- data.frame(
c_BUP_NI_crime = rgamma(n_sim, shape = df_crime_costs["shape", "BUP"], scale = df_crime_costs["scale", "BUP"]),
c_BUP_INJ_crime = c_BUP_NI_crime
)

a <- 3
b <- 4
c <- 5
df_6 <- data.frame(a, b, c)

df_x <- data.frame()
for(i in 1:5){
  x <- i
  
  #df_x <- data.frame(matrix(ncol = 1, nrow = 5))
  df_x <- rbind(df_x, x)
}

df_y = data.frame()

# Defining a for loop with 30 iterations

for (i in 1:30) {
  
  output = c(i^3+3, i^2+2, i+1)
  
  # Using rbind() to append the output of one iteration to the dataframe
  df_y = rbind(df_y, output)
}


## CALIBRATION TEST ##
#### Load packages and functions ####
library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(ggridges) # specialized ridge plots
library(tidyverse)
library(lhs)
library(IMIS)

# To-do: Move functions into R package for OUD model
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/calibration_functions.R")

# Load model inputs #
l_params_all <- load_all_params(file.init = "data/init_params.csv",
                                file.init_dist = "data/init_dist.csv", # calibrate on full trial sample (x% bup; x% met)
                                file.mort = "data/all_cause_mortality.csv",
                                file.death_hr = "data/death_hr.csv",
                                file.frailty = "data/frailty.csv",
                                file.weibull_scale = "data/weibull_scale.csv",
                                file.weibull_shape = "data/weibull_shape.csv",
                                file.unconditional = "data/unconditional.csv",
                                file.overdose = "data/overdose.csv", # includes calibration-related parameters
                                file.fentanyl = "data/fentanyl.csv",
                                file.hiv = "data/hiv_sero.csv",
                                file.hcv = "data/hcv_sero.csv",
                                file.costs = "data/costs.csv",
                                file.crime_costs = "data/crime_costs.csv",
                                file.qalys = "data/qalys.csv")

# Load calibration inputs #
v_cali_param_names <- c("'Overdose rate (treatment)'",
                        "'Overdose rate (treatment + opioid)'",
                        "'Overdose rate (active opioid)'", 
                        "'Overdose rate (inactive opioid)'",
                        "'First month mult (treatment)'",
                        "'First month mult (treatment + opioid)'",
                        "'First month mult (active opioid)'",
                        "'Injection mult'",
                        "'Fentanyl mult'",
                        "'Fatal overdose rate'")
v_par1 <- c(n_TX_OD_shape   = l_params_all$n_TX_OD_shape,
            n_TXC_OD_shape  = l_params_all$n_TXC_OD_shape,
            n_REL_OD_shape   = l_params_all$n_REL_OD_shape,
            n_ABS_OD_low    = l_params_all$n_ABS_OD_low,
            n_TX_OD_mult_shape  = l_params_all$n_TX_OD_mult_shape,
            n_TXC_OD_mult_shape = l_params_all$n_TXC_OD_mult_shape,
            n_REL_OD_mult_shape = l_params_all$n_REL_OD_mult_shape,
            n_INJ_OD_mult_shape = l_params_all$n_INJ_OD_mult_shape,
            n_fent_OD_mult_shape = l_params_all$n_fent_OD_mult_shape,
            n_fatal_OD_shape    = l_params_all$n_fatal_OD_shape) # lower bound estimate for each param
v_par2 <- c(n_TX_OD_scale   = l_params_all$n_TX_OD_scale,
            n_TXC_OD_scale  = l_params_all$n_TXC_OD_scale,
            n_REL_OD_scale   = l_params_all$n_REL_OD_scale,
            n_ABS_OD_high    = l_params_all$n_ABS_OD_high,
            n_TX_OD_mult_scale  = l_params_all$n_TX_OD_mult_scale,
            n_TXC_OD_mult_scale = l_params_all$n_TXC_OD_mult_scale,
            n_REL_OD_mult_scale = l_params_all$n_REL_OD_mult_scale,
            n_INJ_OD_mult_scale = l_params_all$n_INJ_OD_mult_scale,
            n_fent_OD_mult_scale = l_params_all$n_fent_OD_mult_scale,
            n_fatal_OD_scale    = l_params_all$n_fatal_OD_scale)

#### Load calibration targets ####
l_cali_targets <- list(ODF = read.csv(file = "data/cali_target_odf.csv", header = TRUE),
                       ODN = read.csv(file = "data/cali_target_odn.csv", header = TRUE))

# Max calibration periods
n_cali_max_per <- max(c(l_cali_targets$ODF$Time, l_cali_targets$ODN$Time))

#test <- calibration_out(v_params_calib = v_params_calib, l_params_all = l_params_all) 1-EXP(-M21)

test2 <- sample.prior(n_samp = 100)

n_t <- 720
p_fent_exp_base <- 0.50
n_fent_growth_rate <- 0.0131
v_fent_exp_rate <- rep(0, n_t)
#test <- -log(1- p_fent_exp_base)

for(i in 2:n_t){
  v_fent_exp_rate[1] <- -log(1- p_fent_exp_base)
  v_fent_exp_rate[i] <- v_fent_exp[i-1] + n_fent_growth_rate
}
v_fent_exp_prob <- 1 - exp(-v_fent_exp_rate)
v_fent_exp_prob

n_t <- 10
test <- rep(0, n_t)
for(i in 1){
  test[i] <- 2
}
for(i in 2:6){
  test[i] <- 4
}
for(i in 7:n_t){
  test[i] <- 6
}

v_fent_exp_prob <- c(0.2, 0.4, 0.5)
v_fent_exp_prob[1]

# Calibration log-liklihood test
a <- c(0.0245, 0.0300, 0.0275)
b <- c(0.210, 0.210, 0.260)
model_res <- list(a = a, 
                  b = b)

#for(j in 1:n_samp) { # j=1
#  jj <- tryCatch( { 
    ###   Run model for parameter set "v_params" ###
    l_model_res <- calibration_out(v_params_calib = v_params[j, ], 
                                   l_params_all = l_params_all)
    
    ###  Calculate log-likelihood of model outputs to targets  ###
    ## TARGET 1: Fatal overdoses ("fatal_overdose")
    ## Normal log-likelihood  
    v_llik[1, "Fatal Overdoses"] <- sum(dnorm(x = l_cali_targets$ODF$pe,
                                              mean = model_res$a,
                                              sd = l_cali_targets$ODF$se,
                                              log = T))
    ## TARGET 2: Non-fatal overdoses ("overdose")
    ## Normal log-likelihood
    v_llik[1, "Overdoses"] <- sum(dnorm(x = l_cali_targets$ODN$pe,
                                        mean = model_res$a,
                                        sd = l_cali_targets$ODN$se,
                                        log = T))
    ## can give different targets different weights
    v_weights <- c(1, 0.75)
    #v_weights <- rep(1, n_target) # currently weight fatal overdoses 1:1 to overall overdoses
    ## weighted sum
    v_llik_overall[j] <- v_llik[j, ] %*% v_weights
#  }, error = function(e) NA) 
#  if(is.na(jj)) { v_llik_overall <- -Inf }
#} ## End loop over sampled parameter sets
## return GOF
    x <- l_cali_targets$ODF$pe
    test <- model_res$a
    sd <- l_cali_targets$ODF$se
    odf_weight <- l_cali_targets$ODF$weight
    odn_weight <- l_cali_targets$ODN$weight
    
try <-    dnorm(x = l_cali_targets$ODF$pe,
              mean = model_res$a,
              sd = l_cali_targets$ODF$se,
              log = T)
try2 <-    sum(dnorm(x = l_cali_targets$ODF$pe,
                mean = model_res$a,
                sd = l_cali_targets$ODF$se,
                log = T) * l_cali_targets$ODN$weight)

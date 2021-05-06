#############
#TEST MODULE#
#############
dat <- as.matrix(read.csv("C:/Users/Benjamin/Desktop/Book1.csv", header = FALSE))
death <- as.matrix(read.csv("C:/Users/Benjamin/Desktop/Book2.csv", header = FALSE))

a_test <- array(0, dim = c(5, 5, 10),
                dimnames = list())
for (i in 1:10){
  a_test[, , i] <- dat[, ]
}

a_test_trace <- a_test_trace_d <- array(0, dim = c(11, 5, 11))
m_test <- array(0, dim = c(5, 10))
m_test <- death

m_test_death <- 1 - m_test

v_init <- c(1, 0, 0, 0, 0)
a_test_trace[1, , 1] <- v_init

# All model time periods
for(i in 2:10){
  for(j in 1:(i - 1)){
    m_sojourn <- a_test[, , j] * m_test[, i - 1]
    v_current_state <- as.vector(a_test_trace[i - 1, , j])
    v_same_state <- as.vector(v_current_state * diag(m_sojourn))
    a_test_trace[i, ,j + 1] <- v_same_state 
    diag(m_sojourn) <- 0
    v_new_state <- as.vector(v_current_state %*% m_sojourn)
    a_test_trace[i,,1] <- v_new_state + a_test_trace[i,,1]
  }
}

m_M_trace <- array(0, dim = c(11, 5))
for (i in 1:10){
  m_M_trace[i, ] <- rowSums(a_test_trace[i, ,])
}
write.csv(m_M_trace,"C:/Users/Benjamin/Desktop/test_trace.csv", row.names = TRUE)

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




v_test <- rep(1, 721)
v_test2 <- rep(0.5, 721)
one_year <- sum(v_test[1:12])
five_year <- sum(v_test[1:60])

v_test3 <- (v_test * 0.6) + (v_test2 * 0.4)


a_TDP <- array(0, dim = c(144, 144, 3120))
a_TDP2 <- array(0, dim = c(144, 144, 3120))
a_TDP3 <- array(0, dim = c(3120, 144, 3120))

l_dim_s  <- list() # list of health states
# Base health states
BASE <- l_dim_s[[1]] <- c("MET", "BUP", "ABS", "REL", "ODN", "ODF")
# Injection/non-injection stratification
INJECT <- l_dim_s[[2]] <- c("NI", "INJ")
# Episodes (1-3)
EP <-  l_dim_s[[3]] <- c("1", "2", "3")
# HIV/HCV status
SERO <- l_dim_s[[4]] <- c("NEG", "HIV", "HCV", "COI")

n_BASE <- length(BASE)
n_INJECT <- length(INJECT)
n_EP <- length(EP)
n_SERO <- length(SERO)

te <- n_BASE * n_INJECT * n_EP
te


params_updated1 <- c(p_BUP_OD_NI  = l_params_all$p_BUP_OD_NI, 
                     p_MET_OD_NI  = l_params_all$p_MET_OD_NI, 
                     p_REL_OD_NI  = l_params_all$p_REL_OD_NI, 
                     p_ABS_OD_NI  = l_params_all$p_ABS_OD_NI, 
                     p_BUP_OD_INJ = l_params_all$p_BUP_OD_INJ, 
                     p_MET_OD_INJ = l_params_all$p_MET_OD_INJ, 
                     p_REL_OD_INJ = l_params_all$p_REL_OD_INJ, 
                     p_ABS_OD_INJ = l_params_all$p_ABS_OD_INJ)

params_updated2 <- c(p_BUP_OD_NI = 0.6, 
                    p_MET_OD_NI = 0.6, 
                    p_REL_OD_NI = 0.6, 
                    p_ABS_OD_NI = 0.6, 
                    p_BUP_OD_INJ = 0.6, 
                    p_MET_OD_INJ = 0.6, 
                    p_REL_OD_INJ = 0.6, 
                    p_ABS_OD_INJ = 0.6)

p_BUP_OD_NI2 = 0.6

if (typeof(params_updated)!="list"){
  params_updated <- split(unname(params_updated),names(params_updated)) #convert the named vector to a list
}

est1 <- l_params_all$p_BUP_OD_NI
est2 <- p_BUP_OD_NI2

###   Run model for parameter set "v_params" ###
l_model_res <- calibration_out(v_params_calib = v_params[1, ], 
                               l_params_all = l_params_all)

###  Calculate log-likelihood of model outputs to targets  ###
## TARGET 1: Fatal overdoses ("fatal_overdose")
## Normal log-likelihood  
v_llik[j, "fatal_overdose"] <- sum(dnorm(x = l_cali_targets$ODF$pe,
                                         mean = l_model_res$fatal_overdose,
                                         sd = l_cali_targets$ODF$se,
                                         log = T))

## TARGET 2: Non-fatal overdoses ("overdose")
## Normal log-likelihood
v_llik <- sum(dnorm(x = l_cali_targets$ODN$pe,
                                   mean = l_model_res$overdose,
                                   sd = l_cali_targets$ODN$se,
                                   log = T))




sample.prior <- function(n_samp,
                         v_param_names = c("p_BUP_OD_NI", 
                                           "p_MET_OD_NI", 
                                           "p_REL_OD_NI", 
                                           "p_ABS_OD_NI", 
                                           "p_BUP_OD_INJ", 
                                           "p_MET_OD_INJ", 
                                           "p_REL_OD_INJ", 
                                           "p_ABS_OD_INJ"),
                         v_lb = c(p_BUP_OD_NI_lb = l_params_all$p_BUP_OD_NI_lb, 
                                  p_MET_OD_NI_lb = l_params_all$p_MET_OD_NI_lb, 
                                  p_REL_OD_NI_lb = l_params_all$p_REL_OD_NI_lb, 
                                  p_ABS_OD_NI_lb = l_params_all$p_ABS_OD_NI_lb, 
                                  p_BUP_OD_INJ_lb = l_params_all$p_BUP_OD_INJ_lb, 
                                  p_MET_OD_INJ_lb = l_params_all$p_MET_OD_INJ_lb, 
                                  p_REL_OD_INJ_lb = l_params_all$p_REL_OD_INJ_lb, 
                                  p_ABS_OD_INJ_lb = l_params_all$p_ABS_OD_INJ_lb), # lower bound estimate for each param
                         v_ub = c(p_BUP_OD_NI_ub = l_params_all$p_BUP_OD_NI_ub, 
                                  p_MET_OD_NI_ub = l_params_all$p_MET_OD_NI_ub, 
                                  p_REL_OD_NI_ub = l_params_all$p_REL_OD_NI_ub, 
                                  p_ABS_OD_NI_ub = l_params_all$p_ABS_OD_NI_ub, 
                                  p_BUP_OD_INJ_ub = l_params_all$p_BUP_OD_INJ_ub, 
                                  p_MET_OD_INJ_ub = l_params_all$p_MET_OD_INJ_ub, 
                                  p_REL_OD_INJ_ub = l_params_all$p_REL_OD_INJ_ub, 
                                  p_ABS_OD_INJ_ub = l_params_all$p_ABS_OD_INJ_ub)){ # higher bound estimate for each param
  
  
  #v_ub = c(p_S1S2 = 0.50, hr_S1 = 4.5, hr_S2 = 15)){
  n_param <- length(v_param_names)
  m_lhs_unit   <- lhs::randomLHS(n = n_samp, k = n_param) # random latin hypercube sampling
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  for (i in 1:n_param){ # draw parameters
    m_param_samp[, i] <- qunif(m_lhs_unit[,i],
                               min = v_lb[i],
                               max = v_ub[i])
    # ALTERNATIVE prior using beta (or other) distributions
    # m_param_samp[, i] <- qbeta(m_lhs_unit[,i],
    #                            min = 1,
    #                            max = 1)
  }
  return(as.matrix(m_param_samp))
}

X_all_2 <- sample.prior(100)
X_all_3 <- sample.prior(100)

Sig2_global2 <- cov(X_all_2)
Sig2_global3 <- cov(X_all_3)

parscale2 = sqrt(diag(Sig2_global2))
parscale3 = sqrt(diag(Sig2_global3))



which_var <- 
Weights <- 

Sig2 = cov.wt(X_all[which_var,], wt = Weights[which_var]+1/length(Weights), cor = FALSE, center = X_imp, method = "unbias")$cov



log_lik <- function(v_params){ # User defined
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  
  
  source("R/OPTIMA_00_input_parameter_functions.R")
  source("R/OPTIMA_01_model_setup_functions.R")
  source("R/OPTIMA_02_calibration_functions.R")
  
  v_params <- sample.prior(5)
  
  n_samp <- nrow(v_params)
  v_target_names <- c("Fatal Overdoses", "Overdoses")
  n_target       <- length(v_target_names)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_target) 
  colnames(v_llik) <- v_target_names
  v_llik_overall <- numeric(n_samp)
  for(j in 1:n_samp) { # j=1
    jj <- tryCatch( { 
      
      ###   Run model for parameter set "v_params" ###
      l_model_res <- calibration_out(v_params_calib = v_params[j, ], 
                                     l_params_all = l_params_all)
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      ## TARGET 1: Fatal overdoses ("fatal_overdose")
      ## Normal log-likelihood  
      v_llik[j, 1] <- sum(dnorm(x = l_cali_targets$ODF$pe,
                                               mean = l_model_res$fatal_overdose,
                                               sd = l_cali_targets$ODF$se,
                                               log = T))
      
      ## TARGET 2: Non-fatal overdoses ("overdose")
      ## Normal log-likelihood
      v_llik[j, 2] <- sum(dnorm(x = l_cali_targets$ODN$pe,
                                         mean = l_model_res$overdose,
                                         sd = l_cali_targets$ODN$se,
                                         log = T))
      
      
      ## OVERALL
      ## can give different targets different weights
      # To-do: Confirm this calculation
      v_weights <- rep(1, n_target) # currently weight fatal overdoses 1:1 to overall overdoses
      
      ## weighted sum
      v_llik_overall[j] <- v_llik[j, ] %*% v_weights#})
    }, error = function(e) NA) 
    if(is.na(jj)) { v_llik_overall <- -Inf }
  } ## End loop over sampled parameter sets

test1 <- v_llik
test2 <- v_llik_overall




v_llik[1, "Surv"] <- 1
v_llik[1, "Prev"] <- 2
v_llik[1, "PropSick"] <- 3

v_weights <- rep(1, n_target)
## weighted sum
v_llik_overall[j] <- v_llik[j, ] %*% v_weights


test <- log_lik(v_params_calib)
test2 <- log_prior(v_params_calib)
test2_menzies <- log_prior(v_params_calib)

test3 <- prior(test2)



X_all_2 <- sample.prior(5)
l_prior <- log_prior(X_all_2[5,])
llike <- likelihood(X_all_2[5,])
llike2 <- likelihood(X_all_2[6,])
llike3 <- log_lik(X_all_2[6,])

l_model_res <- calibration_out(v_params_calib = X_all_2[6,], 
                               l_params_all = l_params_all)

like_menzies <- likelihood(sample.prior(10))

like_darth <- likelihood(sample.prior(5))


X_all_2[]
#Sig2_global2 <- cov(X_all_2)

log_lik(X_all_2)
l_likelihood(X_all_2)

check_cali1 <- calibration_out(v_params_calib = X_all_2[1, ], 
                              l_params_all = l_params_all)
check_cali2 <- calibration_out(v_params_calib = X_all_2[2, ], 
                               l_params_all = l_params_all)
check_cali3 <- calibration_out(v_params_calib = X_all_2[3, ], 
                               l_params_all = l_params_all)

check_fit1 <- sum(dnorm(x = l_cali_targets$ODF$pe,
                  mean = check_cali1$fatal_overdose,
                  sd = l_cali_targets$ODF$se,
                  log = T))
check_fit2 <- sum(dnorm(x = l_cali_targets$ODF$pe,
                        mean = check_cali2$fatal_overdose,
                        sd = l_cali_targets$ODF$se,
                        log = T))
check_fit3 <- sum(dnorm(x = l_cali_targets$ODF$pe,
                        mean = check_cali3$fatal_overdose,
                        sd = l_cali_targets$ODF$se,
                        log = T))

check_fit4 <- sum(dnorm(x = l_cali_targets$ODN$pe,
                        mean = check_cali1$overdose,
                        sd = l_cali_targets$ODN$se,
                        log = T))
check_fit5 <- sum(dnorm(x = l_cali_targets$ODN$pe,
                        mean = check_cali2$overdose,
                        sd = l_cali_targets$ODN$se,
                        log = T))
check_fit6 <- sum(dnorm(x = l_cali_targets$ODN$pe,
                        mean = check_cali3$overdose,
                        sd = l_cali_targets$ODN$se,
                        log = T))

prior_all = c(prior_all, prior(X_all_2))		# Calculate the prior densities
like_all = c(like_all, likelihood(X_all_2))

X_all <- sample.prior(5)
n_samp <- nrow(X_all)

v_target_names <- c("Fatal Overdoses", "Overdoses")
n_target       <- length(v_target_names)
v_llik <- matrix(0, nrow = n_samp, ncol = n_target) 
colnames(v_llik) <- v_target_names
v_llik_overall <- numeric(n_samp)

#for(j in 1:n_samp) { # j=1
  #jj <- tryCatch( { 
    
    ###   Run model for parameter set "v_params" ###
    l_model_res_1 <- calibration_out(v_params_calib = X_all[1, ], 
                                   l_params_all = l_params_all)
    l_model_res_2 <- calibration_out(v_params_calib = X_all[2, ], 
                                     l_params_all = l_params_all)
    l_model_res_3 <- calibration_out(v_params_calib = X_all[3, ], 
                                     l_params_all = l_params_all)
    l_model_res_4 <- calibration_out(v_params_calib = X_all[4, ], 
                                     l_params_all = l_params_all)
    l_model_res_5 <- calibration_out(v_params_calib = X_all[5, ], 
                                     l_params_all = l_params_all)
    
    ###  Calculate log-likelihood of model outputs to targets  ###
    ## TARGET 1: Fatal overdoses ("fatal_overdose")
    ## Normal log-likelihood  
    #v_llik[j, "Fatal Overdoses"] <- sum(dnorm(x = l_cali_targets$ODF$pe,
    #                                         mean = l_model_res$fatal_overdose,
    #                                         sd = l_cali_targets$ODF$se,
    #                                         log = T))
    
    v_llik[1, "Fatal Overdoses"] <- sum(dnorm(x = l_cali_targets$ODF$pe,
                                              mean = l_model_res_1$fatal_overdose,
                                              sd = l_cali_targets$ODF$se,
                                              log = T))
    v_llik[2, "Fatal Overdoses"] <- sum(dnorm(x = l_cali_targets$ODF$pe,
                                              mean = l_model_res_2$fatal_overdose,
                                              sd = l_cali_targets$ODF$se,
                                              log = T))
    v_llik[3, "Fatal Overdoses"] <- sum(dnorm(x = l_cali_targets$ODF$pe,
                                              mean = l_model_res_3$fatal_overdose,
                                              sd = l_cali_targets$ODF$se,
                                              log = T))
    v_llik[4, "Fatal Overdoses"] <- sum(dnorm(x = l_cali_targets$ODF$pe,
                                              mean = l_model_res_4$fatal_overdose,
                                              sd = l_cali_targets$ODF$se,
                                              log = T))
    v_llik[5, "Fatal Overdoses"] <- sum(dnorm(x = l_cali_targets$ODF$pe,
                                              mean = l_model_res_5$fatal_overdose,
                                              sd = l_cali_targets$ODF$se,
                                              log = T))
    
    ## TARGET 2: Non-fatal overdoses ("overdose")
    ## Normal log-likelihood
    #v_llik[j, "Overdoses"] <- sum(dnorm(x = l_cali_targets$ODN$pe,
    #                                   mean = l_model_res$overdose,
    #                                   sd = l_cali_targets$ODN$se,
    #                                   log = T))
    
    try_that <- v_llik[1, "Overdoses"] <- sum(dnorm(x = as.numeric(l_cali_targets$ODN$pe),
                                        mean = as.numeric(l_model_res_1$overdose),
                                        sd = as.numeric(l_cali_targets$ODN$se),
                                        log = TRUE))
    #try_this <-               sum(dnorm(x = c(0.16, 0.14, 0.12), mean = c(0.175, 0.153, 0.116), sd = c(0.05, 0.05, 0.05), log = TRUE))
    #try_this2 <-              sum(dnorm(x = c(0.00005, 0.0002, 0.0003), mean = c(0.00407, 0.0071, 0.00934), sd = c(1, 1, 1), log = TRUE)) 
    
    
    #x_p1 <- c(0.16, 0.14, 0.12)
    #x_p2 <- c(l_cali_targets$ODN$pe)
    #x_p3 <- l_cali_targets$ODN$pe
    #mean_p1 <- c(0.175, 0.153, 0.116)
    #mean_p2 <- as.numeric(l_model_res_1$overdose)
    #sd_p1 <- c(0.005, 0.005, 0.005)
    #sd_p2 <- as.numeric(l_cali_targets$ODN$se)
    
    v_llik[2, "Overdoses"] <- sum(dnorm(x = l_cali_targets$ODN$pe,
                                        mean = l_model_res_2$overdose,
                                        sd = l_cali_targets$ODN$se,
                                        log = T))
    v_llik[3, "Overdoses"] <- sum(dnorm(x = l_cali_targets$ODN$pe,
                                        mean = l_model_res_3$overdose,
                                        sd = l_cali_targets$ODN$se,
                                        log = T))
    v_llik[4, "Overdoses"] <- sum(dnorm(x = l_cali_targets$ODN$pe,
                                        mean = l_model_res_4$overdose,
                                        sd = l_cali_targets$ODN$se,
                                        log = T))
    v_llik[5, "Overdoses"] <- sum(dnorm(x = l_cali_targets$ODN$pe,
                                        mean = l_model_res_5$overdose,
                                        sd = l_cali_targets$ODN$se,
                                        log = T))
    
    
    ## OVERALL
    ## can give different targets different weights
    # To-do: Confirm this calculation
    v_weights <- rep(1, n_target) # currently weight fatal overdoses 1:1 to overall overdoses
    
    ## weighted sum
    v_llik_overall[1] <- v_llik[1, ] %*% v_weights
    v_llik_overall[2] <- v_llik[2, ] %*% v_weights
    v_llik_overall[3] <- v_llik[3, ] %*% v_weights
    v_llik_overall[4] <- v_llik[4, ] %*% v_weights
    v_llik_overall[5] <- v_llik[5, ] %*% v_weights

    v_like_test <- exp(v_llik_overall) 
    
    plz_work <- log_lik(sample.prior(5))
    
    #v_llik_overall[j] <- v_llik[j]
  #}, error = function(e) NA) 
  #if(is.na(jj)) { v_llik_overall <- -Inf }
#} ## End loop over sampled parameter sets
    
    X_all <- sample.prior(5)
    plz_work <- log_lik(sample.prior(5))

    
    prior_mez <- sample.prior.lhs(5)
    llik <- rep(0,nrow(prior_mez))

    res_j_1 <- mod(c(as.numeric(prior_mez[1,]),1)) 
    
    butt <- llik[1] <- llik[1] + sum(dbinom(c(25,75,50),500,            res_j_1[["prev"]], log=TRUE)) # prevalence likelihood
    dick <- llik[1] <- llik[1] + dnorm(10,              res_j_1[["surv"]],2/1.96, log=TRUE)             # survival likelihood
    poop <- llik[1] <- llik[1] + dnorm(75000,           res_j_1[["tx"]],  5000/1.96, log=TRUE)         # treatment volume likelihood    
    llik
    
    # llik = -41 approx






p_OD <- function(rate = rate,
                 multiplier = multiplier,
                 first_month = FALSE,
                 fatal = FALSE){
  
  # Probability of successful naloxone use
  p_NX_rev <- (p_witness * p_NX_used * p_NX_success)
  
  # Probability of mortality from overdose accounting for baseline overdose fatality and effectiveness of naloxone
  # Subsets overdose into fatal and non-fatal, conditional on different parameters
  p_fatal_OD_NX <- p_fatal_OD * (1 - p_NX_rev)
  
  #convert rates to probabilities - multiply rates by first month multiplier before converting
  # from_state = monthly overdose rate in starting state
  if (first_month){
    p_base_OD <- 1 - exp(-(rate * multiplier)) # check calculation: monthly rate * month should cancel out as long as rates are monthly
    else{
      p_base_OD <- 1 - exp(-rate)
    }
  }
  
  if (fatal){
    p_OD <- ((p_base_OD * (1 - p_fent_exp)) + (p_fent_OD * (p_fent_exp))) * p_fatal_OD_NX
  } else{
    p_OD <- ((p_base_OD * (1 - p_fent_exp)) + (p_fent_OD * (p_fent_exp))) * (1 - p_fatal_OD_NX)
  }
  return(p_OD)
}

# Module to calculate probability of overdose from states
# Probability of overdose
# Non-injection
p_BUP_ODN_NI  <- p_OD(rate = n_BUP_OD_NI, first_month = FALSE, fatal = FALSE)
p_MET_ODN_NI  <- p_OD(rate = n_MET_OD_NI, first_month = FALSE, fatal = FALSE)
p_REL_ODN_NI  <- p_OD(rate = n_REL_OD_NI, first_month = FALSE, fatal = FALSE)
p_ABS_ODN_NI  <- p_OD(rate = n_ABS_OD_NI, first_month = FALSE, fatal = FALSE)
p_BUP_ODF_NI  <- p_OD(rate = n_BUP_OD_NI, first_month = FALSE, fatal = TRUE)
p_MET_ODF_NI  <- p_OD(rate = n_MET_OD_NI, first_month = FALSE, fatal = TRUE)
p_REL_ODF_NI  <- p_OD(rate = n_REL_OD_NI, first_month = FALSE, fatal = TRUE)
p_ABS_ODF_NI  <- p_OD(rate = n_ABS_OD_NI, first_month = FALSE, fatal = TRUE)

# Injection
p_BUP_ODN_INJ <- p_OD(rate = n_BUP_OD_INJ, first_month = FALSE, fatal = FALSE)
p_MET_ODN_INJ <- p_OD(rate = n_MET_OD_INJ, first_month = FALSE, fatal = FALSE)
p_REL_ODN_INJ <- p_OD(rate = n_REL_OD_INJ, first_month = FALSE, fatal = FALSE)
p_ABS_ODN_INJ <- p_OD(rate = n_ABS_OD_INJ, first_month = FALSE, fatal = FALSE)
p_BUP_ODF_INJ <- p_OD(rate = n_BUP_OD_INJ, first_month = FALSE, fatal = TRUE)
p_MET_ODF_INJ <- p_OD(rate = n_MET_OD_INJ, first_month = FALSE, fatal = TRUE)
p_REL_ODF_INJ <- p_OD(rate = n_REL_OD_INJ, first_month = FALSE, fatal = TRUE)
p_ABS_ODF_INJ <- p_OD(rate = n_ABS_OD_INJ, first_month = FALSE, fatal = TRUE)

# Probability of overdose (first month multiplier)
# Non-injection
p_BUP_ODN_NI_4wk  <- p_OD(rate = n_BUP_OD_NI, multiplier = n_BUP_OD_mult, first_month = TRUE, fatal = FALSE)
p_MET_ODN_NI_4wk  <- p_OD(rate = n_MET_OD_NI, multiplier = n_MET_OD_mult, first_month = TRUE, fatal = FALSE)
p_REL_ODN_NI_4wk  <- p_OD(rate = n_REL_OD_NI, multiplier = n_REL_OD_mult, first_month = TRUE, fatal = FALSE)
p_ABS_ODN_NI_4wk  <- p_OD(rate = n_ABS_OD_NI, multiplier = n_ABS_OD_mult, first_month = TRUE, fatal = FALSE)
p_BUP_ODF_NI_4wk  <- p_OD(rate = n_BUP_OD_NI, multiplier = n_BUP_OD_mult, first_month = TRUE, fatal = TRUE)
p_MET_ODF_NI_4wk  <- p_OD(rate = n_MET_OD_NI, multiplier = n_MET_OD_mult, first_month = TRUE, fatal = TRUE)
p_REL_ODF_NI_4wk  <- p_OD(rate = n_REL_OD_NI, multiplier = n_REL_OD_mult, first_month = TRUE, fatal = TRUE)
p_ABS_ODF_NI_4wk  <- p_OD(rate = n_ABS_OD_NI, multiplier = n_ABS_OD_mult, first_month = TRUE, fatal = TRUE)

# Injection
p_BUP_ODN_INJ_4wk <- p_OD(rate = n_BUP_OD_INJ, multiplier = n_BUP_OD_mult, first_month = TRUE, fatal = FALSE)
p_MET_ODN_INJ_4wk <- p_OD(rate = n_MET_OD_INJ, multiplier = n_MET_OD_mult, first_month = TRUE, fatal = FALSE)
p_REL_ODN_INJ_4wk <- p_OD(rate = n_REL_OD_INJ, multiplier = n_REL_OD_mult, first_month = TRUE, fatal = FALSE)
p_ABS_ODN_INJ_4wk <- p_OD(rate = n_ABS_OD_INJ, multiplier = n_ABS_OD_mult, first_month = TRUE, fatal = FALSE)
p_BUP_ODF_INJ_4wk <- p_OD(rate = n_BUP_OD_INJ, multiplier = n_BUP_OD_mult, first_month = TRUE, fatal = TRUE)
p_MET_ODF_INJ_4wk <- p_OD(rate = n_MET_OD_INJ, multiplier = n_MET_OD_mult, first_month = TRUE, fatal = TRUE)
p_REL_ODF_INJ_4wk <- p_OD(rate = n_REL_OD_INJ, multiplier = n_REL_OD_mult, first_month = TRUE, fatal = TRUE)
p_ABS_ODF_INJ_4wk <- p_OD(rate = n_ABS_OD_INJ, multiplier = n_ABS_OD_mult, first_month = TRUE, fatal = TRUE)

#### Overdose ####
# Includes additional calibration related parameters

# Non-injection
# Mean
n_BUP_OD_NI = df_overdose["pe", "BUP_OD_NI"]
n_MET_OD_NI = df_overdose["pe", "MET_OD_NI"]
n_ABS_OD_NI = df_overdose["pe", "ABS_OD_NI"]
n_REL_OD_NI = df_overdose["pe", "REL_OD_NI"]
# Lower bound
n_BUP_OD_NI_lb = df_overdose["low", "BUP_OD_NI"] 
n_MET_OD_NI_lb = df_overdose["low", "MET_OD_NI"]
n_REL_OD_NI_lb = df_overdose["low", "REL_OD_NI"]
n_ABS_OD_NI_lb = df_overdose["low", "ABS_OD_NI"]
# Upper bound
n_BUP_OD_NI_ub = df_overdose["high", "BUP_OD_NI"]
n_MET_OD_NI_ub = df_overdose["high", "MET_OD_NI"]
n_REL_OD_NI_ub = df_overdose["high", "REL_OD_NI"]
n_ABS_OD_NI_ub = df_overdose["high", "ABS_OD_NI"]

# Injection
# Mean
n_BUP_OD_INJ = df_overdose["pe", "BUP_OD_INJ"]
n_MET_OD_INJ = df_overdose["pe", "MET_OD_INJ"]
n_ABS_OD_INJ = df_overdose["pe", "ABS_OD_INJ"]
n_REL_OD_INJ = df_overdose["pe", "REL_OD_INJ"]
# Lower bound
n_BUP_OD_INJ_lb = df_overdose["low", "BUP_OD_INJ"]
n_MET_OD_INJ_lb = df_overdose["low", "MET_OD_INJ"]
n_REL_OD_INJ_lb = df_overdose["low", "REL_OD_INJ"]
n_ABS_OD_INJ_lb = df_overdose["low", "ABS_OD_INJ"]
# Upper bound
n_BUP_OD_INJ_ub = df_overdose["high", "BUP_OD_INJ"]
n_MET_OD_INJ_ub = df_overdose["high", "MET_OD_INJ"]
n_REL_OD_INJ_ub = df_overdose["high", "REL_OD_INJ"]
n_ABS_OD_INJ_ub = df_overdose["high", "ABS_OD_INJ"]

# Overdose transition multipliers for first 4 weeks of treatment and relapse
n_BUP_OD_mult = df_overdose["pe", "BUP_OD_mult"]
n_MET_OD_mult = df_overdose["pe", "MET_OD_mult"]
n_REL_OD_mult = df_overdose["pe", "REL_OD_mult"]
n_ABS_OD_mult = df_overdose["pe", "ABS_OD_mult"]


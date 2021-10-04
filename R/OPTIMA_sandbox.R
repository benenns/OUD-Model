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

# PSA Test Code
df_UP <- read.csv("C:/Users/Ben/Documents/GitHub/OUD-Model/Data/unconditional.csv", row.names = 1, header = TRUE) # Initial parameter values
m_dirichlet_UP = df_UP * 272 # weight matrix by sample size
v_dirichlet_UP_BUP = m_dirichlet_UP["BUP_NI",]
m_BUP_state = rdirichlet(50, c(v_dirichlet_UP_BUP["BUP_NI", "BUPC_NI"], v_dirichlet_UP_BUP["BUP_NI", "MET_NI"], v_dirichlet_UP_BUP["BUP_NI", "METC_NI"], v_dirichlet_UP_BUP["BUP_NI", "ABS_NI"], v_dirichlet_UP_BUP["BUP_NI", "REL_NI"]))
p_BUP_BUPC_NI = m_BUP_state[,1]

## CALIBRATION TEST ##
source("R/OPTIMA_00_input_parameter_functions.R")
source("R/OPTIMA_01_model_setup_functions.R")
source("R/OPTIMA_02_calibration_functions.R")

l_params_all <- load_all_params(file.init = "data/init_params.csv",
                                file.init_dist = "data/init_dist.csv", # calibrate on full trial sample (x% bup; x% met)
                                file.mort = "data/all_cause_mortality.csv",
                                file.death_hr = "data/death_hr.csv",
                                file.frailty = "data/frailty.csv",
                                file.weibull_scale = "data/weibull_scale.csv",
                                file.weibull_shape = "data/weibull_shape.csv",
                                file.unconditional = "data/unconditional.csv",
                                file.overdose = "data/overdose.csv", # includes calibration-related parameters
                                file.hiv = "data/hiv_sero.csv",
                                file.hcv = "data/hcv_sero.csv",
                                file.costs = "data/costs.csv",
                                file.crime_costs = "data/crime_costs.csv",
                                file.qalys = "data/qalys.csv")

l_cali_targets <- list(ODF = read.csv(file = "data/cali_target_odf.csv", header = TRUE),
                       ODN = read.csv(file = "data/cali_target_odn.csv", header = TRUE))

v_params_calib <- c(n_BUP_OD = l_params_all$n_BUP_OD,
                    n_BUPC_OD = l_params_all$n_BUPC_OD,
                    n_MET_OD = l_params_all$n_MET_OD,
                    n_METC_OD = l_params_all$n_METC_OD,
                    n_REL_OD = l_params_all$n_REL_OD,
                    n_ABS_OD = l_params_all$n_ABS_OD,
                    n_TX_OD_mult = l_params_all$n_TX_OD_mult,
                    n_TXC_OD_mult = l_params_all$n_TXC_OD_mult,
                    n_REL_OD_mult = l_params_all$n_REL_OD_mult,
                    n_INJ_OD_mult = l_params_all$n_INJ_OD_mult,
                    n_fatal_OD = l_params_all$n_fatal_OD)

#test <- calibration_out(v_params_calib = v_params_calib, l_params_all = l_params_all)

l_model_res <- calibration_out(v_params_calib = v_params_calib, 
                               l_params_all = l_params_all)

x <- l_cali_targets$ODF$pe
mean_x <- l_model_res$fatal_overdose
sd_x <- l_cali_targets$ODF$se

y <- l_cali_targets$ODN$pe
mean_y <- l_model_res$overdose
sd_y <- l_cali_targets$ODN$se


v_params <- c(1,2,3,4)
v_lb <- c(0.5,1.5,2.5,3.5)
v_ub <- c(1.5,2.5,3.5,4.5)

n_param <- length(v_params)

n_samp <- 10
lprior <- rep(0, n_samp)
for (i in 1:n_param){
  lprior <- lprior + dunif(v_params[, i],
                           min = v_lb[i],
                           max = v_ub[i], 
                           log = T)
}

g <- dgamma(4, shape = 2, scale = 0.5, log = T)
f <- dgamma(4, shape = 2, scale = 0.5, log = F)


m_lhs_unit   <- lhs::randomLHS(n = 100, k = 9)

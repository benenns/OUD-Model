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

#############
#TEST MODULE#
#############
dat <- as.matrix(read.csv("C:/Users/Benjamin/Desktop/Book1.csv", header = FALSE))

a_test <- array(0, dim = c(5, 5, 10),
                dimnames = list())
for (i in 1:10){
  a_test[, , i] <- dat[, ]
}


a_test_trace <- a_test_trace_death <- array(0, dim = c(11, 5, 11))
m_test <- array(0.99, dim = c(5, 10))
m_test
#m_test_death <- 1 - m_test
#m_death <- array (0, dim = c(11, 1))
v_init <- c(1, 0, 0, 0, 0)
a_test_trace[1, , 1] <- a_test_trace_death [1, , 1] <- v_init

#m_soj <- a_test[ , , 1] * m_test[, 2]
#m_soj
#
#current <- as.vector(a_test_trace[1, , 1])
#current
#
#same <- as.vector(current * diag(m_soj))
#same
#
#a_test_trace[2, , 2] <- same
#a_test_trace
#
#diag(m_soj) <- 0
#
#new <- as.vector(current %*% m_soj)
#new

# All model time periods
for(i in 2:10){
  for(j in 1:(i - 1)){
    m_sojourn <- a_test[, , j] * m_test[, i]
    #m_sojourn_death <- a_test[, , j] * m_test_death[, i]
    
    v_current_state <- as.vector(a_test_trace[i - 1, , j])
    
    v_same_state <- as.vector(v_current_state * diag(m_sojourn))
    #v_same_state_death <- as.vector(v_current_state * diag(m_sojourn_death))
    
    
    a_test_trace[i, ,j + 1] <- v_same_state 
    #a_test_trace_death[i, ,j + 1] <- v_same_state_death 
    
    diag(m_sojourn) <- 0
    #diag(m_sojourn_death) <- 0
    
    v_new_state <- as.vector(v_current_state %*% m_sojourn)
    #v_new_state_death <- as.vector(v_current_state %*% m_sojourn_death)
    
    a_test_trace[i,,1] <- v_new_state + a_test_trace[i,,1]
    #a_test_trace_death[i,,1] <- v_new_state_death + a_test_trace_death[i,,1]
    
    #m_death[, i] <- rowSums(a_test_trace_death[i, , ])
  }
}


g <- array(0, dim = c(11, 5))
for (i in 1:11){
g[i, ] <- rowSums(a_test_trace[i, ,])
}

deaths <- array(0, dim = c(11, 1))
for (i in 1:11){
deaths[i,] <- as.vector(1 - rowSums(g[i,]))
}
deaths

a_test
a_test_trace
m_death
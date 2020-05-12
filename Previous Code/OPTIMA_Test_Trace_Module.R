#############
#TEST MODULE#
#############
dat <- as.matrix(read.csv("C:/Users/Benjamin/Desktop/Book1.csv", header = FALSE))

a_test <- array(0, dim = c(5, 5, 10),
                dimnames = list())
for (i in 1:10){
  a_test[, , i] <- dat[, ]
}


a_test_trace <- array(0, dim = c(11, 5, 11))
m_test <- array(0.99, dim = c(5, 10))

m_test_death <- 1 - m_test

v_init <- c(1, 0, 0, 0, 0)
a_test_trace[1, , 1] <- a_test_trace_d[1, , 1] <- v_init

# All model time periods
for(i in 2:10){
  for(j in 1:(i - 1)){
    m_sojourn <- a_test[, , j] * m_test[, i - 1]
    m_sojourn_d <- a_test[, , j] * m_test_death[, i - 1]
    
    v_current_state <- as.vector(a_test_trace[i - 1, , j])
    #v_current_state_d <- as.vector(a_test_trace_d[i - 1, , j])
    
    v_same_state <- as.vector(v_current_state * diag(m_sojourn))
    v_same_state_d <- as.vector(v_current_state * diag(m_sojourn_d))

    a_test_trace[i, ,j + 1] <- v_same_state 
    a_test_trace_d[i, ,j + 1] <- v_same_state_d 

    
    diag(m_sojourn) <- 0
    diag(m_sojourn_d) <- 0
    
    v_new_state <- as.vector(v_current_state %*% m_sojourn)
    v_new_state_d <- as.vector(v_current_state %*% m_sojourn_d)

    
    a_test_trace[i,,1] <- v_new_state + a_test_trace[i,,1]
    a_test_trace_d[i,,1] <- v_new_state_d + a_test_trace_d[i,,1]
  }
}
m_M_trace <- array(0, dim = c(11, 5))
m_M_trace_d <- array(0, dim = c(11, 5))
for (i in 1:10){
  m_M_trace[i, ] <- rowSums(a_test_trace[i, ,])
  m_M_trace_d[i, ] <- rowSums(a_test_trace_d[i, ,])
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


###################################
##### CREATE TRANSITION ARRAY #####
###################################
n_states_test <- 10
n_t_test <- 10

a_TDP_test <- array(0, dim = c(n_states_test, n_states_test, n_t_test))
m_TDP_test <- array(0.9, dim = c(n_states_test, n_t_test))
m_UP_test <- array(0.111111, dim = c(n_states_test, n_states_test))
m_leave_test <- 1 - m_TDP_test

for (i in 1:10){
  a_TDP_test[, , i] <- m_UP_test * m_leave_test[, i]
}
a_TDP_test

for (i in 1:10){
  for (j in 1:10){
    a_TDP_test[j, j, i] <- m_TDP_test[j, i]
  } 
}
a_TDP_test
rowSums(a_TDP_test)


# HIV
p_sero_BUP1_NI <- 0.005
a_TDP[, 1, ] <- a_TDP[, 1, ] * (1 - p_sero_BUP1_NI)
a_TDP[, 2, ] <- a_TDP[, 2, ] * (1 - p_sero_BUP1_NI)
a_TDP[, 3, ] <- a_TDP[, 3, ] * (1 - p_sero_BUP1_NI)
a_TDP[, 4, ] <- a_TDP[, 4, ] * (1 - p_sero_BUP1_NI)
a_TDP[, 5, ] <- a_TDP[, 5, ] * (1 - p_sero_BUP1_NI)

a_TDP[, 6, ] <- a_TDP[, 6, ] * p_sero_BUP1_NI
a_TDP[, 7, ] <- a_TDP[, 7, ] * p_sero_BUP1_NI
a_TDP[, 8, ] <- a_TDP[, 8, ] * p_sero_BUP1_NI
a_TDP[, 9, ] <- a_TDP[, 9, ] * p_sero_BUP1_NI
a_TDP[, 10, ] <- a_TDP[, 10, ] * p_sero_BUP1_NI

a_TDP


m_check <- array(0, dim = c(n_states, n_t))

for (i in 1:n_t){
m_check[, i] <- rowSums(a_TDP[, , i])
}

a <- array(0, dim = c(n_states, n_t))

a[, 1] <- rowSums(a_TDP[1 & 5, ,1])
a

a <- apply(a_TDP, 1, sum)
a


### Function to check transition array
valid <- (apply(a_TDP, 3, function(x) sum(rowSums(x))) == n_states)
if (!isTRUE(all_equal(as.numeric(sum(valid)), as.numeric(n_t)))) {
  if(err_stop) {
    stop("This is not a valid transition Matrix")
  }
  
  if(verbose){
    warning("This is not a valid transition Matrix")
  } 







p_frailty_BUP_NI_2 <- 0.725
p_weibull_shape_BUP_NI <- 0.613
p_weibull_scale_BUP_NI <- 0.153
#t <- 1
other <- function(t){
temp <- p_weibull_scale_BUP_NI * (((t - 1)^p_weibull_shape_BUP_NI) - (t^p_weibull_shape_BUP_NI))
return(temp)
}
test <- vector()
test2 <- vector()

for (t in 1:10){
test[t] <- as.vector(exp(p_frailty_BUP_NI_2 * p_weibull_scale_BUP_NI * (((t - 1)^p_weibull_shape_BUP_NI) - (t^p_weibull_shape_BUP_NI))))
test2[t] <- as.vector(exp(p_frailty_BUP_NI_2 * other(t)))
}
test
test2

for (i in 1:10){
  m_leave[3, i] <- test[i]
}


b <- array(0)
b

a <- array(1:16, c(2, 2, 2))

rowSums(a[, , 1])

rowSums(a)
colSums(a)

}
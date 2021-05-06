## Code developed for Menzies et al. "Bayesian methods for calibrating
## health policy models: a tutorial" Pharmacoeconomics.

rm(list = ls()) # to clean the workspace

##############################################
################### MODEL ####################
##############################################
# This function estimates outcomes describing epidemiology of a hypothetical disease
# as well as outcomes (life-years, costs) for estimating cost-effectiveness of a policy 
# to expand treatment access to individual's with early stage disease.


mod <- function(par_vector,project_future=FALSE) {
  # par_vector: a vector of model parameters
  # project_future: TRUE/FALSE, whether to project future outcomes for policy comparison
  pop_size   <- 1e6             # population size hard-coded as 1 million
  mu_b       <- 0.015           # background mortality rate hard-coded as 0.015
  mu_e       <- par_vector[1]   # cause-specific mortality rate with early-stage disease
  mu_l       <- par_vector[2]   # cause-specific mortality rate with late-stage disease
  mu_t       <- par_vector[3]   # cause-specific mortality rate on treatment
  p          <- par_vector[4]   # transtion rate from early to late-stage disease
  r_l <- r_e <- par_vector[5]   # rate of uptake onto treatment (r_l = late-stage disease;r_e = early-stage disease)
  rho        <- par_vector[6]   # effective contact rate 
  b          <- par_vector[7]   # fraction of population in at-risk group
  c          <- par_vector[8]   # annual cost of treatment  
  
  ######## Prepare to run model ###################
  n_yrs    <- if(project_future) { 51 } else { 30 }  # no. years to simulate (30 to present, 51 for 20 year analytic horizon)
  sim      <- if(project_future) { 1:2 } else { 1 }  # which scenarios to simulate: 1 = base case, 2 = expanded treatment access
  v_mu     <- c(0,0,mu_e,mu_l,mu_t)+mu_b             # vector of mortality rates
  births   <- pop_size*mu_b*c(1-b,b)                 # calculate birth rate for equilibrium population before epidemic
  init_pop <- pop_size*c(1-b,b-0.001,0.001,0,0,0)    # creates starting vector for population
  trace    <- matrix(NA,12*n_yrs,6)                  # creates a table to store simulation trace
  colnames(trace) <- c("N","S","E","L","T","D")
  results  <- list()                                 # creates a list to store results
  
  ######## Run model ###################
  for(s in sim) {
    P0 <- P1 <- init_pop 
    for(m in 1:(12*n_yrs)) {
      lambda    <- rho*sum(P0[3:4])/sum(P0[2:5]) # calculates force of infection
      P1[1:2]   <- P1[1:2]+births/12             # births
      P1[-6]    <- P1[-6]-P0[-6]*v_mu/12         # deaths: N, S, E, L, T, to D
      P1[6]     <- P1[6]+sum(P0[-6]*v_mu/12)     # deaths:N, S, E, L, T, to D
      P1[2]     <- P1[2]-P0[2]*lambda/12         # infection: S to E
      P1[3]     <- P1[3]+P0[2]*lambda/12         # infection: S to E
      P1[3]     <- P1[3]-P0[3]*p/12              # progression: E to L
      P1[4]     <- P1[4]+P0[3]*p/12              # progression: E to L
      P1[4]     <- P1[4]-P0[4]*r_l/12            # treatment uptake: L to T
      P1[5]     <- P1[5]+P0[4]*r_l/12            # treatment uptake: L to T
      if(s==2 & m>(12*30)) {
        P1[3]   <- P1[3]-P0[3]*r_e/12            # treatment uptake: E to T (scenario 2)
        P1[5]   <- P1[5]+P0[3]*r_e/12            # treatment uptake: E to T (scenario 2)
      }
      trace[m,] <- P0 <- P1                      # fill trace, reset pop vectors
    }
    results[[s]] <- trace                        # save results for each scenario
  }
  
  ######## Report results ###################
  if(project_future==FALSE) {
    ## Return calibration metrics, if project_future = FALSE
    return(list(prev = (rowSums(trace[,3:5])/rowSums(trace[,1:5]))[c(10,20,30)*12],  # Prevalence at 10,20,30 years
                surv = 1/(v_mu[3]+p)+ p/(v_mu[3]+p)*(1/v_mu[4]),                     # HIV survival without treatment
                tx   = trace[30*12,5]                                                # Treatment volume at 30 years
    ) )
  } else {
    ## Policy projections for CE analysis, if project_future = TRUE
    return(list(trace0   = results[[1]],     # Trace without expanded treatment access
                trace1   = results[[2]],     # Trace with expanded treatment access
                inc_LY   = sum(results[[2]][(30*12+1):(51*12),-6]-results[[1]][(30*12+1):(51*12),-6])/12,  # incr. LY lived with expanded tx
                inc_cost = sum(results[[2]][(30*12+1):(51*12),5]-results[[1]][(30*12+1):(51*12),5])*c/12   # incr. cost  with expanded tx                  
    ) )  
  }
}

# Test it 
mod(rep(0.5,8),project_future=F)  # works
mod(rep(0.5,8),project_future=T)  # works


###################################################
############### SAMPLE FROM PRIOR #################
###################################################
# This function returns a sample from the prior parameter distribution,
# with each column a different parameter and each row a different parameter set.  
# sample.prior.srs: this uses a simple random sample of the parameter space.
# sample.prior.lhs: this uses a latin-hypercube sample of the parameter space.

## Simple random sample
sample.prior.srs <- function(n) {
  # n: the number of samples desired
  draws <- data.frame(mu_e  = rlnorm(n,log(0.05)-1/2*0.5^2,0.5),
                      mu_l  = rlnorm(n,log(0.25)-1/2*0.5^2,0.5),
                      mu_t  = rlnorm(n,log(0.025)-1/2*0.5^2,0.5),
                      p     = rlnorm(n,log(0.1)-1/2*0.5^2,0.5),
                      r_l   = rlnorm(n,log(0.5)-1/2*0.5^2,0.5),
                      rho   = rlnorm(n,log(0.5)-1/2*0.5^2,0.5),
                      b     = rbeta(n,2,8),
                      c     = rlnorm(n,log(1000)-1/2*0.2^2,0.2)
  )
  return(as.matrix(draws))
}

## Latin hypercube sample
#  install.packages("lhs")
require(lhs)

sample.prior.lhs <- function(n) {
  # n: the number of samples desired
  draws0 <- randomLHS(n=n,k=8)
  draws  <- data.frame( mu_e  = qlnorm(draws0[,1],log(0.05)-1/2*0.5^2,0.5),
                        mu_l  = qlnorm(draws0[,2],log(0.25)-1/2*0.5^2,0.5),
                        mu_t  = qlnorm(draws0[,3],log(0.025)-1/2*0.5^2,0.5),
                        p     = qlnorm(draws0[,4],log(0.1)-1/2*0.5^2,0.5),
                        r_l   = qlnorm(draws0[,5],log(0.5)-1/2*0.5^2,0.5),
                        rho   = qlnorm(draws0[,6],log(0.5)-1/2*0.5^2,0.5),
                        b     = qbeta(draws0[,7],2,8),
                        c     = qlnorm(draws0[,8],log(1000)-1/2*0.2^2,0.2)
  )
  return(as.matrix(draws))
}

sample.prior <- sample.prior.lhs  # use the lhs version as this is more efficient

sample.prior.2 <- function(n_samp,
                         v_param_names = c("p_S1S2", "hr_S1", "hr_S2"),
                         v_lb = c(p_S1S2 = 0.01, hr_S1 = 1.0, hr_S2 = 5),
                         v_ub = c(p_S1S2 = 0.50, hr_S1 = 4.5, hr_S2 = 15)){
  n_param <- length(v_param_names)
  m_lhs_unit   <- lhs::randomLHS(n = n_samp, k = n_param)
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  for (i in 1:n_param){
    m_param_samp[, i] <- qunif(m_lhs_unit[,i],
                               min = v_lb[i],
                               max = v_ub[i])
  }
  return(m_param_samp)
}

# Test it 
prior <- sample.prior(100)  # works.
prior.2 <- sample.prior.2(100)

cov <- cov(prior)
cov.2 <- cov(prior.2)

diag <- diag(cov)
diag.2 <- diag(cov.2)

length <- length(sqrt(diag))
length.2 <- length(sqrt(diag.2))


###################################################
################### RESULTS UNCALIBRATED ##########
###################################################
# This section uses mod and sample.prior functions to calculate results for the analysis 
# via Monte Carlo simulation, without making use of the calibration data.

# Draw sample from prior
#set.seed(1234)                             # set random seed for reproducibility
#samp <- sample.prior(1e5)                  # 100,000 draws from prior
#
## Generate estimates for inc cost and inc LY via MC simulation (may take a while)
#IncCost <- IncLY <- rep(NA,nrow(samp))     # vectors to collect results
#for(i in 1:nrow(samp)) { # i=1
#  tmp        <- mod(samp[i,],project_future=T)
#  IncLY[i]   <- tmp$inc_LY
#  IncCost[i] <- tmp$inc_cost
#  if(i/100==round(i/100,0)) { 
#    cat('\r',paste(i/nrow(samp)*100,"% done",sep="")) 
#  } 
#}
#
#Uncalib <- list(samp=samp,IncCost=IncCost,IncLY=IncLY)
## save(Uncalib,file="Uncalib_results.rData")
## load("Uncalib_results.rData")
#
#### Results for uncalibrated model
#IncCost <- Uncalib$IncCost
#IncLY   <- Uncalib$IncLY
#
## Distribution of inc costs and inc LY
#hist(IncLY,breaks=100,col="blue3",border=F)
#hist(IncCost,breaks=100,col="red3",border=F)
#
## Summary results
#mean(IncLY/1e3); quantile(IncLY/1e3,c(1,39)/40)                     # 213 (6, 775) thousands
#mean(IncCost/1e6); quantile(IncCost/1e6,c(1,39)/40)                 # 277 (-34, 1235) millions
#mean(IncCost)/mean(IncLY); quantile(IncCost/IncLY,c(1,39)/40)       # 1300 (dominant, 5720)


###################################################
################### LIKELIHOOD ####################
###################################################
# This function calculates the log-likelihood for a parameter set or matrix of parameter sets excluding c,
# the parameter for treatment cost. c is fixed at an arbitrary value (1) as it has no role in the calibration.

l_likelihood <- function(par_vector) {
  # par_vector: a vector (or matrix) of model parameters
  if(is.null(dim(par_vector))) par_vector <- t(par_vector)
  llik <- rep(0,nrow(par_vector))
  for(j in 1:nrow(par_vector)) {
    jj <- tryCatch( {
      
      res_j <- mod(c(as.numeric(par_vector[j,]),1)) 
      
      llik[j] <- llik[j] + sum(dbinom(c(25,75,50),500,            res_j[["prev"]], log=TRUE)) # prevalence likelihood
      llik[j] <- llik[j] + dnorm(10,              res_j[["surv"]],2/1.96, log=TRUE)             # survival likelihood
      llik[j] <- llik[j] + dnorm(75000,           res_j[["tx"]],  5000/1.96, log=TRUE)         # treatment volume likelihood
      
    }, error = function(e) NA)
    if(is.na(jj)) { llik[j] <- -Inf } 
  }
  return(llik)
}

# Test it 
l_likelihood(rbind(rep(0.5,7),rep(0.6,7))) # works


###################################################
##################### PRIOR #######################
###################################################
# This function calculates the log-prior for a parameter set or matrix of parameter sets excluding c,
# the parameter for treatment cost. c is fixed at an arbitrary value (1) as it has no role in the calibration.

l_prior <- function(par_vector) {
  # par_vector: a vector (or matrix) of model parameters (omits c)
  if(is.null(dim(par_vector))) par_vector <- t(par_vector)
  lprior <- rep(0,nrow(par_vector))
  lprior <- lprior+dlnorm(par_vector[,1],log(0.05 )-1/2*0.5^2,0.5,log=TRUE)    # mu_e
  lprior <- lprior+dlnorm(par_vector[,2],log(0.25 )-1/2*0.5^2,0.5,log=TRUE)    # mu_l
  lprior <- lprior+dlnorm(par_vector[,3],log(0.025)-1/2*0.5^2,0.5,log=TRUE)    # mu_t
  lprior <- lprior+dlnorm(par_vector[,4],log(0.1  )-1/2*0.5^2,0.5,log=TRUE)    # p
  lprior <- lprior+dlnorm(par_vector[,5],log(0.5  )-1/2*0.5^2,0.5,log=TRUE)    # r_l
  lprior <- lprior+dlnorm(par_vector[,6],log(0.5  )-1/2*0.5^2,0.5,log=TRUE)    # rho
  lprior <- lprior+dbeta( par_vector[,7],2,8,log=TRUE)                         # b
  return(lprior)
}

# Test it 
l_prior(rbind(rep(0.5,7),rep(0.6,7))) # works


###################################################
################## MAP ESTIMATION #################
###################################################
# This section obtains a single best-fitting parameter set via maximum a posteriori estimation.
# To do so, an optimization algorithm is used to identify the parameter set that maximizes
# the sum of log-prior plus log-likelihood, which is equal to the log posterior (plus a constant which can be ignored).
# Different optimization routines are tried, BFGS works best for this example. 

# Function for log-posterior
l_post <- function(par_vector) {
  return( l_prior(par_vector) + l_likelihood(par_vector) )
}

# Optimize with various methods in optim tool-box
optOut_nm   <- optim(rep(.5,7),l_post              ,control=list(fnscale=-1)); optOut_nm    # est max(l_post) = -43.5
optOut_cg   <- optim(rep(.5,7),l_post,method="CG"  ,control=list(fnscale=-1)); optOut_cg    # est max(l_post) = -12.0
optOut_bfgs <- optim(rep(.5,7),l_post,method="BFGS",control=list(fnscale=-1)); optOut_bfgs  # est max(l_post) = -6.94

# Best fitting parameter set, incl. prior mode for c
mode_c  <- exp(log(1000)-1/2*0.2^2-0.2^2)
par_map <- c(optOut_bfgs$par,mode_c)


###################################################
#################### RESULTS MAP ##################
###################################################
# Results for model with best fitting parameter set identified with MAP

res_map <- mod(par_map,project_future = T)
res_map[["inc_cost"]]                       # inc cost
res_map[["inc_LY"]]                         # inc cost
res_map[["inc_cost"]]/res_map[["inc_LY"]]   # icer


###################################################
###################### SIR ########################
###################################################
# This section obtains a sample from the posterior parameter distribution via SIR.

# First, lets redefine sample.prior to only provide samples for the first seven parameters, excluding c.
sample.prior <- function(n) {
  draws0 <- randomLHS(n=n,k=7)
  draws  <- data.frame( mu_e  = qlnorm(draws0[,1],log(0.05)-1/2*0.5^2,0.5),
                        mu_l  = qlnorm(draws0[,2],log(0.25)-1/2*0.5^2,0.5),
                        mu_t  = qlnorm(draws0[,3],log(0.025)-1/2*0.5^2,0.5),
                        p     = qlnorm(draws0[,4],log(0.1)-1/2*0.5^2,0.5),
                        r_l   = qlnorm(draws0[,5],log(0.5)-1/2*0.5^2,0.5),
                        rho   = qlnorm(draws0[,6],log(0.5)-1/2*0.5^2,0.5),
                        b     = qbeta(draws0[,7],2,8))
  return(as.matrix(draws))
}

# Generate 100,000 samples from prior
n_samp <- 100 #1e5
samp_i <- sample.prior(n_samp)

# Calculate likelihood for each parameter set in samp_i
llik_i <- rep(NA,n_samp)
for(i in 1:n_samp) { 
  llik_i[i] <- l_likelihood(samp_i[i,])
  if(i/100==round(i/100,0)) { 
    cat('\r',paste(i/nrow(samp_i)*100,"% done",sep="")) 
  } 
}

# Calculate weights for resample (i.e. exponentiate the log-likelihood)
# Note: subtracting off the maximum log-likelihood before exponentiating 
# helps avoid numerical under/overflow, which would result in weights of Inf or 0.
wt <- exp(llik_i-max(llik_i)) / sum(exp(llik_i-max(llik_i)))

# Resample from samp_i with wt as sampling weights
id_samp  <- sample.int(n_samp,replace=T,prob=wt)
post_sir <- samp_i[id_samp,]

# Unique parameter sets
length(unique(id_samp)) # 797

# Effective sample size
sum(table(id_samp))^2/sum(table(id_samp)^2) # 88.26

# Max weight
max(table(id_samp))/sum(table(id_samp)) # 0.033

# Could use the parameter sets in post_sir to caliclate results, but lets 
# try a more efficient technique: IMIS.

###################################################
###################### IMIS #######################
###################################################
# This section obtains a sample from the posterior parameter distribution via IMIS
#  install.packages("IMIS")
require(IMIS)

# IMIS needs three functions to be defined: 
#   sample.prior -- draws samples, and we have already created this 

#   prior -- evaluates prior density of a parameter set or sets
prior <- function(par_vector) { 
  exp(l_prior(par_vector))
}

#   likelihood -- evaluates likelihood of a parameter set or sets

likelihood <- function(par_vector) { 
  exp(l_likelihood(par_vector)) 
}

# Run IMIS
set.seed(1234)
imis_res <- IMIS(B=100,B.re=1e4,number_k=400,D=1)  

# Draw samples for parameter c from prior
set.seed(1234)
c_i       <- rlnorm(1e4,log(1000)-1/2*0.2^2,0.2)
post_imis <- cbind(imis_res$resample,c_i)

# Unique parameter sets
length(unique(imis_res$resample[,1])) # 6372

# Effective sample size
sum(table(imis_res$resample[,1]))^2/ sum(table(imis_res$resample[,1])^2) # 4713

# Max weight
max(table(imis_res$resample[,1]))/sum(table(imis_res$resample[,1])) # 8e-04


###################################################
################## RESULTS: IMIS ##################
###################################################

#  Generate inc cost / inc LY using imis_post
IncCostC <- IncLYC <- rep(NA,nrow(post_imis))
for(i in 1:nrow(post_imis)) { # i=1
  tmp <- mod(post_imis[i,],project_future=T)
  IncLYC[i]   <- tmp$inc_LY
  IncCostC[i] <- tmp$inc_cost
  if(i/100==round(i/100,0)) { 
    cat('\r',paste(i/nrow(samp_i)*100,"% done",sep="")) 
  } 
}

Calib <- list(post_imis=post_imis,IncCostC=IncCostC,IncLYC=IncLYC)
# save(Calib,file="Calib_results.rData")
# load("Calib_results.rData")

# Distribution of inc costs and inc LY
hist(IncLYC,breaks=100,col="blue3",border=F)
hist(IncCostC,breaks=100,col="red3",border=F)

# Summary results
mean(IncLYC/1e3); quantile(IncLYC/1e3,c(1,39)/40)                     # 130 (64, 228) thousands
mean(IncCostC/1e6); quantile(IncCostC/1e6,c(1,39)/40)                 # 123 (-4, 312) millions
mean(IncCostC)/mean(IncLYC); quantile(IncCostC/IncLYC,c(1,39)/40)     # 947 (dominant, 2010)



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# The end.  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
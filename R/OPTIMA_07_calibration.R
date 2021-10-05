rm(list = ls()) # to clean the workspace

#### Load packages, data and functions ####
#### Load packages and functions ####
library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(lhs)
library(IMIS)

# To-do: Move functions into R package for OUD model
source("R/OPTIMA_00_input_parameter_functions.R")
source("R/OPTIMA_01_model_setup_functions.R")
source("R/OPTIMA_02_calibration_functions.R")

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
                                file.hiv = "data/hiv_sero.csv",
                                file.hcv = "data/hcv_sero.csv",
                                file.costs = "data/costs.csv",
                                file.crime_costs = "data/crime_costs.csv",
                                file.qalys = "data/qalys.csv")

# Load calibration inputs #
v_cali_param_names <- c("n_TX_OD",
                        "n_TXC_OD",
                        "n_REL_OD", 
                        #"n_ABS_OD",
                        "n_TX_OD_mult",
                        "n_TXC_OD_mult",
                        "n_REL_OD_mult",
                        "n_INJ_OD_mult",
                        "n_fatal_OD")
v_shape <- c(n_TX_OD_shape   = l_params_all$n_TX_OD_shape,
             n_TXC_OD_shape  = l_params_all$n_TXC_OD_shape,
             #n_MET_OD_shape   = l_params_all$n_MET_OD_shape,
             #n_METC_OD_shape  = l_params_all$n_METC_OD_shape,
             n_REL_OD_shape   = l_params_all$n_REL_OD_shape, 
             #n_ABS_OD_shape   = l_params_all$n_ABS_OD_shape, 
             n_TX_OD_mult_shape  = l_params_all$n_TX_OD_mult_shape,
             n_TXC_OD_mult_shape = l_params_all$n_TXC_OD_mult_shape,
             n_REL_OD_mult_shape = l_params_all$n_REL_OD_mult_shape,
             n_INJ_OD_mult_shape = l_params_all$n_INJ_OD_mult_shape,
             n_fatal_OD_shape    = l_params_all$n_fatal_OD_shape) # lower bound estimate for each param
v_scale <- c(n_TX_OD_scale   = l_params_all$n_TX_OD_scale,
             n_TXC_OD_scale  = l_params_all$n_TXC_OD_scale,
             #n_MET_OD_scale   = l_params_all$n_MET_OD_scale,
             #n_METC_OD_scale  = l_params_all$n_METC_OD_scale,
             n_REL_OD_scale   = l_params_all$n_REL_OD_scale, 
             #n_ABS_OD_scale   = l_params_all$n_ABS_OD_scale, 
             n_TX_OD_mult_scale  = l_params_all$n_TX_OD_mult_scale,
             n_TXC_OD_mult_scale = l_params_all$n_TXC_OD_mult_scale,
             n_REL_OD_mult_scale = l_params_all$n_REL_OD_mult_scale,
             n_INJ_OD_mult_scale = l_params_all$n_INJ_OD_mult_scale,
             n_fatal_OD_scale    = l_params_all$n_fatal_OD_scale)

#### Load calibration targets ####
#data("03_calibration_targets")
l_cali_targets <- list(ODF = read.csv(file = "data/cali_target_odf.csv", header = TRUE),
                       ODN = read.csv(file = "data/cali_target_odn.csv", header = TRUE))

#### Visualize targets ####
### TARGET 1: Overdose deaths ("ODF")
plotrix::plotCI(x    = l_cali_targets$ODF$Time, 
                y    = l_cali_targets$ODF$pe, 
                ui   = l_cali_targets$ODF$high,
                li   = l_cali_targets$ODF$low,
                #ylim = c(0, 1), 
                xlab = "Month", ylab = "Fatal Overdoses")

### TARGET 2: Non-fatal overdose ("ODN")
plotrix::plotCI(x    = l_cali_targets$ODN$Time, 
                y    = l_cali_targets$ODN$pe, 
                ui   = l_cali_targets$ODN$high,
                li   = l_cali_targets$ODN$low,
                #ylim = c(0, 1), 
                xlab = "Month", ylab = "Non-fatal Overdoses")


#### Run calibration algorithms ####
# Check that it works
#v_params_calib <- c(p_BUP_OD_NI = l_params_all$p_BUP_OD_NI, 
#              p_MET_OD_NI = l_params_all$p_MET_OD_NI, 
#              p_REL_OD_NI = l_params_all$p_REL_OD_NI, 
#              p_ABS_OD_NI = l_params_all$p_ABS_OD_NI, 
#              p_BUP_OD_INJ = l_params_all$p_BUP_OD_INJ, 
#              p_MET_OD_INJ = l_params_all$p_MET_OD_INJ, 
#              p_REL_OD_INJ = l_params_all$p_REL_OD_INJ, 
#              p_ABS_OD_INJ = l_params_all$p_ABS_OD_INJ)

#test <- calibration_out(v_params_calib = v_params_calib, l_params_all = l_params_all)

#### Specify calibration parameters ####
### Specify seed (for reproducible sequence of random numbers)
set.seed(3730687)

### Number of random samples to obtain from the posterior distribution 
n_resamp <- 1000

### Names and number of input parameters to be calibrated
v_param_names  <- c("n_TX_OD",
                    "n_TXC_OD",
                    #"n_MET_OD",
                    #"n_METC_OD",
                    "n_REL_OD", 
                    #"n_ABS_OD",
                    "n_TX_OD_mult",
                    "n_TXC_OD_mult",
                    "n_REL_OD_mult",
                    "n_INJ_OD_mult",
                    "n_fatal_OD")

#n_param        <- length(v_param_names)

### Vector with range on input search space
v_lb = c(n_TX_OD_lb   = l_params_all$n_TX_OD_lb,
         n_TXC_OD_lb  = l_params_all$n_TXC_OD_lb,
         #n_MET_OD_lb   = l_params_all$n_MET_OD_lb,
         #n_METC_OD_lb  = l_params_all$n_METC_OD_lb,
         n_REL_OD_lb   = l_params_all$n_REL_OD_lb, 
         #n_ABS_OD_lb   = l_params_all$n_ABS_OD_lb, 
         n_TX_OD_mult_lb  = l_params_all$n_TX_OD_mult_lb,
         n_TXC_OD_mult_lb = l_params_all$n_TXC_OD_mult_lb,
         n_REL_OD_mult_lb = l_params_all$n_REL_OD_mult_lb,
         n_INJ_OD_mult_lb = l_params_all$n_INJ_OD_mult_lb,
         n_fatal_OD_lb    = l_params_all$n_fatal_OD_lb) # lower bound estimate for each param

v_ub = c(n_TX_OD_ub   = l_params_all$n_TX_OD_ub,
         n_TXC_OD_ub  = l_params_all$n_TXC_OD_ub,
         #n_MET_OD_ub   = l_params_all$n_MET_OD_ub,
         #n_METC_OD_ub  = l_params_all$n_METC_OD_ub,
         n_REL_OD_ub   = l_params_all$n_REL_OD_ub, 
         #n_ABS_OD_ub   = l_params_all$n_ABS_OD_ub, 
         n_TX_OD_mult_ub  = l_params_all$n_TX_OD_mult_ub,
         n_TXC_OD_mult_ub = l_params_all$n_TXC_OD_mult_ub,
         n_REL_OD_mult_ub = l_params_all$n_REL_OD_mult_ub,
         n_INJ_OD_mult_ub = l_params_all$n_INJ_OD_mult_ub,
         n_fatal_OD_ub    = l_params_all$n_fatal_OD_ub) # higher bound estimate for each param

### Number of calibration targets
v_target_names <- c("Fatal Overdoses", "Non-fatal Overdoses")
n_target       <- length(v_target_names)

#### Run IMIS algorithm ####
# CHECK THIS STEP
l_fit_imis <- IMIS(B = 10,      # n_samp = B*10 (was 100 incremental sample size at each iteration of IMIS)
                   B.re = n_resamp,      # "n_resamp" desired posterior sample size
                   number_k = 10,      # maximum number of iterations in IMIS
                   D = 0)
### Obtain posterior
m_calib_post <- l_fit_imis$resample

#### Exploring posterior distribution ####
#### Summary statistics of posterior distribution ####
### Compute posterior mean
v_calib_post_mean <- colMeans(m_calib_post)

### Compute posterior median and 95% credible interval
m_calib_post_95cr <- matrixStats::colQuantiles(m_calib_post, 
                                               probs = c(0.025, 0.5, 0.975))

### Compute posterior values for draw
v_calib_post      <- exp(log_post(m_calib_post))

### Compute maximum-a-posteriori (MAP) as the mode of the sampled values
v_calib_post_map  <- m_calib_post[which.max(v_calib_post), ]

# Summary statistics
df_posterior_summ <- data.frame(
  Parameter = v_param_names,
  Mean      = v_calib_post_mean,
  m_calib_post_95cr,
  MAP       = v_calib_post_map,
  check.names = FALSE)
df_posterior_summ

### Save summary statistics of posterior distribution
## As .RData
save(df_posterior_summ, 
     file = "outputs/summary_posterior.RData")
## As .csv
write.csv(df_posterior_summ, 
          file = "tables/summary_posterior.csv", 
          row.names = FALSE)

#### Visualization of posterior distribution ####
### Rescale posterior to plot density of plots
v_calib_alpha <- scales::rescale(v_calib_post)

### Plot the 1000 draws from the posterior
png("plots/posterior_distribution_joint.png", 
    width = 8, height = 6, units = 'in', res = 300)
s3d <- scatterplot3d::scatterplot3d(x = m_calib_post[, 1],
                                    y = m_calib_post[, 2],
                                    z = m_calib_post[, 3],
                                    color = scales::alpha("black", v_calib_alpha),
                                    xlim = c(v_lb[1], v_ub[1]), 
                                    ylim = c(v_lb[2], v_ub[2]), 
                                    zlim = c(v_lb[3], v_ub[3]),
                                    xlab = v_param_names[1], 
                                    ylab = v_param_names[2], 
                                    zlab = v_param_names[3])
## Add center of Gaussian components
s3d$points3d(l_fit_imis$center, col = "red", pch = 8)
## Add legend
legend(s3d$xyz.convert(0.05, 1.0, 5), 
       col = c("black", "red"), 
       bg = "white", pch = c(1, 8), yjust = 0, 
       legend = c("Posterior sample", "Center of Gaussian components"), 
       cex = 1.1)
dev.off()

### Plot the 1000 draws from the posterior with marginal histograms
png("plots/posterior_distribution_marginal.png", 
    width = 8, height = 6, units = 'in', res = 300)
psych::pairs.panels(m_calib_post)
dev.off()

#### Store posterior and MAP from IMIS calibration ####
save(m_calib_post,
     v_calib_post_map,
     file = "outputs/imis_output.RData")
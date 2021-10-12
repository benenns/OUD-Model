rm(list = ls()) # to clean the workspace

#### Load packages, data and functions ####
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
                                file.hiv = "data/hiv_sero.csv",
                                file.hcv = "data/hcv_sero.csv",
                                file.costs = "data/costs.csv",
                                file.crime_costs = "data/crime_costs.csv",
                                file.qalys = "data/qalys.csv")

# Load calibration inputs #
v_cali_param_names <- c("n_TX_OD",
                        "n_TXC_OD",
                        "n_REL_OD", 
                        "n_TX_OD_mult",
                        "n_TXC_OD_mult",
                        "n_REL_OD_mult",
                        "n_INJ_OD_mult",
                        "n_fatal_OD")
v_shape <- c(n_TX_OD_shape   = l_params_all$n_TX_OD_shape,
             n_TXC_OD_shape  = l_params_all$n_TXC_OD_shape,
             n_REL_OD_shape   = l_params_all$n_REL_OD_shape, 
             n_TX_OD_mult_shape  = l_params_all$n_TX_OD_mult_shape,
             n_TXC_OD_mult_shape = l_params_all$n_TXC_OD_mult_shape,
             n_REL_OD_mult_shape = l_params_all$n_REL_OD_mult_shape,
             n_INJ_OD_mult_shape = l_params_all$n_INJ_OD_mult_shape,
             n_fatal_OD_shape    = l_params_all$n_fatal_OD_shape) # lower bound estimate for each param
v_scale <- c(n_TX_OD_scale   = l_params_all$n_TX_OD_scale,
             n_TXC_OD_scale  = l_params_all$n_TXC_OD_scale,
             n_REL_OD_scale   = l_params_all$n_REL_OD_scale, 
             n_TX_OD_mult_scale  = l_params_all$n_TX_OD_mult_scale,
             n_TXC_OD_mult_scale = l_params_all$n_TXC_OD_mult_scale,
             n_REL_OD_mult_scale = l_params_all$n_REL_OD_mult_scale,
             n_INJ_OD_mult_scale = l_params_all$n_INJ_OD_mult_scale,
             n_fatal_OD_scale    = l_params_all$n_fatal_OD_scale)

#### Load calibration targets ####
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

#### Specify calibration parameters ####
### Specify seed (for reproducible sequence of random numbers)
set.seed(3730687)

### Number of random samples to obtain from the posterior distribution 
n_resamp <- 2000 # to match number of PSA draws

### Names and number of input parameters to be calibrated
v_param_names  <- c("n_TX_OD",
                    "n_TXC_OD",
                    "n_REL_OD",
                    "n_TX_OD_mult",
                    "n_TXC_OD_mult",
                    "n_REL_OD_mult",
                    "n_INJ_OD_mult",
                    "n_fatal_OD")

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

#### Plot prior vs. posterior distribution for calibration parameters ####
# Draw sample prior
m_calib_prior <- sample.prior(n_samp = 2000)

# Prepare data
df_calib_post  <- as.data.frame(m_calib_post)
df_calib_prior <- as.data.frame(m_calib_prior)
# Base overdose rates
df_calib_post_plot_base  <- gather(df_calib_post, parameter, draw, factor_key = TRUE) %>% filter(parameter == "n_BUP_OD" | parameter == "n_BUPC_OD" | parameter == "n_MET_OD" | parameter == "n_METC_OD" | parameter == "n_REL_OD")
# Overdose rate multipliers
df_calib_post_plot_mult  <- gather(df_calib_post, parameter, draw, factor_key = TRUE) %>% filter(parameter == "n_TX_OD_mult" | parameter == "n_TXC_OD_mult" | parameter == "n_REL_OD_mult" | parameter == "n_INJ_OD_mult")
# Fatal overdose rate
df_calib_post_plot_fatal <- gather(df_calib_post, parameter, draw, factor_key = TRUE) %>% filter(parameter == "n_TX_OD_mult" | parameter == "n_TXC_OD_mult" | parameter == "n_REL_OD_mult" | parameter == "n_INJ_OD_mult")

# Base overdose rates
cali_base_od_prior_post <- ggplot(df_calib_post_plot_base, 
                                  aes(x = draw, y = parameter, fill = stat(x))) +
                                  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
                                  scale_fill_viridis_c(name = "Monthly rate", option = "C") +
                                  labs(title = 'Base Overdose Rates')
pdf("Plots/Calibration/cali_base_od_prior_post.pdf", width = 8, height = 6)
cali_base_od_prior_post
dev.off()

# Overdose rate multipliers

# Fatal overdose rate
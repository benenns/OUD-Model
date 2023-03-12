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
library(grid)
library(gridExtra)
library(lattice)

# To-do: Move functions into R package for OUD model
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/calibration_functions.R")

# Load model inputs #
l_params_all <- load_all_params(file.init = "data/Calibration/init_params.csv",
                                file.init_dist = "data/init_dist.csv", # calibrate on BC OUD cohort data for 2018
                                file.mort = "data/all_cause_mortality.csv",
                                file.death_hr = "data/death_hr.csv",
                                file.frailty = "data/frailty.csv",
                                file.weibull = "data/Modified Model Specification/weibull.csv",
                                file.unconditional = "data/Modified Model Specification/unconditional.csv",
                                file.overdose = "data/overdose.csv", # includes calibration-related parameters
                                file.fentanyl = "data/Calibration/fentanyl.csv",
                                file.naloxone = "data/naloxone.csv",### R&R MODIFICATION ###
                                file.hiv = "data/hiv_sero.csv",
                                file.hcv = "data/hcv_sero.csv",
                                file.costs = "data/Modified Model Specification/costs.csv",
                                file.crime_costs = "data/Modified Model Specification/crime_costs.csv",
                                file.qalys = "data/Modified Model Specification/qalys.csv")

# Load calibration inputs #
v_cali_param_names <- c("'Overdose rate (BNX/MET)'",
                        "'Overdose rate (BNX/MET + opioid)'",
                        "'Overdose rate (opioid use)'", 
                        "'Overdose rate (opioid cessation)'",
                        "'First month mult (BNX/MET + opioid)'",
                        "'Fentanyl mult'",
                        "'Fatal overdose rate'",
                        "'Probability overdose witnessed'",
                        "'COVID-19 reduction witness'")

v_par1 <- c(n_TX_OD_shape       = l_params_all$n_TX_OD_shape,
            n_TXC_OD_shape      = l_params_all$n_TXC_OD_shape,
            n_REL_OD_shape      = l_params_all$n_REL_OD_shape,
            n_ABS_OD_shape        = l_params_all$n_ABS_OD_shape,### R&R MODIFICATION ###
            n_TXC_OD_mult_shape = l_params_all$n_TXC_OD_mult_shape,
            n_fent_OD_mult_shape  = l_params_all$n_fent_OD_mult_shape,### R&R MODIFICATION ###
            n_fatal_OD_shape    = l_params_all$n_fatal_OD_shape,
            p_witness_low           = l_params_all$p_witness_low,
            p_witness_covid_adj_low           = l_params_all$p_witness_covid_adj_low)### R&R MODIFICATION ###

v_par2 <- c(n_TX_OD_scale       = l_params_all$n_TX_OD_scale,
            n_TXC_OD_scale      = l_params_all$n_TXC_OD_scale,
            n_REL_OD_scale      = l_params_all$n_REL_OD_scale,
            n_ABS_OD_scale       = l_params_all$n_ABS_OD_scale,### R&R MODIFICATION ###
            n_TXC_OD_mult_scale = l_params_all$n_TXC_OD_mult_scale,
            n_fent_OD_mult_scale = l_params_all$n_fent_OD_mult_scale,### R&R MODIFICATION ###
            n_fatal_OD_scale    = l_params_all$n_fatal_OD_scale,
            p_witness_high           = l_params_all$p_witness_high,
            p_witness_covid_adj_high = l_params_all$p_witness_covid_adj_high)### R&R MODIFICATION ###

#### Load calibration targets ####
l_cali_targets <- list(ODF = read.csv(file = "data/cali_target_odf.csv", header = TRUE),
                       ODN = read.csv(file = "data/cali_target_odn.csv", header = TRUE))

# Max calibration periods
n_cali_max_per <- max(c(l_cali_targets$ODF$Time, l_cali_targets$ODN$Time))

#### Visualize targets ####
### TARGET 1: Overdose deaths ("ODF")
# plotrix::plotCI(x    = l_cali_targets$ODF$Time, 
#                 y    = l_cali_targets$ODF$pe, 
#                 ui   = l_cali_targets$ODF$high,
#                 li   = l_cali_targets$ODF$low,
#                 xlab = "Month", ylab = "Fatal Overdoses")

### TARGET 2: Non-fatal overdose ("ODN")
# plotrix::plotCI(x    = l_cali_targets$ODN$Time, 
#                 y    = l_cali_targets$ODN$pe, 
#                 ui   = l_cali_targets$ODN$high,
#                 li   = l_cali_targets$ODN$low,
#                 xlab = "Month", ylab = "Non-fatal Overdoses")

#### Specify calibration parameters ####
### Set seed
set.seed(3730687)

### Number of random samples to obtain from the posterior distribution 
n_resamp <- 10000 # to match number of PSA draws

### Number of calibration targets
v_target_names <- c("Fatal Overdoses", "Non-fatal Overdoses")
n_target       <- length(v_target_names)

Sys.time()
#### Run IMIS algorithm ####
l_fit_imis <- IMIS(B = 100,      # n_samp = B*10 (was 100 incremental sample size at each iteration of IMIS)
                   B.re = n_resamp,      # "n_resamp" desired posterior sample size
                   number_k = 500,      # maximum number of iterations in IMIS (originally 10)
                   D = 0) # originally 0
Sys.time()
# Unique parameter sets
n_unique <- length(unique(l_fit_imis$resample[,1])) # 6299
# Effective sample size
n_ess <- round(sum(table(l_fit_imis$resample[,1]))^2/ sum(table(l_fit_imis$resample[,1])^2), 0) # 4568
# Max weight
n_max_wt <- max(table(l_fit_imis$resample[,1]))/sum(table(l_fit_imis$resample[,1])) # 0.0009

# Calibration stats
df_cali_stats <- data.frame(n_unique, n_ess, n_max_wt)
df_cali_stats

### Save calibration stats
## As .RData
save(df_cali_stats, 
     file = "outputs/Calibration/cali_stats.RData")
## As .csv
write.csv(df_cali_stats, 
          file = "outputs/Calibration/cali_stats.csv", 
          row.names = FALSE)

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
Sys.time()
### Compute maximum-a-posteriori (MAP) as the mode of the sampled values
v_calib_post_map  <- m_calib_post[which.max(v_calib_post), ]

# Summary statistics
df_posterior_summ <- data.frame(
  Parameter = v_cali_param_names,
  Mean      = v_calib_post_mean,
  m_calib_post_95cr,
  MAP       = v_calib_post_map,
  check.names = FALSE)
df_posterior_summ

### Save summary statistics of posterior distribution
## As .RData
save(df_posterior_summ, 
     file = "outputs/Calibration/summary_posterior.RData")
## As .csv
write.csv(df_posterior_summ, 
          file = "outputs/Calibration/summary_posterior.csv", 
          row.names = FALSE)

#### Visualization of posterior distribution ####
### Rescale posterior to plot density of plots
v_calib_alpha <- scales::rescale(v_calib_post)

### Plot the 1000 draws from the posterior
png("plots/Calibration/posterior_distribution_joint.png", 
    width = 8, height = 6, units = 'in', res = 300)
s3d <- scatterplot3d::scatterplot3d(x = m_calib_post[, 1],
                                    y = m_calib_post[, 2],
                                    z = m_calib_post[, 3],
                                    color = scales::alpha("black", v_calib_alpha),
                                    xlim = c(v_lb[1], v_ub[1]), 
                                    ylim = c(v_lb[2], v_ub[2]), 
                                    zlim = c(v_lb[3], v_ub[3]),
                                    xlab = v_cali_param_names[1], 
                                    ylab = v_cali_param_names[2], 
                                    zlab = v_cali_param_names[3])
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
png("plots/Calibration/posterior_distribution_marginal.png", 
    width = 8, height = 6, units = 'in', res = 300)
psych::pairs.panels(m_calib_post)
dev.off()

#### Store posterior and MAP from IMIS calibration ####
save(m_calib_post,
     v_calib_post_map,
     file = "outputs/Calibration/imis_output.RData")

#### Plot prior vs. posterior distribution for calibration parameters ####
# Load posterior
imis_output <- load(file = "outputs/Calibration/imis_output.RData")

# Draw sample prior
m_calib_prior <- sample.prior(n_samp = n_resamp)

# Prepare data
df_samp_prior <- melt(cbind(Distribution = "Prior", 
                            as.data.frame(m_calib_prior[1:n_resamp, ])), 
                            variable.name = "Parameter")
df_samp_post_imis  <- melt(cbind(Distribution = "Posterior", 
                                 as.data.frame(m_calib_post[1:n_resamp, ])), 
                                 variable.name = "Parameter")

df_calib_prior_post <- rbind(df_samp_prior, df_samp_post_imis)


df_calib_prior_post$Distribution <- ordered(df_calib_prior_post$Distribution, 
                                            levels = c("Prior", 
                                                       "Posterior"))

df_calib_prior_post$Parameter <- factor(df_calib_prior_post$Parameter,
                                        levels = levels(df_calib_prior_post$Parameter),
                                        ordered = TRUE,
                                        labels = v_cali_param_names)

### Plot priors and IMIS posteriors
# TO-DO: Add vertical lines for prior mean and MAP
prior_v_posterior <- ggplot(df_calib_prior_post, 
                         aes(x = value, y = ..density.., fill = Distribution)) +
  facet_wrap(~Parameter, scales = "free", 
             ncol = 3,
             labeller = label_parsed) +
  #geom_vline(data = df_posterior_summ,
  #           aes(xintercept = value, linetype = Type, color = Type)) +
  #scale_x_continuous(breaks = dampack::number_ticks(5)) +
  scale_color_manual("", values = c("black", "navy blue", "tomato")) +
  geom_density(alpha=0.5) +
  theme_bw(base_size = 16) +
  guides(fill = guide_legend(title = "", order = 1),
         linetype = guide_legend(title = "", order = 2),
         color = guide_legend(title = "", order = 2)) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "white",
                                        color = "white"),
        strip.text = element_text(hjust = 0))
prior_v_posterior

ggsave(prior_v_posterior, 
       filename = "Plots/Calibration/prior-v-posterior.pdf", 
       width = 10, height = 7)
ggsave(prior_v_posterior, 
       filename = "Plots/Calibration/prior-v-posterior.png", 
       width = 10, height = 7)

#### Plot model fit against calibration targets ####
# Run model for n_samp posterior distribution draws
# Output list of fatal and total overdoses at T = 1, T = 2, T = 3
m_model_targets_ODF <- m_model_targets_ODN <- matrix(0, nrow = n_resamp, ncol = 4) 

for(i in 1:n_resamp){
  l_model_target_fit <- calibration_out(v_params_calib = m_calib_post[i, ], 
                                        l_params_all = l_params_all)
  
  m_model_targets_ODF[i, 1] <- l_model_target_fit$fatal_overdose[1]
  m_model_targets_ODF[i, 2] <- l_model_target_fit$fatal_overdose[2]
  m_model_targets_ODF[i, 3] <- l_model_target_fit$fatal_overdose[3]
  ### R&R MODIFICATION ###
  m_model_targets_ODF[i, 4] <- l_model_target_fit$fatal_overdose[4]
  
  m_model_targets_ODN[i, 1] <- l_model_target_fit$overdose[1]
  m_model_targets_ODN[i, 2] <- l_model_target_fit$overdose[2]
  m_model_targets_ODN[i, 3] <- l_model_target_fit$overdose[3]
  ### R&R MODIFICATION ###
  m_model_targets_ODN[i, 4] <- l_model_target_fit$overdose[4]
}
Sys.time()
## As .RData
save(m_model_targets_ODF, 
     file = "outputs/Calibration/model_targets_ODF.RData")
save(m_model_targets_ODN, 
     file = "outputs/Calibration/model_targets_ODN.RData")

# Model outputs
m_model_targets_ODF_stats <- cbind(matrixStats::colQuantiles(m_model_targets_ODF, 
                                                             probs = c(0.025, 0.5, 0.975)),
                               matrixStats::colMeans2(m_model_targets_ODF))

m_model_targets_ODN_stats <- cbind(matrixStats::colQuantiles(m_model_targets_ODN, 
                                                             probs = c(0.025, 0.5, 0.975)),
                                   matrixStats::colMeans2(m_model_targets_ODN))

m_time <- matrix(c(12, 24, 36, 48))
m_pop <- matrix(l_cali_targets$ODF$Pop)
m_model_targets_ODF_fit <- cbind(m_model_targets_ODF_stats, m_time, m_pop)
m_model_targets_ODN_fit <- cbind(m_model_targets_ODN_stats, m_time, m_pop)

df_model_targets_ODF_fit <- m_model_targets_ODF_fit %>% as_tibble() %>% setNames(c("ci_low", "Median", "ci_high", "pe", "Time", "Pop")) %>% mutate(Target = "Model output (95% CI)", 
                                                                                                                                                   Num = pe * Pop,
                                                                                                                                                   low = ci_low * Pop,
                                                                                                                                                   high = ci_high * Pop)
df_model_targets_ODN_fit <- m_model_targets_ODN_fit %>% as_tibble() %>% setNames(c("ci_low", "Median", "ci_high", "pe", "Time", "Pop")) %>% mutate(Target = "Model output (95% CI)",
                                                                                                                                                   Num = pe * Pop,
                                                                                                                                                   low = ci_low * Pop,
                                                                                                                                                   high = ci_high * Pop)

# Targets
df_targets_ODF <- l_cali_targets$ODF %>% as_tibble() %>% mutate(Target = "Cali target (95% CI)",
                                                                low = low * Pop,
                                                                high = high * Pop)
df_targets_ODN <- l_cali_targets$ODN %>% as_tibble() %>% mutate(Target = "Cali target (95% CI)",
                                                                low = low * Pop,
                                                                high = high * Pop)

# Combine
df_fit_ODF <- bind_rows(df_targets_ODF, df_model_targets_ODF_fit)
df_fit_ODN <- bind_rows(df_targets_ODN, df_model_targets_ODN_fit)

# Plot fit vs. targets
# Fatal overdose
p_temp_ODF <- ggplot(df_fit_ODF, aes(x = Time, y = Num, group = Target, color = Target)) + 
              geom_line() +
              geom_point()+
              geom_errorbar(aes(ymin = low, ymax = high), width = .5,
                            position = position_dodge(0.05))

plot_fit_ODF <- p_temp_ODF + labs(title = NULL, x = "Year", y = "Fatal overdoses") +
                             theme_classic() +
                             theme(legend.position="none") + 
                             theme(legend.title = element_blank()) +
                             scale_color_manual(values = c('#999999','#E69F00')) +
                             scale_x_continuous(breaks = c(12, 24, 36, 48),
                                                labels = c("2017", "2018", "2019", "2020"))

# Plot for extracting legend only
plot_fit_ODF_leg <- p_temp_ODF + labs(title = NULL, x = "Year", y = "Fatal overdoses") +
  theme_classic() +
  scale_color_manual(values = c('#999999','#E69F00')) +
  scale_x_continuous(breaks = c(12, 24, 36, 48),
                     labels = c("2017", "2018", "2019", "2020"))

# Non-fatal overdose
p_temp_ODN <- ggplot(df_fit_ODN, aes(x = Time, y = Num, group = Target, color = Target)) + 
              geom_line() +
              geom_point()+
              geom_errorbar(aes(ymin = low, ymax = high), width = .5,
                            position = position_dodge(0.05))

plot_fit_ODN <- p_temp_ODN + labs(title = NULL, x = "Year", y = "Non-fatal overdoses")+
                             theme_classic() +
                             theme(legend.position = "none") + 
                             theme(legend.title = element_blank()) +
                             scale_color_manual(values = c('#999999','#E69F00')) +
                             scale_x_continuous(breaks = c(12, 24, 36, 48),
                                                labels = c("2017", "2018", "2019", "2020"))

plot_fit_ODN
# Code to extract legend from plots
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(plot_fit_ODF_leg)

# Combined
#plot_fit_comb <- grid.arrange(plot_fit_ODF, plot_fit_ODN, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

plot_fit_comb <- grid.arrange(arrangeGrob(plot_fit_ODF, plot_fit_ODN, nrow = 1),
                                          mylegend, nrow = 2, heights = c(6, 1))

plot_fit_comb
# Outputs
ggsave(plot_fit_ODF, 
       filename = "Plots/Calibration/target-fit-ODF.png", 
       width = 4, height = 4)
ggsave(plot_fit_ODN, 
       filename = "Plots/Calibration/target-fit-ODN.png", 
       width = 4, height = 4)
ggsave(plot_fit_comb, 
       filename = "Plots/Calibration/target-fit-comb.png", 
       width = 8, height = 4)

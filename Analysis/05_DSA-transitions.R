rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(data.table)
library(formattable)
library(tidyr)
library(RColorBrewer)

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/ICER_functions.R")

# Load parameters
source("Analysis/00_load_parameters.R")

# Load DSA parameters
###################
### Transitions ###
###################
df_dsa_frailty <- read.csv(file = "data/DSA/frailty.csv", row.names = 1, header = TRUE)
#df_dsa_frailty_TS <- read.csv(file = "data/DSA/Trial Specification/frailty.csv", row.names = 1, header = TRUE)

v_dsa_frailty_episode <- unlist(df_dsa_frailty["pe_episode_frailty",])
v_dsa_frailty_concurrent <- unlist(df_dsa_frailty["pe_concurrent_frailty",])
v_dsa_frailty_inj <- unlist(df_dsa_frailty["pe_inj_frailty",])
v_dsa_frailty_kurz <- unlist(df_dsa_frailty["pe_kurz_frailty",])

### BNX threshold SA ###
df_dsa_threshold_MMS <- read.csv(file = "data/DSA/Modified Model Specification/threshold.csv", row.names = 1, header = TRUE)
# Initialize matrices
v_threshold_names_MMS <- colnames(df_dsa_threshold_MMS)
v_threshold_rownames_MMS <- rownames(df_dsa_threshold_MMS)
m_dsa_threshold_MMS <- array(0, dim = c(nrow(df_dsa_threshold_MMS), length(df_dsa_threshold_MMS)),
                             dimnames = list(v_threshold_rownames_MMS, v_threshold_names_MMS))
## Threshold SA ##
# MMS
for (i in 1:nrow(df_dsa_threshold_MMS)){
  m_dsa_threshold_MMS[i,] <- unlist(df_dsa_threshold_MMS[i,])
}

############################################
#### Deterministic sensitivity analysis ####
############################################

################
### Baseline ###
################
# MMS
l_outcomes_MET_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map)
l_outcomes_BUP_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map)
ICER_MMS <- ICER(outcomes_comp = l_outcomes_MET_MMS, outcomes_int = l_outcomes_BUP_MMS)

###################
### Transitions ###
###################
## Frailty ##
# Episode
l_outcomes_MET_frailty_episode_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_episode)
l_outcomes_BUP_frailty_episode_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_episode)
ICER_frailty_episode_MMS <- ICER(outcomes_comp = l_outcomes_MET_frailty_episode_MMS, outcomes_int = l_outcomes_BUP_frailty_episode_MMS)

# Concurrent opioid use
l_outcomes_MET_frailty_concurrent_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_concurrent)
l_outcomes_BUP_frailty_concurrent_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_concurrent)
ICER_frailty_concurrent_MMS <- ICER(outcomes_comp = l_outcomes_MET_frailty_concurrent_MMS, outcomes_int = l_outcomes_BUP_frailty_concurrent_MMS)

# Injection multiplier
l_outcomes_MET_frailty_inj_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_inj)
l_outcomes_BUP_frailty_inj_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_inj)
ICER_frailty_inj_MMS <- ICER(outcomes_comp = l_outcomes_MET_frailty_inj_MMS, outcomes_int = l_outcomes_BUP_frailty_inj_MMS)

# Kurz results (combined for treatment states)
l_outcomes_MET_frailty_kurz_MMS <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_kurz)
l_outcomes_BUP_frailty_kurz_MMS <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_frailty_kurz)
ICER_frailty_kurz_MMS <- ICER(outcomes_comp = l_outcomes_MET_frailty_kurz_MMS, outcomes_int = l_outcomes_BUP_frailty_kurz_MMS)

###############################
### BNX Retention Threshold ###
###############################

# Initialize lists
l_outcomes_MET_threshold_MMS <- list()
l_outcomes_BUP_threshold_MMS <- list()
l_ICER_threshold_MMS <- list()

## Treatment retention (threshold SA for BNX retention) ##
# MMS
for (i in 1:nrow(m_dsa_threshold_MMS)){
  # +i%
  l_outcomes_MET_threshold_MMS[[i]] <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_MMS[i,])
  l_outcomes_BUP_threshold_MMS[[i]] <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_MMS[i,])
  l_ICER_threshold_MMS[[i]] <- ICER(outcomes_comp = l_outcomes_MET_threshold_MMS[[i]], outcomes_int = l_outcomes_BUP_threshold_MMS[[i]])
}

################
### Baseline ###
################
df_baseline_MMS <- data.frame(ICER_MMS$df_incremental, ICER_MMS$df_icer)

###################
### Transitions ###
###################
# Baseline
df_transitions_baseline_MMS <- data.frame(ICER_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_MMS$df_incremental$n_inc_qalys_TOTAL_life, 
                                          ICER_MMS$df_icer$n_icer_TOTAL_life)

# Costs
v_transitions_frailty_episode_MMS <- c(ICER_frailty_episode_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_frailty_episode_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_transitions_frailty_concurrent_MMS <- c(ICER_frailty_concurrent_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_frailty_concurrent_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_transitions_frailty_inj_MMS <- c(ICER_frailty_inj_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_frailty_inj_MMS$df_incremental$n_inc_costs_TOTAL_life)
v_transitions_frailty_kurz_MMS <- c(ICER_frailty_kurz_MMS$df_incremental$n_inc_costs_TOTAL_life, ICER_frailty_kurz_MMS$df_incremental$n_inc_costs_TOTAL_life)

m_transitions_costs_MMS <- rbind(v_transitions_frailty_episode_MMS, v_transitions_frailty_concurrent_MMS, v_transitions_frailty_inj_MMS, v_transitions_frailty_kurz_MMS)

df_transitions_costs_MMS <- as.data.frame(m_transitions_costs_MMS)
colnames(df_transitions_costs_MMS) <- c("Lower", "Upper")
df_transitions_costs_MMS <- as_data_frame(df_transitions_costs_MMS) %>% mutate(diff = abs(Upper - Lower),
                                                                         base = ICER_MMS$df_incremental$n_inc_costs_TOTAL_life) %>%
  add_column(var_name = c("All episode frailty equal", "No difference for concurrent use", "No difference for injection", "Combined TX from Kurz (2021)"))

# QALYs
v_transitions_frailty_episode_MMS <- c(ICER_frailty_episode_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_frailty_episode_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_transitions_frailty_concurrent_MMS <- c(ICER_frailty_concurrent_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_frailty_concurrent_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_transitions_frailty_inj_MMS <- c(ICER_frailty_inj_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_frailty_inj_MMS$df_incremental$n_inc_qalys_TOTAL_life)
v_transitions_frailty_kurz_MMS <- c(ICER_frailty_kurz_MMS$df_incremental$n_inc_qalys_TOTAL_life, ICER_frailty_kurz_MMS$df_incremental$n_inc_qalys_TOTAL_life)

m_transitions_qalys_MMS <- rbind(v_transitions_frailty_episode_MMS, v_transitions_frailty_concurrent_MMS, v_transitions_frailty_inj_MMS, v_transitions_frailty_kurz_MMS)

df_transitions_qalys_MMS <- as.data.frame(m_transitions_qalys_MMS)
colnames(df_transitions_qalys_MMS) <- c("Lower", "Upper")
df_transitions_qalys_MMS <- as_data_frame(df_transitions_qalys_MMS) %>% mutate(diff = abs(Upper - Lower),
                                                                               base = ICER_MMS$df_incremental$n_inc_qalys_TOTAL_life) %>%
  add_column(var_name = c("All episode frailty equal", "No difference for concurrent use", "No difference for injection", "Combined TX from Kurz (2021)"))

##########################
#### BNX Threshold SA ####
##########################
df_threshold_MMS <- data.frame()
df_threshold_MMS_temp <- data.frame()

v_names <- c("n_inc_costs_TOTAL_life", "n_inc_qalys_TOTAL_life", "n_icer_TOTAL_life")

# MMS
for (i in 1:nrow(m_dsa_threshold_MMS)){
  df_threshold_MMS_temp <- data.frame(l_ICER_threshold_MMS[[i]]$df_incremental$n_inc_costs_TOTAL_life, l_ICER_threshold_MMS[[i]]$df_incremental$n_inc_qalys_TOTAL_life,
                                      l_ICER_threshold_MMS[[i]]$df_icer$n_icer_TOTAL_life)
  
  df_threshold_MMS <- rbind(df_threshold_MMS, df_threshold_MMS_temp)
}

rownames(df_threshold_MMS) <- v_threshold_rownames_MMS
colnames(df_threshold_MMS) <- v_names

df_threshold_MMS <- cbind(v_threshold_rownames_MMS, df_threshold_MMS)

# Save outputs
## As .RData ##
save(df_threshold_MMS, 
     file = "outputs/DSA/Modified Model Specification/df_threshold_MMS.RData")

#########################
#### Tornado Diagram ####
#########################
# Costs
v_order_parameters <- df_transitions_costs_MMS %>% arrange(diff) %>%
  mutate(var_name = factor(x = var_name, levels = var_name)) %>%
  select(var_name) %>% unlist() %>% levels()

# width of columns in plot (value between 0 and 1)
width <- 0.75
# get data frame in shape for ggplot and geom_rect
df.2 <- df_transitions_costs_MMS %>% 
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key = 'type', value = 'output.value', Lower:Upper) %>%
  # just reordering columns
  select(var_name, type, output.value, diff, base) %>%
  # create the columns for geom_rect
  mutate(var_name = factor(var_name, levels = v_order_parameters),
         ymin = pmin(output.value, base),
         ymax = pmax(output.value, base),
         xmin = as.numeric(var_name) - width/2,
         xmax = as.numeric(var_name) + width/2)

# create plot
# (use scale_x_continuous to change labels in y axis to name of parameters)
p_tornado_transitions_costs <- ggplot() + 
  geom_rect(data = df.2, 
            aes(ymax = ymax, ymin = ymin, xmax = xmax, xmin = xmin, fill = type)) +
  theme_bw() + 
  scale_fill_manual(values = c("Upper" = "#5ab4ac",
                               "Lower" = "#d8b365")) +
  theme(axis.title.y = element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  geom_hline(yintercept = df.2$base) +
  scale_x_continuous(breaks = c(1:length(v_order_parameters)), 
                     labels = v_order_parameters) +
  xlab("Parameter") + ylab("Incremental Cost") +
  coord_flip()

png(file = "Plots/DSA/Modified Model Spec/tornado_transitions_costs.png", width = 600, height = 600)
p_tornado_transitions_costs
dev.off()

# QALYs
v_order_parameters <- df_transitions_qalys_MMS %>% arrange(diff) %>%
  mutate(var_name = factor(x = var_name, levels = var_name)) %>%
  select(var_name) %>% unlist() %>% levels()

# width of columns in plot (value between 0 and 1)
width <- 0.75
# get data frame in shape for ggplot and geom_rect
df.2 <- df_transitions_qalys_MMS %>% 
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key = 'type', value = 'output.value', Lower:Upper) %>%
  # just reordering columns
  select(var_name, type, output.value, diff, base) %>%
  # create the columns for geom_rect
  mutate(var_name = factor(var_name, levels = v_order_parameters),
         ymin = pmin(output.value, base),
         ymax = pmax(output.value, base),
         xmin = as.numeric(var_name) - width/2,
         xmax = as.numeric(var_name) + width/2)

# create plot
# (use scale_x_continuous to change labels in y axis to name of parameters)
#png(file = "Plots/DSA/Modified Model Spec/tornado_overdose_qalys.png", width = 960, height = 540)
p_tornado_transitions_qalys <- ggplot() + 
  geom_rect(data = df.2, 
            aes(ymax = ymax, ymin = ymin, xmax = xmax, xmin = xmin, fill = type)) +
  theme_bw() + 
  scale_fill_manual(values = c("Upper" = "#5ab4ac",
                               "Lower" = "#d8b365")) +
  theme(axis.title.y=element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  geom_hline(yintercept = df.2$base) +
  scale_x_continuous(breaks = c(1:length(v_order_parameters)), 
                     labels = v_order_parameters) +
  xlab("Parameter") + ylab("Incremental Cost") +
  coord_flip()
#dev.off()

png(file = "Plots/DSA/Modified Model Spec/tornado_transitions_qalys.png", width = 600, height = 600)
p_tornado_transitions_qalys
dev.off()

#########################
#### Threshold plots ####
#########################

# Prepare data for plotting
df_threshold_MMS <- df_threshold_MMS %>% as_tibble() %>% mutate(perc_increase = as.numeric(v_threshold_rownames_MMS))

# MMS
df_threshold_qalys_MMS <- df_threshold_MMS %>% gather("scenario", "inc_qalys", n_inc_qalys_TOTAL_life, na.rm = FALSE, convert = FALSE) %>%
  select(perc_increase, scenario, inc_qalys)

df_threshold_costs_MMS <- df_threshold_MMS %>% gather("scenario", "inc_costs", n_inc_costs_TOTAL_life, na.rm = FALSE, convert = FALSE) %>%
  select(perc_increase, scenario, inc_costs)

## Threshold plots ##
# MMS
# Incremental QALYs
plot_DSA_qalys_MMS_threshold <- ggplot(df_threshold_qalys_MMS, aes(x = perc_increase, y = inc_qalys, group = scenario)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_line(aes(color = scenario)) +
  geom_hline(yintercept = 0) +
  xlab("Percentage Increase in BNX Episode Duration") + ylab("Incremental QALYs") +
  xlim(0, 200) #+
#ylim(-0.01, 0.025)

#plot_DSA_qalys_MMS_threshold

ggsave(plot_DSA_qalys_MMS_threshold, 
       filename = "Plots/DSA/Threshold SA/DSA-BNX-threshold-qalys-MMS.png", 
       width = 6, height = 6)

# Incremental Costs
plot_DSA_costs_MMS_threshold <- ggplot(df_threshold_costs_MMS, aes(x = perc_increase, y = inc_costs, group = scenario)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_line(aes(color = scenario)) +
  geom_hline(yintercept = 0) +
  #scale_x_continuous(limits = c(0, 200))
  xlab("Percentage Increase in BNX Episode Duration") + ylab("Incremental Costs") +
  xlim(0, 200) #+
# ylim(-0.01, 0.025)

#plot_DSA_qalys_MMS_threshold

ggsave(plot_DSA_costs_MMS_threshold, 
       filename = "Plots/DSA/Threshold SA/DSA-BNX-threshold-costs-MMS.png", 
       width = 6, height = 6)

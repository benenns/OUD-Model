rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(data.table)
library(formattable)
library(tidyr)
library(RColorBrewer)
library(Rmisc)
library(grid)
library(gridExtra)
library(lattice)

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/ICER_functions.R")

# Load parameters
source("Analysis/00_load_parameters.R")

# Load DSA parameters
### BNX threshold SA ###
df_dsa_threshold_CAN_MMS <- read.csv(file = "data/DSA/Modified Model Specification/threshold_CAN.csv", row.names = 1, header = TRUE)
df_dsa_threshold_BC_MMS <- read.csv(file = "data/DSA/Modified Model Specification/threshold_BC.csv", row.names = 1, header = TRUE)
df_dsa_threshold_AB_MMS <- read.csv(file = "data/DSA/Modified Model Specification/threshold_AB.csv", row.names = 1, header = TRUE)
df_dsa_threshold_ON_MMS <- read.csv(file = "data/DSA/Modified Model Specification/threshold_ON.csv", row.names = 1, header = TRUE)
df_dsa_threshold_QC_MMS <- read.csv(file = "data/DSA/Modified Model Specification/threshold_QC.csv", row.names = 1, header = TRUE)

# Initialize matrices
v_threshold_names_MMS <- colnames(df_dsa_threshold_CAN_MMS)
v_threshold_rownames_MMS <- rownames(df_dsa_threshold_CAN_MMS)
m_dsa_threshold_CAN_MMS <- array(0, dim = c(nrow(df_dsa_threshold_CAN_MMS), length(df_dsa_threshold_CAN_MMS)),
                                 dimnames = list(v_threshold_rownames_MMS, v_threshold_names_MMS))
m_dsa_threshold_BC_MMS <- array(0, dim = c(nrow(df_dsa_threshold_BC_MMS), length(df_dsa_threshold_BC_MMS)),
                                 dimnames = list(v_threshold_rownames_MMS, v_threshold_names_MMS))
m_dsa_threshold_AB_MMS <- array(0, dim = c(nrow(df_dsa_threshold_AB_MMS), length(df_dsa_threshold_AB_MMS)),
                                 dimnames = list(v_threshold_rownames_MMS, v_threshold_names_MMS))
m_dsa_threshold_ON_MMS <- array(0, dim = c(nrow(df_dsa_threshold_ON_MMS), length(df_dsa_threshold_ON_MMS)),
                                 dimnames = list(v_threshold_rownames_MMS, v_threshold_names_MMS))
m_dsa_threshold_QC_MMS <- array(0, dim = c(nrow(df_dsa_threshold_QC_MMS), length(df_dsa_threshold_QC_MMS)),
                                 dimnames = list(v_threshold_rownames_MMS, v_threshold_names_MMS))

## Threshold SA ##
for (i in 1:nrow(df_dsa_threshold_CAN_MMS)){
  m_dsa_threshold_CAN_MMS[i,] <- unlist(df_dsa_threshold_CAN_MMS[i,])
}
for (i in 1:nrow(df_dsa_threshold_BC_MMS)){
  m_dsa_threshold_BC_MMS[i,] <- unlist(df_dsa_threshold_BC_MMS[i,])
}
for (i in 1:nrow(df_dsa_threshold_AB_MMS)){
  m_dsa_threshold_AB_MMS[i,] <- unlist(df_dsa_threshold_AB_MMS[i,])
}
for (i in 1:nrow(df_dsa_threshold_ON_MMS)){
  m_dsa_threshold_ON_MMS[i,] <- unlist(df_dsa_threshold_ON_MMS[i,])
}
for (i in 1:nrow(df_dsa_threshold_QC_MMS)){
  m_dsa_threshold_QC_MMS[i,] <- unlist(df_dsa_threshold_QC_MMS[i,])
}

###############################
### BNX Retention Threshold ###
###############################

# Initialize lists
l_outcomes_MET_threshold_CAN_MMS <- l_outcomes_MET_threshold_BC_MMS <- l_outcomes_MET_threshold_AB_MMS <- l_outcomes_MET_threshold_ON_MMS <- l_outcomes_MET_threshold_QC_MMS <- list()
l_outcomes_BUP_threshold_CAN_MMS <- l_outcomes_BUP_threshold_BC_MMS <- l_outcomes_BUP_threshold_AB_MMS <- l_outcomes_BUP_threshold_ON_MMS <- l_outcomes_BUP_threshold_QC_MMS <- list()
l_ICER_threshold_CAN_MMS <- l_ICER_threshold_BC_MMS <- l_ICER_threshold_AB_MMS <- l_ICER_threshold_ON_MMS <- l_ICER_threshold_QC_MMS <- list()

## Treatment retention (threshold SA for BNX retention) ##
# Canada
for (i in 1:nrow(m_dsa_threshold_CAN_MMS)){
  # +i%
  l_outcomes_MET_threshold_CAN_MMS[[i]] <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_CAN_MMS[i,])
  l_outcomes_BUP_threshold_CAN_MMS[[i]] <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_CAN_MMS[i,])
  l_ICER_threshold_CAN_MMS[[i]] <- ICER(outcomes_comp = l_outcomes_MET_threshold_CAN_MMS[[i]], outcomes_int = l_outcomes_BUP_threshold_CAN_MMS[[i]])
}

# BC
for (i in 1:nrow(m_dsa_threshold_BC_MMS)){
  # +i%
  l_outcomes_MET_threshold_BC_MMS[[i]] <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_BC_MMS[i,])
  l_outcomes_BUP_threshold_BC_MMS[[i]] <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_BC_MMS[i,])
  l_ICER_threshold_BC_MMS[[i]] <- ICER(outcomes_comp = l_outcomes_MET_threshold_BC_MMS[[i]], outcomes_int = l_outcomes_BUP_threshold_BC_MMS[[i]])
}

# AB
for (i in 1:nrow(m_dsa_threshold_AB_MMS)){
  # +i%
  l_outcomes_MET_threshold_AB_MMS[[i]] <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_AB_MMS[i,])
  l_outcomes_BUP_threshold_AB_MMS[[i]] <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_AB_MMS[i,])
  l_ICER_threshold_AB_MMS[[i]] <- ICER(outcomes_comp = l_outcomes_MET_threshold_AB_MMS[[i]], outcomes_int = l_outcomes_BUP_threshold_AB_MMS[[i]])
}

# ON
for (i in 1:nrow(m_dsa_threshold_ON_MMS)){
  # +i%
  l_outcomes_MET_threshold_ON_MMS[[i]] <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_ON_MMS[i,])
  l_outcomes_BUP_threshold_ON_MMS[[i]] <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_ON_MMS[i,])
  l_ICER_threshold_ON_MMS[[i]] <- ICER(outcomes_comp = l_outcomes_MET_threshold_ON_MMS[[i]], outcomes_int = l_outcomes_BUP_threshold_ON_MMS[[i]])
}

# QC
for (i in 1:nrow(m_dsa_threshold_QC_MMS)){
  # +i%
  l_outcomes_MET_threshold_QC_MMS[[i]] <- outcomes(l_params_all = l_params_MET_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_QC_MMS[i,])
  l_outcomes_BUP_threshold_QC_MMS[[i]] <- outcomes(l_params_all = l_params_BUP_MMS, v_params_calib = v_calib_post_map, v_params_dsa = m_dsa_threshold_QC_MMS[i,])
  l_ICER_threshold_QC_MMS[[i]] <- ICER(outcomes_comp = l_outcomes_MET_threshold_QC_MMS[[i]], outcomes_int = l_outcomes_BUP_threshold_QC_MMS[[i]])
}

##########################
#### BNX Threshold SA ####
##########################
df_threshold_CAN_MMS <- data.frame()
df_threshold_BC_MMS <- data.frame()
df_threshold_AB_MMS <- data.frame()
df_threshold_ON_MMS <- data.frame()
df_threshold_QC_MMS <- data.frame()
df_threshold_MMS_temp <- data.frame()

v_names <- c("n_inc_costs_TOTAL_life", "n_inc_qalys_TOTAL_life", "n_icer_TOTAL_life")

# Canada
for (i in 1:nrow(m_dsa_threshold_CAN_MMS)){
  df_threshold_MMS_temp <- data.frame(l_ICER_threshold_CAN_MMS[[i]]$df_incremental$n_inc_costs_TOTAL_life, l_ICER_threshold_CAN_MMS[[i]]$df_incremental$n_inc_qalys_TOTAL_life,
                                      l_ICER_threshold_CAN_MMS[[i]]$df_icer$n_icer_TOTAL_life)
  
  df_threshold_CAN_MMS <- rbind(df_threshold_CAN_MMS, df_threshold_MMS_temp)
}

rownames(df_threshold_CAN_MMS) <- v_threshold_rownames_MMS
colnames(df_threshold_CAN_MMS) <- v_names

df_threshold_CAN_MMS <- cbind(v_threshold_rownames_MMS, df_threshold_CAN_MMS)

# BC
for (i in 1:nrow(m_dsa_threshold_BC_MMS)){
  df_threshold_MMS_temp <- data.frame(l_ICER_threshold_BC_MMS[[i]]$df_incremental$n_inc_costs_TOTAL_life, l_ICER_threshold_BC_MMS[[i]]$df_incremental$n_inc_qalys_TOTAL_life,
                                      l_ICER_threshold_BC_MMS[[i]]$df_icer$n_icer_TOTAL_life)
  
  df_threshold_BC_MMS <- rbind(df_threshold_BC_MMS, df_threshold_MMS_temp)
}

rownames(df_threshold_BC_MMS) <- v_threshold_rownames_MMS
colnames(df_threshold_BC_MMS) <- v_names

df_threshold_BC_MMS <- cbind(v_threshold_rownames_MMS, df_threshold_BC_MMS)

# AB
for (i in 1:nrow(m_dsa_threshold_AB_MMS)){
  df_threshold_MMS_temp <- data.frame(l_ICER_threshold_AB_MMS[[i]]$df_incremental$n_inc_costs_TOTAL_life, l_ICER_threshold_AB_MMS[[i]]$df_incremental$n_inc_qalys_TOTAL_life,
                                      l_ICER_threshold_AB_MMS[[i]]$df_icer$n_icer_TOTAL_life)
  
  df_threshold_AB_MMS <- rbind(df_threshold_AB_MMS, df_threshold_MMS_temp)
}

rownames(df_threshold_AB_MMS) <- v_threshold_rownames_MMS
colnames(df_threshold_AB_MMS) <- v_names

df_threshold_AB_MMS <- cbind(v_threshold_rownames_MMS, df_threshold_AB_MMS)

# ON
for (i in 1:nrow(m_dsa_threshold_ON_MMS)){
  df_threshold_MMS_temp <- data.frame(l_ICER_threshold_ON_MMS[[i]]$df_incremental$n_inc_costs_TOTAL_life, l_ICER_threshold_ON_MMS[[i]]$df_incremental$n_inc_qalys_TOTAL_life,
                                      l_ICER_threshold_ON_MMS[[i]]$df_icer$n_icer_TOTAL_life)
  
  df_threshold_ON_MMS <- rbind(df_threshold_ON_MMS, df_threshold_MMS_temp)
}

rownames(df_threshold_ON_MMS) <- v_threshold_rownames_MMS
colnames(df_threshold_ON_MMS) <- v_names

df_threshold_ON_MMS <- cbind(v_threshold_rownames_MMS, df_threshold_ON_MMS)

# QC
for (i in 1:nrow(m_dsa_threshold_QC_MMS)){
  df_threshold_MMS_temp <- data.frame(l_ICER_threshold_QC_MMS[[i]]$df_incremental$n_inc_costs_TOTAL_life, l_ICER_threshold_QC_MMS[[i]]$df_incremental$n_inc_qalys_TOTAL_life,
                                      l_ICER_threshold_QC_MMS[[i]]$df_icer$n_icer_TOTAL_life)
  
  df_threshold_QC_MMS <- rbind(df_threshold_QC_MMS, df_threshold_MMS_temp)
}

rownames(df_threshold_QC_MMS) <- v_threshold_rownames_MMS
colnames(df_threshold_QC_MMS) <- v_names

df_threshold_QC_MMS <- cbind(v_threshold_rownames_MMS, df_threshold_QC_MMS)

# Combine scenarios
df_threshold_MMS <- bind_rows(list(df_threshold_CAN_MMS, df_threshold_BC_MMS, df_threshold_AB_MMS, df_threshold_ON_MMS, df_threshold_QC_MMS), .id = "scenario") %>%
  mutate(prov = if_else(scenario == 1, "Canada",
                if_else(scenario == 2, "BC",
                if_else(scenario == 3, "Alberta",
                if_else(scenario == 4, "Ontario",
                if_else(scenario == 5, "Quebec", ""))))))

# Save outputs
## As .RData ##
save(df_threshold_MMS, 
     file = "outputs/DSA/Modified Model Specification/df_threshold_MMS.RData")

#########################
#### Threshold plots ####
#########################

load(file = "outputs/DSA/Modified Model Specification/df_threshold_MMS.RData")

# Prepare data for plotting
df_threshold_MMS <- df_threshold_MMS %>% mutate(perc_increase = as.numeric(v_threshold_rownames_MMS))

df_threshold_qalys_MMS <- df_threshold_MMS %>% select(perc_increase, n_inc_qalys_TOTAL_life, prov)
df_threshold_costs_MMS <- df_threshold_MMS %>% select(perc_increase, n_inc_costs_TOTAL_life, prov)

## Threshold plots ##
# Set colours
# AB, BC, CAN, ON, QC
#v_threshold_colours <- c()

# MMS
# Incremental QALYs
plot_DSA_qalys_MMS_threshold <- ggplot(df_threshold_qalys_MMS, aes(x = perc_increase, y = n_inc_qalys_TOTAL_life, group = prov)) +
  theme_bw() +
  scale_fill_discrete(name = "Provincial Fentanyl Prevalence") +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  geom_line(aes(color = prov), size = 1) +
  scale_colour_brewer(palette = "RdYlBu") +
  #scale_colour_viridis_d() +
  geom_hline(yintercept = 0) +
  #scale_x_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1, suffix = "x"), limits = c(1, 5)) +
  xlab("Increase in BNX Episode Duration") + ylab("Incremental QALYs (BNX vs. MET)") +
  theme(text = element_text(size = 15))

plot_DSA_qalys_MMS_threshold

ggsave(plot_DSA_qalys_MMS_threshold, 
       filename = "Plots/DSA/Threshold SA/DSA-BNX-threshold-qalys-MMS.png", 
       width = 6, height = 6, dpi = 350)

# Incremental Costs
plot_DSA_costs_MMS_threshold <- ggplot(df_threshold_costs_MMS, aes(x = perc_increase, y = n_inc_costs_TOTAL_life, group = prov)) +
  theme_bw() +
  scale_fill_discrete(name = "Provincial Fentanyl Prevalence") +
  theme(legend.position = "none") +
  geom_line(aes(color = prov), size = 1) +
  scale_colour_brewer(palette = "RdYlBu") +
  #scale_color_manual(values = v_threshold_colours) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1, suffix = "x"), limits = c(1, 5)) +
  scale_y_continuous(labels = scales::dollar_format(scale = .001, suffix = "K")) +
  xlab("Increase in BNX Episode Duration") + ylab("Incremental Costs (BNX vs. MET)*") +
  theme(text = element_text(size = 15))

plot_DSA_costs_MMS_threshold

ggsave(plot_DSA_costs_MMS_threshold, 
       filename = "Plots/DSA/Threshold SA/DSA-BNX-threshold-costs-MMS.png", 
       width = 6, height = 6, dpi = 350)

## Combined plot ##
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(plot_DSA_qalys_MMS_threshold)

#plot_DSA_MMS_threshold_comb <- grid.arrange(plot_DSA_costs_MMS_threshold, plot_DSA_qalys_MMS_threshold, ncol = 2, widths = 4:5)

plot_DSA_MMS_threshold_comb <- grid.arrange(arrangeGrob(plot_DSA_qalys_MMS_threshold + theme(legend.position = "none"), 
                                                        plot_DSA_costs_MMS_threshold + theme(legend.position = "none"), nrow = 1),
                                            mylegend, nrow = 2, heights = c(6, 1))

ggsave(plot_DSA_MMS_threshold_comb,
       filename = "Plots/DSA/Threshold SA/DSA-BNX-threshold-combined-MMS.png",
       width = 8, height = 6, dpi = 350)

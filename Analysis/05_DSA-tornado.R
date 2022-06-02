rm(list = ls()) # to clean the workspace

library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(data.table)
library(formattable)
library(tidyr)
library(RColorBrewer)

## Load DSA output files
# Costs
load(file = "outputs/DSA/Modified Model Specification/df_transitions_costs_MMS.RData")
load(file = "outputs/DSA/Modified Model Specification/df_overdose_costs_MMS.RData")
load(file = "outputs/DSA/Modified Model Specification/df_costs_MMS.RData")

# QALYs
load(file = "outputs/DSA/Modified Model Specification/df_transitions_qalys_MMS.RData")
load(file = "outputs/DSA/Modified Model Specification/df_overdose_qalys_MMS.RData")
load(file = "outputs/DSA/Modified Model Specification/df_qalys_MMS.RData")

## Load deterministic outputs for baseline values
load(file = "outputs/ICER/incremental_det_MMS.RData")

## Load PSA outputs for baseline values
load(file = "outputs/PSA/Modified Model Specification/summary_incremental_PSA_MMS.RData")

# Subset by mean
# Deterministic
base_costs <- df_incremental_MMS$n_inc_costs_TOTAL_life[1]
base_qalys <- df_incremental_MMS$n_inc_qalys_TOTAL_life[1]

# PSA
# tbl_df_PSA_mean <- tbl_df_summary_incremental_MMS %>% select(variable, mean)
# base_costs <- tbl_df_PSA_mean %>% deframe %>% getElement("n_inc_costs_TOTAL_life")
# base_qalys <- tbl_df_PSA_mean %>% deframe %>% getElement("n_inc_qalys_TOTAL_life")

## Combine data frames
# Costs
df_DSA_costs_MMS <- rbind(df_transitions_costs_MMS, df_overdose_costs_MMS, df_costs_MMS) %>% 
  mutate(base = base_costs) %>%
  mutate(diff = ifelse(abs(Upper - Lower) > 0, abs(Upper - Lower), abs(base - Upper)))

# QALYs
df_DSA_qalys_MMS <- rbind(df_transitions_qalys_MMS, df_overdose_qalys_MMS, df_qalys_MMS) %>%
  mutate(base = base_qalys) %>%
  mutate(diff = ifelse(abs(Upper - Lower) > 0, abs(Upper - Lower), abs(base - Upper)))

#########################
#### Tornado diagram ####
#########################
### Costs ###
v_order_parameters <- df_DSA_costs_MMS %>% arrange(diff) %>%
  mutate(var_name = factor(x = var_name, levels = var_name)) %>%
  select(var_name) %>% unlist() %>% levels()

v_cost_labels_lower <- c("test1", "test2", "test3", "test4", "test5", "test6", "test7", "test8", "test9", "test10", "test11", "test12", "test13", "test14", "test15", "test16", "test17")
v_cost_labels_upper <- c("test1", "test2", "test3", "test4", "test5", "test6", "test7", "test8", "test9", "test10", "test11", "test12", "test13", "test14", "test15", "test16", "test17")

#eggcounts <- cbind(df_DSA_costs_MMS[, "var_name"], counts)

# width of columns in plot (value between 0 and 1)
width <- 0.6
# get data frame in shape for ggplot and geom_rect
df.2 <- df_DSA_costs_MMS %>% 
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
p_tornado_costs <- ggplot() + 
  geom_rect(data = df.2, 
            aes(ymax = ymax, ymin = ymin, xmax = xmax, xmin = xmin, fill = type)) +
  geom_text(data = eggcounts, aes(y = 1, label = counts), size = 4) +
  theme_bw() + 
  scale_fill_manual(values = c("Upper" = "midnightblue",
                               "Lower" = "slategray2")) +
  theme(axis.title.y = element_blank(), legend.position = 'none', legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = df.2$base) +
  scale_x_continuous(breaks = c(1:length(v_order_parameters)), 
                     labels = v_order_parameters) +
  xlab("Parameter") + ylab("Incremental Cost") +
  coord_flip()

p_tornado_costs

ggsave("Plots/DSA/Modified Model Spec/tornado_costs.png", p_tornado_costs, height = 5, width = 7, dpi = 320)

### QALYs ###
v_order_parameters <- df_DSA_qalys_MMS %>% arrange(diff) %>%
  mutate(var_name = factor(x = var_name, levels = var_name)) %>%
  select(var_name) %>% unlist() %>% levels()

# width of columns in plot (value between 0 and 1)
width <- 0.6
# get data frame in shape for ggplot and geom_rect
df.2 <- df_DSA_qalys_MMS %>% 
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
p_tornado_qalys <- ggplot() + 
  geom_rect(data = df.2, 
            aes(ymax = ymax, ymin = ymin, xmax = xmax, xmin = xmin, fill = type)) +
  theme_bw() + 
  scale_fill_manual(values = c("Upper" = "midnightblue",
                               "Lower" = "slategray2")) +
  theme(axis.title.y = element_blank(), legend.position = 'none',
        legend.title = element_blank()) + 
  geom_hline(yintercept = df.2$base) +
  scale_x_continuous(breaks = c(1:length(v_order_parameters)), 
                     labels = v_order_parameters) +
  xlab("Parameter") + ylab("Incremental QALYs") +
  coord_flip()

p_tornado_qalys

ggsave("Plots/DSA/Modified Model Spec/tornado_qalys.png", p_tornado_qalys, height = 5, width = 7, dpi = 320)

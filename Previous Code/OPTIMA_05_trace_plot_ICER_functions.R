library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(dampack)  # for CEA and calculate ICERs
library(tidyverse)

#### Create plots ####
# Load parameters
l_params_all <- load_all_params(file.init = "data/init_params.csv",
                                file.mort = "data/all_cause_mortality.csv",
                                file.death_hr = "data/death_hr.csv",
                                file.frailty = "data/frailty.csv",
                                file.weibull_scale = "data/weibull_scale.csv",
                                file.weibull_shape = "data/weibull_shape.csv",
                                file.unconditional = "data/unconditional_TP.csv",
                                file.sero = "data/hiv_sero.csv",
                                file.costs = "data/costs.csv",
                                file.crime_costs = "data/crime_costs.csv",
                                file.qalys = "data/qalys.csv")
# Run model
l_out_markov <- markov_model(l_params_all = l_params_all, err_stop = FALSE, verbose = TRUE)

# Prepare data
df_M_agg_trace <- as.data.frame(l_out_markov$m_M_agg_trace)
df_M_agg_trace$month <- as.numeric(rownames(df_M_agg_trace))
df_M_agg_trace_plot <- df_M_agg_trace %>% gather(state, proportion, "Death", "OD", "REL1", "REL", "BUP1", "BUP", "MET1", "MET", "ABS") # health states to plot
df_M_agg_state_time <- df_M_agg_trace %>% gather(state, proportion, "OD", "REL1", "REL", "BUP1", "BUP", "MET1", "MET", "ABS") # alive health states to plot
df_M_agg_state_time <- df_M_agg_state_time %>% 
  group_by(state) %>% 
  summarise_each(funs(sum), proportion) %>%
  mutate(percentage = round((proportion / sum(proportion)) * 100,1))

# Preserve order for plotting
state_order_trace <- factor(df_M_agg_trace_plot$state, levels = c("Death", "OD", "REL1", "REL", "BUP1", "BUP", "MET1", "MET", "ABS"))
state_order_time  <- factor(df_M_agg_state_time$state, levels = c("OD", "REL1", "REL", "BUP1", "BUP", "MET1", "MET", "ABS"))
state_colours_trace <- c("#d9d9d9", "#d53e4f", "#f46d43", "#fdae61", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd")
state_colours_time <- c("#d53e4f", "#f46d43", "#fdae61", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd")

### Markov trace plots ###
main_states_trace_plot <- ggplot(df_M_agg_trace_plot, aes(x = month, y = proportion, fill = state_order_trace)) + 
  theme_bw() +
  geom_area() +
  scale_fill_manual(values = state_colours_trace)
pdf("Plots/Markov Trace/trace_states.pdf")
main_states_trace_plot
dev.off()

### Time spent in health states ###
main_states_time <- ggplot(df_M_agg_state_time, aes(x = state_order_time, y = proportion, fill = state_order_time)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Health State") + ylab("Time") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = state_colours_time) +
  geom_text(aes(label = percentage), hjust = -0.25, size = 3.5) +
  coord_flip(ylim = c(0, 720))
pdf(file = "Plots/Markov Trace/time_states.pdf", width = 8, height = 4)
main_states_time
dev.off()

#######################################
#### State and Transition Outcomes ####
#######################################
### Cost-effectiveness analysis ###
### Vector of total costs for both strategies
v_ted_cost <- c(n_totcost_UC, n_totcost_Tr)
### Vector of effectiveness for both strategies
v_ted_qaly <- c(n_totqaly_UC, n_totqaly_Tr)

### Calculate incremental cost-effectiveness ratios (ICERs)
df_cea <- calculate_icers(cost = v_ted_cost,
                          effect = v_ted_qaly,
                          strategies = v_names_str)
df_cea
### Create CEA table
table_cea <- df_cea
## Format column names
colnames(table_cea)[2:6] <- c("Costs ($)", "QALYs",
                              "Incremental Costs ($)", "Incremental QALYs",
                              "ICER ($/QALY)") # name the columns
## Format rows
table_cea$`Costs ($)` <- comma(round(table_cea$`Costs ($)`, 0))
table_cea$`Incremental Costs ($)` <- comma(round(table_cea$`Incremental Costs ($)`, 0))
table_cea$QALYs <- round(table_cea$QALYs, 3)
table_cea$`Incremental QALYs` <- round(table_cea$`Incremental QALYs`, 3)
table_cea$`ICER ($/QALY)` <- comma(round(table_cea$`ICER ($/QALY)`, 0))
table_cea
### CEA frontier
plot(df_cea) +
  expand_limits(x = 20.8)
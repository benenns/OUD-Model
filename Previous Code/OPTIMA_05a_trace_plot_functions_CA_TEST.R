library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(dampack)  # for CEA and calculate ICERs
library(tidyverse)

# Call model setup functions
source("R/OPTIMA_00_input_parameter_functions_CA_TEST.R")
source("R/OPTIMA_01_model_setup_functions_CA_TEST.R")

# Load parameters
l_params_all_CA_EB <- load_all_params(file.init = "data/CA_EB_Replication/init_params_CA_EB.csv",
                                      file.init_dist = "data/CA_EB_Replication/init_dist_CA_EB.csv",
                                      file.mort = "data/CA_EB_Replication/all_cause_mortality_CA_EB.csv",
                                      file.death_hr = "data/CA_EB_Replication/death_hr_CA_EB.csv",
                                      file.frailty = "data/CA_EB_Replication/frailty_CA_EB.csv",
                                      file.weibull_scale = "data/CA_EB_Replication/weibull_scale_CA_EB.csv",
                                      file.weibull_shape = "data/CA_EB_Replication/weibull_shape_CA_EB.csv",
                                      file.unconditional = "data/CA_EB_Replication/unconditional_TP_CA_EB.csv",
                                      file.sero = "data/CA_EB_Replication/hiv_sero_CA_EB.csv",
                                      file.costs = "data/CA_EB_Replication/costs_CA_EB.csv",
                                      file.crime_costs = "data/CA_EB_Replication/crime_costs_CA_EB.csv",
                                      file.qalys = "data/CA_EB_Replication/qalys_CA_EB.csv")
# Run model
l_out_markov <- markov_model_CA(l_params_all = l_params_all_CA_EB, err_stop = FALSE, verbose = TRUE)

write.csv(l_out_markov$m_M_agg_trace,"C:/Users/Benjamin/Desktop/trace.csv", row.names = TRUE)
#write.csv(l_out_markov$m_M_agg_trace_death,"C:/Users/Benjamin/Desktop/trace_death.csv", row.names = TRUE)
#write.csv(l_out_markov$m_M_agg_trace_sero,"C:/Users/Benjamin/Desktop/trace_sero.csv", row.names = TRUE)
#### Create plots ####
# Prepare data
df_M_agg_trace <- as.data.frame(l_out_markov$m_M_agg_trace)
df_M_agg_trace$month <- as.numeric(rownames(df_M_agg_trace))
df_M_agg_trace_plot <- df_M_agg_trace %>% gather(state, proportion, "Death", "OD", "REL1", "REL", "MET", "ABS") # health states to plot
df_M_agg_state_time <- df_M_agg_trace %>% gather(state, proportion, "OD", "REL1", "REL", "MET", "ABS") # alive health states to plot
df_M_agg_state_time <- df_M_agg_state_time %>% 
  group_by(state) %>% 
  summarise_each(funs(sum), proportion) %>%
  mutate(percentage = round((proportion / sum(proportion)) * 100,1))

# Preserve order for plotting
state_order_trace <- factor(df_M_agg_trace_plot$state, levels = c("Death", "OD", "REL", "REL1", "ABS", "MET"))
state_order_time  <- factor(df_M_agg_state_time$state, levels = c("OD", "REL1", "REL", "MET", "ABS"))
#state_colours_trace <- c("#d9d9d9", "#d53e4f", "#f46d43", "#fdae61", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd") # colour pallette 1
state_colours_trace2 <- c("#d9d9d9", "#b2182b", "#fddbc7", "#d6604d", "#d9f0d3", "#1b7837") # colour pallette 2

#state_colours_time <- c("#d53e4f", "#f46d43", "#fdae61", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd")
state_colours_time2 <- c("#b2182b", "#fddbc7", "#d6604d", "#d9f0d3", "#1b7837") # colour pallette 2

### Markov trace plots ###
main_states_trace_plot <- ggplot(df_M_agg_trace_plot, aes(x = month, y = proportion, fill = state_order_trace)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Time (months)") + ylab("Proportion in state") +
  geom_area() +
  scale_fill_manual(name = "Health States", values = state_colours_trace2)
pdf("Plots/Markov Trace/trace_states_CA_EB.pdf", width = 8, height = 5)
main_states_trace_plot
dev.off()
png("Plots/Markov Trace/trace_states_CA_EB.png", width = 640, height = 480)
main_states_trace_plot
dev.off()

### Time spent in health states ###
main_states_time <- ggplot(df_M_agg_state_time, aes(x = state_order_time, y = proportion, fill = state_order_time)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Health State") + ylab("Time") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = state_colours_time2) +
  geom_text(aes(label = paste0(round((proportion/12),1)," (",percentage,")", "%")), hjust = -0.25, size = 3.5) +
  coord_flip(ylim = c(0, 720))
pdf(file = "Plots/Markov Trace/time_states_CA_EB.pdf", width = 8, height = 4)
main_states_time
dev.off()
png("Plots/Markov Trace/time_states_CA_EB.png", width = 640, height = 480)
main_states_time
dev.off()
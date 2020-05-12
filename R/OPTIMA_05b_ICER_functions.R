library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(dampack)  # for CEA and calculate ICERs
library(tidyverse)

# Load parameters
l_params_all_CA <- load_all_params(file.init = "data/init_params_CA.csv",
                                file.mort = "data/all_cause_mortality_CA.csv",
                                file.death_hr = "data/death_hr_CA.csv",
                                file.frailty = "data/frailty_CA.csv",
                                file.weibull_scale = "data/weibull_scale_CA.csv",
                                file.weibull_shape = "data/weibull_shape_CA.csv",
                                file.unconditional = "data/unconditional_TP_CA.csv",
                                file.sero = "data/hiv_sero_CA.csv",
                                file.costs = "data/costs_CA.csv",
                                file.crime_costs = "data/crime_costs_CA.csv",
                                file.qalys = "data/qalys_CA.csv")
# Run model
l_out_markov <- markov_model(l_params_all = l_params_all_CA, err_stop = FALSE, verbose = TRUE)

#### Cost-effectiveness analysis ####
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
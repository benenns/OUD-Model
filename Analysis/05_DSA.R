library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(tidyverse)
library(xlsx)

# Call model setup functions
# To-do: Move into package eventually
source("R/input_parameter_functions.R")
source("R/generate_psa_parameter_functions.R")
source("R/model_setup_functions.R")
source("R/ICER_functions.R")

# Load parameters
# Calibrated parameter values
load(file = "outputs/imis_output.RData")

# All parameters
# BNX scenario (primary model definition)
l_params_BUP <- load_all_params(file.init = "data/init_params.csv",
                                file.init_dist = "data/init_dist_bup.csv", # Change initial distributions (100% in BUP)
                                file.mort = "data/all_cause_mortality.csv",
                                file.death_hr = "data/death_hr.csv",
                                file.frailty = "data/frailty.csv",
                                file.weibull_scale = "data/weibull_scale.csv",
                                file.weibull_shape = "data/weibull_shape.csv",
                                file.unconditional = "data/unconditional.csv",
                                file.overdose = "data/overdose.csv",
                                file.hiv = "data/hiv_sero.csv",
                                file.hcv = "data/hcv_sero.csv",
                                file.costs = "data/costs.csv",
                                file.crime_costs = "data/crime_costs.csv",
                                file.qalys = "data/qalys.csv")

# Methadone scenario (primary model definition)
l_params_MET <- load_all_params(file.init = "data/init_params.csv",
                                file.init_dist = "data/init_dist_met.csv", # Change initial distributions (100% in MET)
                                file.mort = "data/all_cause_mortality.csv",
                                file.death_hr = "data/death_hr.csv",
                                file.frailty = "data/frailty.csv",
                                file.weibull_scale = "data/weibull_scale.csv",
                                file.weibull_shape = "data/weibull_shape.csv",
                                file.unconditional = "data/unconditional.csv",
                                file.overdose = "data/overdose.csv",
                                file.hiv = "data/hiv_sero.csv",
                                file.hcv = "data/hcv_sero.csv",
                                file.costs = "data/costs.csv",
                                file.crime_costs = "data/crime_costs.csv",
                                file.qalys = "data/qalys.csv")

#### Deterministic sensitivity analysis scenarios ####
## Single parameter OWSA ##
## Overdose ##
# Witnessed OD
# Low
l_outcomes_MET_witness_low <- outcomes(l_params_all = l_params_MET, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_low)
l_outcomes_BUP_witness_low <- outcomes(l_params_all = l_params_BUP, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_low)
ICER_witness_low <- ICER(outcomes_comp = l_outcomes_MET_witness_low, outcomes_int = l_outcomes_BUP_witness_low)
# High
l_outcomes_MET_witness_high <- outcomes(l_params_all = l_params_MET, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_high)
l_outcomes_BUP_witness_high <- outcomes(l_params_all = l_params_BUP, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_witness_high)
ICER_witness_high <- ICER(outcomes_comp = l_outcomes_MET_witness_high, outcomes_int = l_outcomes_BUP_witness_high)

# Fentanyl prevalence
# Low
l_outcomes_MET_fentanyl_low <- outcomes(l_params_all = l_params_MET, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fentanyl_low)
l_outcomes_BUP_fentanyl_low <- outcomes(l_params_all = l_params_BUP, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_fentanyl_low)
ICER_fentanyl_low <- ICER(outcomes_comp = l_outcomes_MET_fentanyl_low, outcomes_int = l_outcomes_BUP_fentanyl_low)

# Naloxone prevalence
l_outcomes_MET_naloxone <- outcomes(l_params_all = l_params_MET, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_naloxone)
l_outcomes_BUP_naloxone <- outcomes(l_params_all = l_params_BUP, v_params_calib = v_calib_post_map, v_params_dsa = v_dsa_naloxone)
ICER_naloxone <- ICER(outcomes_comp = l_outcomes_MET_naloxone, outcomes_int = l_outcomes_BUP_naloxone)

## Cohort characteristics ##
# Starting age

# Male %

## Outcomes ##
# Costs



# QALYs

#### Tornado Diagram ####
# Prepare data for plotting
df_tornado <- data.frame()

# this is throwing some warnings in my computer, but it is reading the data frame correctly
df <- '
Parameter Lower_Bound Upper_Bound UL_Difference
Parameter01 8074 11181 3108 
Parameter02 8177 11007 2831 
Parameter03 8879 10188 1308 
Parameter04 4358 18697 14339 
Parameter05 9073 10087 1013 
Parameter06 12034 7572 4462 
Parameter07 11357 7933 3423 
Parameter08 9769 9202 567 
Parameter09 8833 10403 1570 
Parameter10 13450 4219 9231 
Parameter11 10691 7915 2776 
Parameter12 10036 8792 1244
' %>% read_table2()

# original value of output
base.value <- 9504

# get order of parameters according to size of intervals
# (I use this to define the ordering of the factors which I then use to define the positions in the plot)
order.parameters <- df %>% arrange(UL_Difference) %>%
  mutate(Parameter=factor(x=Parameter, levels=Parameter)) %>%
  select(Parameter) %>% unlist() %>% levels()

# width of columns in plot (value between 0 and 1)
width <- 0.95

# get data frame in shape for ggplot and geom_rect
df.2 <- df %>% 
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key='type', value='output.value', Lower_Bound:Upper_Bound) %>%
  # just reordering columns
  select(Parameter, type, output.value, UL_Difference) %>%
  # create the columns for geom_rect
  mutate(Parameter=factor(Parameter, levels=order.parameters),
         ymin=pmin(output.value, base.value),
         ymax=pmax(output.value, base.value),
         xmin=as.numeric(Parameter)-width/2,
         xmax=as.numeric(Parameter)+width/2)

# create plot
# (use scale_x_continuous to change labels in y axis to name of parameters)
png(width = 960, height = 540)
ggplot() + 
  geom_rect(data = df.2, 
            aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, fill=type)) +
  theme_bw() + 
  theme(axis.title.y=element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  geom_hline(yintercept = base.value) +
  scale_x_continuous(breaks = c(1:length(order.parameters)), 
                     labels = order.parameters) +
  coord_flip()
dev.off()

## Trial health state versus model health state definition ##

## Province-specific analysis ##
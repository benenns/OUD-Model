rm(list = ls()) # to clean the workspace
source("R/ICER_functions.R")
# Run 1-500

load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-16-500runs/ICER_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-16-500runs/incremental_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-16-500runs/outcomes_BUP_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-16-500runs/outcomes_MET_PSA_MMS.RData")

df_ICER_temp1 <- df_ICER_PSA_MMS_comb
df_incremental_temp1 <- df_incremental_PSA_MMS_comb
df_outcomes_BUP_temp1 <- df_outcomes_BUP_PSA_MMS_comb
df_outcomes_MET_temp1 <- df_outcomes_MET_PSA_MMS_comb

# Run 500-2000

load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-19-1500runs/ICER_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-19-1500runs/incremental_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-19-1500runs/outcomes_BUP_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-19-1500runs/outcomes_MET_PSA_MMS.RData")

df_ICER_temp2 <- df_ICER_PSA_MMS_comb
df_incremental_temp2 <- df_incremental_PSA_MMS_comb
df_outcomes_BUP_temp2 <- df_outcomes_BUP_PSA_MMS_comb
df_outcomes_MET_temp2 <- df_outcomes_MET_PSA_MMS_comb

# Combine
df_ICER_PSA_MMS <- rbind(df_ICER_temp1, df_ICER_temp2)
df_incremental_PSA_MMS <- rbind(df_incremental_temp1, df_incremental_temp2)
df_outcomes_BUP_PSA_MMS <- rbind(df_outcomes_BUP_temp1, df_outcomes_BUP_temp2)
df_outcomes_MET_PSA_MMS <- rbind(df_outcomes_MET_temp1, df_outcomes_MET_temp2)


## TEMP CODE ##
m_outcomes_BUP_PSA_MMS <- as.matrix(df_outcomes_BUP_PSA_MMS)
m_outcomes_MET_PSA_MMS <- as.matrix(df_outcomes_MET_PSA_MMS)

m_incremental_PSA_MMS <- m_outcomes_BUP_PSA_MMS - m_outcomes_MET_PSA_MMS
df_incremental_PSA_MMS <- as.data.frame(m_incremental_PSA_MMS)
names(df_incremental_PSA_MMS) <- c("n_inc_costs_TOTAL_6mo", "n_inc_costs_TOTAL_1yr", "n_inc_costs_TOTAL_5yr", "n_inc_costs_TOTAL_10yr", "n_inc_costs_TOTAL_life", 
                                   "n_inc_costs_HEALTH_SECTOR_6mo", "n_inc_costs_HEALTH_SECTOR_1yr", "n_inc_costs_HEALTH_SECTOR_5yr", "n_inc_costs_HEALTH_SECTOR_10yr", "n_inc_costs_HEALTH_SECTOR_life",
                                   "n_inc_costs_CRIMINAL_6mo", "n_inc_costs_CRIMINAL_1yr", "n_inc_costs_CRIMINAL_5yr", "n_inc_costs_CRIMINAL_10yr", "n_inc_costs_CRIMINAL_life",  
                                   "n_inc_costs_TX_6mo", "n_inc_costs_TX_1yr", "n_inc_costs_TX_5yr", "n_inc_costs_TX_10yr", "n_inc_costs_TX_life", 
                                   "n_inc_costs_HRU_6mo", "n_inc_costs_HRU_1yr", "n_inc_costs_HRU_5yr", "n_inc_costs_HRU_10yr", "n_inc_costs_HRU_life",
                                   "n_inc_qalys_TOTAL_6mo", "n_inc_qalys_TOTAL_1yr", "n_inc_qalys_TOTAL_5yr", "n_inc_qalys_TOTAL_10yr", "n_inc_qalys_TOTAL_life")



# Run 2000-5000

load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-24-3000runs/ICER_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-24-3000runs/incremental_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-24-3000runs/outcomes_BUP_PSA_MMS.RData")
load(file = "outputs/PSA/Modified Model Specification/Blockwise-PSA-May-16/May-24-3000runs/outcomes_MET_PSA_MMS.RData")

df_ICER_temp3 <- df_ICER_PSA_MMS_comb
df_incremental_temp3 <- df_incremental_PSA_MMS_comb
df_outcomes_BUP_temp3 <- df_outcomes_BUP_PSA_MMS_comb
df_outcomes_MET_temp3 <- df_outcomes_MET_PSA_MMS_comb

# Combine
df_ICER_PSA_MMS <- rbind(df_ICER_PSA_MMS, df_ICER_temp3)
df_incremental_PSA_MMS <- rbind(df_incremental_PSA_MMS, df_incremental_temp3)
df_outcomes_BUP_PSA_MMS <- rbind(df_outcomes_BUP_PSA_MMS, df_outcomes_BUP_temp3)
df_outcomes_MET_PSA_MMS <- rbind(df_outcomes_MET_PSA_MMS, df_outcomes_MET_temp3)

### Output results
## As .RData
save(df_outcomes_MET_PSA_MMS, 
     file = "outputs/PSA/Modified Model Specification/outcomes_MET_PSA_MMS.RData")
save(df_outcomes_BUP_PSA_MMS, 
     file = "outputs/PSA/Modified Model Specification/outcomes_BUP_PSA_MMS.RData")
save(df_ICER_PSA_MMS, 
     file = "outputs/PSA/Modified Model Specification/ICER_PSA_MMS.RData")
save(df_incremental_PSA_MMS,
     file = "outputs/PSA/Modified Model Specification/incremental_PSA_MMS.RData")

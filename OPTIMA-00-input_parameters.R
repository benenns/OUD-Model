##################################################################################
## OUD Cohort Model
## Deterministic - Set parameter values for: "OPTIMA_01_model_setup_functions.R"
##################################################################################
###############################
### Initial characteristics ###
###############################

v_initial_BUP <- read.csv("data/01_initial_BUP.csv")
v_initial_MET <- read.csv("data/01_initial_MET.csv")

n_age_init <- 35 # age at baseline
n_age_max <- 95 # maximum age of follow up

n_t <- (n_age_max - n_age_init) * 12 # modeling time horizon in months
p_discount <- 0.03

p_male_BUP <- 0.35
p_male_MET <- 0.40

p_HIV_POS <- 0.05 # % of HIV-positive individuals
p_HIV_ART <- 0.75 # % of HIV-positive on-ART

#########################
### Survival analysis ###
#########################
m_frailty <- array(0.4, dim = c(n_EP, n_BASE, n_INJECT),
                   dimnames = list(EP, BASE, INJECT))
m_weibull_scale <- array(0.7, dim = c(n_BASE, n_INJECT),
                         dimnames = list(BASE, INJECT))
m_weibull_shape <- array(0.7, dim = c(n_BASE, n_INJECT),
                         dimnames = list(BASE, INJECT))


########################## 
### HIV Seroconversion ###
##########################
p_sero <- read.csv("data/01_hiv_sero.csv", header = TRUE)
p_sero

# Only applied to injection
p_sero_BUP1_NI <- p_sero[1, 2]
p_sero_BUP_NI  <- p_sero[2, 2]
p_sero_MET1_NI <- p_sero[3, 2]
p_sero_MET_NI  <- p_sero[4, 2]
p_sero_REL1_NI <- p_sero[5, 2]
p_sero_REL_NI  <- p_sero[6, 2]
p_sero_OD_NI   <- p_sero[7, 2]
p_sero_ABS_NI  <- p_sero[8, 2]

###############################
### Age-dependent mortality ###
###############################
# Baseline mortality by age
lt_can_2018 <- read.csv("data/01_all_cause_mortality.csv")
#v_mortality_age <- lt_can_2018[which(lt_can_2018$Age >= n_age_init & lt_can_2018$Age < n_age_max), ] 

v_r_mort_by_age <- lt_can_2018 %>%
  #subset(Age >= n_age_init & Age < n_age_max) %>%
  #filter(Age >= n_age_init & Age < n_age_max) %>%
  select(Total) %>%
  as.matrix()

# Hazard ratios for death probability
hr_s <- read.csv("data/01_death_hr.csv", header = TRUE)
#hr_s
#hr_BUP1 <- hr_s["BUP1", "HR"] figure out why this doesn't work
hr_BUP1_NI <- hr_s[1, 2]
hr_BUP_NI  <- hr_s[2, 2]
hr_MET1_NI <- hr_s[3, 2]
hr_MET_NI  <- hr_s[4, 2]
hr_REL1_NI <- hr_s[5, 2]
hr_REL_NI  <- hr_s[6, 2]
hr_OD_NI   <- hr_s[7, 2]
hr_ABS_NI  <- hr_s[8, 2]

hr_BUP1_INJ <- hr_s[9, 2]
hr_BUP_INJ  <- hr_s[10, 2]
hr_MET1_INJ <- hr_s[11, 2]
hr_MET_INJ  <- hr_s[12, 2]
hr_REL1_INJ <- hr_s[13, 2]
hr_REL_INJ  <- hr_s[14, 2]
hr_OD_INJ   <- hr_s[15, 2]
hr_ABS_INJ  <- hr_s[16, 2]

hr_ABS_HIV  <- hr_s[17, 2]

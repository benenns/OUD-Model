#' Load mortality data
#'
#' \code{load_mort_data} is used to load age-specific mortality from .csv file 
#' into vector.
#'
#' @param file String with the location and name of the file with mortality 
#' data.
#' @return 
#' A vector with mortality by age.
#' @export
load_mort_params <- function(file.mort = NULL, n_male){
  df_lt_can_2018 <- read.csv(file = file.mort)
  v_r_mort_by_age_male <- df_lt_can_2018 %>% select(Male) %>% as.matrix()
  v_r_mort_by_age_female <- df_lt_can_2018 %>% select(Female) %>% as.matrix()
  v_r_mort_by_age <- (v_r_mort_by_age_male * n_male) + (v_r_mort_by_age_female * (1 - n_male)) # weighted mortality
  return(v_r_mort_by_age)
}

#' Load all parameters
#'
#' \code{load_all_params} loads all parameters for the decision model from multiple sources and creates a list.
#'
#' @param file.init String with the location and name of the file with initial set of parameters
#' @param file.init_dist String with the location and name of the file with initial distributions
#' @param file.mort String with the location and name of the file with mortality data
#' @param file.death_hr String with the location and name of death hazard ratios
#' @param file.frailty String with the location and name of the file with frailty estimates
#' @param file.weibull_scale String with the location and name of the file with weibull scale
#' @param file.weibull_shape String with the location and name of the file with weibull shape
#' @param file.unconditional String with the location and name of the file with empirical destination states
#' @param file.overdose String with the location and name of the file with overdose/fentanyl-related parameters
#' @param file.hiv String with the location and name of the file with HIV seroconversion probabilities
#' @param file.hcv String with the location and name of the file with HCV seroconversion probabilities
#' @param file.costs String with the location and name of the file with costs (excluding crime costs)
#' @param file.crime_costs String with the location and name of the file with age-specific crime costs
#' @param file.qalys String with the location and name of the file with HRQoL weights
#' 
#' @return 
#' A list of all parameters used for the decision model.
#' @export
load_all_params <- function(file.init = NULL,
                            file.init_dist = NULL,
                            file.mort = NULL,
                            file.death_hr = NULL,
                            file.frailty = NULL,
                            file.weibull_scale = NULL,
                            file.weibull_shape = NULL,
                            file.unconditional = NULL,
                            file.overdose = NULL,
                            file.hiv = NULL,
                            file.hcv = NULL,
                            file.costs = NULL,
                            file.crime_costs = NULL,
                            file.qalys = NULL){ # User defined
    
  #Load files of all baseline model parameters
  df_init_params <- read.csv(file = file.init, row.names = 1, header = TRUE) # Initial parameter values
  df_init_dist <- read.csv(file = file.init_dist, row.names = 1, header = TRUE) # Initial parameter values
  df_death_hr <- read.csv(file = file.death_hr, row.names = 1, header = TRUE) # Mortality hazard ratios
  df_frailty <- read.csv(file = file.frailty, row.names = 1, header = TRUE) # Episode frailty params
  df_weibull_scale <- read.csv(file = file.weibull_scale, row.names = 1, header = TRUE) # Weibull scale params
  df_weibull_shape <- read.csv(file = file.weibull_shape, row.names = 1, header = TRUE) # Weibull shape params
  df_UP <- read.csv(file = file.unconditional, row.names = 1, header = TRUE) # Unconditional transition probs
  df_overdose <- read.csv(file = file.overdose, row.names = 1, header = TRUE) # Overdose-fentanyl parameters
  df_hiv <- read.csv(file = file.hiv, row.names = 1, header = TRUE) # HIV seroconversion probs
  df_hcv <- read.csv(file = file.hcv, row.names = 1, header = TRUE) # HCV seroconversion probs
  df_costs <- read.csv(file = file.costs, row.names = 1, header = TRUE) # All costs excluding crime
  df_crime_costs <- read.csv(file = file.crime_costs, header = TRUE) # Age-dependent crime costs
  df_qalys <- read.csv(file = file.qalys, row.names = 1, header = TRUE) # QALYs
  
  l_params_all <- list(
    # Initial parameters
    n_age_init = df_init_params["pe", "age_init"], # age at baseline
    n_age_max = df_init_params["pe", "age_max"], # maximum age of follow up
    n_per = df_init_params["pe", "period_yr"], # periods per year (12-months/52-weeks)
    #n_t = (df_init_params["pe", "age_max"] - df_init_params["pe", "age_init"]) * 52, # modeling time horizon in weeks
    n_discount = df_init_params["pe", "discount"], # discount rate
    n_male = df_init_params["pe", "male_prop"], # % male
    n_INJ = df_init_params["pe", "inj_prop"], # % injection
    n_HIV = df_init_params["pe", "hiv_prop"], # % of HIV-positive individuals
    n_HCV = df_init_params["pe", "hcv_prop"], # % of HCV-positive individuals
    n_COI = df_init_params["pe", "coi_prop"], # % of co-infected individuals
    n_HIV_ART = df_init_params["pe", "art_prop"], # % of HIV-positive on-ART (used to calculate costs)
    n_HCV_DAA = df_init_params["pe", "daa_prop"], # % of HIV-positive on-ART (used to calculate costs)
    
    # Initial state distribution
    v_init_dist = as.vector(df_init_dist["pe", ]),
    
    # Mortality
    v_r_mort_by_age = load_mort_params(file = file.mort, n_male = df_init_params["pe", "male_prop"]), # vector of age-specific mortality
    
    # Hazard ratios for death probability
    hr_BUP_NI  = df_death_hr["pe", "BUP_NI"],
    hr_MET_NI  = df_death_hr["pe", "MET_NI"],
    hr_REL_NI  = df_death_hr["pe", "REL_NI"],
    hr_ODN_NI   = df_death_hr["pe", "ODN_NI"],
    hr_ODF_NI   = df_death_hr["pe", "ODF_NI"],
    hr_ABS_NI  = df_death_hr["pe", "ABS_NI"],
    hr_HIV_NI  = df_death_hr["pe", "HIV_NI"],
    hr_HCV_NI  = df_death_hr["pe", "HCV_NI"],
    hr_COI_NI  = df_death_hr["pe", "COI_NI"],
    
    hr_BUP_INJ  = df_death_hr["pe", "BUP_INJ"],
    hr_MET_INJ  = df_death_hr["pe", "MET_INJ"],
    hr_REL_INJ  = df_death_hr["pe", "REL_INJ"],
    hr_ODN_INJ   = df_death_hr["pe", "ODN_INJ"],
    hr_ODF_INJ   = df_death_hr["pe", "ODF_INJ"],
    hr_ABS_INJ  = df_death_hr["pe", "ABS_INJ"],
    hr_HIV_INJ  = df_death_hr["pe", "HIV_INJ"],
    hr_HCV_INJ  = df_death_hr["pe", "HCV_INJ"],
    hr_COI_INJ  = df_death_hr["pe", "COI_INJ"],
    
    # Survival analysis
    p_frailty_BUP_NI_1 = 1,
    p_frailty_BUP_NI_2 = df_frailty["pe", "BUP_NI_2"],
    p_frailty_BUP_NI_3 = df_frailty["pe", "BUP_NI_3"],
    p_frailty_MET_NI_1 = 1,
    p_frailty_MET_NI_2 = df_frailty["pe", "MET_NI_2"],
    p_frailty_MET_NI_3 = df_frailty["pe", "MET_NI_3"],
    p_frailty_ABS_NI_1 = 1,
    p_frailty_ABS_NI_2 = df_frailty["pe", "ABS_NI_2"],
    p_frailty_ABS_NI_3 = df_frailty["pe", "ABS_NI_3"],
    p_frailty_REL_NI_1 = 1,
    p_frailty_REL_NI_2 = df_frailty["pe", "REL_NI_2"],
    p_frailty_REL_NI_3 = df_frailty["pe", "REL_NI_3"],
    
    p_frailty_BUP_INJ_1 = 1,
    p_frailty_BUP_INJ_2 = df_frailty["pe", "BUP_INJ_2"],
    p_frailty_BUP_INJ_3 = df_frailty["pe", "BUP_INJ_3"],
    p_frailty_MET_INJ_1 = 1,
    p_frailty_MET_INJ_2 = df_frailty["pe", "MET_INJ_2"],
    p_frailty_MET_INJ_3 = df_frailty["pe", "MET_INJ_3"],
    p_frailty_ABS_INJ_1 = 1,
    p_frailty_ABS_INJ_2 = df_frailty["pe", "ABS_INJ_2"],
    p_frailty_ABS_INJ_3 = df_frailty["pe", "ABS_INJ_3"],
    p_frailty_REL_INJ_1 = 1,
    p_frailty_REL_INJ_2 = df_frailty["pe", "REL_INJ_2"],
    p_frailty_REL_INJ_3 = df_frailty["pe", "REL_INJ_3"],
    
    # Load weibull parameters
    # Weibull scale
    p_weibull_scale_BUP_NI = df_weibull_scale["pe", "BUP_NI"],
    p_weibull_scale_MET_NI = df_weibull_scale["pe", "MET_NI"],
    p_weibull_scale_REL_NI = df_weibull_scale["pe", "REL_NI"],
    p_weibull_scale_ABS_NI = df_weibull_scale["pe", "ABS_NI"],
    p_weibull_scale_BUP_INJ = df_weibull_scale["pe", "BUP_INJ"],
    p_weibull_scale_MET_INJ = df_weibull_scale["pe", "MET_INJ"],
    p_weibull_scale_REL_INJ = df_weibull_scale["pe", "REL_INJ"],
    p_weibull_scale_ABS_INJ = df_weibull_scale["pe", "ABS_INJ"],

    # Weibull shape
    p_weibull_shape_BUP_NI = df_weibull_shape["pe", "BUP_NI"],
    p_weibull_shape_MET_NI = df_weibull_shape["pe", "MET_NI"],
    p_weibull_shape_REL_NI = df_weibull_shape["pe", "REL_NI"],
    p_weibull_shape_ABS_NI = df_weibull_shape["pe", "ABS_NI"],
    p_weibull_shape_BUP_INJ = df_weibull_shape["pe", "BUP_INJ"],
    p_weibull_shape_MET_INJ = df_weibull_shape["pe", "MET_INJ"],
    p_weibull_shape_REL_INJ = df_weibull_shape["pe", "REL_INJ"],
    p_weibull_shape_ABS_INJ = df_weibull_shape["pe", "ABS_INJ"],

    # Unconditional transition probabilities
    # Non-Injection
    # From BUP1 & BUP
    p_BUP_MET_NI  = df_UP["BUP_NI", "MET_NI"],
    p_BUP_ABS_NI   = df_UP["BUP_NI", "ABS_NI"],
    p_BUP_REL_NI  = df_UP["BUP_NI", "REL_NI"],
    # From MET1 & MET
    p_MET_BUP_NI  = df_UP["MET_NI", "BUP_NI"],
    p_MET_ABS_NI   = df_UP["MET_NI", "ABS_NI"],
    p_MET_REL_NI  = df_UP["MET_NI", "REL_NI"],
    # From ABS
    p_ABS_REL_NI = df_UP["ABS_NI", "REL_NI"],
    # From REL1 & REL
    p_REL_MET_NI  = df_UP["REL_NI", "MET_NI"],
    p_REL_BUP_NI  = df_UP["REL_NI", "BUP_NI"],
    p_REL_ABS_NI   = df_UP["REL_NI", "ABS_NI"],
    # From OD
    p_ODN_MET_NI  = df_UP["ODN_NI", "MET_NI"],
    p_ODN_BUP_NI  = df_UP["ODN_NI", "BUP_NI"],
    p_ODN_ABS_NI  = df_UP["ODN_NI", "ABS_NI"],
    p_ODN_REL_NI  = df_UP["ODN_NI", "REL_NI"],

    # Inj
    # From BUP1 & BUP
    p_BUP_MET_INJ  = df_UP["BUP_INJ", "MET_INJ"],
    p_BUP_ABS_INJ   = df_UP["BUP_INJ", "ABS_INJ"],
    p_BUP_REL_INJ  = df_UP["BUP_INJ", "REL_INJ"],
    # From MET
    p_MET_BUP_INJ  = df_UP["MET_INJ", "BUP_INJ"],
    p_MET_ABS_INJ   = df_UP["MET_INJ", "ABS_INJ"],
    p_MET_REL_INJ  = df_UP["MET_INJ", "REL_INJ"],
    # From ABS
    p_ABS_REL_INJ = df_UP["ABS_INJ", "REL_INJ"],
    # From REL
    p_REL_MET_INJ  = df_UP["REL_INJ", "MET_INJ"],
    p_REL_BUP_INJ  = df_UP["REL_INJ", "BUP_INJ"],
    p_REL_ABS_INJ   = df_UP["REL_INJ", "ABS_INJ"],
    # From OD
    p_ODN_MET_INJ  = df_UP["ODN_INJ", "MET_INJ"],
    p_ODN_BUP_INJ  = df_UP["ODN_INJ", "BUP_INJ"],
    p_ODN_ABS_INJ  = df_UP["ODN_INJ", "ABS_INJ"],
    p_ODN_REL_INJ  = df_UP["ODN_INJ", "REL_INJ"],
    
    #### Overdose ####
    # Includes additional calibration related parameters
    
    # Non-injection
    # Mean
    p_BUP_OD_NI = df_overdose["pe", "BUP_OD_NI"],
    p_MET_OD_NI = df_overdose["pe", "MET_OD_NI"],
    p_ABS_OD_NI = df_overdose["pe", "ABS_OD_NI"],
    p_REL_OD_NI = df_overdose["pe", "REL_OD_NI"],
    # Lower bound
    p_BUP_OD_NI_lb = df_overdose["low", "BUP_OD_NI"], 
    p_MET_OD_NI_lb = df_overdose["low", "MET_OD_NI"],
    p_REL_OD_NI_lb = df_overdose["low", "REL_OD_NI"],
    p_ABS_OD_NI_lb = df_overdose["low", "ABS_OD_NI"],
    # Upper bound
    p_BUP_OD_NI_ub = df_overdose["high", "BUP_OD_NI"],
    p_MET_OD_NI_ub = df_overdose["high", "MET_OD_NI"],
    p_REL_OD_NI_ub = df_overdose["high", "REL_OD_NI"],
    p_ABS_OD_NI_ub = df_overdose["high", "ABS_OD_NI"],
    
    # Injection
    # Mean
    p_BUP_OD_INJ = df_overdose["pe", "BUP_OD_INJ"],
    p_MET_OD_INJ = df_overdose["pe", "MET_OD_INJ"],
    p_ABS_OD_INJ = df_overdose["pe", "ABS_OD_INJ"],
    p_REL_OD_INJ = df_overdose["pe", "REL_OD_INJ"],
    # Lower bound
    p_BUP_OD_INJ_lb = df_overdose["low", "BUP_OD_INJ"],
    p_MET_OD_INJ_lb = df_overdose["low", "MET_OD_INJ"],
    p_REL_OD_INJ_lb = df_overdose["low", "REL_OD_INJ"],
    p_ABS_OD_INJ_lb = df_overdose["low", "ABS_OD_INJ"],
    # Upper bound
    p_BUP_OD_INJ_ub = df_overdose["high", "BUP_OD_INJ"],
    p_MET_OD_INJ_ub = df_overdose["high", "MET_OD_INJ"],
    p_REL_OD_INJ_ub = df_overdose["high", "REL_OD_INJ"],
    p_ABS_OD_INJ_ub = df_overdose["high", "ABS_OD_INJ"],
    
    # Overdose transition multipliers for first 4 weeks of treatment and relapse
    p_BUP_OD_mult = df_overdose["pe", "BUP_OD_mult"],
    p_MET_OD_mult = df_overdose["pe", "MET_OD_mult"],
    p_REL_OD_mult = df_overdose["pe", "REL_OD_mult"],
    
    # Fatal overdose
    p_fatal_OD = df_overdose["pe", "fatal_OD_prob"],
    
    # Fentanyl
    p_fent_OD = df_overdose["pe", "fent_OD_prob"],
    p_fent_exp = df_overdose["pe", "fent_exp_prob"], # prob of fentanyl exposure
    
    # Naloxone
    p_witness = df_overdose["pe", "witness_prob"],
    p_NX_used = df_overdose["pe", "NX_prob"],
    p_NX_success = df_overdose["pe", "NX_success_prob"],
    
    #### Seroconversion ####
    
    # HIV Seroconversion
    # From negative
    # Non-injection
    p_HIV_BUP_NI  = df_hiv["pe", "HIV_BUP_NI"],
    p_HIV_MET_NI  = df_hiv["pe", "HIV_MET_NI"],
    p_HIV_REL_NI  = df_hiv["pe", "HIV_REL_NI"],
    p_HIV_ODN_NI  = df_hiv["pe", "HIV_REL_NI"],
    p_HIV_ABS_NI  = df_hiv["pe", "HIV_ABS_NI"],
    # Injection
    p_HIV_BUP_INJ  = df_hiv["pe", "HIV_BUP_INJ"],
    p_HIV_MET_INJ  = df_hiv["pe", "HIV_MET_INJ"],
    p_HIV_REL_INJ  = df_hiv["pe", "HIV_REL_INJ"],
    p_HIV_ODN_INJ  = df_hiv["pe", "HIV_REL_INJ"],
    p_HIV_ABS_INJ  = df_hiv["pe", "HIV_ABS_INJ"],
    
    # Co-infection conditional on HCV
    # Non-injection
    p_HCV_HIV_BUP_NI  = df_hiv["pe", "COI_BUP_NI"],
    p_HCV_HIV_MET_NI  = df_hiv["pe", "COI_MET_NI"],
    p_HCV_HIV_REL_NI  = df_hiv["pe", "COI_REL_NI"],
    p_HCV_HIV_ODN_NI  = df_hiv["pe", "COI_REL_NI"],
    p_HCV_HIV_ABS_NI  = df_hiv["pe", "COI_ABS_NI"],
    # Injection
    p_HCV_HIV_BUP_INJ  = df_hiv["pe", "COI_BUP_INJ"],
    p_HCV_HIV_MET_INJ  = df_hiv["pe", "COI_MET_INJ"],
    p_HCV_HIV_REL_INJ  = df_hiv["pe", "COI_REL_INJ"],
    p_HCV_HIV_ODN_INJ  = df_hiv["pe", "COI_REL_INJ"],
    p_HCV_HIV_ABS_INJ  = df_hiv["pe", "COI_ABS_INJ"],
    
    # HCV Seroconversion
    # From negative
    # Non-injection
    p_HCV_BUP_NI  = df_hcv["pe", "HCV_BUP_NI"],
    p_HCV_MET_NI  = df_hcv["pe", "HCV_MET_NI"],
    p_HCV_REL_NI  = df_hcv["pe", "HCV_REL_NI"],
    p_HCV_ODN_NI  = df_hcv["pe", "HCV_REL_NI"],
    p_HCV_ABS_NI  = df_hcv["pe", "HCV_ABS_NI"],
    # Injection
    p_HCV_BUP_INJ  = df_hcv["pe", "HCV_BUP_INJ"],
    p_HCV_MET_INJ  = df_hcv["pe", "HCV_MET_INJ"],
    p_HCV_REL_INJ  = df_hcv["pe", "HCV_REL_INJ"],
    p_HCV_ODN_INJ  = df_hcv["pe", "HCV_REL_INJ"],
    p_HCV_ABS_INJ  = df_hcv["pe", "HCV_ABS_INJ"],
    
    # Co-infection conditional on HIV
    # Non-injection
    p_HIV_HCV_BUP_NI  = df_hcv["pe", "COI_BUP_NI"],
    p_HIV_HCV_MET_NI  = df_hcv["pe", "COI_MET_NI"],
    p_HIV_HCV_REL_NI  = df_hcv["pe", "COI_REL_NI"],
    p_HIV_HCV_ODN_NI  = df_hcv["pe", "COI_REL_NI"],
    p_HIV_HCV_ABS_NI  = df_hcv["pe", "COI_ABS_NI"],
    # Injection
    p_HIV_HCV_BUP_INJ  = df_hcv["pe", "COI_BUP_INJ"],
    p_HIV_HCV_MET_INJ  = df_hcv["pe", "COI_MET_INJ"],
    p_HIV_HCV_REL_INJ  = df_hcv["pe", "COI_REL_INJ"],
    p_HIV_HCV_ODN_INJ  = df_hcv["pe", "COI_REL_INJ"],
    p_HIV_HCV_ABS_INJ  = df_hcv["pe", "COI_ABS_INJ"],

    #### Costs ####
    
    # Treatment Costs
    c_BUP_TX  = df_costs["pe", "BUP_TX"],
    c_MET_TX  = df_costs["pe", "MET_TX"],
    c_ODN_TX  = df_costs["pe", "OD_TX"],
    
    # HRU Costs
    # Modify if age-specific
    c_BUP_NI_HRU = df_costs["pe", "BUP_NI_HRU"],
    c_MET_NI_HRU = df_costs["pe", "MET_NI_HRU"],
    c_REL_NI_HRU = df_costs["pe", "REL_NI_HRU"],
    c_ODN_NI_HRU = df_costs["pe", "OD_NI_HRU"],
    c_ABS_NI_HRU = df_costs["pe", "ABS_NI_HRU"],
    c_BUP_INJ_HRU = df_costs["pe", "BUP_INJ_HRU"], 
    c_MET_INJ_HRU = df_costs["pe", "MET_INJ_HRU"],
    c_REL_INJ_HRU = df_costs["pe", "REL_INJ_HRU"],
    c_ODN_INJ_HRU = df_costs["pe", "OD_INJ_HRU"],
    c_ABS_INJ_HRU = df_costs["pe", "ABS_INJ_HRU"], 

    # HIV Costs
    c_HIV_HRU = df_costs["pe", "HIV_HRU"],
    c_HIV_ART = df_costs["pe", "HIV_ART"],
    
    # Crime Costs
    df_crime_costs = subset(df_crime_costs, type=="pe"),
    # Age-specific
    v_c_BUP_NI_crime = df_crime_costs %>% select(BUP_NI) %>% as.matrix(),
    v_c_MET_NI_crime = df_crime_costs %>% select(MET_NI) %>% as.matrix(),
    v_c_REL_NI_crime = df_crime_costs %>% select(REL_NI) %>% as.matrix(),
    v_c_ODN_NI_crime = df_crime_costs %>% select(ODN_NI) %>% as.matrix(),
    v_c_ABS_NI_crime = df_crime_costs %>% select(ABS_NI) %>% as.matrix(),
    
    v_c_BUP_INJ_crime = df_crime_costs %>% select(BUP_INJ) %>% as.matrix(),
    v_c_MET_INJ_crime = df_crime_costs %>% select(MET_INJ) %>% as.matrix(),
    v_c_REL_INJ_crime = df_crime_costs %>% select(REL_INJ) %>% as.matrix(),
    v_c_ODN_INJ_crime = df_crime_costs %>% select(ODN_INJ) %>% as.matrix(),
    v_c_ABS_INJ_crime = df_crime_costs %>% select(ABS_INJ) %>% as.matrix(), 

    #### QALYs ####
    
    u_BUP_NI_NEG = df_qalys["pe", "BUP_NI_NEG"],
    u_MET_NI_NEG = df_qalys["pe", "MET_NI_NEG"],
    u_REL_NI_NEG = df_qalys["pe", "REL_NI_NEG"],
    u_ODN_NI_NEG = df_qalys["pe", "ODN_NI_NEG"],
    u_ABS_NI_NEG = df_qalys["pe", "ABS_NI_NEG"],
    
    u_BUP_INJ_NEG = df_qalys["pe", "BUP_INJ_NEG"],
    u_MET_INJ_NEG = df_qalys["pe", "MET_INJ_NEG"],
    u_REL_INJ_NEG = df_qalys["pe", "REL_INJ_NEG"],
    u_ODN_INJ_NEG  = df_qalys["pe", "ODN_INJ_NEG"],
    u_ABS_INJ_NEG = df_qalys["pe", "ABS_INJ_NEG"],
    
    u_BUP_NI_POS = df_qalys["pe", "BUP_NI_POS"],
    u_MET_NI_POS = df_qalys["pe", "MET_NI_POS"],
    u_REL_NI_POS = df_qalys["pe", "REL_NI_POS"],
    u_ODN_NI_POS = df_qalys["pe", "ODN_NI_POS"],
    u_ABS_NI_POS = df_qalys["pe", "ABS_NI_POS"],
    
    u_BUP_INJ_POS = df_qalys["pe", "BUP_INJ_POS"],
    u_MET_INJ_POS = df_qalys["pe", "MET_INJ_POS"],
    u_REL_INJ_POS = df_qalys["pe", "REL_INJ_POS"],
    u_ODN_INJ_POS = df_qalys["pe", "ODN_INJ_POS"],
    u_ABS_INJ_POS = df_qalys["pe", "ABS_INJ_POS"],
    
    u_HIV_mult = df_qalys["pe", "HIV_mult"] # Modify if using state-HIV-specific QALYs, also add HCV

    ) # Close list
  return(l_params_all) # Return full parameter list
}

#' Update parameters
#'
#' \code{update_param_list} is used to update list of all parameters with new 
#' values for specific parameters.
#'
#' @param l_params_all List with all parameters of decision model
#' @param params_updated Parameters for which values need to be updated
#' @return 
#' A modifed list with all parameters updated.
#' @export

#update_param_list <- function(l_params_all, params_updated){
#  l_params_all <- modifyList(l_params_all, params_updated) #update values
#  return(l_params_all)
#}
update_param_list <- function(l_params_all, params_updated){
  
  if (typeof(params_updated)!="list"){
    params_updated <- split(unname(params_updated),names(params_updated)) #convert the named vector to a list
  }
  l_params_all <- modifyList(l_params_all, params_updated) #update the values
  return(l_params_all)
}
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
load_mort_params <- function(file.mort = NULL, type = NULL){
  df_lt_can_2018 <- read.csv(file = file.mort)
  v_r_mort_by_age <- df_lt_can_2018 %>%
    select(type) %>%
    as.matrix()
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
#' @param file.sero String with the location and name of the file with seroconversion probabilities
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
                            file.sero = NULL,
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
  df_sero <- read.csv(file = file.sero, row.names = 1, header = TRUE) # Seroconversion probs
  df_costs <- read.csv(file = file.costs, row.names = 1, header = TRUE) # All costs excluding crime
  df_crime_costs <- read.csv(file = file.crime_costs, header = TRUE) # Age-dependent crime costs
  df_qalys <- read.csv(file = file.qalys, row.names = 1, header = TRUE) # QALYs
  
  l_params_all <- list(
    # Initial parameters
    n_age_init = df_init_params["pe", "age_init"], # age at baseline
    n_age_max = df_init_params["pe", "age_max"], # maximum age of follow up
    p_discount = df_init_params["pe", "discount"], # discount rate
    p_male = df_init_params["pe", "male_prop"], # % male
    p_INJ = df_init_params["pe", "inj_prop"], # % injection
    p_HIV_POS = df_init_params["pe", "hiv_prop"], # % of HIV-positive individuals
    p_HIV_ART = df_init_params["pe", "art_prop"], # % of HIV-positive on-ART (used to calculate costs)
    
    # Initial state distribution
    v_init_dist = as.vector(df_init_dist["pe", ]),
    
    # Mortality
    v_r_mort_by_age_NI = load_mort_params(file = file.mort, type = "Total_NI"), # vector of age-specific mortality
    v_r_mort_by_age_INJ = load_mort_params(file = file.mort, type = "Total_INJ"), # vector of age-specific mortality
    
    # Hazard ratios for death probability

    hr_MET_NI  = df_death_hr["pe", "MET_NI"],
    hr_REL1_NI = df_death_hr["pe", "REL1_NI"],
    hr_REL_NI  = df_death_hr["pe", "REL_NI"],
    hr_OD_NI   = df_death_hr["pe", "OD_NI"],
    hr_ABS_NI  = df_death_hr["pe", "ABS_NI"],
    hr_HIV_NI  = df_death_hr["pe", "HIV_NI"],
    hr_HIV_OD_NI = df_death_hr["pe", "HIV_OD_NI"],
    

    hr_MET_INJ  = df_death_hr["pe", "MET_INJ"],
    hr_REL1_INJ = df_death_hr["pe", "REL1_INJ"],
    hr_REL_INJ  = df_death_hr["pe", "REL_INJ"],
    hr_OD_INJ   = df_death_hr["pe", "OD_INJ"],
    hr_ABS_INJ  = df_death_hr["pe", "ABS_INJ"],
    hr_HIV_INJ  = df_death_hr["pe", "HIV_INJ"],
    hr_HIV_OD_INJ = df_death_hr["pe", "HIV_OD_INJ"],

    # Survival analysis

    p_frailty_MET_NI_1 = 1,
    p_frailty_MET_NI_2 = df_frailty["pe", "MET_NI_2"],
    p_frailty_MET_NI_3 = df_frailty["pe", "MET_NI_3"],
    p_frailty_ABS_NI_1 = 1,
    p_frailty_ABS_NI_2 = df_frailty["pe", "ABS_NI_2"],
    p_frailty_ABS_NI_3 = df_frailty["pe", "ABS_NI_3"],
    p_frailty_REL_NI_1 = 1,
    p_frailty_REL_NI_2 = df_frailty["pe", "REL_NI_2"],
    p_frailty_REL_NI_3 = df_frailty["pe", "REL_NI_3"],
    p_frailty_OD_NI_1  = 1,
    p_frailty_OD_NI_2  = df_frailty["pe", "OD_NI_2"],
    p_frailty_OD_NI_3  = df_frailty["pe", "OD_NI_3"],

    p_frailty_MET_INJ_1 = 1,
    p_frailty_MET_INJ_2 = df_frailty["pe", "MET_INJ_2"],
    p_frailty_MET_INJ_3 = df_frailty["pe", "MET_INJ_3"],
    p_frailty_ABS_INJ_1 = 1,
    p_frailty_ABS_INJ_2 = df_frailty["pe", "ABS_INJ_2"],
    p_frailty_ABS_INJ_3 = df_frailty["pe", "ABS_INJ_3"],
    p_frailty_REL_INJ_1 = 1,
    p_frailty_REL_INJ_2 = df_frailty["pe", "REL_INJ_2"],
    p_frailty_REL_INJ_3 = df_frailty["pe", "REL_INJ_3"],
    p_frailty_OD_INJ_1  = 1,
    p_frailty_OD_INJ_2  = df_frailty["pe", "OD_INJ_2"],
    p_frailty_OD_INJ_3  = df_frailty["pe", "OD_INJ_3"],

    # Load weibull parameters
    # Weibull scale
    #p_weibull_scale_BUP_NI = df_weibull_scale["pe", "BUP_NI"],
    p_weibull_scale_MET_NI = df_weibull_scale["pe", "MET_NI"],
    p_weibull_scale_REL_NI = df_weibull_scale["pe", "REL_NI"],
    p_weibull_scale_OD_NI  = df_weibull_scale["pe", "OD_NI"],
    p_weibull_scale_ABS_NI = df_weibull_scale["pe", "ABS_NI"],
    #p_weibull_scale_BUP_INJ = df_weibull_scale["pe", "BUP_INJ"],
    p_weibull_scale_MET_INJ = df_weibull_scale["pe", "MET_INJ"],
    p_weibull_scale_REL_INJ = df_weibull_scale["pe", "REL_INJ"],
    p_weibull_scale_OD_INJ  = df_weibull_scale["pe", "OD_INJ"],
    p_weibull_scale_ABS_INJ = df_weibull_scale["pe", "ABS_INJ"],

    # Weibull shape
    #p_weibull_shape_BUP_NI = df_weibull_shape["pe", "BUP_NI"],
    p_weibull_shape_MET_NI = df_weibull_shape["pe", "MET_NI"],
    p_weibull_shape_REL_NI = df_weibull_shape["pe", "REL_NI"],
    p_weibull_shape_OD_NI  = df_weibull_shape["pe", "OD_NI"],
    p_weibull_shape_ABS_NI = df_weibull_shape["pe", "ABS_NI"],
    #p_weibull_shape_BUP_INJ = df_weibull_shape["pe", "BUP_INJ"],
    p_weibull_shape_MET_INJ = df_weibull_shape["pe", "MET_INJ"],
    p_weibull_shape_REL_INJ = df_weibull_shape["pe", "REL_INJ"],
    p_weibull_shape_OD_INJ  = df_weibull_shape["pe", "OD_INJ"],
    p_weibull_shape_ABS_INJ = df_weibull_shape["pe", "ABS_INJ"],

    # Unconditional transition probabilities
    # Non-Injection
    # From MET1 & MET
    #p_MET1_BUP1_NI = df_UP["MET_NI", "BUP1_NI"],
    #p_MET_BUP1_NI  = df_UP["MET_NI", "BUP1_NI"],
    #p_MET1_ABS_NI  = df_UP["MET_NI", "ABS_NI"],
    p_MET_ABS_NI   = df_UP["MET_NI", "ABS_NI"],
    #p_MET1_REL1_NI = df_UP["MET_NI", "REL1_NI"],
    p_MET_REL1_NI  = df_UP["MET_NI", "REL1_NI"],
    #p_MET1_OD_NI   = df_UP["MET_NI", "OD_NI"],
    p_MET_OD_NI    = df_UP["MET_NI", "OD_NI"],
    # From ABS
    p_ABS_REL1_NI = df_UP["ABS_NI", "REL1_NI"],
    p_ABS_OD_NI   = df_UP["ABS_NI", "OD_NI"],
    # From REL1 & REL
    p_REL1_REL_NI = df_UP["REL1_NI", "REL_NI"],
    p_REL1_MET_NI = df_UP["REL1_NI", "MET_NI"],
    p_REL_MET_NI  = df_UP["REL_NI", "MET_NI"],
    #p_REL1_BUP1_NI = df_UP["REL_NI", "BUP1_NI"],
    #p_REL_BUP1_NI  = df_UP["REL_NI", "BUP1_NI"],
    p_REL1_ABS_NI  = df_UP["REL1_NI", "ABS_NI"],
    p_REL_ABS_NI   = df_UP["REL_NI", "ABS_NI"],
    p_REL1_OD_NI   = df_UP["REL1_NI", "OD_NI"],
    p_REL_OD_NI    = df_UP["REL_NI", "OD_NI"],
    # From OD
    p_OD_MET_NI  = df_UP["OD_NI", "MET_NI"],
    #p_OD_BUP1_NI  = df_UP["OD_NI", "BUP1_NI"],
    p_OD_ABS_NI   = df_UP["OD_NI", "ABS_NI"],
    p_OD_REL1_NI  = df_UP["OD_NI", "REL1_NI"],

    # Inj
    # From BUP1 & BUP
    #p_BUP1_MET1_INJ = df_UP["BUP_INJ", "MET1_INJ"],
    #p_BUP_MET1_INJ  = df_UP["BUP_INJ", "MET1_INJ"],
    #p_BUP1_ABS_INJ  = df_UP["BUP_INJ", "ABS_INJ"],
    #p_BUP_ABS_INJ   = df_UP["BUP_INJ", "ABS_INJ"],
    #p_BUP1_REL1_INJ = df_UP["BUP_INJ", "REL1_INJ"],
    #p_BUP_REL1_INJ  = df_UP["BUP_INJ", "REL1_INJ"],
    #p_BUP1_OD_INJ   = df_UP["BUP_INJ", "OD_INJ"],
    #p_BUP_OD_INJ    = df_UP["BUP_INJ", "OD_INJ"],
    # From MET1 & MET
    #p_MET1_BUP1_INJ = df_UP["MET_INJ", "BUP1_INJ"],
    #p_MET_BUP1_INJ  = df_UP["MET_INJ", "BUP1_INJ"],
    #p_MET1_ABS_INJ  = df_UP["MET_INJ", "ABS_INJ"],
    p_MET_ABS_INJ   = df_UP["MET_INJ", "ABS_INJ"],
    #p_MET1_REL1_INJ = df_UP["MET_INJ", "REL1_INJ"],
    p_MET_REL1_INJ  = df_UP["MET_INJ", "REL1_INJ"],
    #p_MET1_OD_INJ   = df_UP["MET_INJ", "OD_INJ"],
    p_MET_OD_INJ    = df_UP["MET_INJ", "OD_INJ"],
    # From ABS
    p_ABS_REL1_INJ = df_UP["ABS_INJ", "REL1_INJ"],
    p_ABS_OD_INJ   = df_UP["ABS_INJ", "OD_INJ"],
    # From REL1
    p_REL1_REL_INJ = df_UP["REL1_INJ", "REL_INJ"],
    p_REL1_MET_INJ = df_UP["REL1_INJ", "MET_INJ"],
    p_REL_MET_INJ  = df_UP["REL_INJ", "MET_INJ"],
    #p_REL1_BUP1_INJ = df_UP["REL_INJ", "BUP1_INJ"],
    #p_REL_BUP1_INJ  = df_UP["REL_INJ", "BUP1_INJ"],
    p_REL1_ABS_INJ  = df_UP["REL1_INJ", "ABS_INJ"],
    p_REL_ABS_INJ   = df_UP["REL_INJ", "ABS_INJ"],
    p_REL1_OD_INJ   = df_UP["REL1_INJ", "OD_INJ"],
    p_REL_OD_INJ    = df_UP["REL_INJ", "OD_INJ"],
    # From OD
    p_OD_MET_INJ  = df_UP["OD_INJ", "MET_INJ"],
    #p_OD_BUP1_INJ  = df_UP["OD_INJ", "BUP1_INJ"],
    p_OD_ABS_INJ   = df_UP["OD_INJ", "ABS_INJ"],
    p_OD_REL1_INJ  = df_UP["OD_INJ", "REL1_INJ"],

    # HIV Seroconversion
    # Non-injection
    #p_sero_BUP1_NI = df_sero["pe", "BUP_NI"],
    #p_sero_BUP_NI  = df_sero["pe", "BUP_NI"],
    #p_sero_MET1_NI = df_sero["pe", "MET_NI"],
    p_sero_MET_NI  = df_sero["pe", "MET_NI"],
    p_sero_REL1_NI = df_sero["pe", "REL1_NI"],
    p_sero_REL_NI  = df_sero["pe", "REL_NI"],
    p_sero_OD_NI   = df_sero["pe", "OD_NI"],
    p_sero_ABS_NI  = df_sero["pe", "ABS_NI"],
    
    # Injection
    #p_sero_BUP1_INJ = df_sero["pe", "BUP_INJ"],
    #p_sero_BUP_INJ  = df_sero["pe", "BUP_INJ"],
    p_sero_MET1_INJ = df_sero["pe", "MET1_INJ"],
    p_sero_MET_INJ  = df_sero["pe", "MET_INJ"],
    p_sero_REL1_INJ = df_sero["pe", "REL1_INJ"],
    p_sero_REL_INJ  = df_sero["pe", "REL_INJ"],
    p_sero_OD_INJ   = df_sero["pe", "OD_INJ"],
    p_sero_ABS_INJ  = df_sero["pe", "ABS_INJ"],

    # Costs
    # Treatment Costs
    #c_BUP_TX  = df_costs["pe", "BUP_TX"],
    c_MET_TX  = df_costs["pe", "MET_TX"],
    c_OD_TX  = df_costs["pe", "OD_TX"],
    
    # HRU Costs
    # Modify if age-specific
    #c_BUP_NI_HRU = df_costs["pe", "BUP_NI_HRU"],
    c_MET_NI_HRU = df_costs["pe", "MET_NI_HRU"],
    c_REL_NI_HRU = df_costs["pe", "REL_NI_HRU"],
    c_OD_NI_HRU = df_costs["pe", "OD_NI_HRU"],
    c_ABS_NI_HRU = df_costs["pe", "ABS_NI_HRU"],
    #c_BUP_INJ_HRU = df_costs["pe", "BUP_INJ_HRU"], 
    c_MET_INJ_HRU = df_costs["pe", "MET_INJ_HRU"],
    c_REL_INJ_HRU = df_costs["pe", "REL_INJ_HRU"],
    c_OD_INJ_HRU = df_costs["pe", "OD_INJ_HRU"],
    c_ABS_INJ_HRU = df_costs["pe", "ABS_INJ_HRU"], 

    # HIV Costs
    c_HIV_HRU = df_costs["pe", "HIV_HRU"],
    c_HIV_ART = df_costs["pe", "HIV_ART"],
    
    # Crime Costs
    #df_crime_costs = subset(df_crime_costs, type >= "pe_25" & type <= "pe_90"),
    df_crime_costs = subset(df_crime_costs, type=="pe"),
    # Age-specific
    #v_c_BUP_NI_crime = df_crime_costs %>% select(BUP_NI) %>% as.matrix(),
    v_c_MET_NI_crime = df_crime_costs %>% select(MET_NI) %>% as.matrix(),
    v_c_REL_NI_crime = df_crime_costs %>% select(REL_NI) %>% as.matrix(),
    v_c_OD_NI_crime = df_crime_costs %>% select(OD_NI) %>% as.matrix(),
    v_c_ABS_NI_crime = df_crime_costs %>% select(ABS_NI) %>% as.matrix(),
    
    #v_c_BUP_INJ_crime = df_crime_costs %>% select(BUP_INJ) %>% as.matrix(),
    v_c_MET_INJ_crime = df_crime_costs %>% select(MET_INJ) %>% as.matrix(),
    v_c_REL_INJ_crime = df_crime_costs %>% select(REL_INJ) %>% as.matrix(),
    v_c_OD_INJ_crime = df_crime_costs %>% select(OD_INJ) %>% as.matrix(),
    v_c_ABS_INJ_crime = df_crime_costs %>% select(ABS_INJ) %>% as.matrix(), 

    # QALYs
    #u_BUP_NI_NEG = df_qalys["pe", "BUP_NI_NEG"],
    u_MET_NI_NEG = df_qalys["pe", "MET_NI_NEG"],
    u_REL_NI_NEG = df_qalys["pe", "REL_NI_NEG"],
    u_OD_NI_NEG = df_qalys["pe", "OD_NI_NEG"],
    u_ABS_NI_NEG = df_qalys["pe", "ABS_NI_NEG"],
    
    #u_BUP_INJ_NEG = df_qalys["pe", "BUP_INJ_NEG"],
    u_MET_INJ_NEG = df_qalys["pe", "MET_INJ_NEG"],
    u_REL_INJ_NEG = df_qalys["pe", "REL_INJ_NEG"],
    u_OD_INJ_NEG = df_qalys["pe", "OD_INJ_NEG"],
    u_ABS_INJ_NEG = df_qalys["pe", "ABS_INJ_NEG"],
    
    #u_BUP_NI_POS = df_qalys["pe", "BUP_NI_POS"],
    u_MET_NI_POS = df_qalys["pe", "MET_NI_POS"],
    u_REL_NI_POS = df_qalys["pe", "REL_NI_POS"],
    u_OD_NI_POS = df_qalys["pe", "OD_NI_POS"],
    u_ABS_NI_POS = df_qalys["pe", "ABS_NI_POS"],
    
    #u_BUP_INJ_POS = df_qalys["pe", "BUP_INJ_POS"],
    u_MET_INJ_POS = df_qalys["pe", "MET_INJ_POS"],
    u_REL_INJ_POS = df_qalys["pe", "REL_INJ_POS"],
    u_OD_INJ_POS = df_qalys["pe", "OD_INJ_POS"],
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
update_param_list <- function(l_params_all, params_updated){
  l_params_all <- modifyList(l_params_all, params_updated) #update values
  return(l_params_all)
}
#####################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Khan M, Daly-Grafstein D, Erdelyi S, Staples JA. AMA-OD_cco create table of cohort characteristics. 2023 Oct 3. Retrieved from: ***LINK*** 
# cco Table 1 creation
# Author: Nicole Hu 
# Date: 2023-09-13
# Updated: 2023-10-03
#####################################################
# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")

# set working dir
setwd("R:/working/AMA-OD_cco/NH/results/Table A")

################## Combine all variables ############
cohort <- readRDS("cco_cohort.rds")
cohort <- cohort %>%
  left_join(readRDS("cco_demo.rds")) %>% 
  left_join(readRDS("cco_medical_hist.rds")) %>% 
  left_join(readRDS("cco_medication.rds")) %>% 
  mutate(male_sex = ifelse(gender == "M",1,0)) %>% 
  mutate(age_group = case_when(index_age < 30 ~ '18-29',
                               index_age < 50 ~ '30-49',
                               TRUE ~ '50+')) %>% 
  mutate(geq1_num_psyc_hosp_prev_6m = as.numeric(num_psyc_hosp_prev_6m > 0)) %>%
  mutate(geq1_oat = as.numeric(oat > 0)) %>% 
  mutate(geq1_opioid = as.numeric(opioid > 0)) %>% 
  mutate(geq1_benzo = as.numeric(benzo > 0)) %>% 
  mutate(geq1_antipsyc = as.numeric(antipsyc > 0)) %>% 
  mutate(other_medication = ifelse(other < 20, other, "20+")) %>% 
  mutate(num_active_cat = case_when(num_active == 0 ~ "0 or 1",
                                    num_active == 1 ~ "0 or 1",
                                    TRUE ~ "2+")) %>% 
  mutate(sdpr_pay_cat = as.numeric(sdpr_pay > 0)) %>% 
  mutate(urban = case_when(substring(od_event_fsa_16,2,2) == "0" ~ "rural",
                           substring(od_event_fsa_16,2,2) != "0" ~ "urban",
                           TRUE ~ "missing FSA"))

#categorical variables
cohort <- cohort %>% 
  mutate(
    dc_bma = factor(as.numeric(discharge == "bma")),
    dc_wa = factor(as.numeric(discharge == "wa")),
    homeless_hist = factor(homeless_hist),
    alcohol = factor(alcohol),
    drugs_nonopioid = factor(drugs_nonopioid),
    psych = factor(psych),
    cci_cat = factor(cci_cat),
    oat_active_ad = factor(oat_active_ad),
    geq1_benzo = factor(geq1_benzo),
    geq1_opioid = factor(geq1_opioid),
    geq1_antipsyc = factor(geq1_antipsyc),
    case = as.numeric(interval == "od")) 

#reset reference level
cohort$discharge <- relevel(as.factor(cohort$discharge), ref = "no discharge")

#save cohort with +- 6 month and +12 month controls for sensitivity analysis
cohort_lbcontrol <- cohort %>% filter(interval %in% c("od","control_1","control_lb6m","control_lb12m") )
saveRDS(cohort_lbcontrol,"cco_coh_lbcontrol.rds")

#save cohort with 3 controls for sensitivity analysis
cohort_3control <- cohort %>% filter(interval %in% c("od","control_1","control_2","control_3") )
saveRDS(cohort_3control,"cco_coh_3control.rds")

#save cohort with 2 controls for primary analysis
cohort_2control <- cohort %>% filter(interval %in% c("od","control_1","control_2") )
saveRDS(cohort_2control, "cco_coh_all.rds") 

#################### Table 1 #######################
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")
# Vector of variables to summarize
table1_vars <- c(
  # demographic
  "index_age", "age_group","gender","male_sex", "pop_density_class", "urban","sdpr_pay_cat",
  
  # medical hist
  "num_hosp_prev_6m","geq1_hosp_prev_6m","tdays_prev_6m",
  "geq1_num_psyc_hosp_prev_6m", "num_clinic_prev_6m", "geq7_clinic_prev_6m", 
  
  
  # comorbidities
  "homeless_hist",
  "alcohol", "drugs_all", "drugs_opioid","drugs_nonopioid","any_sbstcs",
  "mood","anxiety","dm_nc_or_compl","mi","chf","cvd","copd","dem","renal",
  "cancer","afib","osa", "psych","cihd", "htn","endocarditis","om","hiv",  
  "cci", "cci_cat","ivdu_hist",
  
  # medication
  "num_active","num_active_cat", "oat_active_ad" ,
  "geq1_oat", "geq1_benzo", "geq1_opioid", "geq1_antipsyc", "other_medication",
  
  # other
  "study_year","health_region",
  
  # exposure
  "discharge"
)


## Vector of non-normal variables 
non_norm_vars <- c("index_age", "num_clinic_prev_6m", "num_hosp_prev_6m", "tdays_prev_6m",
                   "cci", "num_active")


## Vector of categorical variables
cat_vars <- c("age_group", "male_sex","gender", "pop_density_class","sdpr_pay_cat",
              "geq1_hosp_prev_6m","geq7_clinic_prev_6m",
              "geq1_num_psyc_hosp_prev_6m",
              
              "homeless_hist",
              "alcohol", "drugs_all", "drugs_opioid","drugs_nonopioid","any_sbstcs",
              "mood","anxiety","dm_nc_or_compl","mi","chf","cvd","copd","dem","renal","hiv", 
              "cancer","afib","osa", "psych","cihd", "htn","endocarditis","om", 
              "cci", "cci_cat","ivdu_hist",
              
              "oat_active_ad","num_active_cat",
              "geq1_oat", "geq1_benzo", "geq1_opioid", "geq1_antipsyc", "other_medication",
              
              "study_year","health_region",
              
              "discharge"
)


# save
tab1a <- CreateTableOne(data = cohort, 
                        vars = table1_vars, factorVars = cat_vars)
tab1a_csv <- print(tab1a, nonnormal = non_norm_vars, showAllLevels = T, formatOptions = list(big.mark=','), printToggle = F)
write.csv(tab1a_csv, file = "AMA-OD_cco - Table 1a v6.csv")


# Group by unpooled pre-exposure 
cohort$interval <- factor(cohort$interval, levels = c("od","control_1",'control_2'))
tab1b <- CreateTableOne(data = cohort, 
                        strata = 'interval', 
                        vars = table1_vars, factorVars = cat_vars)
tab1b_csv <- print(tab1b, nonnormal = non_norm_vars, showAllLevels = T, formatOptions = list(big.mark=','), printToggle = F,test = F)

# Group by pooled pre-exposure
tab1c <- CreateTableOne(data = cohort, 
                        strata = 'interval2', 
                        vars = table1_vars, factorVars = cat_vars)
tab1c_csv <- print(tab1c, nonnormal = non_norm_vars, showAllLevels = T, formatOptions = list(big.mark=','),smd = T, printToggle = F)
tab1c_csv <- tab1c_csv[,c(-1,-3)]

#combine pooled and unpooled together
tab1b_csv <- cbind(tab1b_csv,tab1c_csv)
write.csv(tab1b_csv, file = "AMA-OD_cco - Table 1b v6.csv")







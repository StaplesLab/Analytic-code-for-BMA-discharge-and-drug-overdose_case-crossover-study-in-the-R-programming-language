###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Khan M, Daly-Grafstein D, Erdelyi S, Staples JA. AMA-OD_cco create demographic variables. 2023 Sep 20. Retrieved from: https://github.com/StaplesLab/Analytic-code-for-BMA-discharge-and-drug-overdose_case-crossover-study-in-the-R-programming-language/edit/main/AMA-OD_cco%202a%20Create%20varaible%20-%20demographic%20v1.R
# Create demographic variables
# Author: Xiao (Nicole) Hu 
# Date: 2023-09-08
# Updated: 2023-09-20
##########################################################
# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")

#read in cohort
setwd("R:/working/AMA-OD_cco/NH/results/Table A")
cohort <- readRDS("cco_cohort.rds")

# cohort for sensitivity analyses 
#setwd("R:/working/AMA-OD_cco/NH/results/Table D")
#cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table D/allow prior od/cco_cohort_pod.rds") #cohort including od with prior od
#cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table D/nod_coh.rds") #non overdose cohort
#cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table D/od_date/cco_cohort_date.rds") #cohort with alternative overdose date

#population density classification
fsa_xw <- read_delim("R:/working/XW files/fsa_2016_urbanRural.csv")
pop_density <- cohort %>% 
  left_join(fsa_xw, by = c("od_event_fsa_16" = "fsa")) %>% 
  select(episode_id, date, fsaType) %>% 
  mutate(pop_density_class = case_when(
    fsaType %in% c("popCentreMedium", "popCentreLarge") ~ "Urban",
    fsaType %in% c("ruralArea","popCentreSmall") ~ "Rural",
    TRUE ~ "Missing"
  ))

#health authority
lu_pnet_ha_hsda_lha <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_lu_pnet_ha_hsda_lha.csv.gz")
health_region <- cohort %>% 
  left_join(lu_pnet_ha_hsda_lha, by = c("od_event_lha_18" = "lha_id")) %>% 
  select(episode_id, date,ha_area) %>% 
  mutate(health_region = ifelse(ha_area == "09 Unknown HA" | is.na(ha_area), "Unknown/Out of Province HA",ha_area))

#homelessness (blocker: use ICD10, only available for hosp)


#social assistance payments in the past 6 month
sdpr <- read_csv("R:/DATA/2023-08-15_230952/vw_sdpr.csv")
sdpr_pay <- cohort %>% 
  left_join(sdpr %>% select(moh_study_id, sdpr_date, pay), relationship = "many-to-many") %>% 
  filter(sdpr_date <= date %m-% days(1) &
           sdpr_date >= date %m-% days(182)) %>% 
  group_by(episode_id,date) %>% 
  summarise(sdpr_pay = sum(pay)) %>% 
  ungroup() %>% 
  right_join(cohort) %>% 
  mutate(sdpr_pay = replace_na(sdpr_pay,0)) %>% 
  select(episode_id, date,sdpr_pay)

#drug toxicity
drug_toxicity <- read.csv("R:/working/AMA-OD_coh/NH/results/Table A/drug_toxicity.csv")
drug_toxicity <- cohort %>% 
  mutate(ym = substr(date,1,7)) %>% 
  left_join(drug_toxicity, by = "ym") %>% 
  select(episode_id,date, drug_toxicity = n)

demo <- pop_density %>% 
  full_join(health_region) %>% 
  full_join(sdpr_pay) %>% 
  full_join(drug_toxicity)

saveRDS(demo,"cco_demo.rds")

#saveRDS(demo,"allow prior od/cco_demo_pod.rds")
#sdpr_pay %>% full_join(drug_toxicity) %>% saveRDS("R:/working/AMA-OD_cco/NH/results/Table D/nod_cco_demo.rds")
#sdpr_pay %>% full_join(drug_toxicity) %>% saveRDS("R:/working/AMA-OD_cco/NH/results/Table D/od_date/date_cco_demo.rds")

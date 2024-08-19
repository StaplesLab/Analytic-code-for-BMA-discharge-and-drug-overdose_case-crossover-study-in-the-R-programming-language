################################################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Daly-Grafstein D, Khan M, Erdelyi S, Staples JA. AMA-OD_cco sensitivity analyses case-time control analysis. 2023 Dec 22. Retrieved from: ***LINK*** 
# Case-time-control analysis
# Author: Nicole Hu
# Updated: Dec 22, 2023
################################################################################

source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")
setwd("R:/working/AMA-OD_cco/NH/results/Table D")
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")
cco_adj <- readRDS("R:/working/AMA-OD_cco/NH/results/Table C/cco_adj.rds")
person_xwalk <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_person_xwalk.csv.gz")
cohort <- cohort %>% left_join(person_xwalk %>% select(moh_study_id,matching_id))

get_coef <- function(m,se = F){
  summary <- summary(m) 
  model_summary <- 
    cbind(summary$coefficients[,-1] ,summary$conf.int[,c(3,4)])
  model_summary <- as.data.frame(round(model_summary,6))
  model_summary$coef <- rownames(model_summary)
  var_xw <- readxl::read_xlsx("R:/working/XW files/AMA-OD_cco coef_to_variable.xlsx")
  model_summary <- model_summary %>% left_join(var_xw, by = "coef") %>% 
    select(-type) %>% 
    select(variable = coef, `short description` = variable,`exp(coef)`, `se(coef)`, `robust se`, `Pr(>|z|)`, `lower .95`, `upper .95`) %>% 
    mutate(`exp, CI, p` = paste0(sprintf("%.2f",`exp(coef)`),", ",sprintf("%.2f",`lower .95`), "-",sprintf("%.2f",`upper .95`), ", ", "p ", ifelse(`Pr(>|z|)` < 0.001, "<0.001", paste("=", sprintf("%.3f",`Pr(>|z|)`))))) 
  
  if(se == T){
    model_summary <- model_summary %>% select(variable,`short description`,`exp, CI, p`, `robust se`, se = `se(coef)`)
  }else{
    model_summary <- model_summary %>% select(variable,`short description`,`exp, CI, p`, `robust se`)
  }
  return(model_summary)
}

######################## non_od_cohort #########################################
# get matchings for current cco cohort from person xwalk [each case has four controls]
person_coh <- person_xwalk %>% 
  filter(moh_study_id %in% cohort$moh_study_id) %>% 
  select(matching_id) %>% 
  left_join(person_xwalk) %>% 
  group_by(matching_id) %>% 
  mutate(last_od_date = max(last_od_date,na.rm = T)) %>% 
  ungroup() %>% 
  #make sure controls are alive at overdose date
  filter(death_date >= last_od_date | is.na(death_date)) %>% 
  arrange(matching_id)

# randomly select one control for each case
set.seed(3456)
person_coh_control <- person_coh %>% 
  filter(case_control_flag == 0) %>% 
  left_join(cohort %>% select(matching_id,date)) %>% 
  mutate(index_age = as.period(interval(birth_date, date))@year) %>%
  group_by(moh_study_id) %>% 
  #make sure control is at least 18 years old
  mutate(min_index_age = min(index_age)) %>% 
  filter(min_index_age >= 18) %>%
  group_by(matching_id) %>% 
  sample_n(1) %>% 
  ungroup() 

# check cohort size
person_coh %>% filter(case_control_flag == 1) %>% select(moh_study_id) %>% n_distinct()  
length(unique(cohort$moh_study_id))
# 125 individuals for 130 od cases
no_match <- person_xwalk %>% filter(matching_id %in% (person_coh %>% count(matching_id) %>% filter(n < 2) %>% pull(matching_id))) 


# create non_od_cohort
nod_coh <- person_coh_control %>% 
  select(moh_study_id,gender,birth_date,matching_id) %>% 
  left_join(cohort %>% select(episode_id,date,interval,matching_id),by = "matching_id") %>% 
  mutate(index_age = as.period(interval(birth_date, date))@year) %>% 
  arrange(episode_id)




######################## outcome and exposure ##################################
dad <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_dad.csv.gz")

#exposure
nod_coh <- nod_coh %>% 
  left_join(dad,
                     by = "moh_study_id",
                     relationship = "many-to-many") %>% 
  filter(sep_date >= (date-days(28)) & sep_date < date) %>% 
  group_by(episode_id,date) %>% 
  arrange(desc(sep_date),desc(ad_date)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(episode_id,date,sep_disp) %>%
  distinct() %>% 
  mutate(discharge = case_when(sep_disp %in%  c(6,12,61,62,63,64,65) ~ "bma",
                               TRUE ~ "wa")) %>% 
  right_join(nod_coh) %>% 
  mutate(discharge = replace_na(discharge,"no discharge")) %>% 
  arrange(episode_id, date) %>% 
  mutate(
  dc_bma = factor(as.numeric(discharge == "bma")),
  dc_wa = factor(as.numeric(discharge == "wa")))

#case
nod_coh <- nod_coh %>% 
  mutate(case = as.numeric(interval == "od")) %>% 
  mutate(episode_id = episode_id + 50000000)

saveRDS(nod_coh,"R:/working/AMA-OD_cco/NH/results/Table D/nod_coh.rds")

################################## covariates ##################################
nod_coh <- readRDS("R:/working/AMA-OD_cco/NH/results/Table D/nod_coh.rds")

nod_coh <- nod_coh %>%
  left_join(readRDS("nod_cco_demo.rds")) %>% 
  left_join(readRDS("nod_cco_medical_hist.rds")) %>% 
  left_join(readRDS("nod_cco_medication.rds")) %>% 
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
                                    TRUE ~ "2+"))

nod_coh <- nod_coh %>% 
  mutate(
    homeless_hist = factor(homeless_hist),
    alcohol = factor(alcohol),
    drugs_nonopioid = factor(drugs_nonopioid),
    psych = factor(psych),
    cci_cat = factor(cci_cat),
    oat_active_ad = factor(oat_active_ad),
    geq1_benzo = factor(geq1_benzo),
    geq1_opioid = factor(geq1_opioid),
    geq1_antipsyc = factor(geq1_antipsyc)) 

nod_coh$discharge <- relevel(as.factor(nod_coh$discharge), ref = "no discharge")

####################### model with non-od cohort ###############################
cco_adj_nod <- clogit(case ~ discharge + homeless_hist + 
                         sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                         drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                         geq1_benzo + geq1_antipsyc + drug_toxicity + strata(episode_id), 
                       data = nod_coh,
                       cluster = moh_study_id,
                       method  = "efron",
                       robust = T)

####################### model with interaction term ############################
od_cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")
person_xwalk <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_person_xwalk.csv.gz")
od_cohort <- od_cohort %>% left_join(person_xwalk %>% select(moh_study_id,matching_id))
combined_cohort <- od_cohort %>% 
  select(colnames(nod_coh)) %>% 
  filter(matching_id %in% nod_coh$matching_id) %>% 
  mutate(is_od = 1) %>% 
  rbind(nod_coh %>% mutate(is_od = 0)) %>% 
  mutate(is_od = factor(is_od))

cco_adj_int <- clogit(case ~ discharge + homeless_hist + 
                         sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                         drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                         geq1_benzo + geq1_antipsyc + drug_toxicity + is_od*discharge+
                         strata(episode_id), 
                       data = combined_cohort,
                       cluster = moh_study_id,
                       method  = "efron",
                       robust = T)

#check
exp(log(0.34) + log(6.21)) #primary: 2.08
exp(log(1.02) + log(1.38)) #primary: 1.39

############################ Characteristics ###################################
# library(gtsummary)
# table <- tbl_summary(data = combined_cohort, by = is_od)
# print(table)


# Vector of variables to summarize
table1_vars <- c(
  # demographic
  "index_age", "age_group","gender","male_sex", "sdpr_pay",
  
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
  
  
  # exposure
  "discharge"
)


## Vector of non-normal variables 
non_norm_vars <- c("index_age","sdpr_pay", "num_clinic_prev_6m", "num_hosp_prev_6m", "tdays_prev_6m",
                   "cci", "num_active")


## Vector of categorical variables
cat_vars <- c("age_group", "male_sex","gender",
              "geq1_hosp_prev_6m","geq7_clinic_prev_6m",
              "geq1_num_psyc_hosp_prev_6m",
              
              "homeless_hist",
              "alcohol", "drugs_all", "drugs_opioid","drugs_nonopioid","any_sbstcs",
              "mood","anxiety","dm_nc_or_compl","mi","chf","cvd","copd","dem","renal","hiv", 
              "cancer","afib","osa", "psych","cihd", "htn","endocarditis","om", 
              "cci", "cci_cat","ivdu_hist",
              
              "oat_active_ad","num_active_cat",
              "geq1_oat", "geq1_benzo", "geq1_opioid", "geq1_antipsyc", "other_medication",
              
              
              "discharge"
)


# save
tab1a <- CreateTableOne(data = combined_cohort,
                        strata = "is_od",
                        vars = table1_vars, 
                        factorVars = cat_vars)
tab1a_csv <- print(tab1a, nonnormal = non_norm_vars, showAllLevels = F, formatOptions = list(big.mark=','), printToggle = F,smd = T)
write.csv(tab1a_csv, file = "R:/working/AMA-OD_cco/NH/results/Table D/AMA-OD_cco - Case-time-control Table 1a v2.csv")

tab1a_csv <- read.csv("R:/working/AMA-OD_cco/NH/results/Table D/AMA-OD_cco - Case-time-control Table 1a v2.csv")
write.xlsx(list("Characteristics" = tab1a_csv, 
                "non-od group" = get_coef(cco_adj_nod),
                "combined group" = get_coef(cco_adj_int)),
           "R:/working/AMA-OD_cco/NH/results/Table D/AMA-OD_cco - Case-time-control v2.xlsx")


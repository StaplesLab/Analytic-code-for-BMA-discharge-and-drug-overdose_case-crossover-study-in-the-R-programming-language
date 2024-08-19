###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Daly-Grafstein D, Khan M, Erdelyi S, Staples JA. AMA-OD_cco sensitivity analyses using earliest overdose encounter date as overdose date. 2023 Dec 19. Retrieved from: ***LINK*** 
# Use earliest overdose encounter date as overdose date
# Author: Xiao (Nicole) Hu 
# Date: 2023-11-20
# Updated: 2023-12-19
##########################################################

source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")
setwd("R:/working/AMA-OD_cco/NH/results/Table D/od_date")

# read-in episode and encounter data
ee_xwalk <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_episode_encounter_xwalk.csv.gz")
encounter_master <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_encounter_master.csv.gz")
od_cohort <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_episode_master.csv.gz") %>% 
  filter(od_flag == 1)

#link episode with encounter
od_cohort_encounter <- od_cohort %>% 
  select(episode_id,moh_study_id) %>% 
  left_join(ee_xwalk %>% select(-to_keep)) %>% 
  left_join(encounter_master %>% select(encounter_id,ec_start = start_dt_tm, ec_end = stop_dt_tm, ec_fatal = fatal)) %>% #ec_fatal: if the encounter is fatal
  group_by(episode_id) %>% 
  arrange(episode_id, ec_start, .by_group=TRUE) %>% 
  ungroup()

#use earliest od encounter date as the new overdose date
new_date <- od_cohort_encounter %>% filter(od_flag == 1) %>% group_by(episode_id) %>% summarise(od_date = as.Date(min(ec_start)),od_date_tm = min(ec_start))
od_cohort_encounter <- od_cohort_encounter %>% left_join(new_date)
od_cohort <- od_cohort_encounter %>% select(episode_id,moh_study_id,od_date,od_date_tm) %>% distinct()

#create cohort
dad <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_dad.csv.gz")
person_xwalk <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_person_xwalk.csv.gz")
od_cohort <- od_cohort %>% mutate(control_date_1 = as.Date(od_date) %m-% days(182), #6 month
                                  control_date_2 = as.Date(od_date) %m-% days(364), #12 month
                                  control_date_3 = as.Date(od_date) %m-% days(546)) #18 month

# exclude those <18 years old at any pre-control interval
od_filter <- od_cohort %>% 
  left_join(person_xwalk) %>% 
  mutate(control_age = (as.period(interval(birth_date, control_date_3 - days(28))))@year) %>% 
  mutate(younger_18 = as.numeric(control_age < 18 | control_age > 110))

# exclude those with previous od in 18 months + 30 days prior to index od
prior_od <- od_cohort %>% 
  left_join(od_cohort %>% select(moh_study_id,pre_od_date = od_date),
            by = "moh_study_id", relationship = "many-to-many") %>% 
  # filter if od within 18 months + 30 days before index od
  filter(pre_od_date < od_date & 
           pre_od_date >= (od_date %m-% (days(546)+days(28))))

od_filter <- od_filter %>% 
  mutate(prior_od = as.numeric(episode_id %in% prior_od$episode_id)) 

#exclude those earlier than 2015-01-01 + 18 month + 28 days (2016-07-28)
od_filter <- od_filter %>% 
  mutate(earlier_16 = as.numeric(od_date %m-% (days(546) + days(28))) < as.Date("2015-01-01"))


# analytic data set for cco
od_cco <- od_filter %>% 
  #apply exclusion criteria
  filter(prior_od == 0, younger_18 == 0, earlier_16 == 0)

#total cohort
nrow(od_cco) #N = 27608
od_cco %>% select(moh_study_id) %>% n_distinct() #n = 26614


# make control as a row itself 
od_cco <- od_cco %>% 
  mutate(expose_od = NA, expose_control_1 = NA, expose_control_2 = NA,expose_control_3 = NA) %>% 
  select(moh_study_id, od_date,episode_id,expose_od,expose_control_1,expose_control_2,expose_control_3) %>% 
  mutate(id = seq(1:nrow(od_cco))) %>% 
  pivot_longer(c(expose_od, expose_control_1, expose_control_2,expose_control_3)) %>% 
  mutate(date = case_when(
    name == "expose_od" ~ od_date, 
    name == "expose_control_1" ~ as.Date(od_date) %m-% days(182),
    name == "expose_control_2" ~ as.Date(od_date) %m-% days(364),
    name == "expose_control_3" ~ as.Date(od_date) %m-% days(546))) %>% 
  mutate(name = case_when(
    name == "expose_od" ~ "od", 
    name == "expose_control_1" ~ "control_1",
    name == "expose_control_2" ~ "control_2",
    name == "expose_control_3" ~ "control_3")) %>% 
  rename(interval = name) %>%
  select(-value) %>% 
  mutate(interval2 = ifelse(interval == "od","od","control"))

# add exposure variable: based on most recent discharge within 30 days before od/control dates
od_cco <- od_cco %>% 
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
  right_join(od_cco) %>% 
  mutate(discharge = replace_na(discharge,"no discharge")) %>% 
  arrange(episode_id, date)


# index age
cohort <- od_cco %>% 
  left_join(person_xwalk %>% select(-reference,-matching_id,-case_control_flag,-to_keep,od_death = fatal_od_case)) %>% 
  mutate(index_age = as.period(interval(birth_date, date))@year)

cohort <- cohort %>% filter(interval != "control_3") %>% arrange(episode_id)

saveRDS(cohort, "cco_cohort_date.rds")

######################## covariates ###################################
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table D/od_date/cco_cohort_date.rds")
cohort <- cohort %>%
  left_join(readRDS("date_cco_demo.rds")) %>% 
  left_join(readRDS("date_cco_medical_hist.rds")) %>% 
  left_join(readRDS("date_cco_medication.rds")) %>% 
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
    geq1_antipsyc = factor(geq1_antipsyc)) 
cohort$discharge <- relevel(as.factor(cohort$discharge), ref = "no discharge")

########################## model ###################################
cohort <- cohort %>% mutate(case = as.numeric(interval == "od"))

cco_date <- clogit(case ~ discharge + strata(episode_id), 
                       data = cohort,
                       cluster = moh_study_id,
                       method  = "efron",
                       robust = T)

cco_adj_date <- clogit(case ~ discharge + homeless_hist + 
                         sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                         drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                         geq1_benzo + geq1_antipsyc + drug_toxicity + strata(episode_id), 
                       data = cohort,
                       cluster = moh_study_id,
                       method  = "efron",
                       robust = T)
summary(cco_adj_date)

saveRDS(cco_date,"cco_date.rds")
saveRDS(cco_adj_date,"cco_adj_date.rds")

####################### number of discharges ######################
cohort <- cohort %>% left_join(od_cohort_encounter %>% select(episode_id,od_date_tm) %>% distinct())
cohort_hosp <- cohort %>% 
  filter(interval == "od") %>% 
  left_join(dad,
            by = "moh_study_id",
            relationship = "many-to-many") %>% 
  filter(sep_date >= (date-days(574)) & sep_dt_tm < od_date_tm) %>% 
  select(episode_id,moh_study_id,od_date,sep_date,sep_disp) %>% 
  mutate(dc_bma = as.numeric(sep_disp %in% c(6,12,61,62,63,64,65))) %>% 
  mutate(diff_days = as.numeric(sep_date - od_date)) %>% 
  mutate(diff_weeks = ceiling(diff_days/7)) %>% 
  mutate(diff_months = ceiling(diff_days/30))

cnt_d <- cohort_hosp %>% 
  filter(diff_days >= -28) %>% 
  count(diff_days) %>% 
  left_join(cohort_hosp %>% filter(dc_bma == 1) %>% count(diff_days), by = "diff_days") %>% 
  rename("any_dc" = n.x, "dc_bma" = n.y) %>% 
  mutate("non dc_bma" = any_dc - dc_bma) %>% 
  mutate(pct = round(dc_bma/any_dc*100,2))


primary_cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")
primary_cohort_hosp <- primary_cohort %>% 
  filter(interval == "od") %>% 
  left_join(dad,
            by = "moh_study_id",
            relationship = "many-to-many") %>% 
  filter(sep_date >= (date-days(574)) & sep_date < date) %>% 
  select(episode_id,moh_study_id,od_date,sep_date,sep_disp) %>% 
  mutate(dc_bma = as.numeric(sep_disp %in% c(6,12,61,62,63,64,65))) %>% 
  mutate(diff_days = as.numeric(sep_date - od_date)) %>% 
  mutate(diff_weeks = ceiling(diff_days/7)) %>% 
  mutate(diff_months = ceiling(diff_days/30))
primary_cnt_d <- primary_cohort_hosp %>% 
  filter(diff_days >= -28) %>% 
  count(diff_days) %>% 
  left_join(primary_cohort_hosp %>% filter(dc_bma == 1) %>% count(diff_days), by = "diff_days") %>% 
  rename("any_dc" = n.x, "dc_bma" = n.y) %>% 
  mutate("non dc_bma" = any_dc - dc_bma) %>% 
  mutate(pct = round(dc_bma/any_dc*100,2))


write.xlsx(list("primary" = primary_cnt_d,
                "alternative" = cnt_d),
           "AMA_OD_cco od_date_diff.xlsx")




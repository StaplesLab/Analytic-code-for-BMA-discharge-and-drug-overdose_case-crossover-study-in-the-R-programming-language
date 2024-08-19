###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Daly-Grafstein D, Khan M, Erdelyi S, Staples JA. AMA-OD_cco sensitivity analyses - Create controls for case-crossover cohort. 2024 Feb 21. Retrieved from: ***LINK*** 
# Create controls for the AMA-OD case-crossover cohort in sensitivity analysis
# 2 Controls will be indexed 6 months and 12 months prior
# doesn't exclude od events with od in 18 months + 30 days prior to index od
# exclude od events earlier than July 28, 2016 (otherwise we don't know if they have prior od)
# also exclude those <18 years old at index od
# Author: Nicole Hu 
# Date: 2024-02-21
##########################################################

# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")

setwd("R:/working/AMA-OD_cco/NH/results/Table D")
dad <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_dad.csv.gz")
person_xwalk <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_person_xwalk.csv.gz")

od_cohort <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_episode_master.csv.gz") %>% 
  filter(od_flag == 1) %>% 
  mutate(od_date = as.Date(start_dt_tm)) #checked with od_date dataset

#od_cohort <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_od_date.csv.gz")
od_cohort <- od_cohort %>% mutate(control_date_1 = as.Date(od_date) %m-% days(182), #6 month
                                  control_date_2 = as.Date(od_date) %m-% days(364), #12 month
                                  control_date_3 = as.Date(od_date) %m-% days(546) #18 month
                                   ) %>%  
  select(-death_date,-death_date_source,-to_keep)
#ODCr
nrow(od_cohort) #N = 64156
od_cohort %>% select(moh_study_id) %>% n_distinct() #n = 36679

# exclude those <18 years old at any pre-control interval
od_filter <- od_cohort %>% 
  left_join(person_xwalk) %>% 
  mutate(control_age = (as.period(interval(birth_date, control_date_3 - days(28))))@year) %>% 
  mutate(younger_18 = as.numeric(control_age < 18 | control_age > 110))


#exclude those earlier than 2015-01-01 + 18 month + 28 days (2016-07-28)
od_filter <- od_filter %>% 
  mutate(earlier_16 = as.numeric(od_date %m-% (days(546) + days(28))) < as.Date("2015-01-01"))


# analytic data set for cco
od_cco <- od_filter %>% 
  #apply exclusion criteria
  filter(younger_18 == 0, earlier_16 == 0)


#total cohort
nrow(od_cco) #N = 48965
od_cco %>% select(moh_study_id) %>% n_distinct() #n = 28313

# make control as a row itself 
od_cco <- od_cco %>% 
  mutate(expose_od = NA, expose_control_1 = NA, expose_control_2 = NA) %>% 
  select(moh_study_id, od_date,episode_id,expose_od,expose_control_1,expose_control_2) %>% 
  mutate(id = seq(1:nrow(od_cco))) %>% 
  pivot_longer(c(expose_od, expose_control_1, expose_control_2)) %>% 
  mutate(date = case_when(
    name == "expose_od" ~ od_date, 
    name == "expose_control_1" ~ as.Date(od_date) %m-% days(182),
    name == "expose_control_2" ~ as.Date(od_date) %m-% days(364))) %>% 
  mutate(name = case_when(
    name == "expose_od" ~ "od", 
    name == "expose_control_1" ~ "control_1",
    name == "expose_control_2" ~ "control_2")) %>% 
  rename(interval = name) %>% #indicate intervals of case and control 1,2,3
  select(-value) %>% 
  mutate(interval2 = ifelse(interval == "od","od","control")) #whether the interval is case or control

# add exposure variable: based on most recent discharge within 30 days before od/control dates
od_cco <- od_cco %>% 
  left_join(dad,
            by = "moh_study_id",
            relationship = "many-to-many") %>% 
  filter(sep_date >= (date-days(28)) & sep_date < date) %>% 
  group_by(episode_id,date) %>% 
  arrange(desc(sep_date),desc(ad_date)) %>% 
  slice_head(n = 1) %>% #most recent discharge
  ungroup() %>% 
  select(episode_id,date,sep_disp) %>%
  distinct() %>% 
  mutate(discharge = case_when(sep_disp %in%  c(6,12,61,62,63,64,65) ~ "bma",
                               TRUE ~ "wa")) %>% 
  right_join(od_cco) %>% #join back to the whole cohort
  mutate(discharge = replace_na(discharge,"no discharge")) %>% 
  arrange(episode_id, date)

############ sex,age,location ##################################################
#sex&age
person_xwalk <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_person_xwalk.csv.gz")
cohort <- od_cco %>% 
  left_join(person_xwalk %>% select(-reference,-matching_id,-case_control_flag,-to_keep,od_death = fatal_od_case)) %>% 
  mutate(index_age = as.period(interval(birth_date, date))@year)
#location
od_event_location <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_od_event_location.csv.gz")
cohort <- cohort %>% 
  left_join(od_event_location %>% select(episode_id,od_event_fsa_16,od_event_lha_18))
#study month& study year
cohort <- cohort %>% 
  mutate(study_year = (year(date)  - year((ymd("2015-01-01"))) + 1)) %>% 
  mutate(study_month = (study_year - 1)*12 + month(date) - month((ymd("2015-01-01"))) + 1)

saveRDS(cohort, "cco_cohort_pod.rds")

######### covariates ############
cohort <- readRDS("cco_cohort_pod.rds")
cohort <- cohort %>%
  left_join(readRDS("cco_demo_pod.rds")) %>% 
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

#save cohort
saveRDS(cohort, "cco_coh_all_pod.rds") 

#################model#################
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table D/allow prior od/cco_coh_all_pod.rds")
unadj_model = clogit(case ~ discharge + strata(episode_id), 
               cohort,
               cluster = moh_study_id,
               method  = "efron",
               robust = T
)
adj_model = clogit(as.formula(paste("case ~", sub(".*~ ","",cco_adj$formula)[3])), 
             data = cohort,
             cluster = moh_study_id,
             method  = "efron",
             robust = T
)
bma_res <- 
  cbind(
    "bma",
    paste0(cohort %>% filter(interval == "od") %>% nrow(), " (", round((length(which(cohort$interval == "od"))*100)/27584,1) ,")"),
    cohort %>% filter(interval == "od",discharge == "bma") %>% nrow(),
    cohort %>% filter(interval == "control_1",discharge == "bma") %>% nrow(),
    cohort %>% filter(interval == "control_2",discharge == "bma") %>% nrow(),
    cohort %>% filter(interval2 == "control",discharge == "bma") %>% nrow(),
    tidy(unadj_model,exponentiate = T, conf.int = T)[1,] %>% 
      select(estimate,conf.low,conf.high,p.value) %>% 
      mutate(unadjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p", ifelse(p.value < 0.001, "<0.001", paste0("=",sprintf("%.3f",p.value))))) %>% 
      select(unadjusted),
    tidy(adj_model,exponentiate = T, conf.int = T)[1,] %>% 
      select(estimate,robust.se,conf.low,conf.high,p.value) %>% 
      mutate(adjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p", ifelse(p.value < 0.001, "<0.001", paste0("=",sprintf("%.3f",p.value))))) %>% 
      select(adjusted,estimate,robust.se,conf.low,conf.high)
  )
colnames(bma_res) <- c("Label","OD number","Exposure in pre-od","Exposure in control 1","Exposure in control 2","Exposure in pooled control","Unadjusted OR","Adjusted OR","estimate","robust.se","conf.low","conf.high")

wa_res <- 
  cbind(
    "wa",
    paste0(cohort %>% filter(interval == "od") %>% nrow(), " (", round((length(which(cohort$interval == "od"))*100)/27584,1) ,")"),
    cohort %>% filter(interval == "od",discharge == "wa") %>% nrow(),
    cohort %>% filter(interval == "control_1",discharge == "wa") %>% nrow(),
    cohort %>% filter(interval == "control_2",discharge == "wa") %>% nrow(),
    cohort %>% filter(interval2 == "control",discharge == "wa") %>% nrow(),
    tidy(unadj_model,exponentiate = T, conf.int = T)[2,] %>% 
      select(estimate,conf.low,conf.high,p.value) %>% 
      mutate(unadjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p", ifelse(p.value < 0.001, "<0.001", paste0("=",sprintf("%.3f",p.value))))) %>% 
      select(unadjusted),
    tidy(adj_model,exponentiate = T, conf.int = T)[2,] %>% 
      select(estimate,robust.se,conf.low,conf.high,p.value) %>% 
      mutate(adjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p", ifelse(p.value < 0.001, "<0.001", paste0("=",sprintf("%.3f",p.value))))) %>% 
      select(adjusted,estimate,robust.se,conf.low,conf.high)
  )
colnames(wa_res) <- c("Label","OD number","Exposure in pre-od","Exposure in control 1","Exposure in control 2","Exposure in pooled control","Unadjusted OR","Adjusted OR","estimate","robust.se","conf.low","conf.high")

bma_res
wa_res


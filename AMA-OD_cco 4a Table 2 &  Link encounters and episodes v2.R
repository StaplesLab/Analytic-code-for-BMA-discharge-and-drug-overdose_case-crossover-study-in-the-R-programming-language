###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Khan M, Daly-Grafstein D, Erdelyi S, Staples JA. AMA-OD_cco link overdose episodes with corresponding encounters. 2023 Oct 10. Retrieved from: ***LINK*** 
# Link overdose episodes with corresponding encounters and create table 2
# Author: Xiao (Nicole) Hu 
# Date: 2023-09-12
# Updated: 2023-10-10
##########################################################

# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")
setwd("R:/working/AMA-OD_cco/NH/results/Table B")
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")


#link episodes with encounters
ee_xwalk <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_episode_encounter_xwalk.csv.gz")
encounter_master <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_encounter_master.csv.gz")
od_date <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_od_date.csv.gz")

cohort_encounter <- cohort %>% 
  filter(interval == "od") %>%
  mutate(frst_od_age = (as.period(interval(birth_date, first_od_date)))@year) %>% 
  select(episode_id,od_date,moh_study_id,index_age,frst_od_age,od_death) %>% #od_death: if they died from od
  left_join(od_date %>% select(moh_study_id,od_date,od_fatal)) %>% #od_fatal: if the od episode is od
  left_join(ee_xwalk %>% select(-to_keep)) %>% 
  left_join(encounter_master %>% select(encounter_id,ec_start = start_dt_tm, ec_end = stop_dt_tm, ec_fatal = fatal)) %>% #ec_fatal: if the encounter is fatal
  group_by(episode_id) %>% 
  arrange(episode_id, ec_start, .by_group=TRUE) %>% 
  ungroup()

#proportion from each source
cohort_encounter %>% 
  filter(od_fatal == 1) %>% 
  group_by(source) %>% 
  summarise(encounter_n = n(),episode_n = n_distinct(episode_id)) %>% 
  mutate(pct = episode_n/length(unique(cohort_encounter$episode_id))) %>% 
  arrange(desc(pct))

## od death number is higher than fatal od events, because some fatal od events are excluded due to the "no prior od" criteria
cohort_encounter %>% filter(od_death == 1) %>% select(moh_study_id) %>% n_distinct()
cohort_encounter %>% filter(od_fatal == 1) %>% select(moh_study_id) %>% n_distinct()
cohort_encounter %>% filter(ec_fatal == 1) %>% select(moh_study_id) %>% n_distinct()

saveRDS(cohort_encounter,"cohort_encounter.rds")

#################### outcomes ##########################
#location
cohort_encounter <- readRDS("cohort_encounter.rds")
dpic <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_dpic.csv.gz")

location <- cohort_encounter %>% 
  filter(source == "DPIC") %>% 
  left_join(dpic %>% select(dpic_id,exp_site), by = c("source_id" = "dpic_id")) %>% 
  select(episode_id,exp_site) %>% 
  distinct()

#length of stay
dad <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_dad.csv.gz")
los <- cohort_encounter %>% 
  filter(source == "DAD") %>% 
  left_join(dad %>% select(dad_id,ad_date,sep_date), by = c("source_id" = "dad_id")) %>% 
  group_by(episode_id) %>% 
  summarise(
    ad_date = min(ad_date),
    sep_date = max(sep_date),
    los = as.duration(as.period(ad_date %--% sep_date))/ddays(1))

#outcome count number
cohort_encounter <- 
    cohort_encounter %>% mutate(hospitalized = as.numeric(source == "DAD"),
                                ed = as.numeric(source == "NACRS"),
                                bcehs = as.numeric(source =="EPCR" | source == "PCR"))
cohort_od_episode <- cohort_encounter %>% 
  group_by(episode_id) %>% 
  summarise(od_fatal = as.numeric(any(od_fatal == 1)),
            hospitalized = as.numeric(any(hospitalized == 1)),
            ed = as.numeric(any(ed == 1)),
            bcehs = as.numeric(any(bcehs == 1)))
cohort_od_episode <- cohort_od_episode %>% 
  mutate(outcome = case_when(od_fatal == 1 ~ "Fatal", #if fatal then fatal
                             hospitalized == 1 ~ "Taken to hospital", #if non-fatal and from DAD then hospitalized
                             ed == 1 ~ "Taken to ED", #if non-fatal and non-hospitalized and from NACRS then taken to ed
                             bcehs == 1 ~ "Treated on-site (EHS)", #if non-fatal and non-hospitalized and non-ed and from BCEHS then treated on-site
                             TRUE ~ "Treated on-site (EHS)"))  #there are only a few in other categories, consider them as treated on-site
cohort_od_episode %>% count(outcome) %>% arrange(desc(n)) %>% mutate(pct = round(n/nrow(cohort_od_episode)*100,2))

#first substance age
mha_client <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_mha_client.csv.gz")

frst_sbstc_age <- 
  mha_client %>% 
  select(moh_study_id,frst_alcl_use_age_in_yrs,frst_mrjn_use_age_in_yrs,frst_oth_use_age_in_yrs) %>% #n = 84226
  filter(moh_study_id %in% cohort$moh_study_id) %>%  #nrow = 20397
  mutate(across(-1,function(x) as.numeric(gsub("year\\w*",'',x)))) %>%
  pivot_longer(cols = c('frst_alcl_use_age_in_yrs','frst_mrjn_use_age_in_yrs','frst_oth_use_age_in_yrs'),
               names_to = 'sbstc_source',
               values_to = 'frst_sbstc_age') %>%
  filter(!is.na(frst_sbstc_age)) %>% 
  group_by(moh_study_id) %>%
  #the earliest age among three substance use columns and among multiple records of a person
  summarise(min_frst_sbstc_age = min(frst_sbstc_age, na.rm = T)) %>%  #nrow = 4606 
  ungroup() %>% 
  left_join(cohort %>% select(moh_study_id,episode_id),
             by = "moh_study_id") %>% 
  select(episode_id,min_frst_sbstc_age) %>% 
  distinct()
  

########## table 2 dataset ########################### 
tab2_df <- cohort %>% 
  filter(interval == "od") %>% 
  select(episode_id,age_group,male_sex,pop_density_class) %>% 
  left_join(cohort_od_episode) %>%
  left_join(location) %>% 
  replace(is.na(.), "Missing") %>% 
  left_join(frst_sbstc_age) %>% 
  left_join(los) %>% 
  mutate(los_missing = as.numeric(is.na(los))) %>% 
  mutate(frst_age_missing = as.numeric(is.na(min_frst_sbstc_age))) %>% 
  mutate(on_site = as.numeric(hospitalized == 0 & ed == 0))
saveRDS(tab2_df,"tab2_df.rds")

########## table 2 #################################
tab2_df <- readRDS("tab2_df.rds")
table2_vars <- c("age_group","male_sex", "pop_density_class","outcome","los","los_missing","exp_site","min_frst_sbstc_age","frst_age_missing","hospitalized","ed","bcehs","on_site")
non_norm_vars <- c("los","min_frst_sbstc_age")
cat_vars <- c("age_group","male_sex", "pop_density_class","outcome","exp_site","los_missing","frst_age_missing","hospitalized","ed","bcehs","on_site")

tab2a <- CreateTableOne(data = tab2_df, 
                        vars = table2_vars, factorVars = cat_vars)
tab2a_csv <- print(tab2a, nonnormal = non_norm_vars, showAllLevels = T, formatOptions = list(big.mark=','), printToggle = F,test = F)
write.csv(tab2a_csv, file = "Table 2a v2.csv")

# Group by exposure category
tab2b <- CreateTableOne(data = tab2_df,
                        strata = "od_fatal",
                        vars = table2_vars, factorVars = cat_vars)
tab2b_csv <- print(tab2b, nonnormal = non_norm_vars, showAllLevels = T, formatOptions = list(big.mark=','), smd=TRUE, printToggle = F)
write.csv(tab2b_csv, file = "Table 2b v3.csv")

tab2a <- read.csv("Table 2a v2.csv") %>% mutate(X = replace(X, X == "",NA)) %>% fill(X, .direction = "down")
tab2b <- read.csv("Table 2b v3.csv") %>%  mutate(X = replace(X, X == "",NA)) %>% fill(X, .direction = "down")

res <- tab2a %>% 
  left_join(tab2b, by = c("X","level")) 

colnames(res) <- c("Characteristics","Level","Fatal/non-fatal","Non-fatal","Fatal","p","test","SMD")
write.csv(res, file = "AMA-OD_cco Table 2 v4.csv")

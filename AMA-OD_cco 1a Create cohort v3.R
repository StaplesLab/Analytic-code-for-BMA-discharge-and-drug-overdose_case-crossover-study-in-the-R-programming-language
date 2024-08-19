###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Khan M, Daly-Grafstein D, Erdelyi S, Staples JA. AMA-OD_cco create controls for the AMA-OD case-crossover cohort. 2023 Dec 17. Retrieved from: ***LINK*** 
# Create controls for the AMA-OD case-crossover cohort
# Author: Xiao (Nicole) Hu 
# Date: 2023-09-08
# Updated: 2023-12-17
##########################################################

# 2 Controls will be indexed 6 months and 12 months prior
# exclude od events with od in 18 months + 30 days prior to index od
# exclude od events earlier than July 28, 2016 (otherwise we don't know if they have prior od)
# also exclude those <18 years old at index od
##########################################################

# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")

setwd("R:/working/AMA-OD_cco/NH/results/Table A")
dad <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_dad.csv.gz")
person_xwalk <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_person_xwalk.csv.gz")

od_cohort <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_episode_master.csv.gz") %>% 
  filter(od_flag == 1) %>% 
  mutate(od_date = as.Date(start_dt_tm)) #checked with od_date dataset

#od_cohort <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_od_date.csv.gz")
od_cohort <- od_cohort %>% mutate(control_date_1 = as.Date(od_date) %m-% days(182), #6 month
                                  control_date_2 = as.Date(od_date) %m-% days(364), #12 month
                                  control_date_3 = as.Date(od_date) %m-% days(546), #18 month
                                  control_date_lb6m = as.Date(od_date) %m+% days(182), #6 month after od
                                  control_date_lb12m = as.Date(od_date) %m+% days(364)) %>%  #12 month od od
              select(-death_date,-death_date_source,-to_keep)
#ODCr
nrow(od_cohort) #N = 64156
od_cohort %>% select(moh_study_id) %>% n_distinct() #n = 36679

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
nrow(od_cco) #N = 27584
od_cco %>% select(moh_study_id) %>% n_distinct() #n = 26589

# make control as a row itself 
od_cco <- od_cco %>% 
  mutate(expose_od = NA, expose_control_1 = NA, expose_control_2 = NA,expose_control_3 = NA,expose_control_lb6m = NA,expose_control_lb12m = NA) %>% 
  select(moh_study_id, od_date,episode_id,expose_od,expose_control_1,expose_control_2,expose_control_3,expose_control_lb6m,expose_control_lb12m) %>% 
  mutate(id = seq(1:nrow(od_cco))) %>% 
  pivot_longer(c(expose_od, expose_control_1, expose_control_2,expose_control_3,expose_control_lb6m,expose_control_lb12m)) %>% 
  mutate(date = case_when(
      name == "expose_od" ~ od_date, 
      name == "expose_control_1" ~ as.Date(od_date) %m-% days(182),
      name == "expose_control_2" ~ as.Date(od_date) %m-% days(364),
      name == "expose_control_3" ~ as.Date(od_date) %m-% days(546),
      name == "expose_control_lb6m" ~ as.Date(od_date) %m+% days(182),
      name == "expose_control_lb12m" ~ as.Date(od_date) %m+% days(364))) %>% 
  mutate(name = case_when(
      name == "expose_od" ~ "od", 
      name == "expose_control_1" ~ "control_1",
      name == "expose_control_2" ~ "control_2",
      name == "expose_control_3" ~ "control_3",
      name == "expose_control_lb6m" ~ "control_lb6m",
      name == "expose_control_lb12m" ~ "control_lb12m")) %>% 
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

############# flow chart ######################################################

#exclusion
exclusion_tab <- as_tibble(
  rbind(
    cbind("Overdose events with prior overdoses in the past 82 weeks",
          nrow(od_filter %>% filter(prior_od == 1) %>% select(moh_study_id) %>% distinct()), 
          sum(od_filter$prior_od),round(mean(od_filter$prior_od)*100,1)),#24930
    cbind("Overdose events earlier than July 28, 2016",
          nrow(od_filter %>% filter(earlier_16 == 1) %>% select(moh_study_id) %>% distinct()),
          sum(od_filter$earlier_16),round(mean(od_filter$earlier_16)*100,1)),#13122
    cbind("Patient age < 18 years at any control interval",
          nrow(od_filter %>% filter(younger_18 == 1) %>% select(moh_study_id) %>% distinct()),
          sum(od_filter$younger_18),round(mean(od_filter$younger_18)*100,1)),#2726
    cbind("Total exclusion",
          nrow(od_filter %>% filter(prior_od == 1|younger_18 == 1|earlier_16 == 1) %>% select(moh_study_id) %>% distinct()),
          nrow(od_filter %>% filter(prior_od == 1|younger_18 == 1|earlier_16 == 1)),
          round(nrow(od_filter %>% filter(prior_od == 1|younger_18 == 1|earlier_16 == 1))/nrow(od_cohort)*100,1)),
    cbind("Original cohort",
          od_cohort %>% select(moh_study_id) %>% n_distinct(), 
          nrow(od_cohort),NA),
    cbind("Total cohort",
          nrow(od_filter %>% filter(prior_od == 0, younger_18 == 0, earlier_16 == 0) %>% select(moh_study_id) %>% distinct()),
          nrow(od_filter %>% filter(prior_od == 0, younger_18 == 0, earlier_16 == 0)),NA)
  )
) %>% rename("Exclusion criteria" = V1,"n" = V2, "N" = V3, "Excluded percentage" = V4)

#flow chart
discharge_sum <- od_cco %>% 
  group_by(interval,discharge) %>% 
  summarise( n = n_distinct(moh_study_id),N = n()) %>% 
  mutate(pct = round(N/length(unique(od_cco$episode_id))*100,2))

#number of eligible overdoses per unique individual
cohort %>% 
  filter(interval == "od") %>% 
  count(moh_study_id) %>% 
  count(n) %>% 
  mutate(pct = (nn/nrow(cohort %>% 
                         filter(interval == "od") %>% select(moh_study_id) %>% distinct())) *100)

#number of unique hospitals
t <- cohort %>% select(episode_id,date,moh_study_id) %>% 
  left_join(dad,
            by = "moh_study_id",
            relationship = "many-to-many") %>% 
  filter(sep_date >= (date-days(28)) & sep_date < date) %>% 
  group_by(episode_id,date) %>% 
  arrange(desc(sep_date),desc(ad_date)) %>% 
  slice_head(n = 1) %>% #most recent discharge
  ungroup() %>% 
  select(episode_id,date,hosp,sep_disp) %>%
  mutate(discharge = case_when(sep_disp %in%  c(6,12,61,62,63,64,65) ~ "bma",
                               TRUE ~ "wa"))

t %>% filter(discharge == "wa") %>% select(hosp) %>% distinct() %>% nrow()
t %>% filter(discharge == "bma") %>% select(hosp) %>% distinct() %>% nrow()
t %>% select(hosp) %>% distinct() %>% nrow()

#save
write.xlsx(list("flow chart" = discharge_sum, 
                "exclusion" = exclusion_tab), "AMA-OD_cco - flow chart v4.xlsx") 

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

saveRDS(cohort, "cco_cohort.rds")



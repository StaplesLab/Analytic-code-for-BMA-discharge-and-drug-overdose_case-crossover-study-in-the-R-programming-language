###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Daly-Grafstein D, Khan M, Erdelyi S, Staples JA. AMA-OD_cco sensitivity analyses - categorized exposure variable. 2024 Feb 21. Retrieved from: ***LINK*** 
# Sensitivity analysis
# Codes exposure as: 
# - Physician-advised DC after psych hospitalization
# - BMA after psych hospitalization
# - Physician-advised DC after non-psych hospitalization
# - BMA after non-psych hospitalization
# - no discharge in interval (referent)

# Author: Nicole Hu 
# Date: 2024-02-21
##########################################################

source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")
setwd("R:/working/AMA-OD_cco/NH/results/Table D")
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")
cco_adj <- readRDS("R:/working/AMA-OD_cco/NH/results/Table C/cco_adj.rds")
dad <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_dad.csv.gz")


cohort <- cohort %>% left_join(dad,
                           by = "moh_study_id",
                           relationship = "many-to-many") %>% 
  filter(sep_date >= (date-days(28)) & sep_date < date) %>% 
  group_by(episode_id,date) %>% 
  arrange(desc(sep_date),desc(ad_date)) %>% 
  slice_head(n = 1) %>% #most recent discharge
  ungroup() %>% 
  select(episode_id,date,diagx_1) %>%
  distinct() %>% 
  mutate(psych_hosp = as.numeric(str_detect(diagx_1, 'F'))) %>% 
  right_join(cohort) %>% #join back to the whole cohort
  mutate(discharge_long = case_when(dc_bma == 1 & psych_hosp == 1 ~ "bma_psych",
                                    dc_wa == 1 & psych_hosp == 1 ~ "wa_psych",
                                    dc_bma == 1 & psych_hosp == 0 ~ "bma_nonpsych",
                                    dc_wa == 1 & psych_hosp == 0 ~ "wa_nonpsych",
                                    TRUE ~ "no discharge")) %>% 
  arrange(episode_id, date)
cohort$discharge_long <- relevel(as.factor(cohort$discharge_long), ref = "no discharge")

#
cnt <- cohort %>% count(discharge_long,interval)

#
unadj = clogit(case ~ discharge_long + strata(episode_id), 
               data = cohort,
               cluster = moh_study_id,
               method  = "efron",
               robust = T)
adj = clogit(as.formula(gsub(paste("discharge","+"),"", paste("case ~ discharge_long", sub(".*~ ","",cco_adj$formula)[3]))), 
             data = cohort,
             cluster = moh_study_id,
             method  = "efron",
             robust = T)

unadj_coef <- tidy(unadj,exponentiate = T, conf.int = T) %>% 
  select(term,estimate,conf.low,conf.high,p.value) %>% 
  mutate(unadjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p", ifelse(p.value < 0.001, "<0.001", paste0("=",sprintf("%.3f",p.value))))) %>% 
  select(term,unadjusted)


adj_coef <- tidy(adj,exponentiate = T, conf.int = T) %>% 
  filter(str_detect(term,"wa")
         | str_detect(term,"bma")) %>% 
  select(term,estimate,conf.low,conf.high,p.value) %>% 
  mutate(adjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p", ifelse(p.value < 0.001, "<0.001", paste0("=",sprintf("%.3f",p.value))))) %>% 
  select(term,adjusted)

write.xlsx(list("Number of exposure" = cnt,
                "Unadjusted OR" = unadj_coef,
                "Adjsuted OR" = adj_coef), "R:/working/AMA-OD_cco/NH/results/Table D/AMA-OD_cco 5level_exposure.xlsx")


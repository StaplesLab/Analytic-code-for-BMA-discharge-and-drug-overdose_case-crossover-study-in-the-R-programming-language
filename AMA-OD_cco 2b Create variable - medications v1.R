###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Khan M, Daly-Grafstein D, Erdelyi S, Staples JA. AMA-OD_cco create medication variables. 2023 Sep 8. Retrieved from: https://github.com/StaplesLab/Analytic-code-for-BMA-discharge-and-drug-overdose_case-crossover-study-in-the-R-programming-language/edit/main/AMA-OD_cco%202b%20Create%20variable%20-%20medications%20v1.R 
# Create medication variables
# Author: Xiao (Nicole) Hu 
# Date: 2023-09-08
# Updated: 2023-09-08
##########################################################
# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")

#read in cohort
setwd("R:/working/AMA-OD_cco/NH/results/Table A")
cohort <- readRDS("cco_cohort.rds")

# cohort for sensitivity analyses 
setwd("R:/working/AMA-OD_cco/NH/results/Table D/allow prior od")
#cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table D/cco_cohort_pod.rds") #cohort including od with prior od
#cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table D/nod_coh.rds") #non overdose cohort
#cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table D/od_date/cco_cohort_date.rds") #cohort with alternative overdose date


#read in
f <- function(x,pos){
  subset(x, 
         select = c(
           'moh_study_id',
           'din_pin',
           "dispensing_claims_srv_date",
           "dispensing_claims_qty",
           "dispensing_claims_days_sply",
           "dispensing_claims_sply_end_date"),
         moh_study_id %in% cohort$moh_study_id
  )
}

read_chunk <- function(chunk) {
  df <- read_csv_chunked(chunk,
                         callback =  DataFrameCallback$new(f),
                         chunk_size = 1000000) 
  return(df)
}

path_rx <- "R:/DATA/2023-05-09_190103/POC2023-04-05/vw_pnet_rx_group.csv.gz"
pnet_rx_group <- read_chunk(path_rx) 
#saveRDS(pnet_rx_group,"R:/working/AMA-OD_cco/NH/results/Table D/cohort_pnet_rx.rds") 


path_other <- "R:/DATA/2023-05-09_190103/POC2023-04-05/vw_pnet_other.csv.gz"
pnet_other <- read_chunk(path_other)
#saveRDS(pnet_other, "cohort_pnet_other.rds")



pnet <- pnet_rx_group %>% bind_rows(pnet_other) 
#saveRDS(pnet, "cohort_pnet.rds")
#saveRDS(pnet, "R:/working/AMA-OD_cco/NH/results/Table D/allow prior od/nod-cohort_pnet.rds")
#saveRDS(pnet, "R:/working/AMA-OD_cco/NH/results/Table D/od_date/date-cohort_pnet.rds")
#saveRDS(pnet, "R:/working/AMA-OD_cco/NH/results/Table D/od_date/pod-cohort_pnet.rds")

pnet <- readRDS("cohort_pnet.rds")

#number of prescription medications at baseline 
pharm_active <- pnet %>%   
  inner_join(cohort %>% select(episode_id, moh_study_id, date),
             by = "moh_study_id",
             relationship = "many-to-many") %>% 
  filter((date > dispensing_claims_srv_date) & 
           (date <= (dispensing_claims_srv_date + dispensing_claims_days_sply))) %>% 
  select(episode_id, date, din_pin) %>%
  distinct() %>%
  group_by(episode_id,date) %>%
  summarise(num_active = n()) %>%
  ungroup() %>%
  right_join(cohort %>% select(episode_id, date),
             by = c("episode_id","date")) %>%
  mutate(num_active = replace_na(num_active, 0)) %>%
  select(episode_id, date, num_active)

saveRDS(pharm_active,"ttl_meds_ad.rds")

#active OAT prescription at index date
lu_pnet_dinpin_rx_group <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_lu_pnet_dinpin_rx_group.csv.gz")
oat_active <- pnet %>% 
  inner_join(lu_pnet_dinpin_rx_group,
             by = "din_pin") %>% 
  #get all oat prescriptions
  filter(oat == 1) %>%
  #join with the cohort to get everyone's oat prescriptions
  inner_join(cohort %>% select(episode_id,moh_study_id, date),
             by = "moh_study_id",
             relationship = "many-to-many") %>% 
  #include active prescriptions at admission date
  filter((date >= dispensing_claims_srv_date) & 
           (date <= (dispensing_claims_srv_date + dispensing_claims_days_sply + 1))) %>%
  #only care if there is a oat prescription active at index date
  select(episode_id, date, oat) %>%
  distinct() %>% 
  right_join(cohort %>% select(episode_id, date),
             by = c("episode_id","date")) %>% 
  mutate(oat_active_ad = as.numeric(!is.na(oat))) %>% 
  select(episode_id, date, oat_active_ad)

saveRDS(oat_active, "oat_active_ad.rds") #update: 0626


#selected prescription medication filled in the 90d prior to index date
lu_pnet_dinpin_rx_group <- lu_pnet_dinpin_rx_group %>% 
  mutate(antipsyc = as.numeric(grepl("antipsychotic", ahfs_3, ignore.case = T))) %>%
  mutate(opioid = as.numeric(!is.na(opioid_category))) %>% 
  rowwise() %>% 
  mutate(other = as.numeric(sum(c(oat, benzo, opioid, antipsyc, z_drug),na.rm = T) == 0))

sub_meds_90ad <- pnet %>%   
  inner_join(cohort %>% select(episode_id,moh_study_id, date),
             by = "moh_study_id",
             relationship = "many-to-many") %>% 
  filter((date > dispensing_claims_srv_date) & 
           (date <= (dispensing_claims_srv_date + 90))) %>% 
  select(episode_id, date, din_pin) %>% 
  distinct() %>% 
  mutate(oat = as.numeric(din_pin %in% as.numeric(lu_pnet_dinpin_rx_group %>% filter(oat==1) %>% pull(din_pin))),
         benzo = as.numeric(din_pin %in% as.numeric(lu_pnet_dinpin_rx_group %>% filter(benzo==1 | z_drug==1) %>% pull(din_pin))),
         opioid = as.numeric(din_pin %in% as.numeric(lu_pnet_dinpin_rx_group %>% filter(opioid == 1) %>% filter(oat == 0 | is.na(oat)) %>% pull(din_pin))),
         antipsyc = as.numeric(din_pin %in% as.numeric(lu_pnet_dinpin_rx_group %>% filter(antipsyc==1) %>% pull(din_pin))),
         other = as.numeric(din_pin %in% as.numeric(lu_pnet_dinpin_rx_group %>% filter(other==1) %>% pull(din_pin))|
                              !(din_pin %in% as.numeric(lu_pnet_dinpin_rx_group %>% pull(din_pin))))) %>% 
  select(-din_pin) %>% 
  group_by(episode_id, date) %>% 
  summarise(across(everything(), sum)) %>% 
  ungroup() %>% 
  right_join(cohort %>% select(episode_id,date),
             by = c("episode_id","date")) %>% 
  replace(is.na(.), 0)

saveRDS(sub_meds_90ad, "sub_meds_90ad.rds") 

#oat in prior 6 months
oat_6m <- pnet %>%   
  inner_join(cohort %>% select(episode_id,moh_study_id, date),
             by = "moh_study_id",
             relationship = "many-to-many") %>% 
  filter(dispensing_claims_srv_date < date & 
           dispensing_claims_srv_date >= date %m-% days(182)) %>% 
  select(episode_id, date, din_pin) %>% 
  distinct() %>% 
  mutate(oat_6m = as.numeric(din_pin %in% as.numeric(lu_pnet_dinpin_rx_group %>% filter(oat==1) %>% pull(din_pin)))) %>% 
  select(-din_pin) %>% 
  group_by(episode_id, date) %>% 
  summarise(across(everything(), sum)) %>% 
  ungroup() %>% 
  right_join(cohort %>% select(episode_id,date),
             by = c("episode_id","date")) %>% 
  replace(is.na(.), 0) %>% 
  mutate(geq1_oat_6m = as.numeric(oat_6m > 0))

saveRDS(oat_6m,"oat_6m.rds")


#save results
medication <- readRDS("sub_meds_90ad.rds") %>% 
  full_join(readRDS("oat_active_ad.rds")) %>% 
  full_join(readRDS("ttl_meds_ad.rds")) %>% 
  #saveRDS("date_cco_medication.rds") %>% 
  #saveRDS("nod_cco_medication.rds")
  saveRDS("cco_medication.rds")
  

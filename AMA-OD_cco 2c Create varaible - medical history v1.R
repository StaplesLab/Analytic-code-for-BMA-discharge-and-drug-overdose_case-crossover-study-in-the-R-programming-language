###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Khan M, Daly-Grafstein D, Erdelyi S, Staples JA. AMA-OD_cco create medical history variables. 2023 Sep 20. Retrieved from: ***LINK*** 
# Create medical history variables
# Author: Xiao (Nicole) Hu 
# Date: 2023-09-08
# Updated: 2023-09-08
##########################################################
# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")

setwd("R:/working/AMA-OD_cco/NH/results/Table A")
cohort <- readRDS("cco_cohort.rds")

#setwd("R:/working/AMA-OD_cco/NH/results/Table D")
#cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table D/nod_coh.rds") #non overdose cohort
#cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table D/od_date/cco_cohort_date.rds") #cohort with alternative overdose date

################## DAD ###############################
dad <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_dad.csv.gz")
# get all hospitalizations of cohort
hosp <- dad %>% 
  filter(moh_study_id %in% cohort$moh_study_id) %>% 
  select(moh_study_id, dad_id, ad_date, starts_with("sep"), t_days, starts_with("diagx"),
         # Also select hospto, hospfrom, admit variables
         hosp_to,hosp_from, admit, hosp)

# creat los
hosp <- hosp %>%
  mutate(sep_date = ymd(sep_date), 
         ad_date = ymd(ad_date), 
         tlos = as.period(ad_date %--% sep_date),
         los = as.duration(tlos) / ddays(1)) %>% 
  select(-tlos)

saveRDS(hosp,"cohort_hosp.rds")

hosp <- readRDS("cohort_hosp.rds")
# prior hospitalization
prior_hosp <- hosp %>% 
  inner_join(., cohort, by = c("moh_study_id"),relationship = "many-to-many") %>%
  # if the event happened in the interval of [date - 365, date - 1]
  filter(sep_date <= date %m-% days(1) &
           sep_date >= date %m-% days(182)) %>%
  group_by(episode_id,date) %>% 
  summarise(num_hosp_prev_6m=n(), tdays_prev_6m=sum(los)) %>% 
  ungroup() %>% 
  right_join(cohort, by = c("episode_id","date")) %>%
  mutate(num_hosp_prev_6m = replace_na(num_hosp_prev_6m, 0),
         tdays_prev_6m = replace_na(tdays_prev_6m, 0),
         geq1_hosp_prev_6m = ifelse(num_hosp_prev_6m > 0, 1, 0)) %>%
  select(episode_id,date, num_hosp_prev_6m, geq1_hosp_prev_6m, tdays_prev_6m)
saveRDS(prior_hosp,"prior_hosp.rds") 

# Psychiatric hospitalizations in the prior 6 months
prior_psyc_hosp <- cohort %>% 
  #pychiatric hospitalizations are diagx_1 starts with "F"
  left_join(hosp %>% filter((str_detect(diagx_1, 'F'))), 
            by = "moh_study_id", relationship = "many-to-many") %>%
  filter(sep_date <= date %m-% days(1) &
           sep_date >= date %m-% days(182)) %>%
  group_by(episode_id,date) %>% 
  summarise(num_psyc_hosp_prev_6m=n()) %>%
  right_join(cohort, by = c("episode_id","date")) %>%
  mutate(num_psyc_hosp_prev_6m = replace_na(num_psyc_hosp_prev_6m, 0),
         geq1_psyc_hosp_prev_6m = ifelse(num_psyc_hosp_prev_6m > 0, 1, 0)) %>%
  select(episode_id,date, num_psyc_hosp_prev_6m,geq1_psyc_hosp_prev_6m) %>% 
  ungroup()
saveRDS(prior_psyc_hosp, "prior_psyc_hosp.rds")

# DAD Comorbidities
icd10_ca <- readxl::read_xlsx("R:/working/XW files/Staples Lab - PMH XW ICD-10-CA - 2023-07-25.xlsx")
icd10_ca <- icd10_ca %>% mutate(pmh_name = ifelse(pmh_name == "psychosis","psychotic",pmh_name ))
selected_coms <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/selected_coms.rds")

# For the LHS, duplicate the variables used for the join; rename other variables (eg DIAGX_all, TDAYS) to make clear their provenance
hosp <- hosp %>% 
  unite("diagx_all", contains("diagx"), na.rm = TRUE, sep = "_, _") %>% 
  mutate(diagx_all = paste0("_", diagx_all, "_")) %>% 
  dplyr::select(moh_study_id, pmh_sepdate = sep_date, pmh_sepdate_merge = sep_date, pmh_diagx_all = diagx_all, pmh_tdays = los) %>% 
  setDT()

# For the RHS, create the dates used for the join and duplicate these variables
coh_ids <- cohort %>%
  mutate(date_m1d = date %m-% days(1),
         date_m6m = date %m-% days(182)) %>% 
  dplyr::select(episode_id,date,moh_study_id, start = date_m6m, end = date_m1d) %>%
  setDT()

# Complete a right_join (ie. keep all rows from coh_ids)
pmh_dad <- hosp[coh_ids,
                on = .(moh_study_id == moh_study_id,
                       pmh_sepdate_merge >= start,
                       pmh_sepdate_merge <= end)] %>%
  dplyr::select(episode_id,date, pmh_sepdate, pmh_diagx_all, pmh_tdays)


# Map through the list of selected comorbidities
pmh_dad <- map(selected_coms, ~ if_else(str_detect(pmh_dad$pmh_diagx_all, 
                                                   icd10_ca %>% 
                                                     filter(pmh_name == .x) %>% 
                                                     pull(icd10_code) %>% 
                                                     unique() %>% 
                                                     paste0(., collapse = "|")), 
                                        1, 
                                        0)) %>% 
  set_names(paste0(selected_coms, ".dad")) %>% 
  bind_cols(pmh_dad, .)


pmh_dad <-   
  pmh_dad %>%
  mutate(dm_nc_or_compl.dad = case_when(dm_nc.dad == 1 | dm_compl.dad == 1 ~ 1, TRUE ~ 0),
         mild_mod_sev_liver.dad = case_when(liver_mild.dad == 1 | liver_modsev.dad == 1 ~ 1, TRUE ~ 0),
         cancer_or_mc.dad = case_when(cancer.dad == 1 | mets.dad == 1 ~ 1, TRUE ~ 0),
         any_sbstcs.dad = case_when(alcohol.dad == 1 | drugs_all.dad == 1 | drugs_opioid.dad == 1 | drugs_nonopioid.dad == 1 ~ 1, TRUE ~ 0)) %>% 
  #if multiple comorbidities of a given group in 1 ADDATE, just count as 1
  group_by(episode_id,date,pmh_sepdate) %>% summarise(across(.cols = ends_with("dad"),
                                                                              function(x) {ifelse(sum(x) >= 1, 1, 0)})) %>%
  ungroup() %>%
  # sum each comorbidity for each patient
  group_by(episode_id,date) %>%
  # capture the number of visits for a given comorbidity
  summarise(across(.cols = ends_with("dad"),
                   .fns = ~sum(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  replace(is.na(.), 0) #diag = NA when there is no prior hosp

saveRDS(pmh_dad, "pmh_dad.rds")

################## MSP ###############################
# read in cohort msp
f <- function(x,pos){
  subset(x, 
         select = c(moh_study_id, msp_id, serv_date,serv_loc, 
                    fitm,clm_spec,diag_cd_1,diag_cd_2,diag_cd_3),
         moh_study_id %in% cohort$moh_study_id)
}

read_chunk <- function(chunk) {
  df <- read_csv_chunked(chunk,
                         callback =  DataFrameCallback$new(f),
                         chunk_size = 1000000) 
  return(df)
}

path <- "R:/DATA/2023-05-09_190103/POC2023-04-05/vw_msp.csv.gz"
cohort_msp <- read_chunk(path) 
#saveRDS(cohort_msp, "cohort_msp.rds")
#saveRDS(cohort_msp, "cohort_msp_pod.rds")


#Physician clinic visits in prior 6m
cohort_msp <- readRDS("cohort_msp.rds")
MD_spec_XW <- readxl::read_xlsx("R:/working/XW files/Staples Lab - MD specialties included for 'prior clinic visits' - 2023-08-08.xlsx",sheet = "Sheet1")
MD_spec_XW$Code <- as.numeric(MD_spec_XW$Code)
prior_clinic <- cohort_msp  %>%
  filter(serv_loc == "A") %>%
  left_join(MD_spec_XW, by = c("clm_spec" = "Code")) %>% 
  filter(incl_for_outpt_visits == 1) %>% 
  # msp table rows represent billing items, so only count unique dates as one visit
  select(moh_study_id, serv_date) %>% 
  distinct() %>% 
  inner_join(cohort, by = "moh_study_id",relationship = "many-to-many") %>%
  filter(as.Date(serv_date) >= date %m-% days(182)  & 
           as.Date(serv_date) <= date %m-% days(1) ) %>% 
  group_by(episode_id, date) %>%
  summarise(num_clinic_prev_6m = n()) %>%
  ungroup() %>%
  right_join(cohort, by = c("episode_id", "date")) %>%
  mutate(num_clinic_prev_6m = replace_na(num_clinic_prev_6m, 0),
         geq7_clinic_prev_6m = ifelse(num_clinic_prev_6m >= 7, 1,0)) %>%
  select(episode_id, date, num_clinic_prev_6m, geq7_clinic_prev_6m)
saveRDS(prior_clinic, "prior_clinic.rds") 

# median(prior_clinic$num_clinic_prev_6m)
# hist(prior_clinic$num_clinic_prev_6m)
# prior_clinic %>% count(geq7_clinic_prev_6m) %>% mutate(p = n/nrow(cohort))


#MSP Comorbidities
icd9_cm <- readxl::read_xlsx("R:/working/XW files/Staples Lab - PMH XW ICD-9-CM - 2023-07-25.xlsx")
# For the LHS, duplicate the variables used for the join; rename other variables (eg DIAGX_all, TDAYS) to make clear their provenance
msp <- cohort_msp %>% 
  unite("icd_all", contains("diag_cd"), na.rm = TRUE, sep = "_, _") %>% 
  mutate(icd_all = paste0("_", icd_all, "_")) %>% 
  select(moh_study_id, pmh_servdate = serv_date, pmh_servdate_merge = serv_date, pmh_icd_all = icd_all) %>% 
  setDT()

# For the RHS, create the dates used for the join and duplicate these variables
coh_ids <- cohort %>%
  mutate(date_m1d = date %m-% days(1),
         date_m6m = date %m-% days(182)) %>% 
  select(episode_id, date, moh_study_id, start = date_m6m, end = date_m1d) %>%
  setDT()

# Complete a right_join (ie. keep all rows from coh_ids)
pmh_msp <- msp[coh_ids,
               on = .(moh_study_id == moh_study_id,
                      pmh_servdate_merge >= start,
                      pmh_servdate_merge <= end)] %>%
  select(episode_id, date, pmh_servdate, pmh_icd_all)

# For multiple visits at the same serv_date, combine all icd-codes 
pmh_msp <- pmh_msp %>% 
  distinct() %>% 
  group_by(episode_id, date, pmh_servdate) %>% 
  summarise(pmh_icd_all = str_c(pmh_icd_all, collapse = ", ")) %>% 
  ungroup()

# Map through the list of selected comorbidities
pmh_msp <- map(selected_coms, ~ if_else(str_detect(pmh_msp$pmh_icd_all, 
                                                   icd9_cm %>% 
                                                     filter(pmh_name == .x) %>% 
                                                     pull(icd9_code) %>% 
                                                     unique() %>% 
                                                     paste0("_",.,"_", collapse = "|")), #when icd9_code = 3040, it wouldn't be detected as 304
                                        1, 
                                        0)) %>% 
  set_names(paste0(selected_coms, ".msp")) %>% 
  bind_cols(pmh_msp, .)


pmh_msp <-   
  pmh_msp %>%
  mutate(dm_nc_or_compl.msp = case_when(dm_nc.msp == 1 | dm_compl.msp == 1 ~ 1, TRUE ~ 0),
         mild_mod_sev_liver.msp = case_when(liver_mild.msp == 1 | liver_modsev.msp == 1 ~ 1, TRUE ~ 0),
         cancer_or_mc.msp = case_when(cancer.msp == 1 | mets.msp == 1 ~ 1, TRUE ~ 0),
         any_sbstcs.msp = case_when(alcohol.msp == 1 | drugs_all.msp == 1 | drugs_opioid.msp == 1 | drugs_nonopioid.msp == 1 ~ 1, TRUE ~ 0)) %>% 
  # sum each comorbidity for each patient
  group_by(episode_id, date) %>%
  # capture the number of visits for a given comorbidity
  summarise(across(.cols = ends_with("msp"),
                   .fns = ~sum(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  replace(is.na(.), 0)


saveRDS(pmh_msp,"pmh_msp.rds") 

################# Charlson's score ##############################################
pmh_all <- readRDS("pmh_dad.rds") %>% 
  full_join(readRDS("pmh_msp.rds")) 

# Create function to identify comorbidities in the 6m prior to admission (1 hosp or 2 physician services)
get_comorbidity_core <- function(com_list, df) {
  for (i in 1:length(com_list)) {
    name <- com_list[i]
    df[name] <- as.numeric(df[paste0(name,'.dad')] >= 1 | df[paste0(name,'.msp')] >= 2)} 
  return(df)}

pmh_all <- pmh_all %>% 
  get_comorbidity_core(com_list = gsub(".dad","",colnames(readRDS("pmh_dad.rds"))[-c(1:2)]), .) %>%
  mutate(
    # Charlson comorbidity score
    cci = (mi*1 + chf*1 + pvd*1 + cvd*1 + dem*1 + copd*1 + rheum*1 + pud*1 + liver_mild*1 + dm_nc*1 + 
             dm_compl*2 + paraplegia*2 + renal*2 + cancer*2 + liver_modsev*3 + mets*6 + hiv*6),
    cci_cat = if_else(cci>=2, 1, 0) %>% factor(., levels = c(0, 1))) 

saveRDS(pmh_all, "pmh_all.rds")

################## homelessness history #################
mha_dsm_5 <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_mha_dsm_5.csv.gz")
homeless_mha <- 
  cohort %>% 
  select(episode_id,moh_study_id,date) %>% 
  left_join(mha_dsm_5,
            by = "moh_study_id",
            relationship = "many-to-many") %>% 
  filter(as.Date(enrl_assmt_date) >= date %m-% days(182) &
           as.Date(enrl_assmt_date) <= date %m-% days(1)) %>% 
  mutate(homeless = as.numeric(if_any(paste0("dsm_5_enrl_dx_cd_", 1:3), ~ . %in% c("Z59.0","Z59.1")))) %>% 
  group_by(episode_id, date) %>% 
  summarise(homeless.mha = sum(homeless)) %>% 
  right_join(cohort, by = c("episode_id", "date"))  %>% 
  mutate(homeless.mha = replace_na(homeless.mha, 0)) %>% 
  ungroup() %>% 
  select(episode_id, date,homeless.mha)

homeless_hist <- readRDS("pmh_all.rds") %>% select(episode_id,date,homeless.dad) %>% 
  full_join(homeless_mha) %>% 
  mutate(homeless_hist = as.numeric(homeless.dad >= 1 | homeless.mha >= 1)) %>% 
  select(episode_id, date, homeless_hist)

saveRDS(homeless_hist, "homeless_hist.rds")

################### iv_drug use #########################
ivdu_hist <- readRDS("pmh_all.rds") %>% select(episode_id,date,ivdu) %>% 
  left_join(readRDS("oat_6m.rds") %>% select(episode_id,date,geq1_oat_6m)) %>% 
  mutate(ivdu_hist = as.numeric(ivdu == 1 | geq1_oat_6m == 1)) %>% 
  distinct() %>% 
  select(episode_id, date, ivdu_hist)
saveRDS(ivdu_hist, "ivdu_hist.rds")


#save results
medical_hist <- readRDS("prior_hosp.rds") %>% 
  full_join(readRDS("prior_psyc_hosp.rds")) %>% 
  full_join(readRDS("prior_clinic.rds")) %>% 
  full_join(readRDS("pmh_all.rds") %>% select(- ends_with("dad"),- ends_with("msp"), -ivdu, -homeless)) %>% 
  full_join(readRDS("ivdu_hist.rds")) %>% 
  full_join(readRDS("homeless_hist.rds"))

saveRDS(medical_hist,"cco_medical_hist.rds")
#saveRDS(medical_hist,"date_cco_medical_hist.rds")
#saveRDS(medical_hist,"nod_cco_medical_hist.rds")

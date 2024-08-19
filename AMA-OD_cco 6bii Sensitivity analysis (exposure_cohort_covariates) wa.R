###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X,  Yu Y, Daly-Grafstein D, Khan M, Erdelyi S, Staples JA. AMA-OD_cco sensitivity analyses - Part IV. 2023 Nov 16. Retrieved from: ***LINK*** 
# Sensitivity analysis (2)
# alternate exposure
# alternate cohort
# alternate covariates
# add interaction term: male_sex * dc_bma
# analysis unit: one random overdose for each individual
# analysis unit: earliest eligible overdose for each individual
# Author: Xiao (Nicole) Hu 
# Date: 2023-09-29
# Updated: 2023-11-16
##########################################################
# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")
setwd("R:/working/AMA-OD_cco/NH/results/Table D")
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")
cco_adj <- readRDS("R:/working/AMA-OD_cco/NH/results/Table C/cco_adj.rds")
dad <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_dad.csv.gz")

#get number of overdose, number of exposure, odds ratio of DC-BMA
get_pasted_col <- function(df,unadj_model, adj_model, label){
  res <- 
    cbind(
      label,
      paste0(df %>% filter(interval == "od") %>% nrow(), " (", round((length(which(df$interval == "od"))*100)/27584,1) ,")"),
      df %>% filter(interval == "od",discharge == "wa") %>% nrow(),
      df %>% filter(interval == "control_1",discharge == "wa") %>% nrow(),
      df %>% filter(interval == "control_2",discharge == "wa") %>% nrow(),
      df %>% filter(interval2 == "control",discharge == "wa") %>% nrow(),
      tidy(unadj_model,exponentiate = T, conf.int = T) %>% 
        filter(str_detect(term,"wa")) %>% 
        select(estimate,conf.low,conf.high,p.value) %>% 
        mutate(unadjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p = ", ifelse(p.value < 0.001, "<0.001", sprintf("%.3f",p.value)))) %>% 
        select(unadjusted),
      tidy(adj_model,exponentiate = T, conf.int = T) %>% 
        filter(str_detect(term,"wa")) %>% 
        select(estimate,conf.low,conf.high,p.value) %>% 
        mutate(adjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p = ", ifelse(p.value < 0.001, "<0.001", sprintf("%.3f",p.value)))) %>% 
        select(adjusted)
    )
  colnames(res) <- c("Label","OD number","Exposure in pre-od","Exposure in control 1","Exposure in control 2","Exposure in pooled control","Unadjusted OR","Adjusted OR")
  return(res)
}

############# cohort = od+hosp in pre-od and pre-control, control = t0 - 1y, and only adjusted for dc_wa ############
hosp_cohort_1y <- cohort %>% 
  filter(interval != "control_1") %>% 
  left_join(dad, by = "moh_study_id", relationship = "many-to-many") %>% 
  filter(sep_date >= date-days(28) & sep_date < date) %>% 
  select(episode_id,date) %>%
  distinct() %>%
  count(episode_id) %>% 
  filter(n == 2) %>% 
  select(-n) %>% 
  left_join(cohort %>% filter(interval != "control_1"))

unadj = clogit(case ~ dc_wa + strata(episode_id), 
               hosp_cohort_1y,
               cluster = moh_study_id,
               method  = "efron",
               robust = T)
adj = clogit(as.formula(gsub(paste("discharge","+"),"", paste("case ~ dc_wa", sub(".*~ ","",cco_adj$formula)[3]))), 
             data = hosp_cohort_1y,
             cluster = moh_study_id,
             method  = "efron",
             robust = T)
res_hosp_c_1y <- get_pasted_col(hosp_cohort_1y,unadj,adj,"hosp_cohort_1y")

#cohort = od+hosp in pre-od and pre-control, control = t0 - 6m
hosp_cohort_6m <- cohort %>% 
  filter(interval != "control_2") %>% 
  left_join(dad, by = "moh_study_id", relationship = "many-to-many") %>% 
  filter(sep_date >= date-days(28) & sep_date < date) %>% 
  select(episode_id,date) %>%
  distinct() %>%
  count(episode_id) %>% 
  filter(n == 2) %>% 
  select(-n) %>% 
  left_join(cohort %>% filter(interval != "control_2"))
unadj = clogit(case ~ dc_wa + strata(episode_id), 
               hosp_cohort_6m,
               cluster = moh_study_id,
               method  = "efron",
               robust = T)
adj = clogit(as.formula(gsub(paste("discharge","+"),"", paste("case ~ dc_wa", sub(".*~ ","",cco_adj$formula)[3]))), 
             data = hosp_cohort_6m,
             cluster = moh_study_id,
             method  = "efron",
             robust = T)
res_hosp_c_6m <- get_pasted_col(hosp_cohort_6m,unadj,adj,"hosp_cohort_6m")


##################### use dc_wa as exposure #########################
unadj = clogit(case ~ dc_wa + strata(episode_id),
               cohort,
               cluster = moh_study_id,
               method  = "efron",
               robust = T
)
adj = clogit(as.formula(gsub(paste("discharge","+"),"", paste("case ~ dc_wa", sub(".*~ ","",cco_adj$formula)[3]))),
             data = cohort,
             cluster = moh_study_id,
             method  = "efron",
             robust = T)

res_dcwa <- get_pasted_col(cohort,unadj,adj,"dc_wa_nohosp")


################# no comorbidities #############################################
unadj = clogit(case ~ discharge + strata(episode_id), 
               cohort,
               cluster = moh_study_id,
               method  = "efron",
               robust = T
)
adj = clogit(case ~ discharge  +  homeless_hist + 
               sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m  +  cci_cat +
               oat_active_ad + num_active + 
               geq1_benzo + geq1_antipsyc + drug_toxicity +
               strata(episode_id), 
             data = cohort,
             cluster = moh_study_id,
             method  = "efron",
             robust = T)

res_nocom <- get_pasted_col(cohort,unadj,adj,"nocom")


#################### with drugs_opioid #########################################
cco_adj_oud <- clogit(case ~ discharge + homeless_hist + 
                        sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                        drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                        geq1_benzo + geq1_antipsyc + drug_toxicity + drugs_opioid +
                        strata(episode_id), 
                      data = cohort,
                      cluster = moh_study_id,
                      method  = "efron",
                      robust = T)
res_oud <- get_pasted_col(cohort,unadj,cco_adj_oud,"with_oud")


##################### sex interaction ######################################

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
    mutate(`exp, CI, p` = paste0(sprintf("%.2f",`exp(coef)`),", ",sprintf("%.2f",`lower .95`), "-",sprintf("%.2f",`upper .95`), ", ", "p", ifelse(`Pr(>|z|)` < 0.001, "<0.001", paste0("=", sprintf("%.3f",`Pr(>|z|)`))))) 
  if(se == T){
    model_summary <- model_summary %>% select(variable,`short description`,`exp, CI, p`, `robust se`, se = `se(coef)`)
  }else{
    model_summary <- model_summary %>% select(variable,`short description`,`exp, CI, p`, `robust se`)
  }
  return(model_summary)
}


sex_inter <- clogit(case ~ discharge*male_sex + homeless_hist + 
                      sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                      drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                      geq1_benzo + geq1_antipsyc + drug_toxicity +
                      strata(episode_id), 
                    data = cohort,
                    cluster = moh_study_id,
                    method  = "efron",
                    robust = T)


##################### one od per individual ######################################

#randomly select one od per individual
set.seed(1)
random_coh <- cohort %>% 
  filter(interval == "od") %>% 
  group_by(moh_study_id) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(episode_id) %>% 
  left_join(cohort)

unadj_r = clogit(case ~ discharge + strata(episode_id), 
                 random_coh)

adj_r = clogit(case ~ discharge+sdpr_pay+num_hosp_prev_6m+num_clinic_prev_6m+
                 oat_active_ad+num_active+geq1_benzo+
                 geq1_antipsyc+drug_toxicity+strata(episode_id), 
               data = random_coh)

res_random <- get_pasted_col(random_coh,unadj_r,adj_r,"random")

#select first eligible od
first_coh <- cohort %>% 
  filter(interval == "od") %>% 
  group_by(moh_study_id) %>% 
  arrange(date) %>% 
  slice_head(n=1) %>% 
  ungroup() %>% 
  select(episode_id) %>% 
  left_join(cohort)

unadj_f = clogit(case ~ discharge + strata(episode_id), 
                 data = first_coh)
adj_f = clogit(case ~ discharge+sdpr_pay+num_hosp_prev_6m+num_clinic_prev_6m+
                 oat_active_ad+num_active+geq1_benzo+
                 geq1_antipsyc+drug_toxicity+strata(episode_id), 
               data = first_coh)

res_firstod <- get_pasted_col(first_coh,unadj_f,adj_f,"first_od")


################## table ###################

#combine results
od_res <- readRDS("sensitivity_analyses_res_part1_wa.rds") %>% select(-estimate,-robust.se,-conf.low,-conf.high)
od_res <- rbind(od_res,
                #res_hosp,res_nontransfer,
                #res_dcbma,
                res_dcwa,
                res_hosp_c_6m,res_hosp_c_1y,res_random,res_firstod,res_nocom,res_oud)
od_res <- rbind(od_res,get_pasted_col(readRDS("R:/working/AMA-OD_cco/NH/results/Table D/od_date/cco_cohort_date.rds"),
                                      readRDS("R:/working/AMA-OD_cco/NH/results/Table D/od_date/cco_date.rds"),
                                      readRDS("R:/working/AMA-OD_cco/NH/results/Table D/od_date/cco_adj_date.rds"),
                                      "od_date"))
od_res



#format
names <- tibble(Label = c("exposure interval", "exp_2d", "exp_3d","exp_7d","exp_14d", "exp_21d","exp_28d","exp_56d","exp_91d","exp_182d","",
                          "control number", "control3", "control2", "control1_6m","control1_1y", "",
                          "outcome","fatal_od","nonfatal_od","cohort_pm6mo","cohort_p6_12mo","cohort_p12mo", "",
                          "exposure","dc_bma_nohosp","dc_wa_nohosp","",
                          "cohort","hosp_cohort_6m","hosp_cohort_1y","random","first_od","",
                          "covariate","nocom","with_oud","",
                          "date","od_date",""),
                lab_show = c("Alternate exposure interval length", "2 days", "3 days","7 days","14 days", "21 days","28 days (primary)","8 weeks","13 weeks","26 weeks","",
                             "Alternate period between overdose and control dates", "3 control dates = t0-6m, t0-1y, t0-1.5y", "2 control dates = t0-6m, t0-1y",
                             "1 control date = t0-6m","1 control date = t0-1y", "",
                             "Alternate outcome","Outcome = fatal overdose","Outcome = non-fatal overdose",
                             "Outcome = non-fatal overdose (control = t0-6mo and t0+6m)","Outcome = non-fatal overdose (control = t0+6mo and t0+12m)","Outcome = non-fatal overdose (control = t0+12m)","",
                             "Alternate exposure","Exposure = discharge BMA","Exposure = discharge WA","",
                             "Alternate cohort","Cohort = discharged within pre-overdose interval and pre-control(t0-6m) interval",
                             "Cohort = discharged within pre-overdose interval and pre-control(t0-1y) interval",
                             "Cohort = random od for each individual","Cohort = first eligible of for each individual","",
                             "Alternate adjustment","Not adjusted by co-morbidities","Adjusted by OUD as a comorbidity","",
                             "Alternate overdose date","Earliest overdose encounter",""))
bold <- c("exposure interval","control number","outcome","exposure","cohort","covariate","date")
names <- names %>%  mutate(fontface = ifelse(Label %in% bold, "bold","plain"))
res <- names %>% left_join(od_res)



#save 
wb <- loadWorkbook("R:/working/AMA-OD_cco/NH/results/Table D/AMA-OD_cco - Table 4 v4 WA.xlsx")
#sensitivity analyses
addWorksheet(wb, "Table 4 WA1")
writeData(wb,"Table 4 WA1", res %>% select(-fontface,-Label),rowNames = F)
fontface <- createStyle(fontName = "Calibri",fontSize = 12, textDecoration = "bold")
cells_to_style <- expand.grid(rows = which(res$fontface == "bold")+1, cols = 1:ncol(res))
addStyle(wb,"Table 4 WA1", fontface, rows = cells_to_style$rows,cols = cells_to_style$cols)

#save
saveWorkbook(wb, "R:/working/AMA-OD_cco/NH/results/Table D/AMA-OD_cco - Table 4 v4 WA.xlsx",overwrite = T) 

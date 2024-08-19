###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Daly-Grafstein D, Yu Y, Erdelyi S, Staples JA. AMA-OD_cco Conditional logistic regression & bootstrapped resamples for variable selection validation. 2023 Dec 22. Retrieved from: ***LINK*** 
# Conditional logistic regression 
# Use Bootstrap resamples to validate variable selection
# Author: Daniel (adapted for AMA-OD by Nicole)
##########################################################

source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")
setwd("R:/working/AMA-OD_cco/NH/results/Table C")
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")

# Bootstrap
cco_unadj <- clogit(case ~ discharge+strata(episode_id), cohort)
cco_adj_full <- clogit(case ~ discharge+homeless_hist+sdpr_pay+num_hosp_prev_6m+num_clinic_prev_6m+
                         alcohol+drugs_nonopioid+psych+cci_cat+oat_active_ad+num_active+geq1_benzo+geq1_opioid+
                         geq1_antipsyc+drug_toxicity+strata(episode_id), cohort)
cco_adj_full_r <- clogit(case ~ discharge+homeless_hist+sdpr_pay+num_hosp_prev_6m+num_clinic_prev_6m+
                         alcohol+drugs_nonopioid+psych+cci_cat+oat_active_ad+num_active+geq1_benzo+geq1_opioid+
                         geq1_antipsyc+drug_toxicity+strata(episode_id), cohort,
                         cluster = moh_study_id,
                         method  = "efron",
                         robust = T)
cco_adj <- clogit(case ~ discharge + homeless_hist + 
                    sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                    drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                    geq1_benzo + geq1_antipsyc + drug_toxicity + strata(episode_id), cohort)
cco_adj_r <- clogit(case ~ discharge + homeless_hist + 
                      sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                      drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                      geq1_benzo + geq1_antipsyc + drug_toxicity + strata(episode_id), 
                    data = cohort,
                    cluster = moh_study_id,
                    method  = "efron",
                    robust = T)
#summary

model_sum <- function(m){
  summary <- summary(m) 
  model_summary <- 
    cbind(summary$coefficients ,summary$conf.int[,c(3,4)])
  model_summary <- as.data.frame(round(model_summary,6))
  model_summary$term <- rownames(model_summary)
  model_summary <- model_summary %>% 
    select(term, est = coef, or = `exp(coef)`, se = `robust se`, `Pr(>|z|)`, `lower .95`, `upper .95`) 
  return(model_summary)
}

get_coef <- function(m){
  summary <- summary(m) 
  model_summary <- 
    cbind(summary$coefficients[,-1] ,summary$conf.int[,c(3,4)])
  model_summary <- as.data.frame(round(model_summary,6))
  model_summary$term <- rownames(model_summary)
  model_summary <- model_summary %>% 
    select(variable = term, `exp(coef)`, `robust se`, `Pr(>|z|)`, `lower .95`, `upper .95`) %>% 
    mutate(`exp, CI, p` = paste0(round(`exp(coef)`,2),", ",round(`lower .95`,2), "-",round(`upper .95`,2), ", ", "p ", ifelse(`Pr(>|z|)` < 0.001, "< 0.001", paste("= ", round(`Pr(>|z|)`,3))))) %>% 
    select(variable,`exp, CI, p`, `robust se`)
  return(model_summary)
}

# create dataframe to store results of bootstrap variable selection
boot_metrics <- tibble(predictors = model_sum(cco_adj_full_r) %>% pull(term),
                       global_est = model_sum(cco_adj_full_r) %>% pull(est),
                       global_se = model_sum(cco_adj_full_r) %>% pull(se))

boot_metrics <- boot_metrics %>% left_join(model_sum(cco_adj_r) %>% select(term, selected_est = est, selected_se = se),
                                           by = c("predictors" = "term"))
# specify number of bootstrap resamples and initiate matrices to store results
bootnum <- 1000
boot_est <- boot_se <- matrix(0, ncol = length(model_sum(cco_adj_full_r) %>% pull(term)), nrow = bootnum,
                              dimnames = list(NULL, model_sum(cco_adj_full_r) %>% pull(term)))

# repeatedly take bootstrap resamples
# resample the matched set instead of each individual row
# for each resample, perform variable selection via backward elimination
# save chosen variables and their estimates/SE, for variables not chosen set to 0
episode_list <- unique(cohort$episode_id)
cohort <- cohort %>% mutate(index = 1:nrow(cohort))
set.seed(1)
for(i in 1:bootnum){
  # running counter
  print(i)
  od_boot_epi <- sample(episode_list, replace= T) %>% as_tibble()
  od_boot_epi <- od_boot_epi %>% left_join(cohort %>% select(episode_id,index), by = c("value" = "episode_id"), relationship = "many-to-many")
  od_boot <- cohort[od_boot_epi$index,]
  boot_mod <- step(clogit(case ~ discharge+homeless_hist+sdpr_pay+num_hosp_prev_6m+num_clinic_prev_6m+
                         alcohol+drugs_nonopioid+psych+cci_cat+oat_active_ad+num_active+geq1_benzo+geq1_opioid+
                         geq1_antipsyc+drug_toxicity+strata(episode_id),
                          data = od_boot),
                   scope = list(lower = formula(cco_unadj),
                                upper = formula(cco_adj_full)),
                   direction = "backward", trace = 0)
  boot_mod <- clogit(as.formula(paste("case ~", sub(".*~ ","",boot_mod$formula)[3])),
                     data = od_boot,
                     cluster = moh_study_id,
                     method  = "efron",
                     robust = T)
  boot_est[i, names(coef(boot_mod))] <- model_sum(boot_mod)[,"est"]
  boot_se[i, names(coef(boot_mod))] <- model_sum(boot_mod)[, "se"]
}
saveRDS(boot_est,"boot_est.rds")
saveRDS(boot_se,"boot_se.rds")


# calculated bootstrap metrics as in Heinze et al. (2018) Table 5 - ref 14 of STA-MVC - Appendix v29 [Heinze G, Wallisch C, Dunkler D. Variable selection - A review and recommendations for the practicing statistician. Biom J. 2018 May;60(3):431-449.]
boot_metrics <- boot_metrics %>% 
  mutate(boot_inclusion = apply((boot_est != 0)*1, 2, function(x) sum(x)/length(x) * 100),
         rmsdratio = apply(((t(boot_est) - global_est)^2), 1, function(x) sqrt(mean(x)))/global_se,
         boot_meanratio = apply(boot_est, 2, mean)/global_est,
         boot_rel_bias = (boot_meanratio / (boot_inclusion / 100) - 1)*100,
         boot_median = apply(boot_est, 2, median),
         boot_median_OR = exp(boot_median),
         boot_025_quan = apply(boot_est, 2, function(x) quantile(x, 0.025)),
         boot_975_quan = apply(boot_est, 2, function(x) quantile(x, 0.975)),
         boot_025_quan_OR = exp(boot_025_quan),
         boot_975_quan_OR = exp(boot_975_quan))




# model frequency metrics as in Heinze et al. (2018) Table 6
boot_01 <- cbind(((boot_est != 0)*1)[, 
                                     boot_metrics$predictors[order(boot_metrics$boot_inclusion, decreasing = T)]], 
                 count = rep(1, times = bootnum))

model_freq = tibble(aggregate(count ~., data = boot_01, sum))
model_freq <- model_freq %>% 
  mutate(percent = count / bootnum *100) %>% 
  arrange(desc(percent)) %>% 
  mutate(cum_percent = cumsum(percent))

# create column denoting which variables not included in model
model_freq[,"vars_elim"] <- apply(model_freq[,c(2:(ncol(model_freq) -3))], 1, 
                                  function(x) paste(names(x[x==0]), collapse= " & "))

# select cols, create column denoting which is the final selected model
elim_vars <- names(coef(cco_adj_full_r))[!names(coef(cco_adj_full_r)) %in% names(coef(cco_adj_r))]
if(length(elim_vars) == 0){
  model_freq <- model_freq %>% 
    select(vars_elim, count, percent, cum_percent) %>% 
    mutate(selected_model = as.numeric(vars_elim == ""))
} else{
  model_freq <- model_freq %>% 
    select(vars_elim, count, percent, cum_percent) %>% 
    mutate(selected_model = as.numeric(setequal(sort(trimws(unlist(str_split(vars_elim, "&")))), sort(elim_vars))))
}

#final present table
final_res <- boot_metrics %>%
  select(predictors,boot_inclusion,boot_025_quan_OR,boot_975_quan_OR) %>% 
  left_join(get_coef(cco_adj_r), by = c("predictors" = "variable" )) %>% 
  select(predictors, `exp, CI, p`, `robust se`, boot_inclusion,boot_025_quan_OR,boot_975_quan_OR)
  

res <- list(
     "Selected Model" = get_coef(cco_adj_r), 
     "Global Model" = get_coef(cco_adj_full_r),
     "Bootstrap Variable Metrics" = boot_metrics,
     "Model Selection Metrics" = model_freq,
     "Final Table" = final_res)
##########################################################
# Save
write.xlsx(res, file = "AMA-OD_cco - coef_bootstrap v3.xlsx")









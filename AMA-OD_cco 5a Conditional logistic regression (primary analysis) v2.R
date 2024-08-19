###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Khan M, Daly-Grafstein D, Erdelyi S, Staples JA. AMA-OD_cco Conditional logistic regression. 2023 Dec 18. Retrieved from: ***LINK*** 
# Conditional logistic regression
# Backward selection
# OR result table
# Population attributable risk
# Author: Xiao Nicole Hu 
# Date: 2023-09-20
# Updated: 2023-12-18
##########################################################
get_pasted_col <- function(df,unadj_model, adj_model, label){
  res <- 
  cbind(
    label,
    df %>% filter(interval == "od") %>% nrow(),
    df %>% filter(interval == "control_1",dc_bma == 1) %>% nrow(),
    df %>% filter(interval == "control_2",dc_bma == 1) %>% nrow(),
    df %>% filter(interval2 == "control",dc_bma == 1) %>% nrow(),
    tidy(unadj_model,exponentiate = T, conf.int = T)[1,] %>% 
      select(estimate,conf.low,conf.high,p.value) %>% 
      mutate(unadjusted = paste0(round(estimate,2),", ",round(conf.low,2), "-",round(conf.high,2), ", ", "p = ", ifelse(p.value < 0.001, "< 0.001", round(p.value,3)))) %>% 
      select(unadjusted),
    tidy(adj_model,exponentiate = T, conf.int = T)[1,] %>% 
      select(estimate,conf.low,conf.high,p.value) %>% 
      mutate(adjusted = paste0(round(estimate,2),", ",round(conf.low,2), "-",round(conf.high,2), ", ", "p = ", ifelse(p.value < 0.001, "< 0.001", round(p.value,3)))) %>% 
      select(adjusted)
  )
  colnames(res) <- c("Label","OD number","BMA in control 1","BMA in control 2","BMA in pooled control","Unadjusted OR","Adjusted OR")
  return(res)
}
# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")
setwd("R:/working/AMA-OD_cco/NH/results/Table C")
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")


# Unadjusted effects
cco_unadj <- clogit(case ~ discharge+strata(episode_id), cohort)

# Adjusted effects
cco_adj_full <- clogit(case ~ discharge+homeless_hist+sdpr_pay+num_hosp_prev_6m+num_clinic_prev_6m+
                  alcohol+drugs_nonopioid+psych+cci_cat+oat_active_ad+num_active+geq1_benzo+geq1_opioid+
                  geq1_antipsyc+drug_toxicity+strata(episode_id), cohort)
vif(cco_adj_full)

#only geq1_opioid is removed from BW selection
step_AIC <-  step(cco_adj_full,
                  scope = list(lower = formula(cco_unadj),
                             upper = formula(cco_adj_full)),
                  direction = "backward") 

cco_adj <- clogit(case ~ discharge  + homeless_hist + 
                    sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                    drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                    geq1_benzo + geq1_antipsyc + drug_toxicity + strata(episode_id), cohort)



cco_adj_r <- clogit(case ~ discharge  + homeless_hist + 
                      sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                      drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                      geq1_benzo + geq1_antipsyc + drug_toxicity + strata(episode_id), 
                  data = cohort,
                  cluster = moh_study_id,
                  method  = "efron",
                  robust = T)

saveRDS(cco_adj_r,"cco_adj.rds")

############################ Get coef ##########################################
get_pasted_col(cohort,cco_unadj,cco_adj,"cco")

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


cco_adj_male <- clogit(case ~ discharge + homeless_hist + 
                      sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                      drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                      geq1_benzo + geq1_antipsyc + drug_toxicity + strata(episode_id), 
                    data = cohort %>% filter(gender == "M"),
                    cluster = moh_study_id,
                    method  = "efron",
                    robust = T)

cco_adj_female <- clogit(case ~ discharge  + homeless_hist + 
                         sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                         drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                         geq1_benzo + geq1_antipsyc + drug_toxicity + strata(episode_id), 
                       data = cohort %>% filter(gender == "F"),
                       cluster = moh_study_id,
                       method  = "efron",
                       robust = T)

###################### population attributable risk ############################
#multiplying proportion of overdoses with dc_bma in the pre-od interval by (OR - 1)/OR

#unadjusted OR
cco_unadj <- clogit(case ~ discharge+strata(episode_id), cohort)
or_bma <- (tidy(cco_unadj,exponentiate = T, conf.int = T) %>% select(-term,-std.error,-statistic,-p.value))[1,]
or_wa <- (tidy(cco_unadj,exponentiate = T, conf.int = T) %>% select(-term,-std.error,-statistic,-p.value))[2,]


#proportion of overdoses with dc_bma in the pre-od interval
n_bma <- cohort %>% filter(interval == "od", dc_bma == 1) %>% nrow()
total <- cohort %>% filter(interval == "od") %>% nrow()
prop_bma <- n_bma/total
 

#attributable risk
prop_bma*(or_bma-1)/or_bma*100


#proportion of overdoses with dc_wa in the pre-od interval
n_wa <- cohort %>% filter(interval == "od", dc_wa == 1) %>% nrow()
prop_wa <- n_wa/total
 

#attributable risk
prop_wa*(or_wa-1)/or_wa*100

#res
PAR <- rbind(cbind(or = paste0(sprintf("%.2f",or_bma$estimate),", ",sprintf("%.2f",or_bma$conf.low), "-",sprintf("%.2f",or_bma$conf.high)),
                   n = n_bma,total,prop_bma*(or_bma-1)/or_bma*100),
             cbind(or = paste0(sprintf("%.2f",or_wa$estimate),", ",sprintf("%.2f",or_wa$conf.low), "-",sprintf("%.2f",or_wa$conf.high)),
                   n = n_wa,total,prop_wa*(or_wa-1)/or_wa*100))
PAR$PAR <- c("bma","wa")
colnames(PAR) <- c("unadjusted OR","n","total","PAR estimate","conf.low","conf.high", "group")
PAR


res_coef <- list(
  "All cohort" = get_coef(cco_adj_r),
  "Male" = get_coef(cco_adj_male),
  "Female" = get_coef(cco_adj_female),
  "PAR" = PAR)


write.xlsx(res_coef, file = "AMA-OD_cco - coef v3.1.xlsx")




######################## model assumption ######################################

#residual
dev_res <- residuals(cco_adj,type = "deviance")
plot(fitted(cco_adj), dev_res)
abline(h = 0)

anova(cco_adj_full,cco_adj)#keep similar fitness as the full model



#check linearity between continuous variable and log odds
cco_adj <- clogit(case ~ discharge + homeless_hist + 
                    sdpr_pay + num_hosp_prev_6m + num_clinic_prev_6m + alcohol + 
                    drugs_nonopioid + psych + cci_cat + oat_active_ad + num_active + 
                    geq1_benzo + geq1_antipsyc + drug_toxicity + strata(episode_id), cohort)

cohort$predict <- predict(cco_adj,type = "lp")

ggplot(cohort, aes(x = sdpr_pay, y = predict)) + 
  geom_point()+
  geom_smooth()

ggplot(cohort, aes(x = num_hosp_prev_6m, y = predict)) + 
  geom_point()+
  geom_smooth()

ggplot(cohort, aes(x = num_clinic_prev_6m, y = predict)) + 
  geom_point()+
  geom_smooth()

ggplot(cohort, aes(x = num_active, y = predict)) + 
  geom_point()+
  geom_smooth()

ggplot(cohort, aes(x = drug_toxicity, y = predict)) + 
  geom_point()+
  geom_smooth()

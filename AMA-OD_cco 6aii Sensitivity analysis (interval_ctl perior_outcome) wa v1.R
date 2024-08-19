###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X,  Yu Y,Daly-Grafstein D, Khan M, Erdelyi S, Staples JA. AMA-OD_cco sensitivity analyses for physician-approved discharges. 2023 Sep 29. Retrieved from: ***LINK*** 
# Sensitivity analysis (1) [DC-WA]
# alternate exposure interval length
# alternate period
# alternate outcome
# Author: Nicole Hu 
# Date: 2023-09-29
# Updated: 2023-09-29
##########################################################
# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")
setwd("R:/working/AMA-OD_cco/NH/results/Table D")
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")
cco_adj <- readRDS("R:/working/AMA-OD_cco/NH/results/Table C/cco_adj.rds")

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
        mutate(unadjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p", ifelse(p.value < 0.001, "<0.001", paste0("=",sprintf("%.3f",p.value))))) %>% 
        select(unadjusted),
      tidy(adj_model,exponentiate = T, conf.int = T) %>% 
        filter(str_detect(term,"wa")) %>% 
        select(estimate,robust.se,conf.low,conf.high,p.value) %>% 
        mutate(adjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p", ifelse(p.value < 0.001, "<0.001", paste0("=",sprintf("%.3f",p.value))))) %>% 
        select(adjusted,estimate,robust.se,conf.low,conf.high)
    )
  colnames(res) <- c("Label","OD number","Exposure in pre-od","Exposure in control 1","Exposure in control 2","Exposure in pooled control","Unadjusted OR","Adjusted OR","estimate","robust.se","conf.low","conf.high")
  return(res)
}

#fit model in alternative datasets
sub_fit <- function(df,label){
  unadj = clogit(case ~ discharge + strata(episode_id), 
                 df,
                 cluster = moh_study_id,
                 method  = "efron",
                 robust = T
  )
  adj = clogit(as.formula(paste("case ~", sub(".*~ ","",cco_adj$formula)[3])), 
               data = df,
               cluster = moh_study_id,
               method  = "efron",
               robust = T
  )
  res = get_pasted_col(df,unadj,adj,label)
  return(res)
}

##################### alternate exposure interval length #######################
dad <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_dad.csv.gz")

#get discharge in pre-exposure interval with different length n
get_ndays_exp <- function(df,n){
  cohort <- df %>% 
    select(-discharge,-dc_bma,-dc_wa,-sep_disp) %>% 
    left_join(dad,
              by = "moh_study_id",
              relationship = "many-to-many") %>% 
    filter(sep_date >= (date-days(n)) & sep_date < date) %>% 
    group_by(episode_id,date) %>% 
    arrange(desc(sep_date),desc(ad_date)) %>% 
    slice_head(n = 1) %>% 
    ungroup() %>% 
    select(episode_id,date,sep_disp) %>%
    distinct() %>% 
    mutate(discharge = case_when(sep_disp %in%  c(6,12,61,62,63,64,65) ~ "bma",
                                 TRUE ~ "wa")) %>% 
    right_join(df %>% select(-discharge,-dc_bma,-dc_wa,-sep_disp)) %>% 
    mutate(discharge = replace_na(discharge,"no discharge")) %>% 
    mutate(dc_bma = as.numeric(discharge == "bma"),
           dc_wa = as.numeric(discharge == "wa")) %>% 
    arrange(episode_id, date)
  cohort$discharge <- relevel(as.factor(cohort$discharge), ref = "no discharge")
  return(cohort)
}




########################## alternate period ####################################
#3 controls
cohort_3c <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_3control.rds")
#1 control (6 months)
cohort_1c_6m <- cohort %>% filter(interval != "control_2")
#1 control (1 year)
cohort_1c_1y <- cohort %>% filter(interval != "control_1")

########################### alternate outcome ##################################
#fatal od
od_date <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_od_date.csv.gz")
fatal_cohort <- cohort %>% 
  left_join(od_date %>% select(moh_study_id,od_date,od_fatal), 
            by = c("moh_study_id", "od_date")) %>% 
  filter(od_fatal == 1) 

#non-fatal od
nonfatal_cohort <- cohort %>% 
  left_join(od_date %>% select(moh_study_id,od_date,od_fatal), 
            by = c("moh_study_id", "od_date")) %>% 
  filter(od_fatal == 0) 

########################### non-fatal and lookback period ##################################
cohort_lb <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_lbcontrol.rds")
eligible_od <- cohort %>% filter(interval == "od", date %m+% days(364) <= as.Date("2019-12-31")) %>% select(episode_id)
#t0 +- 6 month control and nonfatal
cohort_pm6mo <- eligible_od %>% 
  left_join(cohort_lb) %>% 
  filter(interval %in% c("od","control_lb6m","control_1")) %>% 
  left_join(od_date %>% select(moh_study_id,od_date,od_fatal), 
            by = c("moh_study_id", "od_date")) %>% 
  filter(od_fatal == 0)

#t0 +6 mo and + 12mo control and nonfatal
cohort_p6_12mo <- eligible_od %>% 
  left_join(cohort_lb) %>% 
  filter(interval %in% c("od","control_lb6m","control_lb12m")) %>% 
  left_join(od_date %>% select(moh_study_id,od_date,od_fatal), 
            by = c("moh_study_id", "od_date")) %>% 
  filter(od_fatal == 0) 

#t0 +12mo
cohort_p6_12mo <- eligible_od %>% 
  left_join(cohort_lb) %>% 
  filter(interval %in% c("od","control_lb12m")) %>% 
  left_join(od_date %>% select(moh_study_id,od_date,od_fatal), 
            by = c("moh_study_id", "od_date")) %>% 
  filter(od_fatal == 0) 


############## result #####################################################
data <- list(exp_2d = get_ndays_exp(cohort,2),
             exp_3d = get_ndays_exp(cohort,3),
             exp_7d = get_ndays_exp(cohort,7),
             exp_14d = get_ndays_exp(cohort,14),
             exp_21d = get_ndays_exp(cohort,21),
             exp_28d = get_ndays_exp(cohort,28),
             exp_56d = get_ndays_exp(cohort,56),
             exp_91d = get_ndays_exp(cohort,91),
             exp_182d = get_ndays_exp(cohort,182),
             control3 = cohort_3c,
             control2 = cohort,
             control1_6m = cohort_1c_6m,
             control1_1y = cohort_1c_1y,
             fatal_od = fatal_cohort,
             nonfatal_od = nonfatal_cohort,
             cohort_pm6mo = cohort_pm6mo,
             cohort_p6_12mo = cohort_p6_12mo,
             cohort_pod = readRDS("R:/working/AMA-OD_cco/NH/results/Table D/allow prior od/cco_coh_all_pod.rds")
)

labels <- names(data)

od_res <- numeric(0)
for (l in labels){
  od_res <- rbind(od_res,sub_fit(data[[l]],l))
}
od_res
saveRDS(od_res,"sensitivity_analyses_res_part1_wa.rds")

#################### plot with different interval length ##################

od_res <- readRDS("sensitivity_analyses_res_part1_wa.rds") %>% select(Label,estimate,robust.se,conf.low,conf.high) %>% head(9)
#ests$estimate <- as.numeric(ests$estimate)
od_res$label <- c(2,3,7,14,21,28,56,91,182)

plot.ts <- ggplot(od_res,aes(x = label, y = estimate)) +
  geom_point(aes(size = (1/robust.se)^2), shape = "square") +
  #geom_smooth(formula = y~x)+
  geom_segment(aes(x = label, xend = label, y = conf.low, yend = conf.high)) +
  geom_hline(yintercept = 1, linetype = "dotted", size = 0.25)+ 
  scale_y_continuous(trans = "log",
                     breaks = c(0, 0.7, 1, 1.5, 2, 2.5),
                     labels = c(0, 0.7, 1, 1.5, 2, 2.5)) +
  scale_x_continuous(#trans = "log10",
    breaks = c(2,7,14,21,28,56,91,182),
    labels = c(2,7,14,21,28,56,91,182)) +
  scale_size_continuous(guide = "none") +
  xlab("Lookback interval length (days)") +
  ylab("Adjusted odds ratio") +
  theme_stapleslab
plot.ts
ggsave("R:/working/AMA-OD_cco/NH/results/fig/cco_interval_length_forest_plot_wa v1.tiff", 
       width=14, height=8, units="in", dpi=600, compression = "lzw",
       plot = plot.ts)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/cco_interval_length_forest_plot_wa v1.pdf", 
       width=14, height=8, units="in", dpi=600, 
       plot = plot.ts)



plot.ts.sub <- ggplot(od_res[c(-1,-2),],aes(x = label, y = estimate)) +
  geom_point(aes(size = (1/robust.se)^2), shape = "square") +
  #geom_smooth(formula = y~x)+
  geom_segment(aes(x = label, xend = label, y = conf.low, yend = conf.high)) +
  geom_hline(yintercept = 1, linetype = "dotted", size = 0.25)+ 
  scale_y_continuous(trans = "log",
                     breaks = c(0, 0.7, 1, 1.5, 2, 2.5),
                     labels = c(0, 0.7, 1, 1.5, 2, 2.5)) +
  scale_x_continuous(#trans = "log10",
    breaks = c(7,14,21,28,56,91,182),
    labels = c(7,14,21,28,56,91,182)) +
  scale_size_continuous(guide = "none") +
  xlab("Lookback interval length (days)") +
  ylab("Adjusted odds ratio") +
  theme_stapleslab
plot.ts.sub

ggsave("R:/working/AMA-OD_cco/NH/results/fig/cco_interval_length_forest_plot_remove2&3_wa v1.tiff", 
       width=14, height=8, units="in", dpi=600, compression = "lzw",
       plot = plot.ts.sub)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/cco_interval_length_forest_plot_remove2&3_wa v1.pdf", 
       width=14, height=8, units="in", dpi=600, 
       plot = plot.ts.sub)


###################Combined plot #############
od_res_wa <- readRDS("sensitivity_analyses_res_part1_wa.rds") %>% select(Label,estimate,robust.se,conf.low,conf.high) %>% head(9) %>% mutate(dc = "wa")
od_res_bma <- readRDS("sensitivity_analyses_res_part1.rds") %>% select(Label,estimate,robust.se,conf.low,conf.high) %>% head(9) %>% mutate(dc = "bma")
od_res_com <- rbind(od_res_wa,od_res_bma)

#ests$estimate <- as.numeric(ests$estimate)
od_res_com$label <- rep(c(2,3,7,14,21,28,56,91,182),2)
dw  = 0.3
od_res_com <- od_res_com %>% 
  mutate(xdodge = ifelse(dc == "bma", as.numeric(label) - dw,  as.numeric(label) + dw))
plot.ts <- ggplot(od_res_com,aes(x = xdodge, y = estimate, color = dc)) +
  geom_point(aes(size = (1/robust.se)^2), shape = "square") +
  #geom_smooth(formula = y~x)+
  geom_segment(aes(x = xdodge, xend = xdodge, y = conf.low, yend = conf.high), position = position_dodge(width = dw)) +
  geom_hline(yintercept = 1, linetype = "dotted", size = 0.25)+ 
  scale_y_continuous(trans = "log",
                     breaks = c(0, 0.7, 1, 1.5, 2, 2.5),
                     labels = c(0, 0.7, 1, 1.5, 2, 2.5)) +
  scale_x_continuous(#trans = "log10",
    breaks = c(2,7,14,21,28,56,91,182),
    labels = c(2,7,14,21,28,56,91,182)) +
  scale_size_continuous(guide = "none") +
  xlab("Lookback interval length (days)") +
  ylab("Adjusted odds ratio") +
  scale_color_manual(values = c("wa" ="dodgerblue4" ,"bma" ="red4" ))+
  theme_stapleslab+
  theme(legend.position = "none")
plot.ts
ggsave("R:/working/AMA-OD_cco/NH/results/fig/cco_interval_length_forest_plot_all v2.tiff", 
       width=14, height=8, units="in", dpi=600, compression = "lzw",
       plot = plot.ts)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/cco_interval_length_forest_plot_all v2.pdf", 
       width=14, height=8, units="in", dpi=600, 
       plot = plot.ts)



plot.ts.sub <- ggplot(od_res_com %>% filter(label > 3),aes(x = xdodge, y = estimate,color = dc)) +
  geom_point(aes(size = (1/robust.se)^2), shape = "square") +
  #geom_smooth(formula = y~x)+
  geom_segment(aes(x = xdodge, xend = xdodge, y = conf.low, yend = conf.high)) +
  geom_hline(yintercept = 1, linetype = "dotted", size = 0.25)+ 
  scale_y_continuous(trans = "log",
                     breaks = c(0, 0.7, 1, 1.5, 2, 2.5),
                     labels = c(0, 0.7, 1, 1.5, 2, 2.5)) +
  scale_x_continuous(#trans = "log10",
    breaks = c(7,14,21,28,56,91,182),
    labels = c(7,14,21,28,56,91,182)) +
  scale_size_continuous(guide = "none") +
  xlab("Lookback interval length (days)") +
  ylab("Adjusted odds ratio") +
  scale_color_manual(values = c("wa" ="dodgerblue4" ,"bma" ="red4" ))+
  theme_stapleslab+
  theme(legend.position = "none")
plot.ts.sub

ggsave("R:/working/AMA-OD_cco/NH/results/fig/cco_interval_length_forest_plot_remove2&3_all v2.tiff", 
       width=14, height=8, units="in", dpi=600, compression = "lzw",
       plot = plot.ts.sub)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/cco_interval_length_forest_plot_remove2&3_all v2.pdf", 
       width=14, height=8, units="in", dpi=600, 
       plot = plot.ts.sub)


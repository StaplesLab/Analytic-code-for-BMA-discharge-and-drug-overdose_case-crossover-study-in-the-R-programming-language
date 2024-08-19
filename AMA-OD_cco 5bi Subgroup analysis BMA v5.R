###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Khan M, Daly-Grafstein D, Erdelyi S, Staples JA. AMA-OD_cco subgroup analysis for BMA discharges. 2023 Sep 29. Retrieved from: ***LINK*** 
# Subgroup analysis
# Author: Xiao Nicole Hu 
# Date: 2023-09-25
# Updated: 2023-09-29
##########################################################
get_pasted_col <- function(df,unadj_model, adj_model, label){
  res <- 
    cbind(
      label,
      paste0(df %>% filter(interval == "od") %>% nrow(), " (", round((length(which(df$interval == "od"))*100)/27584,1) ,")"),
      df %>% filter(interval == "od",dc_bma == 1) %>% nrow(),
      df %>% filter(interval == "control_1",dc_bma == 1) %>% nrow(),
      df %>% filter(interval == "control_2",dc_bma == 1) %>% nrow(),
      df %>% filter(interval2 == "control",dc_bma == 1) %>% nrow(),
      tidy(unadj_model,exponentiate = T, conf.int = T)[1,] %>% 
        select(estimate,conf.low,conf.high,p.value) %>% 
        mutate(unadjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p", ifelse(p.value < 0.001, "<0.001", paste0("=",sprintf("%.3f",p.value))))) %>% 
        select(unadjusted),
      tidy(adj_model,exponentiate = T, conf.int = T)[1,] %>% 
        select(estimate,robust.se,conf.low,conf.high,p.value) %>% 
        mutate(adjusted = paste0(sprintf("%.2f",estimate),", ",sprintf("%.2f",conf.low), "-",sprintf("%.2f",conf.high), ", ", "p", ifelse(p.value < 0.001, "<0.001", paste0("=",sprintf("%.3f",p.value))))) %>% 
        select(adjusted,estimate,robust.se,conf.low,conf.high)
    )
  colnames(res) <- c("Label","OD number","BMA in pre-od","BMA in control 1","BMA in control 2","BMA in pooled control","Unadjusted OR","Adjusted OR","estimate","robust.se","conf.low","conf.high")
  return(res)
}

sub_fit <- function(df,label,v,plot = F){
  unadj = clogit(case ~ discharge + strata(episode_id), 
                 df,
                 cluster = moh_study_id,
                 method  = "efron",
                 robust = T
  )
  adj = clogit(as.formula(gsub(paste(v,"+"),"", paste("case ~", sub(".*~ ","",cco_adj$formula)[3]))), 
               data = df,
               cluster = moh_study_id,
               method  = "efron",
               robust = T
  )
  res = get_pasted_col(df,unadj,adj,label)
  return(res)
}


# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")
setwd("R:/working/AMA-OD_cco/NH/results/Table C")
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")
cco_adj <- readRDS("cco_adj.rds")



################################# subgroups ###################################
data <- list(primary = cohort,
             male = cohort %>% filter(gender == "M"),
             female = cohort %>% filter(gender == "F"),
             age18 = cohort %>% filter(interval == "od", age_group == "18-29") %>% select(episode_id) %>% left_join(cohort),
             age30 = cohort %>% filter(interval == "od", age_group == "30-49") %>% select(episode_id) %>% left_join(cohort),
             age50 = cohort %>% filter(interval == "od", age_group == "50+") %>% select(episode_id) %>% left_join(cohort),
             rural = cohort %>% filter(interval == "od", pop_density_class == "Rural") %>% select(episode_id) %>% left_join(cohort),
             urban = cohort %>% filter(interval == "od", pop_density_class == "Urban") %>% select(episode_id) %>% left_join(cohort),
             missing_res = cohort %>% filter(interval == "od", pop_density_class == "Missing") %>% select(episode_id) %>% left_join(cohort),
             iv_drug1 = cohort %>% filter(interval == "od", ivdu_hist == 1) %>% select(episode_id) %>% left_join(cohort),
             iv_drug0 = cohort %>% filter(interval == "od", ivdu_hist == 0) %>% select(episode_id) %>% left_join(cohort),
             act_num0 = cohort %>% filter(interval == "od", num_active_cat == "0 or 1") %>% select(episode_id) %>% left_join(cohort),
             act_num2 = cohort %>% filter(interval == "od", num_active_cat == "2+") %>% select(episode_id) %>% left_join(cohort),
             act_oat1 = cohort %>% filter(interval == "od", oat_active_ad == 1) %>% select(episode_id) %>% left_join(cohort),
             act_oat0 = cohort %>% filter(interval == "od", oat_active_ad == 0) %>% select(episode_id) %>% left_join(cohort),
             med_oat1 = cohort %>% filter(interval == "od", geq1_oat == 1) %>% select(episode_id) %>% left_join(cohort),
             med_oat0 = cohort %>% filter(interval == "od", geq1_oat == 0) %>% select(episode_id) %>% left_join(cohort),
             med_opioid1 = cohort %>% filter(interval == "od", geq1_opioid == 1) %>% select(episode_id) %>% left_join(cohort),
             med_opioid0 = cohort %>% filter(interval == "od", geq1_opioid == 0) %>% select(episode_id) %>% left_join(cohort),
             any_sbstc1 = cohort %>% filter(interval == "od", any_sbstcs == 1) %>% select(episode_id) %>% left_join(cohort),
             any_sbstc0 = cohort %>% filter(interval == "od", any_sbstcs == 0) %>% select(episode_id) %>% left_join(cohort),
             opioid1 = cohort %>% filter(interval == "od", drugs_opioid == 1) %>% select(episode_id) %>% left_join(cohort),
             opioid0 = cohort %>% filter(interval == "od", drugs_opioid == 0) %>% select(episode_id) %>% left_join(cohort),
             alcohol1 = cohort %>% filter(interval == "od", alcohol == 1) %>% select(episode_id) %>% left_join(cohort),
             alcohol0 = cohort %>% filter(interval == "od", alcohol == 0) %>% select(episode_id) %>% left_join(cohort),
             nonopioid1 = cohort %>% filter(interval == "od", drugs_nonopioid == 1) %>% select(episode_id) %>% left_join(cohort),
             nonopioid0 = cohort %>% filter(interval == "od", drugs_nonopioid == 0) %>% select(episode_id) %>% left_join(cohort),
             depression1 = cohort %>% filter(interval == "od", mood == 1) %>% select(episode_id) %>% left_join(cohort),
             depression0 = cohort %>% filter(interval == "od", mood == 0) %>% select(episode_id) %>% left_join(cohort),
             anxiety1 = cohort %>% filter(interval == "od", anxiety == 1) %>% select(episode_id) %>% left_join(cohort),
             anxiety0 = cohort %>% filter(interval == "od", anxiety == 0) %>% select(episode_id) %>% left_join(cohort),
             psych1 = cohort %>% filter(interval == "od", psych == 1) %>% select(episode_id) %>% left_join(cohort),
             psych0 = cohort %>% filter(interval == "od", psych == 0) %>% select(episode_id) %>% left_join(cohort),
             homeless1 = cohort %>% filter(interval == "od", homeless_hist == 1) %>% select(episode_id) %>% left_join(cohort),
             homeless0 = cohort %>% filter(interval == "od", homeless_hist == 0) %>% select(episode_id) %>% left_join(cohort)
             )

formula <- list(primary = "", male = "", female = "",age18 = "",age30 = "", age50 = "", rural = "",urban = "",missing_res = "", iv_drug1 = "", iv_drug0 = "", 
                act_oat1 = "oat_active_ad", act_oat0 = "oat_active_ad", med_opioid1 = "",med_opioid0 = "",any_sbstc1 = "", any_sbstc0 = "alcohol \\+ drugs_nonopioid",
                opioid1 = "",opioid0 = "",
                alcohol1 = "alcohol",alcohol0 = "alcohol", nonopioid1 = "drugs_nonopioid", nonopioid0 ="drugs_nonopioid", depression1 = "psych", depression0 = "",
                anxiety1 = "psych", anxiety0 = "",psych1 = "psych" ,psych0 = "psych",homeless1 = "homeless_hist",homeless0 = "homeless_hist")

labels <- names(data)

od_res <- numeric(0)
for (l in labels){
  od_res <- rbind(od_res,sub_fit(data[[l]],l,formula[[l]]))
}

od_res
saveRDS(od_res,"sub_od_res.rds")


###################### forest plot ###########################################
plot_res <- readRDS("sub_od_res.rds") %>% select(Label,estimate,robust.se,conf.low,conf.high)

lable_name <- read.xlsx("label_name.xlsx")
lable_name <- lable_name %>% replace(is.na(.), "")
plot_res <-  lable_name %>% 
  left_join(plot_res,by = c("label" = "Label") ) %>% 
  mutate(lab_break = paste0(1:nrow(lable_name), label)) %>% 
  mutate(lab_break = factor(lab_break, levels = rev(lab_break))) %>% 
  filter(include_in_plot == 1)

#saveRDS(plot_res,"plot_res.rds")
plot_res <- readRDS("plot_res.rds")
plot_subgroup <- ggplot(plot_res) +
  geom_point(aes(x = lab_break, y = estimate, size = 1/robust.se,color = subgroup_name), shape = 15) +
  geom_segment(aes(x = lab_break, xend = lab_break, y = conf.low, yend = conf.high,color = subgroup_name)) +
  geom_hline(yintercept = 1, lty = 2) +
  coord_flip() +
  scale_x_discrete(breaks = plot_res$lab_break,
                   labels = plot_res$lab_show) +
  scale_y_continuous(trans = "log10",
                     limits = c(0.5,10),
                     breaks = c(0.5,1,2,5,10),
                     labels = c(0.5,1,2,5,10)) +
  scale_colour_viridis_d(direction = -1,begin = 0, end = 0.8)+
  xlab(NULL) +
  ylab("Adjusted OR (95% CI)") +
  theme_stapleslab+
  #theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(margin = margin(0,0,0,0),
                                   hjust = 0,
                                   face = plot_res$fontface,
                                   size = 13))+
  guides(size = "none",color = "none")

plot_subgroup

ggsave("R:/working/AMA-OD_cco/NH/results/fig/cco_forest v3.tiff", 
       width=10, height=nrow(plot_res)*0.2 + 0.25, units="in", dpi=600, compression = "lzw", 
       plot = plot_subgroup)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/cco_forest v3.pdf", 
       width=10, height=nrow(plot_res)*0.2 + 0.25, units="in", dpi=600, 
       plot = plot_subgroup)


########################### table ############################################

table_res <- read.xlsx("label_name.xlsx") %>% select(label,lab_show,fontface) %>% 
  left_join(readRDS("sub_od_res.rds") %>% select(-estimate,-robust.se,-conf.low,-conf.high), by = c("label" = "Label")) %>% 
  filter(!is.na(label)) %>% 
  select(-label)
  
 
wb <- loadWorkbook("R:/working/AMA-OD_cco/NH/results/Table C/AMA-OD_cco - Table 3 v4.xlsx")
addWorksheet(wb, "Table 3 v4.1")
writeData(wb,"Table 3 v4.1", table_res %>% select(-fontface),rowNames = F)
fontface <- createStyle(fontName = "Calibri",fontSize = 12, textDecoration = "bold")
cells_to_style <- expand.grid(rows = which(table_res$fontface == "bold")+1, cols = 1:ncol(table_res))
addStyle(wb,"Table 3 v4.1", fontface, rows = cells_to_style$rows,cols = cells_to_style$cols)
saveWorkbook(wb, "R:/working/AMA-OD_cco/NH/results/Table C/AMA-OD_cco - Table 3 v4.1.xlsx",overwrite = T) 


###################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Hu X, Yu Y, Daly-Grafstein D, Khan M, Erdelyi S, Staples JA. AMA-OD_cco descriptive plots. 2023 Dec 21. Retrieved from: ***LINK*** 
# Descriptive plots
# number of hospitalizations in X months before overdose
# number of DC-BMAs in X months before overdose (range: 0-18 months)
# number of hospitalizations in X days before overdose (range: 2-28 days)
# number of DC-BMAs in X days before overdose (range: 2-28 days)
# monthly overdose risk berfore and after discharge
# monthly overdose risk berfore and after dc-bma
###################################


# libraries
source("R:/working/AMA-OD_coh/NH/code/AMA-OD - 0 Packages and working directory.R")
setwd("R:/working/AMA-OD_cco/NH/results/Table C")
cohort <- readRDS("R:/working/AMA-OD_cco/NH/results/Table A/cco_coh_all.rds")
dad <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_dad.csv.gz")

#get all hospitalizations in the past 82 weeks
cohort_hosp <- cohort %>% select(-sep_disp) %>% 
  filter(interval == "od") %>% 
  left_join(dad,
            by = "moh_study_id",
            relationship = "many-to-many") %>% 
  filter(sep_date >= (date-days(574)) & sep_date < date) %>% 
  select(episode_id,moh_study_id,od_date,sep_date,sep_disp) %>% 
  mutate(dc_bma = as.numeric(sep_disp %in% c(6,12,61,62,63,64,65))) %>% 
  mutate(diff_days = as.numeric(sep_date - od_date)) %>% 
  mutate(diff_weeks = ceiling(diff_days/7)) %>% 
  mutate(diff_months = ceiling(diff_days/30))

################# By months ####################################################
#get count of hospitalization in each month
cnt <- cohort_hosp %>% 
  count(diff_months) %>% 
  left_join(cohort_hosp %>% filter(dc_bma == 1) %>% count(diff_months), by = "diff_months") %>% 
  rename("any_dc" = n.x, "dc_bma" = n.y) %>% 
  mutate("non dc_bma" = any_dc - dc_bma) %>% 
  mutate(pct = round(dc_bma/any_dc*100,2)) %>% 
  mutate(interval = case_when(diff_months == 0 ~ "pre-od",
                              diff_months > -6  ~ "control",
                              diff_months == -6  ~ "pre_control",
                              diff_months > -12 ~ "control",
                              diff_months == -12 ~ "pre_control",
                              diff_months > -18 ~ "control",
                              TRUE ~ "pre_control"))


#make a long format data for bar plot
cnt_long <- reshape2::melt(cnt, id.vars = c("diff_months","pct","interval"), measure.vars = c("non dc_bma","dc_bma")) %>% 
  rename(Discharge = variable)

#DC-BMA
month_bar_bma <- cnt_long %>%  
  filter(diff_months != -19, Discharge == "dc_bma") %>% #sample size is too small
  ggplot(aes(x = diff_months, y = value, fill = interval)) +
  geom_bar(stat = "identity")+
  scale_x_continuous(breaks = c(0,-6,-12,-18),
                     labels = c("overdose",6,12,18)) +
  scale_fill_manual(values = c("dodgerblue4","red4","yellow4"))+
  labs(x = "Month intervals prior to overdose",
       y = "Number of BMA discharges") +
  theme_stapleslab+
  theme(legend.position = "none")
month_bar_bma
ggsave("R:/working/AMA-OD_cco/NH/results/fig/bma_before_od v2.tiff", 
       width=10, height=10, units="in", dpi=600, compression = "lzw", 
       plot = month_bar_bma)

ggsave("R:/working/AMA-OD_cco/NH/results/fig/bma_before_od v2.pdf", 
       width=10, height=10, units="in", dpi=600, 
       plot = month_bar_bma)

#Any discharge
month_bar_anydc <- cnt %>%  
  filter(diff_months != -19) %>% 
  ggplot(aes(x = diff_months, y = any_dc,fill = interval)) +
  geom_bar(stat = "identity")+
  scale_x_continuous(breaks = c(0,-6,-12,-18),
                     labels = c("overdose",6,12,18)) +
  scale_fill_manual(values = c("dodgerblue4","red4","yellow4"))+
  labs(x = "Month intervals prior to overdose",
       y = "Number of any discharges") +
  theme_stapleslab+
  theme(legend.position = "none")
month_bar_anydc

ggsave("R:/working/AMA-OD_cco/NH/results/fig/anydc_before_od v2.tiff", 
       width=10, height=10, units="in", dpi=600, compression = "lzw", 
       plot = month_bar_anydc)

ggsave("R:/working/AMA-OD_cco/NH/results/fig/anydc_before_od v2.pdf", 
       width=10, height=10, units="in", dpi=600, 
       plot = month_bar_anydc)

#BMA proportion
# cnt_long %>%  
#   filter(diff_months != -19,Discharge == "dc_bma") %>% 
#   ggplot(aes(x = diff_months, y = pct)) +
#   geom_line(aes(y = pct), color = "black", size = 1) +
#   geom_point(aes(y = pct), color = "black",size = 1)+
#   labs(x = "Months before index overdose",
#        y = "Discharge BMA percentage") +
#   theme_stapleslab+
#   theme(legend.position = "none")

#combined
month_combined_bma <- cnt_long %>%  
  filter(diff_months != -19,Discharge == "dc_bma") %>% 
  ggplot(aes(x = diff_months, y = value, fill = interval)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = c(0,-6,-12,-18),
                     labels = c("overdose",6,12,18)) +
  scale_fill_manual(values = c("dodgerblue4","red4","yellow4"))+
  geom_line(aes(y = pct*max(cnt_long$value)/100, group = 1), color = "black", size = 1) +
  geom_point(aes(y = pct*max(cnt_long$value)/100), color = "black",size = 1)+
  scale_y_continuous(
    name = "Number of BMA discharges",
    sec.axis = sec_axis(~.*100/max(cnt_long$value),name = "Discharge BMA percentage")
  )+
  labs(x = "Month intervals prior to overdose") +
  theme_stapleslab+
  theme(legend.position = "none") 
month_combined_bma
ggsave("R:/working/AMA-OD_cco/NH/results/fig/bma_before_od_pct_and_n v2.tiff", 
       width=10, height=10, units="in", dpi=600, compression = "lzw", 
       plot = month_combined_bma)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/bma_before_od_pct_and_n v2.pdf", 
       width=10, height=10, units="in", dpi=600, 
       plot = month_combined_bma)

################# By days ######################################################
cnt_d <- cohort_hosp %>% 
  filter(diff_days >= -28) %>% 
  count(diff_days) %>% 
  left_join(cohort_hosp %>% filter(dc_bma == 1) %>% count(diff_days), by = "diff_days") %>% 
  rename("any_dc" = n.x, "dc_bma" = n.y) %>% 
  mutate("non dc_bma" = any_dc - dc_bma) %>% 
  mutate(pct = round(dc_bma/any_dc*100,2))


cnt_long_d <- reshape2::melt(cnt_d, id.vars = c("diff_days","pct"), measure.vars = c("non dc_bma","dc_bma")) %>% 
  rename(Discharge = variable)
cnt_long_d <- cnt_long_d %>%  filter(diff_days != -1)

#PLOT

#DC-BMA
day_bar_bma <- cnt_long_d %>%  
  filter(Discharge == "dc_bma") %>% 
  ggplot(aes(x = diff_days, y = value)) +
  geom_bar(stat = "identity",fill = "dodgerblue4")+
  scale_x_continuous(breaks = c(-2,-7,-14,-21,-28),
                     labels = c(2,7,14,21,28))+
  labs(x = "Days prior to overdose",
       y = "Number of BMA discharges") +
  theme_stapleslab+
  theme(legend.position = "none")
day_bar_bma
ggsave("R:/working/AMA-OD_cco/NH/results/fig/bma_before_od_28d v2.tiff", 
       width=10, height=10, units="in", dpi=600, compression = "lzw", 
       plot = day_bar_bma)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/bma_before_od_28d v2.pdf", 
       width=10, height=10, units="in", dpi=600, 
       plot = day_bar_bma)

#Any discharge
day_bar_anydc <- cnt_d %>%  
  filter(diff_days != -1) %>% 
  ggplot(aes(x = diff_days, y = any_dc)) +
  geom_bar(stat = "identity",fill = "dodgerblue4")+
  scale_x_continuous(breaks = c(-2,-7,-14,-21,-28),
                     labels = c(2,7,14,21,28))+
  labs(x = "Days prior to overdose",
       y = "Number of BMA discharges") +
  theme_stapleslab+
  theme(legend.position = "none")
day_bar_anydc
ggsave("R:/working/AMA-OD_cco/NH/results/fig/anydc_before_od_28d v2.tiff", 
       width=10, height=10, units="in", dpi=600, compression = "lzw", 
       plot = day_bar_anydc)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/anydc_before_od_28d v2.pdf", 
       width=10, height=10, units="in", dpi=600, 
       plot = day_bar_anydc)


#combined
day_combined_bma <- cnt_long_d %>%  
  filter(Discharge == "dc_bma") %>% 
  ggplot(aes(x = diff_days, y = value)) +
  geom_bar(stat = "identity",fill = "dodgerblue4") +
  geom_line(aes(y = pct*max(cnt_long_d$value)/100, group = 1), color = "black", size = 1) +
  geom_point(aes(y = pct*max(cnt_long_d$value)/100), color = "black",size = 1)+
  scale_x_continuous(breaks = c(-2,-7,-14,-21,-28),
                     labels = c(2,7,14,21,28))+
  scale_y_continuous(
    name = "Number of BMA discharges",
    sec.axis = sec_axis(~.*100/max(cnt_long_d$value),name = "Discharge BMA percentage")
  )+
  labs(x = "Days prior to overdose") +
  theme_stapleslab+
  theme(legend.position = "none",
        axis.title.y.right = element_text(vjust= +1.75)) 
day_combined_bma
ggsave("R:/working/AMA-OD_cco/NH/results/fig/bma_before_od_pct_and_n_28d v2.tiff", 
       width=10, height=10, units="in", dpi=600, compression = "lzw", 
       plot = day_combined_bma)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/bma_before_od_pct_and_n_28d v2.pdf", 
       width=10, height=10, units="in", dpi=600, 
       plot = day_combined_bma)


######### monthly overdose risk berfore and after dc-bma ######################
#y axis: how many dc-bma has a overdose in the index month, in the first month after dc-bma, in the second month after dc-bma
cohort = readRDS("R:/working/AMA-OD_coh/NH/results/Table B/cohort_cox.rds")
od_date <- read_delim("R:/DATA/2023-05-09_190103/POC2023-04-05/vw_od_date.csv.gz")

#overdose for each individual who discharged between 2016 and 2018
month_od <- cohort %>% 
  select(-od_date) %>% 
  filter(year(sep_date) > 2015 & year(sep_date) < 2019) %>% #115091 rows
  left_join(od_date %>% select(moh_study_id,od_date)) %>% 
  select(dad_id,moh_study_id,sep_date,ad_date, od_date,dc) %>% 
  mutate(day_diff  = as.numeric(od_date - sep_date),
         ad_day_diff = as.numeric(od_date - ad_date),
         month_diff = case_when((ad_day_diff == 0 | ad_day_diff == -1) ~ 0, #index hospitalization is an overdose hospitalization
                                day_diff > 0 ~ ceiling(day_diff/30),
                                day_diff < 0 ~ (ceiling(day_diff/30) - 1))) %>% 
  filter(month_diff >= -12 & month_diff <= 12) #nrow = 12588

#frequency of difference in months
month_od_cn <- month_od %>% count(dc,month_diff)

#bar plot or overdose before/after all discharge
month_od_plot <- month_od_cn %>%
  ggplot(aes(x = month_diff,y = n,fill = dc)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("dodgerblue4","red4"))+
  scale_x_continuous(breaks = c(-12,-9,-6,-3,0,3,6,9,12),
                     labels = c("-12 \n months", "-9 \n months", "-6 \n months", "-3 \n months",
                                "index \n hospitalization",
                                "3 \n months", "6 \n months", "9 \n months", "12 \n months"))+
  labs(x = "Time between index discharge and overdose (months)",
       y = "Number of overdoses") +
  theme_stapleslab+
  theme(legend.position = "none")
month_od_plot

ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od v2.tiff", 
       width=14, height=10, units="in", dpi=600, compression = "lzw",
       plot = month_od_plot)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od v2.pdf", 
       width=14, height=10, units="in", dpi=600, 
       plot = month_od_plot)

#bar plot or overdose before/after discharge BMA
month_od_bma_plot <- month_od_cn %>% 
  filter(dc == "ama") %>% 
  ggplot(aes(x = month_diff,y = n)) +
  geom_bar(stat = "identity",fill = "red4")+
  scale_x_continuous(breaks = c(-12,-9,-6,-3,0,3,6,9,12),
                     labels = c("-12 \n months", "-9 \n months", "-6 \n months", "-3 \n months",
                                "index \n hospitalization",
                                "3 \n months", "6 \n months", "9 \n months", "12 \n months"))+
  labs(x = "Time between index discharge and overdose (months)",
       y = "Number of overdoses") +
  theme_stapleslab+
  theme(legend.position = "none")

month_od_bma_plot
ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od_bma v2.tiff", 
       width=14, height=10, units="in", dpi=600, compression = "lzw",
       plot = month_od_bma_plot)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od_bma v2.pdf", 
       width=14, height=10, units="in", dpi=600, 
       plot = month_od_bma_plot)

#proportion of overdoses before/after DC-BMA
prop_line <- month_od_cn %>% filter(dc == "ama") %>% rename(bma_n = n) %>% select(-dc) %>% 
  left_join(month_od_cn %>% filter(dc == "wa") %>% rename(wa_n = n) %>% select(-dc)) %>% 
  mutate(prop = bma_n/(bma_n+wa_n)*100) %>% 
  ggplot(aes(x = month_diff,y = prop)) +
  geom_point()+
  geom_line()+
  ylim(10,40)+
  scale_x_continuous(breaks = c(-12,-9,-6,-3,0,3,6,9,12),
                     labels = c("-12 \n months", "-9 \n months", "-6 \n months", "-3 \n months",
                                "index \n hospitalization",
                                "3 \n months", "6 \n months", "9 \n months", "12 \n months"))+
  
  labs(x = "Time between index discharge and overdose (months)",
       y = "Proportion of overdoses before or after DC-BMA (%)") +
  theme_stapleslab+
  theme(legend.position = "none")
prop_line
ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od_bma_prop v2.tiff", 
       width=14, height=10, units="in", dpi=600, compression = "lzw",
       plot = prop_line)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od_bma_prop v2.pdf", 
       width=14, height=10, units="in", dpi=600, 
       plot = prop_line)

#bar plot or overdoses before/after all discharge, removing overdoses followed by death within 1 year
month_od_nodeath <- cohort %>% 
  select(-od_date) %>% 
  filter(year(sep_date) > 2015 & year(sep_date) < 2019) %>% #115091 rows
  left_join(od_date %>% select(moh_study_id,od_date)) %>% 
  filter(is.na(death_date) | death_date >= od_date + days(365)) %>% 
  select(dad_id,moh_study_id,sep_date,ad_date, od_date,dc) %>% 
  mutate(day_diff  = as.numeric(od_date - sep_date),
         ad_day_diff = as.numeric(od_date - ad_date),
         month_diff = case_when((ad_day_diff == 0 | ad_day_diff == -1) ~ 0, #index hospitalization is an overdose hospitalization
                                day_diff > 0 ~ ceiling(day_diff/30),
                                day_diff < 0 ~ (ceiling(day_diff/30) - 1))) %>% 
  filter(month_diff >= -12 & month_diff <= 12)

month_od_cn_nodeath <- month_od_nodeath %>% count(dc,month_diff)
month_od_nodeath_plot <- month_od_cn_nodeath %>% 
  ggplot(aes(x = month_diff,y = n,fill = dc)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("dodgerblue4","red4"))+
  scale_x_continuous(breaks = c(-12,-9,-6,-3,0,3,6,9,12),
                     labels = c("-12 \n months", "-9 \n months", "-6 \n months", "-3 \n months",
                                "index \n hospitalization",
                                "3 \n months", "6 \n months", "9 \n months", "12 \n months"))+
  labs(x = "Time between index discharge and overdose (months)",
       y = "Number of overdoses") +
  theme_stapleslab+
  theme(legend.position = "none")
month_od_nodeath_plot

ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od_nodeath v2.tiff", 
       width=14, height=10, units="in", dpi=600, compression = "lzw",
       plot = month_od_nodeath_plot)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od_nodeath v2.pdf", 
       width=14, height=10, units="in", dpi=600, 
       plot = month_od_nodeath_plot)

#bar plot or overdoses before/after DC_BMA, removing overdoses followed by death within 1 year
month_od_bma_nodeath_plot <- month_od_cn_nodeath %>% 
  filter(dc == "ama") %>% 
  ggplot(aes(x = month_diff,y = n)) +
  geom_bar(stat = "identity",fill = "red4")+
  scale_x_continuous(breaks = c(-12,-9,-6,-3,0,3,6,9,12),
                     labels = c("-12 \n months", "-9 \n months", "-6 \n months", "-3 \n months",
                                "index \n hospitalization",
                                "3 \n months", "6 \n months", "9 \n months", "12 \n months"))+
  labs(x = "Time between index discharge and overdose (months)",
       y = "Number of overdoses") +
  theme_stapleslab+
  theme(legend.position = "none")
month_od_bma_nodeath_plot

ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od_nodeath_bma v2.tiff", 
       width=14, height=10, units="in", dpi=600, compression = "lzw",
       plot = month_od_bma_nodeath_plot)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od_nodeath_bma v2.pdf", 
       width=14, height=10, units="in", dpi=600, 
       plot = month_od_bma_nodeath_plot)



#proportion of overdoses before/after DC-BMA, removing overdoses followed by death within 1 year
prop_nodeath <- month_od_cn_nodeath %>% filter(dc == "ama") %>% rename(bma_n = n) %>% select(-dc) %>% 
  left_join(month_od_cn_nodeath %>% filter(dc == "wa") %>% rename(wa_n = n) %>% select(-dc)) %>% 
  mutate(prop = bma_n/(bma_n+wa_n)*100) %>% 
  ggplot(aes(x = month_diff,y = prop)) +
  geom_point()+
  geom_line()+
  ylim(10,40)+
  scale_x_continuous(breaks = c(-12,-9,-6,-3,0,3,6,9,12),
                     labels = c("-12 \n months", "-9 \n months", "-6 \n months", "-3 \n months",
                                "index \n hospitalization",
                                "3 \n months", "6 \n months", "9 \n months", "12 \n months"))+
  
  labs(x = "Time between index discharge and overdose (months)",
       y = "Proportion of overdoses before or after DC-BMA (%)") +
  theme_stapleslab+
  theme(legend.position = "none")
prop_nodeath
ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od_nodeath_bma_prop v2.tiff", 
       width=14, height=10, units="in", dpi=600, compression = "lzw",
       plot = prop_nodeath)
ggsave("R:/working/AMA-OD_cco/NH/results/fig/time_btw_sep_od_nodeath_bma_prop v2.pdf", 
       width=14, height=10, units="in", dpi=600, 
       plot = prop_nodeath)



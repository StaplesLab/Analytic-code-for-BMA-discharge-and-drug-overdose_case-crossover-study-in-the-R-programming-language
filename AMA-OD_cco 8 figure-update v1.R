###########################################################
# This work is licensed under CC BY-NC-SA 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
# Suggested citation: Khan M, Hu X, Yu Y, Daly-Grafstein D, Erdelyi S, Staples JA. AMA-OD_cco sensitivity analyses forest plot update. 2024 May 29. Retrieved from: ***LINK*** 
# Update Figure (#4) forest plot of sensitivity analysis - alternate lookback interval
# Author: Mayesha Khan adapted from Nicole's code [S:/AMA-OD (lab)/6 Manuscripts/2 AMA-OD_cco/Code/
#                                                         AMA-OD_cco 6aii Sensitivity analysis (interval_ctl perior_outcome) wa v1.R] 
# Date: 2024-05-29
# Updated: 2024-
# Additional file needed: S:\AMA-OD (lab)\6 Manuscripts\2 AMA-OD_cco\Code\AMA-OD_cco 8 sensitivity-analysis-fig4.xlsx
##########################################################

# libraries
library(tidyverse)
library(readxl)


##########################################################

theme_stapleslab <- theme_bw() +
  theme(
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "grey30"),
    panel.border = element_rect(size = 1.25, colour = "grey30", fill = NA),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, colour = "grey30"),
    axis.ticks = element_line(size = 0.5, colour = "grey30"),
    axis.title.x = element_text(vjust= -0.5),
    axis.title.y = element_text(vjust= +1.75)#,
#    plot.margin = margin(t=20,r=20,b=20,l=20)
)


##########################################################
# read in results
adj_models <- read_excel("S:/AMA-OD (lab)/6 Manuscripts/2 AMA-OD_cco/Code/AMA-OD_cco 8 sensitivity-analysis-fig4.xlsx")

# get (robust) SE
adj_models <- adj_models %>% 
  mutate(robust.se = (log(uci) - log(aor))/qnorm(0.975))
  
dw = 0.5
adj_models <- adj_models %>% 
  mutate(xdodge = ifelse(dc == "bma", as.numeric(lookbk) - dw,  as.numeric(lookbk) + dw))

plot.ts.sub <- ggplot(adj_models,aes(x = xdodge, y = aor, color = dc)) +
  geom_point(aes(size = (1/robust.se)^2), shape = "square") +
  geom_segment(aes(x = xdodge, xend = xdodge, y = lci, yend = uci), position = position_dodge(width = dw))  + 
   geom_hline(yintercept = 1, size = 0.25)+ 
  scale_y_continuous(trans = "log",
                     breaks = c(0.5, 1.0, 2.0, 3.5),  #orig: c(0, 0.7, 1, 1.5, 2, 2.5), # jas: c(0.5, 1.0, 2.0, 4.0)
                     labels = c(0.5, 1.0, 2.0, 3.5)) +
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


ggsave("C:/Users/maykhan/Desktop/Staples lab/AMA-OD/2 cco/cco_interval_length_forest_plot_all v3.tiff", 
       width=14, height=8, units="in", dpi=600, compression = "lzw",
       plot = plot.ts.sub)
ggsave("C:/Users/maykhan/Desktop/Staples lab/AMA-OD/2 cco/cco_interval_length_forest_plot_all v3.pdf", 
       width=14, height=8, units="in", dpi=600, 
       plot = plot.ts.sub,
       device = cairo_pdf)

###############################################################################
##                                                                           ##
## 03: Descriptive Analysis Before APC Bounding and APC-I Estimation         ##
##                                                                           ##
## Inputs:                                                                   ##
##   1) data/cps.76.20.recoded.wt.voted.rds                                  ##
##   2) data/cps.76.20.recoded.wt.all.rds                                    ##
##   3) rawdata/mcdonald_vep_turnout_1789_present.csv                        ##
##      (from US Elections Project; see https://www.electproject.org/)       ##
##                                                                           ##
## Notes:                                                                    ##
##   - Generate descriptive tables/figures (specifically, Figures 1, 2, 3,   ##
##     S1(a), and Table 1)                                                   ##
##                                                                           ##
## Final Data Output:                                                        ##
##   - None                                                                  ##
##                                                                           ##
## Final Table/Visual Outputs:                                               ##
##   1) output/tex/apc_freq.tex (Table S1)                                   ##
##   2) output/figure/heatmap_apc.pdf (Figure 1)                             ##
##   3) output/figure/non_resp_by_year.pdf (Figure S1(a))                    ##
##   4) output/figure/turnout_trends.pdf (Figure S1(b))                      ##
##   5) output/figure/C_by_A.pdf (Figure 2)                                  ##
##   6) output/figure/C_by_P.pdf (Figure 3)                                  ##
##                                                                           ##
###############################################################################




#=======================================
# Preliminaries
#=======================================

rm(list=ls())

##########
##########    LIBRARIES
##########

# packages
packages <- c("tidyverse", "survey", "srvyr", "here", "ggplot2", "scales", "xtable", "viridis", "viridisLite", "RColorBrewer", "pals"  )

# install if necessary
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN') 
  }
}

# library
for (package in packages) {
  library(package, character.only=T)
}




#=======================================
# Call-in data and set-up for analysis
#=======================================

# call-in data ------------------------------------------------------------
message('Loading data: `cps.76.20.recoded.wt.voted.rds` ...')

here("data", "cps.76.20.recoded.wt.voted.rds") %>% readRDS() -> cps.76.20 
here("data", "cps.76.20.recoded.wt.all.rds") %>% readRDS() -> cps.76.20.all 

message('Loading data: Actual VEP Turnout estimates ...')

vep.turnout <- read.csv(here("rawdata", "mcdonald_vep_turnout_1789_present.csv"))  %>% 
  filter(P >= 1976)
vep.turnout$P <- as.factor(vep.turnout$P)
vep.turnout$vep_turnout <- vep.turnout$vep_turnout 

#=======================================
# Descriptive Analysis
#=======================================

# A x P cell Ns ------------------------------------------------------------
print(xtable(table(cps.76.20$P, cps.76.20$A), type = "latex", file =  here("output", "tex", "apc_freq.tex"))) # latex

# headtmap of APC and turnout ----------------------------------------------
cps.76.20 %>% 
  mutate(voted = voted * 100) %>%  # rescale the outcome
  as_survey(weights = c(wt.miss.rs)) %>%
  group_by(P, A) %>% 
  summarise(tr = survey_mean(voted)) %>%
  # summarise(tr = mean(voted)) %>%
  ggplot(aes(x = P, y = A)) +
  geom_tile(aes(fill = tr), colour="white", size=0.25)  + 
  labs(x="\nElection Year",y="Age\n", title="" , size =12) +
  guides(fill=guide_legend(title="Turnout %", reverse = TRUE)) +
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0)) +
  geom_text(aes(label = format(round(tr, 1), color="black")),  size = 3.3)   +
  # scale_fill_viridis(option="rainbow", trans = 'reverse', alpha = 0.8) +
  # scale_fill_manual(   guide = guide_legend(reverse = TRUE)) +
  scale_fill_gradientn(colors = alpha(c(brewer.pal(9, "Spectral")), alpha = 1)) +
  # scale_fill_gradient(low="white", high="blue") +
  # scale_fill_distiller(palette = "RdBu",trans = "reverse") +
  theme_grey(base_size=10) +
  theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.5),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    axis.text = element_text(size = 8, color = "black"),
    axis.title= element_text(size= 11, face="bold"))

ggsave(path = here("output", "figures"), filename = "heatmap_apc.pdf", dpi = 1000, device = cairo_pdf, width = 8, height = 3, units = "in")


# Percent of non-response by year ----------------------------------------------

cps.76.20.all %>% 
  as_survey(weights = c(wt)) %>%
  group_by(P) %>% 
  summarise(pr.miss = survey_mean(is.na(voted)) * 100) 


cps.76.20.all %>% 
  as_survey(weights = c(wt)) %>%
  group_by(P) %>% 
  summarise(pr.miss = survey_mean(is.na(voted)) *100 ) %>% 
  ggplot(., aes(x = P, y = pr.miss, group = 1)) +
  geom_line() +
  geom_point()  +
  theme_bw() +   
  labs(title = "", x = "\nElection Year", y = "Proportion of Non-response\n") +
  ylim(0,20) +
  theme(
      legend.text = element_text(size = 11),
      axis.ticks=element_line(size=0.4),
      axis.text = element_text(size = 9, color = "black"),
      axis.title= element_text(size= 12)) 
           
ggsave(path = here("output", "figures"), filename = "non_resp_by_year.pdf", dpi = 1000, device = cairo_pdf, width =13, height =8, units = "cm")


# Marginal turnout rate by period -----------------------------------------

# using iptw
turnout.by.p.wt <- cps.76.20 %>% 
  mutate(voted = voted * 100) %>%  # rescale the outcome
  as_survey(weights = c(wt.miss.rs)) %>%
  group_by(P) %>% 
  summarise(pr.voted.iptw = survey_mean(voted)) 

# using base weights
turnout.by.p.base <- cps.76.20 %>% 
  mutate(voted = voted * 100) %>%  # rescale the outcome
  as_survey(weights = c(wt)) %>%
  group_by(P) %>% 
  summarise(pr.voted.base = survey_mean(voted)) 


# plot
tr.rates <- turnout.by.p.wt %>%
  full_join(., turnout.by.p.base, by = "P")  %>% 
  full_join(., vep.turnout, by = "P") %>% 
  select(pr.voted.iptw, pr.voted.base, vep_turnout, P) %>% 
  reshape2::melt(., id.vars = "P", value.anem = "value", variable.name = "type")

ggplot(data=tr.rates, aes(x=P, y=value, group = type, color = type)) +
  geom_line(aes(linetype=type, color=type), size= 0.8) +
  geom_point(aes(shape = type), size=1.5) +
  scale_color_brewer(palette = "Set1", label = c("CPS (IPTW)", "CPS", "VEP")) +
  scale_linetype_manual(values=c("dotdash", "solid", "longdash")) +
  theme_bw() +   
  labs(title = "", x = "\nElection Year", y = "Turnout %\n") +
  theme(legend.title = element_blank()) + 
  scale_size_discrete() +
  guides(linetype = FALSE, shape = FALSE) +
  theme(legend.position=c(0.3, 0.9), legend.direction = "horizontal",
        legend.background=element_blank(),
        legend.key = element_blank(),
        text = element_text(family = ),
        legend.text = element_text(size = 11),
        axis.ticks=element_line(size=0.4),
        axis.text = element_text(size = 9, color = "black"),
        axis.title= element_text(size = 12)) +
  ylim(40,90) 

ggsave(path = here("output", "figures"), filename = "turnout_trends.pdf", dpi = 1000, device = cairo_pdf, width =13, height =8, units = "cm")

# difference between VEP and CPS
tr.rates[1:12,3] - tr.rates[25:36,3] 



# Cohort Plot ------------------------------------------------------------

# cohort data using 10-year bins for age groups
cohort <- cps.76.20 %>% 
  mutate(age10 = cut_width(age, 10, closed = "left", boundary = 18)) 
cohort$age10 <- recode(cohort$age10, "[68,78]" = "[68,78)")  
table(cohort$age10, cohort$age, useNA = "ifany")
cohort$age10 <- factor(paste0(substr(cohort$age10, 2 ,3),"-", as.integer(substr(cohort$age10, 5,6)) - 1))

cohort %>%
  mutate(voted = voted * 100) %>%  # rescale the outcome
  as_survey(weights = c(wt.miss)) %>%
  group_by(age10, C) %>%
  summarise(pr.tr = survey_mean(voted))



cohort %>% 
  mutate(voted = voted * 100) %>%  # rescale the outcome
  as_survey(weights = c(wt.miss)) %>%
  group_by(A, C) %>% 
  summarise(pr.tr = survey_mean(voted)) %>% 
  ggplot(., aes(C, pr.tr, group = A, color = A, shape = A, linetype = A)) +
  labs(title = "", x = "\nBirth Year", y = "Turnout %\n", color = "Age Groups") +
  geom_point(size = 2) +
  geom_line(linewidth = 0.8) +
  theme_bw() +
  scale_shape_manual(values=c(15,16,7,18,19,8,15,16,17,4,19,20,5,1,19)) +
  scale_linetype_manual(values=c(rep(c(1,2),7), 1)) +
  scale_x_discrete(breaks=c("1899-1902", "1911-1914", "1923-1926", "1935-1938", "1947-1950", 
                            "1959-1962", "1971-1974", "1983-1986", "1995-1998")) +
  theme(
    legend.text = element_text(size = 11),
    axis.ticks=element_line(size=0.4),
    axis.text = element_text(size = 11, color = "black"),
    axis.title= element_text(size = 13)) +
  scale_colour_manual(values=c('#a50026','#d73027','#f46d43','#fdae61','#fee010','#8e0152','#c51b7d','#4393c3', '#4575b4', '#2166ac','#313695', 
                               '#b8e186','#7fbc41','#4d9221','#276419')) +
  # scale_colour_brewer(palette = "Spectral") +
  labs(color  = "Age Groups", linetype = "Age Groups", shape = "Age Groups", nrow = 2) +
  guides(
    color = guide_legend(ncol = 2),
    shape = guide_legend(ncol = 2,),
    linetype = guide_legend(ncol = 2)) +  
     scale_y_continuous(limits = c(30, 90), breaks = seq(30, 90, by = 10)) 

 ggsave(path = here("output", "figures"), filename = "C_by_A.pdf", dpi = 1000, device = cairo_pdf, width = 25, height =8, units = "cm")

# C x P plot
 cohort %>%
   mutate(voted = voted * 100) %>%  # rescale the outcome
   as_survey(weights = c(wt.miss)) %>%
   group_by(P, C) %>%
   summarise(pr.tr = survey_mean(voted)) %>%
   ggplot(., aes(C, pr.tr, group = P, color = P, shape = P, linetype = P)) +
   labs(title = "", x = "\nBirth Year", y = "Turnout %\n", color = "Election Year") +
   geom_point(size = 2) +
   geom_line(linewidth = 0.8) +
   theme_bw() +
   scale_shape_manual(values=c(15,16,1,17,5,15,6,16,17,4,19,2)) +
   scale_linetype_manual(values=c(rep(c(1,2),6))) +
   scale_x_discrete(breaks=c("1899-1902", "1911-1914", "1923-1926", "1935-1938", "1947-1950", 
                             "1959-1962", "1971-1974", "1983-1986", "1995-1998")) +
   theme(
     legend.text = element_text(size = 11),
     axis.ticks=element_line(size=0.4),
     axis.text = element_text(size = 11, color = "black"),
     axis.title= element_text(size = 13)) +
   # scale_colour_viridis_d(option = "turbo") +
   scale_colour_manual(values=c('#762a83','#9970ab','#c2a5cf','#e7d4e8','#d6604d','#f4a582','#a6dba0','#5aae61','#1b7837','#74add1','#4575b4', '#313695')) +
   scale_y_continuous(limits = c(30, 90), breaks = seq(30, 90, by = 10)) +
   labs(color  = "Election Year", linetype = "Election Year", shape = "Election Year", nrow = 2) +
   guides(
     color = guide_legend(ncol = 2),
     shape = guide_legend(ncol = 2,),
     linetype = guide_legend(ncol = 2))

ggsave(path = here("output", "figures"), filename = "C_by_P.pdf", dpi = 1000, device = cairo_pdf, width = 25, height =8, units = "cm")





# P x C plot
# cohort %>%
#   mutate(voted = voted * 100) %>%  # rescale the outcome
#   as_survey(weights = c(wt.miss)) %>%
#   group_by(P, C) %>%
#   summarise(pr.tr = survey_mean(voted)) %>%
#   ggplot(., aes(P, pr.tr, group = C, color = C)) +
#   labs(title = "", x = "\nBirth Year", y = "Turnout Rate\n", color = "Cohort") +
#   geom_point() +
#   geom_line() +
#   theme_bw() +
#   # scale_x_discrete(breaks=c("1898-1901", "1910-1913", "1922-1925", "1934-1937", "1946-1949",
#   #                           "1958-1961", "1970-1973", "1982-1985", "1994-1997")) +
#   theme(
#     legend.text = element_text(size = 10),
#     axis.ticks=element_line(size=0.4),
#     axis.text = element_text(size = 9, color = "black"),
#     axis.title= element_text(size = 12)) +
#   scale_colour_viridis_d(option = "magma") +
#   scale_y_continuous(limits = c(30, 100), breaks = seq(30, 90, by = 10))
# 
# df <- data.frame(name = c("New-Deal Generation", "Boomers"),
#                  start = c("[1918,1922)", "[1946,1950)" ),
#                  end = c("[1930,1934)", "[1962,1966)"),
#                  C = c("ND", "BB"),
#                  stringsAsFactors = FALSE)
# 
# cohort %>%
#   group_by(age10, C) %>%
#   summarise(pr.tr = mean(voted)) %>%
#   ggplot(., aes(C, pr.tr, group = age10, color = age10)) +
#   labs(title = "", x = "\nCohort", y = "Turnout Rate\n", color = "Age Groups") +
#   geom_point() +
#   geom_line() +
#   theme_bw() +
#   scale_x_discrete(breaks=c("[1898,1902)", "[1910,1914)", "[1922,1926)", "[1934,1938)", "[1946,1950)",
#                             "[1958,1962)", "[1970,1974)", "[1982,1986)", "[1994,1998)")) +
#   geom_rect(data=df, aes(ymin= 0.4, ymax=0.9, xmin=start ,xmax= end, fill = C), colour= "white", size=0.5, alpha=0.2) +
#   scale_fill_manual(values=c("ND" = "red", "BB" = "blue"))
# 
# 
# cps.76.20 %>%
#   filter((cohort >= 1918 & cohort <= 1933) | (cohort >= 1946 & cohort <= 1961)) %>%
#   group_by(A, C) %>%
#   summarise(pr.tr = mean(voted)) %>%
#   ggplot(., aes(A, pr.tr, group = C, color = C)) +
#   # labs(title = "", x = "\nCohort", y = "Turnout Rate\n", color = "Age Groups") +
#   geom_point() +
#   geom_line() +
#   theme_bw() +
#   # scale_x_discrete(breaks=c("[1898,1902)", "[1910,1914)", "[1922,1926)", "[1934,1938)", "[1946,1950)",
#   #                           "[1958,1962)", "[1970,1974)", "[1982,1986)", "[1994,1998)")) +
#   theme(
#     legend.text = element_text(size = 10),
#     axis.ticks=element_line(size=0.4),
#     axis.text = element_text(size = 9, color = "black"),
#     axis.title= element_text(size = 12)) +
#   scale_colour_viridis_d(option = "viridis") +
#   ylim(0.4,0.9)
# 
# ggsave(path =here("output", "figures"), filename = "C_by_A.pdf", dpi = 1000, device = cairo_pdf, width = 25, height =8, units = "cm")


#=======================================
# Print Session Info
#=======================================


cat("\n\n---- Session Info ----\n")
print(sessionInfo())

### END OF CODE ###
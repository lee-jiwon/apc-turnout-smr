###############################################################################
##                                                                           ##
## 04: APC Bounding Analysis and APC-I Model Estimation                      ##
##                                                                           ##   
## Inputs:                                                                   ##
##   1) data/cps.76.20.recoded.wt.rds                                        ##
##                                                                           ##
## Notes:                                                                    ##
##   - Estimate APC bounding analysis and APC-I model                        ##
##     using linear probability models                                       ##
##                                                                           ##
## Final Data Output:                                                        ##
##   - None                                                                  ##
##                                                                           ##
## Final Table/Visual Outputs:                                               ##
##   1) output/tex/apci_all.tex (Table S3)                                   ##
##   2) output/tex/bounding_A_nonlin.tex (Table S7)                          ##
##   3) output/tex/bounding_P_nonlin.tex (Table S7)                          ##
##   4) output/tex/bounding_C_nonlin.tex (Table S7)                          ##
##   5) output/tex/bounding_A.tex (Table S3)                                 ##
##   6) output/tex/bounding_P.tex (Table S3)                                 ##
##   7) output/tex/bounding_C.tex (Table S3)                                 ##
##   8) output/figures/nonlinear_apc.pdf (Figure 4)                          ##
##   9) output/figures/total_apc_a12_with_api.pdf (Figure 5(a))              ##
##  10) output/figures/total_apc_a12_c1_with_api.pdf (Figure 5(b))           ##
##  11) output/figures/total_apc_a12_c1_with_api.pdf (Figure 5(c))           ##
##  12) output/figures/intra_cohort_apci.pdf (Figure S3)                     ##
##                                                                           ##
###############################################################################


 
#=======================================
# Preliminaries
#=======================================

# clear  console
rm(list=ls())


# packages
packages <- c("tidyverse", "survey", "srvyr", "weightedcontrasts", 
              "ggplot2", "ggpubr", "xtable", "sandwich", "APCI", "here")

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
# Data and Functions
#=======================================


# call-in data ------------------------------------------------------------
message('Loading data: `cps.76.20.recoded.wt.rds` ...')
cps.76.20 <- readRDS(here("data", "cps.76.20.recoded.wt.voted.rds")) 

# call-in functions -------------------------------------------------------
source(here("script", "00_functions.R"))


#=======================================
# Estimate APC-I model
#=======================================

# Estimate the model
apci.model <- get_apci_identity("voted", "A", "P", "C", "wt.miss", cps.76.20)

# Define the ranges
A.range <- seq(18, 74, 4)
P.range <- seq(1976, 2020, 4)
C.range <- seq(1899, 2002, 4)

# Combine the estimates into data frames
comb.est <- function(range, est_col) {
  as.data.frame(cbind(range, est_col[ ,1], est_col[ ,2], est_col[ ,3]))
}

apci.A.est <- comb.est(A.range, apci.model$age_est[, c("A.main", "A.main.ci.lb", "A.main.ci.ub")])
apci.P.est <- comb.est(P.range, apci.model$period_est[, c("P.main", "P.main.ci.lb", "P.main.ci.ub")])
apci.C.est <- comb.est(C.range, apci.model$cohort_est[, c("C.dev", "C.dev.ci.lb", "C.dev.ci.ub")])
apci.C.slope <- comb.est(C.range, apci.model$cohort_slope_st[, c("C.slope", "C.slope.ci.lb", "C.slope.ci.ub")])


#------------------------------------------------
# export the estimates to latex 
# -----------------------------------------------

A <- cbind(cbind(levels(cps.76.20$A), 
                 paste0(round(apci.model$age_est[,1], 1)," (", round(-(apci.model$age_est[,2] - apci.model$age_est[,1])/(1.96), 1), ")" )))

P <- cbind(cbind(levels(cps.76.20$P), 
                 paste0(round(apci.model$period_est[,1], 1)," (", round(-(apci.model$period_est[,2] - apci.model$period_est[,1])/(1.96), 1), ")" )))

C <- cbind(cbind(levels(cps.76.20$C), 
                 paste0(round(apci.model$cohort_est[,1], 1)," (", round(-(apci.model$cohort_est[,2] - apci.model$cohort_est[,1])/(1.96), 1), ")" )))

print(xtable(rbind(A,P,C), type = "latex", file = here("output", "text", "apci_all.tex"), digits = 1), include.rownames=FALSE)


#=======================================
# APC Bounding Analysis
#=======================================

# polynomial contrast coding
cps.76.20 <- tibble(cps.76.20)
  
contrasts(cps.76.20$A) <- contr.poly.weighted(cps.76.20$A, width = 4)
contrasts(cps.76.20$P) <- contr.poly.weighted(cps.76.20$P, width = 4)
contrasts(cps.76.20$C) <- contr.poly.weighted(cps.76.20$C, width = 4)
contrasts(cps.76.20$P)[, 1] <- 0 # zero out the period linear effect

# declare survey design for sampling weight 
cps.76.20$voted.pr <- cps.76.20$voted * 100 # re-scale the outcome
survey.design <- svydesign(id = ~1, weights = ~wt.miss, data = cps.76.20)  

# estimate the model: linear-probability model
model <- svyglm(voted.pr ~ A + C + P, design = survey.design, family = "gaussian")

# get deviations from linear trend and age/cohort linear effects
A.nonlin <- get_deviations_w_intercept(model, "^A", contrasts(cps.76.20$A)) 
C.nonlin <- get_deviations_w_intercept(model, "^C", contrasts(cps.76.20$C)) 

extract_coefs <- c(0, coef(model)[grepl("^P", names(coef(model)))])
P.nonlin <-  coef(model)[1] + contrasts(cps.76.20$P)  %*% extract_coefs

# check the outputs
A.nonlin 
P.nonlin
C.nonlin

# export the non-linear effect estimates to latex
print(xtable((data.frame(A.nonlin)), type = "latex", file = here("output", "tex", "bounding_A_nonlin.tex"), digits = 1))
print(xtable((data.frame(P.nonlin)), type = "latex", file = here("output", "tex", "bounding_P_nonlin.tex"), digits = 1))
print(xtable((data.frame(C.nonlin)), type = "latex", file = here("output", "tex", "bounding_C_nonlin.tex"), digits = 1))

#------------------------------------------------
# Derive the bounds 
#-----------------------------------------------

# the sums and differences between the linear terms
theta1 <- coef(model)["A.L"] 
theta2 <- coef(model)["C.L"] 

# Assumptions A1 and A2: turnout increases monotonically from 18 to 61 onwards and turnout declines monotonically from 64 onwards
min.alpha.a12 <- -min(diff(A.nonlin[1:11]))/4
max.alpha.a12 <- -max(diff(A.nonlin[14:15]))/4
min.pi.a12 <- theta1 - max.alpha.a12
max.pi.a12 <- theta1 - min.alpha.a12
min.gamma.a12 <- theta2 - theta1 + min.alpha.a12
max.gamma.a12 <- theta2 - theta1 + max.alpha.a12

# check the bounds
message('<Identification regegions based on Assumptions A1 + A2>')
print.boounds(c(min.alpha.a12, max.alpha.a12, min.pi.a12, 
                max.pi.a12, min.gamma.a12, max.gamma.a12))


# Assumption c1: linear cohort effect is non-negative
min.gamma.c1 <- 0
max.gamma.c1 <- Inf
min.alpha.c1 <- theta1 - theta2 + min.gamma.c1
max.alpha.c1 <- theta1 - theta2 + max.gamma.c1
min.pi.c1 <- theta2 - max.gamma.c1
max.pi.c1 <- theta2 - min.gamma.c1

# check the bounds
message('<Identification regegions based on Assumptions C1>')
print.boounds(c(min.alpha.c1, max.alpha.c1, min.pi.c1, max.pi.c1,min.gamma.c1, max.gamma.c1 ))


# Assumption C2: Cohort effects of new-deal generation >= cohort effects of baby-boomers
wt_N <- count(cps.76.20, C, wt = wt.miss) # weighted counts of each cohort
wt.nd <- wt_N[4:7,2]  / sum(wt_N[4:7,2]) # weighted counts of new-deal generation
wt.bb <- wt_N[13:17,2]  / sum(wt_N[13:17,2]) # weighted counts of baby-boomers

weighted.avg.nd <- wt.nd[1,1] + 4 * wt.nd[2,1] + 8 * wt.nd[3,1] + 12 * wt.nd[4,1] 
weighted.avg.bb <- 36 * wt.bb[1,1] + 40 * wt.bb[2,1] + 44 * wt.bb[3,1] + 48 * wt.bb[4,1] + 52 * wt.bb[5,1]

max.gamma.c2 <- (wt.bb[1,1] * C.nonlin[13] + wt.bb[2,1] * C.nonlin[14] + wt.bb[3,1] * C.nonlin[15] + wt.bb[4,1] * C.nonlin[16] +  wt.bb[5,1] * C.nonlin[17] -
  wt.nd[1,1] * C.nonlin[4] - wt.nd[2,1] * C.nonlin[5] - wt.nd[3,1] * C.nonlin[6] - wt.nd[4,1] * C.nonlin[7])  / (weighted.avg.nd - weighted.avg.bb)
min.gamma.c2 <- -Inf

min.alpha.c2 <- theta1 - theta2 + min.gamma.c2
max.alpha.c2 <- theta1 - theta2 + max.gamma.c2
min.pi.c2 <- theta2 - max.gamma.c2
max.pi.c2 <- theta2 - min.gamma.c2

# check the bounds
message('<Identification regegions based on Assumptions C2>')
print.boounds(c(min.alpha.c2, max.alpha.c2, min.pi.c2, max.pi.c2,min.gamma.c2, max.gamma.c2 ))

#------------------------------------------------
# Identification regions of the total effects 
#-----------------------------------------------

# polynomial contrast coding
contrasts(cps.76.20$P) <- contr.poly.weighted(cps.76.20$P, width = 4)

# Assumptions A1 and A2 only
# drive the region of the total effects
cohort.lb.a12 <- C.nonlin + contrasts(cps.76.20$C)[, 1] * min.gamma.a12
cohort.ub.a12 <- C.nonlin + contrasts(cps.76.20$C)[, 1] * max.gamma.a12
age.lb.a12 <- A.nonlin + contrasts(cps.76.20$A)[, 1] * min.alpha.a12
age.ub.a12 <- A.nonlin + contrasts(cps.76.20$A)[, 1] * max.alpha.a12 
period.lb.a12 <- P.nonlin + contrasts(cps.76.20$P)[, 1] * min.pi.a12 
period.ub.a12 <- P.nonlin + contrasts(cps.76.20$P)[, 1] * max.pi.a12 

# Assumptions A1, A2, and C1
# assigning upper/lower bounds
min.alpha.a12.c1 <- max(min.alpha.a12, min.alpha.c1)
max.alpha.a12.c1 <- min(max.alpha.a12, max.alpha.c1)
min.pi.a12.c1 <- max(min.pi.a12, min.pi.c1)
max.pi.a12.c1 <- min(max.pi.a12, max.pi.c1)  
min.gamma.a12.c1 <- max(min.gamma.a12, min.gamma.c1)
max.gamma.a12.c1 <- min(max.gamma.a12, max.gamma.c1)

# check the bounds
message('<Identification regegions based on Assumptions A1 + A2 + C1>')
print.boounds(c(min.alpha.a12.c1, max.alpha.a12.c1, min.pi.a12.c1, max.pi.a12.c1, min.gamma.a12.c1, max.gamma.a12.c1))

# drive the region of the total effects
cohort.lb.a12.c1 <- C.nonlin + contrasts(cps.76.20$C)[, 1] * min.gamma.a12.c1
cohort.ub.a12.c1 <- C.nonlin + contrasts(cps.76.20$C)[, 1] * max.gamma.a12.c1
age.lb.a12.c1 <- A.nonlin + contrasts(cps.76.20$A)[, 1] * min.alpha.a12.c1
age.ub.a12.c1 <- A.nonlin + contrasts(cps.76.20$A)[, 1] * max.alpha.a12.c1
period.lb.a12.c1 <- P.nonlin + contrasts(cps.76.20$P)[, 1] * min.pi.a12.c1
period.ub.a12.c1 <- P.nonlin + contrasts(cps.76.20$P)[, 1] * max.pi.a12.c1


# Assumptions A1, A2, and C2
# assigning upper/lower bounds
min.alpha.a12.c2 <- max(min.alpha.a12, min.alpha.c2)
max.alpha.a12.c2 <- min(max.alpha.a12, max.alpha.c2)
min.pi.a12.c2 <- max(min.pi.a12, min.pi.c2)
max.pi.a12.c2 <- min(max.pi.a12, max.pi.c2)
min.gamma.a12.c2 <- max(min.gamma.a12,  min.gamma.c2)
max.gamma.a12.c2 <- min(max.gamma.a12,  max.gamma.c2)

# check the bounds
message('<Identification regegions based on Assumptions A1 + A2 + C2>')
print.boounds(c(min.alpha.a12.c2, max.alpha.a12.c2, min.pi.a12.c2, 
                max.pi.a12.c2, min.gamma.a12.c2, max.gamma.a12.c2))


# drive the region of the total effects
cohort.lb.a12.c2 <- C.nonlin + contrasts(cps.76.20$C)[, 1] * min.gamma.a12.c2 
cohort.ub.a12.c2  <- C.nonlin + contrasts(cps.76.20$C)[, 1] * max.gamma.a12.c2 
age.lb.a12.c2  <- A.nonlin + contrasts(cps.76.20$A)[, 1] * min.alpha.a12.c2 
age.ub.a12.c2  <- A.nonlin + contrasts(cps.76.20$A)[, 1] * max.alpha.a12.c2 
period.lb.a12.c2  <- P.nonlin + contrasts(cps.76.20$P)[, 1] * min.pi.a12.c2 
period.ub.a12.c2  <- P.nonlin + contrasts(cps.76.20$P)[, 1] * max.pi.a12.c2 

# Assumptions A1, A2, C1 and C2
# assigning upper/lower bounds
min.alpha.a12.c12 <- max(min.alpha.a12, min.alpha.c1, min.alpha.c2)
max.alpha.a12.c12 <- min(max.alpha.a12, max.alpha.c1, max.alpha.c2)
min.pi.a12.c12 <- max(min.pi.a12, min.pi.c1, min.pi.c2)
max.pi.a12.c12 <- min(max.pi.a12, max.pi.c1, max.pi.c2)
min.gamma.a12.c12 <- max(min.gamma.a12, min.gamma.c1,  min.gamma.c2)
max.gamma.a12.c12 <- min(max.gamma.a12, max.gamma.c1, max.gamma.c2)

# check the bounds
message('<Identification regegions based on Assumptions A1 + A2 + C1 + C2>')
print.boounds(c(min.alpha.a12.c12, max.alpha.a12.c12, min.pi.a12.c12, 
                max.pi.a12.c12, min.gamma.a12.c12, max.gamma.a12.c12))

# drive the region of the total effects
cohort.lb.a12.c12 <- C.nonlin + contrasts(cps.76.20$C)[, 1] * min.gamma.a12.c12 
cohort.ub.a12.c12  <- C.nonlin + contrasts(cps.76.20$C)[, 1] * max.gamma.a12.c12 
age.lb.a12.c12  <- A.nonlin + contrasts(cps.76.20$A)[, 1] * min.alpha.a12.c12 
age.ub.a12.c12  <- A.nonlin + contrasts(cps.76.20$A)[, 1] * max.alpha.a12.c12 
period.lb.a12.c12  <- P.nonlin + contrasts(cps.76.20$P)[, 1] * min.pi.a12.c12 
period.ub.a12.c12  <- P.nonlin + contrasts(cps.76.20$P)[, 1] * max.pi.a12.c12 


#=======================================
# Export the estimates to latex
#=======================================

# A helper function
lb_ub <- function(ub , lb) {
  lb_temp <- apply(cbind(ub, lb),1,min,na.rm=TRUE)
  ub_temp <- apply(cbind(ub, lb),1,max,na.rm=TRUE)
  bounds <- paste0("[", substr(round(lb_temp,1), 1 ,4), ", ", substr(round(ub_temp,1), 1,4), "]")
  # bounds <- cbind(rownames(data.frame(ub)), bounds)
  return(bounds)
}

A.a12 <- lb_ub(age.ub.a12, age.lb.a12)
A.a12.c1 <- lb_ub(age.ub.a12.c1, age.lb.a12.c1)
A.a12.c2 <- lb_ub(age.ub.a12.c2, age.lb.a12.c2)
A.a12.c12 <- lb_ub(age.ub.a12.c12, age.lb.a12.c12)

P.a12 <- lb_ub(period.ub.a12, period.lb.a12)
P.a12.c1 <- lb_ub(period.ub.a12.c1, period.lb.a12.c1)
P.a12.c2 <- lb_ub(period.ub.a12.c2, period.lb.a12.c2)
P.a12.c12 <- lb_ub(period.ub.a12.c12, period.lb.a12.c12)

C.a12 <- lb_ub(cohort.ub.a12, cohort.lb.a12)
C.a12.c1 <- lb_ub(cohort.ub.a12.c1, cohort.lb.a12.c1)
C.a12.c2 <- lb_ub(cohort.ub.a12.c2, cohort.lb.a12.c2)
C.a12.c12 <- lb_ub(cohort.ub.a12.c12, cohort.lb.a12.c12)

print(xtable(cbind(rownames(data.frame(age.ub.a12)), A.a12, A.a12.c1, A.a12.c2), type = "latex",
             file = here("output", "tex", "bounding_A.tex"), digits = 1), include.rownames=FALSE)
print(xtable(cbind(rownames(data.frame(period.ub.a12)), P.a12, P.a12.c1, P.a12.c2), type = "latex",
             file = here("output", "tex", "bounding_P.tex"), digits = 1), include.rownames=FALSE)
print(xtable(cbind(rownames(data.frame(cohort.ub.a12)), C.a12, C.a12.c1, C.a12.c2), type = "latex",
             file =  here("output", "tex", "bounding_C.tex"), digits = 1), include.rownames=FALSE)



#=======================================
# Plot the results
#=======================================

# estimates into df
df.age.b <- as.data.frame(cbind(A.range, age.ub.a12, age.lb.a12, A.nonlin))  
df.period..b <-  as.data.frame(cbind(P.range, period.ub.a12, period.lb.a12, P.nonlin))
df.cohort.b <-  as.data.frame(cbind(C.range, cohort.ub.a12, cohort.lb.a12, C.nonlin))

# Non-linear effects ------------------------------------------------------

# data for rectangle to indicate identifying information
A.rect <- data.frame(xmin = c(54, 70), xmax = c(58, 74), ymin = c(-Inf, -Inf), ymax = c(Inf, Inf), group = c("A1", "A2"))
C.rect <- data.frame(xmin = c(1911, 1947, 1899), xmax = c(1926, 1966, 2001),  
                     ymin = c(-Inf, -Inf,-Inf), ymax = c(Inf, Inf,Inf), 
                     color = c("c1", "c2", "c3"),
                     group = c("A1", "A1", "A2"))

A.nolin.g <- ggplot(df.age.b, aes(df.age.b[,1], A.nonlin)) + 
  geom_line(color = "black") +
  geom_point() +
  theme_bw() + 
  labs(title = "Nonlinear Age Estimates", x = "Age", y = "Turnout %") +
  scale_x_continuous(limits = c(18, 76), breaks = c(20, 30, 40, 50, 60, 70 , 74)) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 9, color = "black"),
        axis.title= element_text(size= 11),
        plot.title = element_text(size=11)) +
  geom_rect(data=A.rect, inherit.aes = F, 
            mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=group)) + 
  scale_fill_manual(values = alpha(c("blue3", "red3"), 0.15)) +
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "black", alpha=0.5) +
  # geom_text(data=A.rect, aes(x=c(56, 72), y=c(45, 45), label = c("A1", "A2")), color = c("blue", "red"), size=4) +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) 

P.nolin.g <- ggplot(df.period..b, aes(df.period..b[,1], P.nonlin)) + 
  geom_line(color = "black") +
  geom_point(data = df.period..b, aes(df.period..b[,1], P.nonlin)) +
  theme_bw() + 
  labs(title = "Nonlinear Period Estimates", x = "Election Year", y = "") +
  scale_x_continuous(limits = c(1976, 2020), breaks = c(seq(1980, 2020, by =10))) +
  theme( plot.title = element_text(size=11),  
         axis.text = element_text(size = 9, color = "black"),
         axis.title= element_text(size= 11)) +
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "black", alpha=0.5) +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90))  

C.nolin.g <- ggplot(df.cohort.b, aes(df.cohort.b[,1], C.nonlin)) + 
  geom_line(color = "black") +
  geom_point(data = df.cohort.b, aes(df.cohort.b[,1], C.nonlin)) +
  theme_bw() + labs(title = "Nonlinear Cohort Estimates", x = "Birth Year", y = "") +
  scale_x_continuous(limits = c(1899, 2001), breaks = c(seq(1900, 2000, by =20))) +
  theme( plot.title = element_text(size=11), axis.text = element_text(size = 9, color = "black"), axis.title= element_text(size= 11), legend.position = "none") +
  geom_rect(data=C.rect, inherit.aes = F, 
            mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill=color)) +  
  scale_fill_manual(values = alpha(c("c1" = "blue3", "c2" = "blue3", "c3" = "grey" ), 0.15 )) +
  # geom_text(data=C.rect, aes(x=c(1920, 1956, 1990), y=c(45, 45, 45), label = c("C2", "C2", "C1")), color = c("blue", "blue", "grey"), size=4) +
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "black", alpha=0.5) +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) 

ggarrange(A.nolin.g, P.nolin.g, C.nolin.g,  ncol = 3, nrow = 1)

ggsave(path = here("output", "figures"), filename = "nonlinear_apc.pdf", dpi = 1000, device = cairo_pdf, width = 22, height = 7, units = "cm")



# Total Estimates -----------------------------------------------------------
### Assumptions A1 and A2 
# 
# A.total.a12 <- ggplot(df.age.b, aes(df.age.b[,1], age.ub.a12)) + 
#   # geom_point(aes(df.age.b[,1], (age.lb.a12 + age.ub.a12)/2)) +
#   # geom_line(aes(df.age.b[,1], (age.lb.a12 + age.ub.a12)/2)) +
#   scale_x_continuous(limits = c(18, 76), breaks = c(20, 30, 40, 50, 60, 70)) +
#   geom_ribbon(aes(x=df.age.b[,1], ymax=age.lb.a12, ymin=age.ub.a12), fill= "blue", alpha=.2) +
#   theme_bw() + labs(title = "Total Age Estimates", x = "Age", y = "Turnout %") +
#   geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "blue", alpha=0.5) +
#   scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) +
#   theme( plot.title = element_text(size=11),  
#          axis.text = element_text(size = 9, color = "black"),
#          axis.title= element_text(size= 11))  
# # geom_text(x=40, y=65, label="APC-I", color = "red") +
# # geom_text(x=70, y=85, label="Bounding Analysis", color = "blue")
# 
# P.total.a12 <- ggplot(df.period..b, aes(df.period..b[,1], period.ub.a12)) + 
#   scale_x_continuous(limits = c(1976, 2020), breaks = c(seq(1980, 2020, by =10))) +
#   geom_ribbon(aes(x=df.period..b[,1], ymax=period.lb.a12, ymin=period.ub.a12), fill="blue", alpha=.2) +
#   theme_bw() + labs(title = "Total Period Estimates", x = "Election Year", y = "") +
#   geom_hline(yintercept=model$coefficients[1], linetype="dashed", color =  "blue", alpha=0.5) +
#   theme( plot.title = element_text(size=11),  
#          axis.text = element_text(size = 9, color = "black"),
#          axis.title= element_text(size= 11)) +
#   scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) 
# 
# C.total.a12 <- ggplot(df.cohort.b, aes(df.cohort.b[,1], cohort.ub.a12)) + 
#   scale_x_continuous(limits = c(1899, 2000), breaks = c(seq(1900, 2000, by =20))) +
#   geom_ribbon(aes(x=apci.C.est[,1], ymax=apci.C.est[,4], ymin=apci.C.est[,3]), fill= "red", alpha=.5) +
#   geom_ribbon(aes(x=df.cohort.b[,1], ymax=cohort.lb.a12, ymin=cohort.ub.a12), fill="blue", alpha=.2) +
#   theme_bw() + labs(title = "Total Cohort Estimates", x = "Birth Year", y = "") +
#   geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "blue", alpha=0.5) +  
#   theme( plot.title = element_text(size=11),  
#          axis.text = element_text(size = 9, color = "black"),
#          axis.title= element_text(size= 11)) +
#   scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) 
# 
# ggarrange(A.total.a12, P.total.a12, C.total.a12,  ncol = 3, nrow = 1)
# ggsave(path = here("output", "figures"), filename = "total_apc_a12.pdf", dpi = 1000, device = cairo_pdf, width = 22, height = 7, units = "cm")

### Assumptions A1 and A2

A.total.a12 <- ggplot(df.age.b, aes(df.age.b[,1], age.ub.a12)) + 
  geom_point(data = apci.A.est, aes(apci.A.est[,1], apci.A.est[,2]), color = "red") +
  geom_line(data=apci.A.est, aes(apci.A.est[,1], apci.A.est[,2]), color="red") +
  geom_ribbon(data=apci.A.est, aes(x=apci.A.est[,1], ymax=apci.A.est[,4], ymin=apci.A.est[,3]), fill= "red", alpha=.3) +
  scale_x_continuous(limits = c(18, 76), breaks = c(20, 30, 40, 50, 60, 70)) +
  geom_ribbon(aes(x=df.age.b[,1], ymax=age.lb.a12, ymin=age.ub.a12), fill= "blue", alpha=.2) +
  theme_bw() + labs(title = "Total Age Estimates", x = "Age", y = "Turnout %") +
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "blue", alpha=0.5) +
  geom_hline(yintercept=apci.model$intercept, linetype="dashed", color = "red", alpha=0.5) +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) +
  theme( plot.title = element_text(size=11),  
         axis.text = element_text(size = 9, color = "black"),
         axis.title= element_text(size= 11))  


P.total.a12 <- ggplot(df.period..b, aes(df.period..b[,1], period.ub.a12)) + 
  geom_point(data = apci.P.est, aes(apci.P.est[,1], apci.P.est[,2]), color = "red") +
  geom_line(data=apci.P.est, aes(apci.P.est[,1], apci.P.est[,2]), color="red") + 
  geom_ribbon(data= apci.P.est, aes(x=apci.P.est[,1], ymax=apci.P.est[,4], ymin=apci.P.est[,3]), fill= "red", alpha=.3) +
  scale_x_continuous(limits = c(1976, 2020), breaks = c(seq(1980, 2020, by =10))) +
  geom_ribbon(aes(x=df.period..b[,1], ymax=period.lb.a12, ymin=period.ub.a12), fill="blue", alpha=.2) +
  theme_bw() + labs(title = "Total Period Estimates", x = "Election Year", y = "") +
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "blue", alpha=0.5) +
  geom_hline(yintercept=apci.model$intercept, linetype="dashed", color = "red", alpha=0.5) +
  theme( plot.title = element_text(size=11),  
         axis.text = element_text(size = 9, color = "black"),
         axis.title= element_text(size= 11)) +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) 

C.total.a12 <- ggplot(df.cohort.b, aes(df.cohort.b[,1], cohort.ub.a12)) + 
  geom_point(data = apci.C.est, aes(apci.C.est[,1], apci.C.est[,2]), color = "red") +
  geom_line(data=apci.C.est, aes(apci.C.est[,1], apci.C.est[,2]), color="red") + 
  geom_ribbon(data= apci.C.est, aes(x=apci.C.est[,1], ymax=apci.C.est[,4], ymin=apci.C.est[,3]), fill= "red", alpha=.3) +
  scale_x_continuous(limits = c(1899, 2000), breaks = c(seq(1900, 2000, by =20))) +
  geom_ribbon(aes(x=df.cohort.b[,1], ymax=cohort.lb.a12, ymin=cohort.ub.a12), fill="blue", alpha=.2) +
  theme_bw() + labs(title = "Total Cohort Estimates", x = "Birth Year", y = "") +
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "blue", alpha=0.5) +  
  geom_hline(yintercept=apci.model$intercept, linetype="dashed", color = "red", alpha=0.5) +
  theme( plot.title = element_text(size=11),  
         axis.text = element_text(size = 9, color = "black"),
         axis.title= element_text(size= 11)) +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) 

ggarrange(A.total.a12, P.total.a12, C.total.a12,  ncol = 3, nrow = 1)
ggsave(path = here("output", "figures"), filename = "total_apc_a12_with_api.pdf", dpi = 1000, device = cairo_pdf, width = 22, height = 7, units = "cm")


### Assumption A1, A2, and C1
A.total.a12.c1 <- ggplot(df.age.b, aes(df.age.b[,1], age.ub.a12.c1)) + 
  geom_point(data = apci.A.est, aes(apci.A.est[,1], apci.A.est[,2]), color = "red") +
  geom_line(data=apci.A.est, aes(apci.A.est[,1], apci.A.est[,2]), color="red") +
  geom_ribbon(data=apci.A.est, aes(x=apci.A.est[,1], ymax=apci.A.est[,4], ymin=apci.A.est[,3]), fill= "red", alpha=.3) +
  geom_ribbon(aes(x=df.age.b[,1], ymax=age.lb.a12.c1, ymin=age.ub.a12.c1), fill= "blue", alpha=.2) +
  scale_x_continuous(limits = c(18, 76), breaks = c(20, 30, 40, 50, 60, 70)) +
  theme_bw() + labs(title = "Total Age Estimates", x = "Age", y = "Turnout %") +
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "blue", alpha=0.5) +
  geom_hline(yintercept=apci.model$intercept, linetype="dashed", color = "red", alpha=0.5) +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90))  +
  theme( plot.title = element_text(size=11),  
         axis.text = element_text(size = 9, color = "black"),
         axis.title= element_text(size= 11)) 

P.total.a12.c1 <- ggplot(df.period..b, aes(df.period..b[,1], period.ub.a12.c1)) + 
  geom_point(data = apci.P.est, aes(apci.P.est[,1], apci.P.est[,2]), color = "red") +
  geom_line(data=apci.P.est, aes(apci.P.est[,1], apci.P.est[,2]), color="red") +
  geom_ribbon(data= apci.P.est, aes(x=apci.P.est[,1], ymax=apci.P.est[,4], ymin=apci.P.est[,3]), fill= "red", alpha=.3) +
  geom_ribbon(aes(x=df.period..b[,1], ymax=period.lb.a12.c1, ymin=period.ub.a12.c1), fill="blue", alpha=.2) +
  scale_x_continuous(limits = c(1976, 2020), breaks = c(seq(1980, 2020, by =10))) +
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "blue", alpha=0.5) +   
  theme_bw() + labs(title = "Total Period Estimates", x = "Election Year", y = "") +
  geom_hline(yintercept=apci.model$intercept, linetype="dashed", color = "red", alpha=0.5) +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) +
  theme( plot.title = element_text(size=11),  
         axis.text = element_text(size = 9, color = "black"),
         axis.title= element_text(size= 11)) 

C.total.a12.c1 <- ggplot(df.cohort.b, aes(df.cohort.b[,1], cohort.ub.a12.c1)) + 
  geom_point(data = apci.C.est, aes(apci.C.est[,1], apci.C.est[,2]), color = "red") +
  geom_line(data=apci.C.est, aes(apci.C.est[,1], apci.C.est[,2]), color="red") +
  geom_ribbon(aes(x=df.cohort.b[,1], ymax=cohort.lb.a12.c1, ymin=cohort.ub.a12.c1), fill="blue", alpha=.2) +
  geom_ribbon(data= apci.C.est, aes(x=apci.C.est[,1], ymax=apci.C.est[,4], ymin=apci.C.est[,3]), fill= "red", alpha=.3) +
  theme_bw() + labs(title = "Total Cohort Estimates", x = "Birth Year", y = "") +
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "blue", alpha=0.5) +   
  scale_x_continuous(limits = c(1899, 2000), breaks = c(seq(1900, 2000, by =20))) +
  geom_hline(yintercept=apci.model$intercept, linetype="dashed", color = "red", alpha=0.5) +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) +
  theme( plot.title = element_text(size=11),  
         axis.text = element_text(size = 9, color = "black"),
         axis.title= element_text(size= 11)) 

ggarrange(A.total.a12.c1, P.total.a12.c1, C.total.a12.c1,  ncol = 3, nrow = 1)
ggsave(path = here("output", "figures"), filename = "total_apc_a12_c1_with_api.pdf", dpi = 1000, 
       device = cairo_pdf, width = 22, height = 7, units = "cm")


### Assumption A1, A2, and C2
A.total.a12.c2 <- ggplot(df.age.b, aes(df.age.b[,1], age.ub.a12.c2)) + 
  geom_point(data = apci.A.est, aes(apci.A.est[,1], apci.A.est[,2]), color = "red") +
  geom_line(data=apci.A.est, aes(apci.A.est[,1], apci.A.est[,2]), color="red") +
  geom_ribbon(data=apci.A.est, aes(x=apci.A.est[,1], ymax=apci.A.est[,4], ymin=apci.A.est[,3]), fill= "red", alpha=.3) +
  geom_ribbon(aes(x=df.age.b[,1], ymax=age.lb.a12.c2, ymin=age.ub.a12.c2), fill= "blue", alpha=.2) +
  scale_x_continuous(limits = c(18, 76), breaks = c(20, 30, 40, 50, 60, 70)) +
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "blue", alpha=0.5) + 
  theme_bw() + labs(title = "Total Age Estimates", x = "Age", y = "Turnout %") +
  geom_hline(yintercept=apci.model$intercept, linetype="dashed", color = "red", alpha=0.5) +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) +
  theme( plot.title = element_text(size=11),  
         axis.text = element_text(size = 9, color = "black"),
         axis.title= element_text(size= 11)) 

P.total.a12.c2 <- ggplot(df.period..b, aes(df.period..b[,1], period.ub.a12.c2)) + 
  geom_point(data = apci.P.est, aes(apci.P.est[,1], apci.P.est[,2]), color = "red") +
  geom_line(data=apci.P.est, aes(apci.P.est[,1], apci.P.est[,2]), color="red") +  
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "blue", alpha=0.5) + 
  geom_ribbon(data= apci.P.est, aes(x=apci.P.est[,1], ymax=apci.P.est[,4], ymin=apci.P.est[,3]), fill= "red", alpha=.3) +
  geom_ribbon(aes(x=df.period..b[,1], ymax=period.lb.a12.c2, ymin=period.ub.a12.c2), fill="blue", alpha=.2) +
  scale_x_continuous(limits = c(1976, 2020), breaks = c(seq(1980, 2020, by =10))) +
  theme_bw() + labs(title = "Total Period Estimates", x = "Election Year", y = "") +
  geom_hline(yintercept=apci.model$intercept, linetype="dashed", color = "red", alpha=0.5) +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) +
  theme( plot.title = element_text(size=11),  
         axis.text = element_text(size = 9, color = "black"),
         axis.title= element_text(size= 11)) 

C.total.a12.c2 <- ggplot(df.cohort.b, aes(df.cohort.b[,1], cohort.ub.a12.c2)) + 
  geom_point(data = apci.C.est, aes(apci.C.est[,1], apci.C.est[,2]), color = "red") +
  geom_line(data=apci.C.est, aes(apci.C.est[,1], apci.C.est[,2]), color="red") +
  geom_hline(yintercept=model$coefficients[1], linetype="dashed", color = "blue", alpha=0.5) + 
  geom_ribbon(aes(x=df.cohort.b[,1], ymax=cohort.lb.a12.c2, ymin=cohort.ub.a12.c2), fill="blue", alpha=.2) +
  geom_ribbon(data= apci.C.est, aes(x=apci.C.est[,1], ymax=apci.C.est[,4], ymin=apci.C.est[,3]), fill= "red", alpha=.3) +
  theme_bw() + labs(title = "Total Cohort Estimates", x = "Birth  Year", y = "") +
  scale_y_continuous(breaks=seq(30, 90, 10), limits=c(30, 90)) +
  scale_x_continuous(limits = c(1899, 2000), breaks = c(seq(1900, 2000, by =20))) +
  geom_hline(yintercept=apci.model$intercept, linetype="dashed", color = "red", alpha=0.5) +
  theme( plot.title = element_text(size=11),  
         axis.text = element_text(size = 9, color = "black"),
         axis.title= element_text(size= 11)) 

ggarrange(A.total.a12.c2, P.total.a12.c2, C.total.a12.c2,  ncol = 3, nrow = 1)
ggsave(path = here("output", "figures"), filename = "total_apc_a12_c2_with_apci.pdf", 
       dpi = 1000, device = cairo_pdf, width = 22, height = 7, units = "cm")



# APC-I Model Intra-cohort deviations
C.slope <- ggplot(data = apci.C.slope, mapping = aes(x = apci.C.slope[,1], y = apci.C.slope[,2], ymin = apci.C.slope[,3], ymax = apci.C.slope[,4])) +
  geom_pointrange(size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1, alpha = 0.3) +
  scale_x_continuous(limits = c(1902, 1996), breaks = seq(1900, 1996, 10)) +
  scale_y_continuous(breaks = seq(-15, 15, 5), limits = c(-15, 15)) +
  theme_bw() + labs(title = "APC-I: Intra-Cohort Variation Estimates", x = "Birth Year", y = "") +
  theme(plot.title = element_text(size = 11),  
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 11))

ggsave(path = here("output", "figures"), filename = "intra_cohort_apci.pdf", dpi = 1000, device = cairo_pdf, width = 22, height = 7, units = "cm")


#=======================================
# Print Session Info
#=======================================

cat("\n\n---- Session Info ----\n")
print(sessionInfo())

### END OF CODE ###
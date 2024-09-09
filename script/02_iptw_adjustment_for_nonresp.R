###############################################################################
##                                                                           ##
## 02: IPTW Adjustment for Non-Response in CPS                               ##
##                                                                           ##  
## Inputs:                                                                   ##
##   1) data/cps.76.20.recoded.rds                                           ##
##                                                                           ##
## Notes:                                                                    ##
##   - Models election-year-specific non-response for the turnout question   ##
##   - Extracts missing probabilities from the models to construct           ##
##     missingness-corrected weights                                         ##
##                                                                           ##
## Final Data Output:                                                        ##
##   1) cps.76.20.recoded.wt.all.rds (All respondents)                       ##
##   2) cps.76.20.recoded.wt.voted.rds (Restricted to non-missing Rs)        ##
##                                                                           ##
###############################################################################



#=======================================
# Preliminaries
#=======================================


##########
##########    LIBRARIES
##########

# packages
packages <- c("tidyverse", "survey", "here")

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

# library(cobalt) # to check covariate balance after weighting


#=======================================
# Call-in data and set-up for analysis
#=======================================


# call-in data ------------------------------------------------------------
here("data", "cps.76.20.recoded.rds") %>% readRDS() -> cps.76.20
arrange(cps.76.20, P, id) # sort data by year and id
  
# Variable availability by year -------------------------------------------
table(cps.76.20$metro, cps.76.20$P, useNA = "ifany")
table(cps.76.20$faminc.f, cps.76.20$P, useNA = "ifany")
table(cps.76.20$inttype, cps.76.20$P, useNA = "ifany")
table(cps.76.20$race, cps.76.20$P, useNA = "ifany")
table(cps.76.20$racehisp, cps.76.20$P, useNA = "ifany")
table(cps.76.20$hispanic, cps.76.20$P, useNA = "ifany")
table(cps.76.20$race, cps.76.20$hispanic, useNA = "ifany")
table(cps.76.20$statefips, cps.76.20$P, useNA = "ifany") # state not comparable for year 1976
table(cps.76.20$sex, cps.76.20$P, useNA = "ifany")
table(cps.76.20$marst, cps.76.20$P, useNA = "ifany")
table(cps.76.20$famsize, cps.76.20$P, useNA = "ifany")
table(cps.76.20$nchild, cps.76.20$P, useNA = "ifany")
table(cps.76.20$empstat, cps.76.20$P, useNA = "ifany")
table(cps.76.20$occ.two.digit, useNA = "ifany")
table(cps.76.20$occ.two.digit, cps.76.20$P, useNA = "ifany")
table(cps.76.20$citizen, cps.76.20$P, useNA = "ifany") # only from 1996 onwards
table(cps.76.20$educ, cps.76.20$P, useNA = "ifany") 
table(cps.76.20$A, cps.76.20$P, useNA = "ifany") 
table(cps.76.20$voted, cps.76.20$P, useNA = "ifany")
table(cps.76.20$voted, cps.76.20$citizen, useNA = "ifany")
# table(cps.76.20$voteres, cps.76.20$P, useNA = "ifany") # not available for 1984 and 1988 
# table(cps.76.20$voteresp, cps.76.20$P, useNA = "ifany") 
# table(cps.76.20$voteres, cps.76.20$P, useNA = "ifany") 
# table(cps.76.20$voteres, cps.76.20$votemiss, useNA = "ifany")
# table(cps.76.20$bpl, cps.76.20$year, useNA = "ifany") # only from 1996 onwards

# model specifications ----------------------------------------------------

y <- "votemiss"

x.allyears <- c(
  "inttype", "racehisp", "statefips", "sex", "marst", "famsize", 
  "nchild", "occ.two.digit", "educ", "age", "I(age^2)")

x.1976 <- c()
x.1980 <- c("metro")
x.1984 <- c("metro",  "faminc.f")
x.1988 <- c("metro",  "faminc.f")
x.1992 <- c("metro",  "faminc.f") 
x.1996 <- c("metro",  "faminc.f", "citizen")
x.2000 <- c("metro",  "faminc.f", "citizen")
x.2004 <- c("metro", "faminc.f",  "citizen")
x.2008 <- c("metro", "faminc.f",  "citizen")
x.2012 <- c("metro",  "faminc.f",  "citizen")
x.2016 <- c("metro", "faminc.f",  "citizen")
x.2020 <- c("metro",  "faminc.f",  "citizen")


# loop to generate period-specific specification formulas

P.list <- seq(1976, 2020, by = 4)

for(i in P.list) {
  X <- append(x.allyears, get(paste("x.", i, sep = "")))
  assign(paste("model.", i, sep = ""), as.formula(paste(y, paste(X, collapse = " + "), sep = " ~ ")))
}


# estimate logit models for missing values --------------------------------

for(i in P.list) {
  
  cat('\n', "Missing weight estimation for year", i, '\n')
  
  survey.design <- svydesign(id = ~1, weights = ~wt, data =  cps.76.20[cps.76.20$P == i,])   # declare survey design
  logit.model <- paste("est.logit", i, sep = ".") ; pr.miss <- paste("pr.miss", i, sep = ".") # common model and variable labels
  
  assign(logit.model, svyglm(get(paste("model", i, sep = ".")),  design = survey.design, family = "quasibinomial")) # logit model for missingness
  cat('\n', "Logit model for", i, '\n')
  print(summary(get(logit.model))) # check the estimates
  
  cat('\n', "Checking for no perfect prediction problems........", '\n')
  print(ifelse(
      nrow(subset(cps.76.20, P == i)) ==  length(fitted(get(logit.model))), 
        "All estimated", "Some excluded" )) # check that all observations have been estimated

  assign(pr.miss, predict(get(logit.model), type = "response")) # predicted probabilities
  cat('\n', "Summary of the probabilities before trimming", '\n')
  print(summary(get(pr.miss))) # check the predicted probabilities
  
  lb <-  quantile(get(pr.miss), 0.01); ub <-  quantile(get(pr.miss), 0.99) # upper and lower bounds
  unname(lb) ; unname(ub) # check the upper and lower bounds
  assign(pr.miss, ifelse(get(pr.miss) < lb & !is.na(get(pr.miss)), lb,
                         ifelse(get(pr.miss) > ub & !is.na(get(pr.miss)), ub, get(pr.miss)))) # trim extreme weights (to upper and lower 0.005 percentiles)
  cat('\n', "Summary of the probabilities after trimming", '\n')
  print(summary(get(pr.miss))) # check the trimmed predicted probabilities
  
  print(warnings()) # view warnings
}


# create IPTW for missingness ---------------------------------------------

# combine all year-specific estimates into a single column
pr.miss <- c(as.vector(pr.miss.1976), as.vector(pr.miss.1980), as.vector(pr.miss.1984), as.vector(pr.miss.1988), 
             as.vector(pr.miss.1992), as.vector(pr.miss.1996), as.vector(pr.miss.2000), as.vector(pr.miss.2004), 
             as.vector(pr.miss.2008), as.vector(pr.miss.2012), as.vector(pr.miss.2016), as.vector(pr.miss.2020))

# estimate the iptws (it is noly generated for those who are not missing)
cps.76.20$wt.miss <- ifelse(cps.76.20$votemiss == 1, NA, cps.76.20$wt/(1-pr.miss))

# check that they were implemented properly
cps.76.20 %>% 
  group_by(votemiss) %>% 
  summarize(
    n = n(),
    wt_mean = mean(wt.miss, na.rm = FALSE),
    wt_miss = sum(is.na(wt.miss)))

# check the balance of variables after weighting --------------------------

# covars_list <- as.formula("age+faminc+educ+race+hispanic+marst+inttype+nchild+famsize")

for(i in P.list) {
 
  cat('\n', "Covariate Balance pre and post-weighting for year", i, '\n')
  
  ds.base.all <- svydesign(
    id = ~1, weights = ~wt, data = subset(cps.76.20, P == i))   # declare survey design
  ds.base.nomiss <- svydesign(
    id = ~1, weights = ~wt, data = subset(cps.76.20, P == i & votemiss == 0))  # declare survey design
  ds.iptw.nomiss <- svydesign(
    id = ~1, weights = ~wt.miss, data = subset(cps.76.20, P == i & votemiss == 0))   # declare survey design
  
  base.all <- svymean(~age+faminc.f+educ+race+hispanic+marst+inttype+nchild+famsize, ds.base.all)
  base.nonmiss <- svymean(~age+faminc.f+educ+race+hispanic+marst+inttype+nchild+famsize, ds.base.nomiss)
  ipwt <- svymean(~age+faminc.f+educ+race+hispanic+marst+inttype+nchild+famsize, ds.iptw.nomiss)
  
  cat('\n', "Compare for year", i, '\n')
  print(round(cbind(base.all, base.nonmiss, ipwt), 3))
   
}


# Re-scale by relative N of each year --------------------------------------

message('Re-scaling the weights ...')

# re-scale
wt.miss.rs <- cps.76.20 %>% 
  filter(votemiss == 0) %>% # subset to non-missing
  group_by(P) %>% 
  mutate(wt.miss.sum = sum(wt.miss)) %>% 
  ungroup() %>% 
  mutate(wt.miss.rs = (n() / n_distinct(.$P)) * (wt.miss / wt.miss.sum) ) %>% 
  select(P, id, wt.miss.rs)

cps.76.20 <- full_join(cps.76.20, wt.miss.rs, by = c("P", "id"))

# check that the weights make sense
cps.76.20 %>% 
  filter(votemiss == 0) %>% 
  group_by(P) %>% 
  summarize(mean(wt.miss.rs)) 

cor(cbind(cps.76.20$wt, cps.76.20$wt.miss, cps.76.20$wt.miss.rs ), use = "pairwise.complete.obs")


#=======================================
# Save Data
#=======================================

message('Saving Data ...')

# save all respondents
saveRDS(
  cps.76.20,
  here("data", "cps.76.20.recoded.wt.all.rds")
)

# save only non-missing respondents
saveRDS(
  cps.76.20[cps.76.20$votemiss == 0, ],
  here("data", "cps.76.20.recoded.wt.voted.rds")
)

message('Done!\n')


#=======================================
# Print Session Info
#=======================================

cat("\n\n---- Session Info ----\n")
print(sessionInfo())

### END OF CODE ###
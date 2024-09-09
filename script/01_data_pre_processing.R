###############################################################################
##                                                                           ##
## 01: Data Preprocessing                                                    ##
##                                                                           ##  
## Data Used:                                                                ##
##   1) 1976-2020 CPS November Voter and Registration Supplement             ##
##                                                                           ##
## Raw Data Files:                                                           ##
##   rawdata/cps_vrs_1976_2020.csv                                           ##
##   (Downloadable from IPUMS-CPS)                                           ##
##                                                                           ##
## Notes:                                                                    ##
##   - Imports the 1976-2020 CPS-VRS data                                    ##
##   - Recodes variables for analysis                                        ##
##   - Retains only relevant variables                                       ##
##                                                                           ##
## Final Data Output:                                                        ##
##   1) data/cps.76.20.recoded.rds                                           ##
##                                                                           ##
###############################################################################


#=======================================
# Preliminaries
#=======================================

##########
##########    LIBRARIES
##########



# packages
packages <- c("tidyverse", "data.table", "here")

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

# library(sjlabelled)

# clear console
rm(list = ls())


# List of variables to keep eventually
vars.to.keep <- c(
  "hhid", "id",  "age", "cohort", "A", "P", "C", "metro", "faminc", "faminc.f", "inttype", 
  "race", "statefips", "sex","marst", "famsize", "nchild",  "bpl", "citizen", 
  "hispanic",  "empstat", "occ2010", "occ.emp.comb", "occ.two.digit", "classwkr",  "educ", 
  "diffany",  "voteres", "voteresp", "voted",  "votemiss", "wt", "racehisp", "battlegr.all", "battlegr8", 
  "region", "south")

message('\n\nPreparing CPS 1976-2020 Data --------------------')


#=======================================
# Call-in Data
#=======================================

# call-in data only for those aged between 18 and 79 and born after 1900
message('Loading Cumulative Datafile ...')

# Analytical N construction
# data <- fread(here("rawdata", "cps_vrs_1976_2020.csv")) 
# data <- data %>%  filter(YEAR %in% seq(1976, 2020, by = 4) & AGE >= 18) # 
# data %>% filter(AGE>=18 & CITIZEN %in% c(1,2,3,4,NA)) %>% summarise(n = n()) #  1,160,569 respondents of 18 who should be in universe
# data %>% filter(AGE>=18 & CITIZEN %in% c(1,2,3,4,NA)) %>%  mutate(niu = ifelse(VOTED  ==99,1, 0 )) %>% group_by(niu) %>% summarise(n = n()) #  23579 additional R NIU for unclear reasons hence 1,136,990 IU
# data %>% mutate(niu = ifelse(VOTED  ==99,1, 0 )) %>% 
#   filter(AGE>=18 & CITIZEN %in% c(1,2,3,4, NA) & niu ==0) %>% mutate(nonresp = ifelse(VOTED %in% c(96, 97, 98), 1,0)) %>% group_by(nonresp) %>% summarise(n = n()) # 99163 further dropped for non-response
# data %>% mutate(niu = ifelse(VOTED  ==99,1, 0 ),
#                 nonresp = ifelse(VOTED %in% c(96, 97, 98), 1, 0)) %>% 
#   filter(AGE>=18 & CITIZEN %in% c(1,2,3,4, NA) & niu ==0 & nonresp ==0) %>% mutate(age_lim = ifelse(AGE >= 78, 1, 0)) %>% 
#   group_by(age_lim) %>% summarise(n = n()) # 55,969 further dropped for age 78 or above
# data %>% mutate(niu = ifelse(VOTED  ==99,1, 0),
#                 nonresp = ifelse(VOTED %in% c(96, 97, 98), 1, 0), 
#                 age_lim = ifelse(AGE >= 78, 1, 0)) %>% 
#   filter(AGE>=18 & CITIZEN %in% c(1,2,3,4, NA) & niu ==0 & nonresp ==0 & age_lim == 0) %>% summarise(n = n()) # final N 981,858


data <- fread(here("rawdata", "cps_vrs_1976_2020.csv") ) %>%
  mutate(cohort = YEAR - AGE) %>%  
  filter(
    AGE>=18 & AGE <= 77,
    YEAR %in% seq(1976, 2020, by = 4), # include ages up to 77 given that some years topcode age at 80
    cohort <= 2002,
    cohort >= 1899, # exclude the 2001 cohort to get consistency
    CITIZEN %in% c(1,2,3,4, NA)) %>% # exclude non-citizens from 1996 onwards
  tibble()


#=======================================
# Recoding Variables
#=======================================

### List of Variables to be used and their availability

# age: respondent age (1976-2020)
# year: Survey year (1976-2020)
# metro: Metropolitan central city status (1980-2020)
# faminc: Family income of householder (1982-2020)
# inttype: Interview type (1976-2020)
# race: Race of respondent (1976-2020)
# statefips: State FIPS code (1976-2020)
# sex: Respondent sex (1976-2020)
# marst: Marital status (1976-2020)
# famsize: Number of own family members in hh (1976-2020) 
# nchild: Number of own children in household (1976-2020)
# bpl: Birthplace (1994-2020)
# citizen: Citizenship status (1994-2020)
# hispan: Hispanic origin (1976-2020)
# empstat: Employment status (1976-2020)
# occ2010: Occupation based on 2010 codes (1976-2020)
# voteres: Duration of residence at current address (1976-2020)
# voteresp: Self or proxy respondent for voter supplement (1976-2020)

message('Recode and rename variables ...')

data <- data %>%
  mutate(
    age = AGE,
    cohort = YEAR - AGE,
    A = cut_width(AGE, 4, closed = "left"),
    # C = as.factor(cut_width(cohort, 4, closed = "left")),
    C = cut(cohort, breaks = c(seq(1899, 2003, by =4)), right = F),
    P = as.factor(YEAR),
    metro = case_when(
      METRO == 0 ~ "Not identifiable", 
      METRO == 1 ~ "Non-metro", 
      METRO %in% c(2,3,4) ~ "Metro area",
      METRO == 9 ~ "Missing/unknown"),
    faminc = recode(
      FAMINC,
      `100` = 2500L, `110` = 500L, `111` = 250L, `112` = 750L, 
      `120` = 1500L, `121` = 1250L, `122` = 1750L, 
      `130` = 2500L,`131` = 2250L, `132` = 2750L, 
      `140` = 3500L, `141` = 3250L, `142` = 3750L,
      `150` = 4500L,
      `200` = 6500L, `210` = 6250L, `220` = 5500L, `230` = 7000L, 
      `231` = 6750L, `232` = 6500L, `233` = 7250L, `234` = 7500L,
      `300` = 8750L, `310` = 7750L, `320` = 8250L, `330` = 8750L,
      `340` = 8500L, `350` = 9500L,
      `400` = 12500L, `410` = 10500L, `420` = 11500L, `430` = 11250L,
      `440` = 11000L, `450` = 12500L, `460` = 13500L, `480` = 13500L,
      `490` = 14500L,
      `500` = 17500L, `510` = 15500L, `520` = 16500L, `530` = 17500L,
      `540` = 16250L, `550` = 18250L, `560` = 19000L,
      `600` = 22500L,
      `700` = 37500L, `710` = 27500L, `720` = 32500L, `730` = 37500L,
      `740` = 45000L, 
      `800` = 75000L, `810` = 62500L, `820` = 55000L, `830` = 67500L,
      `840` = 112500L, `841` = 87500L, `842` = 125000L, `843` = 225000L,
      `995` = 99999999L, `996` = 99999999L,
      `997` = 99999999L, `999` = 99999999L
      ),
    faminc.f = as.factor(faminc),
    inttype = recode(
      INTTYPE,
      `0` = "Non-interview",`1` = "In person",  `2` = "Phone",`3` = "Other", `4` = "Other"
    ),
    race = case_when(
      RACE == 100 ~ "White",
      RACE %in% c(200, 801, 805, 806, 807, 810, 811, 814, 816, 818,820, 830) ~ "Black",
      RACE %in% c(300, 802) ~ "Native Indian",
      RACE %in% c(650, 651, 652, 802, 803, 804, 808, 809, 812, 813, 815, 817, 819) ~ "AAPI",
      RACE == 700 ~ "Other",
      RACE == 999 ~ "Missing"
    ),
    statefips = as.factor(STATEFIP),
    battlegr.all = as.factor(ifelse(STATEFIP %in% c(4, 12, 23, 26, 27, 31, 32, 33, 37, 39, 42, 48, 51, 55), 1 , 0)),
    battlegr8 = as.factor(ifelse(STATEFIP %in% c(4, 12, 13, 26, 27, 37, 42, 55), 1 , 0)),
    sex = as.factor(SEX),
    marst = as.factor(ifelse(MARST == 9, as.character(NA), as.factor(MARST))),
    famsize = ifelse(FAMSIZE < 6, FAMSIZE, 6),
    nchild = ifelse(NCHILD < 4, NCHILD, 4),    
    bpl = as.factor(ifelse(BPL == 9900, 1, ifelse(!is.na(BPL), 0, NA))),
    citizen = as.factor(CITIZEN),
    hispanic = ifelse(HISPAN == 0, 0, ifelse(HISPAN !=901 & HISPAN != 902, 1, "Missing")),
    empstat = case_when(
      EMPSTAT == 0 ~ "Missing",
      EMPSTAT %in% c(1, 10, 12) ~ "Employed",
      EMPSTAT %in% c(20, 21, 22) ~ "Unemployed",
      EMPSTAT %in% c(30, 31, 32, 33, 34, 35, 36) ~ "Not in Labor Force"      
    ),
    occ2010 = ifelse(OCC2010 != 9999, as.factor(OCC2010), as.character(NA)),
    occ.emp.comb = ifelse(!is.na(occ2010), occ2010, ifelse(empstat != "Employed", "NIL or unemployed", "Missing")),
    occ.two.digit = case_when(
      OCC2010 >= 10 & OCC2010 <= 430 ~ "Managment in Business, Science, Arts",
      OCC2010 >= 500 & OCC2010 <= 730 ~ "Business Operations Specialists",
      OCC2010 >= 800 & OCC2010 <= 950 ~ "Financial Spcialists",
      OCC2010 >= 1000 & OCC2010 <= 1240 ~ "Computer and Mathematical",
      OCC2010 >= 1300 & OCC2010 <= 1540 ~ "Architecture and Engineering",
      OCC2010 >= 1550 & OCC2010 <= 1560 ~ "Technicians",
      OCC2010 >= 1600 & OCC2010 <= 1980 ~ "Life, Physical, and Social Science",
      OCC2010 >= 2000 & OCC2010 <= 2060 ~ "Community and Social Services ",
      OCC2010 >= 2100 & OCC2010 <= 2150 ~ "Legal",
      OCC2010 >= 2200 & OCC2010 <= 2550 ~ "Education, Training, and Library",
      OCC2010 >= 2600 & OCC2010 <= 2920 ~ "Arts, Design, Entertainment, Sports, and Media",
      OCC2010 >= 3000 & OCC2010 <= 3540 ~ "Healthcare Practitioners and Technicians",
      OCC2010 >= 3600 & OCC2010 <= 3650 ~ "Healthcare Support",
      OCC2010 >= 3700 & OCC2010 <= 3950 ~ "Protective Service",
      OCC2010 >= 4000 & OCC2010 <= 4150 ~ "Food Preparation and Serving",
      OCC2010 >= 4200 & OCC2010 <= 4250 ~ "Building and Grounds Cleaning and Maintenance",
      OCC2010 >= 4300 & OCC2010 <= 4650 ~ "Personal Care and Service",
      OCC2010 >= 4700 & OCC2010 <= 4965 ~ "Sales and Related",
      OCC2010 >= 5000 & OCC2010 <= 5940 ~ "Office and Administrative Support",
      OCC2010 >= 6005 & OCC2010 <= 6130 ~ "Farming, Fisheries, and Forestry",
      OCC2010 >= 6200 & OCC2010 <= 6765 ~ "Construction",
      OCC2010 >= 6800 & OCC2010 <= 6940 ~ "Extraction",
      OCC2010 >= 7000 & OCC2010 <= 7630 ~ "Installation, Maintenance, and Repair",
      OCC2010 >= 7700 & OCC2010 <= 8965 ~ "Production",
      OCC2010 >= 9000 & OCC2010 <= 9830 ~ "Transportation and Material Moving and Military",
      OCC2010 == 9999 & empstat == "Unemployed"  ~ "Unemployed",
      OCC2010 == 9999 & empstat %in% c("Employed", "Missing", "Not in Labor Force") ~ "Not in Labor Force or Missing"
    ),
    classwkr = case_when(
      CLASSWKR == 0 ~ as.character(NA),
      CLASSWKR %in% c(10, 13, 14) ~ "Self-employed",     
      CLASSWKR >= 20 & CLASSWKR <= 28 ~ "Works for salary", 
      CLASSWKR == 29 ~ "Unpaid family work",
      CLASSWKR == 99 ~ as.character(NA)
    ),
    educ = case_when(
      EDUC <= 73 ~ "HS or lss",
      EDUC >= 80 & EDUC <= 100 ~ "Some college",
      EDUC >= 110 & EDUC <= 122 ~ "BA",
      EDUC >= 123 & EDUC <= 900 ~ "Graduate",
      EDUC == 999 ~ as.character(NA)
    ),
    diffany = case_when(
      DIFFANY == 0 ~ "No",
      DIFFANY == 1 ~ "Yes",
      DIFFANY == NA ~ as.character(NA)
    ),
    voteres = case_when(
      VOTERES >= 10 & VOTERES <= 13 ~ "One year or less",
      VOTERES == 20 ~ "One to two years",
      VOTERES >= 30 & VOTERES <= 32 ~ "Three to five years",
      VOTERES >= 33 & VOTERES <= 900 ~ "More than five years",
      VOTERES %in% c(996, 997, 998, 999) ~ "Missing"
    ),
    voteresp = case_when(
      VOTERESP == 1 ~ "Self",
      VOTERESP == 2 ~ "Proxy",
      VOTERESP == 9 ~ "Not in universe"
    ),
    voted = ifelse(VOTED == 1, 0, ifelse(VOTED == 2, 1, NA)),  
    votemiss = case_when(
      VOTED %in% c(1,2) ~ 0,
      VOTED %in% c(96, 97, 98, 99) ~ 1
    ),
    wt = VOSUPPWT
  ) %>%
  mutate(
  region = case_when(
    STATEFIP %in% c(9, 23, 25, 33, 34, 36, 42, 44, 50) ~ "Northeast", 
    STATEFIP %in% c(17, 18, 19, 20, 26, 27, 29, 31, 38, 39, 46, 55) ~ "Midwest",
    STATEFIP %in% c(1, 5, 10, 11, 12, 13, 21, 22, 24, 28, 37, 40, 45, 47, 48, 51, 54)~ "South", 
    STATEFIP %in% c(2, 4, 6, 8, 16, 35, 30, 49, 32, 56, 15, 41, 53) ~ "West"
  ),
  south = case_when(
    STATEFIP %in% c(1, 5, 12, 13, 21, 22, 28, 37, 45, 47, 48, 51) ~ "Southern",
    TRUE ~ "non-Southern"
  ),
  racehisp = 
    case_when(
      race == "White" & hispanic != "1" ~ "Non-hispanic White",
      race == "White" & hispanic == "1" ~  "Hispanic",
      race == "Black" ~ "Black",
      race %in% c("AAPI", "Multiracial", "Native Indian", "Other")  & hispanic == "1" ~ "Hispanic",
      race == "AAPI" & hispanic %in% c("0", "Missing") ~ "AAPI",  
      race == "Multiracial" & hispanic %in% c("0", "Missing") ~ "Multiracial",
      race == "Native Indian" & hispanic %in% c("0", "Missing") ~ "Native Indian",
      race == "Other"  & hispanic %in% c("0", "Missing") ~ "Other")) %>% 
  arrange(P, CPSID, PERNUM) %>%  # create a simpler unique identifier 
  group_by(P, CPSID) %>% 
  mutate(hhid = cur_group_id(),
         id = hhid * 100 + PERNUM)

# table(data$P, data$hispanic, useNA = "ifany")
# table(data$racehisp, data$race , useNA = "ifany")
# table(data$racehisp, data$hispanic , useNA = "ifany")
# table(data$YEAR, data$RACE, useNA = "ifany")
# table(data$YEAR, data$racehisp, useNA = "ifany")

# check whether the individual id is unique
ifelse(n_distinct(data$hhid) == n_distinct(data$CPSID), "Household Ids are unique", "Household Ids are not unique") 
ifelse(n_distinct(data$id) == nrow(data), "Ids are unique", "Ids are not unique") # check whether the individual id is unique


#=======================================
# Assign value labels
#=======================================

# data_r %>% 
#   set_label(, labels = c("White", "Black", "Native Indian",))

# re-define APC labels
data$A <- factor(paste0(substr(data$A, 2 ,3),"-", as.integer(substr(data$A, 5,6)) - 1))
data$C <- factor(paste0(substr(data$C, 2 ,5),"-", as.integer(substr(data$C, 7,10)) - 1))


#=======================================
# Save the datset
#=======================================

message('Saving Data ...')

saveRDS(
  data[vars.to.keep],
  here("data",
    'cps.76.20.recoded.rds'
  )
)

message('Done!\n')


#=======================================
# Print Session Info
#=======================================

cat("\n\n---- Session Info ----\n")
print(sessionInfo())


### END OF CODE ###



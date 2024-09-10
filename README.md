# Replication Materials for "Age, Period, and Cohort Analysis with Bounding and Interactions."

## Overview

This repository provides replication materials for:

* Lee, Jiwon. Forthcoming. "Age, Period, and Cohort Analysis with Bounding and Interactions." _Sociological Methods & Research_


## Roadmap of Code

The content below is copied from the file `roadmap-for-apc-turnout.txt`, which can be found in the `docs` folder of the repository.

### Notes

1. The following data files (.csv) must be downloaded from [HERE](https://www.dropbox.com/scl/fo/21plc08xo0zc59nypuwdf/APVZdluFNZiFq8ASfmsQ4j0?rlkey=uni22opsbafz3fhlw2qg77tec&st=mbcnviuk&dl=0) and placed in the `/rawdata` folder before executing any R code scripts:

   - `cps_vrs_1976_2020.csv` (1976-2020 Current Population Survey Voter and Registration Supplements; CPS-VRS)
   - `mcdonald_vep_turnout_1789_present.csv` (National Election Turnout Rates, 1789-2020)

The 1976-2020 CPS-VRS datasets are provided by the Integrated Public Use Microdata Series ([IPUMS](https://cps.ipums.org/cps/)). These datasets include the following variables (IPUMS labels), though only a subset is ultimately used in the analysis:

"YEAR", "SERIAL", "MONTH", "HWTFINL", "CPSID", "STATEFIP", "METRO", "METAREA", "COUNTY", "STATECENSUS", "CBSASZ", "MSACMSZ", "METFIPS", "INDIVIDCC", "FAMINC", 
"INTTYPE", "PERNUM", "WTFINL", "CPSIDP", "RELATE", "AGE", "SEX", "RACE", "MARST","FAMSIZE", "NCHILD", "BPL", "YRIMMIG", "CITIZEN", "MBPL", "FBPL", "NATIVITY", 
"HISPAN", "EMPSTAT", "LABFORCE", "OCC", "OCC2010", "OCC1990", "IND1990", "OCC1950", "CLASSWKR", "WKSTAT", "EDUC", "DIFFANY", "HOURWAGE", "PAIDHOUR", 
"UNION", "EARNWEEK", "UHRSWORKORG", "WKSWORKORG", "ELIGORG", "OTPAY", "VOWHYNOT", "VOYNOTREG", "VOTEHOW", "VOTEWHEN", "VOREGHOW", "VOREG95", "VOTERES", "VOTERESP",
"VOTED", "VOREG", "VOSUPPWT"
        
The national election turnout rates are obtained from the US Elections Project, managed by Michael P. McDonald at the University of Florida (see [US Elections Project](https://www.electproject.org/election-data/voter-turnout-data)).


2. The following is a master file that calls all individual scripts needed for the analysis. Running this file will execute all the necessary steps in order to produce all final outputs:

   `script/run-for-all-analysis.R`



## Roadmap for the Code

### 0. `script/00_functions.R`

This script defines various helper functions that are used throughout the analysis.  

- Some functions were originally written by [Benjamin Elbers](https://htmlpreview.github.io/?https://github.com/elbersb/weightedcontrasts/blob/master/doc/holford1983.html).
- Other functions are based on the APCI package developed by Jiahui Xu and Liying Luo (see Xu, J., & Luo, L. (2022). *APCI: an R and Stata package for visualizing and analyzing age-period-cohort data.* The R Journal, 14(2), 77).

Both sets of functions have been slightly modified to suit the specific needs of this project.

- **Takes in**:  
  None

- **Calls**:  
  None

- **Yields**:  
  None

---

### 1. `script/01_data_pre_processing.R`

Imports the 1976-2020 CPS-VRS data, recodes variables for analysis, and keeps only the relevant ones.

- **Takes in**:  
  `rawdata/cps_vrs_1976_2020.csv`
- **Calls**:  
  NONE
- **Yields**:  
  `data/cps.76.20.recoded.rds`

---

### 2. `script/02_iptw_adjustment_for_nonresp.R`

Models election-year-specific non-response for the turnout question and extracts missing probabilities from the models to construct missingness-corrected weights.

- **Takes in**:  
  - `data/cps.76.20.recoded.rds`
- **Calls**:  
  NONE
- **Yields**:  
  - `data/cps.76.20.recoded.wt.all.rds` (Includes all respondents)  
  - `data/cps.76.20.recoded.wt.voted.rds` (Restricted to non-missing respondents)

---

### 3. `script/03_descriptives.R`

Generates descriptive tables/figures (specifically, Figures 1, 2, 3, S1(a), and Table 1 in the paper).

- **Takes in**:  
  - `data/cps.76.20.recoded.wt.voted.rds`  
  - `data/cps.76.20.recoded.wt.all.rds`  
  - `rawdata/mcdonald_vep_turnout_1789_present.csv` (from the US Elections Project; see [here](https://www.electproject.org/))
- **Calls**:  
  NONE
- **Yields**:  
  - `output/tex/apc_freq.tex` (Table S1)  
  - `output/figure/heatmap_apc.pdf` (Figure 1)  
  - `output/figure/non_resp_by_year.pdf` (Figure S1(a))  
  - `output/figure/turnout_trends.pdf` (Figure S1(b))  
  - `output/figure/C_by_A.pdf` (Figure 2)  
  - `output/figure/C_by_P.pdf` (Figure 3)

---

### 4. `script/04_apc_bounding_apci_analysis.R`

Estimates APC bounding analysis and APC-I model using linear probability models.

- **Takes in**:  
  - `data/cps.76.20.recoded.wt.rds`
- **Calls**:  
  - `script/00_functions.R`
- **Yields**:  
  - `output/tex/apci_all.tex` (Table S3)  
  - `output/tex/bounding_A_nonlin.tex` (Table S7)  
  - `output/tex/bounding_P_nonlin.tex` (Table S7)  
  - `output/tex/bounding_C_nonlin.tex` (Table S7)  
  - `output/tex/bounding_A.tex` (Table S3)  
  - `output/tex/bounding_P.tex` (Table S3)  
  - `output/tex/bounding_C.tex` (Table S3)  
  - `output/figures/nonlinear_apc.pdf` (Figure 4)  
  - `output/figures/total_apc_a12_with_api.pdf` (Figure 5(a))  
  - `output/figures/total_apc_a12_c1_with_api.pdf` (Figure 5(b))  
  - `output/figures/total_apc_a12_c1_with_api.pdf` (Figure 5(c))  
  - `output/figures/intra_cohort_apci.pdf` (Figure S3)

---

### 5. `script/05_apc_bounding_apci_analysis_logit.R`

Estimates APC bounding analysis and APC-I model using logistic regression models (otherwise identical to `04_apc_bounding_apci_analysis.R`).

- **Takes in**:  
  - `data/cps.76.20.recoded.wt.rds`
- **Calls**:  
  - `script/00_functions.R`
- **Yields**:  
  - `output/figures/nonlinear_apc_logit.pdf`  
  - `output/figures/total_apc_a12_with_api_logit.pdf` (Figure S2(a))  
  - `output/figures/total_apc_a12_c1_with_api_logit.pdf` (Figure S2(b))  
  - `output/figures/total_apc_a12_c2_with_apci_logit.pdf` (Figure S2(c))


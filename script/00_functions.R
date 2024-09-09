###############################################################################
###                                                                         ###
### Function Definitions                                                    ###
###                                                                         ###
###############################################################################

# define helper functions -------------------------------------------------
# Originally written by Benjamin Elbers 
# (source: https://htmlpreview.github.io/?https://github.com/elbersb/weightedcontrasts/blob/master/doc/holford1983.html),
# with slight modifications for this project.

# define helper functions
get_deviations_w_intercept <- function(model, set, contrasts) {
  extract_coefs <- coef(model)[grepl(set, names(coef(model)))]
  extract_coefs <- extract_coefs[2:length(extract_coefs)]
  extract_coefs[is.na(extract_coefs)] <- 0
  deviations <- contrasts[, 2:(1 + length(extract_coefs))] %*% extract_coefs
  coef(model)[1] + deviations[, 1]
}

get_deviations <- function(model, set, contrasts) {
  extract_coefs <- coef(model)[grepl(set, names(coef(model)))]
  extract_coefs <- extract_coefs[2:length(extract_coefs)]
  deviations <- contrasts[, 2:(1 + length(extract_coefs)), drop = FALSE] %*% extract_coefs
  deviations[, 1]
}


get_total_effect <- function(contrasts, linear_coef, nonlinear_coefs) {
  linear_contrast <- contrasts[, 1]
  nonlinear_contrasts <- contrasts[, 2:ncol(contrasts)]
  linear_contrast * linear_coef + nonlinear_contrasts %*% nonlinear_coefs
}

get_linear_contrast <- function(vector) {
  1:n_distinct(vector) - 1/2 * n_distinct(vector) - 1/2
}


# inverse logit function --------------------------------------------------

invlogit <- function(x){
  invlogit <-  exp(x)/(1+exp(x))
  return(invlogit)
}


# misc. helper functions --------------------------------------------------
# Functions originally developed by Xu and Luo for the -APCI- package 
# (see Xu, J., & Luo, L. (2022). APCI: an R and Stata package for visualizing and analyzing age-period-cohort data. The R Journal, 14(2), 77),
# with slight modifications to incorporate intercepts into the estimates

# function to display the bounds of the age, period, and cohort effects
print.boounds <- function(bounds) {
  tibble(age.lower = bounds[1], age.upper = bounds[2], period.lower = bounds[3], 
         period.upper = bounds[4], cohort.lower = bounds[5], cohort.upper = bounds[6])
}



# apci cohort deviation with intercept ------------------------------------

cohortdeviation_w_intercept <- function (A, P, C, model = temp6, weight = "wt", covariate)
{
  r6 = model$coefficients[stringr::str_detect(names(model$coefficients), 
                                              "acc|pcc|(Intercept)")]
  r6se = summary(model)$coef[stringr::str_detect(names(model$coefficients), 
                                                 "acc|pcc|(Intercept)"), "Std. Error"]
  r6p = summary(model)$coef[stringr::str_detect(names(model$coefficients), 
                                                "acc|pcc|(Intercept)"), "Pr(>|t|)"]
  T = array(rep(0, A * P * (A - 1) * (P - 1)), dim = c(A * 
                                                         P, (A - 1) * (P - 1)))
  ind1 = A * 1:(P - 1)
  ind2 = (A * (P - 1) + 1):(A * P - 1)
  ind3 = A * P
  ind = c(ind1, ind2, ind3)
  newind = 1:(A * P)
  newind = newind[-ind]
  T[newind, ] = diag((A - 1) * (P - 1))
  T[ind1, ] = -diag(P - 1)[, rep(1:(P - 1), each = A - 1)]
  T[ind2, ] = -diag(A - 1)[, rep(1:(A - 1), P - 1)]
  T[ind3, ] = rep(1, (A - 1) * (P - 1))
  T <- cbind(rep(1, nrow(T)), T) # bind the column of 1 for the intercept
  row_ind <- stringr::str_detect(rownames(vcov(model)),  "(Intercept)|^acc([0-9])*:pcc([0-9])*$") 
  col_ind <- stringr::str_detect(colnames(vcov(model)), "(Intercept)|^acc([0-9])*:pcc([0-9])*$")
  row_ind_r6 <- stringr::str_detect(names(r6), "(Intercept)|^acc([0-9])*:pcc([0-9])*$")
  col_ind_r6 <- stringr::str_detect(names(r6), "(Intercept)|^acc([0-9])*:pcc([0-9])*$")
  iatemp = vcov(model)[row_ind, col_ind]
  iavcov = T %*% iatemp %*% t(T)
  df = model$df.residual
  iaesti = as.vector(T %*% r6[row_ind_r6])
  iase = sqrt(diag(iavcov))
  iap = pt(-abs(iaesti/iase), df) * 2
  cindex <- sapply(1:P, function(j) {
    seq((A + j - 1), j, -1)
  })
  sig = rep("   ", (A * P))
  sig[iap < 0.05] = "*  "
  sig[iap < 0.01] = "** "
  sig[iap < 0.001] = "***"
  iasig = sig
  cohortindex = as.vector(cindex)
  ia = as.data.frame(cbind(iaesti, iase, iap, iasig, cohortindex))
  cohortint <- sapply(1:C, function(k) {
    O = sum(cindex == k)
    k1 = rep(1/O, O)
    k2 = rep(0, A * P)
    k2[cindex == k] = k1
    contresti = k2 %*% iaesti
    contrse = sqrt(t(k2) %*% iavcov %*% k2)
    t = contresti/contrse
    if (t > 0) {
      p = 2 * pt(t, df, lower.tail = F)
    }
    else {
      p = 2 * pt(t, df, lower.tail = T)
    }
    sig <- "   "
    if (p < 0.05) {
      sig <- "*  "
    }
    if (p < 0.01) {
      sig <- "** "
    }
    if (p < 0.001) {
      sig <- "***"
    }
    c(contresti, contrse, t, p, sig)
  }) %>% t %>% as.data.frame %>% `colnames<-`(c("cohort_average", 
                                                "cohort_average_se", "cohort_average_t", 
                                                "cohort_average_p", "sig"))
  cohortint$cohort_group = seq(1, C)
  cohortint = cohortint[, c("cohort_group", "cohort_average", 
                            "cohort_average_se", "cohort_average_t", 
                            "cohort_average_p", "sig")]
  poly = 1
  cohortslope <- sapply((poly + 1):(C - poly), function(k) {
    o = sum(cindex == k)
    k1 = contr.poly(o)
    k2 = rep(0, A * P)
    k2[cindex == k] = k1[, poly]
    contresti = k2 %*% iaesti
    contrse = sqrt(t(k2) %*% iavcov %*% k2)
    t = contresti/contrse
    if (t > 0) {
      p = 2 * pt(t, df, lower.tail = F)
    }
    else {
      p = 2 * pt(t, df, lower.tail = T)
    }
    sig <- "   "
    if (p < 0.05) {
      sig <- "*  "
    }
    if (p < 0.01) {
      sig <- "** "
    }
    if (p < 0.001) {
      sig <- "***"
    }
    c(contresti, contrse, t, p, sig)
  }) %>% t %>% as.data.frame %>% `colnames<-`(c("cohort_slope", 
                                                "cohort_slope_se", "cohort_slope_t", "cohort_slope_p", 
                                                "sig"))
  cohortslope <- rbind(NA, cohortslope, NA)
  cohortslope$cohort_group = seq(1, C)
  cohortslope = cohortslope[, c("cohort_group", "cohort_slope", 
                                "cohort_slope_se", "cohort_slope_t", "cohort_slope_p", 
                                "sig")]
  list(cohort_average = cohortint, cohort_slope = cohortslope, 
       int_matrix = ia)
}


maineffect_w_intercept <- function (A, P, C, model = temp6) 
{
  
  r6 = model$coefficients[stringr::str_detect(names(model$coefficients), 
                                              "acc|pcc|(Intercept)")]
  
  r6se = summary(model)$coef[stringr::str_detect(names(model$coefficients), 
                                                 "acc|pcc|(Intercept)"), "Std. Error"]
  r6p = summary(model)$coef[stringr::str_detect(names(model$coefficients), 
                                                "acc|pcc|(Intercept)"), "Pr(>|t|)"]
  fullae = array(rep(0, A), dim = c(A, 1))
  fullae
  fullas = array(rep(0, A), dim = c(A, 1))
  fullas
  S1 = array(rep(0, A * (A - 1)), dim = c(A, (A - 1)))
  ind = A * 1:(A - 1)
  newind = 1:(A * (A - 1))
  newind = newind[-ind]
  newind
  S1[newind] = diag(A - 1)
  S1[ind] = rep(-1, (A - 1))
  S1 <- cbind(rep(1, nrow(S1)), S1) # bind the column of 1 for the intercept
  fullae = as.vector(S1 %*% model$coef[stringr::str_detect(names(model$coef), 
                                                           "(Intercept)|^acc([0-9])*$")])
  row_ind <- stringr::str_detect(rownames(vcov(model)), "(Intercept)|^acc([0-9])*$")
  col_ind <- stringr::str_detect(colnames(vcov(model)), "(Intercept)|^acc([0-9])*$")
  fullas = sqrt(diag(S1 %*% vcov(model)[row_ind, col_ind] %*% 
                       t(S1)))
  fullat = fullae/fullas
  fullap = pt(-abs(fullat), df.residual(model)) * 2
  sig = rep("   ", A)
  sig[fullap < 0.05] = "*  "
  sig[fullap < 0.01] = "** "
  sig[fullap < 0.001] = "***"
  fullasig = sig
  fulla = cbind(fullae, fullas, fullap, fullasig)
  fullpe = array(rep(0, P), dim = c(P, 1))
  fullps = array(rep(0, P), dim = c(P, 1))
  S2 = array(rep(0, P * (P - 1)), dim = c(P, (P - 1)))
  ind = P * 1:(P - 1)
  newind = 1:(P * (P - 1))
  newind = newind[-ind]
  S2[newind] = diag(P - 1)
  S2[ind] = rep(-1, (P - 1))
  S2 <- cbind(rep(1, nrow(S2)), S2) # bind the column of 1 for the intercept
  fullpe = as.vector(S2 %*% model$coef[stringr::str_detect(names(model$coef), 
                                                           "(Intercept)|^pcc([0-9])*$")])
  row_ind <- stringr::str_detect(rownames(vcov(model)), "(Intercept)|^pcc([0-9])*$")
  col_ind <- stringr::str_detect(colnames(vcov(model)), "(Intercept)|^pcc([0-9])*$")
  fullps = sqrt(diag(S2 %*% vcov(model)[row_ind, col_ind] %*% 
                       t(S2)))
  fullpt = fullpe/fullps
  fullpp = pt(-abs(fullpt), df.residual(model)) * 2
  sig = rep("   ", P)
  sig[fullpp < 0.05] = "*  "
  sig[fullpp < 0.01] = "** "
  sig[fullpp < 0.001] = "***"
  fullpsig = sig
  fullp = cbind(fullpe, fullps, fullpp, fullpsig)
  inte = as.vector(r6[1])
  intse = r6se[1]
  intp = r6p[1]
  intsig = rep("   ", 1)
  intsig[r6p[1] < 0.05] = "*  "
  intsig[r6p[1] < 0.01] = "** "
  intsig[r6p[1] < 0.001] = "***"
  fullint = cbind(inte, intse, intp, intsig)
  maineff = rbind(fullint, fulla, fullp)
  colnames(maineff) = c("estimate", "se", "p", 
                        "sig")
  rownames(maineff) = c()
  age_results <- maineff[2:(A + 1), ]
  colnames(age_results) = c("age_estimate", "age_se", 
                            "age_p", "sig")
  rownames(age_results) = c()
  period_results <- maineff[(A + 2):(A + P + 1), ]
  colnames(period_results) = c("period_estimate", "period_se", 
                               "period_p", "sig")
  rownames(period_results) = c()
  list(intercept = maineff[1, ], age_effect = cbind(age_group = seq(1, 
                                                                    A), age_results), period_effect = cbind(period_group = seq(1, 
                                                                                                                               P), period_results))
}

# helper function to get apc_i estimates for lpm model
get_apci_identity <- function(y, a, p, c, wt, data){
  
  apc_i <- temp_model(data = data, outcome = paste(y), age = paste(a), period = paste(p), 
                      cohort = paste(c), weight = paste(wt), family = "gaussian")
  ap <- maineffect_w_intercept(A = apc_i$A, P = apc_i$P, C = apc_i$C, model = apc_i$model)  
  c <- cohortdeviation_w_intercept(A =  apc_i$A, P = apc_i$P, C = apc_i$C, model = apc_i$model)
  
  A.main <- as.double(ap$age_effect[,2]) * 100
  P.main <-  as.double(ap$period_effect[,2]) * 100
  C.dev <- as.double(c$cohort_average[,2]) * 100
  C.slope <-  as.double(c$cohort_slope[,2]) * 100
  
  A.main.ci.lb <- A.main - 1.96 * as.double(ap$age_effect[,3]) * 100
  A.main.ci.ub <- A.main + 1.96 * as.double(ap$age_effect[,3]) * 100
  P.main.ci.lb <- P.main - 1.96 * as.double(ap$period_effect[,3]) * 100
  P.main.ci.ub <- P.main + 1.96 * as.double(ap$period_effect[,3]) * 100
  C.dev.ci.lb <- C.dev - 1.96 * as.double(c$cohort_average[,3]) * 100
  C.dev.ci.ub <- C.dev + 1.96 * as.double(c$cohort_average[,3]) * 100
  C.slope.ci.lb <- C.slope - 1.96 * as.double(c$cohort_slope[,3]) * 100
  C.slope.ci.ub <- C.slope + 1.96 * as.double(c$cohort_slope[,3]) * 100
  
  int <- apc_i$model$coefficients[1] * 100
  A <- cbind(A.main, A.main.ci.lb, A.main.ci.ub)
  P <- cbind(P.main, P.main.ci.lb, P.main.ci.ub)
  C.avg <- cbind(C.dev, C.dev.ci.lb, C.dev.ci.ub)
  C.slope <- cbind(C.slope, C.slope.ci.lb, C.slope.ci.ub)
  
  output <- list(intercept = int, age_est = A, period_est = P, cohort_est = C.avg, cohort_slope_st = C.slope)
  return(output)
}

# helper function to get apc_i estimates for logit model
get_apci_logit <- function(y, a, p, c, wt, data){
  
  apc_i <- temp_model(data = data, outcome = paste(y), age = paste(a), period = paste(p), 
                      cohort = paste(c), weight = paste(wt), family = "quasibinomial")
  ap <- maineffect_w_intercept(A = apc_i$A, P = apc_i$P, C = apc_i$C, model = apc_i$model)  
  c <- cohortdeviation_w_intercept(A =  apc_i$A, P = apc_i$P, C = apc_i$C, model = apc_i$model)
  
  # extract AP effects in predicted probabilities
  int <- invlogit(apc_i$model$coefficients[1]) * 100
  A <- invlogit(as.double(ap$age_effect[,2])) * 100
  P <- invlogit(as.double(ap$period_effect[,2])) * 100
  C.avg <- invlogit(as.double(c$cohort_average[,2])) * 100
  C.slope <- invlogit(as.double(c$cohort_slope[,2])) * 100
  
  output <- list(intercept = int, age_est = A, period_est = P, cohort_est = C.avg, cohort_slope_est = C.slope)
  return(output)
  
}




### END OF CODE ###
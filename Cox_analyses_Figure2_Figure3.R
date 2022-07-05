################################################################
# 20-Year Benefit from Adjuvant Goserelin and Tamoxifen in
# Premenopausal Breast Cancer Patients in a Controlled Randomized
# Clinical Trial (STO-5)
# Author: Annelie Johansson
# Date: 1 July 2021
################################################################

library(openxlsx)
library(dplyr)
library(survival)

# Load data ----
# ER-positive patients in the STO-5 trial
load("Data/STO5_trial.RData")
nrow(STO5) # 584

table(STO5$treatment, useNA = "always")
#  C    G    T  G+T <NA> 
# 145  155  135  149    0 

# function
get_res <- function(fit){
  temp <- summary(fit)$conf.int
  ind <- grep("var", rownames(temp))
  HR_95CI <- formatC(temp[ind, c(1,3,4)], format = "f", digits = 2)
  HR_95CI <- c(paste0(HR_95CI[1], " (", HR_95CI[2], "-", HR_95CI[3], ")"))
  return(HR_95CI)
}

# Adjustments:
# age, randomization year (YR_STO5), tumor size, tumor grade, PR status, HER2 status, Ki67 status,
# lymph node status/chemotherapy/radiotherapy (variable "Group"), type of surgery.

#### ALL PATIENTS ####
res_all <- NULL

# * Goserelin vs Control ----
data_temp <- subset(STO5, (treatment == "G" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit))

# * Tamoxifen vs  Control ----
data_temp <- subset(STO5, (treatment == "T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit))

# * Gos+Tam vs Control ----
data_temp <- subset(STO5, (treatment == "G+T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit))

#### GENOMIC LOW RISK PATIENTS ####
STO5_low <- subset(STO5, MammaPrint == "Low")
res_low <- NULL

# * Goserelin vs Control ----
data_temp <- subset(STO5_low, (treatment == "G" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit))

# * Tamoxifen vs  Control ----
data_temp <- subset(STO5_low, (treatment == "T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit))

# * Gos+Tam vs Control ----
data_temp <- subset(STO5_low, (treatment == "G+T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit))

#### GENOMIC HIGH RISK PATIENTS ####
STO5_high <- subset(STO5, MammaPrint == "High")
# very few grade 1, so we put together grade 1 and 2
STO5_high$grade <- factor(STO5_high$grade , levels = c("1", "2", "3"), labels = c("1-2", "1-2", "3"))
res_high <- NULL

# * Goserelin vs Control ----
data_temp <- subset(STO5_high, (treatment == "G" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))

# * Tamoxifen vs  Control ----
data_temp <- subset(STO5_high, (treatment == "T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))

# * Gos+Tam vs Control ----
data_temp <- subset(STO5_high, (treatment == "G+T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
data_temp <- data_temp[-which(data_temp$size == "Unknown"),] # doesn't converge otherwise!
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))


# get results
res_final <- cbind(res_all, res_low, res_high)
rownames(res_final) <- c("G_vs_C", "T_vs_C", "G+T_vs_C")
colnames(res_final) <- c("All Patients", "Genomic Low Risk Patients", "Genomic High Risk Patients")

#          All Patients      
# G_vs_C   "0.49 (0.32-0.75)" *
# T_vs_C   "0.57 (0.38-0.87)" *
# G+T_vs_C "0.63 (0.42-0.94)" *
# Genomic Low Risk Patients
# G_vs_C   "0.75 (0.39-1.46)"       
# T_vs_C   "0.24 (0.10-0.60)" *     
# G+T_vs_C "0.63 (0.33-1.23)"       
# Genomic High Risk Patients
# G_vs_C   "0.24 (0.10-0.54)" *       
# T_vs_C   "0.87 (0.41-1.85)"        
# G+T_vs_C "0.69 (0.35-1.35)"

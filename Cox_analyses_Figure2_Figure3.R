################################################################
# Endocrine therapy & the 70-gene risk signature (MammaPrint)
# in premenopausal ER-positive patients.
# Author: Annelie Johansson
# Date: 8 December 2021
################################################################

library(openxlsx)
library(dplyr)
library(survival)

# Load data ----

setwd("/Volumes/Annelie Encrypted 2/STO25/Zoladex_STO5/")
load("Data_temp/STO5_prepared.RData") # * UPDATE
STO <- STO5
nrow(STO) # 584

table(STO$treatment, useNA = "always")
#  C    G    T  G+T <NA> 
# 145  155  135  149    0 

# function
get_res <- function(fit, data = data_temp){
  temp <- summary(fit)$conf.int
  ind <- grep("var", rownames(temp))
  HR_95CI <- formatC(temp[ind, c(1,3,4)], format = "f", digits = 2)
  HR_95CI <- c(paste0(HR_95CI[1], " (", HR_95CI[2], "-", HR_95CI[3], ")"))
  tr <- which(data$var == levels(data$var)[2])
  ref <- which(data$var == levels(data$var)[1])
  tr_res <- c(length(tr), length(which(data$MetBC20yr[tr] == 1)), HR_95CI)
  ref_res <- c(length(ref), length(which(data$MetBC20yr[ref] == 1)), "1.00 (Ref)")
  res <- rbind(tr_res, ref_res)
  return(res)
}


#### ALL PATIENTS ####
res_all <- NULL

# * Goserelin vs Control ----
data_temp <- subset(STO, (treatment == "G" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit)[1,])

# * Tamoxifen vs  Control ----
data_temp <- subset(STO, (treatment == "T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
#data_temp <- data_temp[-which(data_temp$grade == "98"),]
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit)[1,]) # 0.57 (0.38-0.87)
# doesn't converge for grade 98, but doesn't change results to keep or not

# * Gos+Tam vs Control ----
data_temp <- subset(STO, (treatment == "G+T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
# data_temp <- data_temp[-which(data_temp$HER2status == "Unknown"),]
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit)) # 0.63 (0.42-0.94)
# doesn't converge for HER2 unknown, but doesn't change results to keep or not

#### GENOMIC LOW RISK ####
STO_low <- subset(STO, MammaPrint == "Low")
res_low <- NULL

# * Goserelin vs Control ----
data_temp <- subset(STO_low, (treatment == "G" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit)[1,])

# * Tamoxifen vs  Control ----
data_temp <- subset(STO_low, (treatment == "T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
# data_temp <- data_temp[-which(data_temp$HER2status == "Unknown"),]
# data_temp <- data_temp[-which(data_temp$size == "Unknown"),]
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit)[1,]) # 0.24 (0.10-0.60)
# does not converge for size and HER2 Unknown, excluding gives same estimates

# * Gos+Tam vs Control ----
data_temp <- subset(STO_low, (treatment == "G+T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
# data_temp <- data_temp[-which(data_temp$HER2status == "Unknown"),]
# data_temp <- data_temp[-which(data_temp$size == "Unknown"),]
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit)) # 0.63 (0.33-1.23)
# does not converge for size and HER2 Unknown, excluding gives same estimates

#### GENOMIC HIGH RISK ####
STO_high <- subset(STO, MammaPrint == "High")
# very few grade 1, so we put together grade 1 and 2
STO_high$grade <- factor(STO_high$grade , levels = c("1", "2", "3"), labels = c("1-2", "1-2", "3"))
res_high <- NULL

# * Goserelin vs Control ----
data_temp <- subset(STO_high, (treatment == "G" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit)[1,])

# * Tamoxifen vs  Control ----
data_temp <- subset(STO_high, (treatment == "T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit)[1,])

# * Gos+Tam vs Control ----
data_temp <- subset(STO_high, (treatment == "G+T" | treatment == "C"))
data_temp$var <- droplevels(data_temp$treatment)
data_temp <- data_temp[-which(data_temp$size == "Unknown"),] # doesn't converge otherwise!
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))

# get results
res_final <- rbind(res_all, res_low, res_high)

rownames(res_final) <- rep(c("G_vs_C", "T_vs_C", "G+T_vs_C", "Ref"), 3)
colnames(res_final) <- c("Patients", "DR", "HR_95CI")

# write.table(res_final[1:4,], file = "Tables/cox_analyses_Figure2.txt", sep = "\t")
# write.table(res_final[5:12,], file = "Tables/cox_analyses_Figure3.txt", sep = "\t")

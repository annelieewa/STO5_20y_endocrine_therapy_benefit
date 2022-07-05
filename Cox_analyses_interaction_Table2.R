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
#  145  155  135  149    0 

#### Interaction test between Goserelin and Tamoxifen ####

# function
get_res <- function(fit ){
  temp <- summary(fit)$conf.int
  ind <- grep("var", rownames(temp))
  HR_95CI <- formatC(temp[ind, c(1,3,4)], format = "f", digits = 2)
  HR_95CI <- c(paste0(HR_95CI[1], " (", HR_95CI[2], "-", HR_95CI[3], ")"))
  return(HR_95CI)
}

# Create Gos and Tam variables in STO5
STO5$Gos <- factor(ifelse(STO5$treatment == "G" | STO5$treatment == "G+T", "Gos+", "Gos-"))
STO5$Tam <- factor(ifelse(STO5$treatment == "T" | STO5$treatment == "G+T", "Tam+", "Tam-"))
table(STO5$Gos, STO5$Tam, useNA = "always")
#        Tam- Tam+ <NA>
# Gos-  145  135    0
# Gos+  155  149    0
# <NA>    0    0    0

#### ALL PATIENTS ####
res_all <- NULL

### Effect of Goserelin ----

# * Effect in patients treated with Tamoxifen ----
# G+T vs T
data_temp <- subset(STO5, Tam == "Tam+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos-Tam+"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit))

# * Effect in patients NOT treated with Tamoxifen ----
# G vs C
data_temp <- subset(STO5, Tam == "Tam-")
data_temp$var <- factor(ifelse(data_temp$treatment == "G", "Gos+Tam-", "Gos-Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit))

### Effect of Tamoxifen ----

# * Effect in patients treated with Goserelin ----
# G+T vs G
data_temp <- subset(STO5, Gos == "Gos+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos+Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit))

# * Effect in patients NOT treated with Goserelin ----
# T vs C
data_temp <- subset(STO5, Gos == "Gos-")
data_temp$var <- factor(ifelse(data_temp$treatment == "T", "Gos-Tam+", "Gos-Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit))

# * Interaction Goserelin and Tamoxifen ----
data_temp <- STO5
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ Tam*Gos + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, formatC(summary(fit)$coefficients["TamTam+:GosGos+",5], format = "f", digits = 3))


#### GENOMIC LOW RISK PATIENTS ####
STO5_low <- subset(STO5, MammaPrint == "Low")
STO5_low$grade <- droplevels(STO5_low$grade) # none 98
res_low <- NULL

### Effect of Goserelin ----

# * Effect in patients treated with Tamoxifen ----
# G+T vs T
data_temp <- subset(STO5_low, Tam == "Tam+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos-Tam+"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit))

# * Effect in patients NOT treated with Tamoxifen ----
# G vs C
data_temp <- subset(STO5_low, Tam == "Tam-")
data_temp$var <- factor(ifelse(data_temp$treatment == "G", "Gos+Tam-", "Gos-Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit))

### Effect of Tamoxifen ----

# * Effect in patients treated with Goserelin ----
# G+T vs G
data_temp <- subset(STO5_low, Gos == "Gos+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos+Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit))

# * Effect in patients NOT treated with Goserelin ----
# T vs C
data_temp <- subset(STO5_low, Gos == "Gos-")
data_temp$var <- factor(ifelse(data_temp$treatment == "T", "Gos-Tam+", "Gos-Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit))

# * Interaction Goserelin and Tamoxifen ----
data_temp <- STO5_low
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ Tam*Gos + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, formatC(summary(fit)$coefficients["TamTam+:GosGos+",5], format = "f", digits = 3))


#### GENOMIC HIGH RISK PATIENTS ####
STO5_high <- subset(STO5, MammaPrint == "High")
# very few grade 1, so we put together grade 1 and 2
STO5_high$grade <- factor(STO5_high$grade , levels = c("1", "2", "3"), labels = c("1-2", "1-2", "3"))
res_high <- NULL

### Effect of Goserelin ----

# * Effect in patients treated with Tamoxifen ----
# G+T vs T
data_temp <- subset(STO5_high, Tam == "Tam+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos-Tam+"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))

# * Effect in patients NOT treated with Tamoxifen ----
# G vs C
data_temp <- subset(STO5_high, Tam == "Tam-")
data_temp$var <- factor(ifelse(data_temp$treatment == "G", "Gos+Tam-", "Gos-Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))

### Effect of Tamoxifen ----

# * Effect in patients treated with Goserelin ----
# G+T vs G
data_temp <- subset(STO5_high, Gos == "Gos+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos+Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))

# * Effect in patients NOT treated with Goserelin ----
# T vs C
data_temp <- subset(STO5_high, Gos == "Gos-")
data_temp$var <- factor(ifelse(data_temp$treatment == "T", "Gos-Tam+", "Gos-Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))

# * Interaction Goserelin and Tamoxifen ----
data_temp <- STO5_high
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ Tam*Gos + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, formatC(summary(fit)$coefficients["TamTam+:GosGos+",5], format = "f", digits = 3))


# get results
res_final <- cbind(res_all, res_low, res_high)
rownames(res_final) <- c("Gos_in_Tam", "Gos_without_Tam", "Tam_in_Gos", "Tam_without_Gos", "Interaction_G_T")
colnames(res_final) <- c("All Patients", "Genomic Low Risk Patients", "Genomic High Risk Patients")

#                 All Patients      
# Gos_in_Tam      "1.01 (0.65-1.58)"
# Gos_without_Tam "0.49 (0.32-0.75)" *
# Tam_in_Gos      "1.27 (0.82-1.95)"
# Tam_without_Gos "0.57 (0.38-0.87)" *
# Interaction_G_T "0.016"           
# Genomic Low Risk Patients
# Gos_in_Tam      "2.20 (0.99-4.89)"       
# Gos_without_Tam "0.75 (0.39-1.46)"       
# Tam_in_Gos      "0.85 (0.45-1.62)"       
# Tam_without_Gos "0.24 (0.10-0.60)" *       
# Interaction_G_T "0.080"                  
# Genomic High Risk Patients
# Gos_in_Tam      "0.65 (0.30-1.40)"        
# Gos_without_Tam "0.24 (0.10-0.54)" *      
# Tam_in_Gos      "3.36 (1.39-8.07)" *        
# Tam_without_Gos "0.87 (0.41-1.85)"        
# Interaction_G_T "0.006" 


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

#### Interaction test between Goserelin and Tamoxifen ####

table(STO$treatment, useNA = "always")
#  C    G    T  G+T <NA> 
#  145  155  135  149    0 

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

# Create Gos and Tam variables in STO
STO$Gos <- factor(ifelse(STO$treatment == "G" | STO$treatment == "G+T", "Gos+", "Gos-"))
STO$Tam <- factor(ifelse(STO$treatment == "T" | STO$treatment == "G+T", "Tam+", "Tam-"))
table(STO$Gos, STO$Tam, useNA = "always")
#        Tam- Tam+ <NA>
# Gos-  145  135    0
# Gos+  155  149    0
# <NA>    0    0    0

#### ALL PATIENTS ####
res_all <- NULL

### Effect of Goserelin ----

# * Effect in patients treated with Tamoxifen ----
# G+T vs T
data_temp <- subset(STO, Tam == "Tam+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos-Tam+"))
# data_temp <- data_temp[-which(data_temp$grade == "98"),]
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit)) # 1.01 (0.65-1.58)
# does not converge for grade 98; excluding gives same estimates

# * Effect in patients NOT treated with Tamoxifen ----
# G vs C
data_temp <- subset(STO, Tam == "Tam-")
data_temp$var <- factor(ifelse(data_temp$treatment == "G", "Gos+Tam-", "Gos-Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit))

### Effect of Tamoxifen ----

# * Effect in patients treated with Goserelin ----
# G+T vs G
data_temp <- subset(STO, Gos == "Gos+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos+Tam-"))
# data_temp <- data_temp[-which(data_temp$HER2status == "Unknown"),]
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit)) # 1.27 (0.82-1.95)
# does not converge for HER2 Unknown; excluding gives same estimates

# * Effect in patients NOT treated with Goserelin ----
# T vs C
data_temp <- subset(STO, Gos == "Gos-")
data_temp$var <- factor(ifelse(data_temp$treatment == "T", "Gos-Tam+", "Gos-Tam-"))
# data_temp <- data_temp[-which(data_temp$grade == "98"),]
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_all <- rbind(res_all, get_res(fit)) # 0.57 (0.38-0.87)
# does not converge for grade 98; excluding gives same estimates

# * Interaction Goserelin and Tamoxifen ----
data_temp <- STO
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ Tam*Gos + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
# data_temp <- data_temp[-which(data_temp$grade == "98"),]
res_all <- rbind(res_all, c("", "", formatC(summary(fit)$coefficients["TamTam+:GosGos+",5], format = "f", digits = 3)))
# does not converge for grade 98; excluding gives same estimates


#### LOW RISK-PATIENTS ####

STO_low <- subset(STO, MammaPrint == "Low")
STO_low$grade <- droplevels(STO_low$grade) # none 98

res_low <- NULL

### Effect of Goserelin ----

# * Effect in patients treated with Tamoxifen ----
# G+T vs T
data_temp <- subset(STO_low, Tam == "Tam+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos-Tam+"))
# data_temp <- data_temp[-which(data_temp$size == "Unknown"),]
# data_temp <- data_temp[-which(data_temp$HER2status == "Unknown"),]
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit)) # 2.20 (0.99-4.89)
# does not converge for size and HER2 Unknown; excluding gives same estimates

# * Effect in patients NOT treated with Tamoxifen ----
# G vs C
data_temp <- subset(STO_low, Tam == "Tam-")
data_temp$var <- factor(ifelse(data_temp$treatment == "G", "Gos+Tam-", "Gos-Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit))

### Effect of Tamoxifen ----

# * Effect in patients treated with Goserelin ----
# G+T vs G
data_temp <- subset(STO_low, Gos == "Gos+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos+Tam-"))
# data_temp <- data_temp[-which(data_temp$Ki67status == "Unknown"),]
# data_temp <- data_temp[-which(data_temp$HER2status == "Unknown"),]
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit)) # 0.85 (0.45-1.62)
# does not converge for Ki-67 and HER2 Unknown; excluding gives same estimates

# * Effect in patients NOT treated with Goserelin ----
# T vs C
data_temp <- subset(STO_low, Gos == "Gos-")
data_temp$var <- factor(ifelse(data_temp$treatment == "T", "Gos-Tam+", "Gos-Tam-"))
# data_temp <- data_temp[-which(data_temp$size == "Unknown"),]
# data_temp <- data_temp[-which(data_temp$HER2status == "Unknown"),]
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, get_res(fit)) # 0.24 (0.10-0.60)
# does not converge for size and HER2 Unknown; excluding gives same estimates

# * Interaction Goserelin and Tamoxifen ----
data_temp <- STO_low
# data_temp <- data_temp[-which(data_temp$HER2status == "Unknown"),]
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ Tam*Gos + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_low <- rbind(res_low, c("", "", formatC(summary(fit)$coefficients["TamTam+:GosGos+",5], format = "f", digits = 3)))
# does not converge for HER2 Unknown; excluding gives same estimates


#### HIGH RISK-PATIENTS ####
STO_high <- subset(STO, MammaPrint == "High")
STO_high$grade <- factor(STO_high$grade , levels = c("1", "2", "3"), labels = c("1-2", "1-2", "3"))

res_high <- NULL

### Effect of Goserelin ----

# * Effect in patients treated with Tamoxifen ----
# G+T vs T
data_temp <- subset(STO_high, Tam == "Tam+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos-Tam+"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))

# * Effect in patients NOT treated with Tamoxifen ----
# G vs C
data_temp <- subset(STO_high, Tam == "Tam-")
data_temp$var <- factor(ifelse(data_temp$treatment == "G", "Gos+Tam-", "Gos-Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))

### Effect of Tamoxifen ----

# * Effect in patients treated with Goserelin ----
# G+T vs G
data_temp <- subset(STO_high, Gos == "Gos+")
data_temp$var <- factor(ifelse(data_temp$treatment == "G+T", "Gos+Tam+", "Gos+Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))

# * Effect in patients NOT treated with Goserelin ----
# T vs C
data_temp <- subset(STO_high, Gos == "Gos-")
data_temp$var <- factor(ifelse(data_temp$treatment == "T", "Gos-Tam+", "Gos-Tam-"))
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ var + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, get_res(fit))

# * Interaction Goserelin and Tamoxifen ----
data_temp <- STO_high
fit <- coxph(Surv(Met_20yr, MetBC20yr) ~ Tam*Gos + YR_STO5 + age + size + grade + PRstatus + HER2status + Ki67status + Group + surgery, data = data_temp)
res_high <- rbind(res_high, c("", "", formatC(summary(fit)$coefficients["TamTam+:GosGos+",5], format = "f", digits = 3)))


res_final <- cbind(res_all, res_low, res_high)
rownames(res_final) <- c(rep( "Gos_in_tam",2), rep("Gos_wo_Tam",2),
                         rep("Tam_in_Gos",2), rep("Tam_wo_Gos",2), "interaction_G_T")

# write.table(res_final, file = "Tables/table2_interaction.txt", sep = "\t")

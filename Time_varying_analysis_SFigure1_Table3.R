################################################################
# 20-Year Benefit from Adjuvant Goserelin and Tamoxifen in
# Premenopausal Breast Cancer Patients in a Controlled Randomized
# Clinical Trial (STO-5)
# Author: Annelie Johansson
# Date: 1 July 2021
################################################################

library(survival)
library(survminer)
library(rstpm2) 
library(dplyr)

# * Load data ----
# ER-positive patients in the STO-5 trial
load("Data/STO5_trial.RData")
nrow(STO5) # 584

# * Exclude patients with missing data ----
STO5_use <- subset(STO5, size != "Unknown" & grade != "98" & PRstatus != "Unknown"
                  & HER2status != "Unknown" & Ki67status != "Unknown")
nrow(STO5_use) # 563
STO5_use$size <- droplevels(STO5_use$size)
STO5_use$grade <- droplevels(STO5_use$grade)
STO5_use$PRstatus <- droplevels(STO5_use$PRstatus)
STO5_use$HER2status <- droplevels(STO5_use$HER2status)
STO5_use$Ki67status <- droplevels(STO5_use$Ki67status)

# * Select genomic low-risk patients ---- 
STO5_low <- subset(STO5_use, MammaPrint == "Low")
nrow(STO5_low) # 295

# * Select genomic high-risk patients ----
STO5_high <- subset(STO5_use, MammaPrint == "High")
nrow(STO5_high) # 154

#### Functions for graphs ####

#### Graph function A ----
func_graphA <- function(model, trial_arm, max, linecol, sec, title){
  
  hazard_treated <- predict(model, type = "hazard", newdata = cbind.data.frame(Met_20yr = sec, df1), var = "treatment", se.fit = TRUE) 
  hazard_untreated <- predict(model, type = "hazard", newdata = cbind.data.frame(Met_20yr = sec, df0), var = "treatment", se.fit = TRUE) 
  df <- data.frame(x = c(sec, sec),
                   est = c(hazard_treated$Estimate, hazard_untreated$Estimate),
                   treat = c(rep(trial_arm, length(sec)), rep("Control", length(sec))))
  df$treat <- factor(df$treat, levels = c(trial_arm, "Control"))
  
  graphA <- ggplot() +
    geom_line(data = df, aes(x = x, y = est, group = treat, linetype = treat, color = treat), size = 1) +
    xlab("Years Since Randomization") +
    ylab("Estimated Hazard Rates for DRFI")  +
    ggtitle(title) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 20.1), breaks = seq(0, 20, 5)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max)) +
    scale_color_manual(values = c(linecol, "black")) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       axis.text = element_text(size = 12),
                       axis.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.title = element_blank(),
                       text = element_text(size = 10),
                       legend.position = c(0.4,0.6)) #,
                       #plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 
  return(graphA)
}


#### Graph function B ----
func_graphB <- function(model, trial_arm, max, linecol, sec, title){
  
  HR <- predict(model, type = "hr", newdata = cbind.data.frame(Met_20yr = sec, df0), var = "treatment", se.fit = TRUE) 
  HR$upper[which(HR$upper > max)] <- max
  
  df <- data.frame(x = c(sec, sec), group = c(rep(trial_arm, length(sec)), rep("Control", length(sec))), HR_est = c(HR$Estimate, rep(1, length(sec))),
                   HR_low = c(HR$lower, rep(1, length(sec))), HR_upp = c(HR$upper, rep(1, length(sec))))
  df$group <- factor(df$group, levels = c(trial_arm, "Control"))
  df_treated <- data.frame(x = sec, HR_est = HR$Estimate, HR_low = HR$lower, HR_upp = HR$upper)
  
  events_0 <- model@data$Met_20yr[which(model@data$MetBC20yr == 1 & model@data$treatment == 0)]
  events_1 <- model@data$Met_20yr[which(model@data$MetBC20yr == 1 & model@data$treatment == 1)]
  events_0 <- events_0[order(events_0)]
  events_1 <- events_1[order(events_1)]
  if( any(events_0 %in% events_1) ){
    t <- which(events_0 %in% events_1)
    events_0[t] <- events_0[t]+0.05
  }
  
  df_treated_sig <- df_treated[which(df_treated$HR_upp <= 1),]
  
  graph <- ggplot() +
    geom_line(data = df, aes(x = x, y = HR_est, group = group, linetype = group, color = group), size = 0.8) +
    geom_ribbon(data = df_treated, aes(x = x, ymin = HR_low, ymax = HR_upp), alpha = 0.5, fill = "grey") +
    geom_ribbon(data = df_treated_sig, aes(x = x, ymin = HR_low, ymax = HR_upp), alpha = 0.5, fill = "gray70") +
    ggtitle(title) +
    xlab("Years Since Randomization") +
    ylab("Estimated HR for DRFI") +
    scale_color_manual(values = c(linecol, "black")) +
    scale_x_continuous(expand = c(0, 0),  limits = c(0, 20.1), breaks = seq(0, 20, 5)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,max), breaks = seq(0, max, 0.5)) +
    theme_bw() + theme(panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"),
                       axis.text = element_text(size = 12),
                       axis.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.title = element_blank(),
                       legend.position = c(0.32,0.72),
                       text = element_text(size = 10)) + 
    geom_segment(aes(x = events_0, y = 0, xend = events_0, yend = 0.1), col = "black") + 
    geom_segment(aes(x = events_1, y = 0, xend = events_1, yend = 0.1), col = linecol)
  return(graph)
}

########################################
# Time-varying analysis                #
# using flexible parametric modelling  # 
########################################

sec <- seq(0.05, 20, 0.1)

#### Tamoxifen in patients of low genomic risk ####
df0 <- data.frame(treatment = 0, age = "(44,50]", size = "T1c", grade = "2", Group = "STO5_R1",
                  PRstatus = "Positive", HER2status = "Negative", Ki67status = "Low", surgery = "breast_conserving")
df1 <- data.frame(treatment = 1, age = "(44,50]", size = "T1c", grade = "2", Group = "STO5_R1",
                  PRstatus = "Positive", HER2status = "Negative", Ki67status = "Low", surgery = "breast_conserving")

STO5_temp <- subset(STO5_low, treatment == "C" | treatment == "T")
STO5_temp$treatment <- as.character(STO5_temp$treatment)
STO5_temp$treatment <- ifelse(STO5_temp$treatment == "C", 0, 1)

model <- stpm2(Surv(Met_20yr, MetBC20yr) ~ treatment + factor(age) + factor(Group) +
                 factor(size) + factor(grade) + factor(PRstatus) + factor(HER2status) +
                 factor(Ki67status) + factor(surgery),
               data = STO5_temp, df = 2, tvc=list(treatment = 1), stata.stpm2.compatible=TRUE)
pred <- predict(model, type = "hr", newdata = cbind.data.frame(df0, Met_20yr = c(5,10,15,20)), var = "treatment", se.fit = TRUE) 
rownames(pred) <- paste0("Tamoxifen_", c(5,10,15,20), "_Low")

## graphs
gA_Tam_low <- func_graphA(model, "Tamoxifen", 0.069, "deeppink", sec, "Patients of genomic low risk")
gB_Tam_low <- func_graphB(model, "Tamoxifen", 2, "deeppink", sec, "Patients of genomic low risk")


#### Goserelin in patients of high genomic risk ####
df0 <- data.frame(treatment = 0, age = "(44,50]", size = "T1c", grade = "2", Group = "STO5_R1",
                  PRstatus = "Positive", HER2status = "Negative", Ki67status = "Low", surgery = "breast_conserving")
df1 <- data.frame(treatment = 1, age = "(44,50]", size = "T1c", grade = "2", Group = "STO5_R1",
                  PRstatus = "Positive", HER2status = "Negative", Ki67status = "Low", surgery = "breast_conserving")
if(length(levels(STO5_high$grade)) != 2){
  STO5_high$grade <- factor(as.character(STO5_high$grade), levels = c("1", "2", "3"), labels = c("2", "2", "3"))
}

STO5_temp <- subset(STO5_high, treatment == "C" | treatment == "G")
STO5_temp$treatment <- as.character(STO5_temp$treatment)
STO5_temp$treatment <- ifelse(STO5_temp$treatment == "C", 0, 1)

model <- stpm2(Surv(Met_20yr, MetBC20yr) ~ treatment + factor(age) + factor(Group) +
                 factor(size) + factor(grade) + factor(PRstatus) + factor(HER2status) +
                 factor(Ki67status) + factor(surgery),
               data = STO5_temp, df = 2, tvc=list(treatment = 1), stata.stpm2.compatible=TRUE)
pred_temp <- predict(model, type="hr", newdata = cbind.data.frame(df0, Met_20yr = c(5,10,15,20)), var = "treatment", se.fit = TRUE) 
rownames(pred_temp) <- paste0("Goserelin_", c(5,10,15,20), "_High")
pred <- rbind(pred, pred_temp)

## graphs
gA_Gos_high <- func_graphA(model, "Goserelin", 0.069, "dodgerblue3", sec, "Patients of genomic high risk")
gB_Gos_high <- func_graphB(model, "Goserelin", 2, "dodgerblue3", sec, "Patients of genomic high risk")


### Table3 ----
pred$HR <- paste0(formatC(pred[,1], format="f", digits=2), "_(",
                  formatC(pred[,2], format="f", digits=2), "-",
                  formatC(pred[,3], format="f", digits=2), ")")
table3 <- cbind(c(5,10,15,20,5,10,15,20), pred$HR)
rownames(table3) <- rownames(pred)
colnames(table3) <- c("Year", "HR_95CI")

#                   Year HR_95CI           
# Tamoxifen_5_Low   "5"  "0.24_(0.08-0.74)" *
# Tamoxifen_10_Low  "10" "0.24_(0.09-0.60)" *
# Tamoxifen_15_Low  "15" "0.23_(0.07-0.75)" *
# Tamoxifen_20_Low  "20" "0.23_(0.06-0.92)" *
# Goserelin_5_High  "5"  "0.26_(0.11-0.61)" *
# Goserelin_10_High "10" "0.29_(0.06-1.35)"
# Goserelin_15_High "15" "0.33_(0.03-3.12)"
# Goserelin_20_High "20" "0.33_(0.03-3.87)"


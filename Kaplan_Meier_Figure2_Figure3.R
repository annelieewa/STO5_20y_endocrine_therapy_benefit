################################################################
# Endocrine therapy & the 70-gene risk signature (MammaPrint)
# in premenopausal ER-positive patients.
# Script to create Kaplan-Meier graphs
# Author: Annelie Johansson
# Date: 8 December 2021
################################################################

library(survival)
library(survminer)

setwd("/Volumes/Annelie Encrypted 2/STO25/Zoladex_STO5/")

# Load data ----
load("Data_temp/STO5_prepared.RData") # * UPDATE
STO <- STO5
nrow(STO) # 584

#### FIGURE 2 ####

#### Function for KM graphs ####
km_plot <- function(event, data, palette, linetype, legends, title){
  
  km <- survfit(Surv(Met_20yr, MetBC20yr) ~ treatment, data = data)
  ylab = "DRFI, Proportion %"

  if(surv_pvalue(km, data = data)$pval < 0.001) { pval <- "<0.001"
  }else{ pval <- paste0("=", formatC(surv_pvalue(km, data = data)$pval, format = "f", digits = 3)) }
  
  kmplot <- ggsurvplot(km, data = data, conf.int = FALSE,
                       pval = paste0("Log-rank P", pval),
                       pval.coord = c(13, 0.50),
                       risk.table = TRUE,
                       censor.shape = "",
                       title = title,
                       tables.height = 0.20,
                       tables.theme = theme_cleantable(),
                       ggtheme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.key = element_rect(fill = "white"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")),
                       risk.table.y.text.col = T,
                       risk.table.y.text = T,
                       palette = palette, # input 
                       linetype = linetype, # input 
                       xlim = c(0, 20.8),
                       break.x.by = 5,
                       xlab = "Years from randomization",
                       ylab = ylab, # input 
                       legend.title = "", # input 
                       legend.labs = legends, # input 
                       font.main = c(15, "plain", "black"),
                       font.x = c(15, "plain", "black"),
                       font.y = c(15, "plain", "black"),
                       font.tickslab = c(15, "plain", "black"),
                       font.legend = c(15, "plain", "black"),
                       legend = c(0.80,0.95),
                       fontsize = 5,
                       pval.size = 5,
                       axes.offset = FALSE)
  return(kmplot)
}

## A. Goserelin in all patients ---- 
STO_temp <- STO
STO_temp$treatment <- factor(STO_temp$treatment, levels = c("G", "C"))
palette <- c("dodgerblue3", "black")
linetype <- c("solid","dotted")
legends <- c("Goserelin", "Control")
title <- ""

gA_DRFI <- km_plot(event = "DRFI", data = STO_temp, palette, linetype, legends, title = title)
gA_DRFI

## B. Tamoxifen in all patients ---- 
STO_temp <- STO
STO_temp$treatment <- factor(STO_temp$treatment, levels = c("T", "C"))
palette <- c("deeppink", "black")
linetype <- c("twodash", "dotted")
legends <- c("Tamoxifen", "Control")
title <- ""

gB_DRFI <- km_plot(event = "DRFI", data = STO_temp, palette, linetype, legends, title = title)
# gB_DRFI

## C. Gos+Tam in all patients ---- 
STO_temp <- STO
STO_temp$treatment <- factor(STO_temp$treatment, levels = c("G+T", "C"))
palette <- c("darkcyan", "black")
linetype <- c("dotdash", "dotted")
legends <- c("Gos+Tam", "Control")
title <- ""

gC_DRFI <- km_plot(event = "DRFI", data = STO_temp, palette, linetype, legends, title = title)
# gC_DRFI

pdf(file="Figures/Figure2_original.pdf", width = 18, height = 6, onefile = FALSE)
arrange_ggsurvplots(list(gA_DRFI, gB_DRFI, gC_DRFI), ncol = 3, nrow = 1)
dev.off()


#### FIGURE 3 #####

#### Function for KM graphs ####
km_plot <- function(event, data, palette, linetype, legends, title){
  
  km <- survfit(Surv(Met_20yr, MetBC20yr) ~ treatment, data = data)
  ylab = "DRFI, Proportion %"
  
  if(surv_pvalue(km, data = data)$pval < 0.001) { pval <- "<0.001"
  }else{ pval <- paste0("=", formatC(surv_pvalue(km, data = data)$pval, format = "f", digits = 3)) }
  
  kmplot <- ggsurvplot(km, data = data, conf.int = FALSE,
                       pval = paste0("Log-rank P", pval),
                       pval.coord = c(13, 0.60),
                       risk.table = TRUE,
                       censor.shape = "",
                       title = title,
                       tables.height = 0.20,
                       tables.theme = theme_cleantable(),
                       ggtheme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.key = element_rect(fill = "white"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")),
                       risk.table.y.text.col = T,
                       risk.table.y.text = T,
                       palette = palette, # input 
                       linetype = linetype, # input 
                       xlim = c(0, 20.5),
                       break.x.by = 5,
                       xlab = "Years from randomization",
                       ylab = ylab, # input 
                       legend.title = "", # input 
                       legend.labs = legends, # input 
                       font.main = c(15, "plain", "black"),
                       font.x = c(15, "plain", "black"),
                       font.y = c(15, "plain", "black"),
                       font.tickslab = c(15, "plain", "black"),
                       font.legend = c(15, "plain", "black"),
                       legend = c(0.80,0.95),
                       fontsize = 5,
                       pval.size = 5,
                       axes.offset = FALSE)
  return(kmplot)
}


## A. Goserelin in low risk ----
STO_temp <- subset(STO, MammaPrint == "Low")
STO_temp$treatment <- factor(STO_temp$treatment, levels = c("G", "C"))
palette <- c("dodgerblue3",  "black")
linetype <- c("solid", "dotted")
legends <- c("Goserelin", "Control")
title <- "Patients of low genomic risk"

gA_DRFI <- km_plot(event = "DRFI", data = STO_temp, palette, linetype, legends, title = title)
# gA_DRFI

## B. Tamoxifen in low risk ----
STO_temp <- subset(STO, MammaPrint == "Low")
STO_temp$treatment <- factor(STO_temp$treatment, levels = c("T", "C"))
palette <- c("deeppink", "black")
linetype <- c("twodash", "dotted")
legends <- c("Tamoxifen", "Control")
title <- "Patients of low genomic risk"

gB_DRFI <- km_plot(event = "DRFI", data = STO_temp, palette, linetype, legends, title = title)
# gB_DRFI

## C. Gos+Tam in low risk ----
STO_temp <- subset(STO, MammaPrint == "Low")
STO_temp$treatment <- factor(STO_temp$treatment, levels = c("G+T", "C"))
palette <- c("darkcyan", "black")
linetype <- c("dotdash", "dotted")
legends <- c("Gos+Tam", "Control")
title <- "Patients of low genomic risk"

gC_DRFI <- km_plot(event = "DRFI", data = STO_temp, palette, linetype, legends, title = title)
# gC_DRFI

## D. Goserelin in high risk ----
STO_temp <- subset(STO, MammaPrint == "High")
STO_temp$treatment <- factor(STO_temp$treatment, levels = c("G", "C"))
palette <- c("dodgerblue3",  "black")
linetype <- c("solid", "dotted")
legends <- c("Goserelin", "Control")
title <- "Patients of high genomic risk"

gD_DRFI <- km_plot(event = "DRFI", data = STO_temp, palette, linetype, legends, title = title)
# gD_DRFI

## E. Tamoxifen in high risk ----
STO_temp <- subset(STO, MammaPrint == "High")
STO_temp$treatment <- factor(STO_temp$treatment, levels = c("T", "C"))
palette <- c("deeppink", "black")
linetype <- c("twodash", "dotted")
legends <- c("Tamoxifen", "Control")
title <- "Patients of high genomic risk"

gE_DRFI <- km_plot(event = "DRFI", data = STO_temp, palette, linetype, legends, title = title)
# gE_DRFI

## F. Goserelin + Tamoxifen in high risk ----
STO_temp <- subset(STO, MammaPrint == "High")
STO_temp$treatment <- factor(STO_temp$treatment, levels = c("G+T", "C"))
palette <- c("darkcyan", "black")
linetype <- c("dotdash", "dotted")
legends <- c("Gos+Tam", "Control")
title <- "Patients of high genomic risk"

gF_DRFI <- km_plot(event = "DRFI", data = STO_temp, palette, linetype, legends, title = title)
# gF_DRFI


pdf(file="Figures/Figure3_original.pdf", width = 18, height = 12, onefile = FALSE)
arrange_ggsurvplots(list(gA_DRFI, gD_DRFI,
                         gB_DRFI, gE_DRFI,
                         gC_DRFI, gF_DRFI), ncol = 3, nrow = 2)
dev.off()










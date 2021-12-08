################################################################
# 20-Year Endocrine Therapy Benefit in Premenopausal Breast Cancer
# Time-varying analyses
# Figure 4 and Supplementary Table S3
# Author: Annelie Johansson
# Date: 8 December 2021
################################################################

library(survival)
library(survminer)
library(rstpm2) 
library(dplyr)

setwd("/Volumes/Annelie Encrypted 2/STO25/Zoladex_STO5/")

# * Load data ----
load("Data_temp/STO5_prepared.RData") # * UPDATE!
nrow(STO5) # 584

# * Exclude patients with missing data ----
STO_use <- subset(STO5, size != "Unknown" & grade != "98" & PRstatus != "Unknown"
                  & HER2status != "Unknown" & Ki67status != "Unknown")
nrow(STO_use) # 563
STO_use$size <- droplevels(STO_use$size)
STO_use$grade <- droplevels(STO_use$grade)
STO_use$PRstatus <- droplevels(STO_use$PRstatus)
STO_use$HER2status <- droplevels(STO_use$HER2status)
STO_use$Ki67status <- droplevels(STO_use$Ki67status)

# * Select low-Risk patients ---- 
STO_low <- subset(STO_use, MammaPrint == "Low")
nrow(STO_low) # 295

# * Select high-Risk patients ----
STO_high <- subset(STO_use, MammaPrint == "High")
nrow(STO_high) # 154

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


#### Multiple plot function ----
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout_ If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}

#########################
# Time-varying analysis #
#########################

sec <- seq(0.05, 20, 0.1)

#### Tamoxifen in patients of low genomic risk ####
df0 <- data.frame(treatment = 0, age = "(44,50]", size = "T1c", grade = "2", Group = "STO5_R1",
                  PRstatus = "Positive", HER2status = "Negative", Ki67status = "Low", surgery = "breast_conserving")
df1 <- data.frame(treatment = 1, age = "(44,50]", size = "T1c", grade = "2", Group = "STO5_R1",
                  PRstatus = "Positive", HER2status = "Negative", Ki67status = "Low", surgery = "breast_conserving")

STO_temp <- subset(STO_low, treatment == "C" | treatment == "T")
STO_temp$treatment <- as.character(STO_temp$treatment)
STO_temp$treatment <- ifelse(STO_temp$treatment == "C", 0, 1)

model <- stpm2(Surv(Met_20yr, MetBC20yr) ~ treatment + factor(age) + factor(Group) +
                 factor(size) + factor(grade) + factor(PRstatus) + factor(HER2status)
               + factor(Ki67status) + factor(surgery),
               data = STO_temp, df = 2, tvc=list(treatment = 1), stata.stpm2.compatible=TRUE)
pred <- predict(model, type="hr", newdata = cbind.data.frame(df0, Met_20yr = c(5,10,15,20)), var = "treatment", se.fit = TRUE) 
rownames(pred) <- paste0("Tamoxifen_", c(5,10,15,20), "_Low")

# extra for Tamoxifen: which years is tamoxifen treatment significant?
pred_extra_tam <- predict(model, type="hr", newdata = cbind.data.frame(df0, Met_20yr = c(1:20)), var = "treatment", se.fit = TRUE) 
rownames(pred_extra_tam) <- paste0("Tamoxifen_", c(1:20), "_Low")
pred_extra_tam
#                  Estimate     lower     upper
# Tamoxifen_1_Low  0.2545691 0.02747576 2.3586397
# Tamoxifen_2_Low  0.2499455 0.04467098 1.3985083
# Tamoxifen_3_Low  0.2471330 0.05883880 1.0380004
# Tamoxifen_4_Low  0.2450533 0.07071033 0.8492554
# Tamoxifen_5_Low  0.2433598 0.08044658 0.7361905
# Tamoxifen_6_Low  0.2418929 0.08799060 0.6649819
# Tamoxifen_7_Low  0.2405619 0.09318756 0.6210060
# Tamoxifen_8_Low  0.2393064 0.09586020 0.5974072
# Tamoxifen_9_Low  0.2381010 0.09589383 0.5911962
# Tamoxifen_10_Low 0.2369927 0.09364558 0.5997670
# Tamoxifen_11_Low 0.2359809 0.08988737 0.6195199
# Tamoxifen_12_Low 0.2350624 0.08535999 0.6473094
# Tamoxifen_13_Low 0.2342346 0.08063238 0.6804444
# Tamoxifen_14_Low 0.2334951 0.07607640 0.7166476
# Tamoxifen_15_Low 0.2328406 0.07190023 0.7540271
# Tamoxifen_16_Low 0.2322665 0.06819664 0.7910615
# Tamoxifen_17_Low 0.2317675 0.06498545 0.8265880
# Tamoxifen_18_Low 0.2313373 0.06224458 0.8597849
# Tamoxifen_19_Low 0.2309690 0.05993066 0.8901400
# Tamoxifen_20_Low 0.2306543 0.05798362 0.9175245
## years 4-20

## graphs
gA_Tam_low <- func_graphA(model, "Tamoxifen", 0.069, "deeppink", sec, "Patients of genomic low risk")
gB_Tam_low <- func_graphB(model, "Tamoxifen", 2, "deeppink", sec, "Patients of genomic low risk")


#### Goserelin in patients of high genomic risk ####
df0 <- data.frame(treatment = 0, age = "(44,50]", size = "T1c", grade = "2", Group = "STO5_R1",
                  PRstatus = "Positive", HER2status = "Negative", Ki67status = "Low", surgery = "breast_conserving")
df1 <- data.frame(treatment = 1, age = "(44,50]", size = "T1c", grade = "2", Group = "STO5_R1",
                  PRstatus = "Positive", HER2status = "Negative", Ki67status = "Low", surgery = "breast_conserving")
if(length(levels(STO_high$grade)) != 2){
  STO_high$grade <- factor(as.character(STO_high$grade), levels = c("1", "2", "3"), labels = c("2", "2", "3"))
}

# * Goserelin ----
STO_temp <- subset(STO_high, treatment == "C" | treatment == "G")
STO_temp$treatment <- as.character(STO_temp$treatment)
STO_temp$treatment <- ifelse(STO_temp$treatment == "C", 0, 1)

model <- stpm2(Surv(Met_20yr, MetBC20yr) ~ treatment + factor(age) + factor(Group) +
                 factor(size) + factor(grade) + factor(PRstatus) + factor(HER2status)
               + factor(Ki67status) + factor(surgery),
               data = STO_temp, df = 2, tvc=list(treatment = 1), stata.stpm2.compatible=TRUE)
pred_temp <- predict(model, type="hr", newdata = cbind.data.frame(df0, Met_20yr = c(5,10,15,20)), var = "treatment", se.fit = TRUE) 
rownames(pred_temp) <- paste0("Goserelin_", c(5,10,15,20), "_High")
pred <- rbind(pred, pred_temp)

# extra for Goserelin: which years is goserelin treatment significant?
pred_extra_gos <- predict(model, type="hr", newdata = cbind.data.frame(df0, Met_20yr = c(1:20)), var = "treatment", se.fit = TRUE) 
rownames(pred_extra_gos) <- paste0("Goserelin_", c(1:20), "_High")
pred_extra_gos
#                   Estimate      lower     upper
# Goserelin_1_High  0.2273361 0.07087146 0.7292315
# Goserelin_2_High  0.2370949 0.09083389 0.6188660
# Goserelin_3_High  0.2445253 0.10315988 0.5796110
# Goserelin_4_High  0.2515241 0.10946144 0.5779604
# Goserelin_5_High  0.2582809 0.10949096 0.6092649
# Goserelin_6_High  0.2651027 0.10436596 0.6733941
# Goserelin_7_High  0.2721377 0.09573466 0.7735856
# Goserelin_8_High  0.2794296 0.08530731 0.9152895
# Goserelin_9_High  0.2869323 0.07449429 1.1051873
# Goserelin_10_High 0.2945157 0.06428766 1.3492406
# Goserelin_11_High 0.3019749 0.05528479 1.6494384
# Goserelin_12_High 0.3090515 0.04776457 1.9996590
# Goserelin_13_High 0.3154692 0.04177663 2.3822130
# Goserelin_14_High 0.3209779 0.03722638 2.7675747
# Goserelin_15_High 0.3253951 0.03394530 3.1191942
# Goserelin_16_High 0.3286336 0.03174117 3.4025230
# Goserelin_17_High 0.3307061 0.03042801 3.5942716
# Goserelin_18_High 0.3318160 0.02977269 3.6980825
# Goserelin_19_High 0.3327448 0.02924195 3.7863107
# Goserelin_20_High 0.3336284 0.02874639 3.8720648
# years 1-8

## graphs
gA_Gos_high <- func_graphA(model, "Goserelin", 0.069, "dodgerblue3", sec, "Patients of genomic high risk")
gB_Gos_high <- func_graphB(model, "Goserelin", 2, "dodgerblue3", sec, "Patients of genomic high risk")


### Save graphs ----

pdf(file = "Figures/Suppl_Figure1_original.pdf", width = 10, height = 8)
multiplot(gA_Tam_low, gB_Tam_low,
          gA_Gos_high, gB_Gos_high,
          layout = matrix(c(1,2,3,4), nrow=2, byrow=TRUE))
dev.off()

### Table data ----
pred$HR <- paste0(formatC(pred[,1], format="f", digits=2), "_(",
                  formatC(pred[,2], format="f", digits=2), "-",
                  formatC(pred[,3], format="f", digits=2), ")")
tab_manuscript <- cbind(c(5,10,15,20,5,10,15,20), pred$HR)
rownames(tab_manuscript) <- rownames(pred)
colnames(tab_manuscript) <- c("Year", "HR_95CI")

# write.table(tab_manuscript, file = "Tables/table_time_varying_Table3.txt", sep ="\t")


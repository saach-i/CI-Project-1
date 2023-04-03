# Associate ploidy with outcome
# Survival data from https://xenabrowser.net/datapages/?dataset=TCGA-OV.survival.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
setwd("/Users/saachisachdev/GitHub/CI Project 1/CNV analysis")

surv_data <- read.delim("TCGA-OV.survival.tsv")
head(surv_data)

# ClustNum contains ploidy information
surv_data$sample
rownames(ClustNum) # need to format the sample names and remove the A or B

# matching the TCGA ploidy samples with the TCGA CNV samples
# replace the . with -
surv_data$sample <- gsub(surv_data$sample, pattern="-", replace=".")
surv_data$sample

# function to remove the A or B
substrRight <- function(x, n){
  substr(x, 1, n)}

surv_data$ID <- substrRight(surv_data$sample, n=15) # duplicate samples present 
surv_data$ID # removed the letters at the end

# deleting duplicate samples from surv_data$ID
dim(surv_data)
surv_data <- surv_data[!duplicated(surv_data$ID), ]
dim(surv_data)

rownames(surv_data) <- surv_data$ID

# removing excess columns containing patient IDs
surv_data <- surv_data[,-1] # removes samples column
surv_data <- surv_data[,-2] 
dim(surv_data)

# Identifying shared samples between survival data and ploidy data 
tcga_ploidy_surv <- merge(surv_data, ClustNum, by = "row.names", all = FALSE)
dim(tcga_ploidy_surv)
head(tcga_ploidy_surv)

# Counting the number of deceased vs surviving patients 
library("dplyr")
tcga_ploidy_surv %>%
  count(OS == "1") # 298 patients deceased, 171 patients alive
              # 1 = Deceased, 0 = Survived 

# Survival analysis 
library(survival)
tcga_ploidy_os <- Surv(tcga_ploidy_surv$OS.time, tcga_ploidy_surv$OS)
summary(coxph(Surv(tcga_ploidy_surv$OS.time, tcga_ploidy_surv$OS) ~ ClustNum, data=tcga_ploidy_surv)) # p value is not signficant 


# Week 4.2&5 analysis 
# make a survfit object to use for the ggkm function
survfit_tcga <- survfit(Surv(tcga_ploidy_surv$OS.time, tcga_ploidy_surv$OS) ~ ClustNum, data = tcga_ploidy_surv)

library(devtools)
install_github("michaelway/ggkm")
library(ggkm)

ggkm(survfit_tcga, table = TRUE, risk_title ="Numbers at risk (n=469)", xlabs = "Time-to-event (days)", ylabs=" Overall Survival (OS) Probability", pval = TRUE,
     pval_threshold = 0.001, ystratalabs = c("Diploid", "Polyploid"), ystrataname = "Ploidy", adj_table_title = -0.1,   adj_y_axis_label = -12.5,
 yscale = 100,  shape = ".", legend = TRUE,  point_size = 3)

# need to change to months? publication examples have days


# Overall survival of the cohort 
plot(survfit(tcga_ploidy_os ~ 0), 
     xlab = "time (days)", 
     ylab = "proportion surviving", 
     mark.time = TRUE)

# Generate an object to classify by ClustNum/Ploidy
#ploidy <- rep(NA, nrow(tcga_ploidy_surv))
#ploidy[tcga_ploidy_surv$ClustNum == "1"] <- 0
#ploidy[tcga_ploidy_surv$ClustNum == "2"] <- 1

#ploidy # checking this is filled in 
#survdiff(tcga_ploidy_os ~ ploidy)
#plot(survfit(tcga_ploidy_os ~ ploidy), 
     xlab = "time (days)", 
     ylab = "proportion surviving")


# Reading in the clinical data - this dataset includes information on the progression free survival 
clin_data <- read.csv("/Users/saachisachdev/GitHub/CI Project 1/CNV analysis/TCGA_ov_clinical.csv", row.names = 1)
colnames(clin_data) 
dim(clin_data)

# Identifying shared samples between clin data and ploidy data 
tcga_ploidy_clin <- merge(clin_data, ClustNum, by = "row.names", all = FALSE)
dim(tcga_ploidy_clin) # 470 samples left 
colnames(tcga_ploidy_clin) # contains clustnum

OS_tcga <- survfit(Surv(tcga_ploidy_clin$Overall.survival..days., tcga_ploidy_clin$OS.event) ~ ClustNum, data = tcga_ploidy_clin)
ggkm(OS_tcga, table = TRUE, xlabs = "Time-to-event (days)", ylabs="Overall Survival (OS) probability", pval = TRUE, ystratalabs = c("Diploid", "Polyploid"), ystrataname = "Ploidy", shape = ".", legend = TRUE,  point_size = 3)



ggsurvplot(
  OS_tcga,
  data = tcga_ploidy_clin,
  size = 1,                 # change line size
  pval = TRUE,              # Add p-value
  risk.table = FALSE,        # Add risk table
  legend.labs =
    c("Diploid", "Polyploid"),    # Change legend labels
  risk.table.height = 0.25,
  pval.size	= 7,
  xlab = "Time-to-event (days)",
  ylab = "Survival probablity (OS)",
  legend.title = "Ploidy",
  risk.table.fontsize = 6, fontsize =13, # Useful to change when you have multiple groups
  ggtheme = theme_classic(base_size = 29)      # Change ggplot2 theme
)

#############################
# stage 3 and 4 OS
# stage 3 and 4
tcga_ploidy_clin_stage34 <- tcga_ploidy_clin[!(tcga_ploidy_clin$Stage == "2"),]
tcga_ploidy_clin_stage34 <- tcga_ploidy_clin_stage34 %>% 
  filter(! Row.names == "TCGA.59.2348.01")

tcga_ploidy_clin_pfs_stage <- Surv(tcga_ploidy_clin_stage34$Progression.free.survival..days., tcga_ploidy_clin_stage34$PFS.event)
survfit_tcga_stage <- survfit(Surv(tcga_ploidy_clin_stage34$Progression.free.survival..days., tcga_ploidy_clin_stage34$PFS.event) ~ ClustNum, data=tcga_ploidy_clin_stage34) 
ggkm(survfit_tcga_stage, table = TRUE, xlabs = "Time-to-event", ylabs="Progression Free Survival (PFS) probability", pval = TRUE, ystratalabs = c("Diploid", "Polyploid"), ystrataname = "Ploidy", shape = ".", legend = TRUE,  point_size = 3)

ggsurvfit

tcga_ploidy_clin_os_stage <- Surv(tcga_ploidy_clin_stage34$Overall.survival..days., tcga_ploidy_clin_stage34$OS.event)
survfit_tcga_stage <- survfit(Surv(tcga_ploidy_clin_stage34$Overall.survival..days., tcga_ploidy_clin_stage34$OS.event) ~ ClustNum, data=tcga_ploidy_clin_stage34) 
ggkm(survfit_tcga_stage, table = TRUE, xlabs = "Time-to-event", ylabs="Overall Survival (OS) probability", pval = TRUE, ystratalabs = c("Diploid", "Polyploid"), ystrataname = "Ploidy", shape = ".", legend = TRUE,  point_size = 3)


# pfs
# Removing NAs/Filtering for the samples who have the progression free survival data 
library(tidyr)

tcga_ploidy_clin_pfs <- tcga_ploidy_clin %>% 
  drop_na("Progression.free.survival..days.")

dim(tcga_ploidy_clin)
dim(tcga_ploidy_clin_pfs) # from 470 samples to 396

# removing the extra patient who had an extended PFS in stage 3 and 4 analysis
tcga_ploidy_clin_pfs$Row.names == "TCGA.59.2348.01"

tcga_ploidy_clin_stage34 <- tcga_ploidy_clin_stage34 %>% 
  filter(! Row.names == "TCGA.59.2348.01")

tcga_ploidy_clin_stage34$Progression.free.survival..days.
tcga_ploidy_clin_pfs_stage <- Surv(tcga_ploidy_clin_stage34$Progression.free.survival..days., tcga_ploidy_clin_stage34$PFS.event)
survfit_tcga_stage <- survfit(Surv(tcga_ploidy_clin_stage34$Progression.free.survival..days., tcga_ploidy_clin_stage34$PFS.event) ~ ClustNum, data=tcga_ploidy_clin_stage34) 

ggkm(survfit_tcga_stage, table = TRUE, xlabs = "Time-to-event (days)", ylabs="Progression Free Survival (PFS) probability", pval = TRUE, ystratalabs = c("Diploid", "Polyploid"), ystrataname = "Ploidy", shape = ".", legend = TRUE,  point_size = 3)


#Creating the survival object for progression free survival - use the PFS event (progressed or not) instead of OS 
library(survival)
library(ggplot)
tcga_ploidy_clin_pfs_surv <- Surv(tcga_ploidy_clin_pfs$Progression.free.survival..days., tcga_ploidy_clin_pfs$PFS.event)
summary(coxph(Surv(tcga_ploidy_clin_pfs$Progression.free.survival..days., tcga_ploidy_clin_pfs$PFS.event) ~ ClustNum, data=tcga_ploidy_clin_pfs)) # p value not significant 
library("ggkm")
survfit_tcga_pfs <- survfit(Surv(tcga_ploidy_clin_pfs$Progression.free.survival..days., tcga_ploidy_clin_pfs$PFS.event)~ClustNum, data = tcga_ploidy_clin_pfs)

# do one pfs curve with the patient and one with them removed to show both results are signficant 
ggkm(survfit_tcga_pfs, table = TRUE, xlabs = "Time-to-event (days)", ylabs="Progression Free Survival (PFS) Probability ", risk_title = "Numbers at risk (n=396)", pval = TRUE,
     pval_threshold = 0.001, ystratalabs = c("Diploid", "Polyploid"), ystrataname = "Ploidy", adj_table_title = -0.1,  pvalpos = c(500, 0.25),  adj_y_axis_label = -12.5,
     yscale = 100,  shape = ".", legend = TRUE,  point_size = 3)

# curve without that patient 
colnames(tcga_ploidy_clin_pfs)
#dim(tcga_ploidy_clin_pfs$Row.names[(tcga_ploidy_clin_pfs$Progression.free.survival..days == "5480")]
# remove patient TCGA.59.2348.01 from the data
dim(tcga_ploidy_clin_pfs)

tcga_ploidy_clin_pfs$Row.names == "TCGA.59.2348.01"

tcga_ploidy_clin_pfs_2 <- tcga_ploidy_clin_pfs %>% 
  filter(! Row.names == "TCGA.59.2348.01")

dim(tcga_ploidy_clin_pfs_2) # removed that patient now repeat survival curve 

survfit_tcga_pfs_2 <- survfit(Surv(tcga_ploidy_clin_pfs_2$Progression.free.survival..days., tcga_ploidy_clin_pfs_2$PFS.event)~ClustNum, data = tcga_ploidy_clin_pfs_2)

ggkm(survfit_tcga_pfs_2, table = TRUE,  xlabs = "Time-to-event (days)", ylabs="Progression Free Survival (PFS) Probability ", risk_title = "Numbers at risk (n=395)", pval = TRUE,
     pval_threshold = 0.05, ystratalabs = c("Diploid", "Polyploid"), ystrataname = "Ploidy", adj_table_title = -0.1,  pvalpos = c(500, 0.25),  adj_y_axis_label = -12.5,
     yscale = 100, legend = TRUE,  point_size = 10, size.table = 4,  size.plot = 25, size.summary = 500) 

ggsurvplot(survfit_tcga_pfs_2, risk.table = TRUE,  size =1, xlabs = "Time-to-event (days)", ylabs="Progression Free Survival (PFS) Probability ", risk_title = "Numbers at risk (n=395)", pval = TRUE,
     pval_threshold = 0.05, ystratalabs = c("Diploid", "Polyploid"), ystrataname = "Ploidy", adj_table_title = -0.1,  pvalpos = c(500, 0.25),  adj_y_axis_label = -12.5,
     yscale = 100, legend = TRUE,  point_size = 10, size.table = 4,  size.plot = 25, size.summary = 500) 


ggsurvplot(
  survfit_tcga_pfs_2,
  data = tcga_ploidy_clin_pfs_2,
  size = 1,                 # change line size
  pval = TRUE,              # Add p-value
  risk.table = FALSE,        # Add risk table
  legend.labs =
    c("Diploid", "Polyploid"),    # Change legend labels
  risk.table.height = 0.25,
  pval.size	= 7,
  legend.title = "Ploidy",
  xlab = "Time-to-event (days)",
  ylab = "Survival probablity (PFS)",
  risk.table.fontsize = 6, fontsize =13, # Useful to change when you have multiple groups
  ggtheme = theme_classic(base_size = 29)      # Change ggplot2 theme
)

library(survminer)


# Checking the rest of the clinical data 
#Stage  
fisher.test(tcga_ploidy_clin$ClustNum, tcga_ploidy_clin$Stage) # p value 
table(tcga_ploidy_clin$ClustNum, tcga_ploidy_clin$Stage)

# Box plot facetted by "dose"
p <- ggboxplot(tcga_ploidy_clin, x = "ClustNum", y = "age_at_initial_pathologic_diagnosis",
               color = "ClustNum", palette = "jco",
               add = "jitter",
               facet.by = "Stage", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format")



#tcga_stage_late <- subset(tcga_ploidy_clin, subset = Stage %in% c(3,4))
#tcga_stage_late

#tcga_stage_early <- subset(tcga_ploidy_clin, subset = Stage %in% c(1,2))
#tcga_stage_early

#fisher.test(tcga_stage_late$ClustNum, tcga_stage_late$Stage) 
#table(tcga_stage_late$ClustNum, tcga_stage_late$Stage) 

#fisher.test(tcga_stage_early$ClustNum, tcga_stage_early$Stage) 
#table(tcga_stage_late$ClustNum, tcga_stage_late$Stage) 

table(tcga_ploidy_clin$ClustNum, cut(tcga_ploidy_clin$age_at_initial_pathologic_diagnosis,c(0,30, 60, 90))) # for each stage, more patients are younger (0-60years)

tcga_stage_2 <- subset(tcga_ploidy_clin, subset = Stage %in% 2)
dim(tcga_stage_2)
tcga_stage_2$ClustNum <- as.numeric(tcga_stage_2$ClustNum)
boxplot(tcga_stage_2$age_at_initial_pathologic_diagnosis ~ tcga_stage_2$ClustNum, main = "
        Stage 2 - Ploidy association with age")
t.test(tcga_stage_2$age_at_initial_pathologic_diagnosis ~ tcga_stage_2$ClustNum) # p-value = 0.3061
tcga_stage_2$ClustNum <- as.factor(tcga_stage_2$ClustNum) 

library(ggpubr)


# violin plot
stat.test.stage2 <- tcga_stage_2 %>%
  t_test(age_at_initial_pathologic_diagnosis ~ ClustNum) %>%
  add_significance()

violin_stage2 <- ggplot(tcga_stage_2, aes(x=ClustNum, y=age_at_initial_pathologic_diagnosis))+ 
  labs(y = "Age at initial pathologic diagnosis", x = "Ploidy")  + 
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  theme(axis.text.x = element_text(size = 40))  +
  geom_line(linewidth=1) +
  stat_pvalue_manual(stat.test.stage2, y.position = 110, label = "p.signif", label.size = 11) +
  geom_violin(trim=FALSE, lwd=1, aes(fill = ClustNum)) + geom_boxplot(width=0.1, fill = "white", lwd=1) + theme_classic() 

violin_stage2 + scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + theme_classic(base_size = 28) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))





#age



library(dplyr)
tcga_ploidy_clin%>%
  group_by(ClustNum)%>% 
  summarise(Mean=mean(age_at_initial_pathologic_diagnosis), Max=max(age_at_initial_pathologic_diagnosis), Min=min(age_at_initial_pathologic_diagnosis),)

#



tcga_ploidy_clin <- tcga_ploidy_clin %>% 
  drop_na("Stage")


violin_age <- ggplot(tcga_ploidy_clin, aes(x=ClustNum, y=age_at_initial_pathologic_diagnosis)) + 
  geom_violin(trim=FALSE, lwd=1, aes(fill = ClustNum)) +
  labs(title = "Association between age and ploidy", y = "Age at initial pathologic diagnosis", x = "Ploidy")  + 
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  theme(axis.text.x = element_text(size = 40))  +
  geom_line(size=1) +
  stat_pvalue_manual(stat.test.age, y.position = 100, label = "p.signif", label.size = 7) + 
  geom_boxplot(width=0.1, fill = "white", lwd=1)


bp <- violin_age + scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + theme_classic(base_size = 14)
bp


# Box plot facetted by "stage"

p <- ggviolin(tcga_ploidy_clin, x = "ClustNum", y = "age_at_initial_pathologic_diagnosis",
               color = "ClustNum", palette = "jco",
               facet.by = "Stage", short.panel.labs = FALSE)

# Use only p.format as label. Remove method name.
p + stat_pvalue_manual(stat.test.age, y.position = 100, label = "p.signif", label.size = 7) 

stat.test.age <- tcga_ploidy_clin %>%
  t_test(age_at_initial_pathologic_diagnosis ~ ClustNum) %>%
  add_significance()


violin_age <- ggplot(tcga_ploidy_clin, aes(x=ClustNum, y=age_at_initial_pathologic_diagnosis)) + 
  geom_violin(trim=FALSE, lwd=1, aes(fill = ClustNum)) +
  labs(y = "Age at diagnosis", x = "Ploidy")  + 
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  theme(axis.text.x = element_text(size = 40))  +
  geom_line(size=1) +
  stat_pvalue_manual(stat.test.age, y.position = 100, label = "p", label.size = 13) + 
  geom_boxplot(width=0.1, fill = "white", lwd=1)


bp <- violin_age + scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + theme_classic(base_size = 33)
bp + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


##########################################################################



tcga_stage_3 <- subset(tcga_ploidy_clin, subset = Stage %in% 3)
dim(tcga_stage_3)
boxplot(tcga_stage_3$age_at_initial_pathologic_diagnosis ~ tcga_stage_3$ClustNum, main = "
        Stage 3 - Ploidy association with age")
t.test(tcga_stage_3$age_at_initial_pathologic_diagnosis ~ tcga_stage_3$ClustNum) # p-value = 2.525e-07

stat.test.stage3 <- tcga_stage_3 %>%
  t_test(age_at_initial_pathologic_diagnosis ~ ClustNum) %>%
  add_significance()

violin_age <- ggplot(tcga_stage_3, aes(x=ClustNum, y=age_at_initial_pathologic_diagnosis)) + 
  geom_violin(trim=FALSE, lwd=1, aes(fill = ClustNum)) +
  labs(y = "Age at initial pathologic diagnosis", x = "Ploidy")  + 
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  theme(axis.text.x = element_text(size = 40))  +
  geom_line(linewidth=1) +
  stat_pvalue_manual(stat.test.stage3, y.position = 105, label = "p.signif", label.size = 11) + 
  geom_boxplot(width=0.1, fill = "white", lwd=1)


bp <- violin_age + scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + theme_classic(base_size = 33)
bp + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


tcga_stage_4 <- subset(tcga_ploidy_clin, subset = Stage %in% 4)
dim(tcga_stage_4)
boxplot(tcga_stage_4$age_at_initial_pathologic_diagnosis ~ tcga_stage_4$ClustNum, main = "
        Stage 4 - Ploidy association with age")
t.test(tcga_stage_4$age_at_initial_pathologic_diagnosis ~ tcga_stage_4$ClustNum) # p-value = 0.002106

stat.test.stage4 <- tcga_stage_4 %>%
  t_test(age_at_initial_pathologic_diagnosis ~ ClustNum) %>%
  add_significance()

violin_stage4 <- ggplot(tcga_stage_4, aes(x=ClustNum, y=age_at_initial_pathologic_diagnosis)) + 
  geom_violin(trim=FALSE, lwd=1, aes(fill = ClustNum)) +
  labs(y = "Age at initial pathologic diagnosis", x = "Ploidy")  + 
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  theme(axis.text.x = element_text(size = 40))  +
  geom_line(size=1) +
  stat_pvalue_manual(stat.test.stage4, y.position = 105, label = "p.signif", label.size = 11) + 
  geom_boxplot(width=0.1, fill = "white", lwd=1)


bp <- violin_stage4 + scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + theme_classic(base_size = 33)
bp + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))





#Surgery outcome
summary(coxph(Surv(tcga_ploidy_clin$Progression.free.survival..days., tcga_ploidy_clin$PFS.event) ~ ClustNum+Stage+Surgery.outcome, data=tcga_ploidy_clin))
table(tcga_ploidy_clin$ClustNum,tcga_ploidy_clin$Surgery.outcome)
fisher.test(tcga_ploidy_clin$ClustNum,tcga_ploidy_clin$Surgery.outcome)

summary(coxph(Surv(tcga_ploidy_clin$Progression.free.survival..days., tcga_ploidy_clin$PFS.event) ~ ClustNum+Stage+Surgery.outcome+age_at_initial_pathologic_diagnosis, data=tcga_ploidy_clin))

# Age -  significant. 

table(tcga_ploidy_clin$ClustNum,cut(tcga_ploidy_clin$age_at_initial_pathologic_diagnosis,c(0,60,1000))) 
fisher.test(tcga_ploidy_clin$ClustNum,cut(tcga_ploidy_clin$age_at_initial_pathologic_diagnosis,c(0,60,1000))) # p-value = 2.098e-06

#Statistical test
library(tidyverse)
library(ggplot2)
library(ggprism)
library(patchwork)
library(magrittr)
library(rstatix)
library(ggpubr)
tcga_ploidy_clin$ClustNum <- as.factor(tcga_ploidy_clin$ClustNum) #need to do this for violin plot (apparently) code worked without

stat.test.age <- tcga_ploidy_clin %>%
  t_test(age_at_initial_pathologic_diagnosis ~ ClustNum) %>%
  add_significance()
# Boxplot
ggplot(tcga_ploidy_clin, aes(x = ClustNum, y = age_at_initial_pathologic_diagnosis)) + 
                         labs(y = "Ploidy", x = "Age at initial pathologic diagnosis") +
                           stat_pvalue_manual(stat.test.age, y.position = 100, label = "p.signif", label.size = 7) + 
                      geom_boxplot() 
+ theme_classic()

# Basic violin plot
                         violin_age <- ggplot(tcga_ploidy_clin, aes(x=ClustNum, y=age_at_initial_pathologic_diagnosis)) + 
   geom_violin(trim=FALSE, lwd=1, aes(fill = ClustNum)) +
labs(title = "Association between age and ploidy", y = "Age at initial pathologic diagnosis", x = "Ploidy")  + 
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  theme(axis.text.x = element_text(size = 40))  +
                           geom_line(size=1) +
  stat_pvalue_manual(stat.test.age, y.position = 100, label = "p.signif", label.size = 7) + 
  geom_boxplot(width=0.1, fill = "white", lwd=1) 

violin_age +scale_fill_manual(values=c("#1e9d77", "#7570b3"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + theme_classic(base_size = 14)



# Molecular sub type - significant 
fisher.test(tcga_ploidy_clin$ClustNum,tcga_ploidy_clin$Molecular.subtype) # p-value = 0.0009479
mol_subtype <- table(tcga_ploidy_clin$ClustNum,tcga_ploidy_clin$Molecular.subtype)
mol_subtype <- as.data.frame(mol_subtype)


#non stacked barplot
#ggplot(mol_subtype, aes(x=Var2, y=Freq, fill=factor(Var1,labels=c("Diploid","Polyploid")))) +
  labs(title = "Association of age and ploidy", y = "Number of HGSOC patients", x = "Molecular subtype") +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values=c('#1e9d77','#7570b3'), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) +
  theme_classic()


# Stacked barplot with multiple groups
ggplot(data=mol_subtype, aes(x=Var1, y=Freq, fill=factor(Var2,labels=c("Differentiated","Immunoreactive", "Mesenchymal", "Proliferative")))) +
  geom_bar(stat="identity",  color="black")+
  labs(title = "Association of molecular subtype and ploidy", y = "Number of HGSOC patients", x = "Ploidy") +
  scale_fill_manual(values=c('#1e9d77','#7570b3', "red", "green"), guide = guide_legend(title = "Molecular subtype", labels=c("Differentiated","Immunoreactive", "Mesenchymal", "Proliferative"))) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
theme_classic() 

install.packages("ggrepel")

# Stacked barplot with multiple groups
ggplot(data=mol_subtype, aes(x=Var1, y=Freq, fill=factor(Var2,labels=c("Differentiated","Immunoreactive", "Mesenchymal", "Proliferative")))) +
  geom_bar(position = "fill", stat="identity",  color="black", lwd = 0.85) +
  labs(y = "Proportion", x = "Ploidy") +
  scale_fill_manual(values=c('#1e9d77','#7570b3', "orange", "light blue"), guide = guide_legend(title = "Molecular subtype", labels=c("Differentiated","Immunoreactive", "Mesenchymal", "Proliferative"))) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels=c("Diploid", "Polyploid")) +
  theme_classic(base_size = 32) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

library(ggpubr)

stat.test.stage2 <- tcga_stage_2 %>%
  t_test(age_at_initial_pathologic_diagnosis ~ ClustNum) %>%
  add_significance()


# Primary therapy outcome - group together Remission/Response, Progressive Disease, Stable Disease
table(tcga_ploidy_clin$ClustNum, tcga_ploidy_clin$primary_therapy_outcome_success)

library(stringi)
tcga_ploidy_clin$primary_therapy_outcome_numbered <- stri_replace_all_regex(tcga_ploidy_clin$primary_therapy_outcome_success,
                                  pattern=c('Partial Remission/Response', 'Stable Disease', 'Progressive Disease', "Complete Remission/Response"),
                                  replacement=c('2', '2', '2', '1'),
                                  vectorize=FALSE)

tcga_ploidy_clin$primary_therapy_outcome_numbered
tcga_ploidy_clin$primary_therapy_outcome_numbered <- as.numeric(tcga_ploidy_clin$primary_therapy_outcome_numbered)

# Duplicating the data set so I don't permanently delete patients 
tcga_ploidy_clin_primary <- tcga_ploidy_clin

tcga_ploidy_clin_primary$primary_therapy_outcome_numbered

library(tidyr)
tcga_ploidy_clin_primary <- tcga_ploidy_clin_primary %>% 
  drop_na("primary_therapy_outcome_numbered")

dim(tcga_ploidy_clin_primary)
tcga_ploidy_clin_primary$primary_therapy_outcome_numbered

fisher.test(tcga_ploidy_clin_primary$ClustNum, tcga_ploidy_clin_primary$primary_therapy_outcome_numbered) # not significant 


## Associate KRas amplification with diploid and polyploid (using copy number data)
# Copy number - gene level data 

colnames(tcga_ploidy_cna)
tcga_ploidy_cna[, "KRAS"]

tcga_kras_ploidy <-  as.data.frame(tcga_ploidy_cna$KRAS)
rownames(tcga_kras_ploidy) <- tcga_ploidy_cna$Row.names

head(tcga_kras_ploidy)

#Assigning an amplified or loss value based on positive or negative 
library(data.table)

tcga_kras_ploidy <- data.table(tcga_kras_ploidy)

colnames(tcga_kras_ploidy) # 1 is loss, 2 is amplified  #threshold should be 0.5
tcga_kras_ploidy$Status <- tcga_kras_ploidy[tcga_ploidy_cna$KRAS < 1 , Label := "Unamplified"]
tcga_kras_ploidy$Status <- tcga_kras_ploidy[tcga_ploidy_cna$KRAS > 1 , Label := "Amplified"]
tcga_kras_ploidy

tcga_kras_ploidy <- as.data.frame(tcga_kras_ploidy)
rownames(tcga_kras_ploidy) <- tcga_ploidy_cna$Row.names

tcga_kras_ploidy <- merge(tcga_kras_ploidy, ClustNum, by = "row.names", all = FALSE)
tcga_kras_ploidy
rownames(tcga_kras_ploidy) <- tcga_ploidy_cna$Row.names
tcga_kras_ploidy

tcga_kras_ploidy <- tcga_kras_ploidy[,-1]
tcga_kras_ploidy <- tcga_kras_ploidy[,-1]
tcga_kras_ploidy

library(tidyr)
tcga_kras_ploidy <- tcga_kras_ploidy %>% 
  drop_na("Label")

dim(tcga_kras_ploidy)
tcga_kras_ploidy

class(tcga_kras_ploidy)
class(tcga_kras_ploidy$Label)

table(tcga_kras_ploidy$Label, tcga_kras_ploidy$ClustNum) # table for ploidy vs amplified/loss
fisher.test(tcga_kras_ploidy$Label, tcga_kras_ploidy$ClustNum) # fisher test for ploidy vs amplified/loss # p-value = 0.09

# tcga_kras_ploidy$Label <- as.character(tcga_kras_ploidy$Label) # needs to be a character vector for gsub

# 1 is loss, 2 is amplified - not working anymore
# tcga_kras_ploidy$Status <- gsub(pattern = "Loss", replace= "1", 
                               # gsub(pattern = "Amplified", replace="2", x = tcga_kras_ploidy$Label))

# tcga_kras_ploidy$Status

### Investigating KRAS vs overall survival 
dim(tcga_kras_ploidy)
dim(clin_data)

tcga_kras_clin <- merge(clin_data, tcga_kras_ploidy, by = "row.names", all = FALSE)
dim(tcga_kras_clin)
tcga_kras_clin

library(survival)

Surv(tcga_kras_clin$Overall.survival..days., tcga_kras_clin$OS.event)
summary(coxph(Surv(tcga_kras_clin$Overall.survival..days., tcga_kras_clin$OS.event) ~ Label, data = tcga_kras_clin)) 
summary(coxph(Surv(tcga_kras_clin$Progression.free.survival..days., tcga_kras_clin$PFS.event) ~ Label, data = tcga_kras_clin)) 

survfit_tcga_kras <- survfit(Surv(Overall.survival..days., OS.event) ~Label, data = tcga_kras_clin)
ggkm(survfit_tcga_kras, table = TRUE, xlabs = "Time-to-event", ylabs="Survival", pval=TRUE)






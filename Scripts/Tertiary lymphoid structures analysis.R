# TLS analysis

library(readr)
TLS_TCGA <- read_rds("../TLSMCPclinical_ov_TCGA.RDS")
TLS_TCGA <- data.frame(TLS_TCGA[,-1], row.names=TLS_TCGA[,1])
TLS_ploidy <- merge(ClustNum, TLS_TCGA, by = 0, all = FALSE)

TLS_ploidy_clin <- merge(TLS_ploidy, clin_data, by.x = "Row.names", by.y = 0, all = FALSE)

# split into the stages 
tls_stage_2 <- subset(TLS_ploidy_clin, subset = Stage.y %in% 2)
tls_stage_3 <- subset(TLS_ploidy_clin, subset = Stage.y %in% 3)
tls_stage_4 <- subset(TLS_ploidy_clin, subset = Stage.y %in% 4)

library(ggpubr)
library(rstatix)
library(ggplot2)

# stage 2 - violin plot
tls_stat_stage2 <- tls_stage_2 %>%
  t_test(TLS.x ~ ClustNum) %>%
  add_significance()

violin_stage2 <- ggplot(tls_stage_2, aes(x=ClustNum, y=TLS.x))+ 
  labs(y = "TLS level", x = "Ploidy", title = "Stage 2")  + 
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  theme(axis.text.x = element_text(size = 40))  +
  geom_line(linewidth=1) +
  stat_pvalue_manual(tls_stat_stage2, y.position = 2, label = "p.signif", label.size = 9) +
  geom_violin(trim=FALSE, lwd=1, aes(fill = ClustNum)) + geom_boxplot(width=0.1, fill = "white", lwd=1) + theme_classic() 

violin_stage2 + scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + theme_classic(base_size = 29) +   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# stage 3 - violin plot
tls_stat_stage3 <- tls_stage_3 %>%
  t_test(TLS.x ~ ClustNum) %>%
  add_significance()

violin_stage3 <- ggplot(tls_stage_3, aes(x=ClustNum, y=TLS.x))+ 
  labs(y = "TLS level", x = "Ploidy", title = "Stage 3")  + 
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  theme(axis.text.x = element_text(size = 40))  +
  geom_line(linewidth=1) +
  stat_pvalue_manual(tls_stat_stage3, y.position = 2, label = "p.signif", label.size = 9) +
  geom_violin(trim=FALSE, lwd=1, aes(fill = ClustNum)) + geom_boxplot(width=0.1, fill = "white", lwd=1) + theme_classic() 

violin_stage3 + scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + theme_classic(base_size = 29) +   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# stage 4 - violin plot
tls_stat_stage4 <- tls_stage_4 %>%
  t_test(TLS.x ~ ClustNum) %>%
  add_significance()

violin_stage4 <- ggplot(tls_stage_4, aes(x=ClustNum, y=TLS.x))+ 
  labs(y = "TLS level", x = "Ploidy", title = "Stage 4")  + 
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  theme(axis.text.x = element_text(size = 40))  +
  geom_line(linewidth=1) +
  stat_pvalue_manual(tls_stat_stage4, y.position = 2, label = "p.signif", label.size = 13) +
  geom_violin(trim=FALSE, lwd=1, aes(fill = ClustNum)) + geom_boxplot(width=0.1, fill = "white", lwd=1) + theme_classic() 

violin_stage4 + scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + theme_classic(base_size = 29) +   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# all stages - boxplot
stat.test <- TLS_ploidy_clin %>%
  t_test(TLS.x ~ ClustNum) %>%
  add_significance()
stat.test

ggplot(data = TLS_ploidy_clin, aes(x = ClustNum, y = TLS.x)) + 
  labs(y = "TLS level", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(size=0.8) +
  geom_boxplot(width=0.7, lwd=1, aes(fill = ClustNum))  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat.test, y.position = 1.75, label = "p.signif", label.size = 15) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 29)


# TLS survival analysis 
library(survival)
library(survminer)
library(data.table)

summary(TLS_ploidy_clin$TLS.x) # 3rd quartile 0.9351, median 0.8080, mean 0.8243
TLS_ploidy_clin$Row.names == "TCGA.59.2348.01"
TLS_ploidy_clin <- TLS_ploidy_clin %>% 
  filter(! Row.names == "TCGA.59.2348.01")

# first split into diploid and polyploid 
TLS_ploidy_clin_diploid <- subset(TLS_ploidy_clin, subset = ClustNum %in% 1)
TLS_ploidy_clin_polyploid <- subset(TLS_ploidy_clin, subset = ClustNum %in% 2)

# diploid tumour group analysis 
# using whole group summary statistics to group high and low
TLS_ploidy_clin_diploid <- as.data.table(TLS_ploidy_clin_diploid)
TLS_ploidy_clin_diploid$TLS_level <- TLS_ploidy_clin_diploid[TLS_ploidy_clin_diploid$TLS.x < 0.9351 , TLS_level := "Normal"]
TLS_ploidy_clin_diploid$TLS_level <- TLS_ploidy_clin_diploid[TLS_ploidy_clin_diploid$TLS.x > 0.9351, TLS_level := "High"]
TLS_ploidy_clin_diploid <- as.data.frame(TLS_ploidy_clin_diploid)
TLS_ploidy_clin_diploid$TLS_level

# Diploid - OS 
TLS_survfit <- survfit(Surv(TLS_ploidy_clin_diploid$Overall.survival..days..x, TLS_ploidy_clin_diploid$OS.event.x) ~ TLS_level, data=TLS_ploidy_clin_diploid) # p value is not signficant 

ggsurvplot(
  TLS_survfit,
  data = TLS_ploidy_clin_diploid,
  size = 1,               
  pval = TRUE,             
  legend.labs =
    c("High", "Normal"),    
  risk.table.height = 0.25,
  pval.size	= 7,
  xlab = "Time-to-event (days)",
  ylab = "Survival probablity (OS)",
  legend.title = "TLS level",
  risk.table.fontsize = 6, fontsize =13, 
  ggtheme = theme_classic(base_size = 29)    
)

# Diploid - PFS 
TLS_survfit_PFS <- survfit(Surv(TLS_ploidy_clin_diploid$Progression.free.survival..days..x, TLS_ploidy_clin_diploid$PFS.event.x) ~ TLS_level, data=TLS_ploidy_clin_diploid) # p value is not signficant 
ggkm(TLS_survfit_PFS, table = TRUE, xlabs = "Time-to-event (days)", ylabs="Progression Free Survival Probability ", risk_title = "Numbers at risk", pval = TRUE,
     pval_threshold = 0.05, ystrataname = "TLS level", ystratalabs = c("High", "Normal"), adj_table_title = -0.1,  pvalpos = c(500, 0.25),  adj_y_axis_label = -12.5,
     yscale = 100,  shape = ".", legend = TRUE,  point_size = 3)

ggsurvplot(
  TLS_survfit_PFS,
  data = TLS_ploidy_clin_diploid,
  size = 1,                
  pval = TRUE,             
  legend.labs =
    c("High", "Low"),   
  risk.table.height = 0.25,
  pval.size	= 7,
  xlab = "Time-to-event (days)",
  ylab = "Survival probablity (PFS)",
  legend.title = "TLS level",
  risk.table.fontsize = 6, fontsize =13, 
  ggtheme = theme_classic(base_size = 29)     
)

# polyploid tumour group analysis 
# Polyploid - OS 
# using whole group summary statistics to group high and low
TLS_ploidy_clin_polyploid <- as.data.table(TLS_ploidy_clin_polyploid)
TLS_ploidy_clin_polyploid$TLS_level <- TLS_ploidy_clin_polyploid[TLS_ploidy_clin_polyploid$TLS.x < 0.9351 , TLS_level := "Normal"]
TLS_ploidy_clin_polyploid$TLS_level <- TLS_ploidy_clin_polyploid[TLS_ploidy_clin_polyploid$TLS.x > 0.9351, TLS_level := "High"]
TLS_ploidy_clin_polyploid <- as.data.frame(TLS_ploidy_clin_polyploid)
TLS_ploidy_clin_polyploid$TLS_level

# Polyploid - OS 
TLS_survfit_poly <- survfit(Surv(TLS_ploidy_clin_polyploid$Overall.survival..days..x, TLS_ploidy_clin_polyploid$OS.event.x) ~ TLS_level, data=TLS_ploidy_clin_polyploid) # p value is not signficant 
ggkm(TLS_survfit_poly, table = TRUE, xlabs = "Time-to-event (days)", ylabs="Overall Survival Probability (Polyploid)", risk_title = "Numbers at risk", pval = TRUE,
     pval_threshold = 0.05, ystrataname = "TLS level", ystratalabs = c("High", "Low"), adj_table_title = -0.1,  pvalpos = c(500, 0.25),  adj_y_axis_label = -12.5,
     yscale = 100,  shape = ".", legend = TRUE,  point_size = 3)

ggsurvplot(
  TLS_survfit_poly,
  data = TLS_ploidy_clin_polyploid,
  size = 1,                 # change line size
  pval = TRUE,              # Add p-value
  legend.labs =
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.25,
  pval.size	= 7,
  xlab = "Time-to-event (days)",
  ylab = "Survival probablity (OS)",
  legend.title = "TLS level",
  risk.table.fontsize = 6, fontsize =13, # Useful to change when you have multiple groups
  ggtheme = theme_classic(base_size = 29)      # Change ggplot2 theme
)

# Polyploid - PFS 
TLS_survfit_PFS_poly <- survfit(Surv(TLS_ploidy_clin_polyploid$Progression.free.survival..days..x, TLS_ploidy_clin_polyploid$PFS.event.x) ~ TLS_level, data=TLS_ploidy_clin_polyploid) # p value is not signficant 
ggkm(TLS_survfit_PFS_poly, table = TRUE, xlabs = "Time-to-event (days)", ylabs="Progression Free Survival Probability ", risk_title = "Numbers at risk", pval = TRUE,
     pval_threshold = 0.05, ystrataname = "TLS level", ystratalabs = c("High", "Low"), adj_table_title = -0.1,  pvalpos = c(500, 0.25),  adj_y_axis_label = -12.5,
     yscale = 100,  shape = ".", legend = TRUE,  point_size = 3)

ggsurvplot(
  TLS_survfit_PFS_poly,
  data = TLS_ploidy_clin_polyploid,
  size = 1,                 # change line size
  pval = TRUE,              # Add p-value
  legend.labs =
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.25,
  pval.size	= 7,
  xlab = "Time-to-event (days)",
  ylab = "Survival probablity (PFS)",
  legend.title = "TLS level",
  pval.coord = c(0, 0.04),
  risk.table.fontsize = 6, fontsize =13, # Useful to change when you have multiple groups
  ggtheme = theme_classic(base_size = 29)      # Change ggplot2 theme
)



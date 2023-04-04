# CORE analyis

library(CORE)
CORE_results_tcga_amp1new$coreTable
df <- CORE_results_tcga_amp1new$coreTable
grange <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL, start.field="start", end.field=c("end", "stop"), strand.field="strand", starts.in.df.are.0based=FALSE)
grange <- nonOverlappingGR(gr = grange, by = "score", decreasing = TRUE, verbose = FALSE)

seg_tcga_cnv_del_core <- findOverlaps(CORE_results_tcga_del_1_nonoverlap, seg_tcga_cnv_grange, maxgap=-1, type=c("any"),   select=c("all"))
seg_tcga_cnv_amp_core <- findOverlaps(grange, seg_tcga_cnv_grange, maxgap=-1, type=c("any"),   select=c("all"))

 # after overlap
length(seg_tcga_cnv_amp_core)
length(seg_tcga_cnv_del_core)

seg_tcga_del_core_df <- as.data.frame(seg_tcga_cnv_del_core)
seg_tcga_amp_core_df <- as.data.frame(seg_tcga_cnv_amp_core)

# patient segmentation data
seg_tcga_cnv <- read.delim("../TCGA.OV.sampleMap_SNP6_nocnv_genomicSegment", sep="\t")
rownames(seg_tcga_cnv) <- as.numeric(rownames(seg_tcga_cnv))

# merging segmentation patient data with core results
seg_core_del_id<- merge(seg_tcga_cnv, seg_tcga_del_core_df, by.x = 0, by.y = "subjectHits")
seg_core_del_id <- seg_core_del_id[, -1]
seg_core_amp_id<- merge(seg_tcga_cnv, seg_tcga_amp_core_df, by.x = 0, by.y = "subjectHits")
seg_core_amp_id <- seg_core_amp_id[, -1]

#subsetting copy number values based on deletion and gain thresholds 
seg_core_del_id_filtered <- subset(seg_core_del_id, value < -0.4150375)
seg_core_amp_id_filtered <- subset(seg_core_amp_id, value > 0.5849625) 

#adding the frequency of cores for each patient into a new dataframe
seg_id_freq_del <- as.data.frame(table(seg_core_del_id_filtered$sample))
seg_id_freq_amp <- as.data.frame(table(seg_core_amp_id_filtered$sample))

#Change ID names so that they match (, -), in order to merge with ClustNum
seg_id_freq_del$Var1 <- gsub(seg_id_freq_del$Var1, pattern="-", replace=".")
seg_id_freq_amp$Var1 <- gsub(seg_id_freq_amp$Var1, pattern="-", replace=".")

# merging with ClustNum (contains ploidy classification for each patient)
# deletions first
tcga_ploidy_core_del <- merge(seg_id_freq_del, ClustNum, by.x = "Var1", by.y = 0, all = FALSE)
tcga_ploidy_core_amp <- merge(seg_id_freq_amp, ClustNum, by.x = "Var1", by.y = 0, all = FALSE)

library(ggplot2)
library(ggprism)
library(patchwork)
library(magrittr)
library(rstatix)
library(ggpubr)

# All cohort deletions
# Statistical test
stat.test_all_deletions <- tcga_ploidy_core_del %>%
  t_test(Freq~ClustNum) %>%
  add_significance()

#Boxplot
ggplot(data = tcga_ploidy_core_del, aes(x = ClustNum, y = Freq)) + 
  labs(y = "Mean recurrent deletions per tumour", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(linewidth=0.8) +
  geom_boxplot(width=0.7, lwd=1, aes(fill = ClustNum))  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat.test_all_deletions, y.position = 300, label = "p.signif", label.size = 9) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 29)

# All cohort amplifications
# Statistical test
stat.test_all_amp <- tcga_ploidy_core_amp %>%
  t_test(Freq~ClustNum) %>%
  add_significance()

#Boxplot
 ggplot(data = tcga_ploidy_core_amp, aes(x = ClustNum, y = Freq)) + 
  labs(y = "Mean recurrent gains per tumour", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(size=0.8) +
  geom_boxplot(width=0.7, lwd=1, aes(fill = ClustNum))  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat.test_all_amp, y.position = 120, label = "p.signif", label.size = 9) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 29)


# Stage - early vs late 
#merging with the clinical data 
tcga_core_clin_amp <- merge(tcga_ploidy_core_amp, clin_data, by.x = "Var1",  by.y = 0, all = FALSE)
tcga_del_early_amp <- subset(tcga_core_clin_amp, subset = Stage %in% 2)
tcga_del_amp_stage3 <- subset(tcga_core_clin_amp, subset = Stage %in% 3)
tcga_del_amp_stage4 <- subset(tcga_core_clin_amp, subset = Stage %in% 4)


library(dplyr)

# Statistical test
# Stage 2 amplifications
stat.test_early_amp <- tcga_del_early_amp %>%
  t_test(Freq~ClustNum) %>%
  add_significance()

# Boxplot
tcga_del_early_amp$ClustNum <- as.factor(tcga_del_early_amp$ClustNum)
ggplot(data = tcga_del_early_amp, aes(x = ClustNum, y = Freq)) + 
  labs(y = "Mean recurrent gains per tumour", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(size=0.8) +
  geom_boxplot(width=0.7, lwd=1, aes(fill = ClustNum))  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat.test_early_amp, y.position = 80, label = "p.signif", label.size = 10) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 28)

# Stage 3 amplifications
stat_del_amp_stage3 <- tcga_del_amp_stage3 %>%
  t_test(Freq~ClustNum) %>%
  add_significance()

ggplot(data = tcga_del_amp_stage3, aes(x = ClustNum, y = Freq)) + 
  labs(y = "Mean recurrent gains per tumour", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(size=0.8) +
  geom_boxplot(width=0.7, lwd=1, aes(fill = ClustNum))  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat_del_amp_stage3, y.position = 120, label = "p.signif", label.size = 13) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 28)

# Stage 4 amplifications
stat_del_amp_stage4 <- tcga_del_amp_stage4 %>%
  t_test(Freq~ClustNum) %>%
  add_significance()

ggplot(data = tcga_del_amp_stage4, aes(x = ClustNum, y = Freq)) + 
  labs(y = "Mean recurrent gains per tumour", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(size=0.8) +
  geom_boxplot(width=0.7, lwd=1, aes(fill = ClustNum))  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat_del_amp_stage4, y.position = 110, label = "p.signif", label.size = 13) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 28)

# Deletions for each stage

tcga_core_clin_del <- merge(tcga_ploidy_core_del, clin_data, by.x = "Var1",  by.y = 0, all = FALSE)
tcga_del_early_del <- subset(tcga_core_clin_del, subset = Stage %in% 2)
tcga_del_del_stage3 <- subset(tcga_core_clin_del, subset = Stage %in% 3)
tcga_del_del_stage4 <- subset(tcga_core_clin_del, subset = Stage %in% 4)

# Stage 2 deletions
#Statistical test
stat.test_early_del <- tcga_del_early_del %>%
  t_test(Freq~ClustNum) %>%
  add_significance()

# Boxplot 
bxp <- ggboxplot(tcga_del_early_del, x = "ClustNum", y = "Freq", xlab = "Ploidy", ylab = "Mean recurrent deletions CNAs/COREs per tumour", title = "Deletions - Early HGSOC stage")
stat.test_early_del_1 <- stat.test_early_del %>% add_xy_position(x = "ClustNum")
bxp + 
  stat_pvalue_manual(stat.test_early_del_1, y.position = 270, label = "p", label.size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  scale_fill_brewer(palette="RdBu") + theme_minimal()

tcga_del_early_del$ClustNum <- as.factor(tcga_del_early_del$ClustNum)

ggplot(data = tcga_del_early_del, aes(x = ClustNum, y = Freq)) + 
  labs(y = "Mean recurrent deletions per tumour", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(size=0.8) +
  geom_boxplot(width=0.7, lwd=1, aes(fill = ClustNum))  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat.test_early_del, y.position = 250, label = "p.signif", label.size = 10) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 28)

# Stage 3 deletions
tcga_del_del_stage3%>%
  group_by(ClustNum)%>% 
  summarise(Mean=mean(Freq), Median=median(Freq))

stat_del_stage3 <- tcga_del_del_stage3 %>%
  t_test(Freq~ClustNum) %>%
  add_significance()

ggplot(data = tcga_del_del_stage3, aes(x = ClustNum, y = Freq)) + 
  labs(y = "Mean recurrent deletions per tumour", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(size=0.8) +
  geom_boxplot(width=0.7, lwd=1, aes(fill = ClustNum))  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat_del_stage3, y.position = 270, label = "p.signif", label.size = 13) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 28)


library(dplyr)
tcga_del_del_stage4%>%
  group_by(ClustNum)%>% 
  summarise(Mean=mean(Freq), Median=median(Freq))

stat_del_stage4 <- tcga_del_del_stage4 %>%
  t_test(Freq~ClustNum) %>%
  add_significance()

ggplot(data = tcga_del_del_stage4, aes(x = ClustNum, y = Freq)) + 
  labs(y = "Mean recurrent deletions per tumour", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(size=0.8) +
  geom_boxplot(width=0.7, lwd=1, aes(fill = ClustNum))  +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat_del_stage4, y.position = 260, label = "p.signif", label.size = 13) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 28)














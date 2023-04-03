# Mutational analysis 
mut <- read.delim("/Users/saachisachdev/GitHub/CI Project 1/CNV analysis/mutations_ov_tcga.txt", fill = T)
mut[, 2]<-gsub(mut[, 2], pattern="-", replace=".")
mut_ploidy <- merge(ClustNum, mut, by.x =0, by.y = "SAMPLE_ID", all = FALSE)
rownames(mut_ploidy) <- mut_ploidy$Row.names
mut_ploidy1 <- mut_ploidy

library(dplyr)
library("tidyr")

#all WT are 1, mutations are 2 
colnames(mut_ploidy1)
mut_ploidy1 <- mut_ploidy1 %>%
  mutate(GABRA6 = case_when(GABRA6 == 'WT' ~ '1',
                         TRUE ~ '2'))

mut_ploidy1 <- mut_ploidy1 %>%
  mutate(BRCA1 = case_when(BRCA1 == 'WT' ~ '1',
                            TRUE ~ '2'))

mut_ploidy1 <- mut_ploidy1 %>%
  mutate(BRCA2 = case_when(BRCA2 == 'WT' ~ '1',
                           TRUE ~ '2'))
mut_ploidy1 <- mut_ploidy1 %>%
  mutate(CSMD3 = case_when(CSMD3 == 'WT' ~ '1',
                           TRUE ~ '2'))

mut_ploidy1 <- mut_ploidy1 %>%
  mutate(NF1 = case_when(NF1 == 'WT' ~ '1',
                           TRUE ~ '2'))
mut_ploidy1 <- mut_ploidy1 %>%
  mutate(CDK12 = case_when(CDK12 == 'WT' ~ '1',
                         TRUE ~ '2'))

mut_ploidy1 <- mut_ploidy1 %>%
  mutate(FAT3 = case_when(FAT3 == 'WT' ~ '1',
                           TRUE ~ '2'))

mut_ploidy1 <- mut_ploidy1 %>%
  mutate(RB1 = case_when(RB1 == 'WT' ~ '1',
                           TRUE ~ '2'))
mut_ploidy1 <- mut_ploidy1 %>%
  mutate(TP53 = case_when(TP53 == 'WT' ~ '1',
                           TRUE ~ '2'))


mut_ploidy1 <- subset(mut_ploidy1, select = -c(STUDY_ID, Row.names))
mut_ploidy2 <- mutate_all(mut_ploidy1, function(x) as.numeric(as.character(x)))

BRCA1 <- chisq.test(table(mut_ploidy2$ClustNum, mut_ploidy2$BRCA1))
BRCA1 <- BRCA1$p.value

BRCA2 <- chisq.test(table(mut_ploidy2$ClustNum, mut_ploidy2$BRCA2))
BRCA2 <- BRCA2$p.value

CSMD3 <- chisq.test(table(mut_ploidy2$ClustNum, mut_ploidy2$CSMD3))
CSMD3 <- CSMD3$p.value

NF1 <- chisq.test(table(mut_ploidy2$ClustNum, mut_ploidy2$NF1))
NF1 <- NF1$p.value

CDK12 <- chisq.test(table(mut_ploidy2$ClustNum, mut_ploidy2$CDK12))
CDK12 <- CDK12$p.value

FAT3 <- chisq.test(table(mut_ploidy2$ClustNum, mut_ploidy2$FAT3))
FAT3 <- FAT3$p.value

GABRA6 <- chisq.test(table(mut_ploidy2$ClustNum, mut_ploidy2$GABRA6))
GABRA6 <- GABRA6$p.value

RB1 <- chisq.test(table(mut_ploidy2$ClustNum, mut_ploidy2$RB1))
RB1 <- RB1$p.value

TP53 <- chisq.test(table(mut_ploidy2$ClustNum, mut_ploidy2$TP53))
TP53 <- TP53$p.value

pvalues <- c(BRCA1 , BRCA2 , CSMD3 , NF1,CDK12, FAT3,  GABRA6,RB1,  TP53)
adjusted <- p.adjust(pvalues,method="fdr")
chi_results = data.frame(BRCA1 = BRCA1, BRCA2 = BRCA2, CSMD3 =CSMD3, NF1=NF1,CDK12=CDK12, FAT3=FAT3,  GABRA6=GABRA6,RB1=RB1,  TP53=TP53)
rownames(chi_results) <- "p.value"
chi_results1 <- rbind(chi_results, adjusted)
rownames(chi_results1) <- c("p.value", "FDR adjusted")

write.csv(chi_results1, "/Users/saachisachdev/GitHub/CI Project 1/CNV analysis/chi_results.csv")


# Mutation plot 
library(reshape2)
mut_ploidy1$Row.names <- rownames(mut_ploidy1)
x2 <- melt(mut_ploidy1,id.var="Row.names",measure.var=c("BRCA1","BRCA2","CSMD3","NF1","CDK12","FAT3","GABRA6","RB1","TP53"))
x3 <- merge(ClustNum, x2, by.x = 0, by.y = "Row.names", all = FALSE)
 
 exp.heatmap <- ggplot(x3, aes(x=reorder(Row.names,ClustNum), y=variable, fill = value)) +
   geom_raster(hjust = 0.9, vjust = 0.9) +
   scale_fill_manual(values=c("grey", "red4")) +
   facet_grid(~ ClustNum, switch = "x", scales = "free_x", space = "free_x") +
   scale_x_discrete(labels=c("Diploid","Polyploid")) +
   theme_classic(base_size = 18)
 
 exp.heatmap + geom_hline(yintercept = 0.5 + 0:35, colour = "white", size = 3) 
 
 
 






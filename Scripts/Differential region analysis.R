# Limma - differential region analysis 
library(limma)
Limma_SigCN <- read.csv("/Users/saachisachdev/GitHub/CI Project 1/CNV analysis/limma_signficantCN_TCGA", header = TRUE)
Limma_SigCN[1:3, 1:3]
ploidy_tcga_ovc$ARRAY
ploidy_tcga_ovc$ClustNum <- as.numeric(ploidy_tcga_ovc$ClustNum)

Limma_ClustNum <- as.data.frame(ploidy_tcga_ovc$ClustNum, ploidy_tcga_ovc$ARRAY)
row.names(Limma_SigCN)=Limma_SigCN[,1]
Limma_SigCN_t <- t(Limma_SigCN[,10:421])

# merge this dataset with ploidy/ClustNum and group into diploid/polyploid
Limma_Combined <- merge(Limma_ClustNum, Limma_SigCN_t, by.x = "row.names", by.y = "row.names", all = FALSE)
rownames(Limma_Combined) <- Limma_Combined[,"Row.names"]
Limma_Combined <- Limma_Combined[,-1]
Limma_Combined[1:4, 1:4]
dim(Limma_Combined)

design_cna <- model.matrix(~ Limma_Combined$`ploidy_tcga_ovc$ClustNum`)
dim(design_cna)
dim(Limma_Combined)

Limma <- Limma_SigCN[,10:421]
dim(Limma)

fit_3 <- lmFit(t(Limma_Combined[,2:100]), design=cbind(intercept=1,group=Limma_Combined$`ploidy_tcga_ovc$ClustNum`))
fit_3 <- eBayes(fit_3)
dge_cna_adjp <- topTable(fit_3
         , coef = 2, n =1000)
         
ID_info <- Limma_SigCN[,1:6]

dge_cna_adjp_df <- merge(dge_cna_adjp, ID_info, by= 0)
rownames(dge_cna_adjp_df) <- dge_cna_adjp_df$Row.names
dge_cna_adjp_df <- dge_cna_adjp_df[, -1]

dge_cna_adjp_df <- dge_cna_adjp_df[order(dge_cna_adjp_df$adj.P.Val), ]

write.csv(dge_cna_adjp_df, "/Users/saachisachdev/GitHub/CI Project 1/CNV analysis/dge_cna_adjp_df.csv", row.names=TRUE)

# volcano plot
library(EnhancedVolcano)
Limma_Combined <- Limma_Combined %>% rename_with(make.names) # remove spaces in colnames 
Limma_Combined %>% group_by(ploidy_tcga_ovc.ClustNum) %>%
  summarize(first=quantile(Deletion.Peak.48...CN.values,probs=0.25),
                           second=quantile(Deletion.Peak.48...CN.values,probs=0.5),
                                          third=quantile(Deletion.Peak.48...CN.values,probs=0.75))

library(tidyverse)
library(ggpubr)
library(rstatix)

dge_cna_adjp_df <- dge_cna_adjp_df[order(dge_cna_adjp_df$Unique.Name), ]
ID_AMP <- dge_cna_adjp_df[1:45, 7]
ID_DEL <- dge_cna_adjp_df[46:99, 7]

# key for deletion and gain colours
keyvals.colour <- ifelse(rownames(dge_cna_adjp_df) %in% ID_DEL, "#1A58A8",
                         ifelse(rownames(dge_cna_adjp_df) %in% ID_AMP, "#E33092", "#8F8F8F"))

keyvals.colour[is.na(keyvals.colour)] <- 'black'
  names(keyvals.colour)[keyvals.colour == "#8F8F8F"] <- 'Not Signficant'
  names(keyvals.colour)[keyvals.colour == "#1A58A8"] <- 'Deletions'
  names(keyvals.colour)[keyvals.colour == "#E33092"] <- 'Gains'


  # most signficant regions and their chromosome descriptor
dge_cna_adjp_df_FC <- dge_cna_adjp_df[order(dge_cna_adjp_df$logFC), ]
FC <- dge_cna_adjp_df_FC$Descriptor[1:3]
FC_2 <- dge_cna_adjp_df_FC$Descriptor[98:99]

dge_cna_adjp_df_P <- dge_cna_adjp_df[order(dge_cna_adjp_df$adj.P.Val), ]
P <- dge_cna_adjp_df_P$Descriptor[1:4]
labels <- c(P, FC_2, FC)

EnhancedVolcano(dge_cna_adjp_df,
                  title = "Differentially expressed regions",
                  subtitle = "HGSOC Diploid and Polyploid",
                  lab = dge_cna_adjp_df$Descriptor,
                  selectLab = labels,
                  x = 'logFC',
                  y = 'adj.P.Val',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  pCutoff = 0.25,
                  FCcutoff = -0.15,
                  xlim = c(-0.2, 0.5),
                  ylim = c(0, 16),
                  pointSize = 2.9,
                  labSize = 3,
                  labCol = 'black',
                  labFace = 'bold',
                  boxedLabels = TRUE,
                  colAlpha = 4/5,
                  legendPosition = 'bottom',
                  legendLabSize = 14,
                  legendIconSize = 4.0,
                  drawConnectors = TRUE,
                arrowheads = FALSE,
                  widthConnectors = .5,
                  colCustom = keyvals.colour,
                  colConnectors = 'black', 
                  caption = bquote(~Log[2]~ "fold change cutoff, -0.15; FDR adjusted p-value cutoff, 0.25"),
  )
  
  
  
# Further analysis of specific regions
# Deletion.Peak.48 analysis
stat.test <- Limma_Combined %>%
  t_test(Deletion.Peak.48...CN.values ~ ploidy_tcga_ovc.ClustNum) %>%
  add_significance()
stat.test

ggplot(data = Limma_Combined, aes(x = ploidy_tcga_ovc.ClustNum, y = Deletion.Peak.48...CN.values)) + 
  labs(y = "Copy number value 
       (17q11.2)", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(linewidth=0.8) + geom_boxplot(width=0.7, lwd=1, aes(fill = ploidy_tcga_ovc.ClustNum)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat.test, y.position = 1, label = "p.signif", label.size = 15) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 29)

# Amplification.Peak.28 analysis
stat.test <- Limma_Combined %>%
  t_test(Amplification.Peak.28...CN.values ~ ploidy_tcga_ovc.ClustNum) %>%
  add_significance()
stat.test

ggplot(data = Limma_Combined, aes(x = ploidy_tcga_ovc.ClustNum, y = Amplification.Peak.28...CN.values)) + 
  labs(y = "Copy number value (12p11.21)", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(linewidth=0.8) + geom_boxplot(width=0.7, lwd=1, aes(fill = ploidy_tcga_ovc.ClustNum)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat.test, y.position = 3.6, label = "p.signif", label.size = 15) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 29)


# Amplification.Peak.19 analysis
stat.test <- Limma_Combined %>%
  t_test(Amplification.Peak.19...CN.values ~ ploidy_tcga_ovc.ClustNum) %>%
  add_significance()
stat.test

ggplot(data = Limma_Combined, aes(x = ploidy_tcga_ovc.ClustNum, y = Amplification.Peak.19...CN.values)) + 
  labs(y = "Copy number value 
       (8q24.12)", x = "Ploidy") +
  scale_x_discrete(labels=c("Diploid","Polyploid")) +
  geom_line(linewidth=0.8) + geom_boxplot(width=0.7, lwd=1, aes(fill = ploidy_tcga_ovc.ClustNum)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  stat_pvalue_manual(stat.test, y.position = 4, label = "p.signif", label.size = 15) +
  scale_fill_manual(values=c("#96c8a2", "#b19cd9"), guide = guide_legend(title = "Ploidy", labels=c("Diploid","Polyploid"))) + 
  theme_classic(base_size = 29)




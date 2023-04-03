# load the ABSOLUTE data file containing absolute ploidy values
ploidy_tcga <- read.delim("/Users/saachisachdev/GitHub/CI Project 1/CNV analysis/TCGA Ploidy data.txt", sep="\t")

#subsetting into a new df  so that only tcga ovarian cancer data remains
ploidy_tcga_ovc<-subset(ploidy_tcga, COHORT=="TCGA Ovarian Cancer") #need to use a double equal sign here

#changing row names to the TCGA sample names 
ploidy_tcga_ovc <- data.frame(ploidy_tcga_ovc, row.names = 2)
ploidy_tcga_ovc$ploidy<-as.numeric(ploidy_tcga_ovc$ploidy) #changing the ploidy values from characters to numeric

#density plot
d <- density(as.numeric(na.omit(ploidy_tcga_ovc$ploidy))) #returns the density data
plot(d, main = "TCGA Ovarian Cancer: genome ploidy density")
abline(v = 2.45)


library(ggpubr)
library(cowplot)

ggplot(ploidy_tcga_ovc, mapping = aes(x = ploidy)) + 
  labs(x = "Ploidy", y = "Density") +
  geom_histogram(aes(y = after_stat(density)),
                 colour = 1, fill = "white") +
  geom_vline(xintercept=2.45,linetype=2) +
  xlim(1, 8) + 
  geom_density(lwd = 0.65, colour = 1,
               fill = 2, alpha = 0.25) + theme_classic(base_size=32)

# Separation based on the absolute ploidy value of 2.45 
#new column overwritting ClustNum
ploidy_tcga_ovc$ClustNum[ploidy_tcga_ovc$ploidy > 2.45] <- "2" # polyploid
ploidy_tcga_ovc$ClustNum[ploidy_tcga_ovc$ploidy < 2.45] <- "1" # diploid
ploidy_tcga_ovc$ClustNum[ploidy_tcga_ovc$ploidy == 2.45] <- "1" # diploid
ploidy_tcga_ovc <- na.omit(ploidy_tcga_ovc)

# sum number of diploid and polyploid samples
table(ploidy_tcga_ovc$ClustNum) # 192 diploid, 278 polyploid 
# creating a new dataframe for subsequent analyses
ClustNum <- data.frame(ClustNum = ploidy_tcga_ovc$ClustNum)
rownames(ploidy_tcga_ovc)<-gsub(rownames(ploidy_tcga_ovc), pattern="-", replace=".")

#changing the row names to remove the extra letters (in order to match with order dataframes)
#creating a function to use
substrRight <- function(x, n){
  substr(x, 1, n)}

rownames(ploidy_tcga_ovc) <- substrRight(rownames(ploidy_tcga_ovc), n=15)
rownames(ClustNum) <- rownames(ploidy_tcga_ovc)



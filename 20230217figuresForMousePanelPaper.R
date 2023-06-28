### script is used to for figures and results section found in mouse panel paper

### figure 1 graphic is just outline of the process
### made specifically with powerpoint and illustrator


### Figure 2 should contain conversion stats and genotyping analysis
### genotyping script here /mnt/DATA6/kevin_recovery/scripts/20190305snpMousePanel.R
### and /mnt/DATA6/kevin_recovery/scripts/20190624updatedSnpReduction.R

library(reshape2)
library(gridExtra)
#library(tidyverse)
library(FactoMineR)
library(VariantAnnotation)
library(DescTools)
snpMcaTable <- read.table(file = "/mnt/DATA6/kevin_recovery/ComprehensiveMouse/20190423_SnpTableMca.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
snpMcaTable.red <- snpMcaTable[,c(1,2,4,5,10,13)]


missingData <- union(which(is.na(snpMcaTable.red$X129S1.SvImJ)),union(which(is.na(snpMcaTable.red$FVB.NJ)),union(which(is.na(snpMcaTable.red$C3H.HeJ)),
                                                                                                                 union(which(is.na(snpMcaTable.red$C57BL.6J)), which(is.na(snpMcaTable.red$BALB.cJ))))))


snpMcaTable.red2 <- snpMcaTable.red[-missingData,]
snpMcaTable.red2 <- snpMcaTable.red2[-which(duplicated(snpMcaTable.red2$RS)),]
snpMcaTable.red2 <- snpMcaTable.red2[-which(is.na(snpMcaTable.red2$RS)),]
snpMcaTableMat <- snpMcaTable.red2[,2:ncol(snpMcaTable.red2)]
rownames(snpMcaTableMat) <- snpMcaTable.red2$RS

snp_mca <- MCA(t(snpMcaTableMat), graph = FALSE, ncp = 5)
summary(snp_mca)


snp_contrib_df <- data.frame(snp_mca$var$contrib)

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20230228mca1v2.pdf", useDingbats = FALSE)
plot.MCA(snp_mca, axes = c(1,2), invisible = "var", xlim = c(-1.5,1.5), ylim = c(-1.5,1.5))
dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20230228mca1v3.pdf", useDingbats = FALSE)
plot.MCA(snp_mca, axes = c(1,3), invisible = "var", xlim = c(-1.5,1.5), ylim = c(-1.5,1.5))
dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20230228mca1v4.pdf", useDingbats = FALSE)
plot.MCA(snp_mca, axes = c(1,4), invisible = "var", xlim = c(-1.5,1.5), ylim = c(-1.5,1.5))
dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20230228mca2v3.pdf", useDingbats = FALSE)
plot.MCA(snp_mca, axes = c(2,3), invisible = "var", xlim = c(-1.5,1.5), ylim = c(-1.5,1.5))
dev.off()



### from the look of the things Dim 1 + 2 separates into 3 groups 129, C57Bl6, and (FVB + BALB + C3H); Dim 3 specifically separates the 3 group cluster
### seems like 2 is the least important .. but we'll see if this changes with different combinations

### will next look at contribution - want at most 300

snp_contrib_df <- snp_mca$var$contrib
dim_1_100 <- sort(snp_contrib_df[,1], decreasing = TRUE)[1:100]
dim_2_100 <- sort(snp_contrib_df[,2], decreasing = TRUE)[1:100]
dim_3_100 <- sort(snp_contrib_df[,3], decreasing = TRUE)[1:100]

dim_1_50 <- sort(snp_contrib_df[,1], decreasing = TRUE)[1:50]
dim_2_50 <- sort(snp_contrib_df[,2], decreasing = TRUE)[1:50]
dim_3_50 <- sort(snp_contrib_df[,3], decreasing = TRUE)[1:50]

testList <- c(names(dim_1_100), names(dim_2_100), names(dim_3_100))
testList <- str_remove(testList, "_.*")

testList2 <- c(names(dim_1_50), names(dim_2_50), names(dim_3_50))
testList2 <- str_remove(testList2, "_.*")

testList3 <- c(names(dim_1_50), names(dim_2_50), names(dim_3_100))
testList3 <- str_remove(testList3, "_.*")


### making graphs for admixture results
### labels are taken from pop files

sampleDf <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/vcfs/20190624/1633/plink.fam", stringsAsFactors = FALSE)
sampleNames <- sampleDf$V2
sampleNames <- str_remove(sampleNames, "_sorted.bam")
admixture1000 <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/vcfs/20190624/1633/plink.5.Q")
admixture500 <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/vcfs/20190624/300/myplinkBed.5.Q")

colnames(admixture1000) <- c("129", "BALB", "C3H", "FVB001", "C57BL/J")
colnames(admixture500) <- c("129", "BALB", "C3H", "FVB001", "C57BL/J")


admixture1000_small <- admixture1000
admixture1000_small <- data.frame(apply(admixture1000_small, 2, function(x) signif(x, 2)), check.names = FALSE)
admixture1000_small$sample <- sampleNames
admixture1000_small_melt  <- reshape2::melt(admixture1000_small, id.vars = "sample", variable.name = "strain")
admixture1000_small_melt$strain<- as.character(admixture1000_small_melt$strain)
colnames(admixture1000_small_melt)[2:3] <- c("strain", "fraction")
fulladmixture1000 <- admixture1000_small_melt

admixture500_small <- admixture500
admixture500_small <- data.frame(apply(admixture500_small, 2, function(x) signif(x, 2)), check.names = FALSE)
admixture500_small$sample <- sampleNames
admixture500_small_melt  <- reshape2::melt(admixture500_small, id.vars = "sample", variable.name = "strain")
admixture500_small_melt$strain<- as.character(admixture500_small_melt$strain)
colnames(admixture500_small_melt)[2:3] <- c("strain", "fraction")
fulladmixture500 <- admixture500_small_melt

colPalette <- c("#009E73","#D55E00", "#8b0000", "#800080", "#ffff00")

keepSamples <- c("129_6NJ_5050", "129_BALB_5050", "129_C3H_5050", "129_FVB_5050",
                 "129_SRA_5050", "6NJ_FVB_5050", "BALB_6NJ_5050",
                 "BALB_C3H_5050", "BALB_FVB_5050", "BALB_SRA_5050","C3H_6NJ_5050",
                 "C3H_FVB_5050", "C3H_SRA_5050", "FVB_SRA_5050",
                 "129S1_SvImJ", "BALB_cJ", "C3H_HeJ", "FVB_NJ", "SRA_C57BL_6J_v3")

fulladmixture1000 <- fulladmixture1000[which(fulladmixture1000$sample %in% keepSamples),]
fulladmixture500 <- fulladmixture500[which(fulladmixture500$sample %in% keepSamples),]


a <- ggplot(fulladmixture1000, aes(x = sample, y = fraction, fill=strain)) +
  geom_bar(position="fill", stat="identity", color = "black") + 
  scale_fill_manual(values=colPalette) + theme_bw() + ggtitle("1000 SNPs") + 
  theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5))


b <- ggplot(fulladmixture500, aes(x = sample, y = fraction, fill=strain)) +
  geom_bar(position="fill", stat="identity", color = "black") + 
  scale_fill_manual(values=colPalette) + theme_bw() +  ggtitle("500 SNPs") +
  theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5))

grid.arrange(a,b, nrow = 2)


# colLabels <- c("129S1_SvImJ", "BALB", "C3H_HeJ", "FVB", "C57Bl_6J")

### structure plot code is in brahm. make these two and make sure to color code the strains correctly between graphs
### figure two should also have the results for 

### not structure plot
# makeAdmixturePlot <- function(x, y){
#   rownames(x) <- sampleNames
#   colnames(x) <- colLabels
#   test <- melt(as.matrix(x))
#   test$Var1 <- as.character(test$Var1)
#   test$Var2 <- as.character(test$Var2)
#   
#   graphVar <- ggplot(test, aes(Var2, Var1)) +
#     geom_tile(aes(fill = value)) + 
#     geom_text(aes(label = round(value, 3))) +
#     scale_fill_gradient(low = "white",high = "red") + 
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#     ylab("Samples") + xlab("Main Strains") + ggtitle(y)
#   return(graphVar)
# }
# 
# graph1600 <- makeAdmixturePlot(admixture1600, "Admixture 1600")
# graph300 <- makeAdmixturePlot(admixture300, "Admixture 300")
# graph150 <- makeAdmixturePlot(admixture150, "Admixture 150")
# graph200 <- makeAdmixturePlot(admixture200, "Admixture 200")
# 
# 
# pdf("/mnt/DATA4/kevhu/ComprehensiveMouse/20190625admixtureResults.pdf",useDingbats = TRUE, height = 7, width = 21)
# grid.arrange(graph1600, graph300, graph150, graph200, nrow = 1, ncol=4)
# dev.off()


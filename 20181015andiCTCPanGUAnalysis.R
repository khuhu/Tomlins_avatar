### best way to start looking at Andi's data is to look at the combined coverage file 
### first I ran Dan's pipeline on some of Andi's dataset with WBC normals and CTC plus PanGU normals
### using this to see what is normally expected 
library(stringr)
library(reshape2)
source("/mnt/DATA4/kevhu/scripts/fastReadFile.R")

### trying to perform poisson regression. idea of x being samples and y being the amplicons is the thinking that the signal from these samples are independent of sampling method
### so the signal/amplicon is dependent on these samples
### overall looking at the modified z-score with log values biases reducing the lowly expressed amplicoms which is what we're looking for


panGuAnno <- read.table("/home/kevhu/data/bedFiles/WG_IAD127899.20170720.designed.forscript.GC.bed", sep = "\t",
                        stringsAsFactors = FALSE)
panGuAnno$V9 <- paste0("AMP_", seq_along(rownames(panGuAnno)))

allRegSamps <- faster.readfile("/mnt/DATA4/kevhu/PanGuPanelAnalysis/all_regCombinedCov.txt",0)
allRegSamps <- allRegSamps[,-c(246,247)]
allRegSamps2 <- data.frame(lapply(allRegSamps[,2:ncol(allRegSamps)], as.numeric), stringsAsFactors = FALSE)
rownames(allRegSamps2) <- allRegSamps$AmpliconId
allRegSamps2 <- data.frame(t(allRegSamps2),stringsAsFactors = FALSE)

a <- sample(1:2900, size = 1)
boxplot(log10((allRegSamps2[,a:(a+200)] +1)), outline = FALSE, xaxt="n", main = "Count distribution amplicons 1173-1374")


medExp <- apply(log10(allRegSamps2 + 1), 2,median)
med <- median(medExp)
modifiedZ <- 0.6745*(medExp - med)/mad(medExp)

hist(modifiedZ)
boxplot(modifiedZ)
summary(modifiedZ)
length(which(abs(modifiedZ) > 2.58))

badAmps <- names(modifiedZ)[which(abs(modifiedZ) > 2.58)]
goodAmps <- names(modifiedZ)[-which(colnames(allRegSamps2) %in% badAmps)]

pdf(file = "/mnt/DATA4/kevhu/PanGuPanelAnalysis/20181016goodVsBadAmps.pdf", useDingbats = FALSE, onefile = TRUE)
boxplot((allRegSamps2[,c(which(colnames(allRegSamps2) %in% goodAmps[300:400]),
                         which(colnames(allRegSamps2) %in% badAmps))]), outline = FALSE, xaxt="n",
        main = "Good versus bad amplicons")
dev.off()


### looking at CTCs ... if t-test doesn't work, maybe then ill use poisson regression on each sample and look at the difference in residuals
### two-sample t for the welch doesn't look so good if our only worry is a difference from 0, we can just just one-sample t-test (Welch)

dataDf <- read.table("/mnt/DATA4/kevhu/PanGuPanelAnalysis/20181107combinedCov.txt", sep = ",",
                     header = FALSE, stringsAsFactors = FALSE)
colnames(dataDf) <- dataDf[1,]
dataDf <- dataDf[-1,]
dataDf <- dataDf[,-which(grepl("Con", colnames(dataDf)))]

andiAnalysis <- read.xlsx2("/mnt/DATA4/kevhu/PanGuPanelAnalysis/AZ_CTC_GoodAmpliconsV1.xlsx", sheetIndex = 1, header = FALSE, stringsAsFactors = FALSE)
colnames(andiAnalysis) <- andiAnalysis[1,]
andiAnalysis <- andiAnalysis[-1,]
andiSamps <- str_replace_all(colnames(andiAnalysis)[7:ncol(andiAnalysis)],"-","_")
which(andiSamps %in% test)
andiSamps[c(1, 14,34,35,36)] <-  c("646_CTC_R8_Reseq", "674_CTC_R5_Reseq","655_CTC_R4_PanGU", "651_CTC_R5_PanGU", "647_CTC_R3_PanGU")
test <- andiSamps[-which(andiSamps %in% colnames(dataDf))]
andiSamps[-which(andiSamps %in% colnames(dataDf))]
dataDf.reduced <- dataDf[,c(1,which(colnames(dataDf) %in% andiSamps))]


ctcAmps <- dataDf$AmpliconId
dataDf2 <- data.frame(lapply(dataDf[,2:ncol(dataDf.reduced)], as.numeric), stringsAsFactors = FALSE)
rownames(dataDf2) <- ctcAmps
dataDf2 <- dataDf2[-which(rownames(dataDf2) %in% badAmps),]
colnames(dataDf2) <- colnames(dataDf.reduced[2:ncol(dataDf.reduced)])


normalWBC <- which(grepl("WBC", colnames(dataDf2)))

#welchResults <- data.frame("AMP","p-value",stringsAsFactors = FALSE)
#for(i in seq_along(rownames(dataDf2))){
#  dummyTest <- t.test(x = dataDf2[i,group1], var.equal = FALSE, alternative = "greater")
#  tmp <- c(rownames(dataDf2)[i],dummyTest$p.value)
#  welchResults <- rbind(welchResults, tmp)
#}
#colnames(welchResults) <- welchResults[1,]
#welchResults <- welchResults[-1,]

#cutoff <- 0.001/nrow(welchResults)
#CTC_badAmps <- welchResults$AMP[(which(welchResults$`p-value` <  cutoff))]

#dfCTC <- dataDf2[,group1]
#dfReg <- dataDf2[,group2]

#par(mfrow=c(2,1))
#boxplot(log2(t(dfCTC[which(rownames(dfCTC) %in% CTC_badAmps[1:100]),]) + 1), las = 2, cex.axis = 0.5)
#boxplot(log2(t(dfReg[which(rownames(dfCTC) %in% CTC_badAmps[1:100]),]) + 1), las = 2, cex.axis = 0.5)
                
#which(colnames(allRegSamps2) %in% badAmps))]), outline = FALSE, xaxt="n",
#        main = "Good versus bad amplicons")



### trying median instead again like previously. doesnt work -> use poisson regression
dfWBC <- log10(dataDf2[,normalWBC]+1)

ctcMedExp <- apply(dfWBC, 1, median)
ctcModZ <- 0.6745*(ctcMedExp - median(ctcMedExp))/mad(ctcMedExp)
length(which(abs(modifiedZ) > 1.9))
hist(ctcModZ)
goodAmps <- names(which(abs(modifiedZ) < 1.9))


boxplot((t(dfWBC[which(rownames(dfWBC) %in% goodAmps),])[,1:100] + 1),
        las = 2, cex.axis = 0.5, ylim=c(0,50))

### one-sided t-test

testResults <- data.frame("AMP","p-value",stringsAsFactors = FALSE)
for(i in seq_along(rownames(dfWBC))){
  dummyTest <- t.test(x = log10(dfWBC[i,] + 1), alternative = "greater")
  tmp <- c(rownames(dfWBC)[i],dummyTest$p.value)
  testResults <- rbind(testResults, tmp)
}
colnames(testResults) <- testResults[1,]
testResults <- testResults[-1,]
testResults$`p-value`[which(testResults$`p-value` == NaN)] <- 0
testResults$`p-value` <- as.numeric(testResults$`p-value`)

cutoff <- 0.01/nrow(testResults)
length(which(abs(testResults$`p-value`) < cutoff))

goodAmpsWbc <- dfWBC[which(abs(testResults$`p-value`) < cutoff & abs(testResults$`p-value`) > 0),]
badAmpsWbc <- dfWBC[which(abs(testResults$`p-value`) > cutoff | abs(testResults$`p-value`) == 0),]


par(mfrow=c(2,1))
boxplot(log2(t(badAmpsWbc[100:200,]) + 1), las = 2, cex.axis = 0.5, main = "Bad WBC Amps", ylim = c(0,12))
boxplot(log2(t(goodAmpsWbc[100:200,]) + 1), las = 2, cex.axis = 0.5, main = "Good WBC Amps", ylim = c(0,12))
summary(apply(goodAmpsWbc, 1, median))

panGuAnno$V8[which(panGuAnno$V9 %in% rownames(badAmpsWbc))]
panGuAnno$V8[which(panGuAnno$V9 %in% rownames(goodAmpsWbc))]

ctcBedFile <- panGuAnno[which(panGuAnno$V9 %in% rownames(goodAmpsWbc)),]
ctcBedFile <- ctcBedFile[,1:8]
andiCtcBed <- panGuAnno[which(panGuAnno$V4 %in% andiAnalysis$region_id),]
andiCtcBed <- andiCtcBed[,1:8]

write.table(ctcBedFile, file = "/mnt/DATA4/kevhu/PanGuPanelAnalysis/20181107ctcBedfileCutoff.bed",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(andiCtcBed, file = "/mnt/DATA4/kevhu/PanGuPanelAnalysis/20181107andiCtcBedfile.bed",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



### making index and normal files

### normal files
allIdx <- read.table("/mnt/DATA4/kevhu/PanGuPanelAnalysis/20181107CtcIdx.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
ctcNormalTxt <- allIdx[which(allIdx$V1 %in% colnames(dfWBC)),]
#ctcNormalTxt$V1[1:10] <- str_replace(ctcNormalTxt$V1[1:10], "_", "-")

write.table(ctcNormalTxt, file = "/home/kevhu/data/normals/PanGU/ampliconCovIdx.PanGU_CTC_normals.n12.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(ctcNormalTxt$V1, file = "/home/kevhu/data/normals/PanGU/ampliconCovIdx.PanGU_CTC_normals.n12.IDlist.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


### idx files
experimentalIdx <- colnames(dataDf2)[-which(grepl("WBC",colnames(dataDf2)))]
experimentalIdx2 <- allIdx[which(allIdx$V1 %in% experimentalIdx),]
#experimentalIdx2$V1 <- str_replace(experimentalIdx2$V1,"_","-")

experimentalIdx2 <- rbind(experimentalIdx2, ctcNormalTxt)

write.table(experimentalIdx2, file = "/mnt/DATA4/kevhu/PanGuPanelAnalysis/20181107testCtcIdx.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



### making heatmaps for the different bed file calls

panGuAnno <- read.table("/home/kevhu/data/bedFiles/WG_IAD127899.20170720.designed.forscript.GC.bed", sep = "\t",
                        stringsAsFactors = FALSE)
panGuAnno$V9 <- paste0("AMP_", seq_along(rownames(panGuAnno)))


andi4pool <- read.table("/mnt/DATA4/kevhu/PanGuPanelAnalysis/20181127cnRatio4PoolsNormal.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(andi4pool) <- str_remove(colnames(andi4pool), "X")
andi4pool2 <- data.frame(lapply(andi4pool[,2:ncol(andi4pool)], function(x) as.numeric(as.character(x))))
rownames(andi4pool2) <- andi4pool$Gene
colnames(andi4pool2) <- colnames(andi4pool)[2:ncol(andi4pool)]
andi4pool2 <- as.matrix(log2(andi4pool2))


andiNoPool <- read.table("/mnt/DATA4/kevhu/PanGuPanelAnalysis/20181205cnGeneRatioNoPool6Amp82Gene.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(andiNoPool) <- str_remove(colnames(andiNoPool), "X")
andiNoPool2 <- data.frame(lapply(andiNoPool[,2:ncol(andiNoPool)], function(x) as.numeric(as.character(x))))
rownames(andiNoPool2) <- andiNoPool$Gene
colnames(andiNoPool2) <- colnames(andiNoPool)[2:ncol(andiNoPool)]
andiNoPool2 <- as.matrix(log2(andiNoPool2))

colnames(andiNoPool2)[37:48]
finalOrder <- c(1,3,4,5,6,7,2,8,9,10,11,12) + 36

andiNoPool3 <- cbind(andiNoPool2[,1:36] , andiNoPool2[,finalOrder])

mybedString <- read.table("/mnt/DATA4/kevhu/PanGuPanelAnalysis/20181108myBedmoreStringent.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(mybedString) <- mybedString[1,]
mybedString <- mybedString[-1,]
mybedString2 <- data.frame(lapply(mybedString[,2:ncol(mybedString)], function(x) as.numeric(as.character(x))))
rownames(mybedString2) <- mybedString$Gene
colnames(mybedString2) <- colnames(mybedString)[2:ncol(mybedString)]
mybedString2  <- as.matrix(log2(mybedString2))
mybedString3  <- mybedString2
mybedString3[which(mybedString3 > -1 & mybedString3 < 0.5849625)] <- 0

mybed <- read.table("/mnt/DATA4/kevhu/PanGuPanelAnalysis/20181108myBedLessStringestGeneEst.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(mybed) <- mybed[1,]
mybed <- mybed[-1,]
mybed2 <- data.frame(lapply(mybed[,2:ncol(mybed)], function(x) as.numeric(as.character(x))))
rownames(mybed2) <- mybed$Gene
colnames(mybed2) <- colnames(mybed)[2:ncol(mybed)]
mybed2 <- as.matrix(log2(mybed2))




heatMapCol <- colorRampPalette(c("blue","white","red"))(1000)
quantile.range <- quantile(as.matrix(andi4pool2), probs = seq(0, 1, 0.01))
colors.breaks2 <- seq(-3,3,6/1000)


sampleNumber <- str_replace(colnames(andi4pool2), "_.*","")
sampleNumber[37:50] <- c( rep("WBC", 10), "Pool","Pool","ArtiPool","ArtiPool")
annoCol4Pools<- data.frame("Sample.Type" = c(rep("CTC", 36), rep("WBC", 14)),
                           "Sample.Number" = sampleNumber)


sampleNumber2 <- str_replace(colnames(andiNoPool2), "_.*","")
sampleNumber2[37:48] <- rep(c("Pool",rep("WBC", 5)),2)
annoColNoPools <- data.frame("Sample.Type" = c(rep("CTC", 36), rep("WBC", 12)),
                              "Sample.Number" = sampleNumber2)


annoTmp <- panGuAnno[match(rownames(andiNoPool3), panGuAnno$V8),c(1,8)]
annoRowNoPools <- data.frame("Chromosome" = annoTmp$V1)




rownames(annoCol4Pools) <- colnames(andi4pool2) 
rownames(annoColNoPools) <- colnames(andiNoPool3)
rownames(annoRowNoPools) <- rownames(andiNoPool3)

anno.colors <- list("Sample.Type"  = c("CTC" = "#8B0000", "WBC" = "#EEC591"), "Sample.Number" = c("646" = "purple4", "647" = "khaki1","650" = "indianred3","651" = "orangered1", "655" = "peru", "663" = "yellow1",
                                                                                                  "669" = "tomato2", "673" = "turquoise1", "674" = "lightgoldenrod1","675" = "steelblue4",
                                                                                                  "676" = "violet", "ArtiPool" = "slateblue4","Pool" = "violetred4", "WBC" = "wheat2"))

anno.colorsNoPool <- list("Sample.Type"  = c("CTC" = "#8B0000", "WBC" = "#EEC591"), "Sample.Number" = c("646" = "purple4", "647" = "khaki1","650" = "indianred3","651" = "orangered1", "655" = "peru", "663" = "yellow1",
                                                                                                  "669" = "tomato2", "673" = "turquoise1", "674" = "lightgoldenrod1","675" = "steelblue4",
                                                                                                  "676" = "violet", "WBC" = "wheat2","Pool" = "violetred4"),
                          "Chromosome" = c("chr1" = "grey0", "chr2" = "gray84", "chr3" = "grey0", "chr4" = "gray84", "chr5" = "grey0",
                                           "chr6" = "gray84", "chr7" = "grey0", "chr8" = "gray84", "chr9" = "grey0", "chr10" = "gray84",
                                           "chr11" = "grey0", "chr12" = "gray84", "chr13" = "grey0", "chr14" = "gray84", "chr15" = "grey0",
                                           "chr16" = "gray84", "chr17" = "grey0", "chr18" = "gray84", "chr19" = "grey0", "chr20" = "gray84",
                                           "chr21" = "grey0", "chr22" = "gray84", "chrX" = "grey0"))


andi4pool2[which(andi4pool2 > 3)] <- 3
andi4pool2[which(andi4pool2 < -3)] <- -3
pheatmap(andi4pool2, color = heatMapCol,
         cluster_cols = TRUE, cluster_rows = TRUE,fontsize = 6,
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         breaks = colors.breaks2, annotation_colors = anno.colors, 
         annotation_col = annoCol4Pools, border_color = "NA",
         cellwidth = 15, cellheight = 6, treeheight_row = 20,
         treeheight_col = 100)

andiNoPool3[which(andiNoPool3 > 3)] <- 3
andiNoPool3[which(andiNoPool3 < -3)] <- -3

### next two lines are only for thresholds
#andiNoPool3[which(andiNoPool3 < log2(1.5) & andiNoPool3  > log2(.75))] <- 0
andiNoPool3[which(andiNoPool3 < log2(1.5) & andiNoPool3  > log2(.5))] <- 0

pheatmap(t(andiNoPool3), color = heatMapCol,
         cluster_cols = FALSE, cluster_rows = TRUE,fontsize = 6,
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         breaks = colors.breaks2, annotation_colors = anno.colorsNoPool, 
         annotation_row = annoColNoPools, border_color = "grey98", annotation_col = annoRowNoPools, 
         cellwidth = 10, cellheight = 4, treeheight_row = 50,
         treeheight_col = 100)



mybed2[which(mybed2 > 3)] <- 3
mybed2[which(mybed2 < -3)] <- -3
pheatmap(mybed2, color = heatMapCol,
         cluster_cols = TRUE, cluster_rows = TRUE,fontsize = 5,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         breaks = colors.breaks2, annotation_colors = anno.colors, 
         annotation_col = annoCol, border_color = "NA",
         cellwidth = 15, cellheight = 5, treeheight_row = 20,
         treeheight_col = 20)


mybedString2[which(mybedString2 > 3)] <- 3
mybedString2[which(mybedString2 < -3)] <- -3
pheatmap(mybedString2, color = heatMapCol,
         cluster_cols = TRUE, cluster_rows = TRUE,fontsize = 5,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         breaks = colors.breaks2, annotation_colors = anno.colors, 
         annotation_col = annoCol, border_color = "NA",
         cellwidth = 15, cellheight = 5, treeheight_row = 20,
         treeheight_col = 20)

boxplot(andi2, las = 2)
boxplot(andi3, las = 2)
boxplot(mybedString2, las = 2)



andi3[which(andi3 > 3)] <- 3
andi3[which(andi3 < -3)] <- -3
pheatmap(andi3, color = heatMapCol,
         cluster_cols = TRUE, cluster_rows = TRUE,fontsize = 5,
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         breaks = colors.breaks2, annotation_colors = anno.colors, 
         annotation_col = annoCol, border_color = "gray95",
         cellwidth = 15, cellheight = 5, treeheight_row = 20,
         treeheight_col = 20)


mybedString3[which(mybedString3 > 3)] <- 3
mybedString3[which(mybedString3 < -3)] <- -3
pheatmap(mybedString3, color = heatMapCol,
         cluster_cols = TRUE, cluster_rows = TRUE,fontsize = 5,
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         breaks = colors.breaks2, annotation_colors = anno.colors, 
         annotation_col = annoCol, border_color = "gray95",
         cellwidth = 15, cellheight = 5, treeheight_row = 20,
         treeheight_col = 20)

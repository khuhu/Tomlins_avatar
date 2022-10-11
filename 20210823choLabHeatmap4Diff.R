firstUpper <- function(gene){
  firstLetter <- toupper(substr(gene, start = 1, stop = 1))
  restOfGene <- tolower(substr(gene, start = 2, stop = nchar(gene)))
  res <- paste0(firstLetter, restOfGene)
  return(res)
}

library(stringr)
library(RTNsurvival)

tcDf <- read.table("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

### KC-01 to 05 is in mixe run with frearon samples

statsTab1 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-138-MG_cho_20210621_354_343/plugin_out/coverageAnalysis_out.668/Auto_MG_cho_20210621_eros_343.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

statsTab2 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/plugin_out/coverageAnalysis_out.681/Auto_MG_cho_202106_3TS_eros_349.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

statsTab3 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-142-MG_cho_20210701_357_353/plugin_out/coverageAnalysis_out.693/Auto_MG_cho_20210701_eros_353.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)


statsTab4 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/plugin_out/coverageAnalysis_out.577/Auto_MG_EFD4_BBN_eros_304.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

combinedCalls1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-138-MG_cho_20210621_354_343/combinedCalls.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)

combinedCalls2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/combinedCalls.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)

combinedCalls3 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-142-MG_cho_20210701_357_353/combinedCalls.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)

combinedCalls4 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/combinedCalls.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)


allStats <- rbind(statsTab1, statsTab2, statsTab3, statsTab4)
allStats$On.Target <- as.numeric(str_remove(allStats$On.Target, "%"))
allStats$Uniformity <- as.numeric(str_remove(allStats$Uniformity, "%"))
badSamps <- allStats$Sample.Name[which(allStats$On.Target < 80 | allStats$Mean.Depth < 100 | allStats$Uniformity < 80)]

cnGeneDf1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-138-MG_cho_20210621_354_343/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
cnGeneDf1 <- cnGeneDf1[,1:49]
if(length(which(colnames(cnGeneDf1) %in% badSamps)) > 0){
  cnGeneDf1 <- cnGeneDf1[,-which(colnames(cnGeneDf1) %in% badSamps)]
}

colnames(cnGeneDf1) <- str_remove(colnames(cnGeneDf1), "_.*")
colnames(cnGeneDf1) <- str_remove(colnames(cnGeneDf1), "[[:punct:]]")
colnames(cnGeneDf1) <- str_replace_all(colnames(cnGeneDf1), " ", "")
colnames(cnGeneDf1) <- str_remove(colnames(cnGeneDf1), "O")


cnGeneDf2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
cnGeneDf2 <- cnGeneDf2[,1:32]
if(length(which(colnames(cnGeneDf2) %in% badSamps)) > 0){
  cnGeneDf2 <- cnGeneDf2[,-which(colnames(cnGeneDf2) %in% badSamps)]
}

colnames(cnGeneDf2) <- str_remove(colnames(cnGeneDf2), "X.*")
colnames(cnGeneDf2) <- str_remove(colnames(cnGeneDf2), "[[:punct:]]")
colnames(cnGeneDf2) <- str_replace_all(colnames(cnGeneDf2), " ", "")
colnames(cnGeneDf2) <- str_remove(colnames(cnGeneDf2), "O")
colnames(cnGeneDf2) <- str_remove(colnames(cnGeneDf2), "_MG.*")

cnGeneDf3 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-142-MG_cho_20210701_357_353/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
cnGeneDf3 <- cnGeneDf3[,1:56]

if(length(which(colnames(cnGeneDf3) %in% badSamps)) > 0){
  cnGeneDf3 <- cnGeneDf3[,-which(colnames(cnGeneDf3) %in% badSamps)]
}

colnames(cnGeneDf3) <- str_remove(colnames(cnGeneDf3), "_X.*")
colnames(cnGeneDf3) <- str_remove(colnames(cnGeneDf3), "[[:punct:]]")
colnames(cnGeneDf3) <- str_replace_all(colnames(cnGeneDf3), " ", "")
colnames(cnGeneDf3) <- str_remove(colnames(cnGeneDf3), "O")
colnames(cnGeneDf3) <- str_remove(colnames(cnGeneDf3), "_MG.*")

cnGeneDf4 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_bigBed/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
cnGeneDf4 <- cnGeneDf4[,1:33]

if(length(which(colnames(cnGeneDf4) %in% badSamps)) > 0){
  cnGeneDf4 <- cnGeneDf4[,-which(colnames(cnGeneDf4) %in% badSamps)]
}

colnames(cnGeneDf4) <- str_remove(colnames(cnGeneDf4), "_X.*")
colnames(cnGeneDf4) <- str_remove(colnames(cnGeneDf4), "[[:punct:]]")
colnames(cnGeneDf4) <- str_replace_all(colnames(cnGeneDf4), " ", "")
colnames(cnGeneDf4) <- str_remove(colnames(cnGeneDf4), "O")
colnames(cnGeneDf4) <- str_remove(colnames(cnGeneDf4), "_MG.*")

combinedCnDf <- cbind(cnGeneDf1, cnGeneDf2[,2:ncol(cnGeneDf2)],
                      cnGeneDf3[,2:ncol(cnGeneDf3)], cnGeneDf4[,2:ncol(cnGeneDf4)])

colnames(combinedCnDf)[which(colnames(combinedCnDf) == "12167")] <- "12167met"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "12401RT2")] <- "12401RT"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "14433MT")] <- "14433met"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "14433LTS")] <- "14433LT"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "2027LTS")] <- "2027LT"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "14109RT")] <- "14109LT"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "2016RT")] <- "2016LT"
# colnames(combinedCnDf)[which(colnames(combinedCnDf) == "KC01_MG")] <- "KC01"
# colnames(combinedCnDf)[which(colnames(combinedCnDf) == "KC02_MG")] <- "KC02"
# colnames(combinedCnDf)[which(colnames(combinedCnDf) == "KC03_MG")] <- "KC03"
# colnames(combinedCnDf)[which(colnames(combinedCnDf) == "KC04_MG")] <- "KC04"
# colnames(combinedCnDf)[which(colnames(combinedCnDf) == "KC05_MG")] <- "KC05"

allSampNames <- colnames(combinedCnDf)

choCombinedCalls <- rbind(combinedCalls1, combinedCalls2, combinedCalls3, combinedCalls4)
choCombinedCalls$Sample <- str_remove(choCombinedCalls$Sample , "_X.*")
choCombinedCalls$Sample <- str_remove(choCombinedCalls$Sample, "[[:punct:]]")
choCombinedCalls$Sample <- str_replace_all(choCombinedCalls$Sample, " ", "")
choCombinedCalls$Sample <- str_remove(choCombinedCalls$Sample, "O")
choCombinedCalls$Sample <- str_remove(choCombinedCalls$Sample, "_.*")
choCombinedCalls$Sample <- str_remove_all(choCombinedCalls$Sample, " ")
choCombinedCalls$Sample <- str_remove(choCombinedCalls$Sample , "X.*")
choCombinedCalls$Sample[which(choCombinedCalls$Sample == "12167MT")] <- "12167met"
choCombinedCalls$Sample[which(choCombinedCalls$Sample == "12401RT2")] <- "12401RT"
choCombinedCalls$Sample[which(choCombinedCalls$Sample == "14433MT")] <- "14433met"
choCombinedCalls$Sample[which(choCombinedCalls$Sample == "14433LTS")] <- "14433LT"
choCombinedCalls$Sample[which(choCombinedCalls$Sample == "2027LTS")] <- "2027LT"
choCombinedCalls$Sample[which(choCombinedCalls$Sample == "14109RT")] <- "14109LT"
choCombinedCalls$Sample[which(choCombinedCalls$Sample == "2016RT")] <- "2016LT"
# choCombinedCalls$Sample[which(choCombinedCalls$Sample == "KC01_MG")] <- "KC01"
# choCombinedCalls$Sample[which(choCombinedCalls$Sample == "KC02_MG")] <- "KC02"
# choCombinedCalls$Sample[which(choCombinedCalls$Sample == "KC03_MG")] <- "KC03"
# choCombinedCalls$Sample[which(choCombinedCalls$Sample == "KC04_MG")] <- "KC04"
# choCombinedCalls$Sample[which(choCombinedCalls$Sample == "KC05_MG")] <- "KC05"
choCombinedCalls <- choCombinedCalls[-grep("Del",choCombinedCalls$Gene),]
choCombinedCalls <- choCombinedCalls[which(choCombinedCalls$Sample %in% allSampNames),]
allSampNames[-which(allSampNames %in% choCombinedCalls$Sample)]

geneNames <- cnGeneDf1$Gene
geneNames <- geneNames[-grep("Del", geneNames)]
choSigMat <- NULL
for (i in unique(choCombinedCalls$Sample)) {
  tmpDf <- choCombinedCalls[which(choCombinedCalls$Sample == i),]
  tmpDf <- tmpDf[match(geneNames,tmpDf$Gene),]
  tmpVector <- ifelse(10^tmpDf$Log10QValue < 0.05, 1, 0)
  choSigMat <- rbind(choSigMat, tmpVector)
}
rownames(choSigMat) <- unique(choCombinedCalls$Sample)
colnames(choSigMat) <- geneNames
choSigMat <- t(choSigMat)

annoTab1 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20210823choLabAging.xlsx")
annoTab1$ID <- str_remove_all(annoTab1$ID, " ")
annoTab1$ID <- str_remove_all(annoTab1$ID, "O")

annoTab1_1 <- annoTab1[which(annoTab1$Study == "Aging (aged)"),]
annoTab1_2 <- annoTab1[which(annoTab1$Study == "Aging (control)"),]


anno1_1Idx <- match(annoTab1_1$ID, colnames(combinedCnDf))
combinedCnAnno1_1 <- combinedCnDf[ , anno1_1Idx]
combinedCnAnno1_1 <- cbind("Gene" = combinedCnDf$Gene, combinedCnAnno1_1)

delRowsAnno1_1 <- combinedCnAnno1_1[grep("Del", combinedCnAnno1_1$Gene), 2:ncol(combinedCnAnno1_1)]
combinedCnAnno1_1 <- combinedCnAnno1_1[-grep("Del", combinedCnAnno1_1$Gene),]
combinedCnDf$Gene <- firstUpper(combinedCnDf$Gene)

### moved the transposition and tc here
combinedCnAnno1_1 <- t(combinedCnAnno1_1[,-1])
matchingTc <- tcDf$tc[match(tolower(rownames(combinedCnAnno1_1)), tcDf$sample)]
for (i in which(!is.na(matchingTc))) {
  tmpVector <- combinedCnAnno1_1[i,]
  tmpVector[tmpVector > 1] <- tmpVector[tmpVector > 1] / matchingTc[i]
  tmpVector[tmpVector < 1] <- tmpVector[tmpVector < 1] * matchingTc[i]
  combinedCnAnno1_1[i,] <- tmpVector
}


#combinedCnAnno1_1[1:nrow(combinedCnAnno1_1),] <- lapply(combinedCnAnno1_1[1:nrow(combinedCnAnno1_1),], function(x) log2(x))
#rownames(combinedCnAnno1_1) <- combinedCnDf$Gene[-grep("del", combinedCnDf$Gene)]
combinedCnAnno1_1 <- log2(combinedCnAnno1_1)
colnames(combinedCnAnno1_1) <- combinedCnDf$Gene[-grep("del", combinedCnDf$Gene)]

#combinedCnAnno1_1 <- t(combinedCnAnno1_1)
#matchingTc <- tcDf$tc[match(tolower(rownames(combinedCnAnno1_1)), tcDf$sample)]
#combinedCnAnno1_1[which(!is.na(matchingTc)),] <- sweep(combinedCnAnno1_1[which(!is.na(matchingTc)),], 1,
#                                                 matchingTc[which(!is.na(matchingTc))], "/")


combinedCnAnno1_1[combinedCnAnno1_1 > log2(2/3) & combinedCnAnno1_1 < log2(4/3)] <- 0
#combinedCnAnno1_1[abs(combinedCnAnno1_1) < 0.2] <- 0
combinedCnAnno1_1 <- t(combinedCnAnno1_1)
combinedCnAnno1_1[combinedCnAnno1_1 < -3] <- -3
combinedCnAnno1_1[combinedCnAnno1_1 > 3] <- 3

# anno1_1TcIdx <- match(colnames(choSigMat), colnames(combinedCnAnno1_1), )
anno1_1TcIdx <- match(colnames(combinedCnAnno1_1), colnames(choSigMat))
#anno1_1TcIdx <- anno1_1TcIdx[-which(is.na(anno1_1TcIdx))]
choSigMatAnno1_1 <- choSigMat[,anno1_1TcIdx]
combinedCnAnno1_1[choSigMatAnno1_1 == 0] <- 0

annoTab1_1graph <- data.frame("Label" = annoTab1_1$Study, "Type" = annoTab1_1$Type,
                             stringsAsFactors = FALSE)
rownames(annoTab1_1graph) <- annoTab1_1$ID
annoTab1_1graph <- annoTab1_1graph[which(rownames(annoTab1_1graph) %in% colnames(combinedCnAnno1_1)),]

annoCol1 <- list("Label" = c("Aging (aged)" = "#556b2f", "Aging (control)" = "black"),
                "Type"  = c("HGSC" = "red", "MMMT" = "dodgerblue"))

heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-3,3,6/1000)

combinedCnAnno1_1 <- t(combinedCnAnno1_1)

### 20211004: adding semi-supervised clustering for the two groups - can't

heatmap_Anno1_1genes <- pheatmap(mat = combinedCnAnno1_1, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                             breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                             border_color = "black", annotation_row = annoTab1_1graph, annotation_colors = annoCol1)

rownames(delRowsAnno1_1) <- combinedCnDf$Gene[grep("del", combinedCnDf$Gene)]
delRowsAnno1_1 <- t(log2(delRowsAnno1_1))
# delRowsAnno1[which(!is.na(matchingTc)),] <- sweep(delRowsAnno1[which(!is.na(matchingTc)),], 1,
#                                                      matchingTc[which(!is.na(matchingTc))], "/")
# delRowsAnno1[delRowsAnno1 > log2(2/3) & delRowsAnno1 < log2(4/3)] <- 0
# delRowsAnno1[delRowsAnno1 < -3] <- -3
# delRowsAnno1[delRowsAnno1 > 3] <- 3

heatmap_Anno1_1dels <- pheatmap(delRowsAnno1_1[heatmap_Anno1_1genes$tree_row$order,], cluster_rows = FALSE, cluster_cols = FALSE,
                               color = heatMapCol,breaks = colors.breaks, fontsize = 5, cellwidth = 5,
                               cellheight = 10,silent = FALSE, border_color = "black")


heatMapCol_tc <- colorRampPalette(c("#FFFFFF","#FFA500"))(100)
colors.breaks_tc <- seq(0,1,1/100)
tc_labels <- matchingTc
tc_labels<- t(tc_labels)
heatmap_Anno1_1tc <- pheatmap(tc_labels[heatmap_Anno1_1genes$tree_row$order], cluster_rows = FALSE, cluster_cols = FALSE,
                             color = heatMapCol_tc, breaks = colors.breaks_tc, fontsize = 5, cellwidth = 5,
                             cellheight = 10,silent = FALSE, border_color = "black")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20211004choAnno1_1Heatmap_gene.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_Anno1_1genes
dev.off()


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20211004choAnno1_1Heatmap_del.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_Anno1_1dels
dev.off()

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20211004choAnno1_1Heatmap_tc.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_Anno1_1tc
dev.off()

###
###
### second split
anno1_2Idx <- match(annoTab1_2$ID, colnames(combinedCnDf))
anno1_2Idx <- anno1_2Idx[-which(is.na(anno1_2Idx))]
combinedCnAnno1_2 <- combinedCnDf[ , anno1_2Idx]
combinedCnAnno1_2 <- cbind("Gene" = combinedCnDf$Gene, combinedCnAnno1_2)

delRowsAnno1_2 <- combinedCnAnno1_2[grep("del", combinedCnAnno1_2$Gene),
                                    2:ncol(combinedCnAnno1_2)]
combinedCnAnno1_2 <- combinedCnAnno1_2[-grep("del", combinedCnAnno1_2$Gene),]
combinedCnDf$Gene <- firstUpper(combinedCnDf$Gene)

combinedCnAnno1_2 <- t(combinedCnAnno1_2[,-1])
matchingTc <- tcDf$tc[match(tolower(rownames(combinedCnAnno1_2)), tcDf$sample)]
for (i in which(!is.na(matchingTc))) {
  tmpVector <- combinedCnAnno1_2[i,]
  tmpVector[tmpVector > 1] <- tmpVector[tmpVector > 1] / matchingTc[i]
  tmpVector[tmpVector < 1] <- tmpVector[tmpVector < 1] * matchingTc[i]
  combinedCnAnno1_2[i,] <- tmpVector
}


combinedCnAnno1_2 <- log2(combinedCnAnno1_2)
colnames(combinedCnAnno1_2) <- combinedCnDf$Gene[-grep("del", combinedCnDf$Gene)]

# combinedCnAnno1_2 <- combinedCnAnno1_2[,-1]
# combinedCnAnno1_2[1:nrow(combinedCnAnno1_2),] <- lapply(combinedCnAnno1_2[1:nrow(combinedCnAnno1_2),], function(x) log2(x))
# rownames(combinedCnAnno1_2) <- combinedCnDf$Gene[-grep("del", combinedCnDf$Gene)]
# 
# combinedCnAnno1_2 <- t(combinedCnAnno1_2)
# matchingTc <- tcDf$tc[match(tolower(rownames(combinedCnAnno1_2)), tcDf$sample)]
# combinedCnAnno1_2[which(!is.na(matchingTc)),] <- sweep(combinedCnAnno1_2[which(!is.na(matchingTc)),], 1,
#                                                        matchingTc[which(!is.na(matchingTc))], "/")
combinedCnAnno1_2[combinedCnAnno1_2 > log2(2/3) & combinedCnAnno1_2 < log2(4/3)] <- 0
#combinedCnAnno1_2[abs(combinedCnAnno1_2) < 0.2] <- 0
combinedCnAnno1_2 <- t(combinedCnAnno1_2)
combinedCnAnno1_2[combinedCnAnno1_2 < -3] <- -3
combinedCnAnno1_2[combinedCnAnno1_2 > 3] <- 3

anno1_2TcIdx <- match(colnames(combinedCnAnno1_2), colnames(choSigMat))
#anno1_2TcIdx <- match(colnames(choSigMat),colnames(combinedCnAnno1_2))
#anno1_2TcIdx <- anno1_2TcIdx[-which(is.na(anno1_2TcIdx))]
choSigMatAnno1_2 <- choSigMat[,anno1_2TcIdx]
combinedCnAnno1_2[choSigMatAnno1_2 == 0] <- 0

annoTab1_2graph <- data.frame("Label" = annoTab1_2$Study, "Type" = annoTab1_2$Type,
                              stringsAsFactors = FALSE)
rownames(annoTab1_2graph) <- annoTab1_2$ID
annoTab1_2graph <- annoTab1_2graph[which(rownames(annoTab1_2graph) %in% colnames(combinedCnAnno1_2)),]

annoCol1 <- list("Label" = c("Aging (aged)" = "#556b2f", "Aging (control)" = "#191970"),
                 "Type"  = c("HGSC" = "red", "MMMT" = "dodgerblue"))

heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-3,3,6/1000)

combinedCnAnno1_2 <- t(combinedCnAnno1_2)

heatmap_Anno1_2genes <- pheatmap(mat = combinedCnAnno1_2, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                                 breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                                 border_color = "black", annotation_row = annoTab1_2graph, annotation_colors = annoCol1)

rownames(delRowsAnno1_2) <- combinedCnDf$Gene[grep("del", combinedCnDf$Gene)]
delRowsAnno1_2 <- t(log2(delRowsAnno1_2))

heatmap_Anno1_2dels <- pheatmap(delRowsAnno1_2[heatmap_Anno1_2genes$tree_row$order,], cluster_rows = FALSE, cluster_cols = FALSE,
                                color = heatMapCol,breaks = colors.breaks, fontsize = 5, cellwidth = 5,
                                cellheight = 10,silent = FALSE, border_color = "black")


heatMapCol_tc <- colorRampPalette(c("#FFFFFF","#FFA500"))(100)
colors.breaks_tc <- seq(0,1,1/100)
tc_labels <- matchingTc
tc_labels<- t(tc_labels)
heatmap_Anno1_2tc <- pheatmap(tc_labels[heatmap_Anno1_2genes$tree_row$order], cluster_rows = FALSE, cluster_cols = FALSE,
                              color = heatMapCol_tc, breaks = colors.breaks_tc, fontsize = 5, cellwidth = 5,
                              cellheight = 10,silent = FALSE, border_color = "black")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20211004choAnno1_2Heatmap_gene.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_Anno1_2genes
dev.off()


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20211004choAnno1_2Heatmap_del.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_Anno1_2dels
dev.off()

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20211004choAnno1_2Heatmap_tc.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_Anno1_2tc
dev.off()


###
###
###

annoTab2 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20210823choLabPreReg.xlsx")
annoTab2$ID <- str_remove_all(annoTab2$ID, " ")
annoTab2$ID <- str_remove_all(annoTab2$ID, "O")


anno2Idx <- match(annoTab2$ID, colnames(combinedCnDf))
combinedCnanno2 <- combinedCnDf[ , anno2Idx]
combinedCnanno2 <- cbind("Gene" = combinedCnDf$Gene, combinedCnanno2)

delRowsanno2 <- combinedCnanno2[grep("del", combinedCnanno2$Gene), 2:ncol(combinedCnanno2)]
combinedCnanno2 <- combinedCnanno2[-grep("del", combinedCnanno2$Gene),]
combinedCnDf$Gene <- firstUpper(combinedCnDf$Gene)
combinedCnanno2 <- combinedCnanno2[,-1]
combinedCnanno2[1:nrow(combinedCnanno2),] <- lapply(combinedCnanno2[1:nrow(combinedCnanno2),], function(x) log2(x))
rownames(combinedCnanno2) <- combinedCnDf$Gene[-grep("del", combinedCnDf$Gene)]

combinedCnanno2 <- t(combinedCnanno2)
matchingTc <- tcDf$tc[match(tolower(rownames(combinedCnanno2)), tcDf$sample)]
combinedCnanno2[which(!is.na(matchingTc)),] <- sweep(combinedCnanno2[which(!is.na(matchingTc)),], 1,
                                                     matchingTc[which(!is.na(matchingTc))], "/")
combinedCnanno2[combinedCnanno2 > log2(2/3) & combinedCnanno2 < log2(4/3)] <- 0
combinedCnanno2 <- t(combinedCnanno2)
combinedCnanno2[combinedCnanno2 < -3] <- -3
combinedCnanno2[combinedCnanno2 > 3] <- 3

anno2TcIdx <- match(colnames(choSigMat),colnames(combinedCnanno2))
anno2TcIdx <- anno2TcIdx[-which(is.na(anno2TcIdx))]
choSigMatanno2 <- choSigMat[,anno2TcIdx]
combinedCnanno2[choSigMatanno2 == 0] <- 0

annoTab2_graph <- data.frame("Label" = annoTab2$Study, "Type" = annoTab2$Type,
                             stringsAsFactors = FALSE)
rownames(annoTab2_graph) <- annoTab2$ID
annoTab2_graph <- annoTab2_graph[which(rownames(annoTab2_graph) %in% colnames(combinedCnanno2)),]

annoCol2 <- list("Label" = c("Pregnancy (high parity)" = "#556b2f", "Pregnancy (control)" = "black"),
                 "Type"  = c("HGSC" = "red", "MMMT" = "dodgerblue"))

heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-3,3,6/1000)

combinedCnanno2 <- t(combinedCnanno2)
heatmap_anno2_genes <- pheatmap(mat = combinedCnanno2, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                                breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                                border_color = "black", annotation_row = annoTab2_graph, annotation_colors = annoCol2)


rownames(delRowsanno2) <- combinedCnDf$Gene[grep("del", combinedCnDf$Gene)]
delRowsanno2 <- t(delRowsanno2)
heatmap_anno2_dels <- pheatmap(log2(delRowsanno2)[heatmap_anno2_genes$tree_row$order,], cluster_rows = FALSE, cluster_cols = FALSE,
                               color = heatMapCol,breaks = colors.breaks, fontsize = 5, cellwidth = 5,
                               cellheight = 10,silent = FALSE, border_color = "black")


heatMapCol_tc <- colorRampPalette(c("#FFFFFF","#FFA500"))(100)
colors.breaks_tc <- seq(0,1,1/100)
tc_labels <- matchingTc
tc_labels<- t(tc_labels)
heatmap_anno2_tc <- pheatmap(tc_labels[heatmap_anno2_genes$tree_row$order], cluster_rows = FALSE, cluster_cols = FALSE,
                             color = heatMapCol_tc, breaks = colors.breaks_tc, fontsize = 5, cellwidth = 5,
                             cellheight = 10,silent = FALSE, border_color = "black")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210823choAnno2Heatmap_gene.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_anno2_genes
dev.off()


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210823choAnno2Heatmap_del.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_anno2_dels
dev.off()

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210823choAnno2Heatmap_tc.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_anno2_tc
dev.off()

###
###
###

annoTab3 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20210823choLabPregNf.xlsx")
annoTab3$ID <- str_remove_all(annoTab3$ID, " ")
annoTab3$ID <- str_remove_all(annoTab3$ID, "O")


anno3Idx <- match(annoTab3$ID, colnames(combinedCnDf))
anno3Idx <- anno3Idx[-which(is.na(anno3Idx))]
combinedCnanno3 <- combinedCnDf[ , anno3Idx]
combinedCnanno3 <- cbind("Gene" = combinedCnDf$Gene, combinedCnanno3)

delRowsanno3 <- combinedCnanno3[grep("del", combinedCnanno3$Gene), 2:ncol(combinedCnanno3)]
combinedCnanno3 <- combinedCnanno3[-grep("del", combinedCnanno3$Gene),]
combinedCnDf$Gene <- firstUpper(combinedCnDf$Gene)
combinedCnanno3 <- combinedCnanno3[,-1]
combinedCnanno3[1:nrow(combinedCnanno3),] <- lapply(combinedCnanno3[1:nrow(combinedCnanno3),], function(x) log2(x))
rownames(combinedCnanno3) <- combinedCnDf$Gene[-grep("del", combinedCnDf$Gene)]

combinedCnanno3 <- t(combinedCnanno3)
matchingTc <- tcDf$tc[match(tolower(rownames(combinedCnanno3)), tcDf$sample)]
combinedCnanno3[which(!is.na(matchingTc)),] <- sweep(combinedCnanno3[which(!is.na(matchingTc)),], 1,
                                                     matchingTc[which(!is.na(matchingTc))], "/")
combinedCnanno3[combinedCnanno3 > log2(2/3) & combinedCnanno3 < log2(4/3)] <- 0
combinedCnanno3 <- t(combinedCnanno3)
combinedCnanno3[combinedCnanno3 < -3] <- -3
combinedCnanno3[combinedCnanno3 > 3] <- 3

anno3TcIdx <- match(colnames(choSigMat),colnames(combinedCnanno3))
anno3TcIdx <- anno3TcIdx[-which(is.na(anno3TcIdx))]
choSigMatanno3 <- choSigMat[,anno3TcIdx]
combinedCnanno3[choSigMatanno3 == 0] <- 0

annoTab3_graph <- data.frame("Label" = annoTab3$Study, "Type" = annoTab3$Type,
                             stringsAsFactors = FALSE)
rownames(annoTab3_graph) <- annoTab3$ID
annoTab3_graph <- annoTab3_graph[which(rownames(annoTab3_graph) %in% colnames(combinedCnanno3)),]

annoCol3 <- list("Label" = c("Pregnancy (high parity_Nf1fl/+)" = "#556b2f", "Pregnancy (control_Nf1fl/+)" = "black"),
                 "Type"  = c("HGSC" = "red", "MMMT" = "dodgerblue", "eHGSC" = "yellow"))

heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-3,3,6/1000)

combinedCnanno3 <- t(combinedCnanno3)
heatmap_anno3_genes <- pheatmap(mat = combinedCnanno3, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                                breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                                border_color = "black", annotation_row = annoTab3_graph, annotation_colors = annoCol3)


rownames(delRowsanno3) <- combinedCnDf$Gene[grep("del", combinedCnDf$Gene)]
delRowsanno3 <- t(delRowsanno3)
heatmap_anno3_dels <- pheatmap(log2(delRowsanno3)[heatmap_anno3_genes$tree_row$order,], cluster_rows = FALSE, cluster_cols = FALSE,
                               color = heatMapCol,breaks = colors.breaks, fontsize = 5, cellwidth = 5,
                               cellheight = 10,silent = FALSE, border_color = "black")


heatMapCol_tc <- colorRampPalette(c("#FFFFFF","#FFA500"))(100)
colors.breaks_tc <- seq(0,1,1/100)
tc_labels <- matchingTc
tc_labels<- t(tc_labels)
heatmap_anno3_tc <- pheatmap(tc_labels[heatmap_anno3_genes$tree_row$order], cluster_rows = FALSE, cluster_cols = FALSE,
                             color = heatMapCol_tc, breaks = colors.breaks_tc, fontsize = 5, cellwidth = 5,
                             cellheight = 10,silent = FALSE, border_color = "black")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210823choAnno3Heatmap_gene.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_anno3_genes
dev.off()


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210823choAnno3Heatmap_del.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_anno3_dels
dev.off()

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210823choAnno3Heatmap_tc.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_anno3_tc
dev.off()


###
###
###

annoTab4 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20210823choLabSim.xlsx")
annoTab4$ID <- str_remove_all(annoTab4$ID, " ")
annoTab4$ID <- str_remove_all(annoTab4$ID, "O")


anno4Idx <- match(annoTab4$ID, colnames(combinedCnDf))
combinedCnanno4 <- combinedCnDf[ , anno4Idx]
combinedCnanno4 <- cbind("Gene" = combinedCnDf$Gene, combinedCnanno4)

delRowsanno4 <- combinedCnanno4[grep("del", combinedCnanno4$Gene), 2:ncol(combinedCnanno4)]
combinedCnanno4 <- combinedCnanno4[-grep("del", combinedCnanno4$Gene),]
combinedCnDf$Gene <- firstUpper(combinedCnDf$Gene)
combinedCnanno4 <- combinedCnanno4[,-1]
combinedCnanno4[1:nrow(combinedCnanno4),] <- lapply(combinedCnanno4[1:nrow(combinedCnanno4),], function(x) log2(x))
rownames(combinedCnanno4) <- combinedCnDf$Gene[-grep("del", combinedCnDf$Gene)]

combinedCnanno4 <- t(combinedCnanno4)

matchingTc <- as.numeric(1 - delRowsanno4[1,])
# matchingTc <- tcDf$tc[match(tolower(rownames(combinedCnanno4)), tcDf$sample)]
combinedCnanno4[which(!is.na(matchingTc)),] <- sweep(combinedCnanno4[which(!is.na(matchingTc)),], 1,
                                                     matchingTc[which(!is.na(matchingTc))], "/")
combinedCnanno4[combinedCnanno4 > log2(2/3) & combinedCnanno4 < log2(4/3)] <- 0
combinedCnanno4 <- t(combinedCnanno4)
combinedCnanno4[combinedCnanno4 < -3] <- -3
combinedCnanno4[combinedCnanno4 > 3] <- 3

anno4TcIdx <- match(colnames(choSigMat),colnames(combinedCnanno4))
anno4TcIdx <- anno4TcIdx[-which(is.na(anno4TcIdx))]
choSigMatanno4 <- choSigMat[,anno4TcIdx]
combinedCnanno4[choSigMatanno4 == 0] <- 0

annoTab4_graph <- data.frame("Label" = annoTab4$Study, "Type" = annoTab4$Type,
                             stringsAsFactors = FALSE)
rownames(annoTab4_graph) <- annoTab4$ID
annoTab4_graph <- annoTab4_graph[which(rownames(annoTab4_graph) %in% colnames(combinedCnanno4)),]

annoCol4 <- list("Label" = c("High-dose simvastatin chow" = "#556b2f",
                             "Low-dose simvastatin chow" = "yellow",
                             "Control" = "black"),
                 "Type"  = c("HGSC" = "red", "SBRC" = "dodgerblue"))

heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-3,3,6/1000)

combinedCnanno4 <- t(combinedCnanno4)
heatmap_anno4_genes <- pheatmap(mat = combinedCnanno4, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                                breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                                border_color = "black", annotation_row = annoTab4_graph, annotation_colors = annoCol4)


rownames(delRowsanno4) <- combinedCnDf$Gene[grep("del", combinedCnDf$Gene)]
delRowsanno4 <- t(delRowsanno4)
heatmap_anno4_dels <- pheatmap(log2(delRowsanno4)[heatmap_anno4_genes$tree_row$order,], cluster_rows = FALSE, cluster_cols = FALSE,
                               color = heatMapCol,breaks = colors.breaks, fontsize = 5, cellwidth = 5,
                               cellheight = 10,silent = FALSE, border_color = "black")


heatMapCol_tc <- colorRampPalette(c("#FFFFFF","#FFA500"))(100)
colors.breaks_tc <- seq(0,1,1/100)
tc_labels <- matchingTc
tc_labels<- t(tc_labels)
heatmap_anno4_tc <- pheatmap(tc_labels[heatmap_anno4_genes$tree_row$order], cluster_rows = FALSE, cluster_cols = FALSE,
                             color = heatMapCol_tc, breaks = colors.breaks_tc, fontsize = 5, cellwidth = 5,
                             cellheight = 10,silent = FALSE, border_color = "black")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210823choanno4Heatmap_gene.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_anno4_genes
dev.off()


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210823choanno4Heatmap_del.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_anno4_dels
dev.off()

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210823choanno4Heatmap_tc.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_anno4_tc
dev.off()

### gene tables for xiaomang

write.table(combinedCnAnno1, "/mnt/DATA5/tmp/kev/misc/20210825agingGeneMatrix.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(delRowsAnno1, "/mnt/DATA5/tmp/kev/misc/20210825agingDeletionMatrix.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(combinedCnanno2, "/mnt/DATA5/tmp/kev/misc/20210825PregnancyGeneMatrix.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(delRowsanno2, "/mnt/DATA5/tmp/kev/misc/20210825PregnancyDeletionMatrix.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(combinedCnanno3, "/mnt/DATA5/tmp/kev/misc/20210825PregnancyNf1GeneMatrix.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(delRowsanno3, "/mnt/DATA5/tmp/kev/misc/20210825PregnancyNf1DeletionMatrix.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(combinedCnanno4, "/mnt/DATA5/tmp/kev/misc/20210825SimvastatenGeneMatrix.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(delRowsanno4, "/mnt/DATA5/tmp/kev/misc/20210825SimvastatenDeletionMatrix.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


load("/mnt/DATA5/tmp/kev/misc/20221004ascastSeg.Robj")


logRseg <- ascat.bc$Tumor_LogR_segmented
markers <- ascat.bc$SNPpos
i <- unique(markers$Chr)[1]
res <- NULL
for (i in unique(markers$Chr)) {
  tmpMarker <- markers[which(markers == i),]
  tmpMarker2 <- tmpMarker[1:(nrow(tmpMarker)-1), ]
  tmpMarker2$end <- tmpMarker$Position[2:nrow(tmpMarker)] - 1
  res <- rbind(res, tmpMarker2)
}

logRseg_red <- logRseg[which(rownames(logRseg) %in% rownames(res)),]

logRsegDf <- data.frame("chr" = res$Chr, "start" = res$Position, "end" = res$end,
                        "markers" = rownames(logRseg_red), logRseg_red, check.names = FALSE)



source("/home/kevhu/scripts/20210802syntenyFunctions.R")

library(GenomicRanges)
library(stringr)
library(irr)
library(DescTools)

cat_to_cont <- function(df){
  df <- str_replace_all(df, "none", "0")
  df <- str_replace_all(df, "gain", "1")
  df <- str_replace_all(df, "loss", "-1")
  df <- as.numeric(df)
  return(df)
}

snpRes05 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/snpVsPanel/20220503snpRes05V3.xlsx")

allPloidyCalls <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220517allHgscPloidySnpComp_driverZero.xlsx")


panelV3 <- allPloidyCalls
colnames(panelV3)[2] <- "V3"


eros <- c("/mnt/DATA6/mouseData/copynumber/")
listOfDirectories <- c("Auto_user_AUS5-138-MG_cho_20210621_354_343", "Auto_user_AUS5-76-MG_test1_255_185",
                       "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                       "Auto_user_AUS5-120-MG_EFD4_BBN_334_304")


allNoshadSeg <- read.table("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genes/absoluteRes.txt", sep = "\t",
                           stringsAsFactors = FALSE, header = TRUE, fill = TRUE)


panelV3$sample[which(panelV3$sample == "2027lte")] <- "mg4"
panelV3$sample[which(panelV3$sample == "13085lt")] <- "mg15"
panelV3$sample[which(panelV3$sample == "14399rt")] <- "mg20"
panelV3$sample[which(panelV3$sample == "14656peritonealmt")] <- "14656peritnealmt"
panelV3$sample[which(panelV3$sample == "133576rt")] <- "13576rt"
panelV3$sample[which(panelV3$sample == "14150rt")] <- "14150lt"
panelV3$sample[which(panelV3$sample == "14154rot")] <- "14154lt"


snpRes05$sample[which(snpRes05$sample == "2027lte")] <- "mg4"
snpRes05$sample[which(snpRes05$sample == "13085lt")] <- "mg15"
snpRes05$sample[which(snpRes05$sample == "14399rt")] <- "mg20"
snpRes05$sample[which(snpRes05$sample == "14656peritonealmt")] <- "14656peritnealmt"
snpRes05$sample[which(snpRes05$sample == "133576rt")] <- "13576rt"
snpRes05$sample[which(snpRes05$sample == "14150rt")] <- "14150lt"
snpRes05$sample[which(snpRes05$sample == "14154rot")] <- "14154lt"

panelV3 <- panelV3[match(snpRes05$sample, panelV3$sample),]


### need gene level calls
listOfDir <- paste0("/mnt/DATA6/mouseData/copynumber/", c("Reanalysis_AUS5-76-MG_test1_217", "Auto_user_AUS5-120-MG_EFD4_BBN_334_304",
                                                         "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349","Auto_user_AUS5-138-MG_cho_20210621_354_343",
                                                          "Auto_user_AUS5-142-MG_cho_20210701_357_353"))

### formatting weird for dir 

i <- listOfDir[1]
allGeneTable <- NULL
for (i in seq_along(listOfDir)) {
  setwd(listOfDir[i])
  tmpTable <- read.table("cnMatrix_gene.txt", stringsAsFactors = FALSE, 
                         header = TRUE, row.names = NULL,
                         check.names = FALSE)
  if (i == 1) {
    allGeneTable <- tmpTable
  } else{
    allGeneTable <- cbind(allGeneTable, tmpTable[, 2:ncol(tmpTable)])
  }
}


allGeneTable2 <- allGeneTable
allGeneTable2[, 2:ncol(allGeneTable2)] <- log2(allGeneTable2[, 2:ncol(allGeneTable2)])
geneMatrix <- allGeneTable2[, 2:ncol(allGeneTable2)]
geneMatrix[abs(geneMatrix) < 0.2] <- 0

mouseBed <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed", sep = "\t", 
                       header = FALSE, stringsAsFactors = FALSE)


geneNames <- firstUpper(allGeneTable2$Gene)
allGeneTable2$Gene <- geneNames

allGeneTable2



unwantedColumns <- c("AmpliconId", "AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene", "NumGC", "Length", "GC", "TotalPool", "Weights")
allAmpTable <- NULL
for (i in seq_along(listOfDir)) {
  setwd(listOfDir[i])
  tmpTable <- read.table("cnAmplicon_matrix.txt", stringsAsFactors = FALSE, 
                         header = TRUE, row.names = NULL,
                         check.names = FALSE, sep = "\t")
  tmpTable <- tmpTable[, -which(colnames(tmpTable) %in% unwantedColumns)]
  if (i == 1) {
    allAmpTable <- tmpTable
  } else{
    allAmpTable<- cbind(allAmpTable, tmpTable[, 2:ncol(tmpTable)])
  }
}


genomeInfo <- read.table("cnAmplicon_matrix.txt", stringsAsFactors = FALSE, 
                       header = TRUE, row.names = NULL,
                       check.names = FALSE, sep = "\t")
genomeInfo  <- genomeInfo[, which(colnames(genomeInfo) %in% unwantedColumns)]

geneBed <- NULL
i <- geneNames[1]
for (i in geneNames) {
  tmpBed <- mouseBed[which(mouseBed$V8 == i), ]
  if (nrow(tmpBed) < 1) {
    print(i)
    next()
  }
  tmpVec <- c(tmpBed$V1[1], min(tmpBed$V2), max(tmpBed$V3), tmpBed$V8[1])
  geneBed <- rbind(geneBed, tmpVec)
}

geneBed <- data.frame(geneBed)
rownames(geneBed) <- NULL

geneBed <- geneBed[-grep(paste("Rb1", "Brca1", "Nf1", "Trp53", sep = "|"), geneBed$X4),]

### remove commonly altered genes in ov mice



geneMatrix_df <- data.frame("Genes" = allGeneTable2$Gene, geneMatrix, check.names = FALSE)
geneMatrix_df <- geneMatrix_df[which(geneMatrix_df$Genes %in% geneBed$X4),]

geneMatrix_df2 <- data.frame("Chromosome" = geneBed$X1, "Start" = geneBed$X2,
                             "End" = geneBed$X3, geneMatrix_df, check.names = FALSE)

samepleNames <- colnames(geneMatrix_df2)
samepleNames <- nameStripper(samepleNames)


samepleNamesSnp <- colnames(allAmpTable)
samepleNamesSnp <- nameStripper(samepleNamesSnp)



# logRsegDf$chr <- paste0("chr", logRsegDf$chr)
snpGrange <- GRanges(seqnames = paste0("chr", logRsegDf$chr), 
                     IRanges(logRsegDf$start, logRsegDf$end))

ngsGrange <- GRanges(seqnames = geneMatrix_df2$Chromosome,
                     IRanges(start = as.numeric(geneMatrix_df2$Start), 
                             end = as.numeric(geneMatrix_df2$End)))

ngsAmpGrange <- GRanges(seqnames = paste0("chr", genomeInfo$ChromNum),
                     IRanges(start = as.numeric(genomeInfo$StartPos), 
                             end = as.numeric(genomeInfo$EndPos)))



# i <- colnames(logRseg_red)[1]

colnamesStripped <- str_remove(nameStripper(colnames(logRseg_red)), "\\-")
colnamesStripped[14] <- "13085rt"

finalCompRes <- NULL
namesCol <- NULL
for (i in seq_along(colnames(logRseg_red))) {
  geneSignal <- 2^unlist(geneMatrix_df2[which(samepleNames %in% colnamesStripped[i])])
  if (length(geneSignal) == 0) {
    print(colnamesStripped[i])
    next()
  }
  tmpSnp <- 2^unlist(logRseg_red[, colnames(logRseg_red)[i]])
  # tmpSnp[abs(tmpSnp) < 0.2] <- 0
  tmpResSnp <- NULL
  for (j in seq_along(ngsGrange)) {
    tmpOv <- findOverlaps(ngsGrange[j], snpGrange)
    tmpResSnp <- c(tmpResSnp, mean(tmpSnp[subjectHits(tmpOv)]))
  }
  # tmpVec <- list(colnames(logRseg_red)[i], tmpResSnp/geneSignal)
  # finalCompRes <- rbind(finalCompRes, tmpVec)
  finalCompRes <- cbind(finalCompRes, tmpResSnp/geneSignal)
  namesCol <- c(namesCol, colnames(logRseg_red)[i])
}

colnames(finalCompRes) <- namesCol

finalCompRes2 <- finalCompRes
rownames(finalCompRes2) <- NULL
finalCompRes2 <- data.frame(finalCompRes2, check.names = FALSE)
# finalCompRes2$Genes <- geneMatrix_df$Genes
geneStats <- apply(finalCompRes2, 1, function(x) quantile(x, seq(0,1, 0.05)))

geneStatsDf <- data.frame(geneStats)
colnames(geneStatsDf) <- geneMatrix_df$Genes


### snp

finalSnp <- NULL
namesCol2 <- NULL
for (i in seq_along(colnames(logRseg_red))) {
  geneSignal <- unlist(allAmpTable[which(samepleNamesSnp %in% colnamesStripped[i])])
  if (length(geneSignal) == 0) {
    print(colnamesStripped[i])
    next()
  }
  tmpSnp <- 2^unlist(logRseg_red[, colnames(logRseg_red)[i]])
  # tmpSnp[abs(tmpSnp) < 0.2] <- 0
  tmpResSnp <- NULL
  for (j in seq_along(ngsAmpGrange)) {
    tmpOv <- findOverlaps(ngsAmpGrange[j], snpGrange)
    tmpResSnp <- c(tmpResSnp, mean(tmpSnp[subjectHits(tmpOv)]))
  }
  # tmpVec <- list(colnames(logRseg_red)[i], tmpResSnp/geneSignal)
  # finalCompRes <- rbind(finalCompRes, tmpVec)
  finalSnp <- cbind(finalSnp, tmpResSnp/geneSignal)
  namesCol2 <- c(namesCol2, colnames(logRseg_red)[i])
}

colnames(finalSnp) <- namesCol2
finalSnp2 <- finalSnp
finalSnp2
### it seems the deciles per gene are independent of the number of amplicons used to detect a gene - noise is sample based
### think on this and look to just do by probe
sampleCv <- apply(finalCompRes2, 2, function(x) sd(unlist(x))/mean(unlist(x)))
geneCv <- apply(finalCompRes2, 1, function(x) sd(unlist(x))/mean(unlist(x)))

sampleZ <- apply(finalCompRes2, 2, function(x) (unlist(x) - mean(unlist(x)))/sd(unlist(x)))
geneZ <- apply(finalCompRes2, 1, function(x) (unlist(x) - mean(unlist(x)))/sd(unlist(x)))

colnames(geneZ) <- geneMatrix_df$Genes

ampCounts <- NULL
for (i in geneMatrix_df$Genes) {
  tmpMBed <- mouseBed[which(mouseBed$V8 == i),]
  tmpCounts <- nrow(tmpMBed)
  names(tmpCounts) <- i
  ampCounts <- c(ampCounts, tmpCounts)
}

gcContent <- NULL
for (i in colnames(geneZ)) {
  tmpGc <- mean(genomeInfo$GC[which(genomeInfo$Gene == i)])
  gcContent <- c(gcContent, tmpGc)
}

names(gcContent) <- colnames(geneZ)
meanGene <- apply(finalCompRes2, 1, mean)
names(meanGene) <- colnames(geneZ)

meanGene[which(meanGene < 0.9)]
meanGene[which(meanGene > 1.1)]

gcContent[which(names(gcContent) %in% names(meanGene[which(meanGene < 0.9)]))]
gcContent[which(names(gcContent) %in% names(meanGene[which(meanGene > 1.1)]))]


geneCvSnp <- apply(finalSnp2, 1, function(x) sd(x)/mean(x))

par(mfrow=c(2,2))

par(cex.axis = 0.5)
boxplot(geneStatsDf, las = 2, ylim = c(0, 2), main = "cnr ratio ngs/snp",
        xlab = "genes", ylab = "cnr(ngs)/cnr(snp)")
axis(2, at = seq(0, 2, 0.1), las = 2)

plot(ampCounts, geneCv,  main = "coefficient of var. vs # amplicons",
     xlab = "coefficient of var per gene", ylab = "amplicon count")

par(cex.axis = 0.5)
boxplot(t(geneZ), las = 2, ylim = c(-3, 3), main = "z-score of genes between samples",
        xlab = "z-score", ylab = "sample")
axis(2, at = seq(-3, 3, 0.5), las = 2)

plot(gcContent, meanGene, main = "GC content vs meanGene", 
     xlab = "GC content", ylab = "mean gene ratio of cnr (ngs/snp)")

### looking more in depth into SNps


boxplot(t(finalSnp2), xaxt="n", ylim = c(0,2))
boxplot(geneCvSnp, ylim = c(0,1))
quantile(geneCvSnp, seq(0,1, 0.05), na.rm = TRUE)
genomeInfo[which(is.na(as.numeric(geneCvSnp))),]
genomeInfo[which(geneCvSnp > 0.3),]
geneCvSnp[which(geneCvSnp > 0.3)]
rownames(finalSnp2) <- NULL
allProbes <- as.vector(finalSnp2)
allProbes <- allProbes[-which(is.na(allProbes))]
meanAll <- mean(allProbes)
sdAll <- sd(allProbes)
snpZscore <- apply(finalSnp2, 1, function(x) (x - meanAll)/sdAll)

boxplot(snpZscore, ylim = c(0,5), xaxt="n")
boxplot(snpZscore, xaxt="n")

genomeInfo$cv <- geneCvSnp

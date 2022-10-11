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



geneBed <- NULL
i <- geneNames[1]
for (i in geneNames) {
  tmpBed <- mouseBed[which(mouseBed$V8 == i), ]
  tmpVec <- c(tmpBed$V1[1], min(tmpBed$V2), max(tmpBed$V3), tmpBed$V8[1])
  geneBed <- rbind(geneBed, tmpVec)
}

geneBed <- data.frame(geneBed)
geneBed <- geneBed[-which(is.na(geneBed$X1)),]
rownames(geneBed) <- NULL

geneMatrix_df <- data.frame("Genes" = allGeneTable2$Gene, geneMatrix, check.names = FALSE)
geneMatrix_df <- geneMatrix_df[which(geneMatrix_df$Genes %in% geneBed$X4),]

geneMatrix_df2 <- data.frame("Chromosome" = geneBed$X1, "Start" = geneBed$X2,
                             "End" = geneBed$X3, geneMatrix_df, check.names = FALSE)

samepleNames <- colnames(geneMatrix_df2)
samepleNames <- nameStripper(samepleNames)


# logRsegDf$chr <- paste0("chr", logRsegDf$chr)
snpGrange <- GRanges(seqnames = logRsegDf$chr, 
                     IRanges(logRsegDf$start, logRsegDf$end))

ngsGrange <- GRanges(seqnames = geneMatrix_df2$Chromosome,
                     IRanges(start = as.numeric(geneMatrix_df2$Start), 
                             end = as.numeric(geneMatrix_df2$End)))

i <- colnames(logRseg_red)[1]

colnamesStripped <- str_remove(nameStripper(colnames(logRseg_red)), "\\-")
colnamesStripped[14] <- "13085rt"

finalCompRes <- NULL
for (i in seq_along(colnames(logRseg_red))) {
  geneSignal <- unlist(geneMatrix_df2[which(samepleNames %in% colnamesStripped[i])])
  tmpSnp <- unlist(logRseg_red[, colnames(logRseg_red)[i]])
  # tmpSnp[abs(tmpSnp) < 0.2] <- 0
  tmpResSnp <- NULL
  for (j in seq_along(ngsGrange)) {
    tmpOv <- findOverlaps(ngsGrange[j], snpGrange)
    tmpResSnp <- c(tmpResSnp, mean(tmpSnp[subjectHits(tmpOv)]))
  }
  tmpVec <- c(colnames(logRseg_red)[i], tmpResSnp/geneSignal)
  finalCompRes <- rbind(finalCompRes, tmpVec)
}

which(is.na(as.numeric(finalCompRes[,2])))

finalCompRes2 <- finalCompRes[-which(is.na(as.numeric(finalCompRes[,2]))), ]
finalCompRes2[, 2:ncol(finalCompRes2)] <- as.numeric(finalCompRes2[, 2:ncol(finalCompRes2)])
finalCompRes3 <- t(finalCompRes2)
colnames(finalCompRes3) <- finalCompRes3[1,]
finalCompRes3  <- finalCompRes3[-1,]
rownames(finalCompRes3) <- NULL
finalCompRes3 <- data.frame(finalCompRes3)
finalCompRes3 <- lapply(finalCompRes3, function(x) unlist(as.numeric(x)))
boxplot(finalCompRes3, las = 2, ylim = (c(-1, 1)))

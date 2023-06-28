### idea is to calculate RMSE per sample 
### should pick out 5 samples with clean profiles and diploid + relatively high estimated tumor content
### per sample calculate rmse by and group them by gamma

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

# allPloidyCalls <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220517allHgscPloidySnpComp.xlsx")
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

allNoshadSeg_filt <- allNoshadSeg[grep(paste0(paste(panelV3$sample, sprintf("%0.2f", panelV3$purityV2),
                                                    panelV3$ploidyV2, sep = "_"), collapse = "|"), allNoshadSeg$sample),]
allNoshadSeg_filt$sC <- as.numeric(allNoshadSeg_filt$sC)
allNoshadSeg_filt[,2:3] <- lapply(allNoshadSeg_filt[,2:3], as.numeric)
dupVec <- paste0(allNoshadSeg_filt$sample, allNoshadSeg_filt$chr, allNoshadSeg_filt$start, allNoshadSeg_filt$end)
allNoshadSeg_filt <- allNoshadSeg_filt[-which(duplicated(dupVec)), ]


### the five samples chosen ....
### make sure top use round function or cutoff to determine only clonal events when comparing

# sampleCompList <- c("mg4", "13576rt", "1628lt", "15723lt", "13179lt", "12643lt", "15676rt", "1685rt")
# samplePloidy <- c(2, 1.6, 2, 2, 2, 2, 2, 2)

sampleCompList <- c("mg4", "1628lt", "15723lt", "13179lt", "12643lt", "15676rt", "1685rt")
samplePloidy <- c(2, 2, 2, 2, 2, 2, 2)


gofDf <- read.table("/mnt/DATA5/tmp/kev/misc/gammaOptimalResAll.txt", sep = "\t", stringsAsFactors = FALSE,
                    header = TRUE)
gofDf <- gofDf[grep(paste0("psi", 2, "$"), gofDf$sample), ]

gam15 <- read.table("/mnt/DATA5/tmp/kev/misc/allSegResults_0.15.txt", sep = "\t",
                    stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
gam15 <- gam15[grep(paste0("psi", 2, "$"), gam15$sample), ]
gam15 <- gam15[grep(paste0(sampleCompList, collapse = "|"), gam15$sample), ]
gam15 <- gam15[which(gam15$chr %in% 1:19), ]
gam15$chr <- paste0("chr", gam15$chr)
gam15$totalSc <- as.numeric(gam15$nMajor) + as.numeric(gam15$nMinor)


gam20 <- read.table("/mnt/DATA5/tmp/kev/misc/allSegResults_0.2.txt", sep = "\t",
                    stringsAsFactors = FALSE, header = TRUE,  fill = TRUE)
gam20 <- gam20[grep(paste0("psi", 2, "$"), gam20$sample), ]
gam20 <- gam20[grep(paste0(sampleCompList, collapse = "|"), gam20$sample), ]
gam20 <- gam20[which(gam20$chr %in% 1:19), ]
gam20$chr <- paste0("chr", gam20$chr)
gam20$totalSc <- as.numeric(gam20$nMajor) + as.numeric(gam20$nMinor)


gam25 <- read.table("/mnt/DATA5/tmp/kev/misc/allSegResults_0.25.txt", sep = "\t",
                    stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
gam25 <- gam25[grep(paste0("psi", 2, "$"), gam25$sample), ]
gam25 <- gam25[grep(paste0(sampleCompList, collapse = "|"), gam25$sample), ]
gam25 <- gam25[which(gam25$chr %in% 1:19), ]
gam25$chr <- paste0("chr", gam25$chr)
gam25$totalSc <- as.numeric(gam25$nMajor) + as.numeric(gam25$nMinor)


gam30 <- read.table("/mnt/DATA5/tmp/kev/misc/allSegResults_0.3.txt", sep = "\t",
                    stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
gam30 <- gam30[grep(paste0("psi", 2, "$"), gam30$sample), ]
gam30 <- gam30[grep(paste0(sampleCompList, collapse = "|"), gam30$sample), ]
gam30 <- gam30[which(gam30$chr %in% 1:19), ]
gam30$chr <- paste0("chr", gam30$chr)
gam30$totalSc <- as.numeric(gam30$nMajor) + as.numeric(gam30$nMinor)

gam35 <- read.table("/mnt/DATA5/tmp/kev/misc/allSegResults_0.35.txt", sep = "\t",
                    stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
gam35 <- gam35[grep(paste0("psi", 2, "$"), gam35$sample), ]
gam35 <- gam35[grep(paste0(sampleCompList, collapse = "|"), gam35$sample), ]
gam35 <- gam35[which(gam35$chr %in% 1:19), ]
gam35$chr <- paste0("chr", gam35$chr)
gam35$totalSc <- as.numeric(gam35$nMajor) + as.numeric(gam35$nMinor)


gam40 <- read.table("/mnt/DATA5/tmp/kev/misc/allSegResults_0.4.txt", sep = "\t",
                    stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
gam40 <- gam40[grep(paste0("psi", 2, "$"), gam40$sample), ]
gam40 <- gam40[grep(paste0(sampleCompList, collapse = "|"), gam40$sample), ]
gam40 <- gam40[which(gam40$chr %in% 1:19), ]
gam40$chr <- paste0("chr", gam40$chr)
gam40$totalSc <- as.numeric(gam40$nMajor) + as.numeric(gam40$nMinor)


gam45 <- read.table("/mnt/DATA5/tmp/kev/misc/allSegResults_0.45.txt", sep = "\t",
                    stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
gam45 <- gam45[grep(paste0("psi", 2, "$"), gam45$sample), ]
gam45 <- gam45[grep(paste0(sampleCompList, collapse = "|"), gam45$sample), ]
gam45 <- gam45[which(gam45$chr %in% 1:19), ]
gam45$chr <- paste0("chr", gam45$chr)
gam45$totalSc <- as.numeric(gam45$nMajor) + as.numeric(gam45$nMinor)


gam50 <- read.table("/mnt/DATA5/tmp/kev/misc/allSegResults_0.5.txt", sep = "\t",
                    stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
gam50 <- gam50[grep(paste0("psi", 2, "$"), gam50$sample), ]
gam50 <- gam50[grep(paste0(sampleCompList, collapse = "|"), gam50$sample), ]
gam50 <- gam50[which(gam50$chr %in% 1:19), ]
gam50$chr <- paste0("chr", gam50$chr)
gam50$totalSc <- as.numeric(gam50$nMajor) + as.numeric(gam50$nMinor)


gam55 <- read.table("/mnt/DATA5/tmp/kev/misc/allSegResults_0.55.txt", sep = "\t",
                    stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
gam55 <- gam55[grep(paste0("psi", 2, "$"), gam55$sample), ]
gam55 <- gam55[grep(paste0(sampleCompList, collapse = "|"), gam55$sample), ]
gam55 <- gam55[which(gam55$chr %in% 1:19), ]
gam55$chr <- paste0("chr", gam55$chr)
gam55$totalSc <- as.numeric(gam55$nMajor) + as.numeric(gam55$nMinor)


gam60 <- read.table("/mnt/DATA5/tmp/kev/misc/allSegResults_0.6.txt", sep = "\t",
                    stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
gam60 <- gam60[grep(paste0("psi", 2, "$"), gam60$sample), ]
gam60 <- gam60[grep(paste0(sampleCompList, collapse = "|"), gam60$sample), ]
gam60 <- gam60[which(gam60$chr %in% 1:19), ]
gam60$chr <- paste0("chr", gam60$chr)
gam60$totalSc <- as.numeric(gam60$nMajor) + as.numeric(gam60$nMinor)


allNoshadSeg_filt2 <- allNoshadSeg_filt[grep(paste0(sampleCompList, collapse = "|"), allNoshadSeg_filt$sample), ]
allNoshadSeg_filt2$roundSc <- round(allNoshadSeg_filt2$sC)


### redid below to iterate over each gamamdf 

# resTable <- NULL
# i <- unique(allNoshadSeg_filt2$sample)[1]
# for (i in unique(allNoshadSeg_filt2$sample)) {
#   tmpSample <- str_remove(i, "\\_.*")
#   
#   ### get 1Mb bins for panel with corresponding value
#   tmpDf <- allNoshadSeg_filt2[which(allNoshadSeg_filt2$sample == i), ]
#   tmpDf <- tmpDf[which(tmpDf$chr %in% paste0("chr", 1:19)),]
#   tmpGR <- GRanges(seqnames = tmpDf$chr, IRanges(start = tmpDf$start, end = tmpDf$end))
#   tmpBins <- unlist(tile(x = tmpGR, width = 1e6))
#   
#   ### rough way of getting rid of dupe bins .... the bins may have different levels
#   ### of contribution so it's not the best ... only ~ 30 bins have more than 1
#   binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmpGR))))
#   tmpHits <- subjectHits(findOverlaps(tmpBins, tmpGR))
#   tmpHits <- tmpHits[-binDups]
#   tmpBinSc <- tmpDf$roundSc[tmpHits]
#   
#   ### iterate over different gammas
#   ### get values from bins create before
#   ### length should be same. create rmse for it
#   
#   tmp15 <- gam15[grep(tmpSample, gam15$sample), ]
#   dupVec <- paste0(tmp15$sample, tmp15$chr, tmp15$startpos, tmp15$endpos)
#   if (length(which(duplicated(dupVec))) > 1) {
#     tmp15 <- tmp15[-which(duplicated(dupVec)),]
#   }
#   j <- unique(tmp15$sample)[1]
#   for (j in unique(tmp15$sample)) {
#     tmp15_2 <- tmp15[which(tmp15$sample == j),]
#     tmp15_2_grange <- GRanges(tmp15_2$chr, IRanges(as.numeric(tmp15_2$startpos), as.numeric(tmp15_2$endpos)))
#     binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmp15_2_grange))))
#     tmpHits2 <- subjectHits(findOverlaps(tmpBins, tmp15_2_grange))
#     tmpHits2 <- tmpHits2[-binDups]
#     
#     ascatBins <- tmp15_2$totalSc[tmpHits2]
#     
#     resTable <- rbind(resTable, list(j, tmpSample, 0.15, signif(RMSE(ascatBins, tmpBinSc), 3)))
#     
#   }
#   
#   tmp20 <- gam20[grep(tmpSample, gam20$sample), ]
#   dupVec <- paste0(tmp20$sample, tmp20$chr, tmp20$startpos, tmp20$endpos)
#   if (length(which(duplicated(dupVec))) > 1) {
#     tmp20 <- tmp20[-which(duplicated(dupVec)),]
#   }
#   j <- unique(tmp20$sample)[1]
#   for (j in unique(tmp20$sample)) {
#     tmp20_2 <- tmp20[which(tmp20$sample == j),]
#     tmp20_2_grange <- GRanges(tmp20_2$chr, IRanges(as.numeric(tmp20_2$startpos), as.numeric(tmp20_2$endpos)))
#     binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmp20_2_grange))))
#     tmpHits2 <- subjectHits(findOverlaps(tmpBins, tmp20_2_grange))
#     tmpHits2 <- tmpHits2[-binDups]
#     
#     ascatBins <- tmp20_2$totalSc[tmpHits2]
#     
#     resTable <- rbind(resTable, list(j, tmpSample, 0.20, signif(RMSE(ascatBins, tmpBinSc), 3)))
#     
#   }
#   
#   tmp25 <- gam25[grep(tmpSample, gam25$sample), ]
#   dupVec <- paste0(tmp25$sample, tmp25$chr, tmp25$startpos, tmp25$endpos)
#   if (length(which(duplicated(dupVec))) > 1) {
#     tmp25 <- tmp25[-which(duplicated(dupVec)),]
#   }
#   j <- unique(tmp25$sample)[1]
#   for (j in unique(tmp25$sample)) {
#     tmp25_2 <- tmp25[which(tmp25$sample == j),]
#     tmp25_2_grange <- GRanges(tmp25_2$chr, IRanges(as.numeric(tmp25_2$startpos), as.numeric(tmp25_2$endpos)))
#     binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmp25_2_grange))))
#     tmpHits2 <- subjectHits(findOverlaps(tmpBins, tmp25_2_grange))
#     tmpHits2 <- tmpHits2[-binDups]
#     
#     ascatBins <- tmp25_2$totalSc[tmpHits2]
#     
#     resTable <- rbind(resTable, list(j, tmpSample, 0.25, signif(RMSE(ascatBins, tmpBinSc), 3)))
#     
#   }
#   
#   tmp30 <- gam30[grep(tmpSample, gam30$sample), ]
#   dupVec <- paste0(tmp30$sample, tmp30$chr, tmp30$startpos, tmp30$endpos)
#   if (length(which(duplicated(dupVec))) > 1) {
#     tmp30 <- tmp30[-which(duplicated(dupVec)),]
#   }
#   j <- unique(tmp30$sample)[1]
#   for (j in unique(tmp30$sample)) {
#     tmp30_2 <- tmp30[which(tmp30$sample == j),]
#     tmp30_2_grange <- GRanges(tmp30_2$chr, IRanges(as.numeric(tmp30_2$startpos), as.numeric(tmp30_2$endpos)))
#     binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmp30_2_grange))))
#     tmpHits2 <- subjectHits(findOverlaps(tmpBins, tmp30_2_grange))
#     tmpHits2 <- tmpHits2[-binDups]
#     
#     ascatBins <- tmp30_2$totalSc[tmpHits2]
#     
#     resTable <- rbind(resTable, list(j, tmpSample, 0.30, signif(RMSE(ascatBins, tmpBinSc), 3)))
#     
#   }
#   
#   tmp35 <- gam35[grep(tmpSample, gam35$sample), ]
#   dupVec <- paste0(tmp35$sample, tmp35$chr, tmp35$startpos, tmp35$endpos)
#   if (length(which(duplicated(dupVec))) > 1) {
#     tmp35 <- tmp35[-which(duplicated(dupVec)),]
#   }
#   j <- unique(tmp35$sample)[1]
#   for (j in unique(tmp35$sample)) {
#     tmp35_2 <- tmp35[which(tmp35$sample == j),]
#     tmp35_2_grange <- GRanges(tmp35_2$chr, IRanges(as.numeric(tmp35_2$startpos), as.numeric(tmp35_2$endpos)))
#     binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmp35_2_grange))))
#     tmpHits2 <- subjectHits(findOverlaps(tmpBins, tmp35_2_grange))
#     tmpHits2 <- tmpHits2[-binDups]
#     
#     ascatBins <- tmp35_2$totalSc[tmpHits2]
#     
#     resTable <- rbind(resTable, list(j, tmpSample, 0.35, signif(RMSE(ascatBins, tmpBinSc), 3)))
#     
#   }
#   
#   tmp40 <- gam40[grep(tmpSample, gam40$sample), ]
#   dupVec <- paste0(tmp40$sample, tmp40$chr, tmp40$startpos, tmp40$endpos)
#   if (length(which(duplicated(dupVec))) > 1) {
#     tmp40 <- tmp40[-which(duplicated(dupVec)),]
#   }
#   j <- unique(tmp40$sample)[1]
#   for (j in unique(tmp40$sample)) {
#     tmp40_2 <- tmp40[which(tmp40$sample == j),]
#     tmp40_2_grange <- GRanges(tmp40_2$chr, IRanges(as.numeric(tmp40_2$startpos), as.numeric(tmp40_2$endpos)))
#     binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmp40_2_grange))))
#     tmpHits2 <- subjectHits(findOverlaps(tmpBins, tmp40_2_grange))
#     tmpHits2 <- tmpHits2[-binDups]
#     
#     ascatBins <- tmp40_2$totalSc[tmpHits2]
#     
#     resTable <- rbind(resTable, list(j, tmpSample, 0.40, signif(RMSE(ascatBins, tmpBinSc), 3)))
#     
#   }
#   
#   tmp45 <- gam45[grep(tmpSample, gam45$sample), ]
#   dupVec <- paste0(tmp45$sample, tmp45$chr, tmp45$startpos, tmp45$endpos)
#   if (length(which(duplicated(dupVec))) > 1) {
#     tmp45 <- tmp45[-which(duplicated(dupVec)),]
#   }
#   j <- unique(tmp45$sample)[1]
#   for (j in unique(tmp45$sample)) {
#     tmp45_2 <- tmp45[which(tmp45$sample == j),]
#     tmp45_2_grange <- GRanges(tmp45_2$chr, IRanges(as.numeric(tmp45_2$startpos), as.numeric(tmp45_2$endpos)))
#     binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmp45_2_grange))))
#     tmpHits2 <- subjectHits(findOverlaps(tmpBins, tmp45_2_grange))
#     tmpHits2 <- tmpHits2[-binDups]
#     
#     ascatBins <- tmp45_2$totalSc[tmpHits2]
#     
#     resTable <- rbind(resTable, list(j, tmpSample, 0.45, signif(RMSE(ascatBins, tmpBinSc), 3)))
#     
#   }
#   
#   tmp50 <- gam50[grep(tmpSample, gam50$sample), ]
#   dupVec <- paste0(tmp50$sample, tmp50$chr, tmp50$startpos, tmp50$endpos)
#   if (length(which(duplicated(dupVec))) > 1) {
#     tmp50 <- tmp50[-which(duplicated(dupVec)),]
#   }
#   j <- unique(tmp50$sample)[1]
#   for (j in unique(tmp50$sample)) {
#     tmp50_2 <- tmp50[which(tmp50$sample == j),]
#     tmp50_2_grange <- GRanges(tmp50_2$chr, IRanges(as.numeric(tmp50_2$startpos), as.numeric(tmp50_2$endpos)))
#     binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmp50_2_grange))))
#     tmpHits2 <- subjectHits(findOverlaps(tmpBins, tmp50_2_grange))
#     tmpHits2 <- tmpHits2[-binDups]
#     
#     ascatBins <- tmp50_2$totalSc[tmpHits2]
#     
#     resTable <- rbind(resTable, list(j, tmpSample, 0.50, signif(RMSE(ascatBins, tmpBinSc), 3)))
#     
#   }
#   
#   tmp55 <- gam55[grep(tmpSample, gam55$sample), ]
#   dupVec <- paste0(tmp55$sample, tmp55$chr, tmp55$startpos, tmp55$endpos)
#   if (length(which(duplicated(dupVec))) > 1) {
#     tmp55 <- tmp55[-which(duplicated(dupVec)),]
#   }
#   j <- unique(tmp55$sample)[1]
#   for (j in unique(tmp55$sample)) {
#     tmp55_2 <- tmp55[which(tmp55$sample == j),]
#     tmp55_2_grange <- GRanges(tmp55_2$chr, IRanges(as.numeric(tmp55_2$startpos), as.numeric(tmp55_2$endpos)))
#     binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmp55_2_grange))))
#     tmpHits2 <- subjectHits(findOverlaps(tmpBins, tmp55_2_grange))
#     tmpHits2 <- tmpHits2[-binDups]
#     
#     ascatBins <- tmp55_2$totalSc[tmpHits2]
#     
#     resTable <- rbind(resTable, list(j, tmpSample, 0.55, signif(RMSE(ascatBins, tmpBinSc), 3)))
#     
#   }
#   
#   tmp60 <- gam60[grep(tmpSample, gam60$sample), ]
#   dupVec <- paste0(tmp60$sample, tmp60$chr, tmp60$startpos, tmp60$endpos)
#   if (length(which(duplicated(dupVec))) > 1) {
#     tmp60 <- tmp60[-which(duplicated(dupVec)),]
#   }
#   j <- unique(tmp60$sample)[1]
#   for (j in unique(tmp60$sample)) {
#     tmp60_2 <- tmp60[which(tmp60$sample == j),]
#     tmp60_2_grange <- GRanges(tmp60_2$chr, IRanges(as.numeric(tmp60_2$startpos), as.numeric(tmp60_2$endpos)))
#     binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmp60_2_grange))))
#     tmpHits2 <- subjectHits(findOverlaps(tmpBins, tmp60_2_grange))
#     tmpHits2 <- tmpHits2[-binDups]
#     
#     ascatBins <- tmp60_2$totalSc[tmpHits2]
#     
#     resTable <- rbind(resTable, list(j, tmpSample, 0.60, signif(RMSE(ascatBins, tmpBinSc), 3)))
#     
#   }
# }
# 


resTable <- NULL
i <- unique(allNoshadSeg_filt2$sample)[1]
for (i in unique(allNoshadSeg_filt2$sample)) {
  tmpSample <- str_remove(i, "\\_.*")
  
  ### get 1Mb bins for panel with corresponding value
  tmpDf <- allNoshadSeg_filt2[which(allNoshadSeg_filt2$sample == i), ]
  tmpDf <- tmpDf[which(tmpDf$chr %in% paste0("chr", 1:19)),]
  tmpGR <- GRanges(seqnames = tmpDf$chr, IRanges(start = tmpDf$start, end = tmpDf$end))
  tmpBins <- unlist(tile(x = tmpGR, width = 1e6))
  
  ### rough way of getting rid of dupe bins .... the bins may have different levels
  ### of contribution so it's not the best ... only ~ 30 bins have more than 1
  binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmpGR))))
  tmpHits <- subjectHits(findOverlaps(tmpBins, tmpGR))
  tmpHits <- tmpHits[-binDups]
  tmpBinSc <- tmpDf$roundSc[tmpHits]
  
  ### iterate over different gammas
  ### get values from bins create before
  ### length should be same. create rmse for it
  
  
  for (k in seq(15, 60, 5)) {
    tmpGam <- eval(parse(text = paste0("gam", k)))
    
    tmpGam2 <- tmpGam[grep(tmpSample, tmpGam$sample), ]
    dupVec <- paste0(tmpGam2$sample, tmpGam2$chr, tmpGam2$startpos, tmpGam2$endpos)
    if (length(which(duplicated(dupVec))) > 1) {
      tmpGam2 <- tmpGam2[-which(duplicated(dupVec)),]
    }
    
    tmpGammaDf <- tmpGam2
    for (j in unique(tmpGammaDf$sample)) {
      tmpGammaDf_2 <- tmpGammaDf[which(tmpGammaDf$sample == j),]
      tmpGammaDf_2_grange <- GRanges(tmpGammaDf_2$chr, IRanges(as.numeric(tmpGammaDf_2$startpos), as.numeric(tmpGammaDf_2$endpos)))
      binDups <- which(duplicated(queryHits(findOverlaps(tmpBins, tmpGammaDf_2_grange))))
      tmpHits2 <- subjectHits(findOverlaps(tmpBins, tmpGammaDf_2_grange))
      tmpHits2 <- tmpHits2[-binDups]
      
      ascatBins <- tmpGammaDf_2$totalSc[tmpHits2]
      
      tmpGof <- gofDf$goodnessOfFit[intersect(grep(j, gofDf$sample), grep(paste0(k/100, "$"), gofDf$gamma))]
      if(length(tmpGof) == 0){
        tmpGof <- NA
      }
      # gofDf[intersect(grep("1628lt_rho0.88_psi2", gofDf$sample), grep(paste0(0.2, "$"), gofDf$gamma)),]
 
      resTable <- rbind(resTable, list(j, tmpSample, k/100, signif(RMSE(ascatBins, tmpBinSc), 3), tmpGof))
    }
  }
}



resTable
resTableDf <- data.frame(resTable, stringsAsFactors = FALSE)
colnames(resTableDf) <- c("name", "sample", "gamma", "RMSE", "gof")
resTableDf <- lapply(resTableDf, unlist)
resTableDf <- data.frame(resTableDf, stringsAsFactors = FALSE)

a <- ggplot(resTableDf, aes(x = gamma, y = RMSE, color = sample)) + geom_point() + ggtitle("RMSE vs gamma: 1Mb bin absolute CN comparison") +
  theme(plot.title = element_text(hjust = 0.5))

b <- ggplot(resTableDf, aes(x = gamma, y = gof, color = sample)) + geom_point() + ggtitle("Goodness of fit vs gamma") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(a, b, ncol = 1)

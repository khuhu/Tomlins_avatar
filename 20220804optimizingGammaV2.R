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
# allPloidyCalls <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220517allHgscPloidySnpComp_driverZero.xlsx")
allPloidyCalls <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220719allHgscPloidyDriverGene.xlsx")

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

allNoshadSeg_filt <- allNoshadSeg[grep(paste0(paste(panelV3$sample, sprintf("%0.2f", panelV3$purity),
                                                    panelV3$ploidy, sep = "_"), collapse = "|"), allNoshadSeg$sample),]
allNoshadSeg_filt$sC <- as.numeric(allNoshadSeg_filt$sC)
allNoshadSeg_filt[,2:3] <- lapply(allNoshadSeg_filt[,2:3], as.numeric)
dupVec <- paste0(allNoshadSeg_filt$sample, allNoshadSeg_filt$chr, allNoshadSeg_filt$start, allNoshadSeg_filt$end)
allNoshadSeg_filt <- allNoshadSeg_filt[-which(duplicated(dupVec)), ]


### the five samples chosen ....
### make sure top use round function or cutoff to determine only clonal events when comparing

# sampleCompList <- c("mg4", "13576rt", "1628lt", "15723lt", "13179lt", "12643lt", "15676rt", "1685rt")
# samplePloidy <- c(2, 1.6, 2, 2, 2, 2, 2, 2)


# sampleCompList <- c("mg4", "1628lt", "15723lt", "13179lt", "12643lt", "15676rt", "1685rt")
# samplePloidy <- c(2, 2, 2, 2, 2, 2, 2, 1.6)

sampleCompList <- c("12641lt", "133576rt", "kc07", "kc10", "6786lt", "kc01", "2027_LT-E", "14399RT",
                    "14120rt", "2611lt", "3695rt", "2285lt", "12342lt", "kc06", "6826rt", "14943lt")
samplePloidy <- c(1.6, 1.6, 2.4, 3, 1.8, 2.8, 1.8, 1.8,
                  2.4, 4, 2.6, 2.6, 2, 1.8, 2, 3.8)


sampleCompList2 <- c("12641lt", "13576rt", "kc07", "kc10", "6786lt", "kc01", "mg4", "mg20",
                    "14120rt", "2611lt", "3695rt", "2285lt", "12342lt", "kc06", "6826rt", "14943lt")


samplePloidy2 <- c(2, 2, 3, 3, 2, 3, 2, 2,
                  3, 4, 3, 3, 2, 2, 2, 4)

### 14943lt has really odd fitting SNP profile and also 6826rt - not great

### need to find set of samples I'm confident in first

gofDf <- read.table("/mnt/DATA5/tmp/kev/misc/testingGammaDir/gammaOptimalResAll.txt", sep = "\t", stringsAsFactors = FALSE,
                    header = TRUE)


# gam15 <- read.table("/mnt/DATA5/tmp/kev/misc/allSegResults_0.15.txt", sep = "\t",
#                     stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
# gam15 <- gam15[grep(paste0("psi", 2, "$"), gam15$sample), ]
# gam15 <- gam15[grep(paste0(sampleCompList, collapse = "|"), gam15$sample), ]
# gam15 <- gam15[which(gam15$chr %in% 1:19), ]
# gam15$chr <- paste0("chr", gam15$chr)
# gam15$totalSc <- as.numeric(gam15$nMajor) + as.numeric(gam15$nMinor)

### creating for loop to get all gammas values
gamDir <- "/mnt/DATA5/tmp/kev/misc/testingGammaDir/"
setwd(gamDir)
gamFilles <- system("ls *allSeg*", intern = TRUE)
gamFiles2 <- paste0(gamDir, gamFilles)

seqList <- seq(0.05, 1, 0.05)
i <- 1
for (i in seq_along(gamFiles2)) {
  tmpTable <- read.table(gamFiles2[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
  tmpTable <- tmpTable[which(tmpTable$chr %in% 1:19), ]
  tmpTable$chr <- paste0("chr", tmpTable$chr)
  tmpTable$totalSc <- as.numeric(tmpTable$nMajor) + as.numeric(tmpTable$nMinor)
  
  tmpIdx <- NULL
  j <- 1 
  for (j in seq_along(sampleCompList)) {
    tmpRes <- which(grepl(paste0("^", sampleCompList[j], ".*psi", samplePloidy[j], "$"), tmpTable$sample))
    tmpIdx <- c(tmpIdx, tmpRes)
  }
  
  tmpTable$sample <- str_replace_all(tmpTable$sample, "133576rt", "13576rt")
  tmpTable$sample <- str_replace_all(tmpTable$sample, "2027_LT-E", "mg4")
  tmpTable$sample <- str_replace_all(tmpTable$sample, "14399RT", "mg20")
  
  "133576rt"
  
  tmpTable <- tmpTable[tmpIdx,]
  
  assign(paste0("gam", seqList[i]), tmpTable)
}



allNoshadSeg_filt2 <- allNoshadSeg_filt[grep(paste0(sampleCompList2, collapse = "|"), allNoshadSeg_filt$sample), ]
allNoshadSeg_filt2$roundSc <- round(allNoshadSeg_filt2$sC)

sampleDf <- data.frame("sample" = sampleCompList, "sample2" = sampleCompList2,
                       "ploidy" = samplePloidy, "ploidy2" = samplePloidy2,
                       stringsAsFactors = FALSE)

allNoshadSeg_filt2$sample2 <- str_remove(allNoshadSeg_filt2$sample, "\\_.*")
allNoshadSeg_filt2$gLsC <- allNoshadSeg_filt2$roundSc - sampleDf$ploidy2[match(allNoshadSeg_filt2$sample2, sampleDf$sample2)]


gofDf$sample <- str_replace_all(gofDf$sample, "133576rt", "13576rt")
gofDf$sample <- str_replace_all(gofDf$sample, "2027_LT-E", "mg4")
gofDf$sample <- str_replace_all(gofDf$sample, "14399RT", "mg20")

resTable <- NULL
i <- unique(allNoshadSeg_filt2$sample)[1]
for (i in unique(allNoshadSeg_filt2$sample)) {
  tmpSample <- str_remove(i, "\\_.*")
  
  print(tmpSample)
  ### get 1Mb bins for panel with corresponding value
  tmpDf <- allNoshadSeg_filt2[which(allNoshadSeg_filt2$sample == i), ]
  tmpDf <- tmpDf[which(tmpDf$chr %in% paste0("chr", 1:19)),]
  tmpGR <- GRanges(seqnames = tmpDf$chr, IRanges(start = tmpDf$start, end = tmpDf$end))
  tmpBins <- unlist(tile(x = tmpGR, width = 1e4))
  
  tmpNgsQuery <- queryHits(findOverlaps(tmpBins, tmpGR))
  tmpNgsSubject <- subjectHits(findOverlaps(tmpBins, tmpGR))
  
  
  registerDoParallel(28)
  res <- foreach(y = 1:length(tmpBins), .combine='c') %dopar% {
    tmpNgsRes <- which(tmpNgsQuery == y)
    if (length(tmpNgsRes) == 1) {
      tmpRes <- tmpDf$gLsC[tmpNgsSubject[tmpNgsRes]]
    } else if (length(tmpNgsRes) > 1){
      tmpRes <- mean(tmpDf$gLsC[tmpNgsSubject[tmpNgsRes]])
    } else {
      tmpRes <- 0
    }
    tmpRes
  }
  
  stopImplicitCluster()
  tmpBinSc <- res
  
  
  
  
  ### iterate over different gammas
  ### get values from bins create before
  ### length should be same. create rmse for it
  
  k <- 1
  registerDoParallel(20)
  resPar <- foreach(k = seq(0.05, 1, 0.05), .combine='rbind') %dopar% {
    # for (k in seq(0.05, 1, 0.05)) {
    tmpGam <- eval(parse(text = paste0("gam", k)))
    
    tmpGam2 <- tmpGam[grep(tmpSample, tmpGam$sample), ]
    dupVec <- paste0(tmpGam2$sample, tmpGam2$chr, tmpGam2$startpos, tmpGam2$endpos)
    if (length(which(duplicated(dupVec))) > 1) {
      tmpGam2 <- tmpGam2[-which(duplicated(dupVec)),]
    }
    
    tmpGammaDf <- tmpGam2
    tmpGammaDf$sample2 <- str_remove(tmpGammaDf$sample, "\\_.*")
    tmpGammaDf$gLsC <- tmpGammaDf$totalSc - sampleDf$ploidy2[match(tmpGammaDf$sample2, sampleDf$sample2)]
    
    
    j <- unique(tmpGammaDf$sample)[1]
    tmpParRes <- NULL
    for (j in unique(tmpGammaDf$sample)) {
      tmpGammaDf_2 <- tmpGammaDf[which(tmpGammaDf$sample == j),]
      tmpGammaDf_2_grange <- GRanges(tmpGammaDf_2$chr, IRanges(as.numeric(tmpGammaDf_2$startpos), as.numeric(tmpGammaDf_2$endpos)))
      
      tmpAscatQuery <- queryHits(findOverlaps(tmpBins, tmpGammaDf_2_grange))
      tmpAscatSubject <- subjectHits(findOverlaps(tmpBins, tmpGammaDf_2_grange))
      
      
      # registerDoParallel(28)
      res2 <- rep(0, length(tmpBins))
      
      
      for(i in 1:length(tmpBins)) {
        tmpAscatRes <- which(tmpAscatQuery == i)
        if (length(tmpAscatRes) == 1) {
          res2[i] <- tmpGammaDf_2$gLsC[tmpAscatSubject[tmpAscatRes]]
        } else if (length(tmpAscatRes) > 1){
          res2[i] <- mean(tmpGammaDf_2$gLsC[tmpAscatSubject[tmpAscatRes]])
        } else {
          next()
        }
      }
      
      
      ascatBins <- res2
      
      # tmpGof <- gofDf$goodnessOfFit[intersect(grep(paste0(j, "$), gofDf$sample), grep(paste0(k, "$"), gofDf$gamma))]
      tmpGof <- gofDf$goodnessOfFit[which(grepl(paste0(j, "$"), gofDf$sample) & grepl(paste0(k, "$"), gofDf$gamma))]
      
      if (length(tmpGof) > 1) {
        tmpGof <- tmpGof[2]
      } else if(length(tmpGof) == 0){
        tmpGof <- NA
      }
      
      # gofDf[intersect(grep("1628lt_rho0.88_psi2", gofDf$sample), grep(paste0(0.2, "$"), gofDf$gamma)),]
      # resTable <- rbind(resTable, list(j, tmpSample, k, signif(RMSE(ascatBins2, tmpBinSc2), 3), tmpGof))
      tmpParRes <- rbind(tmpParRes, list(j, tmpSample, k, signif(RMSE(ascatBins, tmpBinSc), 3), tmpGof))
    }
    tmpParRes
  }
  stopImplicitCluster()
  resTable <- rbind(resTable, resPar)
}


### need the multiple combine so that I can return which indices where missing to remove from RMSE comparison
### can't use just zero because it'll remove things that were truly zero, not bins that were not detected in ASCAT


resTable
resTableDf <- data.frame(resTable, stringsAsFactors = FALSE)
colnames(resTableDf) <- c("name", "sample", "gamma", "RMSE", "gof")
resTableDf <- lapply(resTableDf, unlist)
resTableDf <- data.frame(resTableDf, stringsAsFactors = FALSE)

a <- ggplot(resTableDf, aes(x = gamma, y = RMSE, color = sample)) + 
  geom_point(alpha = 0.2) + ggtitle("RMSE vs gamma: 10kb bin absolute CN comparison") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5))

b <- ggplot(resTableDf, aes(x = gamma, y = gof, color = sample)) + 
  geom_point(alpha = 0.2) + ggtitle("Goodness of fit vs gamma") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(a, b, ncol = 1)


resTableDf$paste <- paste0(resTableDf$sample, resTableDf$gamma)
resTableDf2 <- NULL
i <- unique(resTableDf$paste)[1]
for (i in unique(resTableDf$paste)) {
  tmpResTable <- resTableDf[which(resTableDf$paste == i), ]
  tmpRMSE <- mean(tmpResTable$RMSE)
  tmpGof <- mean(tmpResTable$gof)
  resTableDf2 <- rbind(resTableDf2, c(tmpResTable$sample[1], tmpResTable$gamma[1],
                                      tmpRMSE, tmpGof))
}

resTableDf2 <- data.frame(resTableDf2)
colnames(resTableDf2) <- c("sample", "gamma", "RMSE", "gof")
resTableDf2[, 3:4] <- lapply(resTableDf2[, 3:4], as.numeric)
resTableDf2$ploidy <- as.character(sampleDf$ploidy2[match(resTableDf2$sample, sampleDf$sample2)])

c <- ggplot(resTableDf2, aes(x = gamma, y = RMSE, color = sample)) + 
  geom_point(alpha = 0.5) + ggtitle("RMSE vs gamma: 10kb bin absolute CN comparison") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5))

d <- ggplot(resTableDf2, aes(x = gamma, y = gof, color = sample)) + 
  geom_point(alpha = 0.5) + ggtitle("Goodness of fit vs gamma") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(c, d, ncol = 1)

c2 <- ggplot(resTableDf2, aes(x = gamma, y = RMSE, color = ploidy)) + 
  geom_point(alpha = 0.25) + ggtitle("RMSE vs gamma: 10kb bin absolute CN comparison") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5))

d2 <- ggplot(resTableDf2, aes(x = gamma, y = gof, color = ploidy)) + 
  geom_point(alpha = 0.25) + ggtitle("Goodness of fit vs gamma") +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100))

grid.arrange(c2, d2, ncol = 1)


for (i in unique(resTableDf2$gamma)) {
  tmp <- resTableDf2[which(resTableDf2$gamma == i), ]
  print(paste("gamma", tmp$gamma[1], "mean rmse", signif(median(tmp$RMSE), digits = 3),
              "sd rmse", signif(sd(tmp$RMSE), digits = 3),
              "gof mean", signif(mean(tmp$gof), digits = 3),
              "gof sd", signif(sd(tmp$gof), digits = 3)))
}

# save.image(file = "/mnt/DATA5/tmp/kev/misc/2022081610kbWindow.Rdata")


### can try labeling based on regular or polyploidy for the graphs
### for now redo the 10kb windows but with only gains and losses



resTableGL <- NULL
i <- unique(allNoshadSeg_filt2$sample)[3]
for (i in unique(allNoshadSeg_filt2$sample)) {
  tmpSample <- str_remove(i, "\\_.*")
  
  print(tmpSample)
  ### get 1Mb bins for panel with corresponding value
  tmpDf <- allNoshadSeg_filt2[which(allNoshadSeg_filt2$sample == i), ]
  tmpDf <- tmpDf[which(tmpDf$chr %in% paste0("chr", 1:19)),]
  tmpGR <- GRanges(seqnames = tmpDf$chr, IRanges(start = tmpDf$start, end = tmpDf$end))
  tmpBins <- unlist(tile(x = tmpGR, width = 1e4))
  
  tmpNgsQuery <- queryHits(findOverlaps(tmpBins, tmpGR))
  tmpNgsSubject <- subjectHits(findOverlaps(tmpBins, tmpGR))
  
  
  registerDoParallel(28)
  res <- foreach(y = 1:length(tmpBins), .combine='c') %dopar% {
    tmpNgsRes <- which(tmpNgsQuery == y)
    if (length(tmpNgsRes) == 1) {
      tmpRes <- tmpDf$gLsC[tmpNgsSubject[tmpNgsRes]]
    } else if (length(tmpNgsRes) > 1){
      tmpRes <- mean(tmpDf$gLsC[tmpNgsSubject[tmpNgsRes]])
    } else {
      tmpRes <- 0
    }
    tmpRes
  }
  
  stopImplicitCluster()
  tmpBinSc <- res
  
  
  
  
  ### iterate over different gammas
  ### get values from bins create before
  ### length should be same. create rmse for it
  
  k <- 1
  registerDoParallel(20)
  resPar <- foreach(k = seq(0.05, 1, 0.05), .combine='rbind') %dopar% {
  # for (k in seq(0.05, 1, 0.05)) {
    tmpGam <- eval(parse(text = paste0("gam", k)))
    
    tmpGam2 <- tmpGam[grep(tmpSample, tmpGam$sample), ]
    dupVec <- paste0(tmpGam2$sample, tmpGam2$chr, tmpGam2$startpos, tmpGam2$endpos)
    if (length(which(duplicated(dupVec))) > 1) {
      tmpGam2 <- tmpGam2[-which(duplicated(dupVec)),]
    }
    
    tmpGammaDf <- tmpGam2
    tmpGammaDf$sample2 <- str_remove(tmpGammaDf$sample, "\\_.*")
    tmpGammaDf$gLsC <- tmpGammaDf$totalSc - sampleDf$ploidy2[match(tmpGammaDf$sample2, sampleDf$sample2)]

    
    j <- unique(tmpGammaDf$sample)[1]
    tmpParRes <- NULL
    for (j in unique(tmpGammaDf$sample)) {
      tmpGammaDf_2 <- tmpGammaDf[which(tmpGammaDf$sample == j),]
      tmpGammaDf_2_grange <- GRanges(tmpGammaDf_2$chr, IRanges(as.numeric(tmpGammaDf_2$startpos), as.numeric(tmpGammaDf_2$endpos)))
      
      tmpAscatQuery <- queryHits(findOverlaps(tmpBins, tmpGammaDf_2_grange))
      tmpAscatSubject <- subjectHits(findOverlaps(tmpBins, tmpGammaDf_2_grange))
      
      
      # registerDoParallel(28)
      res2 <- rep(0, length(tmpBins))

      
      for(i in 1:length(tmpBins)) {
        tmpAscatRes <- which(tmpAscatQuery == i)
        if (length(tmpAscatRes) == 1) {
          res2[i] <- tmpGammaDf_2$gLsC[tmpAscatSubject[tmpAscatRes]]
        } else if (length(tmpAscatRes) > 1){
          res2[i] <- mean(tmpGammaDf_2$gLsC[tmpAscatSubject[tmpAscatRes]])
        } else {
          next()
        }
      }
      
      
      ascatBins <- res2
      
      # tmpGof <- gofDf$goodnessOfFit[intersect(grep(paste0(j, "$), gofDf$sample), grep(paste0(k, "$"), gofDf$gamma))]
      tmpGof <- gofDf$goodnessOfFit[which(grepl(paste0(j, "$"), gofDf$sample) & grepl(paste0(k, "$"), gofDf$gamma))]
      
      if (length(tmpGof) > 1) {
        tmpGof <- tmpGof[2]
      } else if(length(tmpGof) == 0){
        tmpGof <- NA
      }
      
      tmpBinSc2 <- tmpBinSc
      ascatBins2 <- ascatBins
      
      tmpIdx <- c(which(ascatBins2 != 0), which(tmpBinSc2 != 0))
      if (length(which(duplicated(tmpIdx))) > 1) {
        tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
      }
      
      tmpBinSc2 <- tmpBinSc2[tmpIdx]
      ascatBins2 <- ascatBins2[tmpIdx]
      
      # gofDf[intersect(grep("1628lt_rho0.88_psi2", gofDf$sample), grep(paste0(0.2, "$"), gofDf$gamma)),]
      # resTableGL <- rbind(resTableGL, list(j, tmpSample, k, signif(RMSE(ascatBins2, tmpBinSc2), 3), tmpGof))
      tmpParRes <- rbind(tmpParRes, list(j, tmpSample, k, signif(RMSE(ascatBins2, tmpBinSc2), 3), tmpGof))
    }
    tmpParRes
  }
  stopImplicitCluster()
  resTableGL <- rbind(resTableGL, resPar)
}




resTableGL
resTableGLDf <- data.frame(resTableGL, stringsAsFactors = FALSE)
colnames(resTableGLDf) <- c("name", "sample", "gamma", "RMSE", "gof")
resTableGLDf <- lapply(resTableGLDf, unlist)
resTableGLDf <- data.frame(resTableGLDf, stringsAsFactors = FALSE)

e <- ggplot(resTableGLDf, aes(x = gamma, y = RMSE, color = sample)) + 
  geom_point(alpha = 0.2) + ggtitle("RMSE vs gamma: 10kb bin absolute CN comparison") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5))

f <- ggplot(resTableGLDf, aes(x = gamma, y = gof, color = sample)) + 
  geom_point(alpha = 0.2) + ggtitle("Goodness of fit vs gamma") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(e, f, ncol = 1)


resTableGLDf$paste <- paste0(resTableGLDf$sample, resTableGLDf$gamma)
resTableGLDf2 <- NULL
i <- unique(resTableGLDf$paste)[1]
for (i in unique(resTableGLDf$paste)) {
  tmpresTableGL <- resTableGLDf[which(resTableGLDf$paste == i), ]
  tmpRMSE <- mean(tmpresTableGL$RMSE)
  tmpGof <- mean(tmpresTableGL$gof)
  resTableGLDf2 <- rbind(resTableGLDf2, c(tmpresTableGL$sample[1], tmpresTableGL$gamma[1],
                                      tmpRMSE, tmpGof))
}

resTableGLDf2 <- data.frame(resTableGLDf2)
colnames(resTableGLDf2) <- c("sample", "gamma", "RMSE", "gof")
resTableGLDf2[, 3:4] <- lapply(resTableGLDf2[, 3:4], as.numeric)
resTableGLDf2$ploidy <- as.character(sampleDf$ploidy2[match(resTableGLDf2$sample, sampleDf$sample2)])

g <- ggplot(resTableGLDf2, aes(x = gamma, y = RMSE, color = sample)) + 
  geom_point(alpha = 0.5) + ggtitle("RMSE vs gamma: 100kb bin absolute CN comparison: gains/losses") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5))

h <- ggplot(resTableGLDf2, aes(x = gamma, y = gof, color = sample)) + 
  geom_point(alpha = 0.5) + ggtitle("Goodness of fit vs gamma") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(g, h, ncol = 1)


g2 <- ggplot(resTableGLDf2, aes(x = gamma, y = RMSE, color = ploidy)) + 
  geom_point(alpha = 0.25) + ggtitle("RMSE vs gamma: 10kb bin absolute CN comparison: gains/losses") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5)) 

h2 <- ggplot(resTableGLDf2, aes(x = gamma, y = gof, color = ploidy)) + 
  geom_point(alpha = 0.25) + ggtitle("Goodness of fit vs gamma") + scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) + 
  theme(plot.title = element_text(hjust = 0.5)) 

grid.arrange(g2, h2, ncol = 1)

for (i in unique(resTableGLDf2$gamma)) {
  tmp <- resTableGLDf2[which(resTableGLDf2$gamma == i), ]
  print(paste("gamma", tmp$gamma[1], "mean rmse", signif(median(tmp$RMSE), digits = 3),
              "sd rmse", signif(sd(tmp$RMSE), digits = 3),
              "gof mean", signif(mean(tmp$gof), digits = 3),
              "gof sd", signif(sd(tmp$gof), digits = 3)))
}


### doing above but by gene
###
###


bedFile <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD124056_167_Designed.del.nopool.gc.bed",
                      sep = "\t", stringsAsFactors = FALSE, header = FALSE)

bedFile2 <- NULL
for (i in unique(bedFile$V8)) {
  tmpBed <- bedFile[which(bedFile$V8 == i),]
  ampCount <- nrow(tmpBed)
  
  if (ampCount < 5) {
    print(i)
    next()
  } else{
    tmpChr <- tmpBed$V1[1]
    tmpStart <- min(tmpBed$V2)
    tmpEnd <- max(tmpBed$V3)
    tmpGene <- i
    bedFile2 <- rbind(bedFile2, c(tmpChr, tmpStart, tmpEnd, tmpGene))
  }
}

bedFile2 <- data.frame(bedFile2, stringsAsFactors = FALSE)
bedFile2$X2 <- as.numeric(bedFile2$X2)
bedFile2$X3 <- as.numeric(bedFile2$X3)

geneTable <- NULL
i <- unique(allNoshadSeg_filt2$sample)[1]
for (i in unique(allNoshadSeg_filt2$sample)) {
  tmpSample <- str_remove(i, "\\_.*")
  
  print(tmpSample)
  tmpDf <- allNoshadSeg_filt2[which(allNoshadSeg_filt2$sample == i), ]
  tmpDf <- tmpDf[which(tmpDf$chr %in% paste0("chr", 1:19)),]
  tmpGR <- GRanges(seqnames = tmpDf$chr, IRanges(start = tmpDf$start, end = tmpDf$end))
  tmpBins <- GRanges(seqnames = bedFile2$X1, IRanges(start = bedFile2$X2, end = bedFile2$X3))
  
  tmpNgsQuery <- queryHits(findOverlaps(tmpBins, tmpGR))
  tmpNgsSubject <- subjectHits(findOverlaps(tmpBins, tmpGR))
  
  
  registerDoParallel(28)
  res <- foreach(y = 1:length(tmpBins), .combine='c') %dopar% {
    tmpNgsRes <- which(tmpNgsQuery == y)
    if (length(tmpNgsRes) == 1) {
      tmpRes <- tmpDf$gLsC[tmpNgsSubject[tmpNgsRes]]
    } else if (length(tmpNgsRes) > 1){
      tmpRes <- mean(tmpDf$gLsC[tmpNgsSubject[tmpNgsRes]])
    } else {
      tmpRes <- 0
    }
    tmpRes
  }
  
  stopImplicitCluster()
  tmpBinSc <- res
  
  
  
  
  ### iterate over different gammas
  ### get values from bins create before
  ### length should be same. create rmse for it
  
  k <- 1
  registerDoParallel(20)
  resPar <- foreach(k = seq(0.05, 1, 0.05), .combine='rbind') %dopar% {
    # for (k in seq(0.05, 1, 0.05)) {
    tmpGam <- eval(parse(text = paste0("gam", k)))
    
    tmpGam2 <- tmpGam[grep(tmpSample, tmpGam$sample), ]
    dupVec <- paste0(tmpGam2$sample, tmpGam2$chr, tmpGam2$startpos, tmpGam2$endpos)
    if (length(which(duplicated(dupVec))) > 1) {
      tmpGam2 <- tmpGam2[-which(duplicated(dupVec)),]
    }
    
    tmpGammaDf <- tmpGam2
    tmpGammaDf$sample2 <- str_remove(tmpGammaDf$sample, "\\_.*")
    tmpGammaDf$gLsC <- tmpGammaDf$totalSc - sampleDf$ploidy2[match(tmpGammaDf$sample2, sampleDf$sample2)]
    
    
    j <- unique(tmpGammaDf$sample)[1]
    tmpParRes <- NULL
    for (j in unique(tmpGammaDf$sample)) {
      tmpGammaDf_2 <- tmpGammaDf[which(tmpGammaDf$sample == j),]
      tmpGammaDf_2_grange <- GRanges(tmpGammaDf_2$chr, IRanges(as.numeric(tmpGammaDf_2$startpos), as.numeric(tmpGammaDf_2$endpos)))
      
      tmpAscatQuery <- queryHits(findOverlaps(tmpBins, tmpGammaDf_2_grange))
      tmpAscatSubject <- subjectHits(findOverlaps(tmpBins, tmpGammaDf_2_grange))
      
      
      # registerDoParallel(28)
      res2 <- rep(0, length(tmpBins))
      
      
      for(i in 1:length(tmpBins)) {
        tmpAscatRes <- which(tmpAscatQuery == i)
        if (length(tmpAscatRes) == 1) {
          res2[i] <- tmpGammaDf_2$gLsC[tmpAscatSubject[tmpAscatRes]]
        } else if (length(tmpAscatRes) > 1){
          res2[i] <- mean(tmpGammaDf_2$gLsC[tmpAscatSubject[tmpAscatRes]])
        } else {
          next()
        }
      }
      
      
      ascatBins <- res2
      
      # tmpGof <- gofDf$goodnessOfFit[intersect(grep(paste0(j, "$), gofDf$sample), grep(paste0(k, "$"), gofDf$gamma))]
      tmpGof <- gofDf$goodnessOfFit[which(grepl(paste0(j, "$"), gofDf$sample) & grepl(paste0(k, "$"), gofDf$gamma))]
      
      if (length(tmpGof) > 1) {
        tmpGof <- tmpGof[2]
      } else if(length(tmpGof) == 0){
        tmpGof <- NA
      }
      
      # gofDf[intersect(grep("1628lt_rho0.88_psi2", gofDf$sample), grep(paste0(0.2, "$"), gofDf$gamma)),]
      # geneTable <- rbind(geneTable, list(j, tmpSample, k, signif(RMSE(ascatBins2, tmpBinSc2), 3), tmpGof))
      tmpParRes <- rbind(tmpParRes, list(j, tmpSample, k, signif(RMSE(ascatBins, tmpBinSc), 3), tmpGof))
    }
    tmpParRes
  }
  stopImplicitCluster()
  geneTable <- rbind(geneTable, resPar)
}

geneTableDf <- data.frame(geneTable, stringsAsFactors = FALSE)
colnames(geneTableDf) <- c("name", "sample", "gamma", "RMSE", "gof")
geneTableDf <- lapply(geneTableDf, unlist)
geneTableDf <- data.frame(geneTableDf, stringsAsFactors = FALSE)

z <- ggplot(geneTableDf, aes(x = gamma, y = RMSE, color = sample)) + 
  geom_point(alpha = 0.2) + ggtitle("RMSE vs gamma: 10kb bin absolute CN comparison") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5))

y <- ggplot(geneTableDf, aes(x = gamma, y = gof, color = sample)) + 
  geom_point(alpha = 0.2) + ggtitle("Goodness of fit vs gamma") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(z, y, ncol = 1)


geneTableDf$paste <- paste0(geneTableDf$sample, geneTableDf$gamma)
geneTableDf2 <- NULL
i <- unique(geneTableDf$paste)[1]
for (i in unique(geneTableDf$paste)) {
  tmpgeneTable <- geneTableDf[which(geneTableDf$paste == i), ]
  tmpRMSE <- mean(tmpgeneTable$RMSE)
  tmpGof <- mean(tmpgeneTable$gof)
  geneTableDf2 <- rbind(geneTableDf2, c(tmpgeneTable$sample[1], tmpgeneTable$gamma[1],
                                      tmpRMSE, tmpGof))
}

geneTableDf2 <- data.frame(geneTableDf2)
colnames(geneTableDf2) <- c("sample", "gamma", "RMSE", "gof")
geneTableDf2[, 3:4] <- lapply(geneTableDf2[, 3:4], as.numeric)
geneTableDf2$ploidy <- as.character(sampleDf$ploidy2[match(geneTableDf2$sample, sampleDf$sample2)])

w <- ggplot(geneTableDf2, aes(x = gamma, y = RMSE, color = sample)) + 
  geom_point(alpha = 0.5) + ggtitle("RMSE vs gamma: 10kb bin absolute CN comparison") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5))

v <- ggplot(geneTableDf2, aes(x = gamma, y = gof, color = sample)) + 
  geom_point(alpha = 0.5) + ggtitle("Goodness of fit vs gamma") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(w, v, ncol = 1)

w2 <- ggplot(geneTableDf2, aes(x = gamma, y = RMSE, color = ploidy)) + 
  geom_point(alpha = 0.25) + ggtitle("RMSE vs gamma: gene absolute CN comparison") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5))

v2 <- ggplot(geneTableDf2, aes(x = gamma, y = gof, color = ploidy)) + 
  geom_point(alpha = 0.25) + ggtitle("Goodness of fit vs gamma") +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100))

grid.arrange(w2, v2, ncol = 1)




geneTableGL <- NULL
i <- unique(allNoshadSeg_filt2$sample)[3]
for (i in unique(allNoshadSeg_filt2$sample)) {
  tmpSample <- str_remove(i, "\\_.*")
  
  print(tmpSample)
  ### get 1Mb bins for panel with corresponding value
  tmpDf <- allNoshadSeg_filt2[which(allNoshadSeg_filt2$sample == i), ]
  tmpDf <- tmpDf[which(tmpDf$chr %in% paste0("chr", 1:19)),]
  tmpGR <- GRanges(seqnames = tmpDf$chr, IRanges(start = tmpDf$start, end = tmpDf$end))
  tmpBins <- GRanges(seqnames = bedFile2$X1, IRanges(start = bedFile2$X2, end = bedFile2$X3))
  
  tmpNgsQuery <- queryHits(findOverlaps(tmpBins, tmpGR))
  tmpNgsSubject <- subjectHits(findOverlaps(tmpBins, tmpGR))
  
  
  registerDoParallel(28)
  res <- foreach(y = 1:length(tmpBins), .combine='c') %dopar% {
    tmpNgsRes <- which(tmpNgsQuery == y)
    if (length(tmpNgsRes) == 1) {
      tmpRes <- tmpDf$gLsC[tmpNgsSubject[tmpNgsRes]]
    } else if (length(tmpNgsRes) > 1){
      tmpRes <- mean(tmpDf$gLsC[tmpNgsSubject[tmpNgsRes]])
    } else {
      tmpRes <- 0
    }
    tmpRes
  }
  
  stopImplicitCluster()
  tmpBinSc <- res
  
  
  
  
  ### iterate over different gammas
  ### get values from bins create before
  ### length should be same. create rmse for it
  
  k <- 1
  registerDoParallel(20)
  resPar <- foreach(k = seq(0.05, 1, 0.05), .combine='rbind') %dopar% {
    # for (k in seq(0.05, 1, 0.05)) {
    tmpGam <- eval(parse(text = paste0("gam", k)))
    
    tmpGam2 <- tmpGam[grep(tmpSample, tmpGam$sample), ]
    dupVec <- paste0(tmpGam2$sample, tmpGam2$chr, tmpGam2$startpos, tmpGam2$endpos)
    if (length(which(duplicated(dupVec))) > 1) {
      tmpGam2 <- tmpGam2[-which(duplicated(dupVec)),]
    }
    
    tmpGammaDf <- tmpGam2
    tmpGammaDf$sample2 <- str_remove(tmpGammaDf$sample, "\\_.*")
    tmpGammaDf$gLsC <- tmpGammaDf$totalSc - sampleDf$ploidy2[match(tmpGammaDf$sample2, sampleDf$sample2)]
    
    
    j <- unique(tmpGammaDf$sample)[1]
    tmpParRes <- NULL
    for (j in unique(tmpGammaDf$sample)) {
      tmpGammaDf_2 <- tmpGammaDf[which(tmpGammaDf$sample == j),]
      tmpGammaDf_2_grange <- GRanges(tmpGammaDf_2$chr, IRanges(as.numeric(tmpGammaDf_2$startpos), as.numeric(tmpGammaDf_2$endpos)))
      
      tmpAscatQuery <- queryHits(findOverlaps(tmpBins, tmpGammaDf_2_grange))
      tmpAscatSubject <- subjectHits(findOverlaps(tmpBins, tmpGammaDf_2_grange))
      
      
      # registerDoParallel(28)
      res2 <- rep(0, length(tmpBins))
      
      
      for(i in 1:length(tmpBins)) {
        tmpAscatRes <- which(tmpAscatQuery == i)
        if (length(tmpAscatRes) == 1) {
          res2[i] <- tmpGammaDf_2$gLsC[tmpAscatSubject[tmpAscatRes]]
        } else if (length(tmpAscatRes) > 1){
          res2[i] <- mean(tmpGammaDf_2$gLsC[tmpAscatSubject[tmpAscatRes]])
        } else {
          next()
        }
      }
      
      
      ascatBins <- res2
      
      # tmpGof <- gofDf$goodnessOfFit[intersect(grep(paste0(j, "$), gofDf$sample), grep(paste0(k, "$"), gofDf$gamma))]
      tmpGof <- gofDf$goodnessOfFit[which(grepl(paste0(j, "$"), gofDf$sample) & grepl(paste0(k, "$"), gofDf$gamma))]
      
      if (length(tmpGof) > 1) {
        tmpGof <- tmpGof[2]
      } else if(length(tmpGof) == 0){
        tmpGof <- NA
      }
      
      tmpBinSc2 <- tmpBinSc
      ascatBins2 <- ascatBins
      
      tmpIdx <- c(which(ascatBins2 != 0), which(tmpBinSc2 != 0))
      if (length(which(duplicated(tmpIdx))) > 1) {
        tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
      }
      
      tmpBinSc2 <- tmpBinSc2[tmpIdx]
      ascatBins2 <- ascatBins2[tmpIdx]
      
      # gofDf[intersect(grep("1628lt_rho0.88_psi2", gofDf$sample), grep(paste0(0.2, "$"), gofDf$gamma)),]
      # geneTableGL <- rbind(geneTableGL, list(j, tmpSample, k, signif(RMSE(ascatBins2, tmpBinSc2), 3), tmpGof))
      tmpParRes <- rbind(tmpParRes, list(j, tmpSample, k, signif(RMSE(ascatBins2, tmpBinSc2), 3), tmpGof))
    }
    tmpParRes
  }
  stopImplicitCluster()
  geneTableGL <- rbind(geneTableGL, resPar)
}



geneTableGLDf <- data.frame(geneTableGL, stringsAsFactors = FALSE)
colnames(geneTableGLDf) <- c("name", "sample", "gamma", "RMSE", "gof")
geneTableGLDf <- lapply(geneTableGLDf, unlist)
geneTableGLDf <- data.frame(geneTableGLDf, stringsAsFactors = FALSE)

q <- ggplot(geneTableGLDf, aes(x = gamma, y = RMSE, color = sample)) + 
  geom_point(alpha = 0.2) + ggtitle("RMSE vs gamma: 10kb bin absolute CN comparison") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5))

r <- ggplot(geneTableGLDf, aes(x = gamma, y = gof, color = sample)) + 
  geom_point(alpha = 0.2) + ggtitle("Goodness of fit vs gamma") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(q, r, ncol = 1)


geneTableGLDf$paste <- paste0(geneTableGLDf$sample, geneTableGLDf$gamma)
geneTableGLDf2 <- NULL
i <- unique(geneTableGLDf$paste)[1]
for (i in unique(geneTableGLDf$paste)) {
  tmpgeneTableGL <- geneTableGLDf[which(geneTableGLDf$paste == i), ]
  tmpRMSE <- mean(tmpgeneTableGL$RMSE)
  tmpGof <- mean(tmpgeneTableGL$gof)
  geneTableGLDf2 <- rbind(geneTableGLDf2, c(tmpgeneTableGL$sample[1], tmpgeneTableGL$gamma[1],
                                          tmpRMSE, tmpGof))
}

geneTableGLDf2 <- data.frame(geneTableGLDf2)
colnames(geneTableGLDf2) <- c("sample", "gamma", "RMSE", "gof")
geneTableGLDf2[, 3:4] <- lapply(geneTableGLDf2[, 3:4], as.numeric)
geneTableGLDf2$ploidy <- as.character(sampleDf$ploidy2[match(geneTableGLDf2$sample, sampleDf$sample2)])

t <- ggplot(geneTableGLDf2, aes(x = gamma, y = RMSE, color = sample)) + 
  geom_point(alpha = 0.5) + ggtitle("RMSE vs gamma: 100kb bin absolute CN comparison: gains/losses") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5))

u <- ggplot(geneTableGLDf2, aes(x = gamma, y = gof, color = sample)) + 
  geom_point(alpha = 0.5) + ggtitle("Goodness of fit vs gamma") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(t, u, ncol = 1)


t2 <- ggplot(geneTableGLDf2, aes(x = gamma, y = RMSE, color = ploidy)) + 
  geom_point(alpha = 0.25) + ggtitle("RMSE vs gamma: gene absolute CN comparison: gains/losses") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0, 5)) 

u2 <- ggplot(geneTableGLDf2, aes(x = gamma, y = gof, color = ploidy)) + 
  geom_point(alpha = 0.25) + ggtitle("Goodness of fit vs gamma") + scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) + 
  theme(plot.title = element_text(hjust = 0.5)) 

grid.arrange(t2, u2, ncol = 1)

for (i in unique(geneTableGLDf2$gamma)) {
  tmp <- geneTableGLDf2[which(geneTableGLDf2$gamma == i), ]
  print(paste("gamma", tmp$gamma[1], "mean rmse", signif(median(tmp$RMSE), digits = 3),
              "sd rmse", signif(sd(tmp$RMSE), digits = 3),
              "gof mean", signif(mean(tmp$gof), digits = 3),
              "gof sd", signif(sd(tmp$gof), digits = 3)))
}

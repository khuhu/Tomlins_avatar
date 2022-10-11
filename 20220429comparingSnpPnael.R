library(stringr)
library(irr)

panelV3 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/snpVsPanel/20220425panelPloidyV3_comp.xlsx")
# snpRes05 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/snpVsPanel/20220425snpPloidyRes05.xlsx")
snpRes05 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/snpVsPanel/20220503snpRes05V2.xlsx")

panelV3$V3_ploidy_var <- round(panelV3$V3_ploidy_var, digits = 2)


eros <- c("/mnt/DATA6/mouseData/copynumber/")
listOfDirectories <- c("Auto_user_AUS5-138-MG_cho_20210621_354_343", "Auto_user_AUS5-76-MG_test1_255_185",
                       "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                       "Auto_user_AUS5-120-MG_EFD4_BBN_334_304")

allNoshadSeg <- read.table("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15/absoluteRes.txt", sep = "\t",
                           stringsAsFactors = FALSE, header = TRUE, fill = TRUE)

# allNoshadSeg <- read.table("/mnt/DATA5/tmp/kev/noshadSnpFiltHclustGapped15/absoluteRes.txt", sep = "\t",
#                            stringsAsFactors = FALSE, header = TRUE, fill = TRUE)

panelV3$sample[which(panelV3$sample == "2027lte")] <- "mg4"
panelV3$sample[which(panelV3$sample == "13085lt")] <- "mg15"
panelV3$sample[which(panelV3$sample == "14399rt")] <- "mg20"
panelV3$sample[which(panelV3$sample == "14656peritonealmt")] <- "14656peritnealmt"
panelV3$sample[which(panelV3$sample == "133576rt")] <- "13576rt"
panelV3$sample[which(panelV3$sample == "14150rt")] <- "14150lt"
panelV3$sample[which(panelV3$sample == "14154rot")] <- "14154lt"

allNoshadSeg_filt <- allNoshadSeg[grep(paste0(paste(panelV3$sample, panelV3$V3_purity, panelV3$V3_ploidy_var, sep = "_"), collapse = "|"), allNoshadSeg$sample),]
# allNoshadSeg_filt <- allNoshadSeg[grep(paste0(paste(panelV3$sample, panelV3$V4_purity, round(panelV3$V4_ploidy, 2), sep = "_"), collapse = "|"), allNoshadSeg$sample),]
allNoshadSeg_filt$sC <- as.numeric(allNoshadSeg_filt$sC)
# check
# paste(panelV3$sample, panelV3$V3_purity, panelV3$V3_ploidy_var, sep = "_")[-which(paste(panelV3$sample, panelV3$V3_purity, panelV3$V3_ploidy_var, sep = "_") %in% unique(allNoshadSeg_filt$sample))]


# allComps <- paste0(snpRes05$ploidy_int, ":", panelV3$V3)
# table(allComps)

# snp ratios diploid 17/22; triploid 15/20; tetraploid 2/4
# overall match rate is 74%

# 
# snpRes05$sample[which(allComps == "2:4")]
# snpRes05$sample[which(allComps == "4:2")]
# 
# snpRes05$sample[which(allComps == "2:3")]
# snpRes05$sample[which(allComps == "3:2")]
# snpRes05$sample[which(allComps == "3:4")]
# snpRes05$sample[which(allComps == "4:3")]
# 

# res <- NULL
# for (i in unique(allNoshadSeg_filt$sample)) {
#   tmpDf <- allNoshadSeg_filt[which(allNoshadSeg_filt$sample == i),]
#   subCloneCount <- length(which(abs(tmpDf$sC - round(tmpDf$sC)) > 0.2))
#   res <- c(res, subCloneCount)
# }
# 
# 
# matchIdx <- grep(paste0(panelV3$sample[which(panelV3$V3 == snpRes05$ploidy_int)], collapse = "|"), unique(allNoshadSeg_filt$sample))
# nonMatchIdx <- grep(paste0(panelV3$sample[-which(panelV3$V3 == snpRes05$ploidy_int)], collapse = "|"), unique(allNoshadSeg_filt$sample))
# t.test(x = res[matchIdx], y = res[nonMatchIdx])
# 

### can't find reason empirical explanation for discordant calls
### circle back after these next two analyses i.e resolution of calls i.e range of how large
### aneuploidy concordance .... think of proper comparison if ploidies are wrong?
### assume diploid for all 46 samples and then compare changes/aneuploidy loss or 4 ... since even the smallest changes?

### per sample look at gains and losses corr then overall corr

# allAscatSegs <- read.table("/mnt/DATA5/tmp/kev/misc/20220429AscatIllusRes05Seg100.txt", sep = "\t",
#                            stringsAsFactors = FALSE, header = TRUE)

allAscatSegs <- read.table("/mnt/DATA5/tmp/kev/misc/20220512AscatIllusRes05Seg100.txt", sep = "\t",
                           stringsAsFactors = FALSE, header = TRUE)

allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "14399RT", "14399rt")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "2027_LT-E", "2027lte")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "13085ROT\\(L\\)", "13085lt")

allAscatSegs_filt <- allAscatSegs[grep(paste0(paste0(snpRes05$sample, "_rho", snpRes05$purity_2, "_psi", snpRes05$ploidy_var_2), collapse = "|"), allAscatSegs$sample),]
allAscatSegs_filt[, 2:6] <- lapply(allAscatSegs_filt[, 2:6], as.numeric)

# paste0(snpRes05$sample, "_rho", snpRes05$purity, "_psi", snpRes05$ploidy_var)[-which(paste0(snpRes05$sample, "_rho", snpRes05$purity, "_psi", snpRes05$ploidy_var) %in% unique(allAscatSegs_filt$sample))]

# checking aneupoloidy status first

allAscatSegs_filt$totalSc <- as.numeric(allAscatSegs_filt$nMajor) + as.numeric(allAscatSegs_filt$nMinor)
allAscatSegs_filt$chr <- paste0("chr", allAscatSegs_filt$chr)
allAscatSegs_filt <-  allAscatSegs_filt[which(allAscatSegs_filt$chr %in% unique(allNoshadSeg_filt$chr)),]
allAscatSegs_filt$size <- allAscatSegs_filt$endpos - allAscatSegs_filt$startpos

allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "15774rt", "tmp")
allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "15774lt", "15774rt")
allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "tmp", "15774lt")


ascatAneploidy <- matrix(0, nrow = length(unique(snpRes05$sample)), ncol = length(unique(allNoshadSeg_filt$chr)))
rownames(ascatAneploidy) <- unique(snpRes05$sample)
colnames(ascatAneploidy) <- unique(allNoshadSeg_filt$chr)


### >0.20 cutoff from 2018 paper Allison M Taylor for arm-level stuff more or less doubled b/c mouse have no arms like chr14 in humans
for (i in rownames(ascatAneploidy)) {
  tmpDf <- allAscatSegs_filt[grep(i, allAscatSegs_filt$sample),]
  tmpPloidy <- snpRes05$ploidy_int[which(snpRes05$sample == i)]
  tmpDf$totalSc_n <- tmpDf$totalSc - tmpPloidy
  for (j in unique(tmpDf$chr)) {
    tmpChr <- tmpDf[which(tmpDf$chr == j), ]
    chrThres <- sum(tmpChr$size) * 0.6
    tmpAneu <- "none"
    if (sum(tmpChr$size[which(tmpChr$totalSc_n > 0)]) > chrThres) {
      tmpAneu <- "gain"
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < 0)]) > chrThres){
      tmpAneu <- "loss"
    } else{
      tmpAneu <- "none"
    }
    ascatAneploidy[i,j] <- tmpAneu
  }
}



noshadAneploidy <- matrix(0, nrow = length(unique(panelV3$sample)), ncol = length(unique(allNoshadSeg_filt$chr)))
rownames(noshadAneploidy) <- unique(panelV3$sample)
colnames(noshadAneploidy) <- unique(allNoshadSeg_filt$chr)
allNoshadSeg_filt$size <- as.numeric(allNoshadSeg_filt$end.off) - as.numeric(allNoshadSeg_filt$start.off)


### .8 cut off for things to be rounded to 
i <- rownames(noshadAneploidy)[2]
for (i in rownames(noshadAneploidy)) {
  tmpDf <- allNoshadSeg_filt[grep(i, allNoshadSeg_filt$sample),]
  tmpPloidy <- as.numeric(panelV3$V3[which(panelV3$sample == i)])
  tmpDf$totalSc_n <- tmpDf$sC- tmpPloidy
  tmpDf$totalSc_n[which(abs(tmpDf$totalSc_n) < 0.2)] <- 0
  j <- "chr12"
  for (j in unique(tmpDf$chr)) {
    tmpChr <- tmpDf[which(tmpDf$chr == j), ]
    chrThres <- sum(tmpChr$size) * 0.6
    tmpAneu <- "none"
    if (sum(tmpChr$size[which(tmpChr$totalSc_n > 0.8)]) > chrThres) {
      tmpAneu <- "gain"
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < -0.8)]) > chrThres){
      tmpAneu <- "loss"
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n > 0.2)]) > chrThres){
      tmpAneu <- "sub"
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < -0.2)]) > chrThres){
      tmpAneu <- "sub"
    } else{
      tmpAneu <- "none"
    }
    noshadAneploidy[i,j] <- tmpAneu
  }
}

### ask Scott about my assumptions for >= 0.8 being clonal but also 0.5 =< x < 0.8 being subclonal rounded to clonal


### look for overall comparion - replace sub with none


rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "2027lte")] <- "mg4"
rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "13085lt")] <- "mg15"
rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "14399rt")] <- "mg20"
rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "14656peritonealmt")] <- "14656peritnealmt"
rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "133576rt")] <- "13576rt"
rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "14150rt")] <- "14150lt"
rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "14154rot")] <- "14154lt"

clonalCompNoshad <- noshadAneploidy
clonalCompNoshad[clonalCompNoshad == "sub"] <- "none"


### general, gains, losses
### then subclonal respectively

# cat_to_cont <- function(df){
#   df <- str_replace_all(df, "none", "0")
#   df <- str_replace_all(df, "gain", "1")
#   df <- str_replace_all(df, "loss", "-1")
#   df <- as.numeric(df)
#   return(df)
# }


mcConcord <- NULL
ckappConcord <- NULL
concordList <- NULL
corList <- NULL
i <- rownames(ascatAneploidy)[2]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- clonalCompNoshad[i, 1:19]
  concord <- length(which(tmpAscat == tmpNoshad))/19
  concordList <- c(concordList, concord)
  
  tmpKapp <- kappa2(cbind(tmpAscat, tmpNoshad))
  tmpMc <- tryCatch(mcnemar.test(tmpAscat, tmpNoshad), 
                  error = function(e) NULL)
  ckappConcord <- c(ckappConcord, tmpKapp$value)
  if (is.null(tmpMc)) {
    mcConcord <- c(mcConcord, NA)
  } else{
    mcConcord <- c(mcConcord, tmpMc$p.value)
  }
}



mcConcordGains <- NULL
ckappConcordGains <- NULL
concordListGains <- NULL
corGains <- NULL
i <- rownames(ascatAneploidy)[2]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- clonalCompNoshad[i, 1:19]
  
  tmpIdx <- c(which(tmpAscat == "gain"), which(tmpNoshad == "gain"))
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListGains <- c(concordListGains, concord)
    corGains <- c(corGains, NA)
    mcConcordGains <- c(mcConcordGains, NA)
    ckappConcordGains <- c(ckappConcordGains, NA)
    
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListGains <- c(concordListGains, concord)
    
    tmpKapp <- kappa2(cbind(tmpAscat, tmpNoshad))
    ckappConcordGains <- c(ckappConcordGains, tmpKapp$value)
    
    tmpMc <- tryCatch(mcnemar.test(tmpAscat, tmpNoshad), 
                      error = function(e) NULL)
    if (is.null(tmpMc)) {
      mcConcordGains <- c(mcConcordGains, NA)
    } else{
      mcConcordGains <- c(mcConcordGains, tmpMc$p.value)
    }
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListGains <- c(concordListGains, concord)
    
    tmpKapp <- kappa2(cbind(tmpAscat, tmpNoshad))
    ckappConcordGains <- c(ckappConcordGains, tmpKapp$value)
    
    tmpMc <- tryCatch(mcnemar.test(tmpAscat, tmpNoshad), 
                      error = function(e) NULL)
    if (is.null(tmpMc)) {
      mcConcordGains <- c(mcConcordGains, NA)
    } else{
      mcConcordGains <- c(mcConcordGains, tmpMc$p.value)
    }
  }
}


ckappConcordLosses <- NULL
concordListLosses <- NULL
corLosses <- NULL
i <- rownames(ascatAneploidy)[1]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- clonalCompNoshad[i, 1:19]
  
  
  tmpIdx <- c(which(tmpAscat == "loss"), which(tmpNoshad == "loss"))
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListLosses <- c(concordListLosses, concord)
    corLosses <- c(corLosses, NA)
    mcConcordLosses <- c(mcConcordLosses, NA)
    ckappConcordLosses <- c(ckappConcordLosses, NA)
    
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListLosses <- c(concordListLosses, concord)
    
    tmpKapp <- kappa2(cbind(tmpAscat, tmpNoshad))
    ckappConcordLosses <- c(ckappConcordLosses, tmpKapp$value)
    

  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListLosses <- c(concordListLosses, concord)
    
    tmpKapp <- kappa2(cbind(tmpAscat, tmpNoshad))
    ckappConcordLosses <- c(ckappConcordLosses, tmpKapp$value)
  }
}

### same thing but with subclones

ckappConcordSub <- NULL
concordListSub <- NULL
i <- rownames(ascatAneploidy)[2]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- noshadAneploidy[i, 1:19]
  
  
  if (length(which(tmpNoshad == "sub")) > 0) {
    # tmpAscat <- tmpAscat[-which(tmpNoshad == "sub")]
    # tmpNoshad <- tmpNoshad[which(tmpNoshad == "sub")]
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  concord <- length(which(tmpAscat == tmpNoshad))/19
  concordListSub <- c(concordListSub, concord)
  tmpKapp <- kappa2(cbind(tmpAscat, tmpNoshad))
  ckappConcordSub <- c(ckappConcordSub, tmpKapp$value)
}


tmp <- NULL
concordListSubGains <- NULL
i <- rownames(ascatAneploidy)[16]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- noshadAneploidy[i, 1:19]
  
  if (length(which(tmpNoshad == "sub")) > 0) {
    # tmpAscat <- tmpAscat[-which(tmpNoshad == "sub")]
    # tmpNoshad <- tmpNoshad[which(tmpNoshad == "sub")]
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  
  tmpIdx <- c(which(tmpAscat == "gain"), which(tmpNoshad == "gain"))
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListSubGains <- c(concordListSubGains, concord)
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListSubGains  <- c(concordListSubGains , concord)
    
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListSubGains  <- c(concordListSubGains , concord)
  }
  
  tmp <- c(tmp, paste(i, concord))
}


concordListSubLosses <- NULL
i <- rownames(ascatAneploidy)[1]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- noshadAneploidy[i, 1:19]
  
  if (length(which(tmpNoshad == "sub")) > 0) {
    # tmpAscat <- tmpAscat[-which(tmpNoshad == "sub")]
    # tmpNoshad <- tmpNoshad[which(tmpNoshad == "sub")]
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  tmpIdx <- c(which(tmpAscat == "loss"), which(tmpNoshad == "loss"))
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListSubLosses <- c(concordListSubLosses, concord)  
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListSubLosses <- c(concordListSubLosses, concord)
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListSubLosses <- c(concordListSubLosses, concord)
  }
}


### check if factoring correct ploidy changes anything

matchedAscat <- ascatAneploidy[which(snpRes05$ploidy_int == panelV3$V3),]
matchedNoshad <- noshadAneploidy[which(snpRes05$ploidy_int == panelV3$V3),]

clonalCompNoshadM <- matchedNoshad
clonalCompNoshadM[clonalCompNoshadM == "sub"] <- "none"


concordListMatched <- NULL
for (i in rownames(matchedAscat)) {
  tmpAscat <- matchedAscat[i, 1:19]
  tmpNoshad <- clonalCompNoshadM[i, 1:19]
  concord <- length(which(tmpAscat == tmpNoshad))/19
  concordListMatched <- c(concordListMatched, concord)
}


concordListMatchedGains <- NULL
for (i in rownames(matchedAscat)) {
  tmpAscat <- matchedAscat[i, 1:19]
  tmpNoshad <- clonalCompNoshadM[i, 1:19]
  
  tmpIdx <- c(which(tmpAscat == "gain"), which(tmpNoshad == "gain"))
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListMatchedGains <- c(concordListMatchedGains, concord)  
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListMatchedGains <- c(concordListMatchedGains, concord)
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListMatchedGains <- c(concordListMatchedGains, concord)
  }
}


concordListMatchedLosses <- NULL
i <- rownames(matchedAscat)[2]
for (i in rownames(matchedAscat)) {
  tmpAscat <- matchedAscat[i, 1:19]
  tmpNoshad <- clonalCompNoshadM[i, 1:19]
  
  
  tmpIdx <- c(which(tmpAscat == "loss"), which(tmpNoshad == "loss"))
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListMatchedLosses <- c(concordListMatchedLosses, concord)  
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListMatchedLosses <- c(concordListMatchedLosses, concord)
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListMatchedLosses <- c(concordListMatchedLosses, concord)
  }
}



concordListMatchedSub <- NULL
for (i in rownames(matchedAscat)) {
  tmpAscat <- matchedAscat[i, 1:19]
  tmpNoshad <- matchedNoshad[i, 1:19]
  
  if (length(which(tmpNoshad == "sub")) > 0) {
    # tmpAscat <- tmpAscat[-which(tmpNoshad == "sub")]
    # tmpNoshad <- tmpNoshad[which(tmpNoshad == "sub")]
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  concord <- length(which(tmpAscat == tmpNoshad))/19
  concordListMatchedSub <- c(concordListMatchedSub, concord)
}


tmp2 <- NULL
concordListMatchedSubGains <- NULL
i <- rownames(matchedAscat)[9]
for (i in rownames(matchedAscat)) {
  tmpAscat <- matchedAscat[i, 1:19]
  tmpNoshad <- matchedNoshad[i, 1:19]
  
  if (length(which(tmpNoshad == "sub")) > 0) {
    # tmpAscat <- tmpAscat[-which(tmpNoshad == "sub")]
    # tmpNoshad <- tmpNoshad[which(tmpNoshad == "sub")]
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  tmpIdx <- c(which(tmpAscat == "gain"), which(tmpNoshad == "gain"))
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListMatchedSubGains <- c(concordListMatchedSubGains, concord)  
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListMatchedSubGains <- c(concordListMatchedSubGains, concord)
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListMatchedSubGains <- c(concordListMatchedSubGains, concord)
  }
  tmp2 <- c(tmp2, paste(i, concord))
}


concordListMatchedSubLosses <- NULL
for (i in rownames(matchedAscat)) {
  tmpAscat <- matchedAscat[i, 1:19]
  tmpNoshad <- matchedNoshad[i, 1:19]
  
  if (length(which(tmpNoshad == "sub")) > 0) {
    # tmpAscat <- tmpAscat[-which(tmpNoshad == "sub")]
    # tmpNoshad <- tmpNoshad[which(tmpNoshad == "sub")]
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  tmpIdx <- c(which(tmpAscat == "loss"), which(tmpNoshad == "loss"))
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListMatchedSubLosses <- c(concordListMatchedSubLosses, concord)  
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListMatchedSubLosses <- c(concordListMatchedSubLosses, concord)
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- length(which(tmpAscat[tmpIdx] == tmpNoshad[tmpIdx]))/length(tmpIdx)
    concordListMatchedSubLosses <- c(concordListMatchedSubLosses, concord)
  }
}


### doesn't seem to matter - check to see which and why it happens to have no concordance
### tumor content?  i.e difference in the matched TC for most almost large difference of 10% for some
### 15774 lt and rt swap for one set


#### more than half are below 50% concordance ... so look at those samples specifically. could be aneuploidy cutoffs used
matchedNoshadDf <- panelV3[match(rownames(matchedNoshad), panelV3$sample), ]
matchedNoshadDf$matchedSubLosses <- concordListMatchedSubLosses

# wanted to see if the two are in fact swapped 15774 lt and rt
panelV3$sublos <- concordListSubLosses
panelV3$subgain <- concordListSubGains


summary(concordListLosses)
summary(concordListSubLosses)
summary(concordListMatchedLosses)
summary(concordListMatchedSubLosses)

summary(concordListGains)
summary(concordListSubGains)
summary(concordListMatchedGains)
summary(concordListMatchedSubGains)

plot(matchedNoshadDf$V3, concordListMatchedSubLosses)

lossesCount <- apply(matchedAscat, 1, function(x) length(which(x == "loss")))
gainsCount <- apply(matchedAscat, 1, function(x) length(which(x == "gain")))


dev.off()
png(filename = "/mnt/DATA5/tmp/kev/misc/20220509altCountVFactors.png", width = 1300, height = 1000)
par(mfrow=c(2,2))
plot(lossesCount, concordListMatchedSubLosses, main = "# of losses v. concord")
plot(gainsCount, concordListMatchedSubGains, main = "# gain v. concord")
plot(matchedNoshadDf$V3_purity, concordListMatchedSubLosses, main = "purity v. cocord loss")
plot(matchedNoshadDf$V3_purity, concordListMatchedSubGains, main = "purity v. concord gain")
dev.off()

cor(lossesCount[-which(lossesCount == 0)],  concordListMatchedSubLosses[-which(lossesCount == 0)])
cor(gainsCount[-which(gainsCount == 0)],  concordListMatchedSubGains[-which(gainsCount == 0)])
cor(matchedNoshadDf$V3_purity[-which(is.na(concordListMatchedSubLosses))],
    concordListMatchedSubLosses[-which(is.na(concordListMatchedSubLosses))])
cor(matchedNoshadDf$V3_purity[-which(is.na(concordListMatchedSubGains))],
    concordListMatchedSubGains[-which(is.na(concordListMatchedSubGains))])


matchedPloidy <- as.numeric(panelV3$V3[which(panelV3$V3 == snpRes05$ploidy_int)])
concordListMatchedSubLosses_euploid <- concordListMatchedSubLosses[which(matchedPloidy == 2)]
concordListMatchedSubLosses_polyploid <- concordListMatchedSubLosses[which(matchedPloidy != 2)]
concordListMatchedSubGains_euploid <- concordListMatchedSubGains[which(matchedPloidy == 2)]
concordListMatchedSubGains_polyploid <- concordListMatchedSubGains[which(matchedPloidy != 2)]



summary(concordListMatchedSubLosses_euploid)
summary(concordListMatchedSubLosses_polyploid)
summary(concordListMatchedSubGains_euploid)
summary(concordListMatchedSubGains_polyploid)

alterations <- c(rep("euploidLoss", length(concordListMatchedSubLosses_euploid)), rep("euploidGain", length(concordListMatchedSubGains_euploid)),
                 rep("polyploidLoss", length(concordListMatchedSubLosses_polyploid)), rep("polyploidGain", length(concordListMatchedSubGains_polyploid)))
concordance <- c(concordListMatchedSubLosses_euploid, concordListMatchedSubGains_euploid, 
                 concordListMatchedSubLosses_polyploid, concordListMatchedSubGains_polyploid)

ggplotDf <- data.frame("alterations" = alterations, "concordance" = concordance)
library(ggplot2)

dev.off()
png(filename = "/mnt/DATA5/tmp/kev/misc/20220509altPloidySep.png", width = 1300, height = 1000)
ggplot(data = ggplotDf) + geom_boxplot(aes(x = alterations, y  = concordance), color = c("red", "blue", "red", "blue"))
dev.off()


alterations_raw <- c(rep("loss", length(concordListLosses)), rep("gain", length(concordListGains)))
concordance_raw <- c(concordListLosses, concordListGains)
ggplotDf_raw <- data.frame("alterations" = alterations_raw, "concordance" = concordance_raw)

dev.off()
png(filename = "/mnt/DATA5/tmp/kev/misc/20220509altPloidySepRaw.png", width = 1300, height = 1000)
ggplot(data = ggplotDf_raw) + geom_boxplot(aes(x = alterations, y  = concordance), color = c("red", "blue"))
dev.off()


alterations_sub <- c(rep("loss", length(concordListSubLosses)), rep("gain", length(concordListSubGains)))
concordance_sub <- c(concordListSubLosses, concordListSubGains)
ggplotDf_sub <- data.frame("alterations" = alterations_sub, "concordance" = concordance_sub)


dev.off()
png(filename = "/mnt/DATA5/tmp/kev/misc/20220509altPloidySepSub.png", width = 1300, height = 1000)
ggplot(data = ggplotDf_sub) + geom_boxplot(aes(x = alterations, y  = concordance), color = c("red", "blue"))
dev.off()


cListSubLosses_m_e <-  concordListSubLosses[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 2)]
cListSubLosses_m_p <-  concordListSubLosses[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 != 2)]
cListSubLosses_um_e <-  concordListSubLosses[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 2)]
cListSubLosses_um_p <-  concordListSubLosses[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 != 2)]

cListSubGains_m_e <-  concordListSubGains[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 2)]
cListSubGains_m_p <-  concordListSubGains[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 != 2)]
cListSubGains_um_e <-  concordListSubGains[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 2)]
cListSubGains_um_p <-  concordListSubGains[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 != 2)]


alterations <- c(rep("loss_m_e", length(cListSubLosses_m_e)), rep("loss_m_p", length(cListSubLosses_m_p)),
                 rep("loss_um_e", length(cListSubLosses_um_e)), rep("loss_um_p", length(cListSubLosses_um_p)),
                 rep("gain_m_e", length(cListSubGains_m_e)), rep("gain_m_p", length(cListSubGains_m_p)),
                 rep("gain_um_e", length(cListSubGains_um_e)), rep("gain_um_p", length(cListSubGains_um_p)))

concordance <- c(cListSubLosses_m_e, cListSubLosses_m_p, cListSubLosses_um_e, cListSubLosses_um_p,
                 cListSubGains_m_e, cListSubGains_m_p, cListSubGains_um_e, cListSubGains_um_p)

ggplotDf_m_p <- data.frame("alterations" = alterations, "concordance" = concordance)
ggplotDf_m_p$alterations <- factor(ggplotDf_m_p$alterations, levels = unique(ggplotDf_m_p$alterations))


n_fun <- function(x){
  return(data.frame(y = 1.15, label = paste0("n = ",length(x))))
}

ggplot(data = ggplotDf_m_p, aes(x = alterations, y  = concordance)) + geom_boxplot(color = c(rep("blue", 4), rep("red", 4))) + 
  stat_summary(fun.data = n_fun, geom = "text") + scale_y_continuous(breaks = seq(0, 1.2, by = 0.1))


### perfect situation

cListSubLosses_m_e <-  concordListSubLosses[which(snpRes05$ploidy_int == panelV3$V4_matched & panelV3$V4_matched == 2)]
cListSubLosses_m_p <-  concordListSubLosses[which(snpRes05$ploidy_int ==  panelV3$V4_matched & panelV3$V4_matched != 2)]

cListSubGains_m_e <-  concordListSubGains[which(snpRes05$ploidy_int == panelV3$V4_matched & panelV3$V4_matched == 2)]
cListSubGains_m_p <-  concordListSubGains[which(snpRes05$ploidy_int == panelV3$V4_matched & panelV3$V4_matched != 2)]



alterations <- c(rep("loss_m_e", length(cListSubLosses_m_e)), rep("loss_m_p", length(cListSubLosses_m_p)),
                 rep("gain_m_e", length(cListSubGains_m_e)), rep("gain_m_p", length(cListSubGains_m_p)))

concordance <- c(cListSubLosses_m_e, cListSubLosses_m_p, 
                 cListSubGains_m_e, cListSubGains_m_p)

ggplotDf_m_p <- data.frame("alterations" = alterations, "concordance" = concordance)
ggplotDf_m_p$alterations <- factor(ggplotDf_m_p$alterations, levels = unique(ggplotDf_m_p$alterations))


n_fun <- function(x){
  return(data.frame(y = 1.15, label = paste0("n = ",length(x))))
}

dev.off()
pdf(file = "/mnt/DATA5/tmp/kev/snpVsPanel/20220510pseudoMatchBoxplot.pdf")
ggplot(data = ggplotDf_m_p, aes(x = alterations, y  = concordance)) + geom_boxplot(color = c(rep("blue", 2), rep("red", 2))) + 
  stat_summary(fun.data = n_fun, geom = "text") + scale_y_continuous(breaks = seq(0, 1.2, by = 0.1))
dev.off()


### could be strictly triploid that's off or and WGD - why WGD? signal too low i.e filtering out CNR < 0.2 ?
### graph accuracy by purity
### possibly recall ploidy and change lrs based on missing regions - can use knownSeg files from PSCBS
### gain and lost count per sample have no  correlation - there are no real measurements of magnitude 
### i.e it's just gain, loss and none- maybe table both vars and then 



### specific samples that are just bad? switch the swapped samples



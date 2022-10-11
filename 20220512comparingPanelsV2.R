library(stringr)
library(irr)

# panelV3 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/snpVsPanel/20220425panelPloidyV3_comp.xlsx")
# snpRes05 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/snpVsPanel/20220425snpPloidyRes05.xlsx")
# snpRes05 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/snpVsPanel/20220503snpRes05V2.xlsx")
snpRes05 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/snpVsPanel/20220503snpRes05V3.xlsx")

# snpRes05$purity_2[which(snpRes05$sample == "2611lt")] <- 0.53
# snpRes05$ploidy_var_2[which(snpRes05$sample == "2611lt")] <- 2
# 
# snpRes05$ploidy_int[which(snpRes05$sample == "14120rt")] <- 2
# snpRes05$purity_2[which(snpRes05$sample == "14120rt")] <- 0.52
# snpRes05$ploidy_var_2[which(snpRes05$sample == "14120rt")] <- 2.4
# 
# snpRes05$purity_2[which(snpRes05$sample == "15774lt")] <- 0.77
# snpRes05$ploidy_var_2[which(snpRes05$sample == "15774lt")] <- 2.8
# 
# snpRes05$purity_2[which(snpRes05$sample == "15774rt")] <- 0.68
# snpRes05$ploidy_var_2[which(snpRes05$sample == "15774rt")] <- 2.8
# 
# snpRes05$ploidy_int[which(snpRes05$sample == "14656peritnealmt")] <- 2


# panelV3$V3_ploidy_var <- round(panelV3$V3_ploidy_var, digits = 2)

# allPloidyCalls <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220517allHgscPloidy.xlsx")
# allPloidyCalls <- allPloidyCalls[-which(allPloidyCalls$ploidy_int == "NA"), ]

# allPloidyCallsV2 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220701pbsGlobalSmoothSnpComp.xlsx")
# allPloidyCalls$ploidy_int[1:48] <- allPloidyCallsV2$ploidy_int
# allPloidyCalls$ploidyV2[1:48] <- allPloidyCallsV2$ploidy
# allPloidyCalls$purityV2[1:48] <- allPloidyCallsV2$purity

allPloidyCalls <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220517allHgscPloidySnpComp.xlsx")


panelV3 <- allPloidyCalls
# panelV3 <- panelV3[-which(duplicated(panelV3$sample)), ]
# panelV3 <- panelV3[-which(panelV3$sample %in% c("mg4", "mg15", "mg20")), ]
colnames(panelV3)[2] <- "V3"


# panelV3$purityV2[which(panelV3$sample == "2611lt")] <- 0.5
# panelV3$ploidyV2[which(panelV3$sample == "2611lt")] <- 2.02
# 
# panelV3$purityV2[which(panelV3$sample == "14120rt")] <- 0.545
# panelV3$ploidyV2[which(panelV3$sample == "14120rt")] <- 2.02

#second 2285lt and 2007rt and 14943lt"

eros <- c("/mnt/DATA6/mouseData/copynumber/")
listOfDirectories <- c("Auto_user_AUS5-138-MG_cho_20210621_354_343", "Auto_user_AUS5-76-MG_test1_255_185",
                       "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                       "Auto_user_AUS5-120-MG_EFD4_BBN_334_304")

# allNoshadSeg <- read.table("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV2/absoluteRes.txt", sep = "\t",
#                            stringsAsFactors = FALSE, header = TRUE, fill = TRUE)


# allNoshadSeg <- read.table("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genesG/absoluteRes.txt", sep = "\t",
#                            stringsAsFactors = FALSE, header = TRUE, fill = TRUE)

# allNoshadSeg <- read.table("/mnt/DATA5/tmp/kev/noshadSnpFiltHclustGapped15/absoluteRes.txt", sep = "\t",
#                            stringsAsFactors = FALSE, header = TRUE, fill = TRUE)

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


# allNoshadSeg_filt <- allNoshadSeg[grep(paste0(paste(panelV3$sample, panelV3$V3_purity, panelV3$V3_ploidy_var, sep = "_"), collapse = "|"), allNoshadSeg$sample),]
# allNoshadSeg_filt <- allNoshadSeg[grep(paste0(paste(panelV3$sample, panelV3$V4_purity, round(panelV3$V4_ploidy, 2), sep = "_"), collapse = "|"), allNoshadSeg$sample),]
# allNoshadSeg_filt$sC <- as.numeric(allNoshadSeg_filt$sC)

allNoshadSeg_filt <- allNoshadSeg[grep(paste0(paste(panelV3$sample, sprintf("%0.2f", panelV3$purityV2),
                                                    panelV3$ploidyV2, sep = "_"), collapse = "|"), allNoshadSeg$sample),]
allNoshadSeg_filt$sC <- as.numeric(allNoshadSeg_filt$sC)
allNoshadSeg_filt[,2:3] <- lapply(allNoshadSeg_filt[,2:3], as.numeric)


### per sample look at gains and losses corr then overall corr

# allAscatSegs <- read.table("/mnt/DATA5/tmp/kev/misc/20220429AscatIllusRes05Seg100.txt", sep = "\t",
#                            stringsAsFactors = FALSE, header = TRUE)

allAscatSegs <- read.table("/mnt/DATA5/tmp/kev/misc/20220512AscatIllusRes05Seg100.txt", sep = "\t",
                           stringsAsFactors = FALSE, header = TRUE)

allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "14399RT", "14399rt")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "2027_LT-E", "2027lte")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "13085ROT\\(L\\)", "13085lt")

# paste0(snpRes05$sample, "_rho", snpRes05$purity, "_psi", snpRes05$ploidy_var)[-which(paste0(snpRes05$sample, "_rho", snpRes05$purity, "_psi", snpRes05$ploidy_var) %in% unique(allAscatSegs_filt$sample))]

# checking aneupoloidy status first


allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "15774rt", "tmp")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "15774lt", "15774rt")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "tmp", "15774lt")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "2027lte",  "mg4")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "13085lt",  "mg15")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "14399rt",  "mg20")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "14656peritonealmt",  "14656peritnealmt")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "133576rt",  "13576rt")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "14150rt",  "14150lt")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "14154rot",  "14154lt")



allAscatSegs_filt <- allAscatSegs[which(allAscatSegs$sample %in% paste0(snpRes05$sample, "_rho", snpRes05$purity_2, "_psi", snpRes05$ploidy_var_2)),]
allAscatSegs_filt[, 2:6] <- lapply(allAscatSegs_filt[, 2:6], as.numeric)
allAscatSegs_filt$totalSc <- as.numeric(allAscatSegs_filt$nMajor) + as.numeric(allAscatSegs_filt$nMinor)
allAscatSegs_filt$chr <- paste0("chr", allAscatSegs_filt$chr)
allAscatSegs_filt <-  allAscatSegs_filt[which(allAscatSegs_filt$chr %in% unique(allNoshadSeg_filt$chr)),]
allAscatSegs_filt$size <- allAscatSegs_filt$endpos - allAscatSegs_filt$startpos


ascatAneploidy <- matrix(0, nrow = length(unique(snpRes05$sample)), ncol = length(unique(allNoshadSeg_filt$chr)))
rownames(ascatAneploidy) <- unique(snpRes05$sample)
colnames(ascatAneploidy) <- unique(allNoshadSeg_filt$chr)

i <- rownames(ascatAneploidy)[4]
### >0.20 cutoff from 2018 paper Allison M Taylor for arm-level stuff more or less doubled b/c mouse have no arms like chr14 in humans
for (i in rownames(ascatAneploidy)) {
  print(i)
  tmpDf <- allAscatSegs_filt[grep(i, allAscatSegs_filt$sample),]
  tmpPloidy <- snpRes05$ploidy_int[which(snpRes05$sample == i)]
  tmpDf$totalSc_n <- tmpDf$totalSc - tmpPloidy
  for (j in unique(tmpDf$chr)) {
    tmpChr <- tmpDf[which(tmpDf$chr == j), ]
    chrThres <- sum(tmpChr$size) * 0.7
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
  # tmpDf$totalSc_n[which(abs(tmpDf$totalSc_n) < 0.2)] <- 0
  j <- "chr12"
  for (j in unique(tmpDf$chr)) {
    tmpChr <- tmpDf[which(tmpDf$chr == j), ]
    chrThres <- sum(tmpChr$size) * 0.7
    tmpAneu <- "none"
    if (sum(tmpChr$size[which(tmpChr$totalSc_n > 0.8)]) > chrThres) {
      tmpAneu <- "gain"
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < -0.8)]) > chrThres){
      tmpAneu <- "loss"
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < 0.8 & tmpChr$totalSc_n > 0.2)]) > chrThres){
      tmpAneu <- "sub"
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n > -0.8 & tmpChr$totalSc_n < -0.2)]) > chrThres){
      tmpAneu <- "sub"
    } else{
      tmpAneu <- "none"
    }
    noshadAneploidy[i,j] <- tmpAneu
  }
}


### ask Scott about my assumptions for >= 0.8 being clonal but also 0.5 =< x < 0.8 being subclonal rounded to clonal


### look for overall comparion - replace sub with none


# rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "2027lte")] <- "mg4"
# rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "13085lt")] <- "mg15"
# rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "14399rt")] <- "mg20"
# rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "14656peritonealmt")] <- "14656peritnealmt"
# rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "133576rt")] <- "13576rt"
# rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "14150rt")] <- "14150lt"
# rownames(ascatAneploidy)[which(rownames(ascatAneploidy) == "14154rot")] <- "14154lt"


# badIdx <- which(rownames(noshadAneploidy) %in% c("2611lt", "14120rt", "15774lt"))
# noshadAneploidy <- noshadAneploidy[-badIdx,]
# ascatAneploidy <- ascatAneploidy[-badIdx,]

clonalCompNoshad <- noshadAneploidy
clonalCompNoshad[clonalCompNoshad == "sub"] <- "none"

### changing how I do these statistics based on proportion of all - since data fashioned before only works for agreement statistic
### kappa should work only if looking at overall concordance i.e losses, gains and no change

cat_to_cont <- function(df){
  df <- str_replace_all(df, "none", "0")
  df <- str_replace_all(df, "gain", "1")
  df <- str_replace_all(df, "loss", "-1")
  df <- as.numeric(df)
  return(df)
}


ckappConcord <- NULL
concordList <- NULL
i <- rownames(ascatAneploidy)[2]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- clonalCompNoshad[i, 1:19]
  
  concord <- agree(cbind(cat_to_cont(tmpAscat), cat_to_cont(tmpNoshad)))
  concordList <- c(concordList, concord$value)
  
  tmpKapp <- kappa2(cbind(tmpAscat, tmpNoshad))
  ckappConcord <- c(ckappConcord, tmpKapp$value)
}

kapGainVec <- NULL
concordListGains <- NULL
i <- rownames(ascatAneploidy)[2]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- clonalCompNoshad[i, 1:19]
  
  tmpIdx <- c(which(tmpAscat == "gain"), which(tmpNoshad == "gain"))
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListGains <- c(concordListGains, concord)
    
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListGains <- c(concordListGains, concord$value)
    kapGainVec <- rbind(kapGainVec, cbind(tmpAscat[tmpIdx],tmpNoshad[tmpIdx]))
    
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListGains <- c(concordListGains, concord$value)
    kapGainVec <- rbind(kapGainVec, cbind(tmpAscat[tmpIdx],tmpNoshad[tmpIdx]))
  }
}

kapLossVec <- NULL
concordListLosses <- NULL
i <- rownames(ascatAneploidy)[1]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- clonalCompNoshad[i, 1:19]
  
  
  tmpIdx <- c(which(tmpAscat == "loss"), which(tmpNoshad == "loss"))
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListLosses <- c(concordListLosses, concord)
    
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListLosses <- c(concordListLosses, concord$value)
    kapLossVec <- rbind(kapLossVec, cbind(tmpAscat[tmpIdx],tmpNoshad[tmpIdx]))
    
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListLosses <- c(concordListLosses, concord$value)
    kapLossVec <- rbind(kapLossVec, cbind(tmpAscat[tmpIdx],tmpNoshad[tmpIdx]))
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
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  concord <- agree(cbind(cat_to_cont(tmpAscat), cat_to_cont(tmpNoshad)))
  concordListSub <- c(concordListSub, concord$value)
  
  tmpKapp <- kappa2(cbind(tmpAscat, tmpNoshad))
  ckappConcordSub <- c(ckappConcordSub, tmpKapp$value)
}



kapGainSubVec <- NULL
concordListSubGains <- NULL
i <- rownames(ascatAneploidy)[16]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- noshadAneploidy[i, 1:19]
  
  if (length(which(tmpNoshad == "sub")) > 0) {
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  
  tmpIdx <- c(which(tmpAscat == "gain"), which(tmpNoshad == "gain"))
  if (length(tmpIdx) < 1) {
    concordListSubGains <- c(concordListSubGains, NA)
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListSubGains  <- c(concordListSubGains , concord$value)
    
    kapGainSubVec <- rbind(kapGainSubVec, cbind(tmpAscat[tmpIdx],tmpNoshad[tmpIdx]))
    
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListSubGains  <- c(concordListSubGains , concord$value)
    
    kapGainSubVec <- rbind(kapGainSubVec, cbind(tmpAscat[tmpIdx],tmpNoshad[tmpIdx]))
  }
}

kapLossSubVec <- NULL
concordListSubLosses <- NULL
i <- rownames(ascatAneploidy)[41]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- noshadAneploidy[i, 1:19]
  
  if (length(which(tmpNoshad == "sub")) > 0) {
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  tmpIdx <- c(which(tmpAscat == "loss"), which(tmpNoshad == "loss"))
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListSubLosses <- c(concordListSubLosses, NA)  
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListSubLosses <- c(concordListSubLosses, concord$value)
    
    kapLossSubVec <- rbind(kapLossSubVec, cbind(tmpAscat[tmpIdx],tmpNoshad[tmpIdx]))
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListSubLosses <- c(concordListSubLosses, concord$value)
    
    kapLossSubVec <- rbind(kapLossSubVec, cbind(tmpAscat[tmpIdx],tmpNoshad[tmpIdx]))
  }
}
# 
# kappa2(kapGainVec)
# kappa2(kapGainSubVec)
# kappa2(kapLossVec)
# kappa2(kapLossSubVec)

### separate kapp and mcnemar test on each of the changes individually

noshadAneploidySub <- noshadAneploidy
noshadAneploidySub[noshadAneploidySub == "sub"] <- ascatAneploidy[noshadAneploidySub == "sub"]

noshadAneploidyClonal <- noshadAneploidy
noshadAneploidyClonal[noshadAneploidyClonal == "sub"] <- "none"

badIdx <- which(rownames(noshadAneploidy) %in% c("2611lt", "14120rt", "15774lt"))

### do something similar to what i did before, but not separting it by sample


# table(ascatAneploidy[-which(rownames(noshadAneploidy) %in% c("2611lt", "14120rt", "15774lt"))])
# table(noshadAneploidy[-which(rownames(noshadAneploidy) %in% c("2611lt", "14120rt", "15774lt"))])
# table(noshadAneploidySub[-which(rownames(noshadAneploidy) %in% c("2611lt", "14120rt", "15774lt"))])
# 
# tmpAscatCount <- c(81, 213)
# tmpNoshadCount <- c(61, 128)
# tmpNoshadCountSub <- c(84, 201)
# 
# 
# mcnemar.test(as.matrix(cbind(tmpAscatCount, tmpNoshadCount)))
# mcnemar.test(as.matrix(cbind(tmpAscatCount, tmpNoshadCountSub)))
# 


ckappConcordSub_m_e <-  ckappConcordSub[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 2)]
ckappConcordSub_m_p <-  ckappConcordSub[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 != 2)]
ckappConcordSub_um_e <-  ckappConcordSub[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 2)]
ckappConcordSub_um_p <-  ckappConcordSub[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 != 2)]

alterations <- c(rep("ckapp_m_e", length(ckappConcordSub_m_e)), rep("ckapp_m_p", length(ckappConcordSub_m_p)),
                 rep("ckapp_um_e", length(ckappConcordSub_um_e)), rep("ckapp_um_p", length(ckappConcordSub_um_p)))

concordance <- c(ckappConcordSub_m_e, ckappConcordSub_m_p, ckappConcordSub_um_e, ckappConcordSub_um_p)

ggplotDf_ck <- data.frame("alterations" = alterations, "ckap" = concordance)
ggplotDf_ck$alterations <- factor(ggplotDf_ck$alterations, levels = unique(ggplotDf_ck$alterations))


n_fun <- function(x){
  return(data.frame(y = 1.15, label = paste0("n = ",length(x))))
}


# ck_graph <- ggplot(data = ggplotDf_ck, aes(x = alterations, y  = ckap)) + geom_boxplot() + 
#   geom_jitter(color="black", size= 1 , alpha=0.9) +
#   stat_summary(fun.data = n_fun, geom = "text") + scale_y_continuous(breaks = seq(0, 1.2, by = 0.1)) +  
#   ggtitle("cohen's k: 0.8 length 0.8 absCn") + theme(plot.title = element_text(hjust = 0.5)) + 
#   xlab("match + ploidy")


ck_graph <- ggplot(data = ggplotDf_ck, aes(x = alterations, y  = ckap)) + geom_boxplot() + 
  geom_jitter(color="black", size= 1 , alpha=0.9) +
  stat_summary(fun.data = n_fun, geom = "text") + ylim(c(0,1.2)) + 
  ggtitle("cohen's k: 0.8 length 0.8 absCn") + theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("match + ploidy")



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

# pa_graph <- ggplot(data = ggplotDf_m_p, aes(x = alterations, y  = concordance/100)) + geom_boxplot(color = c(rep("blue", 4), rep("red", 4))) + 
#   geom_jitter(color="black", size= 1 , alpha=0.9) +
#   stat_summary(fun.data = n_fun, geom = "text") + scale_y_continuous(breaks = seq(0, 1.2, by = 0.1)) + 
#   ylab("percentile agreement") + xlab("type of alterations") + 
#   ggtitle("percentile agreement: 0.8 length 0.8 absCn") + theme(plot.title = element_text(hjust = 0.5))


pa_graph <- ggplot(data = ggplotDf_m_p, aes(x = alterations, y  = concordance/100)) + geom_boxplot(color = c(rep("blue", 4), rep("red", 4))) + 
  geom_jitter(color="black", size= 1 , alpha=0.9) +
  stat_summary(fun.data = n_fun, geom = "text") + ylim(c(0,1.2)) + 
  ylab("percentile agreement") + xlab("type of alterations") + 
  ggtitle("percentile agreement: 0.8 length 0.8 absCn") + theme(plot.title = element_text(hjust = 0.5))


gridExtra::grid.arrange(ck_graph, pa_graph, ncol = 2)
dev.off()
png("/mnt/DATA5/tmp/kev/snpVsPanel/20220627AneuploidyKappa.png", width = 1800, height = 1200)
gridExtra::grid.arrange(ck_graph, pa_graph, ncol = 2)
dev.off()
### can maybe just combine the stats of matched and non-matched
### huge range for mathced euploid ... why? go look at those specific profiles
### maybe the sub clonal cutoff is why?



### do these discrepancies happen at certain chromosomes?i.e just bad coverage?



### for all samples or each sample from SNP data set get grange of change
### and cross reference that with panel data ...
### i.e what perfect of each large scale change is detected ...

panelV3$kappa <- ckappConcord
panelV3$kappaSub <- ckappConcordSub
panelV3$gainSub <- concordListSubGains
panelV3$lossSub <- concordListSubLosses


### purity seems to only affect quality of losses

cor(panelV3$V3_purity[-which(is.na(panelV3$kappaSub))],
    panelV3$kappaSub[-which(is.na(panelV3$kappaSub))])

cor(panelV3$V3_purity[-which(is.na(panelV3$gainSub))],
    panelV3$gainSub[-which(is.na(panelV3$gainSub))])

cor(panelV3$V3_purity[-which(is.na(panelV3$lossSub))],
    panelV3$lossSub[-which(is.na(panelV3$lossSub))])


### ploidy negatively correlated i.e higher ploidy less concordant
cor(panelV3$V3_ploidy_var[-which(is.na(panelV3$kappaSub))],
    panelV3$kappaSub[-which(is.na(panelV3$kappaSub))])

cor(panelV3$V3_ploidy_var[-which(is.na(panelV3$gainSub))],
    panelV3$gainSub[-which(is.na(panelV3$gainSub))])

cor(panelV3$V3_ploidy_var[-which(is.na(panelV3$lossSub))],
    panelV3$lossSub[-which(is.na(panelV3$lossSub))])


library(GenomicRanges)
ascatGainsAndLosses <- NULL
i <- snpRes05$sample[2]
for (i in snpRes05$sample) {
  
  tmpAneuploidy <- ascatAneploidy[which(rownames(ascatAneploidy) == i),]
  tmpAneuploidy <- tmpAneuploidy[1:19]
  tmpSkipChroms <- names(tmpAneuploidy)[which(tmpAneuploidy != "none")]
  tmpDf <- allAscatSegs_filt[grep(i, allAscatSegs_filt$sample),]
  tmpPloidy <- as.numeric(snpRes05$ploidy_int[which(snpRes05$sample == i)])
  tmpDf$totalSc_n <- tmpDf$totalSc- tmpPloidy
  tmpDf <- tmpDf[-which(tmpDf$chr %in% tmpSkipChroms),]
  ascatGainsAndLosses <- rbind(ascatGainsAndLosses, tmpDf)
}

# ascatGainsAndLosses <- ascatGainsAndLosses[which(abs(ascatGainsAndLosses$totalSc_n) > 0.2),]
# 

noshadGainsAndLosses <- NULL
i <- snpRes05$sample[2]
for (i in snpRes05$sample) {
  
  tmpAneuploidy <- ascatAneploidy[which(rownames(noshadAneploidySub) == i),]
  tmpAneuploidy <- tmpAneuploidy[1:19]
  tmpSkipChroms <- names(tmpAneuploidy)[which(tmpAneuploidy != "none")]
  
  tmpDf <- allNoshadSeg_filt[grep(i, allNoshadSeg_filt$sample),]
  tmpPloidy <- as.numeric(panelV3$V3[which(panelV3$sample == i)])
  tmpDf$totalSc_n <- tmpDf$sC- tmpPloidy
  tmpDf <- tmpDf[-which(tmpDf$chr %in% tmpSkipChroms),]
  noshadGainsAndLosses <- rbind(noshadGainsAndLosses, tmpDf)
}

noshadGainsAndLosses <- noshadGainsAndLosses[which(abs(noshadGainsAndLosses$totalSc_n) > 0.7),]

### i wonder if there are just chromosomes that aren't concordant in terms of aneuploidy ... or ones with lower quality 
### need to do it by sample
ascatGrange <- GRanges(seqnames = ascatGainsAndLosses$chr,
                       IRanges(start = ascatGainsAndLosses$startpos, end = ascatGainsAndLosses$endpos))
noshadGrange <- GRanges(seqnames = noshadGainsAndLosses$chr,
                          IRanges(start = as.numeric(noshadGainsAndLosses$start), end = as.numeric(noshadGainsAndLosses$end)))


panelBed <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed",
                       sep = "\t", header = FALSE)
panelBedGrange <- GRanges(seqnames = panelBed$V1,
                          IRanges(start = panelBed$V2, end = panelBed$V3))


### number overlaps, percentage covered, same direction
# numberOv <- NULL
# percentOv <- NULL
# directionOv <- NULL
# i <- snpRes05$sample[1]
# for (i in snpRes05$sample) {
#   print(i)
#   tmpAscatGrange <- ascatGrange[grep(i, ascatGainsAndLosses$sample)]
#   tmpNoshadGrange <- noshadGrange[grep(i, noshadGainsAndLosses$sample)]
#   tmpOv <- findOverlaps(tmpAscatGrange, tmpNoshadGrange)
#   # tmpOv <- findOverlaps(tmpNoshadGrange, tmpAscatGrange)
#   
#   cnAscat <- sign(ascatGainsAndLosses$totalSc_n[grep(i, ascatGainsAndLosses$sample)])
#   cnNoshad <- sign(noshadGainsAndLosses$totalSc_n[grep(i, noshadGainsAndLosses$sample)])
#   
#   
#   tmpNumber <- rep(0, length(tmpAscatGrange))
#   tmpPercent <- rep(0, length(tmpAscatGrange))
#   tmpDirection <- rep(0, length(tmpAscatGrange))
#   
#   j <- unique(queryHits(tmpOv))[1]
#   for (j in unique(queryHits(tmpOv))) {
#     print(j)
#     tmpAscatGrange2 <- tmpAscatGrange[j]
#     noshadIdx <- subjectHits(tmpOv)[which(queryHits(tmpOv) == j)]
#     tmpNoshadGrange2 <- tmpNoshadGrange[noshadIdx]
#     
#     tmpNumber[j] <- length(tmpNoshadGrange2)
#     tmpInter <- GenomicRanges::intersect(tmpAscatGrange2, tmpNoshadGrange2)
#     tmpPercent[j] <- tmpInter@ranges@width/tmpAscatGrange2@ranges@width
#     
#     tmpDirection2 <- 0
#     for (k in 1:length(tmpNoshadGrange2)) {
#       tmpInter2 <- GenomicRanges::intersect(tmpAscatGrange2, tmpNoshadGrange2[k])
#       tmpDirection2 <- tmpDirection2 + tmpInter2@ranges@width/tmpAscatGrange2@ranges@width * cnAscat[j] * cnNoshad[noshadIdx[k]]
#     }
#     tmpDirection[j] <- tmpDirection2
#     print(warnings())
#   }
#   numberOv <- c(numberOv, tmpNumber)
#   percentOv <- c(percentOv, tmpPercent)
#   directionOv <- c(directionOv, tmpDirection)
# }

### be sure to get graph marker coverage by number of markers in panel



## possibly different b/c of where it is parallelized
### try to parallelize at the sample level i.e at i

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


registerDoParallel(20)

res <- foreach(i = snpRes05$sample, .combine='comb', .multicombine=TRUE, .packages = c('GenomicRanges')) %dopar% {
  tmpAscatGrange <- ascatGrange[grep(i, ascatGainsAndLosses$sample)]
  tmpNoshadGrange <- noshadGrange[grep(i, noshadGainsAndLosses$sample)]
  tmpOv <- findOverlaps(tmpAscatGrange, tmpNoshadGrange)
  cnAscat <- sign(ascatGainsAndLosses$totalSc_n[grep(i, ascatGainsAndLosses$sample)])
  cnNoshad <- sign(noshadGainsAndLosses$totalSc_n[grep(i, noshadGainsAndLosses$sample)])
  
  tmpNumber <- rep(0, length(tmpAscatGrange))
  tmpPercent <- rep(0, length(tmpAscatGrange))
  tmpDirection <- rep(0, length(tmpAscatGrange))
  
  for(j in unique(queryHits(tmpOv))){
    
    tmpAscatGrange2 <- tmpAscatGrange[j]
    noshadIdx <- subjectHits(tmpOv)[which(queryHits(tmpOv) == j)]
    tmpNoshadGrange2 <- tmpNoshadGrange[noshadIdx]
    
    tmpNumber[j] <- length(tmpNoshadGrange2)
    tmpInter <- GenomicRanges::intersect(tmpAscatGrange2, tmpNoshadGrange2)
    tmpPercent[j] <- tmpInter@ranges@width/tmpAscatGrange2@ranges@width
    
    tmpDirection2 <- 0
    for (k in 1:length(tmpNoshadGrange2)) {
      tmpInter2 <- GenomicRanges::intersect(tmpAscatGrange2, tmpNoshadGrange2[k])
      tmpDirection2 <- tmpDirection2 + tmpInter2@ranges@width/tmpAscatGrange2@ranges@width * cnAscat[j] * cnNoshad[noshadIdx[k]]
    }
    tmpDirection[j] <- tmpDirection2
  }
  
  list(tmpNumber, tmpPercent, tmpDirection)
  
}


numberOvPar <- do.call(c, res[[1]])
percentOvPar<- do.call(c, res[[2]])
directionOvPar <- do.call(c, res[[3]])



ascatGainsAndLosses$numberOv <- numberOvPar
ascatGainsAndLosses$percentOv <- percentOvPar
ascatGainsAndLosses$directionOv <- directionOvPar
ascatGainsAndLosses$sizeMb <- ascatGainsAndLosses$size/1e6
ascatGainsAndLosses$sample <- factor(ascatGainsAndLosses$sample, levels = unique(ascatGainsAndLosses$sample))

numMarkers <- NULL
for (i in 1:length(ascatGrange)) {
  res <- length(findOverlaps(ascatGrange[i], panelBedGrange))
  numMarkers <- c(numMarkers, res)
}

ascatGainsAndLosses$numMarks <- numMarkers


ggplot(ascatGainsAndLosses) + geom_boxplot(aes(chr, percentOv))
qual <- ggplot(ascatGainsAndLosses) + geom_boxplot(aes(sample, percentOv)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
size <- ggplot(ascatGainsAndLosses) + geom_boxplot(aes(sample, sizeMb)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
# pdf(file = "/mnt/DATA5/tmp/kev/snpVsPanel/20220516overlapBySample.pdf", width = 12, height = 6)
gridExtra::grid.arrange(qual, size, ncol = 2)
dev.off()

mb_graph <- ggplot(ascatGainsAndLosses) + geom_point(aes(sizeMb, percentOv))
nm_graph <- ggplot(ascatGainsAndLosses) + geom_point(aes(numMarkers, percentOv))

dev.off()
# pdf(file = "/mnt/DATA5/tmp/kev/snpVsPanel/20220516overlapSizeMarkers.pdf", width = 12, height = 6)
gridExtra::grid.arrange(mb_graph, nm_graph, ncol = 2)
dev.off()

ggplot(ascatGainsAndLosses) + geom_point(aes(sizeMb/numMarkers, percentOv))

ascatGainsAndLosses_filt <- ascatGainsAndLosses[which(ascatGainsAndLosses$percentOv == 1),]
cor(x = ascatGainsAndLosses_filt$percentOv, y = ascatGainsAndLosses_filt$numMarks)

# look to see if the 0 overlap areas are subclonal events oareas with low cnr signals
### doing reverse for ascat vs noshad

registerDoParallel(20)

resN <- foreach(i = snpRes05$sample, .combine='comb', .multicombine=TRUE, .packages = c('GenomicRanges')) %dopar% {
  tmpAscatGrange <- ascatGrange[grep(i, ascatGainsAndLosses$sample)]
  tmpNoshadGrange <- noshadGrange[grep(i, noshadGainsAndLosses$sample)]
  
  tmpOv <- findOverlaps(tmpNoshadGrange, tmpAscatGrange)
  
  cnAscat <- sign(ascatGainsAndLosses$totalSc_n[grep(i, ascatGainsAndLosses$sample)])
  cnNoshad <- sign(noshadGainsAndLosses$totalSc_n[grep(i, noshadGainsAndLosses$sample)])
  
  tmpNumber <- rep(0, length(tmpNoshadGrange))
  tmpPercent <- rep(0, length(tmpNoshadGrange))
  tmpDirection <- rep(0, length(tmpNoshadGrange))
  
  for(j in unique(queryHits(tmpOv))){
    tmpNoshadGrange2 <- tmpNoshadGrange[j]
    ascatIdx <- subjectHits(tmpOv)[which(queryHits(tmpOv) == j)]
    tmpAscatGrange2 <- tmpAscatGrange[ascatIdx]
    
    tmpNumber[j] <- length(tmpAscatGrange2)
    tmpInter <- GenomicRanges::intersect(tmpNoshadGrange2, tmpAscatGrange2)
    tmpPercent[j] <- tmpInter@ranges@width/tmpNoshadGrange2@ranges@width
    
    tmpDirection2 <- 0
    for (k in 1:length(tmpAscatGrange2)) {
      tmpInter2 <- GenomicRanges::intersect(tmpNoshadGrange2, tmpAscatGrange2[k])
      tmpDirection2 <- tmpDirection2 + tmpInter2@ranges@width/tmpNoshadGrange2@ranges@width * cnNoshad[j] * cnAscat[ascatIdx[k]]
    }
    tmpDirection[j] <- tmpDirection2
  }
  
  list(tmpNumber, tmpPercent, tmpDirection)
  
}



numberOvParN <- do.call(c, resN[[1]])
percentOvParN <- do.call(c, resN[[2]])
directionOvParN <- do.call(c, resN[[3]])



noshadGainsAndLosses$numberOv <- numberOvParN
noshadGainsAndLosses$percentOv <- percentOvParN
noshadGainsAndLosses$directionOv <- directionOvParN
noshadGainsAndLosses$sizeMb <- noshadGainsAndLosses$size/1e6
noshadGainsAndLosses$sample <- factor(noshadGainsAndLosses$sample, levels = unique(noshadGainsAndLosses$sample))

numMarkersN <- NULL
for (i in 1:length(noshadGrange)) {
  res <- length(findOverlaps(noshadGrange[i], panelBedGrange))
  numMarkersN <- c(numMarkersN, res)
}

noshadGainsAndLosses$numMarks <- numMarkersN


ggplot(noshadGainsAndLosses) + geom_boxplot(aes(chr, percentOv))
qualN <- ggplot(noshadGainsAndLosses) + geom_boxplot(aes(sample, percentOv)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
sizeN <- ggplot(noshadGainsAndLosses) + geom_boxplot(aes(sample, sizeMb)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
# pdf(file = "/mnt/DATA5/tmp/kev/snpVsPanel/20220516overlapBySample.pdf", width = 12, height = 6)
gridExtra::grid.arrange(qualN, sizeN, ncol = 2)
dev.off()


dev.off()
pdf(file = "/mnt/DATA5/tmp/kev/snpVsPanel/20220627overlapBySample.pdf", width = 12, height = 6)
gridExtra::grid.arrange(qualN, sizeN, ncol = 2)
dev.off()

dev.off()
png(file = "/mnt/DATA5/tmp/kev/snpVsPanel/20220627overlapBySample.png", width = 1200, height = 900)
gridExtra::grid.arrange(qualN, sizeN, ncol = 2)
dev.off()



mb_graphN <- ggplot(noshadGainsAndLosses) + geom_point(aes(sizeMb, percentOv))
nm_graphN <- ggplot(noshadGainsAndLosses) + geom_point(aes(numMarkersN, percentOv))

dev.off()
# pdf(file = "/mnt/DATA5/tmp/kev/snpVsPanel/20220516overlapSizeMarkers.pdf", width = 12, height = 6)
gridExtra::grid.arrange(mb_graphN, nm_graphN, ncol = 2)
dev.off()

ggplot(noshadGainsAndLosses) + geom_point(aes(sizeMb/numMarkersN, percentOv))

ascatGainsAndLosses_filt <- ascatGainsAndLosses[which(ascatGainsAndLosses$percentOv == 1),]
cor(x = ascatGainsAndLosses_filt$percentOv, y = ascatGainsAndLosses_filt$numMarks)


length(which(snpRes05$ploidy_int == allPloidyCallsV2$ploidy_int & snpRes05$ploidy_int == 2))/length(which(snpRes05$ploidy_int == 2))
length(which(snpRes05$ploidy_int == allPloidyCallsV2$ploidy_int & snpRes05$ploidy_int == 3))/length(which(snpRes05$ploidy_int == 3))
length(which(snpRes05$ploidy_int == allPloidyCallsV2$ploidy_int & snpRes05$ploidy_int == 4))/length(which(snpRes05$ploidy_int == 4))


### previous version /home/kevhu/scripts/20220512comparingPanelsV2.R

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


allAscatSegs <- read.table("/mnt/DATA5/tmp/kev/misc/20220512AscatIllusRes05Seg100.txt", sep = "\t",
                           stringsAsFactors = FALSE, header = TRUE)

allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "14399RT", "14399rt")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "2027_LT-E", "2027lte")
allAscatSegs$sample <- str_replace_all(allAscatSegs$sample, "13085ROT\\(L\\)", "13085lt")

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
tmpChrList <- paste0("chr", c(1:19))

i <- rownames(ascatAneploidy)[4]
### >0.20 cutoff from 2018 paper Allison M Taylor for arm-level stuff more or less doubled b/c mouse have no arms like chr14 in humans
for (i in rownames(ascatAneploidy)) {
  print(i)
  tmpDf <- allAscatSegs_filt[grep(i, allAscatSegs_filt$sample),]
  tmpPloidy <- snpRes05$ploidy_int[which(snpRes05$sample == i)]
  tmpDf$totalSc_n <- tmpDf$totalSc - tmpPloidy
  j <- "chr3"
  for (j in tmpChrList) {
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


finalRes <- NULL
variVec <- seq(0.05,1, 0.05)
for (z in variVec) {
  for (i in rownames(noshadAneploidy)) {
    tmpDf <- allNoshadSeg_filt[grep(i, allNoshadSeg_filt$sample),]
    tmpPloidy <- as.numeric(panelV3$V3[which(panelV3$sample == i)])
    tmpDf$totalSc_n <- tmpDf$sC- tmpPloidy
    
    for (j in tmpChrList) {
      tmpChr <- tmpDf[which(tmpDf$chr == j), ]
      chrThres <- sum(tmpChr$size) * z
      tmpAneu <- "none"
      if (sum(tmpChr$size[which(tmpChr$totalSc_n > 0.8)]) > chrThres) {
        tmpAneu <- "gain"
      } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < -0.8)]) > chrThres){
        tmpAneu <- "loss"
        # } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < 0.8 & tmpChr$totalSc_n > 0.2)]) > chrThres){
      } else if(sum(tmpChr$size[which(tmpChr$totalSc_n > 0.5)]) > chrThres){
        tmpAneu <- "sub"
        # } else if(sum(tmpChr$size[which(tmpChr$totalSc_n > -0.8 & tmpChr$totalSc_n < -0.2)]) > chrThres){
      } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < -0.5)]) > chrThres){
        tmpAneu <- "sub"
      } else{
        tmpAneu <- "none"
      }
      noshadAneploidy[i,j] <- tmpAneu
    }
  }
  
  tmpChrKapp <- NULL
  tmpChrCcc <- NULL
  for (y in tmpChrList) {
    tmpAscat <- ascatAneploidy[, y]
    tmpNoshad <- noshadAneploidy[, y]
    
    if (length(which(tmpNoshad == "sub")) > 0) {
      tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
    }
    
    tmpKapp <- kappa2(cbind(tmpAscat, tmpNoshad))
    tmpChrKapp <- c(tmpChrKapp, tmpKapp$value)
    
  }
  finalRes <- cbind(finalRes, tmpChrKapp)
}

colnames(finalRes) <- variVec
rownames(finalRes) <- tmpChrList

vectorOfCutoffs <- NULL
vectorOfMax <- apply(finalRes, 1, max)
for (i in seq_along(vectorOfMax)) {
  tmpCutoff <- colnames(finalRes)[which(finalRes[i,] == vectorOfMax[i])]
  vectorOfCutoffs <- c(vectorOfCutoffs, tmpCutoff[1])
}

vectorOfCutoffs <- as.numeric(vectorOfCutoffs)
vectorOfCutoffs[which(vectorOfCutoffs < 0.5)] <- 0.5
### use cutoffs from above per chromosome


noshadAneploidy <- matrix(0, nrow = length(unique(panelV3$sample)), ncol = length(unique(allNoshadSeg_filt$chr)))
rownames(noshadAneploidy) <- unique(panelV3$sample)
colnames(noshadAneploidy) <- unique(allNoshadSeg_filt$chr)
allNoshadSeg_filt$size <- as.numeric(allNoshadSeg_filt$end.off) - as.numeric(allNoshadSeg_filt$start.off)


### use vector of cutoffs
for (i in rownames(noshadAneploidy)) {
  tmpDf <- allNoshadSeg_filt[grep(i, allNoshadSeg_filt$sample),]
  tmpPloidy <- as.numeric(panelV3$V3[which(panelV3$sample == i)])
  tmpDf$totalSc_n <- tmpDf$sC- tmpPloidy
  
  for (j in seq_along(tmpChrList)) {
    tmpChr <- tmpDf[which(tmpDf$chr == tmpChrList[j]), ]
    chrThres <- sum(tmpChr$size) * vectorOfCutoffs[j]
    tmpAneu <- "none"
    if (sum(tmpChr$size[which(tmpChr$totalSc_n > 0.8)]) > chrThres) {
      tmpAneu <- "gain"
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < -0.8)]) > chrThres){
      tmpAneu <- "loss"
      # } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < 0.8 & tmpChr$totalSc_n > 0.1)]) > chrThres){
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n > 0.5)]) > chrThres){
      tmpAneu <- "sub"
      # } else if(sum(tmpChr$size[which(tmpChr$totalSc_n > -0.8 & tmpChr$totalSc_n < -0.1)]) > chrThres){
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < -0.5)]) > chrThres){
      tmpAneu <- "sub"
    } else{
      tmpAneu <- "none"
    }
    noshadAneploidy[i,j] <- tmpAneu
  }
}




### same thing but with subclones



noshadAneploidy_sub <- NULL

ckappConcordSub <- NULL
concordListSub <- NULL
cccList <- NULL
i <- rownames(ascatAneploidy)[2]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- noshadAneploidy[i, 1:19]
  # tmpAscat <- ascatAneploidy[i, ]
  # tmpNoshad <- noshadAneploidy[i, ]
  # 
  if (length(which(tmpNoshad == "sub")) > 0) {
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  concord <- agree(cbind(cat_to_cont(tmpAscat), cat_to_cont(tmpNoshad)))
  concordListSub <- c(concordListSub, concord$value)
  
  tmpKapp <- kappa2(cbind(tmpAscat, tmpNoshad))
  ckappConcordSub <- c(ckappConcordSub, tmpKapp$value)
  
  noshadAneploidy_sub <- rbind(noshadAneploidy_sub, tmpNoshad)
}

rownames(noshadAneploidy_sub) <- rownames(noshadAneploidy)

### trying something new - finding optimal cutoff per chromosomes b/c there are variable
### coverage per chromosome

subGainLength <- NULL

concordListSubGains <- NULL
i <- rownames(ascatAneploidy)[16]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- noshadAneploidy[i, 1:19]
  
  if (length(which(tmpNoshad == "sub")) > 0) {
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  
  tmpIdx <- c(which(tmpAscat == "gain"), which(tmpNoshad == "gain"))
  subGainLength <- c(subGainLength, length(tmpIdx))
  
  if (length(tmpIdx) < 1) {
    concordListSubGains <- c(concordListSubGains, NA)
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListSubGains  <- c(concordListSubGains , concord$value)
    
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListSubGains  <- c(concordListSubGains , concord$value)
  }
}



subLossLength <- NULL

concordListSubLosses <- NULL
i <- rownames(ascatAneploidy)[41]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy[i, 1:19]
  tmpNoshad <- noshadAneploidy[i, 1:19]
  
  if (length(which(tmpNoshad == "sub")) > 0) {
    tmpNoshad[which(tmpNoshad == "sub")] <- tmpAscat[which(tmpNoshad == "sub")]
  }
  
  tmpIdx <- c(which(tmpAscat == "loss"), which(tmpNoshad == "loss"))
  subLossLength <- c(subLossLength, length(tmpIdx))
  
  if (length(tmpIdx) < 1) {
    concord <- NA
    concordListSubLosses <- c(concordListSubLosses, NA)  
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListSubLosses <- c(concordListSubLosses, concord$value)
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    concord <- agree(cbind(cat_to_cont(tmpAscat[tmpIdx]), cat_to_cont(tmpNoshad[tmpIdx])))
    concordListSubLosses <- c(concordListSubLosses, concord$value)
  }
}


### looking at gains and losses detection by chromosome



### make a 4 column data frame (1) sample (2) kappa overall (3) kappa loss (4) kappa gain

allKappDf <- data.frame("sample" = rownames(ascatAneploidy), "subLosses" = concordListSubLosses, "lossCount" = subLossLength,
                        "subGains" = concordListSubGains, "gainCount" = subGainLength)
allKappDf$tc <- panelV3$purityV2[match(panelV3$sample, allKappDf$sample)]
allKappDf$panelPloidy <- panelV3$V3[match(panelV3$sample, allKappDf$sample)]
allKappDf$snpPloidy <- snpRes05$ploidy_int[match(snpRes05$sample, allKappDf$sample)]
allKappDf$match <- ifelse(allKappDf$panelPloidy == allKappDf$snpPloidy, "yes", "no")

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

ggplotDf_ck <- ggplotDf_ck[-which(is.na(ggplotDf_ck$ckap)), ]

ck_graph <- ggplot(data = ggplotDf_ck, aes(x = alterations, y  = ckap)) + geom_boxplot() + 
  geom_jitter(color="black", size= 1 , alpha=0.9) +
  stat_summary(fun.data = n_fun, geom = "text") + ylim(c(0,1.2)) + 
  ggtitle("Aneuploidy cohen's k: NGS(Var) & SNP(0.7)") + theme(plot.title = element_text(hjust = 0.5)) + 
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
ggplotDf_m_p <- ggplotDf_m_p[-which(is.na(ggplotDf_m_p$concordance)), ]

pa_graph <- ggplot(data = ggplotDf_m_p, aes(x = alterations, y  = concordance/100)) + geom_boxplot(color = c(rep("blue", 4), rep("red", 4))) + 
  geom_jitter(color="black", size= 1 , alpha=0.9) +
  stat_summary(fun.data = n_fun, geom = "text") + ylim(c(0,1.2)) + 
  ylab("percentile agreement") + xlab("type of alterations") + 
  ggtitle("Gains/loses % agreement: NGS(Var) & SNP(0.7)") + theme(plot.title = element_text(hjust = 0.5))


gridExtra::grid.arrange(ck_graph, pa_graph, ncol = 2)


dev.off()
# png("/mnt/DATA5/tmp/kev/snpVsPanel/20220714AneuploidyKappa.png", width = 1500, height = 900)
gridExtra::grid.arrange(ck_graph, pa_graph, ncol = 2)
dev.off()




### making it by all matched and all ploidies
###
###

cListSubLosses_m_dip <-  concordListSubLosses[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 2)]
cListSubLosses_m_trip <-  concordListSubLosses[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 3)]
cListSubLosses_m_tetra <-  concordListSubLosses[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 4)]
cListSubLosses_um_dip <-  concordListSubLosses[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 2)]
cListSubLosses_um_trip <-  concordListSubLosses[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 3)]
cListSubLosses_um_tetra <-  concordListSubLosses[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 4)]

cListSubGains_m_dip <-  concordListSubGains[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 2)]
cListSubGains_m_trip <-  concordListSubGains[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 3)]
cListSubGains_m_tetra <-  concordListSubGains[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 4)]
cListSubGains_um_dip <-  concordListSubGains[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 2)]
cListSubGains_um_trip <-  concordListSubGains[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 3)]
cListSubGains_um_tetra <-  concordListSubGains[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 4)]


alterations_all <- c(rep("loss_m_dip", length(cListSubLosses_m_dip)), rep("loss_m_trip", length(cListSubLosses_m_trip)),
                     rep("loss_m_tetra", length(cListSubLosses_m_tetra)),
                     rep("loss_um_dip", length(cListSubLosses_um_dip)),
                     rep("loss_um_trip", length(cListSubLosses_um_trip)), rep("loss_um_tetra", length(cListSubLosses_um_tetra)),
                     
                     rep("gain_m_dip", length(cListSubGains_m_dip)), rep("gain_m_trip", length(cListSubGains_m_trip)),
                     rep("gain_m_tetra", length(cListSubGains_m_tetra)),
                     rep("gain_um_dip", length(cListSubGains_um_dip)), rep("gain_um_trip", length(cListSubGains_um_trip)),
                     rep("gain_um_tetra", length(cListSubGains_um_tetra)))

concordance_all <- c(cListSubLosses_m_dip, cListSubLosses_m_trip, cListSubLosses_m_tetra,
                     cListSubLosses_um_dip, cListSubLosses_um_trip, cListSubLosses_um_tetra,
                     cListSubGains_m_dip, cListSubGains_m_trip, cListSubGains_m_tetra,
                     cListSubGains_um_dip, cListSubGains_um_trip, cListSubGains_um_tetra)


ggplotDf_m_p_all <- data.frame("alterations" = alterations_all, "concordance" = concordance_all)
ggplotDf_m_p_all$alterations <- factor(ggplotDf_m_p_all$alterations, levels = unique(ggplotDf_m_p_all$alterations))

pa_graph_all <- ggplot(data = ggplotDf_m_p_all, aes(x = alterations, y  = concordance/100)) + geom_boxplot(color = c(rep("blue", 6), rep("red", 6))) + 
  geom_jitter(color="black", size= 1 , alpha=0.9) +
  stat_summary(fun.data = n_fun, geom = "text") + ylim(c(0,1.2)) + 
  ylab("percentile agreement") + xlab("type of alterations") + 
  ggtitle("Gains/loses % agreement: NGS(Var) & SNP(0.7)") + theme(plot.title = element_text(hjust = 0.5))




ckappConcordSub_m_dip <-  ckappConcordSub[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 2)]
ckappConcordSub_m_trip <-  ckappConcordSub[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 3)]
ckappConcordSub_m_tetra <-  ckappConcordSub[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 4)]
ckappConcordSub_um_dip <-  ckappConcordSub[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 2)]
ckappConcordSub_um_trip <-  ckappConcordSub[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 3)]
ckappConcordSub_um_tetra <-  ckappConcordSub[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 4)]


alterations_all2 <- c(rep("ckapp_m_dip", length(ckappConcordSub_m_dip)), rep("ckapp_m_trip", length(ckappConcordSub_m_trip)),
                      rep("ckapp_m_tetra", length(ckappConcordSub_m_tetra)),
                      rep("ckapp_um_dip", length(ckappConcordSub_um_dip)), rep("ckapp_um_trip", length(ckappConcordSub_um_trip)),
                      rep("ckapp_um_tetra", length(ckappConcordSub_um_tetra)))

concordance_all2 <- c(ckappConcordSub_m_dip, ckappConcordSub_m_trip, ckappConcordSub_m_tetra,
                      ckappConcordSub_um_dip, ckappConcordSub_um_trip, ckappConcordSub_um_tetra)

ggplotDf_ck_all <- data.frame("alterations" = alterations_all2, "ckap" = concordance_all2)
ggplotDf_ck_all$alterations <- factor(ggplotDf_ck_all$alterations, levels = unique(ggplotDf_ck_all$alterations))




ck_graph_all <- ggplot(data = ggplotDf_ck_all, aes(x = alterations, y  = ckap)) + geom_boxplot() + 
  geom_jitter(color="black", size= 1 , alpha=0.9) +
  stat_summary(fun.data = n_fun, geom = "text") + ylim(c(0,1.2)) + 
  ggtitle("Aneuploidy cohen's k: NGS(Var) & SNP(0.7)") + theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("match + ploidy")

gridExtra::grid.arrange(ck_graph_all, pa_graph_all, ncol = 1)

dev.off()
# png("/mnt/DATA5/tmp/kev/snpVsPanel/20220718AneuploidyByPloidy.png", width = 1500, height = 900)
gridExtra::grid.arrange(ck_graph_all, pa_graph_all, ncol = 1)
dev.off()


### graph by sample too ... and just color the changes ....


### quick logic check match
# which(panelV3$V3 == 2 & panelV3$V3 == snpRes05$ploidy_int)
# which(panelV3$V3 != 2 & panelV3$V3 == snpRes05$ploidy_int)
# which(panelV3$V3 == 2 & panelV3$V3 != snpRes05$ploidy_int)
# which(panelV3$V3 != 2 & panelV3$V3 != snpRes05$ploidy_int)
# 
# which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 2)
# which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 != 2)
# which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 2)
# which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 != 2)


# ck_graph_violin <- ggplot(data = ggplotDf_ck, aes(x = alterations, y  = ckap)) + geom_violin() + geom_boxplot(width = 0.05) +
#   geom_jitter(color="black", size= 1 , alpha=0.9) +
#   stat_summary(fun.data = n_fun, geom = "text") + ylim(c(0,1.2)) + 
#   ggtitle("Aneuploidy cohen's k: NGS(Var) & SNP(0.7)") + theme(plot.title = element_text(hjust = 0.5)) + 
#   xlab("match + ploidy")
# 
# pa_graph_violin <- ggplot(data = ggplotDf_m_p, aes(x = alterations, y  = concordance/100)) + geom_violin() + geom_boxplot(width = 0.05) +
#   geom_jitter(color="black", size= 1 , alpha=0.9) +
#   stat_summary(fun.data = n_fun, geom = "text") + ylim(c(0,1.2)) + 
#   ylab("percentile agreement") + xlab("type of alterations") + 
#   ggtitle("Gains/loses % agreement: NGS(Var) & SNP(0.7)") + theme(plot.title = element_text(hjust = 0.5))
# 
# gridExtra::grid.arrange(ck_graph_violin, pa_graph_violin, ncol = 2)

### calculating large genomic changes

ascatGainsAndLosses <- NULL
for (i in snpRes05$sample) {
  
  # tmpAneuploidy <- ascatAneploidy[which(rownames(ascatAneploidy) == i),]
  # 
  tmpAneuploidy <- noshadAneploidy[which(rownames(noshadAneploidy_sub) == i),]
  tmpAneuploidy <- tmpAneuploidy[1:19]
  tmpSkipChroms <- names(tmpAneuploidy)[which(tmpAneuploidy != "none")]
  
  tmpDf <- allAscatSegs_filt[grep(i, allAscatSegs_filt$sample),]
  tmpPloidy <- as.numeric(snpRes05$ploidy_int[which(snpRes05$sample == i)])
  tmpDf$totalSc_n <- tmpDf$totalSc- tmpPloidy
  
  if (length(tmpSkipChroms) > 0) {
    tmpDf <- tmpDf[-which(tmpDf$chr %in% tmpSkipChroms),]
  }
  
  ascatGainsAndLosses <- rbind(ascatGainsAndLosses, tmpDf)
}

ascatGainsAndLosses <- ascatGainsAndLosses[which(abs(ascatGainsAndLosses$totalSc_n) > 0), ]



noshadGainsAndLosses <- NULL
i <- snpRes05$sample[4]
for (i in snpRes05$sample) {
  
  print(i)
  tmpAneuploidy <- noshadAneploidy[which(rownames(noshadAneploidy_sub) == i),]
  tmpAneuploidy <- tmpAneuploidy[1:19]
  tmpSkipChroms <- names(tmpAneuploidy)[which(tmpAneuploidy != "none")]
  
  tmpDf <- allNoshadSeg_filt[grep(i, allNoshadSeg_filt$sample),]
  tmpPloidy <- as.numeric(panelV3$V3[which(panelV3$sample == i)])
  tmpDf$totalSc_n <- tmpDf$sC- tmpPloidy
  tmpDf$ploidy <- tmpPloidy
  if (length(tmpSkipChroms) > 0) {
    tmpDf <- tmpDf[-which(tmpDf$chr %in% tmpSkipChroms),]
  }
  
  noshadGainsAndLosses <- rbind(noshadGainsAndLosses, tmpDf)
}

noshadGainsAndLosses <- noshadGainsAndLosses[which(abs(noshadGainsAndLosses$totalSc_n) > 0.5),]

### filtering to zero for both noshad and ascat only used in creating min common region
# noshadGainsAndLosses$totalSc_n[which(abs(noshadGainsAndLosses$totalSc_n) < 0.1)] <- 0


ascatGrange <- GRanges(seqnames = ascatGainsAndLosses$chr,
                       IRanges(start = ascatGainsAndLosses$startpos, end = ascatGainsAndLosses$endpos))
noshadGrange <- GRanges(seqnames = noshadGainsAndLosses$chr,
                        IRanges(start = as.numeric(noshadGainsAndLosses$start), end = as.numeric(noshadGainsAndLosses$end)))


panelBed <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed",
                       sep = "\t", header = FALSE)
load("/mnt/DATA5/tmp/kev/misc/20220628badIdx.RData")
panelBed <- panelBed[-badAmps,]
panelBedGrange <- GRanges(seqnames = panelBed$V1,
                          IRanges(start = panelBed$V2, end = panelBed$V3))

panelReduced <- NULL
for (i in unique(panelBed$V8)) {
  tmpBed <- panelBed[which(panelBed$V8 == i), ]
  panelReduced <- rbind(panelReduced, c(tmpBed$V1[1], mean(tmpBed$V2), mean(tmpBed$V3), tmpBed$V8[1]))
}

panelReduced  <- data.frame(panelReduced)
panelReduced[2:3] <- lapply(panelReduced[2:3], function(x) round(as.numeric(x)))

panelBedGrangeRed <- GRanges(seqnames = panelReduced$X1,
                             IRanges(start = panelReduced$X2, end = panelReduced$X3))

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

registerDoParallel(20)

i <- snpRes05$sample[2]
resN <- foreach(i = snpRes05$sample, .combine='comb', .multicombine=TRUE, .packages = c('GenomicRanges')) %dopar% {
  tmpAscatGrange <- ascatGrange[grep(i, ascatGainsAndLosses$sample)]
  tmpNoshadGrange <- noshadGrange[grep(i, noshadGainsAndLosses$sample)]
  
  tmpOv <- findOverlaps(tmpNoshadGrange, tmpAscatGrange)
  
  cnAscat <- sign(ascatGainsAndLosses$totalSc_n[grep(i, ascatGainsAndLosses$sample)])
  cnNoshad <- sign(noshadGainsAndLosses$totalSc_n[grep(i, noshadGainsAndLosses$sample)])
  
  tmpNumber <- rep(0, length(tmpNoshadGrange))
  tmpPercent <- rep(0, length(tmpNoshadGrange))
  tmpDirection <- rep(0, length(tmpNoshadGrange))
  j <- 8
  for(j in unique(queryHits(tmpOv))){
    tmpNoshadGrange2 <- tmpNoshadGrange[j]
    ascatIdx <- subjectHits(tmpOv)[which(queryHits(tmpOv) == j)]
    tmpAscatGrange2 <- tmpAscatGrange[ascatIdx]
    
    tmpNumber[j] <- length(tmpAscatGrange2)
    tmpInter <- GenomicRanges::intersect(tmpNoshadGrange2, tmpAscatGrange2)
    tmpPercent[j] <- sum(tmpInter@ranges@width)/tmpNoshadGrange2@ranges@width
    
    print(sum(tmpInter@ranges@width)/tmpNoshadGrange2@ranges@width)
    
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

# noshadGainsAndLosses$numMarks <- 0
# for (i in 1:length(noshadGrange)) {
#   noshadGainsAndLosses$numMarks[i] <- length(findOverlaps(noshadGrange[i], panelBedGrange))
# }

# noshadGainsAndLosses$numMarksRed <- 0
# for (i in 1:length(noshadGrange)) {
#   noshadGainsAndLosses$numMarksRed[i] <- length(findOverlaps(noshadGrange[i], panelBedGrangeRed))
# }


noshadGainsAndLosses$type <- ifelse(sign(noshadGainsAndLosses$totalSc_n) == 1, "darkred", "darkblue")

# noshadGainsAndLosses_filt <- noshadGainsAndLosses[which(noshadGainsAndLosses$sizeMb < 1), ]
# noshadGainsAndLosses <- noshadGainsAndLosses[which(noshadGainsAndLosses$sizeMb > 5),]

noshadGainsAndLosses5Mb <- noshadGainsAndLosses[which(noshadGainsAndLosses$sizeMb > 5),]

qualN <- ggplot(noshadGainsAndLosses5Mb, aes(sample, percentOv)) + geom_boxplot() + 
  geom_jitter(col = noshadGainsAndLosses5Mb$type, size= 1 , alpha=0.8) +
  ggtitle("Large CNAs (> 5Mb)") + ylab("Percent overlap") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))


# pdf(file = "/mnt/DATA5/tmp/kev/snpVsPanel/20220516overlapBySample.pdf", width = 12, height = 6)
qualN


# noshadGainsAndLosses_filt <- noshadGainsAndLosses[-which(noshadGainsAndLosses$sizeMb < 50), ]
# noshadGainsAndLosses_filt <- noshadGainsAndLosses_filt[-which(noshadGainsAndLosses$sizeMb < 25), ]

# qualN2 <- ggplot(noshadGainsAndLosses_filt, aes(sample, percentOv)) + geom_boxplot() +
#   geom_jitter(col = noshadGainsAndLosses_filt$type, size= 1 , alpha=0.8) +
#   ggtitle("Large CNAs (< 1Mb)") + ylab("Percent overlap") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))
# 
# gridExtra::grid.arrange(qualN, qualN2)

dev.off()
# png("/mnt/DATA5/tmp/kev/snpVsPanel/20220703LargeCnas.png", width = 2000, height = 1000)
qualN
dev.off()

# table(noshadGainsAndLosses_filt$sample)/table(noshadGainsAndLosses$sample)

length(which(panelV3$V3 == 2 & panelV3$V3 == snpRes05$ploidy_int))/length(which(panelV3$V3 == 2))
length(which(panelV3$V3 == 3 & panelV3$V3 == snpRes05$ploidy_int))/length(which(panelV3$V3 == 3))
length(which(panelV3$V3 == 4 & panelV3$V3 == snpRes05$ploidy_int))/length(which(panelV3$V3 == 4))

length(which(snpRes05$ploidy_int == 2 & snpRes05$ploidy_int == panelV3$V3))/length(which(snpRes05$ploidy_int == 2))
length(which(snpRes05$ploidy_int == 3 & snpRes05$ploidy_int == panelV3$V3))/length(which(snpRes05$ploidy_int == 3))
length(which(snpRes05$ploidy_int == 4 & snpRes05$ploidy_int == panelV3$V3))/length(which(snpRes05$ploidy_int == 4))




# ggplot(noshadGainsAndLosses, aes(sizeMb, percentOv)) + geom_point(col = noshadGainsAndLosses$type) + facet_wrap(vars(chr, type))
# ggplot(noshadGainsAndLosses, aes(numMarks, percentOv)) + geom_point(col = noshadGainsAndLosses$type) + facet_wrap(vars(chr, type))
# ggplot(noshadGainsAndLosses, aes(sizeMb/numMarks, percentOv)) + geom_point(col = noshadGainsAndLosses$type) + facet_wrap(vars(chr, type))
# ggplot(noshadGainsAndLosses, aes(numMarks, percentOv)) + geom_point(col = noshadGainsAndLosses$type) + facet_wrap(vars(ploidy, chr))
# 
# 
# ggplot(noshadGainsAndLosses_filt, aes(sizeMb, percentOv)) + geom_point(col = noshadGainsAndLosses_filt$type) + facet_wrap(vars(ploidy, chr))

# noshadGainsAndLosses_gains <- noshadGainsAndLosses[which(noshadGainsAndLosses$type == "darkred"), ]
# noshadGainsAndLosses_losses <- noshadGainsAndLosses[which(noshadGainsAndLosses$type == "darkblue"), ]
# 
# qualN_g <- ggplot(noshadGainsAndLosses_gains, aes(sample, percentOv)) + geom_boxplot() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# qualN_l <- ggplot(noshadGainsAndLosses_losses, aes(sample, percentOv)) + geom_boxplot() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# gridExtra::grid.arrange(qualN_g, qualN_l)
# 


### this minimal region was made using everything ... i.e including aneuploidy calls
### needed equal regions for interpolating breakpoints

# noshadGainsAndLosses_form <- noshadGainsAndLosses[c("sample", "chr", "start", "end",
#                                           "K", "percentOv")]
# colnames(noshadGainsAndLosses_form) <- c("sampleID", "chrom", "start.pos",
#                                     "end.pos", "n.probes", "mean")
# noshadGainsAndLosses_form$n.probes <- NA
# noshadGainsAndLosses_form$chrom <- str_remove(noshadGainsAndLosses_form$chrom, "chr")
# noshadGainsAndLosses_form <- noshadGainsAndLosses_form[which(noshadGainsAndLosses_form$chrom %in% as.character(1:19)), ]
# noshadGainsAndLosses_form$sampleID <- as.character(noshadGainsAndLosses_form$sampleID)
# 
# noshadGainsAndLosses_freq <- getFreqData(noshadGainsAndLosses_form)
# 
# minComReg <- NULL
# for (i in unique(noshadGainsAndLosses_freq$chr)) {
#   tmpDf <- noshadGainsAndLosses_freq[which(noshadGainsAndLosses_freq$chr == i), ]
#   for (j in 1:(nrow(tmpDf) - 1)) {
#     minComReg <- rbind(minComReg, c(tmpDf$chr[j], tmpDf$pos[j], tmpDf$pos[(j + 1)] - 1))
#   }
# }
# minComReg <- data.frame(minComReg, stringsAsFactors = FALSE)
# colnames(minComReg) <- c("chr", "start", "end")

# write.table(minComReg, "/mnt/DATA5/tmp/kev/misc/20220714minimalCommonRegionCNVex.txt", sep = "\t", row.names = FALSE,
#             col.names = TRUE, quote = FALSE)



### make min common region first, then redo call skipping aneuploidies
### 20220624: after reran taking out aneuploidy regions focusing strictly on  large CNAs
### 20220714: include aneuploidy - makes sense b/c it includes everything and FP and TP calls

minComReg <- read.table("/mnt/DATA5/tmp/kev/misc/20220714minimalCommonRegionCNVex.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE)


# minComReg <- read.table("/mnt/DATA5/tmp/kev/misc/20220706minimalCommonRegionCNVex.txt",
#                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

minComReg$meanPercentOv <- 0
minComReg$medPercentOv <- 0
minComReg$countOv <- 0
minComReg$propZero <- 0

# minComReg$meanPercentOvQual <- 0
# minComReg$countOvQual <- 0



minComRegGrange <- GRanges(seqnames = paste0("chr", minComReg$chr), 
                           IRanges(start = as.numeric(minComReg$start), end = as.numeric(minComReg$end)))

# zeroCounts <- NULL
i <- 50
for (i in 1:length(minComRegGrange)) {
  tmpOv <- findOverlaps(minComRegGrange[i], noshadGrange)
  if(length(tmpOv) < 1){
    next()
  } else{
    minComReg$meanPercentOv[i] <- mean(noshadGainsAndLosses$percentOv[subjectHits(tmpOv)])
    minComReg$medPercentOv[i] <- median(noshadGainsAndLosses$percentOv[subjectHits(tmpOv)])
    minComReg$countOv[i] <- length(which(noshadGainsAndLosses$percentOv[subjectHits(tmpOv)] != 0))
    minComReg$propZero[i] <- length(which(noshadGainsAndLosses$percentOv[subjectHits(tmpOv)] == 0))/48
    # zeroCounts <- c(zeroCounts, noshadGainsAndLosses$sample[subjectHits(tmpOv)][which(noshadGainsAndLosses$percentOv[subjectHits(tmpOv)] == 0)])
  }
}


# quantile(table(zeroCounts), seq(0, 1, 0.05))
quantile(minComReg$propZero, seq(0, 1, 0.05))
quantile(minComReg$propZero, seq(0, 1, 0.05))
# quantile(minComReg$meanPercentOv, seq(0, 1, 0.05))
# quantile(minComReg$medPercentOv, seq(0, 1, 0.05))

# lowQualSamps <- names(table(zeroCounts))[which(table(zeroCounts) > 60)]

# zscore <- (table(zeroCounts) - mean(table(zeroCounts)))/sd(table(zeroCounts))
# zeroMad <- 0.6745 * (table(zeroCounts) - median(table(zeroCounts)))/mad(table(zeroCounts))
# lowQualSamps <- names(which(abs(zeroMad) > 2))



# noshadGainsAndLosses_goodQual <- noshadGainsAndLosses[-which(noshadGainsAndLosses$sample %in% lowQualSamps), ]
# # noshadGainsAndLosses_goodQual <- noshadGainsAndLosses
# testGrange2 <- GRanges(seqnames = noshadGainsAndLosses_goodQual$chr,
#                        IRanges(noshadGainsAndLosses_goodQual$start,
#                                noshadGainsAndLosses_goodQual$end))
# 
# for (i in 1:length(minComRegGrange)) {
#   tmpOv <- findOverlaps(minComRegGrange[i], testGrange2)
#   if(length(tmpOv) < 1){
#     next()
#   } else{
#     minComReg$meanPercentOvQual[i] <- mean(noshadGainsAndLosses_goodQual$percentOv[subjectHits(tmpOv)])
#     minComReg$countOvQual[i] <- length(noshadGainsAndLosses_goodQual$percentOv[subjectHits(tmpOv)])
#     }
# }

minComReg5 <- minComReg
# quantile(minComReg5$meanPercentOvQual, seq(0, 1, 0.05))
# quantile(minComReg5$medPercentOv, seq(0, 1, 0.05))
# minComReg5 <- minComReg5[which(minComReg5$meanPercentOvQual > 0.5), ]

quantile(minComReg5$propZero, seq(0, 1, 0.01))

minComReg5 <- minComReg5[which(minComReg5$propZero < 0.05), ]
minComReg5$size <- (as.numeric(minComReg5$end) - as.numeric(minComReg5$start))/1e6

for (i in unique(minComReg5$chr)) {
  print(paste(i, sum(minComReg5$size[which(minComReg5$chr == i)])))
}




minComRegGrange2 <- GRanges(seqnames = paste0("chr", minComReg5$chr), 
                            IRanges(start = as.numeric(minComReg5$start), end = as.numeric(minComReg5$end)))

zeroCounts <- NULL
for (i in 1:length(minComRegGrange2)) {
  tmpOv <- findOverlaps(minComRegGrange2[i], noshadGrange)
  if(length(tmpOv) < 1){
    next()
  } else{
    zeroCounts <- c(zeroCounts, noshadGainsAndLosses$sample[subjectHits(tmpOv)][which(noshadGainsAndLosses$percentOv[subjectHits(tmpOv)] == 0)])
  }
}

zscore <- (table(zeroCounts) - mean(table(zeroCounts)))/sd(table(zeroCounts))
lowQualSamps <- names(which(abs(zscore) > 2))

### next I'll basically only call these minimal regions for gains and losses

noshadGainsAndLosses1Mb <- noshadGainsAndLosses[which(noshadGainsAndLosses$sizeMb > 10), ]
minRegionNoshad <- minComReg5[, 1:3]
uniqSamps <- unique(noshadGainsAndLosses1Mb$sample)

minRegionNoshadGrange <- GRanges(seqnames = paste0("chr", minRegionNoshad$chr),
                                 IRanges(minRegionNoshad$start, minRegionNoshad$end))

emptyRes <- matrix(NA, nrow = nrow(minRegionNoshad), ncol = length(uniqSamps))
colnames(emptyRes) <- uniqSamps
i <- 1
for (i in 1:length(uniqSamps)) {
  tmpDf <- noshadGainsAndLosses1Mb[which(noshadGainsAndLosses1Mb$sample == uniqSamps[i]),]
  tmpGrange <- GRanges(seqnames = tmpDf$chr, IRanges(start = tmpDf$start, end = tmpDf$end))
  
  tmpOv <- findOverlaps(minRegionNoshadGrange, tmpGrange)
  emptyRes[queryHits(tmpOv), i] <- tmpDf$percentOv[subjectHits(tmpOv)] * sign(tmpDf$totalSc_n[subjectHits(tmpOv)])
}

emptyRes_melt <- melt(emptyRes)
emptyRes_melt <- emptyRes_melt[-which(is.na(emptyRes_melt$value)),]
emptyRes_melt <- emptyRes_melt[, 2:3]
colnames(emptyRes_melt) <- c("sample", "percentOv")
emptyRes_melt$sample <- as.character(emptyRes_melt$sample)

emptyRes_melt$sample <- factor(emptyRes_melt$sample, levels = unique(emptyRes_melt$sample))
emptyRes_melt$type <- ifelse(sign(emptyRes_melt$percentOv) == 1, "darkred", "darkblue")
emptyRes_melt$percentOv <- abs(emptyRes_melt$percentOv)

# axisColors <- ifelse(unique(emptyRes_melt$sample) %in% lowQualSamps, "red", "black")
axisColors <- rep("black", length(unique(emptyRes_melt$sample)))
# diffPloidy <- allKappDf$sample[which(allKappDf$cor < 0.5 | allKappDf$pearsonChanges < 0)]
# axisColors[grep(paste0(diffPloidy, collapse = "|"), unique(as.character(emptyRes_melt$sample)))] <- "goldenrod"

qualN2 <- ggplot(emptyRes_melt, aes(sample, percentOv)) + geom_boxplot() + 
  geom_jitter(col = emptyRes_melt$type, size= 1 , alpha=0.8) + ylim(c(0, 1)) +
  ggtitle("Large CNAs (> 10 Mb ; Proportion Zero < 0.05)") + ylab("Percent overlap") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = axisColors), plot.title = element_text(hjust = 0.5))

qualN2

noshadGainsAndLosses2 <- noshadGainsAndLosses
noshadGainsAndLosses2$sample <- factor(noshadGainsAndLosses2$sample, levels = unique(noshadGainsAndLosses2$sample))
qualN <- ggplot(noshadGainsAndLosses2, aes(sample, percentOv)) + geom_boxplot() + 
  geom_jitter(col = noshadGainsAndLosses2$type, size= 1 , alpha=0.8) +
  ggtitle("Large CNAs") + ylab("Percent overlap") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))



# noshadGainsAndLosses_orig <- noshadGainsAndLosses
# noshadGainsAndLosses_orig <- noshadGainsAndLosses_orig[which(noshadGainsAndLosses_orig$sizeMb > 10), ]
# noshadGainsAndLosses_orig$sample <- factor(noshadGainsAndLosses_orig$sample, levels = unique(emptyRes_melt$sample))

# qualN <- ggplot(noshadGainsAndLosses_orig, aes(sample, percentOv)) + geom_boxplot() + 
#   geom_jitter(size= 1 , alpha=0.8) +
#   ggtitle("Large CNAs (> 5Mb)") + ylab("Percent overlap") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = axisColors), plot.title = element_text(hjust = 0.5))

dev.off()
# png("/mnt/DATA5/tmp/kev/snpVsPanel/20220707LargeCnas.png", width = 1500, height = 1000)
png("/mnt/DATA5/tmp/kev/snpVsPanel/20220718LargeCnas.png", width = 1500, height = 1000)
grid.arrange(qualN, qualN2)
dev.off()


dev.off()
# png("/mnt/DATA5/tmp/kev/snpVsPanel/20220707largeCnasBinomial.png", width = 500, height = 800)
ggplot(noshadGainsAndLosses_orig, aes(x = type, y = percentOv)) + geom_violin() + 
  geom_jitter(col = noshadGainsAndLosses_orig$type,size= 1 , alpha=0.9) +
  ylab("Percent overlap") + xlab("Gains(Red), Losses (Blue)") + 
  theme(axis.text.x = element_blank()) + theme_bw()
dev.off()



### testing instead of agreement using pearson correlation
### could also try using bins

noshadAneploidy_pear <- matrix(0, nrow = length(unique(panelV3$sample)), ncol = length(unique(allNoshadSeg_filt$chr)))
rownames(noshadAneploidy_pear) <- unique(panelV3$sample)
colnames(noshadAneploidy_pear) <- unique(allNoshadSeg_filt$chr)

### use vector of cutoffs
i <- rownames(noshadAneploidy_pear)[46]
for (i in rownames(noshadAneploidy_pear)) {
  tmpDf <- allNoshadSeg_filt[grep(i, allNoshadSeg_filt$sample),]
  tmpPloidy <- as.numeric(panelV3$V3[which(panelV3$sample == i)])
  tmpDf$totalSc_n <- tmpDf$sC- tmpPloidy
  
  j <- 3
  for (j in seq_along(tmpChrList)) {
    tmpChr <- tmpDf[which(tmpDf$chr == tmpChrList[j]), ]
    tmpChr
    chrThres <- sum(tmpChr$size) * vectorOfCutoffs[j]
    tmpAneu <- 0
    if (sum(tmpChr$size[which(tmpChr$totalSc_n > 0.2)]) > chrThres) {
      tmpWeights <- tmpChr$size[which(tmpChr$totalSc_n > 0.2)]/sum(tmpChr$size[which(tmpChr$totalSc_n > 0.2)])
      tmpAneu <- weighted.mean(tmpChr$totalSc_n[which(tmpChr$totalSc_n > 0.2)], tmpWeights)
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < -0.2)]) > chrThres){
      tmpWeights <- tmpChr$size[which(tmpChr$totalSc_n < -0.2)]/sum(tmpChr$size[which(tmpChr$totalSc_n < -0.2)])
      tmpAneu <- weighted.mean(tmpChr$totalSc_n[which(tmpChr$totalSc_n < -0.2)])
    } else{
      tmpAneu <- 0
    }
    noshadAneploidy_pear[i,j] <- tmpAneu
  }
}


ascatAneploidy_pear <- matrix(0, nrow = length(unique(snpRes05$sample)), ncol = length(unique(allNoshadSeg_filt$chr)))
rownames(ascatAneploidy_pear) <- unique(snpRes05$sample)
colnames(ascatAneploidy_pear) <- unique(allNoshadSeg_filt$chr)
tmpChrList <- paste0("chr", c(1:19))

for (i in rownames(ascatAneploidy_pear)) {
  tmpDf <- allAscatSegs_filt[grep(i, allAscatSegs_filt$sample),]
  tmpPloidy <- snpRes05$ploidy_int[which(snpRes05$sample == i)]
  tmpDf$totalSc_n <- tmpDf$totalSc - tmpPloidy
  j <- "chr3"
  for (j in tmpChrList) {
    tmpChr <- tmpDf[which(tmpDf$chr == j), ]
    chrThres <- sum(tmpChr$size) * 0.7
    tmpAneu <- 0
    if (sum(tmpChr$size[which(tmpChr$totalSc_n > 0)]) > chrThres) {
      tmpWeights <- tmpChr$size[which(tmpChr$totalSc_n > 0)]/sum(tmpChr$size[which(tmpChr$totalSc_n > 0)])
      tmpAneu <- weighted.mean(tmpChr$totalSc_n[which(tmpChr$totalSc_n > 0)], tmpWeights)
    } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < 0)]) > chrThres){
      tmpWeights <- tmpChr$size[which(tmpChr$totalSc_n < 0)]/sum(tmpChr$size[which(tmpChr$totalSc_n < 0)])
      tmpAneu <- weighted.mean(tmpChr$totalSc_n[which(tmpChr$totalSc_n < 0)], tmpWeights)
    } else{
      tmpAneu <- 0
    }
    ascatAneploidy_pear[i,j] <- tmpAneu
  }
}



pearsonAll <- NULL
i <- rownames(ascatAneploidy)[2]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy_pear[i, 1:19]
  tmpNoshad <- noshadAneploidy_pear[i, 1:19]
  tmpNoshad[which(tmpNoshad > 0.2 & tmpNoshad < 0.8)] <- tmpAscat[which(tmpNoshad > 0.2 & tmpNoshad < 0.8)]
  tmpNoshad[which(tmpNoshad < -0.2 & tmpNoshad > -0.8)] <- tmpAscat[which(tmpNoshad < -0.2 & tmpNoshad > -0.8)]
  pearsonAll <- c(pearsonAll, cor(tmpAscat, tmpNoshad))
}

pearsonChanges <- NULL
i <- rownames(ascatAneploidy)[16]
for (i in rownames(ascatAneploidy)) {
  tmpAscat <- ascatAneploidy_pear[i, 1:19]
  tmpNoshad <- noshadAneploidy_pear[i, 1:19]
  
  tmpNoshad[which(tmpNoshad > 0 & tmpNoshad < 0.8)] <- tmpAscat[which(tmpNoshad > 0 & tmpNoshad < 0.8)]
  tmpNoshad[which(tmpNoshad < 0 & tmpNoshad > -0.8)] <- tmpAscat[which(tmpNoshad < 0 & tmpNoshad > -0.8)]
  tmpIdx <- c(which(tmpAscat != 0), which(tmpNoshad != 0))
  
  if (length(tmpIdx) < 1) {
    pearsonChanges <- c(pearsonChanges, NA)
  } else if (length(which(duplicated(tmpIdx))) < 1) {
    pearsonChanges  <- c(pearsonChanges , cor(tmpAscat[tmpIdx], tmpNoshad[tmpIdx]))
  } else{
    tmpIdx <- tmpIdx[-which(duplicated(tmpIdx))]
    pearsonChanges  <- c(pearsonChanges , cor(tmpAscat[tmpIdx], tmpNoshad[tmpIdx]))
  }
}

allKappDf$pearsonChanges <- pearsonChanges

### pearson anueploidy graphs

pearsonChanges_m_dip <-  pearsonChanges[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 2)]
pearsonChanges_m_trip <-  pearsonChanges[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 3)]
pearsonChanges_m_tetra <-  pearsonChanges[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 4)]
pearsonChanges_um_dip <-  pearsonChanges[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 2)]
pearsonChanges_um_trip <-  pearsonChanges[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 3)]
pearsonChanges_um_tetra <-  pearsonChanges[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 4)]

alterations_all <- c(rep("aneu_m_2n", length(pearsonChanges_m_dip)),
                     rep("aneu_m_3n", length(pearsonChanges_m_trip)),
                     rep("aneu_m_4n", length(pearsonChanges_m_tetra)),
                     rep("aneu_um_2n", length(pearsonChanges_um_dip)),
                     rep("aneu_um_3n", length(pearsonChanges_um_trip)),
                     rep("aneu_um_4n", length(pearsonChanges_um_tetra)))

concordance_all <- c(pearsonChanges_m_dip, pearsonChanges_m_trip, pearsonChanges_m_tetra,
                     pearsonChanges_um_dip, pearsonChanges_um_trip, pearsonChanges_um_tetra)


ggplotDf_m_p_all <- data.frame("alterations" = alterations_all, "concordance" = concordance_all)
ggplotDf_m_p_all$alterations <- factor(ggplotDf_m_p_all$alterations, levels = unique(ggplotDf_m_p_all$alterations))

pa_graph_all <- ggplot(data = ggplotDf_m_p_all, aes(x = alterations, y  = concordance)) + geom_boxplot() +
  geom_jitter(color="black", size= 1 , alpha=0.9) +
  stat_summary(fun.data = n_fun, geom = "text") + ylim(c(0,1.2)) + 
  ylab("percentile agreement") + xlab("type of alterations") + 
  ggtitle("Pearson's R Aneuploidy: NGS(Var) & SNP(0.7)") + theme(plot.title = element_text(hjust = 0.5))


pearsonAll_m_dip <-  pearsonAll[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 2)]
pearsonAll_m_trip <-  pearsonAll[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 3)]
pearsonAll_m_tetra <-  pearsonAll[which(snpRes05$ploidy_int == panelV3$V3 & panelV3$V3 == 4)]
pearsonAll_um_dip <-  pearsonAll[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 2)]
pearsonAll_um_trip <-  pearsonAll[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 3)]
pearsonAll_um_tetra <-  pearsonAll[which(snpRes05$ploidy_int != panelV3$V3 & panelV3$V3 == 4)]


alterations_all2 <- c(rep("all_m_2n", length(pearsonAll_m_dip)),
                      rep("all_m_3n", length(pearsonAll_m_trip)),
                      rep("all_m_4n", length(pearsonAll_m_tetra)),
                      rep("all_um_2n", length(pearsonAll_um_dip)),
                      rep("all_um_3n", length(pearsonAll_um_trip)),
                      rep("all_um_4m", length(pearsonAll_um_tetra)))

concordance_all2 <- c(pearsonAll_m_dip, pearsonAll_m_trip, pearsonAll_m_tetra,
                      pearsonAll_um_dip, pearsonAll_um_trip, pearsonAll_um_tetra)

ggplotDf_ck_all <- data.frame("alterations" = alterations_all2, "ckap" = concordance_all2)
ggplotDf_ck_all$alterations <- factor(ggplotDf_ck_all$alterations, levels = unique(ggplotDf_ck_all$alterations))




ck_graph_all <- ggplot(data = ggplotDf_ck_all, aes(x = alterations, y  = ckap)) + geom_boxplot() + 
  geom_jitter(color="black", size= 1 , alpha=0.9) +
  stat_summary(fun.data = n_fun, geom = "text") + ylim(c(0,1.2)) + 
  ggtitle("All chromosomes Pearson's R: NGS(Var) & SNP(0.7)") + theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("match + ploidy")

gridExtra::grid.arrange(ck_graph_all, pa_graph_all, ncol = 2)

allKappDf$cor <- pearsonAll

dev.off()
png("/mnt/DATA5/tmp/kev/snpVsPanel/20220718AneuploidyByPloidyPearson.png", width = 1000, height = 400)
gridExtra::grid.arrange(ck_graph_all, pa_graph_all, ncol = 2)
dev.off()

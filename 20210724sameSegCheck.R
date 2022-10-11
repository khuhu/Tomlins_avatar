countAneuCna_tcga_test <- function(df, ploidyDf){
  res <- NULL
  # this was genius and got rid of the loop - literally infinitely faster *pats self on back*
  df$newTotalCN <- round2(df$Copynumber - ploidyDf$ploidy[match(df$Sample, ploidyDf$array)], 0)
  
  
  for(i in unique(df$Sample)){
    resMat <- matrix(data = 0, nrow = 1, ncol = 4) # aneuploidy gain, loss: CNA 15mb gain, loss
    sampDf <- df[which(df$Sample == i),]
    human_arm <- read.table("/mnt/DATA5/tmp/kev/misc/20210713human_arm_syn.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = TRUE)
    
    
    ### counting for anueploidy
    for (j in unique(human_arm$chr)) {
      skipVar <- "no"
      df_chr <- sampDf[sampDf$Chromosome == j,]
      df_chr <- df_chr[order(df_chr$Start),]
      minLength80 <- min(human_arm$length[which(human_arm$chromosome == j)]) * 0.8
      minLength <- 15000000
      
      
      # new loop below gets segments of gains or losses
      delIdx <- which(sign(df_chr$newTotalCN) == -1)
      ampIdx <- which(sign(df_chr$newTotalCN) == 1)
      
      # thought experiment - this actually undercounts because of arm
      # level events - overcount cnas, if I don't skip here
      if(length(delIdx) == 0 & length(ampIdx) == 0){
        next()
      }
      
      if (sum(df_chr$length[delIdx]) > minLength80) {
        resMat[1,2] <- resMat[1,2] + 1
        skipVar <- "yes"
        print(paste0("chr", j, "loss:", resMat[1,2]))
        # next()
        # print(paste(i, j))
      } 
      
      if (sum(df_chr$length[ampIdx]) > minLength80) {
        resMat[1,1] <- resMat[1,1] + 1
        skipVar <- "yes"
        print(paste0("chr", j, "gain", resMat[1,1]))
        # next()
        # print(paste(i, j))
      }
      
      if (skipVar == "yes") {
        print(paste0("chrom", j, "skipped"))
        next() 
      }
      
      
      ### cna counting
      
      cn_sign <- sign(df_chr$newTotalCN)
      k <- 1
      idx1 <- 0
      idx2 <- 0
      tmpSign <- 0
      while (k < (length(cn_sign) + 1)) {
        if (cn_sign[k] == 0 & idx1 == 0) { #start + non end string of zeros
          k <- k + 1
          print(1)
          next()
        } else if(cn_sign[k] != 0 & idx1 == 0 & k == length(cn_sign)){ # if there is a single segment at the end
          idx1 <- length(cn_sign)
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          #print(tmpSign <- cn_sign[k])
          if (tmpLength > minLength & tmpLength < minLength80) {
            if (tmpSign == 1) {
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$Chromosome[idx1], ":",df_chr$Start[idx1], "-",df_chr$End[idx2]))
          }
          k <- k + 1
          print(2)
          next()
        } else if(cn_sign[k] != 0 & idx1 == 0){
          idx1 <- k
          tmpSign <- cn_sign[k]
          k <- k + 1
          print(3)
          next()
        } else if(idx1 != 0 & tmpSign == cn_sign[k] & k == length(cn_sign)){
          idx2 <- k
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength & tmpLength <  minLength80) {
            if (tmpSign == 1) {
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              resMat[1,4] <- resMat[1,4] + 1
            }
          }
          k <- k + 1
          print("long-end")
          print(paste0(i,"/","chr", df_chr$Chromosome[idx1], ":",df_chr$Start[idx1], "-",df_chr$End[idx2]))
          next()
        } else if (idx1 != 0 & tmpSign == cn_sign[k] & k != length(cn_sign)){
          idx2 <- k
          k <- k + 1
          print(4)
          next()
        } else if(cn_sign[k] == (tmpSign * -1) & idx2 == 0){ # special case where -1 and 1 are adjacent
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          print(paste0(tmpLength, minLength))
          if (tmpLength > minLength & tmpLength < minLength80) {
            if (tmpSign == 1) {
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$Chromosome[idx1], ":",df_chr$Start[idx1], "-",df_chr$End[idx2]))
          }
          
          if (k == length(cn_sign)) {
            idx1 <- k
            idx2 <- k
            tmpSign <- cn_sign[k]
            if (tmpSign == 1) {
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$Chromosome[idx1], ":",df_chr$Start[idx1], "-",df_chr$End[idx2]))
          }
          
          idx1 <- k
          idx2 <- 0
          tmpSign <- cn_sign[k]
          k <- k + 1
          print("adjacent")
          next()
        } else if(tmpSign != cn_sign[k] & idx2 == 0){
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength & tmpLength < minLength80) {
            if (tmpSign == 1) {
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$Chromosome[idx1], ":",df_chr$Start[idx1], "-",df_chr$End[idx2]))
          }
          idx1 <- 0
          idx2 <- 0
          k <- k + 1
          print(5)
          next()
        } else if(cn_sign[k] == (tmpSign * -1) & idx2 != 0){
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength & tmpLength < minLength80) {
            if (tmpSign == 1) {
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$Chromosome[idx1], ":",df_chr$Start[idx1], "-",df_chr$End[idx2]))
          }
          idx1 <- k
          tmpSign <- cn_sign[k]
          k <- k + 1
          idx2 <- 0
          print("adjacent 2")
          next()
        } else if(tmpSign != cn_sign[k] & idx2 != 0){
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength & tmpLength < minLength80) {
            if (tmpSign == 1) {
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$Chromosome[idx1], ":",df_chr$Start[idx1], "-",df_chr$End[idx2]))
          }
          idx1 <- 0
          idx2 <- 0
          k <- k + 1
          print(6)
          next()
        }
      }
    }
    
    rownames(resMat) <- i
    res <- rbind(res, resMat)
  }
  return(res)
}



separateSegments_mV2 <- function(df){
  res <- NULL
  res2 <- NULL
  df$length <- df$end.pos - df$start.pos
  
  for (i in unique(df$sampleID)) {
    sampDf <- df[which(df$sampleID == i),]
    for (j in unique(sampDf$chrom)) {
      df_chr <- sampDf[sampDf$chrom == j,]
      minLength <- sum(df_chr$length) * 0.8
      aTable <- df_chr
      cTable <- df_chr
      
      cn_sign <- sign(df_chr$mean)
      
      idx1 <- NULL
      idx2 <- NULL
      
      i <- 1
      idx1 <- 0
      idx2 <- 0
      tmpSign <- 0
      while (i < (length(cn_sign) + 1)) {
        if (cn_sign[i] == 0 & idx1 == 0) { #start + non end string of zeros
          i <- i + 1
          print(1)
          next()
        } else if(cn_sign[i] != 0 & idx1 == 0 & i == length(cn_sign)){ # if there is a single segment at the end
          idx1 <- length(cn_sign)
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          i <- i + 1
          print(2)
          next()
        } else if(cn_sign[i] != 0 & idx1 == 0){
          idx1 <- i
          tmpSign <- cn_sign[i]
          i <- i + 1
          print(3)
          next()
        } else if(idx1 != 0 & tmpSign == cn_sign[i] & i == length(cn_sign)){
          idx2 <- i
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          i <- i + 1
          print("long-end")
          next()
        } else if (idx1 != 0 & tmpSign == cn_sign[i] & i != length(cn_sign)){
          idx2 <- i
          i <- i + 1
          print(4)
          next()
        } else if(cn_sign[i] == (tmpSign * -1) & idx2 == 0){ # special case where -1 and 1 are adjacent
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          
          
          idx1 <- i
          tmpSign <- cn_sign[i]
          i <- i + 1
          print("adjacent")
          next()
        } else if(tmpSign != cn_sign[i] & idx2 == 0){
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          idx1 <- 0
          idx2 <- 0
          i <- i + 1
          print(5)
          next()
        } else if(cn_sign[i] == (tmpSign * -1) & idx2 != 0){
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          idx1 <- i + 1
          tmpSign <- cn_sign[i]
          i <- i + 1
          idx2 <- 0
          print("adjacent 2")
          next()
        } else if(tmpSign != cn_sign[i] & idx2 != 0){
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          idx1 <- 0
          i <- i + 1
          print(6)
          next()
        }
      }
      res <- rbind(res, aTable)
      res2 <- rbind(res2, cTable)
    }
  }
  return(list(res, res2))
}


segResults_1 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/segResults.txt",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
segResults_2 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/segResults.txt",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
segResults_3 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/segResults.txt",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
segResults_4 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-76-MG_test1_255_185/segResults.txt",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")


c(unique(segResults_1$ID), unique(segResults_2$ID),
  unique(segResults_3$ID), unique(segResults_4$ID))[duplicated(c(unique(segResults_1$ID), unique(segResults_2$ID),
             unique(segResults_3$ID), unique(segResults_4$ID)))]

zscore_gc_oe_ratios <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/gcCorrectedCounts_matrix.txt", sep = "\t",
                                  stringsAsFactors = FALSE, header = TRUE)

mouseNormal <- c("MG_17X49", "MG_18X50", "MG_23X55", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

all_segRes <- rbind(segResults_1, segResults_2, segResults_3, segResults_4)
all_segRes <- all_segRes[-which(all_segRes$ID %in% mouseNormal),]
all_segRes$ID <- str_remove(all_segRes$ID, "_MG_X.*")
all_segRes$ID <- str_remove(str_remove(all_segRes$ID, "^X"), "_X.*")
all_segRes$ID <- tolower(str_remove_all(str_remove_all(str_remove(all_segRes$ID, "X.*"), "_"), "\\."))
all_segRes$ID  <- str_remove(all_segRes$ID, "o")

annoTable$mouse_id <- str_remove(str_remove(str_remove_all(tolower(annoTable$mouse_id), "_"),"-"), "o")
annoTable$type <- tolower(annoTable$type)
annoTable$type[142] <- "ehgsc"

all_segRes <- all_segRes[which(all_segRes$ID %in% annoTable$mouse_id[which(annoTable$geno %in% c("BPRN", "BPN", "BPP"))]),]


matchingTc <- tcDf$tc[match(tolower(all_segRes$ID), tcDf$sample)]
all_segRes$seg.mean[which(!is.na(matchingTc))] <- all_segRes$seg.mean[which(!is.na(matchingTc))]/matchingTc[which(!is.na(matchingTc))]


segZscores <- calcZscore(all_segRes)
segZfilt <- segZscoresFilt_zeroOut(all_segRes, segZscores)

tmpSeg <- cbind(segZfilt[,1:4], NA, segZfilt[,5])
colnames(tmpSeg) <- c("sampleID","chrom", "start.pos","end.pos", "n.probes", "mean")
tmpSeg$string <- paste0(tmpSeg$sampleID, tmpSeg$chrom, tmpSeg$start.pos, tmpSeg$end.pos)
tmpSeg <- tmpSeg[-which(duplicated(tmpSeg$string)),]
tmpSeg <- tmpSeg[order(tmpSeg$chrom, tmpSeg$start.pos),]


tmpSeg2 <- tmpSeg[which(tmpSeg$sampleID == "15767rt" & tmpSeg$chrom == 11),]
tmpSeg2$mean <- 0
armRes <- separateSegments_mV2(tmpSeg2)

tmpSeg2 <- tmpSeg[which(tmpSeg$sampleID == "2027lts" & tmpSeg$chrom == 11),]
tmpSeg2$mean <- c(1,rep(0,6))
armRes <- separateSegments_mV2(tmpSeg2)

tmpSeg2 <- tmpSeg[which(tmpSeg$sampleID == "2027lts" & tmpSeg$chrom == 11),]
tmpSeg2$mean <- c(1, 1, 1, 0, 0, 0, 0)
armRes <- separateSegments_mV2(tmpSeg2)


tmpSeg2 <- tmpSeg[which(tmpSeg$sampleID == "2027lts" & tmpSeg$chrom == 11),]
tmpSeg2$mean <- c(0, 0, 0, 0, 0, 0, 1)
armRes <- separateSegments_mV2(tmpSeg2)

tmpSeg2 <- tmpSeg[which(tmpSeg$sampleID == "2027lts" & tmpSeg$chrom == 11),]
tmpSeg2$mean <- c(1, 0, 1, 1, 1, 0, 1)
armRes <- separateSegments_mV2(tmpSeg2)


tmpSeg2 <- tmpSeg[which(tmpSeg$sampleID == "2027lts" & tmpSeg$chrom == 11),]
tmpSeg2$mean <- c(1, -1, 1, 1, 1, -1, 1)
armRes <- separateSegments_mV2(tmpSeg2)

tmpSeg2 <- tmpSeg[which(tmpSeg$sampleID == "2027lts" & tmpSeg$chrom == 11),]
tmpSeg2$mean <- c(1, -1, -1, 0, 1, -1, 1)
armRes <- separateSegments_mV2(tmpSeg2)


tmpSeg3 <- tmpSeg[which(tmpSeg$sampleID == "14120rt"),]
armRes <- countAneuCna_mouse(tmpSeg3)


### so far so good. check three more with CNAs - aneuploidy count no problem
### then run with the tcga data - make sure to convert the column names 
gi_count_mouse <- countAneuCna_mouse(tmpSeg)
colnames(gi_count_mouse) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")

finalTable2 <- finalTable2[order(finalTable2$Chromosome, finalTable2$Start),]
gi_count_human <- countAneuCna_tcga(finalTable2, tcga_ploidy)
colnames(gi_count_human) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")

# tmpRes <- countAneuCna_tcga_test(finalTable2[which(finalTable2$Sample == "TCGA-OR-A5J4-01"),], tcga_ploidy)
# 


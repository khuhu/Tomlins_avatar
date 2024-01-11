### general pipeline will be reading in a file from list of tables
### input should also be target mouse model (1) amp and del file to create synteny table (2) allMouseAneu_bprn file - needed for synteny perm
### general code will be for each list of files read in as table, process to create synteny table -> permute with target mouse model
### output should be an amp and del table each with the p-values, expected frequencies
### each iteration of the table should annotate the human cancer it makes a comparison with


### create dummy code for 1 iteraiton and then run to see if I'm missing anything. 
### need to make sure final table has everything I need to look at comparison of synteny between cancer types
### do I use single arm measurements or a sum of .... i.e some type of metric showing the uniqueness of this mouse model profile
### goes back to fidelity of the model - this will be the end of the first set of figures i.e a-d
### everything after will be the rna analysis (e to ?) - should consist at looking at rna profiles of arm changes
### which genes are changed in the direction of the arm changes. 




# listOfFiles <- c("gdac.broadinstitute.org_COADREAD-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/broad_values_by_arm.txt")
# annotationTable <- coad_anno2
# mouseAmp <- allMouseAneu_coadAdenoCar_amp_bed
# mouseDel <- allMouseAneu_coadAdenoCar_del_bed
# tmpBroadMouseMelt2 <- allMouseAneu_coadAdenoCar
tcgaSyntenyPerm <- function(listOfFiles, annotationTable = NULL, mouseAmp, mouseDel, tmpBroadMouseMelt2){
  
  dir <- "/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/"
  
  allCans <- c("BLCA", "BRCA", "CESC", "COADREAD", "ESCA", "GBM",
               "GBMLGG", "HNSC", "KIPAN", "KIRC", "LGG", "LIHC",
               "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD",
               "SARC", "STES", "TGCT", "THYM")
  
  tmpAmpRes <- NULL
  tmpDelRes <- NULL
  for (j in seq_along(listOfFiles)) {
    setwd(dir)
    setwd(listOfFiles[j])
    
    tmpBroad <- read.table("broad_values_by_arm.txt", sep = "\t", stringsAsFactors = FALSE,
                           header = TRUE, check.names = FALSE)
    colnames(tmpBroad)[2:ncol(tmpBroad)] <- gsub("(^.*?-.{3}?)-.*", "\\1",
                                                 colnames(tmpBroad)[2:ncol(tmpBroad)])
    if (allCans[j] == "COADREAD") {
      tmpBroad <- tmpBroad[, c(1, grep(paste0(annotationTable$Sample.ID, collapse = "|"), colnames(tmpBroad)))] 
    }
    tmpBroad <- tmpBroad[-grep("X", tmpBroad$`Chromosome Arm`), ]
    tmpBroadMat <- tmpBroad[, 2:ncol(tmpBroad)]
    tmpBroad2 <- cbind(hg19ArmLocations[, c("chrStripped", "start", "end")], tmpBroadMat)
    tmpBroadMelt <- melt(tmpBroad2, id.vars = c("chrStripped", "start", "end"))
    tmpBroadMelt2 <- data.frame("sampleID" = tmpBroadMelt$variable, "chrom" = tmpBroadMelt$chrStripped,
                                "start.pos" = tmpBroadMelt$start, "end.pos" = tmpBroadMelt$end, 
                                "n.probes" = NA, "mean" = tmpBroadMelt$value)
    tmpBroadMelt2$str <- paste0(tmpBroadMelt2$sampleID, tmpBroadMelt2$chrom, 
                                tmpBroadMelt2$start.pos, tmpBroadMelt2$end.pos)
    tmpBroadMelt2$length <- tmpBroadMelt2$end.pos - tmpBroadMelt2$start.pos
    
    tmpArm_freq <- getFreqData(tmpBroadMelt2)
    tmpArm_freq_res <- ampsDels(tmpArm_freq)
    tmpArm_amp_bed <- reducingFreqBed2(tmpArm_freq_res[[1]])
    tmpArm_del_bed <- reducingFreqBed2(tmpArm_freq_res[[3]])
    
    tmpSynTable <- circosFreq2(mouseAmp, mouseDel, tmpArm_amp_bed,
                               tmpArm_del_bed, plot = FALSE, ref = "human")
    
    ### synteny permutation portion
    
    tmpSynTable$str <- paste(tmpSynTable$h_chr,
                             tmpSynTable$h_start,
                             tmpSynTable$h_end)
    tmpArms <- tmpBroad2
    tmpSynTable2 <- tmpSynTable[-which(duplicated(tmpSynTable$str)), ]
    tmpSynTable2$h_chr <- as.numeric(str_remove(tmpSynTable2$h_chr, "h_chr"))
    tmpSynTable2$m_chr <- as.numeric(str_remove(tmpSynTable2$m_chr, "m_chr"))
    tmpSynTable2 <- tmpSynTable2[order(tmpSynTable2$h_chr, tmpSynTable2$h_start, decreasing = FALSE),]
    tmpArmsGr <- GRanges(seqnames = tmpArms$chrStripped, IRanges(start = tmpArms$start, end = tmpArms$end))
    tmpHumanGrange <- GRanges(seqnames = tmpSynTable2$h_chr, IRanges(start = tmpSynTable2$h_start, end = tmpSynTable2$h_end))
    tmpMouseGrange <- GRanges(seqnames = tmpSynTable2$m_chr, IRanges(start = tmpSynTable2$m_start, end = tmpSynTable2$m_end))
    
    tmpHumanRegionMat <- data.frame("chr" = tmpSynTable2$h_chr, "start" = tmpSynTable2$h_start, "end" = tmpSynTable2$h_end)
    for (i in 4:ncol(tmpArms)) {
      
      tmpIdx <- subjectHits(findOverlaps(query = tmpHumanGrange, subject = tmpArmsGr))
      testQuery <- queryHits(findOverlaps(query = tmpHumanGrange, subject = tmpArmsGr))
      
      
      if (length(which(duplicated(testQuery))) > 0) {
        tmpIdx <-  tmpIdx[-which(duplicated(testQuery))]
        testQuery <- testQuery[-which(duplicated(testQuery))]
      }
      
      ### error probably happens b/c there is no match for one of the queries
      tmpHumanRegionMat[ ,i] <- NA
      tmpHumanRegionMat[ ,i] <- tmpArms[tmpIdx, i]
    }
    tmpHumanRegionMat2 <- tmpHumanRegionMat[ , 4:ncol(tmpHumanRegionMat)]
    
    cl <- makeCluster(25)
    registerDoParallel(cl)
    start <- Sys.time()
    tmpHumanPerm <- foreach(i=1:25,
                           .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                             tmpAmpAll <- NULL
                             tmpDelAll <- NULL
                             tmpAmpFreq <- NULL
                             tmpDelFreq <- NULL
                             j <- 0
                             while (j < 400) {
                               tmp <- apply(tmpHumanRegionMat2, 2, function(x) sample(x, nrow(tmpHumanRegionMat2)))
                               ampFreq <- apply(tmp, 1, function(x) length(which(x > 0.2))/length(x))
                               delFreq <- apply(tmp, 1, function(x) length(which(x < -0.2))/length(x))
                               tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0.2)]) * length(which(x > 0.2))/length(x))
                               tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < -0.2)])) * length(which(x < -0.2))/length(x))
                               tmpAmpAll <- cbind(tmpAmpAll, tmpAmp)
                               tmpDelAll <- cbind(tmpDelAll, tmpDel)
                               tmpAmpFreq <- cbind(tmpAmpFreq, ampFreq)
                               tmpDelFreq <- cbind(tmpDelFreq, delFreq)
                               j <- j + 1
                             }
                             return(list(tmpAmpAll, tmpDelAll, tmpAmpFreq, tmpDelFreq))
                           }
    
    
    stopCluster(cl)
    print( Sys.time() - start )

    tmpMouseRegionMat <- data.frame("chr" = tmpSynTable2$m_chr, "start" = tmpSynTable2$m_start, "end" = tmpSynTable2$m_end)
    for (i in seq_along(unique(tmpBroadMouseMelt2$sampleID))) {
      
      tmp <- tmpBroadMouseMelt2[which(tmpBroadMouseMelt2$sampleID == unique(tmpBroadMouseMelt2$sampleID)[i]),]
      tmpGrange <- GRanges(seqnames = tmp$chrom, IRanges(start = tmp$start.pos, end = tmp$end.pos))
      
      
      tmpIdx <- subjectHits(findOverlaps(query = tmpMouseGrange, subject = tmpGrange))
      testQuery <- queryHits(findOverlaps(query = tmpMouseGrange, subject = tmpGrange))
      
      if (length(which(duplicated(testQuery))) > 0) {
        tmpIdx <-  tmpIdx[-which(duplicated( testQuery))]
        testQuery <- testQuery[-which(duplicated( testQuery))]
      }
      ### error probably happens b/c there is no match for one of the queries
      
      tmpMouseRegionMat[ , i + 3] <- NA
      tmpMouseRegionMat[testQuery, i + 3] <- tmp$mean[tmpIdx]
    }
    
    tmpMouseRegionMat2 <-tmpMouseRegionMat[ , 4:ncol(tmpMouseRegionMat)]
    cl <- makeCluster(25)
    registerDoParallel(cl)
    start <- Sys.time()
    tmpMousePerm <- foreach(i=1:25,
                             .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                               tmpAmpAll <- NULL
                               tmpDelAll <- NULL
                               tmpAmpFreq <- NULL
                               tmpDelFreq <- NULL
                               j <- 0
                               while (j < 400) {
                                 tmp <- apply(tmpMouseRegionMat2, 2, function(x) sample(x, nrow(tmpMouseRegionMat2)))
                                 ampFreq <- apply(tmp, 1, function(x) length(which(x > 0.2))/length(x))
                                 delFreq <- apply(tmp, 1, function(x) length(which(x < -0.2))/length(x))
                                 tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0.2)]) * length(which(x > 0.2))/length(x))
                                 tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < -0.2)])) * length(which(x < -0.2))/length(x))
                                 tmpAmpAll <- cbind(tmpAmpAll, tmpAmp)
                                 tmpDelAll <- cbind(tmpDelAll, tmpDel)
                                 tmpAmpFreq <- cbind(tmpAmpFreq, ampFreq)
                                 tmpDelFreq <- cbind(tmpDelFreq, delFreq)
                                 j <- j + 1
                               }
                               return(list(tmpAmpAll, tmpDelAll, tmpAmpFreq, tmpDelFreq))
                             }
    
    
    stopCluster(cl)
    print( Sys.time() - start )
    
    mousePermAmpFreq <- do.call(cbind, tmpMousePerm[[3]])
    mousePermDelFreq <- do.call(cbind, tmpMousePerm[[4]])
    
    humanPermAmpFreq <- do.call(cbind, tmpHumanPerm[[3]])
    humanPermDelFreq <- do.call(cbind, tmpHumanPerm[[4]])
    
    tmpAmpSynTable <- data.frame("chr" = tmpSynTable2$m_chr, "start" = tmpSynTable2$m_start, "end" = tmpSynTable2$m_end)
    tmpDelSynTable <- data.frame("chr" = tmpSynTable2$m_chr, "start" = tmpSynTable2$m_start, "end" = tmpSynTable2$m_end)
    tmpAmpSynTable$h_freq <- apply(tmpHumanRegionMat2, 1, function(x) length(which(x > 0.2))/length(x))
    tmpAmpSynTable$m_freq <- apply(tmpMouseRegionMat2, 1, function(x) length(which(x > 0.2))/length(x))
    tmpDelSynTable$h_freq <- apply(tmpHumanRegionMat2, 1, function(x) length(which(x < -0.2))/length(x))
    tmpDelSynTable$m_freq <- apply(tmpMouseRegionMat2, 1, function(x) length(which(x < -0.2))/length(x))

    tmpHumanPermAmp2 <- humanPermAmpFreq
    tmpMousePermAmp2 <- mousePermAmpFreq
    tmpHumanPermDel2 <- humanPermDelFreq
    tmpMousePermDel2 <- mousePermDelFreq
    
    tmpPermStatTbl <- NULL
    for (i in 1:396) {
      
      # tmpAmp <- (tmpHumanPermAmp2[i, ] + tmpMousePermAmp2[i, ])/2
      # tmpDel <- (tmpHumanPermDel2[i, ] + tmpMousePermDel2[i, ])/2
      
      tmpAmp <- (tmpHumanPermAmp2[i, ] * tmpMousePermAmp2[i, ])
      tmpDel <- (tmpHumanPermDel2[i, ] * tmpMousePermDel2[i, ])
      
      if (tmpAmpSynTable$h_freq[i] == 0 | tmpAmpSynTable$m_freq[i] == 0) {
        tmpAmpZ <- 0
      } else{
        # tmpAmpZ <- ((tmpAmpSynTable$h_freq[i] + tmpAmpSynTable$m_freq[i]) - mean(tmpAmp, na.rm = TRUE))/sd(tmpAmp, na.rm = TRUE)
        tmpAmpZ <- ((tmpAmpSynTable$h_freq[i] * tmpAmpSynTable$m_freq[i]) - mean(tmpAmp, na.rm = TRUE))/sd(tmpAmp, na.rm = TRUE)
        
      }
      
      if (tmpDelSynTable$h_freq[i] == 0 | tmpDelSynTable$m_freq[i] == 0) {
        tmpDelZ <- 0
      } else{
        # tmpDelZ <- ((tmpDelSynTable$h_freq[i] + tmpDelSynTable$m_freq[i]) -  mean(tmpDel, na.rm = TRUE))/sd(tmpDel, na.rm = TRUE)
        tmpDelZ <- ((tmpDelSynTable$h_freq[i] * tmpDelSynTable$m_freq[i]) -  mean(tmpDel, na.rm = TRUE))/sd(tmpDel, na.rm = TRUE)
        }
      
      tmpPermStatTbl <- rbind(tmpPermStatTbl, data.frame("ampZ" = tmpAmpZ, "delZ" = tmpDelZ))
    }
    
    
    tmpPermStatTbl$ampP <- pnorm(q=tmpPermStatTbl$ampZ, lower.tail=FALSE)
    tmpPermStatTbl$delP <- pnorm(q=tmpPermStatTbl$delZ, lower.tail=FALSE)
    
    tmpResTable_amp <- tmpSynTable2
    tmpResTable_amp$h_freq <- tmpAmpSynTable$h_freq
    tmpResTable_amp$m_freq <- tmpAmpSynTable$m_freq
    tmpResTable_amp$ampP <- tmpPermStatTbl$ampP
    tmpResTable_amp$adjP <- p.adjust(tmpResTable_amp$ampP, method = "BH", n = 2 * length(tmpResTable_amp$ampP))
    tmpResTable_amp$freqStat <- (tmpResTable_amp$h_freq * tmpResTable_amp$m_freq)
    tmpResTable_amp$cancer <- allCans[j]
    
    tmpResTable_del <- tmpSynTable2
    tmpResTable_del$h_freq <- tmpDelSynTable$h_freq
    tmpResTable_del$m_freq <- tmpDelSynTable$m_freq
    tmpResTable_del$delP <- tmpPermStatTbl$delP
    tmpResTable_del$adjP2<- p.adjust(tmpResTable_del$delP, method = "BH", n = 2 * length(tmpResTable_del$delP))
    tmpResTable_del$freqStat <- (tmpResTable_del$h_freq * tmpResTable_del$m_freq)
    tmpResTable_del$cancer <- allCans[j]
    
    tmpAmpRes <- rbind(tmpAmpRes, tmpResTable_amp)
    tmpDelRes <- rbind(tmpDelRes, tmpResTable_del)
  }
  return(list(tmpAmpRes, tmpDelRes))
}


dir <- "/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/"
setwd(dir)
allDirs <- list.dirs()
allDirs <- allDirs[grep("TP", allDirs)]

listOfFiles <- allDirs
allCancerSynteny <- tcgaSyntenyPerm(listOfFiles = allDirs, annotationTable = coad_anno2,
                                    mouseAmp = allMouseAneu_coadAdenoCar_amp_bed,
                                    mouseDel = allMouseAneu_coadAdenoCar_del_bed,
                                    tmpBroadMouseMelt2 = allMouseAneu_coadAdenoCar)

allCancerSynGainDf <- allCancerSynteny[[1]]
allCancerSynLossDf <- allCancerSynteny[[2]]

ggplot(allCancerSynGainDf) + geom_point(aes(x = m_freq, y = -log10(adjP))) + facet_wrap(~cancer)

ggplot(allCancerSynGainDf) + geom_boxplot(aes(x = cancer, y = -log10(adjP))) 

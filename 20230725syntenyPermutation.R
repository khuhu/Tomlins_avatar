### permutation method:
### have two matrices where it's sample (row) by region (col)
### permute every row (sample) of each matrix n times and calculate frequncy of gained/loss for each column (region)
### these set of values represent for each region the background distribution of freqeuncy of region gain/loss



### first make the region by 
hgsc_bprn_arm_allSynTable$str <- paste(hgsc_bprn_arm_allSynTable$h_chr, hgsc_bprn_arm_allSynTable$h_start, hgsc_bprn_arm_allSynTable$h_end)
tp53OvArms <- broadBySampleGistic2
tmpSynTable <- hgsc_bprn_arm_allSynTable[-which(duplicated(hgsc_bprn_arm_allSynTable$str)), ]
tmpSynTable$h_chr <- as.numeric(str_remove(tmpSynTable$h_chr, "h_chr"))
tmpSynTable$m_chr <- as.numeric(str_remove(tmpSynTable$m_chr, "m_chr"))
tmpSynTable <- tmpSynTable[order(tmpSynTable$h_chr, tmpSynTable$h_start, decreasing = FALSE),]
tp53OvArmsGr <- GRanges(seqnames = tp53OvArms$chrStripped, IRanges(start = tp53OvArms$start, end = tp53OvArms$end))
ovSyntenyGrange <- GRanges(seqnames = tmpSynTable$h_chr, IRanges(start = tmpSynTable$h_start, end = tmpSynTable$h_end))
ovMm10SyntenyGrange <- GRanges(seqnames = tmpSynTable$m_chr, IRanges(start = tmpSynTable$m_start, end = tmpSynTable$m_end))
tmpSynTable$stat <- tmpSynTable$h_freq + tmpSynTable$m_freq

humanOvRegionMat <- data.frame("chr" = tmpSynTable$h_chr, "start" = tmpSynTable$h_start, "end" = tmpSynTable$h_end)
for (i in 4:ncol(tp53OvArms)) {
  
  tmpIdx <- subjectHits(findOverlaps(query = ovSyntenyGrange, subject = tp53OvArmsGr))
  testQuery <- queryHits(findOverlaps(query = ovSyntenyGrange, subject = tp53OvArmsGr))
  
  
  if (length(which(duplicated(testQuery))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated(testQuery))]
    testQuery <- testQuery[-which(duplicated(testQuery))]
  }
  
  ### error probably happens b/c there is no match for one of the queries
  humanOvRegionMat[ ,i] <- NA
  humanOvRegionMat[ ,i] <- tp53OvArms[tmpIdx, i]
}



mouseOvRegionMat <- data.frame("chr" = tmpSynTable$m_chr, "start" = tmpSynTable$m_start, "end" = tmpSynTable$m_end)
for (i in seq_along(unique(allMouseAneu_bprn$sampleID))) {
  
  tmp <- allMouseAneu_bprn[which(allMouseAneu_bprn$sampleID == unique(allMouseAneu_bprn$sampleID)[i]),]
  tmpGrange <- GRanges(seqnames = tmp$chrom, IRanges(start = tmp$start.pos, end = tmp$end.pos))
  
  print(unique(allMouseAneu_bprn$sampleID)[i])
  
  tmpIdx <- subjectHits(findOverlaps(query = ovMm10SyntenyGrange, subject = tmpGrange))
  testQuery <- queryHits(findOverlaps(query = ovMm10SyntenyGrange, subject = tmpGrange))
  
  if (length(which(duplicated(testQuery))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated( testQuery))]
    testQuery <- testQuery[-which(duplicated( testQuery))]
  }
  ### error probably happens b/c there is no match for one of the queries
  
  mouseOvRegionMat[ , i + 3] <- NA
  mouseOvRegionMat[testQuery, i + 3] <- tmp$mean[tmpIdx]
}

### 20230808 changing everything to chunks i.e almost as fast as each additional cpu

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


humanOvRegionMat2 <- humanOvRegionMat[ , 4:ncol(humanOvRegionMat)]
cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
synHgscPerm <- foreach(i=1:25,
                    .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                      tmpAmpAll <- NULL
                      tmpDelAll <- NULL
                      tmpAmpFreq <- NULL
                      tmpDelFreq <- NULL
                      j <- 0
                      while (j < 400) {
                        tmp <- apply(humanOvRegionMat2, 2, function(x) sample(x, nrow(humanOvRegionMat2)))
                        ampFreq <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
                        delFreq <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
                        tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
                        tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
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

humanPermAmpStat2 <- do.call(cbind, synHgscPerm[[1]])
humanPermDelStat2 <- do.call(cbind, synHgscPerm[[2]])
humanPermAmp <- do.call(cbind, synHgscPerm[[3]])
humanPermDel <- do.call(cbind, synHgscPerm[[4]])


# start <- Sys.time()
# humanOvRegionMat2 <- humanOvRegionMat[ , 4:ncol(humanOvRegionMat)]
# humanPermAmp <- NULL
# humanPermDel <- NULL
# humanPermAmpStat2 <- NULL
# humanPermDelStat2 <- NULL
# for (i in 1:10000) {
#   tmp <- apply(humanOvRegionMat2, 2, function(x) sample(x, 396))
#   tmpAmp <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
#   tmpDel <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
#   tmpAmpStat2 <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
#   tmpDelStat2 <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
#   
#   humanPermAmp <- cbind(humanPermAmp, tmpAmp)
#   humanPermDel <- cbind(humanPermDel, tmpDel)
#   humanPermAmpStat2 <- cbind(humanPermAmpStat2,  tmpAmpStat2)
#   humanPermDelStat2 <- cbind(humanPermDelStat2, tmpDelStat2)
# }
# print( Sys.time() - start )
# colnames(humanPermDel) <- NULL
# colnames(humanPermAmp) <- NULL
# colnames(humanPermDelStat2) <- NULL
# colnames(humanPermAmpStat2) <- NULL

mouseOvRegionMat2 <- mouseOvRegionMat[ , 4:ncol(mouseOvRegionMat)]
cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
synBprnPerm <- foreach(i=1:25,
                       .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                         tmpAmpAll <- NULL
                         tmpDelAll <- NULL
                         tmpAmpFreq <- NULL
                         tmpDelFreq <- NULL
                         j <- 0
                         while (j < 400) {
                           tmp <- apply(mouseOvRegionMat2, 2, function(x) sample(x, nrow(mouseOvRegionMat2)))
                           ampFreq <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
                           delFreq <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
                           tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
                           tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
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

mousePermAmpStat2 <- do.call(cbind, synBprnPerm[[1]])
mousePermDelStat2 <- do.call(cbind, synBprnPerm[[2]])
mousePermAmp <- do.call(cbind, synBprnPerm[[3]])
mousePermDel <- do.call(cbind, synBprnPerm[[4]])

# start <- Sys.time()
# mouseOvRegionMat2 <- mouseOvRegionMat[ , 4:ncol(mouseOvRegionMat)]
# mousePermAmp <- NULL
# mousePermDel <- NULL
# mousePermAmpStat2 <- NULL
# mousePermDelStat2 <- NULL
# for (i in 1:10000) {
#   tmp <- apply(mouseOvRegionMat2, 2, function(x) sample(x, 396))
#   tmpAmp <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
#   tmpDel <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
#   tmpAmpStat2 <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
#   tmpDelStat2 <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
#   
#   mousePermAmp <- cbind(mousePermAmp, tmpAmp)
#   mousePermDel <- cbind(mousePermDel, tmpDel)
#   mousePermAmpStat2 <- cbind(mousePermAmpStat2, tmpAmpStat2)
#   mousePermDelStat2 <- cbind(mousePermDelStat2, tmpDelStat2)
# }
# print( Sys.time() - start )


ampSynTable <- data.frame("chr" = tmpSynTable$m_chr, "start" = tmpSynTable$m_start, "end" = tmpSynTable$m_end)
delSynTable <- data.frame("chr" = tmpSynTable$m_chr, "start" = tmpSynTable$m_start, "end" = tmpSynTable$m_end)
ampSynTable$h_freq <- apply(humanOvRegionMat2, 1, function(x) length(which(x > 0))/length(x))
ampSynTable$m_freq <- apply(mouseOvRegionMat2, 1, function(x) length(which(x > 0))/length(x))
delSynTable$h_freq <- apply(humanOvRegionMat2, 1, function(x) length(which(x < 0))/length(x))
delSynTable$m_freq <- apply(mouseOvRegionMat2, 1, function(x) length(which(x < 0))/length(x))
ampSynTable$h_stat2 <- apply(humanOvRegionMat2, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
ampSynTable$m_stat2 <- apply(mouseOvRegionMat2, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
delSynTable$h_stat2 <- apply(humanOvRegionMat2, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
delSynTable$m_stat2 <- apply(mouseOvRegionMat2, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))



ampSynTable$h_freq[which(ampSynTable$h_freq < 0.1)] <- 0
ampSynTable$m_freq[which(ampSynTable$m_freq < 0.1)] <- 0
delSynTable$h_freq[which(delSynTable$h_freq < 0.1)] <- 0
delSynTable$m_freq[which(delSynTable$m_freq < 0.1)] <- 0
ampSynTable$h_stat2[which(ampSynTable$h_freq < 0.1)] <- 0
ampSynTable$m_stat2[which(ampSynTable$m_freq < 0.1)] <- 0
delSynTable$h_stat2[which(delSynTable$h_freq < 0.1)] <- 0
delSynTable$m_stat2[which(delSynTable$m_freq < 0.1)] <- 0


humanPermAmp2 <- humanPermAmp
mousePermAmp2 <- mousePermAmp
humanPermDel2 <- humanPermDel
mousePermDel2 <- mousePermDel
humanPermAmpStat2v2 <- humanPermAmpStat2
mousePermAmpStat2v2 <- mousePermAmpStat2
humanPermDelStat2v2 <- humanPermDelStat2
mousePermDelStat2v2 <- mousePermDelStat2

humanPermAmp2[humanPermAmp2 < 0.1] <- 0
mousePermAmp2[mousePermAmp2 < 0.1] <- 0
humanPermDel2[humanPermDel2 < 0.1] <- 0
mousePermDel2[mousePermDel2 < 0.1] <- 0
humanPermAmpStat2v2[humanPermAmp2 < 0.1] <- 0
mousePermAmpStat2v2[mousePermAmp2 < 0.1] <- 0
humanPermDelStat2v2[humanPermDel2 < 0.1] <- 0
mousePermDelStat2v2[mousePermDel2 < 0.1] <- 0


permStatResTbl <- NULL
for (i in 1:396) {
  tmpAmp <- humanPermAmp2[i, ] + mousePermAmp2[i, ]
  tmpDel <- humanPermDel2[i, ] + mousePermDel2[i, ]
  
  tmpAmpStat2 <- humanPermAmpStat2v2[i, ] + mousePermAmpStat2v2[i, ]
  tmpDelStat2 <- humanPermDelStat2v2[i, ] + mousePermDelStat2v2[i, ]
  
  if (ampSynTable$h_freq[i] == 0 | ampSynTable$m_freq[i] == 0) {
    tmpAmpZ <- 0
    tmpAmpZStat2 <- 0
  } else{
    tmpAmpZ <- ((ampSynTable$h_freq[i] + ampSynTable$m_freq[i]) - mean(tmpAmp, na.rm = TRUE))/sd(tmpAmp, na.rm = TRUE)
    tmpAmpZStat2 <- ((ampSynTable$h_stat2[i] + ampSynTable$m_stat2[i]) - mean(tmpAmpStat2, na.rm = TRUE))/sd(tmpAmpStat2, na.rm = TRUE)
  }
  
  if (delSynTable$h_freq[i] == 0 | delSynTable$m_freq[i] == 0) {
    tmpDelZ <- 0
    tmpDelZStat2 <- 0
  } else{
    tmpDelZ <- ((delSynTable$h_freq[i] + delSynTable$m_freq[i]) - mean(tmpAmp, na.rm = TRUE))/sd(tmpAmp, na.rm = TRUE)
    tmpDelZStat2 <- ((delSynTable$h_stat2[i] + delSynTable$m_stat2[i]) - mean(tmpDelStat2, na.rm = TRUE))/sd(tmpDelStat2, na.rm = TRUE)
  }
  
  permStatResTbl <- rbind(permStatResTbl, data.frame("ampZ" = tmpAmpZ, "delZ" = tmpDelZ, "ampStat2Z" = tmpAmpZStat2, "delStat2Z" = tmpDelZStat2))
}



### only looking at frequencies greater than what we expected from permutation test
# permStatResTbl$ampP <- pnorm(q=permStatResTbl$ampZ, lower.tail=FALSE)
# permStatResTbl$delP <- pnorm(q=permStatResTbl$delZ, lower.tail=FALSE)
permStatResTbl$ampP2 <- pnorm(q=permStatResTbl$ampStat2Z, lower.tail=FALSE)
permStatResTbl$delP2 <- pnorm(q=permStatResTbl$delStat2Z, lower.tail=FALSE)

# zeroRegionIdx <- unique(c(which(permStatResTbl$ampP == 0.5), which(permStatResTbl$delP == 0.5)))
# zeroRegionIdx <- zeroRegionIdx[order(zeroRegionIdx)]
# zeroRegionIdx2 <- unique(c(which(permStatResTbl$ampP2 == 0.5), which(permStatResTbl$delP2 == 0.5)))
# zeroRegionIdx2 <- zeroRegionIdx2[order(zeroRegionIdx2)]

# permStatResTbl2 <- permStatResTbl[-zeroRegionIdx, ]
# permStatResTbl2$ampPAdj <- p.adjust(permStatResTbl2$ampP, method = "BH", n = 2 * length(permStatResTbl2$ampP))
# permStatResTbl2$delPAdj <- p.adjust(permStatResTbl2$delP, method = "BH", n = 2 * length(permStatResTbl2$delP))
# 
# tmpSynTable2 <- tmpSynTable[-zeroRegionIdx, ]
# hgsc_bprn_arm_allSynTable_amp <- tmpSynTable2[which(permStatResTbl2$ampPAdj < 0.05), ]
# hgsc_bprn_arm_allSynTable_del <- tmpSynTable2[which(permStatResTbl2$delPAdj < 0.05), ]

# tmpSynTable2 <- tmpSynTable[-zeroRegionIdx, ]
# hgsc_bprn_arm_allSynTable_amp <- tmpSynTable2[which(permStatResTbl2$ampPAdj < 0.05), ]
# hgsc_bprn_arm_allSynTable_del <- tmpSynTable2[which(permStatResTbl2$delPAdj < 0.05), ]


permStatResTbl_amp <- permStatResTbl
permStatResTbl_del <- permStatResTbl
permStatResTbl_amp$ampPAdj2 <- p.adjust(permStatResTbl_amp$ampP2, method = "BH", n = (length(permStatResTbl_amp$ampP2) + length(permStatResTbl_del$delP2)))
permStatResTbl_del$delPAdj2 <- p.adjust(permStatResTbl_del$delP2, method = "BH", n = (length(permStatResTbl_amp$ampP2) + length(permStatResTbl_del$delP2)))

tmpSynTable_amp <- tmpSynTable
tmpSynTable_del <- tmpSynTable
hgsc_bprn_arm_allSynTable_amp <- tmpSynTable_amp[which(permStatResTbl_amp$ampPAdj2 < 0.05), ]
hgsc_bprn_arm_allSynTable_del <- tmpSynTable_del[which(permStatResTbl_del$delPAdj2 < 0.05), ]


cancerGeneCensusGr <- GRanges(seqnames = cancerGeneCensus$chr,
                              IRanges(start = as.numeric(cancerGeneCensus$start), end = as.numeric(cancerGeneCensus$end)))
cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bprn_arm_allSynTable_amp$h_chr, IRanges(hgsc_bprn_arm_allSynTable_amp$h_start, end = hgsc_bprn_arm_allSynTable_amp$h_end))
                                                     ,cancerGeneCensusGr)),]
cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bprn_arm_allSynTable_del$h_chr, IRanges(hgsc_bprn_arm_allSynTable_del$h_start, end = hgsc_bprn_arm_allSynTable_del$h_end))
                                          ,cancerGeneCensusGr)),]


### doing something similar but for other genotypes
###
###

hgsc_bpp_arm_allSynTable$str <- paste(hgsc_bpp_arm_allSynTable$h_chr, hgsc_bpp_arm_allSynTable$h_start, hgsc_bpp_arm_allSynTable$h_end)
tmpSynTable_bpp <- hgsc_bpp_arm_allSynTable[-which(duplicated(hgsc_bpp_arm_allSynTable$str)), ]
tmpSynTable_bpp$h_chr <- as.numeric(str_remove(tmpSynTable_bpp$h_chr, "h_chr"))
tmpSynTable_bpp$m_chr <- as.numeric(str_remove(tmpSynTable_bpp$m_chr, "m_chr"))
tmpSynTable_bpp <- tmpSynTable_bpp[order(tmpSynTable_bpp$h_chr, tmpSynTable_bpp$h_start, decreasing = FALSE),]



mouseOvRegionMat_bpp <- data.frame("chr" = tmpSynTable_bpp$m_chr, "start" = tmpSynTable_bpp$m_start, "end" = tmpSynTable_bpp$m_end)
for (i in seq_along(unique(allMouseAneu_bpp$sampleID))) {
  
  tmp <- allMouseAneu_bpp[which(allMouseAneu_bpp$sampleID == unique(allMouseAneu_bpp$sampleID)[i]),]
  tmpGrange <- GRanges(seqnames = tmp$chrom, IRanges(start = tmp$start.pos, end = tmp$end.pos))
  
  print(unique(allMouseAneu_bpp$sampleID)[i])
  
  tmpIdx <- subjectHits(findOverlaps(query = ovMm10SyntenyGrange, subject = tmpGrange))
  testQuery <- queryHits(findOverlaps(query = ovMm10SyntenyGrange, subject = tmpGrange))
  
  if (length(which(duplicated(testQuery))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated( testQuery))]
    testQuery <- testQuery[-which(duplicated( testQuery))]
  }
  ### error probably happens b/c there is no match for one of the queries
  
  mouseOvRegionMat_bpp[ , i + 3] <- NA
  mouseOvRegionMat_bpp[testQuery, i + 3] <- tmp$mean[tmpIdx]
}


mouseOvRegionMat2_bpp <- mouseOvRegionMat_bpp[ , 4:ncol(mouseOvRegionMat_bpp)]
cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
synBppPerm <- foreach(i=1:25,
                       .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                         tmpAmpAll <- NULL
                         tmpDelAll <- NULL
                         tmpAmpFreq <- NULL
                         tmpDelFreq <- NULL
                         j <- 0
                         while (j < 400) {
                           tmp <- apply(mouseOvRegionMat2_bpp, 2, function(x) sample(x, nrow(mouseOvRegionMat2_bpp)))
                           ampFreq <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
                           delFreq <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
                           tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
                           tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
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

mousePermAmpStat_bpp <- do.call(cbind, synBppPerm[[1]])
mousePermDelStat_bpp <- do.call(cbind, synBppPerm[[2]])
mousePermAmp_bpp <- do.call(cbind, synBprnPerm[[3]])
mousePermDel_bpp <- do.call(cbind, synBprnPerm[[4]])

# start <- Sys.time()
# mouseOvRegionMat2_bpp <- mouseOvRegionMat_bpp[ , 4:ncol(mouseOvRegionMat_bpp)]
# mousePermAmp_bpp <- NULL
# mousePermDel_bpp <- NULL
# mousePermAmpStat_bpp <- NULL
# mousePermDelStat_bpp <- NULL
# for (i in 1:10000) {
#   tmp <- apply(mouseOvRegionMat2_bpp, 2, function(x) sample(x, 396))
#   tmpAmp <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
#   tmpDel <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
#   tmpAmpStat2 <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
#   tmpDelStat2 <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
#   
#   mousePermAmp_bpp <- cbind(mousePermAmp_bpp, tmpAmp)
#   mousePermDel_bpp <- cbind(mousePermDel_bpp, tmpDel)
#   mousePermAmpStat_bpp <- cbind(mousePermAmpStat_bpp, tmpAmpStat2)
#   mousePermDelStat_bpp <- cbind(mousePermDelStat_bpp, tmpDelStat2)
#   
# }
# print( Sys.time() - start )


ampSynTable_bpp <- data.frame("chr" = tmpSynTable_bpp$m_chr, "start" = tmpSynTable_bpp$m_start, "end" = tmpSynTable_bpp$m_end)
delSynTable_bpp <- data.frame("chr" = tmpSynTable_bpp$m_chr, "start" = tmpSynTable_bpp$m_start, "end" = tmpSynTable_bpp$m_end)
ampSynTable_bpp$h_freq <- apply(humanOvRegionMat2, 1, function(x) length(which(x > 0))/length(x))
ampSynTable_bpp$m_freq <- apply(mouseOvRegionMat2_bpp, 1, function(x) length(which(x > 0))/length(x))
delSynTable_bpp$h_freq <- apply(humanOvRegionMat2, 1, function(x) length(which(x < 0))/length(x))
delSynTable_bpp$m_freq <- apply(mouseOvRegionMat2_bpp, 1, function(x) length(which(x < 0))/length(x))
ampSynTable_bpp$h_stat <- apply(humanOvRegionMat2, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
ampSynTable_bpp$m_stat <- apply(mouseOvRegionMat2_bpp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
delSynTable_bpp$h_stat <- apply(humanOvRegionMat2, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
delSynTable_bpp$m_stat <- apply(mouseOvRegionMat2_bpp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))


ampSynTable_bpp$h_freq[which(ampSynTable_bpp$h_freq < 0.1)] <- 0
ampSynTable_bpp$m_freq[which(ampSynTable_bpp$m_freq < 0.1)] <- 0
delSynTable_bpp$h_freq[which(delSynTable_bpp$h_freq < 0.1)] <- 0
delSynTable_bpp$m_freq[which(delSynTable_bpp$m_freq < 0.1)] <- 0
ampSynTable_bpp$h_stat[which(ampSynTable_bpp$h_freq < 0.1)] <- 0
ampSynTable_bpp$m_stat[which(ampSynTable_bpp$m_freq < 0.1)] <- 0
delSynTable_bpp$h_stat[which(delSynTable_bpp$h_freq < 0.1)] <- 0
delSynTable_bpp$m_stat[which(delSynTable_bpp$m_freq < 0.1)] <- 0

humanPermAmp2_bpp <- humanPermAmp
mousePermAmp2_bpp <- mousePermAmp_bpp
humanPermDel2_bpp <- humanPermDel
mousePermDel2_bpp <- mousePermDel_bpp
mousePermAmpStat2_bpp <- mousePermAmpStat_bpp
mousePermDelStat2_bpp <- mousePermDelStat_bpp
humanPermAmpStat2v2 <- humanPermAmpStat2
humanPermDelStat2v2 <- humanPermDelStat2

humanPermAmp2_bpp[humanPermAmp2_bpp < 0.1] <- 0
mousePermAmp2_bpp[mousePermAmp2_bpp < 0.1] <- 0
humanPermDel2_bpp[humanPermDel2_bpp < 0.1] <- 0
mousePermDel2_bpp[mousePermDel2_bpp < 0.1] <- 0
humanPermAmpStat2v2[humanPermAmp2 < 0.1] <- 0
humanPermDelStat2v2[humanPermDel2 < 0.1] <- 0
mousePermAmpStat2_bpp[mousePermAmp2_bpp < 0.1] <- 0
mousePermDelStat2_bpp[mousePermDel2_bpp < 0.1] <- 0


permStatResTbl_bpp <- NULL
for (i in 1:396) {
  tmpAmp <- humanPermAmp2_bpp[i, ] + mousePermAmp2_bpp[i, ]
  tmpDel <- humanPermDel2_bpp[i, ] + mousePermDel2_bpp[i, ]
  
  tmpAmpStat2 <- humanPermAmpStat2v2[i, ] + mousePermAmpStat2_bpp[i, ]
  tmpDelStat2 <- humanPermDelStat2v2[i, ] + mousePermDelStat2_bpp[i, ]
  
  if (ampSynTable_bpp$h_freq[i] == 0 | ampSynTable_bpp$m_freq[i] == 0) {
    tmpAmpZ <- 0
    tmpAmpZStat2 <- 0
  } else{
    tmpAmpZ <- ((ampSynTable_bpp$h_freq[i] + ampSynTable_bpp$m_freq[i]) - mean(tmpAmp, na.rm = TRUE))/sd(tmpAmp, na.rm = TRUE)
    tmpAmpZStat2 <- ((ampSynTable$h_stat2[i] + ampSynTable$m_stat2[i]) - mean(tmpAmpStat2, na.rm = TRUE))/sd(tmpAmpStat2, na.rm = TRUE)
  }
  
  if (delSynTable_bpp$h_freq[i] == 0 | delSynTable_bpp$m_freq[i] == 0) {
    tmpDelZ <- 0
    tmpDelZStat2 <- 0
  } else{
    tmpDelZ <- ((delSynTable_bpp$h_freq[i] + delSynTable_bpp$m_freq[i]) -  mean(tmpDel, na.rm = TRUE))/sd(tmpDel, na.rm = TRUE)
    tmpDelZStat2 <- ((delSynTable$h_stat2[i] + delSynTable$m_stat2[i]) - mean(tmpDelStat2, na.rm = TRUE))/sd(tmpDelStat2, na.rm = TRUE)
  }
  
  permStatResTbl_bpp <- rbind(permStatResTbl_bpp, data.frame("ampZ" = tmpAmpZ, "delZ" = tmpDelZ, "ampStat2Z" = tmpAmpZStat2, "delStat2Z" = tmpDelZStat2))
}


### only looking at frequencies greater than what we expected from permutation test
# permStatResTbl_bpp$ampP <- pnorm(q=permStatResTbl_bpp$ampZ, lower.tail=FALSE)
# permStatResTbl_bpp$delP <- pnorm(q=permStatResTbl_bpp$delZ, lower.tail=FALSE)
# zeroRegionIdx <- unique(c(which(permStatResTbl_bpp$ampP == 0.5), which(permStatResTbl_bpp$delP == 0.5)))
# zeroRegionIdx <- zeroRegionIdx[order(zeroRegionIdx)]
# 
# permStatResTbl2_bpp <- permStatResTbl_bpp[-zeroRegionIdx, ]
# permStatResTbl2_bpp$ampPAdj <- p.adjust(permStatResTbl2_bpp$ampP, method = "BH", n = 2 * length(permStatResTbl2_bpp$ampP))
# permStatResTbl2_bpp$delPAdj <- p.adjust(permStatResTbl2_bpp$delP, method = "BH", n = 2 * length(permStatResTbl2_bpp$delP))
# 
# tmpSynTable2_bpp <- tmpSynTable_bpp[-zeroRegionIdx, ]
# 
# hgsc_bpp_arm_allSynTable_amp <- tmpSynTable2_bpp[which(permStatResTbl2_bpp$ampPAdj < 0.05), ]
# hgsc_bpp_arm_allSynTable_del <- tmpSynTable2_bpp[which(permStatResTbl2_bpp$delPAdj < 0.05), ]

permStatResTbl_bpp$ampP2 <- pnorm(q=permStatResTbl_bpp$ampStat2Z, lower.tail=FALSE)
permStatResTbl_bpp$delP2 <- pnorm(q=permStatResTbl_bpp$delStat2Z, lower.tail=FALSE)

permStatResTbl_bpp_amp <- permStatResTbl_bpp
permStatResTbl_bpp_del <- permStatResTbl_bpp

permStatResTbl_bpp_amp$ampPAdj2 <- p.adjust(permStatResTbl_bpp_amp$ampP2, method = "BH", n = (length(permStatResTbl_bpp_amp$ampP2) + length(permStatResTbl_bpp_del$ampP2)))
permStatResTbl_bpp_del$delPAdj2 <- p.adjust(permStatResTbl_bpp_del$delP2, method = "BH", n = (length(permStatResTbl_bpp_amp$ampP2) + length(permStatResTbl_bpp_del$ampP2)))

tmpSynTable2_bpp_amp <- tmpSynTable_bpp
tmpSynTable2_bpp_del <- tmpSynTable_bpp

hgsc_bpp_arm_allSynTable_amp <- tmpSynTable2_bpp_amp[which(permStatResTbl_bpp_amp$ampPAdj2 < 0.05), ]
hgsc_bpp_arm_allSynTable_del <- tmpSynTable2_bpp_del[which(permStatResTbl_bpp_del$delPAdj2 < 0.05), ]

cancerGeneCensusGr <- GRanges(seqnames = cancerGeneCensus$chr,
                              IRanges(start = as.numeric(cancerGeneCensus$start), end = as.numeric(cancerGeneCensus$end)))
cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bpp_arm_allSynTable_amp$h_chr, IRanges(hgsc_bpp_arm_allSynTable_amp$h_start, end = hgsc_bpp_arm_allSynTable_amp$h_end))
                                          ,cancerGeneCensusGr)),]
cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bpp_arm_allSynTable_del$h_chr, IRanges(hgsc_bpp_arm_allSynTable_del$h_start, end = hgsc_bpp_arm_allSynTable_del$h_end))
                                          ,cancerGeneCensusGr)),]


###  bpn 
###
###

hgsc_bpn_arm_allSynTable$str <- paste(hgsc_bpn_arm_allSynTable$h_chr, hgsc_bpn_arm_allSynTable$h_start, hgsc_bpn_arm_allSynTable$h_end)
tmpSynTable_bpn <- hgsc_bpn_arm_allSynTable[-which(duplicated(hgsc_bpn_arm_allSynTable$str)), ]
tmpSynTable_bpn$h_chr <- as.numeric(str_remove(tmpSynTable_bpn$h_chr, "h_chr"))
tmpSynTable_bpn$m_chr <- as.numeric(str_remove(tmpSynTable_bpn$m_chr, "m_chr"))
tmpSynTable_bpn <- tmpSynTable_bpn[order(tmpSynTable_bpn$h_chr, tmpSynTable_bpn$h_start, decreasing = FALSE),]


mouseOvRegionMat_bpn <- data.frame("chr" = tmpSynTable_bpn$m_chr, "start" = tmpSynTable_bpn$m_start, "end" = tmpSynTable_bpn$m_end)
for (i in seq_along(unique(allMouseAneu_bpn$sampleID))) {
  
  tmp <- allMouseAneu_bpn[which(allMouseAneu_bpn$sampleID == unique(allMouseAneu_bpn$sampleID)[i]),]
  tmpGrange <- GRanges(seqnames = tmp$chrom, IRanges(start = tmp$start.pos, end = tmp$end.pos))
  
  print(unique(allMouseAneu_bprn$sampleID)[i])
  
  tmpIdx <- subjectHits(findOverlaps(query = ovMm10SyntenyGrange, subject = tmpGrange))
  testQuery <- queryHits(findOverlaps(query = ovMm10SyntenyGrange, subject = tmpGrange))
  
  if (length(which(duplicated(testQuery))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated( testQuery))]
    testQuery <- testQuery[-which(duplicated( testQuery))]
  }
  ### error probably happens b/c there is no match for one of the queries
  
  mouseOvRegionMat_bpn[ , i + 3] <- NA
  mouseOvRegionMat_bpn[testQuery, i + 3] <- tmp$mean[tmpIdx]
}

mouseOvRegionMat2_bpn <- mouseOvRegionMat_bpn[ , 4:ncol(mouseOvRegionMat_bpn)]
cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
synBpnPerm <- foreach(i=1:25,
                      .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                        tmpAmpAll <- NULL
                        tmpDelAll <- NULL
                        tmpAmpFreq <- NULL
                        tmpDelFreq <- NULL
                        j <- 0
                        while (j < 400) {
                          tmp <- apply(mouseOvRegionMat2_bpn, 2, function(x) sample(x, nrow(mouseOvRegionMat2_bpn)))
                          ampFreq <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
                          delFreq <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
                          tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
                          tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
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

mousePermAmpStat_bpn <- do.call(cbind, synBpnPerm[[1]])
mousePermDelStat_bpn <- do.call(cbind, synBpnPerm[[2]])
mousePermAmp_bpn <- do.call(cbind, synBpnPerm[[3]])
mousePermDel_bpn <- do.call(cbind, synBpnPerm[[4]])

# start <- Sys.time()
# mouseOvRegionMat2_bpn <- mouseOvRegionMat_bpn[ , 4:ncol(mouseOvRegionMat_bpn)]
# mousePermAmp_bpn <- NULL
# mousePermDel_bpn <- NULL
# mousePermAmpStat_bpn <- NULL
# mousePermDelStat_bpn <- NULL
# for (i in 1:10000) {
#   tmp <- apply(mouseOvRegionMat2_bpn, 2, function(x) sample(x, 396))
#   tmpAmp <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
#   tmpDel <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
#   tmpAmpStat2 <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
#   tmpDelStat2 <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
#   
#   mousePermAmp_bpn <- cbind(mousePermAmp_bpn, tmpAmp)
#   mousePermDel_bpn <- cbind(mousePermDel_bpn, tmpDel)
#   mousePermAmpStat_bpn <- cbind(mousePermAmpStat_bpn, tmpAmpStat2)
#   mousePermDelStat_bpn <- cbind(mousePermDelStat_bpn, tmpDelStat2)
# }
# print( Sys.time() - start )


ampSynTable_bpn <- data.frame("chr" = tmpSynTable_bpn$m_chr, "start" = tmpSynTable_bpn$m_start, "end" = tmpSynTable_bpn$m_end)
delSynTable_bpn <- data.frame("chr" = tmpSynTable_bpn$m_chr, "start" = tmpSynTable_bpn$m_start, "end" = tmpSynTable_bpn$m_end)
ampSynTable_bpn$h_freq <- apply(humanOvRegionMat2, 1, function(x) length(which(x > 0))/length(x))
ampSynTable_bpn$m_freq <- apply(mouseOvRegionMat2_bpn, 1, function(x) length(which(x > 0))/length(x))
delSynTable_bpn$h_freq <- apply(humanOvRegionMat2, 1, function(x) length(which(x < 0))/length(x))
delSynTable_bpn$m_freq <- apply(mouseOvRegionMat2_bpn, 1, function(x) length(which(x < 0))/length(x))
ampSynTable_bpn$h_stat <- apply(humanOvRegionMat2, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
ampSynTable_bpn$m_stat <- apply(mouseOvRegionMat2_bpn, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
delSynTable_bpn$h_stat <- apply(humanOvRegionMat2, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
delSynTable_bpn$m_stat <- apply(mouseOvRegionMat2_bpn, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))

ampSynTable_bpn$h_freq[which(ampSynTable_bpn$h_freq < 0.1)] <- 0
ampSynTable_bpn$m_freq[which(ampSynTable_bpn$m_freq < 0.1)] <- 0
delSynTable_bpn$h_freq[which(delSynTable_bpn$h_freq < 0.1)] <- 0
delSynTable_bpn$m_freq[which(delSynTable_bpn$m_freq < 0.1)] <- 0
ampSynTable_bpn$h_stat[which(ampSynTable_bpn$h_freq < 0.1)] <- 0
ampSynTable_bpn$m_stat[which(ampSynTable_bpn$m_freq < 0.1)] <- 0
delSynTable_bpn$h_stat[which(delSynTable_bpn$h_freq < 0.1)] <- 0
delSynTable_bpn$m_stat[which(delSynTable_bpn$m_freq < 0.1)] <- 0

humanPermAmp2_bpn <- humanPermAmp
mousePermAmp2_bpn <- mousePermAmp_bpn
humanPermDel2_bpn <- humanPermDel
mousePermDel2_bpn <- mousePermDel_bpn
mousePermAmpStat2_bpn <- mousePermAmpStat_bpn
mousePermDelStat2_bpn <- mousePermDelStat_bpn
humanPermAmpStat2v2 <- humanPermAmpStat2
humanPermDelStat2v2 <- humanPermDelStat2

humanPermAmp2_bpn[humanPermAmp2_bpn < 0.1] <- 0
mousePermAmp2_bpn[mousePermAmp2_bpn < 0.1] <- 0
humanPermDel2_bpn[humanPermDel2_bpn < 0.1] <- 0
mousePermDel2_bpn[mousePermDel2_bpn < 0.1] <- 0
humanPermAmpStat2v2[humanPermAmp2 < 0.1] <- 0
humanPermDelStat2v2[humanPermDel2 < 0.1] <- 0
mousePermAmpStat2_bpn[mousePermAmp2_bpn < 0.1] <- 0
mousePermDelStat2_bpn[mousePermDel2_bpn < 0.1] <- 0

permStatResTbl_bpn <- NULL
for (i in 1:396) {
  tmpAmp <- humanPermAmp2_bpn[i, ] + mousePermAmp2_bpn[i, ]
  tmpDel <- humanPermDel2_bpn[i, ] + mousePermDel2_bpn[i, ]
  
  tmpAmpStat2 <- humanPermAmpStat2v2[i, ] + mousePermAmpStat2_bpn[i, ]
  tmpDelStat2 <- humanPermDelStat2v2[i, ] + mousePermDelStat2_bpn[i, ]
  
  if (ampSynTable_bpn$h_freq[i] == 0 | ampSynTable_bpn$m_freq[i] == 0) {
    tmpAmpZ <- 0
    tmpAmpZStat2 <- 0
  } else{
    tmpAmpZ <- ((ampSynTable_bpn$h_freq[i] + ampSynTable_bpn$m_freq[i]) - mean(tmpAmp, na.rm = TRUE))/sd(tmpAmp, na.rm = TRUE)
    tmpAmpZStat2 <- ((ampSynTable$h_stat2[i] + ampSynTable$m_stat2[i]) - mean(tmpAmpStat2, na.rm = TRUE))/sd(tmpAmpStat2, na.rm = TRUE)
  }
  
  if (delSynTable_bpn$h_freq[i] == 0 | delSynTable_bpn$m_freq[i] == 0) {
    tmpDelZ <- 0
    tmpDelZStat2 <- 0
  } else{
    tmpDelZ <- ((delSynTable_bpn$h_freq[i] + delSynTable_bpn$m_freq[i]) -  mean(tmpDel, na.rm = TRUE))/sd(tmpDel, na.rm = TRUE)
    tmpDelZStat2 <- ((delSynTable$h_stat2[i] + delSynTable$m_stat2[i]) - mean(tmpDelStat2, na.rm = TRUE))/sd(tmpDelStat2, na.rm = TRUE)
  }
  
  permStatResTbl_bpn <- rbind(permStatResTbl_bpn, data.frame("ampZ" = tmpAmpZ, "delZ" = tmpDelZ, "ampStat2Z" = tmpAmpZStat2, "delStat2Z" = tmpDelZStat2))
}


permStatResTbl_bpn$ampP2 <- pnorm(q=permStatResTbl_bpn$ampStat2Z, lower.tail=FALSE)
permStatResTbl_bpn$delP2 <- pnorm(q=permStatResTbl_bpn$delStat2Z, lower.tail=FALSE)

permStatResTbl_bpn_amp <- permStatResTbl_bpn[-which(permStatResTbl_bpn$ampP2 == 0.5), ]
permStatResTbl_bpn_del <- permStatResTbl_bpn[-which(permStatResTbl_bpn$delP2 == 0.5), ]

permStatResTbl_bpn_amp$ampPAdj2 <- p.adjust(permStatResTbl_bpn_amp$ampP2, method = "BH", n = (length(permStatResTbl_bpn_amp$ampP2) + length(permStatResTbl_bpn_del$ampP2)))
permStatResTbl_bpn_del$delPAdj2 <- p.adjust(permStatResTbl_bpn_del$delP2, method = "BH", n = (length(permStatResTbl_bpn_amp$ampP2) + length(permStatResTbl_bpn_del$ampP2)))

tmpSynTable2_bpn_amp <- tmpSynTable_bpn[-which(permStatResTbl_bpn$ampP2 == 0.5), ]
tmpSynTable2_bpn_del <- tmpSynTable_bpn[-which(permStatResTbl_bpn$delP2 == 0.5), ]

hgsc_bpn_arm_allSynTable_amp <- tmpSynTable2_bpn_amp [which(permStatResTbl_bpn_amp$ampPAdj2 < 0.05), ]
hgsc_bpn_arm_allSynTable_del <- tmpSynTable2_bpn_del[which(permStatResTbl_bpn_del$delPAdj2 < 0.05), ]



cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bprn_arm_allSynTable_amp$h_chr, IRanges(hgsc_bprn_arm_allSynTable_amp$h_start, end = hgsc_bprn_arm_allSynTable_amp$h_end))
                                          ,cancerGeneCensusGr)),]

cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bpp_arm_allSynTable_amp$h_chr, IRanges(hgsc_bpp_arm_allSynTable_amp$h_start, end = hgsc_bpp_arm_allSynTable_amp$h_end))
                                          ,cancerGeneCensusGr)),]




cancerGeneCensusAll[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bprn_arm_allSynTable_del$h_chr, IRanges(hgsc_bprn_arm_allSynTable_del$h_start, end = hgsc_bprn_arm_allSynTable_del$h_end))
                                          ,cancerGeneCensusAllGr)),]

cancerGeneCensusAll[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bpp_arm_allSynTable_del$h_chr, IRanges(hgsc_bpp_arm_allSynTable_del$h_start, end = hgsc_bpp_arm_allSynTable_del$h_end))
                                             ,cancerGeneCensusAllGr)),]


cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bpn_arm_allSynTable_amp$h_chr, IRanges(hgsc_bpn_arm_allSynTable_amp$h_start, end = hgsc_bpn_arm_allSynTable_amp$h_end))
                                          ,cancerGeneCensusGr)),]

cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bpp_arm_allSynTable_amp$h_chr, IRanges(hgsc_bpp_arm_allSynTable_amp$h_start, end = hgsc_bpp_arm_allSynTable_amp$h_end))
                                          ,cancerGeneCensusGr)),]




### create simple plots of regions found to be significant, possible hallmark pathway overlaps
### do all this but for coad samples
### lastly take RSEM values and then gain and loss status for arms and look at differences in gene expression
### the pathway overlap can be done by this 
### redo the co-occurence graphs, they're off, there should be double the the rows and columns

hgsc_bprn_arm_allSynTable_simple <- hgsc_bprn_arm_allSynTable[which(hgsc_bprn_arm_allSynTable$str %in% c(hgsc_bprn_arm_allSynTable_amp$str, hgsc_bprn_arm_allSynTable_del$str)),]
hgsc_bprn_arm_allSynTable_simple <- hgsc_bprn_arm_allSynTable_simple[-which(hgsc_bprn_arm_allSynTable_simple$h_freq == 0 | hgsc_bprn_arm_allSynTable_simple$m_freq == 0),]
circosFreqSimple(hgsc_bprn_arm_allSynTable_simple, filename = "20230804bprnPermRes")

hgsc_bpp_arm_allSynTable_simple <- hgsc_bpp_arm_allSynTable[which(hgsc_bpp_arm_allSynTable$str %in% c(hgsc_bpp_arm_allSynTable_amp$str, hgsc_bpp_arm_allSynTable_del$str)),]
hgsc_bpp_arm_allSynTable_simple <- hgsc_bpp_arm_allSynTable_simple[-which(hgsc_bpp_arm_allSynTable_simple$h_freq == 0 | hgsc_bpp_arm_allSynTable_simple$m_freq == 0),]
circosFreqSimple(hgsc_bpp_arm_allSynTable_simple, filename = "20230802bppPermRes")

hgsc_bpn_arm_allSynTable_simple <- hgsc_bpn_arm_allSynTable[which(hgsc_bpn_arm_allSynTable$str %in% c(hgsc_bpn_arm_allSynTable_amp$str, hgsc_bpn_arm_allSynTable_del$str)),]
hgsc_bpn_arm_allSynTable_simple <- hgsc_bpn_arm_allSynTable_simple[-which(hgsc_bpn_arm_allSynTable_simple$h_freq == 0 | hgsc_bpn_arm_allSynTable_simple$m_freq == 0),]
circosFreqSimple(hgsc_bpn_arm_allSynTable_simple, filename = "20230802bpnPermRes")

### mouse chromosome 11 in bprn and not bpp - probably b/c in bpp always had two copy knock-out of brca1 and trp53


hgscAlleleTable <- readxl::read_xlsx("/mnt/DATA6/kevin_recovery//choLab/20190926ScottEvalAnno.xlsx")
hgscAlleleTable$Trp53Geno <- "0/0"
hgscAlleleTable$Rb1Geno <- "0/0"
hgscAlleleTable$Nf1Geno <- "0/0"
hgscAlleleTable$Brca1Geno <- "0/0"

for (i in 1:nrow(hgscAlleleTable)) {
  if(hgscAlleleTable$Brca1.1[i] == "Brca1flox" & hgscAlleleTable$Brca1.2[i] == 0 | hgscAlleleTable$Brca1.2[i] == "Brca1flox" & hgscAlleleTable$Brca1.1[i] == 0){
    hgscAlleleTable$Brca1Geno[i] <- "0/1"
  } else if(hgscAlleleTable$Brca1.1[i] == "Brca1flox" & hgscAlleleTable$Brca1.2[i] == "Brca1flox"){
    hgscAlleleTable$Brca1Geno[i] <- "1/1"
  }
  
  if(hgscAlleleTable$Trp53.1[i] == "Trp53flox" & hgscAlleleTable$Trp53.2[i] == 0 | hgscAlleleTable$Trp53.2[i] == "Trp53flox" & hgscAlleleTable$Trp53.1[i] == 0){
    hgscAlleleTable$Trp53Geno[i] <- "0/1"
  } else if(hgscAlleleTable$Trp53.1[i] == 0 & hgscAlleleTable$Trp53.2[i] == "Trp53LSL-R172H"){
    hgscAlleleTable$Trp53Geno[i] <- "0/1"
  } else if(hgscAlleleTable$Trp53.1[i] == "Trp53flox" & hgscAlleleTable$Trp53.2[i] == "Trp53flox"){
    hgscAlleleTable$Trp53Geno[i] <- "1/1"
  }
  
  if(hgscAlleleTable$Nf1.1[i] == "Nf1flox" & hgscAlleleTable$Nf1.2[i] == 0 | hgscAlleleTable$Nf1.2[i] == "Nf1flox" & hgscAlleleTable$Nf1.1[i] == 0){
    hgscAlleleTable$Nf1Geno[i] <- "0/1"
  } else if(hgscAlleleTable$Nf1.1[i] == "Nf1flox" & hgscAlleleTable$Nf1.2[i] == "Nf1flox"){
    hgscAlleleTable$Nf1Geno[i] <- "1/1"
  }
  
  if(hgscAlleleTable$Rb1.1[i] == "Rb1flox" & hgscAlleleTable$Rb1.2[i] == 0 | hgscAlleleTable$Rb1.2[i] == "Rb1flox" & hgscAlleleTable$Rb1.1[i] == 0){
    hgscAlleleTable$Rb1Geno[i] <- "0/1"
  } else if(hgscAlleleTable$Rb1.1[i] == "Rb1flox" & hgscAlleleTable$Rb1.2[i] == "Rb1flox"){
    hgscAlleleTable$Rb1Geno[i] <- "1/1"
  }
}

### coad, separated by tumor type b/c a lot were from same biological sample, not high enough n without
###
###


coad_adenoCar_arm_allSynTable$str <- paste(coad_adenoCar_arm_allSynTable$h_chr,
                                           coad_adenoCar_arm_allSynTable$h_start,
                                           coad_adenoCar_arm_allSynTable$h_end)
apcCoadArms <- broadBySampleGisticCoad2
tmpSynTable_adenoCar <- coad_adenoCar_arm_allSynTable[-which(duplicated(coad_adenoCar_arm_allSynTable$str)), ]
tmpSynTable_adenoCar$h_chr <- as.numeric(str_remove(tmpSynTable_adenoCar$h_chr, "h_chr"))
tmpSynTable_adenoCar$m_chr <- as.numeric(str_remove(tmpSynTable_adenoCar$m_chr, "m_chr"))
tmpSynTable_adenoCar <- tmpSynTable_adenoCar[order(tmpSynTable_adenoCar$h_chr, tmpSynTable_adenoCar$h_start, decreasing = FALSE),]
apcCoadArmsGr <- GRanges(seqnames = apcCoadArms$chrStripped, IRanges(start = apcCoadArms$start, end = apcCoadArms$end))
coadSyntenyGrange <- GRanges(seqnames = tmpSynTable_adenoCar$h_chr, IRanges(start = tmpSynTable_adenoCar$h_start, end = tmpSynTable_adenoCar$h_end))
coadMm10SyntenyGrange <- GRanges(seqnames = tmpSynTable_adenoCar$m_chr, IRanges(start = tmpSynTable_adenoCar$m_start, end = tmpSynTable_adenoCar$m_end))

humanCoadRegionMat <- data.frame("chr" = tmpSynTable_adenoCar$h_chr, "start" = tmpSynTable_adenoCar$h_start, "end" = tmpSynTable_adenoCar$h_end)
for (i in 4:ncol(apcCoadArms)) {
  
  tmpIdx <- subjectHits(findOverlaps(query = coadSyntenyGrange, subject = apcCoadArmsGr))
  testQuery <- queryHits(findOverlaps(query = coadSyntenyGrange, subject = apcCoadArmsGr))
  
  
  if (length(which(duplicated(testQuery))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated(testQuery))]
    testQuery <- testQuery[-which(duplicated(testQuery))]
  }
  
  ### error probably happens b/c there is no match for one of the queries
  humanCoadRegionMat[ ,i] <- NA
  humanCoadRegionMat[ ,i] <- apcCoadArms[tmpIdx, i]
}


humanCoadRegionMat2 <- humanCoadRegionMat[ , 4:ncol(humanCoadRegionMat)]
cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
synCoadPerm <- foreach(i=1:25,
                       .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                         tmpAmpAll <- NULL
                         tmpDelAll <- NULL
                         tmpAmpFreq <- NULL
                         tmpDelFreq <- NULL
                         j <- 0
                         while (j < 400) {
                           tmp <- apply(humanCoadRegionMat2, 2, function(x) sample(x, nrow(humanCoadRegionMat2)))
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

humanPermAmpCoadAdenoCarStat <- do.call(cbind, synCoadPerm[[1]])
humanPermDelCoadAdenoCarStat <- do.call(cbind, synCoadPerm[[2]])
humanPermAmpCoadAdenoCarFreq <- do.call(cbind, synCoadPerm[[3]])
humanPermDelCoadAdenoCarFreq <- do.call(cbind, synCoadPerm[[4]])

# start <- Sys.time()
# humanCoadRegionMat2 <- humanCoadRegionMat[ , 4:ncol(humanCoadRegionMat)]
# humanPermAmpCoadAdenoCar <- NULL
# humanPermDelCoadAdenoCar <- NULL
# humanPermAmpCoadAdenoCarStat <- NULL
# humanPermDelCoadAdenoCarStat <- NULL
# for (i in 1:10000) {
#   tmp <- apply(humanCoadRegionMat2, 2, function(x) sample(x, 396))
#   tmpAmp <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
#   tmpDel <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
#   tmpAmpStat <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
#   tmpDelStat <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
#   
#   humanPermAmpCoadAdenoCar <- cbind(humanPermAmpCoadAdenoCar, tmpAmp)
#   humanPermDelCoadAdenoCar <- cbind(humanPermDelCoadAdenoCar, tmpDel)
#   humanPermAmpCoadAdenoCarStat <- cbind(humanPermAmpCoadAdenoCarStat,  tmpAmpStat)
#   humanPermDelCoadAdenoCarStat <- cbind(humanPermDelCoadAdenoCarStat, tmpDelStat)
# }
# print( Sys.time() - start )


mouseCoadRegionMat_adenoCar <- data.frame("chr" = tmpSynTable_adenoCar$m_chr, "start" = tmpSynTable_adenoCar$m_start, "end" = tmpSynTable_adenoCar$m_end)
for (i in seq_along(unique(allMouseAneu_coadAdenoCar$sampleID))) {
  
  tmp <- allMouseAneu_coadAdenoCar[which(allMouseAneu_coadAdenoCar$sampleID == unique(allMouseAneu_coadAdenoCar$sampleID)[i]),]
  tmpGrange <- GRanges(seqnames = tmp$chrom, IRanges(start = tmp$start.pos, end = tmp$end.pos))
  
  print(unique(allMouseAneu_coadAdenoCar$sampleID)[i])
  
  tmpIdx <- subjectHits(findOverlaps(query = coadMm10SyntenyGrange, subject = tmpGrange))
  testQuery <- queryHits(findOverlaps(query = coadMm10SyntenyGrange, subject = tmpGrange))
  
  if (length(which(duplicated(testQuery))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated( testQuery))]
    testQuery <- testQuery[-which(duplicated( testQuery))]
  }
  ### error probably happens b/c there is no match for one of the queries
  
  mouseCoadRegionMat_adenoCar[ , i + 3] <- NA
  mouseCoadRegionMat_adenoCar[testQuery, i + 3] <- tmp$mean[tmpIdx]
}

mouseCoadRegionMat_adenoCar2 <- mouseCoadRegionMat_adenoCar[ , 4:ncol(mouseCoadRegionMat_adenoCar)]
cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
synMmCoadPerm <- foreach(i=1:25,
                         .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                           tmpAmpAll <- NULL
                           tmpDelAll <- NULL
                           tmpAmpFreq <- NULL
                           tmpDelFreq <- NULL
                           j <- 0
                           while (j < 400) {
                             tmp <- apply(mouseCoadRegionMat_adenoCar2, 2, function(x) sample(x, nrow(mouseCoadRegionMat_adenoCar2)))
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

# mousePermAmpStat_adenoCar <- do.call(cbind, synMmCoadPerm[[1]]) 
# mousePermDelStat_adenoCar <- do.call(cbind, synMmCoadPerm[[2]])
mousePermAmpFreq_adenoCar <- do.call(cbind, synMmCoadPerm[[3]])
mousePermDelFreq_adenoCar <- do.call(cbind, synMmCoadPerm[[4]])

# start <- Sys.time()
# mouseCoadRegionMat_adenoCar2 <-mouseCoadRegionMat_adenoCar[ , 4:ncol(mouseCoadRegionMat_adenoCar)]
# mousePermAmp_adenoCar <- NULL
# mousePermDel_adenoCar <- NULL
# mousePermAmpStat_adenoCar <- NULL
# mousePermDelStat_adenoCar <- NULL
# for (i in 1:10000) {
#   tmp <- apply(mouseCoadRegionMat_adenoCar, 2, function(x) sample(x, 396))
#   tmpAmp <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
#   tmpDel <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
#   tmpAmpStat2 <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
#   tmpDelStat2 <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
#   
#   mousePermAmp_adenoCar <- cbind(mousePermAmp_adenoCar, tmpAmp)
#   mousePermDel_adenoCar <- cbind(mousePermDel_adenoCar, tmpDel)
#   mousePermAmpStat_adenoCar <- cbind(mousePermAmpStat_adenoCar, tmpAmpStat2)
#   mousePermDelStat_adenoCar <- cbind(mousePermDelStat_adenoCar, tmpDelStat2)
#   
# }
# print( Sys.time() - start )


coad_adenoCar_arm_allSynTable$str <- paste(coad_adenoCar_arm_allSynTable$h_chr, coad_adenoCar_arm_allSynTable$h_start, coad_adenoCar_arm_allSynTable$h_end)
tmpSynTable_adenoCar <- coad_adenoCar_arm_allSynTable[-which(duplicated(coad_adenoCar_arm_allSynTable$str)), ]
tmpSynTable_adenoCar$h_chr <- as.numeric(str_remove(tmpSynTable_adenoCar$h_chr, "h_chr"))
tmpSynTable_adenoCar$m_chr <- as.numeric(str_remove(tmpSynTable_adenoCar$m_chr, "m_chr"))
tmpSynTable_adenoCar <- tmpSynTable_adenoCar[order(tmpSynTable_adenoCar$h_chr, tmpSynTable_adenoCar$h_start, decreasing = FALSE),]

ampSynTable_adenoCar <- data.frame("chr" = tmpSynTable_adenoCar$m_chr, "start" = tmpSynTable_adenoCar$m_start, "end" = tmpSynTable_adenoCar$m_end)
delSynTable_adenoCar <- data.frame("chr" = tmpSynTable_adenoCar$m_chr, "start" = tmpSynTable_adenoCar$m_start, "end" = tmpSynTable_adenoCar$m_end)
ampSynTable_adenoCar$h_freq <- apply(humanCoadRegionMat2, 1, function(x) length(which(x > 0.2))/length(x))
ampSynTable_adenoCar$m_freq <- apply(mouseCoadRegionMat_adenoCar2, 1, function(x) length(which(x > 0.2))/length(x))
delSynTable_adenoCar$h_freq <- apply(humanCoadRegionMat2, 1, function(x) length(which(x < -0.2))/length(x))
delSynTable_adenoCar$m_freq <- apply(mouseCoadRegionMat_adenoCar2, 1, function(x) length(which(x < -0.2))/length(x))
# ampSynTable_adenoCar$h_stat <- apply(humanCoadRegionMat2, 1, function(x) mean(x[which(x > 0.2)]) * length(which(x > 0.2))/length(x))
# ampSynTable_adenoCar$m_stat <- apply(mouseCoadRegionMat_adenoCar2, 1, function(x) mean(x[which(x > 0.2)]) * length(which(x > 0.2))/length(x))
# delSynTable_adenoCar$h_stat <- apply(humanCoadRegionMat2, 1, function(x) abs(mean(x[which(x < -0.2)])) * length(which(x < -0.2))/length(x))
# delSynTable_adenoCar$m_stat <- apply(mouseCoadRegionMat_adenoCar2, 1, function(x) abs(mean(x[which(x < -0.2)])) * length(which(x < -0.2))/length(x))

# ampSynTable_adenoCar$h_freq[which(ampSynTable_adenoCar$h_freq < 0.35)] <- 0
# ampSynTable_adenoCar$m_freq[which(ampSynTable_adenoCar$m_freq < 0.35)] <- 0
# delSynTable_adenoCar$h_freq[which(delSynTable_adenoCar$h_freq < 0.35)] <- 0
# delSynTable_adenoCar$m_freq[which(delSynTable_adenoCar$m_freq < 0.35)] <- 0
# ampSynTable_adenoCar$h_stat[which(ampSynTable_adenoCar$h_freq < 0.35 | ampSynTable_adenoCar$m_freq < 0.35)] <- 0
# ampSynTable_adenoCar$m_stat[which(ampSynTable_adenoCar$h_freq < 0.35 | ampSynTable_adenoCar$m_freq < 0.35)] <- 0
# delSynTable_adenoCar$h_stat[which(delSynTable_adenoCar$h_freq < 0.35 | delSynTable_adenoCar$m_freq < 0.35)] <- 0
# delSynTable_adenoCar$m_stat[which(delSynTable_adenoCar$h_freq < 0.35 | delSynTable_adenoCar$m_freq < 0.35)] <- 0

humanPermAmp2_adenoCar <- humanPermAmpCoadAdenoCarFreq
mousePermAmp2_adenoCar <- mousePermAmpFreq_adenoCar
humanPermDel2_adenoCar <- humanPermDelCoadAdenoCarFreq
mousePermDel2_adenoCar <- mousePermDelFreq_adenoCar
mousePermAmpStat2_adenoCar <- mousePermAmpStat_adenoCar
mousePermDelStat2_adenoCar <- mousePermDelStat_adenoCar
humanPermAmpStatCoad2 <- humanPermAmpCoadAdenoCarStat
humanPermDelStatCoad2 <- humanPermDelCoadAdenoCarStat

# humanPermAmp2_adenoCar[humanPermAmp2_adenoCar < 0.35] <- 0
# mousePermAmp2_adenoCar[mousePermAmp2_adenoCar < 0.35] <- 0
# humanPermDel2_adenoCar[humanPermDel2_adenoCar < 0.35] <- 0
# mousePermDel2_adenoCar[mousePermDel2_adenoCar < 0.35] <- 0
# humanPermAmpStatCoad2[humanPermAmp2_adenoCar < 0.35] <- 0
# humanPermDelStatCoad2[humanPermDel2_adenoCar < 0.35] <- 0
# mousePermAmpStat2_adenoCar[mousePermAmp2_adenoCar < 0.35] <- 0
# mousePermDelStat2_adenoCar[mousePermDel2_adenoCar < 0.35] <- 0

permStatResTbl_adenoCar <- NULL
for (i in 1:396) {
  tmpAmp <- (humanPermAmp2_adenoCar[i, ] + mousePermAmp2_adenoCar[i, ])/2
  tmpDel <- (humanPermDel2_adenoCar[i, ] + mousePermDel2_adenoCar[i, ])/2
  
  tmpAmpStat2 <- humanPermAmpStatCoad2[i, ] + mousePermAmpStat2_adenoCar[i, ]
  tmpDelStat2 <- humanPermDelStatCoad2[i, ] + mousePermDelStat2_adenoCar[i, ]
  
  if (ampSynTable_adenoCar$h_freq[i] == 0 | ampSynTable_adenoCar$m_freq[i] == 0) {
    tmpAmpZ <- 0
    tmpAmpZStat2 <- 0
  } else{
    tmpAmpZ <- ((ampSynTable_adenoCar$h_freq[i] + ampSynTable_adenoCar$m_freq[i]) - mean(tmpAmp, na.rm = TRUE))/sd(tmpAmp, na.rm = TRUE)
    tmpAmpZStat2 <- ((ampSynTable_adenoCar$h_stat[i] + ampSynTable_adenoCar$m_stat[i]) - mean(tmpAmpStat2, na.rm = TRUE))/sd(tmpAmpStat2, na.rm = TRUE)
  }
  
  if (delSynTable_adenoCar$h_freq[i] == 0 | delSynTable_adenoCar$m_freq[i] == 0) {
    tmpDelZ <- 0
    tmpDelZStat2 <- 0
  } else{
    tmpDelZ <- ((delSynTable_adenoCar$h_freq[i] + delSynTable_adenoCar$m_freq[i]) -  mean(tmpDel, na.rm = TRUE))/sd(tmpDel, na.rm = TRUE)
    tmpDelZStat2 <- ((delSynTable_adenoCar$h_stat[i] + delSynTable_adenoCar$m_stat[i]) - mean(tmpDelStat2, na.rm = TRUE))/sd(tmpDelStat2, na.rm = TRUE)
  }
  
  ### alternatively I can look at them mouse and human separately and only nominate ones where they're both significant
  tmpAmpZ_hfreq <- (ampSynTable_adenoCar$h_freq[i] - mean(humanPermAmp2_adenoCar[i, ], na.rm = TRUE))/sd(humanPermAmp2_adenoCar[i, ], na.rm = TRUE)
  tmpDelZ_hfreq <- (delSynTable_adenoCar$h_freq[i] - mean(humanPermDel2_adenoCar[i, ], na.rm = TRUE))/sd(humanPermDel2_adenoCar[i, ], na.rm = TRUE)
  
  tmpAmpZ_mfreq <- (ampSynTable_adenoCar$m_freq[i] - mean(mousePermAmp2_adenoCar[i, ], na.rm = TRUE))/sd(mousePermAmp2_adenoCar[i, ], na.rm = TRUE)
  tmpDelZ_mfreq <- (delSynTable_adenoCar$m_freq[i] - mean(mousePermDel2_adenoCar[i, ], na.rm = TRUE))/sd(mousePermDel2_adenoCar[i, ], na.rm = TRUE)
  
  
  
  # permStatResTbl_adenoCar <- rbind(permStatResTbl_adenoCar,
  #                                  data.frame("ampZ" = tmpAmpZ, "delZ" = tmpDelZ))
  
  permStatResTbl_adenoCar <- rbind(permStatResTbl_adenoCar,
                                   data.frame("ampZ_hfreq" = tmpAmpZ_hfreq, "delZ_hfreq" = tmpDelZ_hfreq,
                                              "ampZ_mfreq" = tmpAmpZ_mfreq, "delZ_mfreq" = tmpDelZ_mfreq))
}

permStatResTbl_adenoCar$ampZ_hfreq2 <- pnorm(q=permStatResTbl_adenoCar$ampZ_hfreq, lower.tail=FALSE)
permStatResTbl_adenoCar$delZ_hfreq2 <- pnorm(q=permStatResTbl_adenoCar$delZ_hfreq, lower.tail=FALSE)
permStatResTbl_adenoCar$ampZ_mfreq2 <- pnorm(q=permStatResTbl_adenoCar$ampZ_mfreq, lower.tail=FALSE)
permStatResTbl_adenoCar$delZ_mfreq2 <- pnorm(q=permStatResTbl_adenoCar$delZ_mfreq, lower.tail=FALSE)


# permStatResTbl_adenoCar$ampP2 <- pnorm(q=permStatResTbl_adenoCar$ampStat2Z, lower.tail=FALSE)
# permStatResTbl_adenoCar$delP2 <- pnorm(q=permStatResTbl_adenoCar$delStat2Z, lower.tail=FALSE)
# 
# permStatResTbl_adenoCar$ampP <- pnorm(q=permStatResTbl_adenoCar$ampZ, lower.tail=FALSE)
# permStatResTbl_adenoCar$delP <- pnorm(q=permStatResTbl_adenoCar$delZ, lower.tail=FALSE)

# adenoCar_arm_allSynTable_amp$ampP2 <- permStatResTbl_adenoCar$ampP2
# adenoCar_arm_allSynTable_amp$ampP2 <- permStatResTbl_adenoCar$ampP2
# adenoCar_arm_allSynTable_amp$adjP2 <- p.adjust(adenoCar_arm_allSynTable_amp$ampP2, method = "BH", n = 2 * length(adenoCar_arm_allSynTable_amp$ampP2))
# adenoCar_arm_allSynTable_amp$ampP <- permStatResTbl_adenoCar$ampP
# adenoCar_arm_allSynTable_amp$adjP <- p.adjust(adenoCar_arm_allSynTable_amp$ampP, method = "BH", n = 2 * length(adenoCar_arm_allSynTable_amp$ampP))
# adenoCar_arm_allSynTable_amp$freqStat <- (adenoCar_arm_allSynTable_amp$h_freq + adenoCar_arm_allSynTable_amp$m_freq)/2



adenoCar_arm_allSynTable_amp <- tmpSynTable_adenoCar
adenoCar_arm_allSynTable_amp$h_freq <- ampSynTable_adenoCar$h_freq
adenoCar_arm_allSynTable_amp$m_freq <- ampSynTable_adenoCar$m_freq
adenoCar_arm_allSynTable_amp$hfreqP <- permStatResTbl_adenoCar$ampZ_hfreq2
adenoCar_arm_allSynTable_amp$hfreqAdjP <- p.adjust(adenoCar_arm_allSynTable_amp$hfreqP, method = "BH",
                                                   n = 4 * length(adenoCar_arm_allSynTable_amp$hfreqP))
adenoCar_arm_allSynTable_amp$mfreqP <- permStatResTbl_adenoCar$ampZ_mfreq2
adenoCar_arm_allSynTable_amp$mfreqAdjP <- p.adjust(adenoCar_arm_allSynTable_amp$mfreqP, method = "BH",
                                                   n = 4 * length(adenoCar_arm_allSynTable_amp$mfreqP))
adenoCar_arm_allSynTable_amp$sigCheck <- ifelse(adenoCar_arm_allSynTable_amp$hfreqAdjP < 0.05 & adenoCar_arm_allSynTable_amp$mfreqAdjP < 0.05,
                                                "good", "bad")




# adenoCar_arm_allSynTable_del <- tmpSynTable_adenoCar
# adenoCar_arm_allSynTable_del$h_freq <- delSynTable_adenoCar$h_freq
# adenoCar_arm_allSynTable_del$m_freq <- delSynTable_adenoCar$m_freq
# adenoCar_arm_allSynTable_del$delP2 <- permStatResTbl_adenoCar$delP2
# adenoCar_arm_allSynTable_del$adjP2 <- p.adjust(adenoCar_arm_allSynTable_del$delP2, method = "BH", n = 2 * length(adenoCar_arm_allSynTable_del$delP2))


adenoCar_arm_allSynTable_del <- tmpSynTable_adenoCar
adenoCar_arm_allSynTable_del$h_freq <- delSynTable_adenoCar$h_freq
adenoCar_arm_allSynTable_del$m_freq <- delSynTable_adenoCar$m_freq
adenoCar_arm_allSynTable_del$hfreqP <- permStatResTbl_adenoCar$delZ_hfreq2
adenoCar_arm_allSynTable_del$hfreqAdjP <- p.adjust(adenoCar_arm_allSynTable_del$hfreqP, method = "BH",
                                                   n = 4 * length(adenoCar_arm_allSynTable_del$hfreqP))
adenoCar_arm_allSynTable_del$mfreqP <- permStatResTbl_adenoCar$delZ_mfreq2
adenoCar_arm_allSynTable_del$mfreqAdjP <- p.adjust(adenoCar_arm_allSynTable_del$mfreqP, method = "BH",
                                                   n = 4 * length(adenoCar_arm_allSynTable_del$mfreqP))
adenoCar_arm_allSynTable_del$sigCheck <- ifelse(adenoCar_arm_allSynTable_del$hfreqAdjP < 0.05 & adenoCar_arm_allSynTable_del$mfreqAdjP < 0.05,
                                                "good", "bad")



cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = adenoCar_arm_allSynTable_amp$h_chr, IRanges(adenoCar_arm_allSynTable_amp$h_start, end = adenoCar_arm_allSynTable_amp$h_end))
                                          ,cancerGeneCensusGr)),]

cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = adenoCar_arm_allSynTable_del$h_chr, IRanges(adenoCar_arm_allSynTable_del$h_start, end = adenoCar_arm_allSynTable_del$h_end))
                                          ,cancerGeneCensusGr)),]

### adenoma


coad_adeno_arm_allSynTable$str <- paste(coad_adeno_arm_allSynTable$h_chr, coad_adeno_arm_allSynTable$h_start, coad_adeno_arm_allSynTable$h_end)
tmpSynTable_Adeno <- coad_adeno_arm_allSynTable[-which(duplicated(coad_adeno_arm_allSynTable$str)), ]
tmpSynTable_Adeno$h_chr <- as.numeric(str_remove(tmpSynTable_Adeno$h_chr, "h_chr"))
tmpSynTable_Adeno$m_chr <- as.numeric(str_remove(tmpSynTable_Adeno$m_chr, "m_chr"))
tmpSynTable_Adeno <- tmpSynTable_Adeno[order(tmpSynTable_Adeno$h_chr, tmpSynTable_Adeno$h_start, decreasing = FALSE),]


mouseCoadRegionMat_Adeno <- data.frame("chr" = tmpSynTable_Adeno$m_chr, "start" = tmpSynTable_Adeno$m_start, "end" = tmpSynTable_Adeno$m_end)
for (i in seq_along(unique(allMouseAneu_coadAdeno$sampleID))) {
  
  tmp <- allMouseAneu_coadAdeno[which(allMouseAneu_coadAdeno$sampleID == unique(allMouseAneu_coadAdeno$sampleID)[i]),]
  tmpGrange <- GRanges(seqnames = tmp$chrom, IRanges(start = tmp$start.pos, end = tmp$end.pos))
  
  print(unique(allMouseAneu_coadAdeno$sampleID)[i])
  
  tmpIdx <- subjectHits(findOverlaps(query = coadMm10SyntenyGrange, subject = tmpGrange))
  testQuery <- queryHits(findOverlaps(query = coadMm10SyntenyGrange, subject = tmpGrange))
  
  if (length(which(duplicated(testQuery))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated( testQuery))]
    testQuery <- testQuery[-which(duplicated( testQuery))]
  }
  ### error probably happens b/c there is no match for one of the queries
  
  mouseCoadRegionMat_Adeno[ , i + 3] <- NA
  mouseCoadRegionMat_Adeno[testQuery, i + 3] <- tmp$mean[tmpIdx]
}


start <- Sys.time()
mouseCoadRegionMat_Adeno2 <-mouseCoadRegionMat_Adeno[ , 4:ncol(mouseCoadRegionMat_Adeno)]
mousePermAmp_Adeno <- NULL
mousePermDel_Adeno <- NULL
mousePermAmpStat_Adeno <- NULL
mousePermDelStat_Adeno <- NULL
for (i in 1:10000) {
  tmp <- apply(mouseCoadRegionMat_Adeno, 2, function(x) sample(x, 396))
  tmpAmp <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
  tmpDel <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
  tmpAmpStat2 <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
  tmpDelStat2 <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
  
  mousePermAmp_Adeno <- cbind(mousePermAmp_Adeno, tmpAmp)
  mousePermDel_Adeno <- cbind(mousePermDel_Adeno, tmpDel)
  mousePermAmpStat_Adeno <- cbind(mousePermAmpStat_Adeno, tmpAmpStat2)
  mousePermDelStat_Adeno <- cbind(mousePermDelStat_Adeno, tmpDelStat2)
  
}
print( Sys.time() - start )


ampSynTable_Adeno <- data.frame("chr" = tmpSynTable_Adeno$m_chr, "start" = tmpSynTable_Adeno$m_start, "end" = tmpSynTable_Adeno$m_end)
delSynTable_Adeno <- data.frame("chr" = tmpSynTable_Adeno$m_chr, "start" = tmpSynTable_Adeno$m_start, "end" = tmpSynTable_Adeno$m_end)
ampSynTable_Adeno$h_freq <- apply(humanCoadRegionMat2, 1, function(x) length(which(x > 0))/length(x))
ampSynTable_Adeno$m_freq <- apply(mouseCoadRegionMat_Adeno2, 1, function(x) length(which(x > 0))/length(x))
delSynTable_Adeno$h_freq <- apply(humanCoadRegionMat2, 1, function(x) length(which(x < 0))/length(x))
delSynTable_Adeno$m_freq <- apply(mouseCoadRegionMat_Adeno2, 1, function(x) length(which(x < 0))/length(x))
ampSynTable_Adeno$h_stat <- apply(humanCoadRegionMat2, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
ampSynTable_Adeno$m_stat <- apply(mouseCoadRegionMat_Adeno2, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
delSynTable_Adeno$h_stat <- apply(humanCoadRegionMat2, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
delSynTable_Adeno$m_stat <- apply(mouseCoadRegionMat_Adeno2, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))

ampSynTable_Adeno$h_freq[which(ampSynTable_Adeno$h_freq < 0.1)] <- 0
ampSynTable_Adeno$m_freq[which(ampSynTable_Adeno$m_freq < 0.1)] <- 0
delSynTable_Adeno$h_freq[which(delSynTable_Adeno$h_freq < 0.1)] <- 0
delSynTable_Adeno$m_freq[which(delSynTable_Adeno$m_freq < 0.1)] <- 0
ampSynTable_Adeno$h_stat[which(ampSynTable_Adeno$h_freq < 0.1)] <- 0
ampSynTable_Adeno$m_stat[which(ampSynTable_Adeno$m_freq < 0.1)] <- 0
delSynTable_Adeno$h_stat[which(delSynTable_Adeno$h_freq < 0.1)] <- 0
delSynTable_Adeno$m_stat[which(delSynTable_Adeno$m_freq < 0.1)] <- 0

humanPermAmp2_Adeno <- humanPermAmp
mousePermAmp2_Adeno <- mousePermAmp_Adeno
humanPermDel2_Adeno <- humanPermDel
mousePermDel2_Adeno <- mousePermDel_Adeno
mousePermAmpStat2_Adeno <- mousePermAmpStat_Adeno
mousePermDelStat2_Adeno <- mousePermDelStat_Adeno

humanPermAmp2_Adeno[humanPermAmp2_Adeno < 0.1] <- 0
mousePermAmp2_Adeno[mousePermAmp2_Adeno < 0.1] <- 0
humanPermDel2_Adeno[humanPermDel2_Adeno < 0.1] <- 0
mousePermDel2_Adeno[mousePermDel2_Adeno < 0.1] <- 0
humanPermAmpStat2v2[humanPermAmp2 < 0.1] <- 0
humanPermDelStat2v2[humanPermDel2 < 0.1] <- 0
mousePermAmpStat2_Adeno[mousePermAmp2_Adeno < 0.1] <- 0
mousePermDelStat2_Adeno[mousePermDel2_Adeno < 0.1] <- 0


permStatResTbl_Adeno <- NULL
for (i in 1:396) {
  tmpAmp <- humanPermAmp2_Adeno[i, ] + mousePermAmp2_Adeno[i, ]
  tmpDel <- humanPermDel2_Adeno[i, ] + mousePermDel2_Adeno[i, ]
  
  tmpAmpStat2 <- humanPermAmpStat2v2[i, ] + mousePermAmpStat2_Adeno[i, ]
  tmpDelStat2 <- humanPermDelStat2v2[i, ] + mousePermDelStat2_Adeno[i, ]
  
  if (ampSynTable_Adeno$h_freq[i] == 0 | ampSynTable_Adeno$m_freq[i] == 0) {
    tmpAmpZ <- 0
    tmpAmpZStat2 <- 0
  } else{
    tmpAmpZ <- ((ampSynTable_Adeno$h_freq[i] + ampSynTable_Adeno$m_freq[i]) - mean(tmpAmp, na.rm = TRUE))/sd(tmpAmp, na.rm = TRUE)
    tmpAmpZStat2 <- ((ampSynTable_Adeno$h_stat[i] + ampSynTable_Adeno$m_stat[i]) - mean(tmpAmpStat2, na.rm = TRUE))/sd(tmpAmpStat2, na.rm = TRUE)
  }
  
  if (delSynTable_Adeno$h_freq[i] == 0 | delSynTable_Adeno$m_freq[i] == 0) {
    tmpDelZ <- 0
    tmpDelZStat2 <- 0
  } else{
    tmpDelZ <- ((delSynTable_Adeno$h_freq[i] + delSynTable_Adeno$m_freq[i]) -  mean(tmpDel, na.rm = TRUE))/sd(tmpDel, na.rm = TRUE)
    tmpDelZStat2 <- ((delSynTable_Adeno$h_stat[i] + delSynTable_Adeno$m_stat[i]) - mean(tmpDelStat2, na.rm = TRUE))/sd(tmpDelStat2, na.rm = TRUE)
  }
  
  permStatResTbl_Adeno <- rbind(permStatResTbl_Adeno, data.frame("ampZ" = tmpAmpZ, "delZ" = tmpDelZ, "ampStat2Z" = tmpAmpZStat2, "delStat2Z" = tmpDelZStat2))
}


permStatResTbl_Adeno$ampP2 <- pnorm(q=permStatResTbl_Adeno$ampStat2Z, lower.tail=FALSE)
permStatResTbl_Adeno$delP2 <- pnorm(q=permStatResTbl_Adeno$delStat2Z, lower.tail=FALSE)

permStatResTbl_Adeno_amp <- permStatResTbl_Adeno
permStatResTbl_Adeno_del <- permStatResTbl_Adeno

permStatResTbl_Adeno_amp$ampPAdj2 <- p.adjust(permStatResTbl_Adeno_amp$ampP2, method = "BH", n = (length(permStatResTbl_Adeno_amp$ampP2) + length(permStatResTbl_Adeno_del$ampP2)))
permStatResTbl_Adeno_del$delPAdj2 <- p.adjust(permStatResTbl_Adeno_del$delP2, method = "BH", n = (length(permStatResTbl_Adeno_amp$ampP2) + length(permStatResTbl_Adeno_del$ampP2)))

tmpSynTable2_Adeno_amp <- tmpSynTable_Adeno
tmpSynTable2_Adeno_del <- tmpSynTable_Adeno

coad_adeno_arm_allSynTable_amp <- tmpSynTable2_Adeno_amp[which(permStatResTbl_Adeno_amp$ampPAdj2 < 0.05), ]
coad_adeno_arm_allSynTable_del <- tmpSynTable2_Adeno_del[which(permStatResTbl_Adeno_del$delPAdj2 < 0.05), ]


coad_adenoCar_arm_allSynTable_simple <- coad_adenoCar_arm_allSynTable[which(coad_adenoCar_arm_allSynTable$str %in% c(adenoCar_arm_allSynTable_amp$str, adenoCar_arm_allSynTable_del$str)),]
coad_adenoCar_arm_allSynTable_simple <- coad_adenoCar_arm_allSynTable_simple[-which(coad_adenoCar_arm_allSynTable_simple$h_freq == 0 | coad_adenoCar_arm_allSynTable_simple$m_freq == 0),]
circosFreqSimple(coad_adenoCar_arm_allSynTable_simple, filename = "20230804coadAdenoCarPermRes")

coad_adeno_arm_allSynTable_simple <- coad_adeno_arm_allSynTable[which(coad_adeno_arm_allSynTable$str %in% c(coad_adeno_arm_allSynTable_amp$str, coad_adeno_arm_allSynTable_del$str)),]
coad_adeno_arm_allSynTable_simple <- coad_adeno_arm_allSynTable_simple[-which(coad_adeno_arm_allSynTable_simple$h_freq == 0 | coad_adeno_arm_allSynTable_simple$m_freq == 0),]
circosFreqSimple(coad_adeno_arm_allSynTable_simple, filename = "20230804coadAdenoPermRes")

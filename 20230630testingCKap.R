allCnaStatsDf2 <- NULL
allCnaStatsDf2_ck <- NULL
for (i in c(1, 5, 10, 20)) {
  
  for (j in c(0:4)) {
    ### read in tables after they've gone through bedtools
    tmpDisJoin <- read.table(paste0("/mnt/DATA5/tmp/kev/misc/20230606ngsInter_", j, "_", i, "Mb.bed"),
                             sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    
    colnames(tmpDisJoin) <- c("chrom", "start", "end","sample", "seg.mean", "str", "size")
    
    tmpDisJoinAneuploidy <- matrix(0, nrow = length(unique(sampleNames2)), ncol = length(unique(segResFilt3$chrom)))
    rownames(tmpDisJoinAneuploidy) <- unique(sampleNames2)
    colnames(tmpDisJoinAneuploidy) <- unique(segResFilt3$chrom)
    
    tmp_cna <- NULL
    k <- rownames(tmpDisJoinAneuploidy)[1]
    for (k in rownames(tmpDisJoinAneuploidy)) {
      tmpDf <- tmpDisJoin[grep(k, tmpDisJoin$sample),]
      
      for (l in seq_along(tmpChrList)) {
        tmpChr <- tmpDf[which(tmpDf$chrom == tmpChrList[l]), ]
        tmpChr2 <- tmpChr
        chrThres <- sum(tmpChr$size) * 0.7
        tmpAneu <- "none"
        if (sum(tmpChr$size[which(tmpChr$seg.mean > 0.2)]) > chrThres) {
          tmpAneu <- "gain"
          tmpChr2$seg.mean[which(tmpChr2$seg.mean > 0.2)] <- 0
        } else if(sum(tmpChr$size[which(tmpChr$seg.mean < -0.2)]) > chrThres){
          tmpAneu <- "loss"
          tmpChr2$seg.mean[which(tmpChr2$seg.mean < -0.2)] <- 0
        } else{
          tmpAneu <- "none"
        }
        tmpDisJoinAneuploidy[k, l] <- tmpAneu
        tmp_cna <- rbind(tmp_cna, tmpChr2)
      }
    }
    
    ### commented out to compare good samples to all
    
    tmpGood <- tmpDisJoinAneuploidy[which(rownames(tmpDisJoinAneuploidy) %in% goodSamps), 1:19]
    tmpGoodVec <- as.vector(tmpGood)
    tmpGoodVec <- factor(tmpGoodVec, levels = c("none", "gain", "loss"))
    # 
    goodAscatVec2n <- factor(ascatAneploidy2n[which(rownames(ascatAneploidy2n) %in% goodSamps), ], levels = c("none", "gain", "loss"))
    goodAscatVec3n <- factor(ascatAneploidy3n[which(rownames(ascatAneploidy3n) %in% goodSamps), ], levels = c("none", "gain", "loss"))
    goodAscatVec4n <- factor(ascatAneploidy4n[which(rownames(ascatAneploidy4n) %in% goodSamps), ], levels = c("none", "gain", "loss"))
    
    goodConMat2n_pr <- confusionMatrix(tmpGoodVec, reference = goodAscatVec2n, mode = "prec_recall")
    goodConMat3n_pr <- confusionMatrix(tmpGoodVec, reference = goodAscatVec3n, mode = "prec_recall")
    goodConMat4n_pr <- confusionMatrix(tmpGoodVec, reference = goodAscatVec4n, mode = "prec_recall")
    
    goodConMat2n_ck <- irr::kappa2(cbind(tmpGoodVec, goodAscatVec2n), weight = "squared")
    goodConMat3n_ck <- irr::kappa2(cbind(tmpGoodVec, goodAscatVec3n), weight = "squared")
    goodConMat4n_ck <- irr::kappa2(cbind(tmpGoodVec, goodAscatVec4n), weight = "squared")
    
    ### comment out for all samples
    tmp_cna <-  tmp_cna[which(tmp_cna$sample %in% goodSamps), ]
    
    if (length(which(is.na(tmp_cna$seg.mean))) > 0) {
      tmp_cna <- tmp_cna[-which(is.na(tmp_cna$seg.mean)),]
    }
    
    tmp_cna$seg.mean[which(abs(tmp_cna$seg.mean) < 0.2)] <- 0
    
    tmp_cna2 <- NULL
    k <- unique(tmp_cna$sample)[1]
    for (k in unique(tmp_cna$sample)) {
      tmpDf <- tmp_cna[which(tmp_cna$sample == k), ]
      l <- unique(tmpDf$chrom)[16]
      for (l in unique(tmpDf$chrom)) {
        tmpDfChr <- tmpDf[which(tmpDf$chrom == l),]
        if (any(tmpDfChr$seg.mean != 0)) {
          
          tmpDfChr2 <- tmpDfChr[-which(tmpDfChr$seg.mean == 0), ]
          tmpDfChr <- tmpDfChr[which(tmpDfChr$seg.mean == 0), ]
          
          for (m in unique(tmpDfChr2$seg.mean)) {
            tmpDfChr3 <- tmpDfChr2[which(tmpDfChr2$seg.mean == m),]
            tmpChrSeg <- tmpDfChr3[1, 1:6]
            tmpChrSeg$start <- min(tmpDfChr3$start)
            tmpChrSeg$end <- max(tmpDfChr3$end)
            tmpChrSeg$size <- tmpChrSeg$end - tmpChrSeg$start
            tmpDfChr <- rbind(tmpDfChr, tmpChrSeg)
          }
          tmpDfChr <- tmpDfChr[order(tmpDfChr$start), ]
          tmp_cna2 <- rbind(tmp_cna2, tmpDfChr)
        } else{
          tmp_cna2 <- rbind(tmp_cna2, tmpDfChr)
        }
      }
    }
    
    if (length(which(tmp_cna2$seg.mean != 0)) == 0) {
      next() 
    }
    
    
    ### ascat portion - need tmpGrange 1 to change if I use different ploidy
    for (q in 2:4) {
      if(q == 2){
        tmpGrange <- tmpGrange1
        tmpDfPloidy <- tmpDf2n
      } else if(q == 3){
        tmpGrange <- tmpGrange2
        tmpDfPloidy <- tmpDf3n
      } else if(q == 4){
        tmpGrange <- tmpGrange3
        tmpDfPloidy <- tmpDf4n
      }
      tmpAscatCna <- NULL
      k <- unique(tmpGrange$sample)[1]
      for (k in unique(tmpGrange$sample)) {
        tmpN <- tmpDfPloidy[which(tmpDfPloidy$sample == k),]
        tmpRes <- tmp_cna2[which(tmp_cna2$seg.mean != 0) ,1:4]
        
        if (grepl("133576rt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "13576rt"), ]
        } else if (grepl("14150rt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "14150lt"), ]
        } else if (grepl("14154rot", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "14154lt"), ]
        } else if (grepl("133576rt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "13576rt"), ]
        } else if (grepl("14656peritonealmt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "14656peritnealmt"), ] 
        } else{
          tmpRes <- tmpRes[which(tmpRes$sample == allSampsCorDf2_bestV2$sampleStripped2[which(allSampsCorDf2_bestV2$string5 == k)]), ]
        }
        
        if (nrow(tmpRes) == 0) {
          next()
        }
        
        queryIdx <- 1:nrow(tmpRes)
        tmpGr <- GRanges(seqnames = tmpN$chr, IRanges(start = tmpN$startpos, end = tmpN$endpos))
        tmpGr2 <- GRanges(seqnames = tmpRes$chrom, IRanges(start = tmpRes$start, end = tmpRes$end))
        tmpIdx <- subjectHits(findOverlaps(query = tmpGr2, subject = tmpGr))
        testQuery <- queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))
        
        if (length(which(duplicated(queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))))) > 0) {
          tmpIdx <-  tmpIdx[-which(duplicated(queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))))]
          testQuery <- testQuery[-which(duplicated(queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))))]
        }
        ### error probably happens b/c there is no match for one of the queries
        tmpRes$cn <- NA
        tmpRes$cn[testQuery] <- tmpN$totalSc[tmpIdx]
        tmpRes2 <- tmpRes[, c("sample", "chrom", "start", "end", "cn")]
        tmpAscatCna <- rbind(tmpAscatCna, tmpRes2)
      }
      
      tmpAscatCna$absCn <- tmpAscatCna$cn - q
      tmpAscatCna$type <- "none"
      tmpAscatCna$type[which(tmpAscatCna$absCn > 0)] <- "gain"
      tmpAscatCna$type[which(tmpAscatCna$absCn < 0)] <- "loss"
      tmpAscatCna$size <- tmpAscatCna$end - tmpAscatCna$start
      
      ### for only good samples
      tmpAscatCna <- tmpAscatCna[which(tmpAscatCna$sample %in% goodSamps), ]
      
      
      tmp_cna3 <- tmp_cna2[which(tmp_cna2$seg.mean != 0),]
      if (nrow(tmp_cna3) == 0) {
        next()
      }
      tmp_cna3 <- tmp_cna3[order(match(tmp_cna3$sample, unique(tmpAscatCna$sample))),]
      tmp_cna3$type <- "none"
      tmp_cna3$type[which(tmp_cna3$seg.mean >  0.2)] <- "gain"
      tmp_cna3$type[which(tmp_cna3$seg.mean <  -0.2)] <- "loss"
      
      tmpLess1Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size < 1e6)], levels = c("none", "gain", "loss")),
                                    reference = factor(tmpAscatCna$type[which(tmpAscatCna$size < 1e6)], levels = c("none", "gain", "loss")),
                                    mode = "prec_recall")
      
      tmpGreater1Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 1e6)], levels = c("none", "gain", "loss")),
                                       reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 1e6)], levels = c("none", "gain", "loss")),
                                       mode = "prec_recall")
      
      tmpGreater5Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 5e6)], levels = c("none", "gain", "loss")),
                                       reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 5e6)], levels = c("none", "gain", "loss")),
                                       mode = "prec_recall")
      
      tmpGreater10Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 10e6)], levels = c("none", "gain", "loss")),
                                        reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 10e6)], levels = c("none", "gain", "loss")),
                                        mode = "prec_recall")
      
      tmpGreater20Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 20e6)], levels = c("none", "gain", "loss")),
                                        reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 20e6)], levels = c("none", "gain", "loss")),
                                        mode = "prec_recall")
      
      ### cohen's kappa - treating these as nominal - 
      
      # tmpLess1Mb_ck <- tryCatch(psych::cohen.kappa(cbind(tmpAscatCna$type[which(tmpAscatCna$size < 1e6)],
      #                                                    tmp_cna3$type[which(tmpAscatCna$size < 1e6)])),
      #                           error = function(w) print(NA))
      # 
      # tmpGreater1Mb_ck <- tryCatch(psych::cohen.kappa(cbind(tmpAscatCna$type[which(tmpAscatCna$size > 1e6)],
      #                                                       tmp_cna3$type[which(tmpAscatCna$size > 1e6)])),
      #                              error = function(w) print(NA))
      # 
      # tmpGreater5Mb_ck <- tryCatch(psych::cohen.kappa(cbind(tmpAscatCna$type[which(tmpAscatCna$size > 5e6)],
      #                                                       tmp_cna3$type[which(tmpAscatCna$size > 5e6)])),
      #                              error = function(w) print(NA))
      # 
      # tmpGreater10Mb_ck <-tryCatch(psych::cohen.kappa(cbind(tmpAscatCna$type[which(tmpAscatCna$size > 1e7)],
      #                                                       tmp_cna3$type[which(tmpAscatCna$size > 1e7)])),
      #                              error = function(w) print(NA))
      # 
      # tmpGreater20Mb_ck <- tryCatch(psych::cohen.kappa(cbind(tmpAscatCna$type[which(tmpAscatCna$size > 2e7)],
      #                                                        tmp_cna3$type[which(tmpAscatCna$size > 2e7)])),
      #                               error = function(w) print(NA))
      
      ### switch to quadratic for difference in importance for getting none change to gain or loss vs gain to loss
      
      tmpLess1Mb_ck <- tryCatch(irr::kappa2(cbind(factor(tmp_cna3$type[which(tmpAscatCna$size < 1e6)], levels = c( "gain", "none","loss")),
                                                  factor(tmpAscatCna$type[which(tmpAscatCna$size < 1e6)], levels = c("gain", "none", "loss"))),
                                            weight = "squared"), error = function(w) print(NA))
      
      tmpGreater1Mb_ck <- tryCatch(irr::kappa2(cbind(factor(tmp_cna3$type[which(tmpAscatCna$size > 1e6)], levels = c( "gain", "none","loss")),
                                                     factor(tmpAscatCna$type[which(tmpAscatCna$size > 1e6)], levels = c("gain", "none", "loss"))),
                                               weight = "squared"), error = function(w) print(NA))
      
      tmpGreater5Mb_ck <- tryCatch(irr::kappa2(cbind(factor(tmp_cna3$type[which(tmpAscatCna$size > 5e6)], levels = c( "gain", "none","loss")),
                                                     factor(tmpAscatCna$type[which(tmpAscatCna$size > 5e6)], levels = c("gain", "none", "loss"))),
                                               weight = "squared"), error = function(w) print(NA))
      
      tmpGreater10Mb_ck <- tryCatch(irr::kappa2(cbind(factor(tmp_cna3$type[which(tmpAscatCna$size > 1e7)], levels = c( "gain", "none","loss")),
                                                      factor(tmpAscatCna$type[which(tmpAscatCna$size > 1e7)], levels = c("gain", "none", "loss"))),
                                                weight = "squared"), error = function(w) print(NA))
      
      tmpGreater20Mb_ck <- tryCatch(irr::kappa2(cbind(factor(tmp_cna3$type[which(tmpAscatCna$size > 2e7)], levels = c( "gain", "none","loss")),
                                                      factor(tmpAscatCna$type[which(tmpAscatCna$size > 2e7)], levels = c("gain", "none", "loss"))),
                                                weight = "squared"), error = function(w) print(NA))
      
      
      if (is.na(tmpGreater1Mb_ck)) {
        tmpGreater1Mb_ck <- tmpLess1Mb_ck
        tmpGreater1Mb_ck$value <- NA
      }
      if (is.na(tmpGreater5Mb_ck)) {
        tmpGreater5Mb_ck <- tmpLess1Mb_ck
        tmpGreater5Mb_ck$value <- NA
      }
      if (is.na(tmpGreater10Mb_ck)) {
        tmpGreater10Mb_ck <- tmpLess1Mb_ck
        tmpGreater10Mb_ck$value <- NA
      }
      if (is.na(tmpGreater20Mb_ck)) {
        tmpGreater20Mb_ck <- tmpLess1Mb_ck
        tmpGreater20Mb_ck$value <- NA
      }
      
      
      tmpCnaStatsDf <- NULL
      tmpCnaStatsDf_ck <- NULL
      
      tmpCnaStatsDf <- rbind(tmpCnaStatsDf,
                             data.frame("type" = rep("aneuploidy", 3), "change"  = rownames(eval(parse(text = paste0("goodConMat", q,"n", "_pr")))$byClass),
                                        eval(parse(text = paste0("goodConMat", q, "n", "_pr")))$byClass),
                             data.frame("type" = rep("<1Mb", 3), "change"  = rownames(tmpLess1Mb$byClass), tmpLess1Mb$byClass),
                             data.frame("type" = rep(">1Mb", 3), "change"  = rownames(tmpGreater1Mb$byClass), tmpGreater1Mb$byClass),
                             data.frame("type" = rep(">5Mb", 3), "change"  = rownames(tmpGreater5Mb$byClass), tmpGreater5Mb$byClass),
                             data.frame("type" = rep(">10Mb", 3), "change"  = rownames(tmpGreater10Mb$byClass), tmpGreater10Mb$byClass),
                             data.frame("type" = rep(">20Mb", 3), "change"  = rownames(tmpGreater20Mb$byClass), tmpGreater20Mb$byClass))
      
      tmpCnaStatsDf_ck <- rbind(data.frame("type" = "aneuploidy", "kappa" = eval(parse(text = paste0("goodConMat", q,"n", "_ck")))$value),
                                data.frame("type" = "<1Mb", "kappa" = tmpLess1Mb_ck$value),
                                data.frame("type" = ">1Mb", "kappa" = tmpGreater1Mb_ck$value),
                                data.frame("type" = ">5Mb", "kappa" = tmpGreater5Mb_ck$value),
                                data.frame("type" = ">10Mb", "kappa" = tmpGreater10Mb_ck$value),
                                data.frame("type" = ">20Mb", "kappa" = tmpGreater20Mb_ck$value))
      
      tmpCnaStatsDf$windowSize <- paste0(i, "Mb")
      tmpCnaStatsDf$ploidy <- paste0(q, "n")
      tmpCnaStatsDf$markerFilt <- j
      tmpCnaStatsDf_ck$windowSize <- paste0(i, "Mb")
      tmpCnaStatsDf_ck$ploidy <- paste0(q, "n")
      tmpCnaStatsDf_ck$markerFilt <- j
      
      allCnaStatsDf2 <- rbind(allCnaStatsDf2, tmpCnaStatsDf)
      allCnaStatsDf2_ck <- rbind(allCnaStatsDf2_ck, tmpCnaStatsDf_ck)
      
    }
  }
}

allCnaStatsDf2_goodSamps <- allCnaStatsDf2
allCnaStatsDf2ck_goodSamps <- allCnaStatsDf2_ck


### playing with cohen's kappa
### plus looking differnces in weighted + unweighted and looking at the differrences

rater1 <- c(3,3,3,3) # rater one's ratings
rater2 <- c(3,3,3,3) # rater one's ratings
tmp <- psych::cohen.kappa(cbind(rater1,rater2)) 


psych::cohen.kappa(cbind(factor(tmp_cna3$type[which(tmpAscatCna$size > 10e6)], levels = c( "gain", "none","loss")),
                factor(tmpAscatCna$type[which(tmpAscatCna$size > 10e6)], levels = c("gain", "none", "loss"))))



tmpGreater10Mb_ck <- tryCatch(psych::cohen.kappa(cbind(tmpAscatCna$type[which(tmpAscatCna$size > 1e7)],
                                                      tmp_cna3$type[which(tmpAscatCna$size > 1e7)])),
                             error = function(w) print(NA))


tryCatch(irr::kappa2(cbind(factor(tmp_cna3$type[which(tmpAscatCna$size > 1e7)], levels = c( "gain", "none","loss")),
                           factor(tmpAscatCna$type[which(tmpAscatCna$size > 1e7)], levels = c("gain", "none", "loss"))),
                     weight = "equal"), error = function(w) print(NA))


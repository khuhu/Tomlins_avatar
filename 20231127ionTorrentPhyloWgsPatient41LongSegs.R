### original for patient 34

round2 = function(x, digits = 0) {
  posneg = sign(x)
  z = abs(x)*10^digits
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^digits
  z*posneg
}

expectedIntCalc <- function(p, log2r, pl){
  res <- (2^(log2r)*pl - 2 * (1-p))/p
  print(res)
}


listOfVcfReport <- c("Auto_user_AUS5-48-PR986to991PR10001toPR10002_CCP_221_122",
                     "Auto_user_AUS5-49-PR10003toPR10005PR903_908_930_956_CCP_222_124",
                     "Reanalysis_user_draco-19-Simpa_PR922_PGU_foreign_525_jEBJ_556")
listOfVcfReport <- paste0(listOfVcfReport, "_anno.txt")


listOfSamples <- c("PR10001_20_CCP", "PR10002_20_CCP", "PR10003_20_CCP", "PR10004_20_CCP", "PR10005_20_CCP")

finalVariantList <- read.table("/mnt/DATA5/tmp/kev/misc/Final variant matrix_Pt41.csv",
                               sep = ",", header = TRUE)

listOfVars <- finalVariantList$CALL
listOfTcs <- read.table("/mnt/DATA5/tmp/kev/misc/sriSimpa_tumor_fractions.csv",
                        sep = ",", header = TRUE)
listOfTcs <- listOfTcs[which(listOfTcs$PATIENT == "Pt41"), ]
listOfTcs$sample2 <- paste0(str_remove(listOfTcs$SAMPLE, "-"), "_20_CCP")
# listOfTcs$sample2 <- str_remove(listOfTcs$SAMPLE, "-")
listOfTcs  <- listOfTcs[match(listOfSamples, listOfTcs$sample2),]
listOfTcs2 <- listOfTcs[, c("SAMPLE", "TUMOR_CONTENT_FINAL", "sample2")]
listOfTcs2$tc2 <- listOfTcs2$TUMOR_CONTENT_FINAL
listOfTcs2$tc2[which(listOfTcs2$tc2 == "N/A")] <- 0.5
listOfTcs2$tc2[which(listOfTcs2$tc2 < 0.5)] <- 0.5
listOfTcs2$tc2 <- as.numeric(listOfTcs2$tc2)
listOfTcs2$tc2[which(listOfTcs2$TUMOR_CONTENT_FINAL == "Small")] <- 0.5




### script used to create ssm_data.txt
### create ssm.txt first- mu_r and mu_v are tumor specific, but can just put
### their Illumina values for now

# tmpTable <- read.table(listOfVcfReport[5], sep = "\t", header = TRUE)


phyloWgs_ssm_filt <- function(listOfVcfReport, listOfSamples, listOfVars){
  tmpReport <- NULL
  setwd("/mnt/DATA6/mouseData/reportAnno/")
  for (i in listOfVcfReport) {
    tmp <- read.table(i, sep = "\t", stringsAsFactors = FALSE,
                      header = TRUE)
    ### only needed so far b/c naming conventions
    if (length(which(grepl("-", tmp$Sample))) > 0) {
      tmp$Sample <- paste0(str_remove(tmp$Sample, "-"), "_20_CCP")
    }
    tmpReport <- rbind(tmpReport, tmp)
  }
  tmpReport$Sample <- str_remove(tmpReport$Sample, "_20_PANGU")
  tmpReport2 <- tmpReport[which(tmpReport$Sample %in% listOfSamples), ]
  
  ### next create matrix for each unique mutation
  tmpReport2$str <- paste(tmpReport2$Gene.refGene, tmpReport2$Start, tmpReport2$End,
                          tmpReport2$Ref, tmpReport2$Alt)
  res <- data.frame("id" = paste0("s", 0:(length(unique(tmpReport2$str)) -1)),
                    "gene" = unique(tmpReport2$str), "a" = rep(NA, length(unique(tmpReport2$str))),
                    "d" = rep(NA, length(unique(tmpReport2$str))), "mu_r" = rep(0.999, length(unique(tmpReport2$str))),
                    "mu_v" = rep(0.499, length(unique(tmpReport2$str))),
                    "tmpStr" = unique(tmpReport2$str))
  res$filter_str <- tmpReport2$AAChange.refGene[match(res$tmpStr, tmpReport2$str)]
  res$chr <- str_remove(tmpReport2$Chr[match(res$tmpStr, tmpReport2$str)], "chr")
  res$start <- tmpReport2$Start[match(res$tmpStr, tmpReport2$str)]
  res$end <- tmpReport2$End[match(res$tmpStr, tmpReport2$str)]
  
  ### ony filt
  tmpVar <- str_split(listOfVars, " ")
  res <- res[grep(paste0(unlist(lapply(tmpVar, "[", 1)), collapse = "|"), res$filter_str), ]
  res <- res[grep(paste0(unlist(lapply(tmpVar, "[", 2)), collapse = "|"), res$filter_str), ]
  
  ###  need to calculate and d for each sample and variant
  resA <- list()
  resD <- list()
  for (i in seq_along(listOfSamples)) {
    tmpDf <- tmpReport2[which(tmpReport2$Sample == listOfSamples[i]),]
    ### alternative allele count
    tmpD <- round(tmpDf$FDP[match(res$gene, tmpDf$str)])
    naRef <- round(median(tmpD[-which(is.na(tmpD))]))
    tmpD[is.na(tmpD)] <- naRef
    ### I would need to rerun tvc so I force call variants for all variant positions
    tmpA <-  round(tmpD - tmpDf$FAO[match(res$gene, tmpDf$str)])
    
    
    ### special case for sample with no variants at all
    if (listOfSamples[i] == "PR10003_20_CCP") {
      tmpA[is.na(tmpA)] <- rep(400, length(tmpA))
      tmpD[is.na(tmpD)] <- rep(400, length(tmpD))
    } else{
      tmpA[is.na(tmpA)] <- naRef
    }
    ### reason for this is if no variant A=D
    
    
    
    resA[[i]] <- tmpA
    resD[[i]] <- tmpD
  }
  for (i in 1:nrow(res)) {
    res$a[i] <- paste0(unlist(lapply(resA, "[", i)), collapse = ",")
    res$d[i] <- paste0(unlist(lapply(resD, "[", i)), collapse = ",")
  }
  return(res)
}

### still need pangu variants ... oddly not here
ssmRes <- phyloWgs_ssm_filt(listOfVcfReport = listOfVcfReport,
                            listOfSamples = listOfSamples,
                            listOfVars = listOfVars)



listOfReport <- str_remove(listOfVcfReport, "\\_anno.txt")
### only for patient 34
listOfReport <- listOfReport[1:2]
allSeg <- NULL
for (i in listOfReport) {
  setwd("/mnt/DATA6/mouseData/copynumber/")
  setwd(i)
  tmpTable <- read.table("./segResults.txt", sep = "\t", header = TRUE)
  allSeg <- rbind(allSeg, tmpTable)
}

allSeg <- allSeg[which(allSeg$ID %in% listOfSamples), ]
allSeg <- allSeg[which(abs(allSeg$seg.mean) > 0.32 & allSeg$q.val2 < 0.05), ]
allSeg$width <- (allSeg$loc.end - allSeg$loc.start)/1e6

allSeg$str <- paste(allSeg$ID, allSeg$chrom, allSeg$loc.start, allSeg$loc.end)
# allSeg <- allSeg[-which(duplicated(allSeg$str)), ]
### questions do we include gene level events? i.e segment 2-3 markers


allSeg_filt <- allSeg[which(allSeg$num.mark > 8), ]
### min size of arm for autosomes
allSeg_filt <- allSeg_filt[which(allSeg_filt$width > 11.8), ]

tcVector <- listOfTcs2$tc2[match(allSeg_filt$ID, listOfTcs2$sample2)]

exptectedIntCn <- NULL
for (i in 1:length(allSeg_filt$seg.mean)) {
  exptectedIntCn <- c(exptectedIntCn, expectedIntCalc(tcVector[i], allSeg_filt$seg.mean[i], 2))
}
allSeg_filt$seg.mean2 <- exptectedIntCn


allSeg_filt$seg.mean2_round <- round2(allSeg_filt$seg.mean2)
allSeg_filt2 <- allSeg_filt[, c("ID", "chrom", "loc.start", "loc.end","num.mark", "seg.mean", "seg.mean2", "seg.mean2_round", "width", "str")]
allSeg_filt2$major_cn <- round2(allSeg_filt2$seg.mean2_round/2)
allSeg_filt2$minor_cn <- round2(allSeg_filt2$seg.mean2_round/2)
allSeg_filt2$minor_cn <- ifelse(allSeg_filt2$seg.mean2_round %% 2 == 0, allSeg_filt2$minor_cn, allSeg_filt2$minor_cn-1)


### most of the events > 5 copies are from a single gene; probably bad coverage of amplicon,
### or extremely large area

### need to do disjoin separate along with A and D counting because of disjoined areas
allSeg_filt2Gr <- GRanges(seqnames = allSeg_filt2$chrom,
                          IRanges(start = allSeg_filt2$loc.start, end = allSeg_filt2$loc.end))
allSeg_filt2Gr_disjoin <- disjoin(allSeg_filt2Gr)

ccpGcBed <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-37-PR902-926_CCP_550_210_099/amplicon.GCinput.txt",
                       sep = "\t", stringsAsFactors = FALSE, header = TRUE)
ccpAmpliconGr <- GRanges(seqnames = ccpGcBed$ChromNum,
                         IRanges(start = ccpGcBed$StartPos,
                                 end = ccpGcBed$EndPos))

ccpDisjoinOverlap <- findOverlaps(ccpAmpliconGr, allSeg_filt2Gr_disjoin)

allAmplicons <- NULL
for (i in listOfReport) {
  setwd("/mnt/DATA6/mouseData/copynumber/")
  setwd(i)
  tmpTable <- read.table("./amplicon.combinedCoverage.input.txt", sep = ",", header = TRUE)
  tmpTable <- tmpTable[match(ccpGcBed$AmpliconId, tmpTable$AmpliconId), ]
  if (i == "Auto_user_AUS5-48-PR986to991PR10001toPR10002_CCP_221_122") {
    allAmplicons <- rbind(allAmplicons, tmpTable)
  } else{
    allAmplicons <- cbind(allAmplicons, tmpTable)
  }
}

allAmplicons2 <- allAmplicons[, c(1, which(colnames(allAmplicons) %in% listOfSamples))]


allSeg_filt2Gr_disjoinDf <- data.frame(allSeg_filt2Gr_disjoin)

i <- listOfSamples[1]
allSampsDisjoinReads <- NULL
for (i in listOfSamples) {
  tmpDf <- allAmplicons2[[i]]
  tmpDfAllSeg <- allSeg_filt2[which(allSeg_filt2$ID  == i), ]
  tmpGr <- GRanges(seqnames = tmpDfAllSeg$chrom,
                   IRanges(start = tmpDfAllSeg$loc.start,
                           end = tmpDfAllSeg$loc.end))
  j <- 110
  for (j in 1:length(allSeg_filt2Gr_disjoin)) {
    tmpReads <- sum(tmpDf[which(subjectHits(ccpDisjoinOverlap) == j)])
    tmpOverlap <- findOverlaps(allSeg_filt2Gr_disjoin[j], tmpGr)
    if (length(subjectHits(tmpOverlap)) == 0) {
      ### another large caveat is we assume diploid tumors
      ### so if there is no gain or loss found we assume 1 copy of each
      tmpMajor <- 1
      tmpMinor <- 1
    } else{
      tmpMajor <- tmpDfAllSeg$major_cn[subjectHits(tmpOverlap)]
      tmpMinor <- tmpDfAllSeg$minor_cn[subjectHits(tmpOverlap)]
    }
    allSampsDisjoinReads <- rbind(allSampsDisjoinReads, data.frame("sample" = i, "chr" = allSeg_filt2Gr_disjoinDf$seqnames[j],
                                                                   "start" = allSeg_filt2Gr_disjoinDf$start[j],
                                                                   "end" = allSeg_filt2Gr_disjoinDf$end[j], "total_reads" = tmpReads,
                                                                   "major_cn" = tmpMajor, "minor_cn"= tmpMinor))
  }
}

### split the reads based on cn, a is # of ref reads, d is $ of total reads
### older versions did it incorrectly, the state calculations were incorrect i.e 0,0; 0/1; 1/1 were not correct
### that's why they looked funky
allSampsDisjoinReads$a <- ifelse(allSampsDisjoinReads$major_cn + allSampsDisjoinReads$minor_cn > 2, 
                                 round(allSampsDisjoinReads$minor * (allSampsDisjoinReads$total_reads/(allSampsDisjoinReads$major_cn + allSampsDisjoinReads$minor_cn))),
                                 round(((allSampsDisjoinReads$minor_cn + allSampsDisjoinReads$major_cn) * allSampsDisjoinReads$total_reads)/2))
allSampsDisjoinReads$d <- allSampsDisjoinReads$total_reads

### creating each individual disjoined entry - now I need to condense 

allSampsDisjoinReads$str <- paste(allSampsDisjoinReads$chr, allSampsDisjoinReads$start, allSampsDisjoinReads$end,
                                  allSampsDisjoinReads$major_cn, allSampsDisjoinReads$minor_cn)
allSampsDisjoinReads$str2 <- paste(allSampsDisjoinReads$chr, allSampsDisjoinReads$start, allSampsDisjoinReads$end)

### this is fine
preResCnv <- NULL
preResCnvCoords <- NULL
i <- unique(allSampsDisjoinReads$str2)[3]
for (i in unique(allSampsDisjoinReads$str2)) {
  tmpDf <- allSampsDisjoinReads[which(allSampsDisjoinReads$str2 == i),]
  tmpDf2 <- tmpDf[match(listOfSamples, tmpDf$sample),]
  
  tmpA <- paste0(tmpDf2$a, collapse = ",")
  tmpD <- paste0(tmpDf2$d, collapse = ",")
  
  j <- unique(tmpDf2$str)[2]
  for(j in unique(tmpDf2$str)){
    tmpDf3 <- tmpDf2[which(tmpDf2$str == j), ]
    
    preResCnv <- rbind(preResCnv, data.frame("a" = tmpA, "d" = tmpD,
                                             "physical_cnvs" = paste0("chrom=", tmpDf3$chr[1], ",", 
                                                                      "start=", tmpDf3$start[1], ",",
                                                                      "end=", tmpDf3$end[1], ",",
                                                                      "major_cn=", tmpDf3$major_cn[1],",",
                                                                      "minor_cn=", tmpDf3$minor_cn[1], ",",
                                                                      "cell_prev=", paste0(listOfTcs2$tc2, collapse = "|"))))
    preResCnvCoords <- rbind(preResCnvCoords, data.frame("chr" = tmpDf3$chr[1], "start" = tmpDf3$start[1], "end" = tmpDf3$end[1]))
    
  }
}




preResCnvCoordsGr <- GRanges(seqnames = preResCnvCoords$chr, 
                             IRanges(start = preResCnvCoords$start,
                                     end = preResCnvCoords$end))
preResCnvCoordsGr <- GRanges(seqnames = preResCnvCoords$chr, 
                             IRanges(start = preResCnvCoords$start,
                                     end = preResCnvCoords$end))

### need to add cnv identifier and lastly the ssms, for ssms need to first estimate the major and minor copies there
ssmResGr <- GRanges(seqnames = ssmRes$chr,
                    IRanges(start = ssmRes$start, end = ssmRes$end))
ssmRes_filt <- ssmRes[, c("id", "gene", "a", "d", "mu_r", "mu_v")]
ssmRes_filt$id <- paste0("s", 0:(nrow(ssmRes_filt)-1))

preResCnv$ssms <- ""
ssmOverlap <- findOverlaps(preResCnvCoordsGr, ssmResGr)
preResCnv$ssms[queryHits(ssmOverlap)] <- ssmRes_filt$id[subjectHits(ssmOverlap)]

### more manual part, look at overlapping ssms and assign accordingly cns
### s2, s3, s5 all (1, 0)
### determined by 
preResCnv$cnv <- paste0("c", 0:(nrow(preResCnv) - 1))
preResCnvFinal <- preResCnv[, c("cnv", "a", "d", "ssms", "physical_cnvs")]


# write.table(preResCnvFinal, "/mnt/DATA5/tmp/kev/misc/20231202patient41LongSegPhylowgsCnv.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# 
# write.table(ssmRes_filt, "/mnt/DATA5/tmp/kev/misc/20231202patient41LongSegPhylowgsSsm.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# 




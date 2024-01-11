### original for patient 34
### I think it may be most effective if for cnvs only use large, almost aneuploidy like changes
### and gene based

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



listOfVcfReport <- c("Auto_user_AUS5-37-PR902-926_CCP_550_210_099",
                     "Auto_user_AUS5-38-PR902_905_907to908_923to925_PR927_CCP_211_101",
                     "Auto_user_AUS5-39-P928toPR935_CCP_212_103",
                     "Auto_user_AUS5-40-PR936toPR943_CCP_213_105",
                     "Reanalysis_user_draco-19-Simpa_PR922_PGU_foreign_525_jEBJ_556",
                     "Reanalysis_user_cetus-24-Simpa_PR930_PGU_foreign_513_cjQM_544")
listOfVcfReport <- paste0(listOfVcfReport, "_anno.txt")


listOfSamples <- c("PR922_20_CCP","PR923_20_CCP", "PR924_20_CCP", "PR925_20_CCP",
                   "PR937_20_CCP", "PR938_20_CCP", "PR939_20_CCP",
                   "PR929_20_CCP","PR930_20_CCP")

finalVariantList <- read.table("/mnt/DATA5/tmp/kev/misc/Final variant matrix_Pt33.csv",
                               sep = ",", header = TRUE)

listOfVars <- finalVariantList$CALL
listOfVars[which(listOfVars == "ATM p.L1809fs")] <- "ATM p.T1810Lfs"
listOfTcs <- read.table("/mnt/DATA5/tmp/kev/misc/sriSimpa_tumor_fractions.csv",
                        sep = ",", header = TRUE)
listOfTcs <- listOfTcs[which(listOfTcs$PATIENT == "Pt33"), ]
listOfTcs$sample2 <- paste0(str_remove(listOfTcs$SAMPLE, "-"), "_20_CCP")
listOfTcs  <- listOfTcs[match(listOfSamples, listOfTcs$sample2),]
listOfTcs2 <- listOfTcs[, c("SAMPLE", "TUMOR_CONTENT_FINAL", "sample2")]
listOfTcs2$tc2 <- listOfTcs2$TUMOR_CONTENT_FINAL
listOfTcs2$tc2[which(listOfTcs2$tc2 == "N/A")] <- 0.5
listOfTcs2$tc2[which(listOfTcs2$tc2 < 0.5)] <- 0.5
listOfTcs2$tc2 <- as.numeric(listOfTcs2$tc2)




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
    if (listOfSamples[i] == "PR929_20_CCP") {
      tmpA[is.na(tmpA)] <- rep(908, length(tmpA))
      tmpD[is.na(tmpD)] <- rep(908, length(tmpD))
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
listOfReport <- listOfReport[1:4]
# allSeg <- NULL
# for (i in listOfReport) {
#   setwd("/mnt/DATA6/mouseData/copynumber/")
#   setwd(i)
#   tmpTable <- read.table("./segResults.txt", sep = "\t", header = TRUE)
#   allSeg <- rbind(allSeg, tmpTable)
# }
# 
# allSeg <- allSeg[which(allSeg$ID %in% listOfSamples), ]
# allSeg <- allSeg[which(abs(allSeg$seg.mean) > 0.32 & allSeg$q.val2 < 0.05), ]
# allSeg$width <- (allSeg$loc.end - allSeg$loc.start)/1e6

### questions do we include gene level events? i.e segment 2-3 markers


# allSeg_filt <- allSeg[which(allSeg$num.mark > 8), ]
# vectorTc <- listOfTcs2$tc2[match(allSeg_filt$ID, listOfTcs2$sample2)]
### (1) get rid of log2 transform
### (2)  with tumor/normal ratio get rid of normal signal - multiply by 2
### essentially what we divide by after correcting for panel bias

# allSeg_filt$seg.mean2 <- (2^(allSeg_filt$seg.mean * 2) / vectorTc)
# 
# allSeg_filt$seg.mean2_round <- round2(allSeg_filt$seg.mean2)
# allSeg_filt2 <- allSeg_filt[, c("ID", "chrom", "loc.start", "loc.end","num.mark", "seg.mean", "seg.mean2", "seg.mean2_round", "width")]
# allSeg_filt2$major_cn <- round2(allSeg_filt2$seg.mean2_round/2)
# allSeg_filt2$minor_cn <- round2(allSeg_filt2$seg.mean2_round/2)
# allSeg_filt2$minor_cn <- ifelse(allSeg_filt2$seg.mean2_round %% 2 == 0, allSeg_filt2$minor_cn, allSeg_filt2$minor_cn-1)



# hg19cyto <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genome/hg19/cytoBand.txt", sep = "\t",
#                        col.names = c("chrom","chromStart","chromEnd","name","gieStain"))


# hg19ArmLocations <- NULL
# for (i in unique(hg19cyto$chrom)) {
#   tmp <- hg19cyto[which(hg19cyto$chrom == i),]
#   tmpP <- tmp[grep("p", tmp$name), ]
#   tmpQ <- tmp[grep("q", tmp$name), ]
#   hg19ArmLocations <- rbind(hg19ArmLocations, data.frame("chrom" = rep(i, 2), "arm" = c("p", "q"),
#                                                          "start" = c(min(tmpP$chromStart), min(tmpQ$chromStart)),                                                        "end" = c(max(tmpP$chromEnd), max(tmpQ$chromEnd))))
# }
# 
# hg19ArmLocations$width <- (hg19ArmLocations$end - hg19ArmLocations$start)/1e6

### can see almost all arms > 20Mb 

# allSeg_filt2 <- allSeg_filt2[which(allSeg_filt2$width > 20), ]



### most of the events > 5 copies are from a single gene; probably bad coverage of amplicon,
### or extremely large area
# allSeg_filt2 <- allSeg_filt2[-which(allSeg_filt2$seg.mean2_round > 3),]
# 
# if (length(which(allSeg_filt2$seg.mean2_round == 2))) {
#   allSeg_filt2 <- allSeg_filt2[-which(allSeg_filt2$seg.mean2_round == 2), ]
# }


### further filtering odd cases with multiple overlapping segmentation

# allSeg_filt2$str <- paste(allSeg_filt2$ID, allSeg_filt2$chrom, allSeg_filt2$loc.start, allSeg_filt2$loc.end)
# allSeg_filt2 <- allSeg_filt2[-which(duplicated(allSeg_filt2$str)), ]
# 
# badSegStr <- c("PR923_20_CCP 1 164761908 186649320", "PR923_20_CCP 8 57080732 145740486",
#                "PR924_20_CCP 8 30915961 57080732", "PR924_20_CCP 8 71025806 145742887",
#                "PR930_20_CCP 15 57355973 99500619")
# 
# allSeg_filt2 <- allSeg_filt2[-which(allSeg_filt2$str %in% badSegStr), ]


### get gene based call then disjoin
### complicated because there are 409 gene targets on ccp, so might hyper segment some of the arms
### when disjoining

### so alternative is only gene

listOfReport <- listOfReport[1:4]
allSeg <- NULL
for (i in listOfReport) {
  setwd("/mnt/DATA6/mouseData/copynumber/")
  setwd(i)
  tmpTable <- read.table("./combinedCalls.txt", sep = "\t", header = TRUE)
  if (i == listOfReport[1]) {
    allSeg <- tmpTable
  } else{
    allSeg <- rbind(allSeg, tmpTable)
  }
}

allSeg$str <- paste(allSeg$Gene, allSeg$Sample)
allSeg <- allSeg[-which(duplicated(allSeg$str)), ]

allSeg2 <- allSeg[which(abs(log2(allSeg$CopyNumberRatio)) > 0.32 & allSeg$Log10QValue < -1.30103), ]
allSeg2 <- allSeg2[which(allSeg2$Sample %in% listOfSamples), ]


ccpGcBed <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-37-PR902-926_CCP_550_210_099/amplicon.GCinput.txt",
                       sep = "\t", stringsAsFactors = FALSE, header = TRUE)
ccpGcBedRed <- NULL
i <- unique(ccpGcBed$Gene)[1]
for (i in unique(ccpGcBed$Gene)) {
  tmpDf <- ccpGcBed[which(ccpGcBed$Gene == i), ]
  ccpGcBedRed <- rbind(ccpGcBedRed, data.frame("Gene" = i,
                                               "chrom" = tmpDf$ChromNum[1],
                                               "start" = min(tmpDf$StartPos),
                                               "end" = max(tmpDf$EndPos)))
}

allSeg2$chrom <- ccpGcBedRed$chrom[match(allSeg2$Gene, ccpGcBedRed$Gene)]
allSeg2$loc.start <- ccpGcBedRed$start[match(allSeg2$Gene, ccpGcBedRed$Gene)]
allSeg2$loc.end <- ccpGcBedRed$end[match(allSeg2$Gene, ccpGcBedRed$Gene)]
allSeg2$width <- (allSeg2$loc.end - allSeg2$loc.start)/1e6
allSeg2$str <- paste(allSeg2$chrom, allSeg2$loc.start, allSeg2$loc.end)

tcVector <- listOfTcs2$tc2[match(allSeg2$Sample, listOfTcs2$sample2)]

allSeg2$seg.mean <- log2(allSeg2$CopyNumberRatio)

exptectedIntCn <- NULL
for (i in 1:length(allSeg2$seg.mean)) {
  exptectedIntCn <- c(exptectedIntCn, expectedIntCalc(tcVector[i], allSeg2$seg.mean[i], 2))
}
allSeg2$seg.mean2 <- exptectedIntCn

allSeg2$seg.mean2_round <- round2(allSeg2$seg.mean2)

colnames(allSeg2)[2]  <- "ID"
allSeg2 <- allSeg2[, c("ID", "Gene","chrom", "loc.start", "loc.end", "seg.mean", "seg.mean2", "seg.mean2_round")]
allSeg2$major_cn <- round2(allSeg2$seg.mean2_round/2)
allSeg2$minor_cn <- round2(allSeg2$seg.mean2_round/2)
allSeg2$minor_cn <- ifelse(allSeg2$seg.mean2_round %% 2 == 0, allSeg2$minor_cn, allSeg2$minor_cn-1)
allSeg2$str2 <- paste(allSeg2$chrom, allSeg2$loc.start, allSeg2$loc.end)
allSeg3 <- allSeg2[-which(duplicated(allSeg2$str2)), ]
# allSeg2 <- allSeg2[-which(allSeg2$seg.mean2_round == 2), ]

### need to do disjoin separate along with A and D counting because of disjoined areas
# allSeg2Gr <- GRanges(seqnames = allSeg2$chrom,
#                           IRanges(start = allSeg2$loc.start, end = allSeg2$loc.end))

### got rid of disjoin because I'm measureing just genes

allSeg2Gr_disjoin <- GRanges(seqnames = allSeg3$chrom,
                             IRanges(start = allSeg3$loc.start, end = allSeg3$loc.end))

ccpGcBed <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-37-PR902-926_CCP_550_210_099/amplicon.GCinput.txt",
                       sep = "\t", stringsAsFactors = FALSE, header = TRUE)
ccpAmpliconGr <- GRanges(seqnames = ccpGcBed$ChromNum,
                         IRanges(start = ccpGcBed$StartPos,
                                 end = ccpGcBed$EndPos))

ccpDisjoinOverlap <- findOverlaps(ccpAmpliconGr, allSeg2Gr_disjoin)

allAmplicons <- NULL
for (i in listOfReport) {
  setwd("/mnt/DATA6/mouseData/copynumber/")
  setwd(i)
  tmpTable <- read.table("./amplicon.combinedCoverage.input.txt", sep = ",", header = TRUE)
  tmpTable <- tmpTable[match(ccpGcBed$AmpliconId, tmpTable$AmpliconId), ]
  if (i == "Auto_user_AUS5-37-PR902-926_CCP_550_210_099") {
    allAmplicons <- rbind(allAmplicons, tmpTable)
  } else{
    allAmplicons <- cbind(allAmplicons, tmpTable)
  }
}

allAmplicons2 <- allAmplicons[, c(1, which(colnames(allAmplicons) %in% listOfSamples))]


allSeg2Gr_disjoinDf <- data.frame(allSeg2Gr_disjoin)

i <- listOfSamples[1]
allSampsDisjoinReads <- NULL
for (i in listOfSamples) {
  tmpDf <- allAmplicons2[[i]]
  tmpDfAllSeg <- allSeg2[which(allSeg2$ID  == i), ]
  tmpGr <- GRanges(seqnames = tmpDfAllSeg$chrom,
                   IRanges(start = tmpDfAllSeg$loc.start,
                           end = tmpDfAllSeg$loc.end))
  j <- 110
  for (j in 1:length(allSeg2Gr_disjoin)) {
    tmpReads <- sum(tmpDf[which(subjectHits(ccpDisjoinOverlap) == j)])
    tmpOverlap <- findOverlaps(allSeg2Gr_disjoin[j], tmpGr)
    if (length(subjectHits(tmpOverlap)) == 0) {
      ### another large caveat is we assume diploid tumors
      ### so if there is no gain or loss found we assume 1 copy of each
      tmpMajor <- 1
      tmpMinor <- 1
    } else{
      tmpMajor <- tmpDfAllSeg$major_cn[subjectHits(tmpOverlap)]
      tmpMinor <- tmpDfAllSeg$minor_cn[subjectHits(tmpOverlap)]
    }
    allSampsDisjoinReads <- rbind(allSampsDisjoinReads, data.frame("sample" = i, "chr" = allSeg2Gr_disjoinDf$seqnames[j],
                                                                   "start" = allSeg2Gr_disjoinDf$start[j],
                                                                   "end" = allSeg2Gr_disjoinDf$end[j], "total_reads" = tmpReads,
                                                                   "major_cn" = tmpMajor, "minor_cn"= tmpMinor))
  }
}

### split the reads based on cn, a is # of ref reads, d is $ of total reads

allSampsDisjoinReads$a <- round(allSampsDisjoinReads$minor * (allSampsDisjoinReads$total_reads/(allSampsDisjoinReads$major_cn + allSampsDisjoinReads$minor_cn)))
allSampsDisjoinReads$d <- allSampsDisjoinReads$total_reads
### because its for reads that deleted on all allele for tumor sample
allSampsDisjoinReads$a[which(is.na(allSampsDisjoinReads$a))] <- 0


### creating each individual disjoined entry - now I need to condense 

allSampsDisjoinReads$str <- paste(allSampsDisjoinReads$chr, allSampsDisjoinReads$start, allSampsDisjoinReads$end,
                                  allSampsDisjoinReads$major_cn, allSampsDisjoinReads$minor_cn)
allSampsDisjoinReads$str2 <- paste(allSampsDisjoinReads$chr, allSampsDisjoinReads$start, allSampsDisjoinReads$end)


### creating cnvs file with differing states of the sample positions
###
### 
# preResCnv <- NULL
# preResCnvCoords <- NULL
# i <- unique(allSampsDisjoinReads$str)[110]
# for (i in unique(allSampsDisjoinReads$str)) {
#   tmpDf <- allSampsDisjoinReads[which(allSampsDisjoinReads$str == i),]
#   tmpDf2 <- tmpDf[match(listOfSamples, tmpDf$sample),]
#   tmpDf2$sample <- listOfSamples
#   
#   tmpReads <- round((sum(tmpDf2$a[-which(is.na(tmpDf2$a))]) + sum(tmpDf2$d[-which(is.na(tmpDf2$d))]))/ (2*length(which(!is.na(tmpDf2$a)))))
#   
#   ### get rid of the weird infinity reads - when both ref and alt have 0 copies
#   if (is.infinite(tmpReads)) {
#     next()
#   }
#   tmpDf2$a <- ifelse(is.na(tmpDf2$a), tmpReads, tmpDf2$a)
#   tmpDf2$d <- ifelse(is.na(tmpDf2$d), tmpReads, tmpDf2$d)
#   preResCnv <- rbind(preResCnv, data.frame("a" = paste0(tmpDf2$a, collapse = ","), "d" = paste0(tmpDf2$d, collapse = ","),
#                                            "physical_cnvs" = paste0("chrom=", tmpDf$chr[1], ",", 
#                                                                     "start=", tmpDf$start[1], ",",
#                                                                     "end=", tmpDf$end[1], ",",
#                                                                     "major_cn=", tmpDf$major_cn[1],",",
#                                                                     "minor_cn=", tmpDf$minor_cn[1], ",",
#                                                                     "cell_prev=", paste0(listOfTcs2$tc2, collapse = "|"))))
#   preResCnvCoords <- rbind(preResCnvCoords, data.frame("chr" = tmpDf$chr[1], "start" = tmpDf$start[1], "end" = tmpDf$end[1]))
# }

### doing it by position instead of cn state


preResCnv <- NULL
preResCnvCoords <- NULL
i <- unique(allSampsDisjoinReads$str2)[110]
for (i in unique(allSampsDisjoinReads$str2)) {
  tmpDf <- allSampsDisjoinReads[which(allSampsDisjoinReads$str2 == i),]
  tmpDf2 <- tmpDf[match(listOfSamples, tmpDf$sample),]
  
  tmpA <- paste0(tmpDf2$a, collapse = ",")
  tmpD <- paste0(tmpDf2$d, collapse = ",")
  
  j <- unique(tmpDf2$str)[1]
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
preResCnv$ssms[which(preResCnv$ssm == "s13")] <- "s13,1,2"
preResCnv$ssms[which(preResCnv$ssm == "s10")] <- "s10,0,1"
preResCnv$ssms[which(preResCnv$ssm == "s8")] <- "s8,1,2"
preResCnv$ssms[which(preResCnv$ssm == "s6")] <- "s6,1,1"
preResCnv$ssms[which(preResCnv$ssm == "s4")] <- "s4,1,1"
preResCnv$ssms[which(preResCnv$ssm == "s2")] <- "s2,1,1"

### first time it shows a ssm predates a cnv


preResCnvFinal <- preResCnv[, c("cnv", "a", "d", "ssms", "physical_cnvs")]


# write.table(preResCnvFinal, "/mnt/DATA5/tmp/kev/misc/20231129patient33phylowgsGeneV2Cnv.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# 
# write.table(ssmRes_filt, "/mnt/DATA5/tmp/kev/misc/20231129patient33phylowgsGeneV2Ssm.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)





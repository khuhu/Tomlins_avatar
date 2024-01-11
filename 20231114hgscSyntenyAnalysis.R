### same pipline but for hgsc - from 20231027allSyntenyAnalysisCoad.R

### note: if rmd is made later, then can just split below code into block based on processing mouse, human, graphs, permutation etc

source("/home/kevhu/scripts/20210802syntenyFunctions.R")
source("/home/kevhu/scripts/20230713newSyntenyDfFunctions.R")
source("/home/kevhu/scripts/20230713newSyntenyDfFunctionV2.R")

### q should be ploidy
segResFilt4 <- read.table("/mnt/DATA5/tmp/kev/misc/20230907segRes2.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
segResFilt4$chrom <- paste0("chr", segResFilt4$chrom)
segResFilt4$size <- segResFilt4$loc.end - segResFilt4$loc.start

allNgsDis_all <- disjoin(GRanges(seqnames = segResFilt4$chrom,
                                 IRanges(start = segResFilt4$loc.start, end = segResFilt4$loc.end)))

### all there minus x chromosome I forgot to remove

allNgsDisDf_all <- data.frame(allNgsDis_all)
allNgsDisDf_all <- allNgsDisDf_all[-which(allNgsDisDf_all$width < 2),]

allNgsDisDf_allGr <- GRanges(seqnames = allNgsDisDf_all$seqnames,
                             IRanges(start = allNgsDisDf_all$start, end = allNgsDisDf_all$end))

ngsDisjoin2 <- NULL
for (i in unique(segResFilt4$ID)) {
  print(i)
  tmpNgs <- segResFilt4[which(segResFilt4$ID == i),]
  tmpRes <- allNgsDisDf_all[,1:3]
  tmpRes$sample <- i
  queryIdx <- 1:nrow(tmpRes)
  tmpGrange <- GRanges(seqnames = tmpNgs$chrom, IRanges(start = tmpNgs$loc.start, end = tmpNgs$loc.end))
  tmpIdx <- subjectHits(findOverlaps(query = allNgsDisDf_allGr, subject = tmpGrange))
  testQuery <- queryHits(findOverlaps(query = allNgsDisDf_allGr, subject = tmpGrange))
  
  
  if (length(which(duplicated(queryHits(findOverlaps(query = allNgsDisDf_allGr, subject = tmpGrange))))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated(queryHits(findOverlaps(query = allNgsDisDf_allGr, subject = tmpGrange))))]
    testQuery <- testQuery[-which(duplicated(queryHits(findOverlaps(query = allNgsDisDf_allGr, subject = tmpGrange))))]
  }
  ### error probably happens b/c there is no match for one of the queries
  tmpRes$cn <- NA
  tmpRes$cn[testQuery] <- tmpNgs$seg.mean[tmpIdx]
  tmpRes2 <- tmpRes[, c("sample", "seqnames", "start", "end", "cn")]
  ngsDisjoin2  <- rbind(ngsDisjoin2 , tmpRes2)
}

colnames(ngsDisjoin2)[c(2,5)] <- c("chrom", "seg.mean")
ngsDisjoin2$string <- paste(ngsDisjoin2$chrom, ngsDisjoin2$start, ngsDisjoin2$end)
ngsDisjoin2$size <- ngsDisjoin2$end - ngsDisjoin2$start


tmp_cna <- NULL
tmpDisJoinA <- read.table("/mnt/DATA5/tmp/kev/misc/20230907ngsInter_0_10Mb.bed",
                          sep = "\t", header = FALSE, stringsAsFactors = FALSE)

tmpDisJoin <- read.table(paste0("/mnt/DATA5/tmp/kev/misc/20230907ngsInter_", 3, "_", 10, "Mb.bed"),
                         sep = "\t", header = FALSE, stringsAsFactors = FALSE)

colnames(tmpDisJoin) <- c("chrom", "start", "end","sample", "seg.mean", "str", "size")
colnames(tmpDisJoinA) <- c("chrom", "start", "end","sample", "seg.mean", "str", "size")
tmpDisJoinAneuploidy <- matrix(0, nrow = length(unique(segResFilt4$ID)), ncol = length(unique(segResFilt4$chrom)))
rownames(tmpDisJoinAneuploidy) <- unique(segResFilt4$ID)
colnames(tmpDisJoinAneuploidy) <- unique(segResFilt4$chrom)
tmpCnaAneu <- NULL

k <- unique(segResFilt4$ID)[1]
for (k in rownames(tmpDisJoinAneuploidy)) {
  tmpDf <- tmpDisJoin[grep(k, tmpDisJoin$sample),]
  tmpDfA <- tmpDisJoinA[grep(k, tmpDisJoinA$sample),]
  tmpChrList <- unique(tmpDf$chrom)
  
  for (l in seq_along(tmpChrList)) {
    tmpChr <- tmpDf[which(tmpDf$chrom == tmpChrList[l]), ]
    tmpChr2 <- tmpChr
    tmpChrA <- tmpDfA[which(tmpDfA$chrom == tmpChrList[l]), ]
    tmpChr3 <- tmpChrA
    chrThres <- sum(tmpChr$size) * 0.7
    tmpAneu <- "none"
    if (sum(tmpChrA$size[which(tmpChrA$seg.mean > 0.2)]) > chrThres) {
      tmpAneu <- "gain"
      tmpChr2$seg.mean[which(tmpChr2$seg.mean > 0.2)] <- 0
      tmpChr3$seg.mean[which(tmpChr3$seg.mean < 0.2)] <- 0
    } else if(sum(tmpChrA$size[which(tmpChrA$seg.mean < -0.2)]) > chrThres){
      tmpAneu <- "loss"
      tmpChr2$seg.mean[which(tmpChr2$seg.mean < -0.2)] <- 0
      tmpChr3$seg.mean[which(tmpChr3$seg.mean > -0.2)] <- 0
    } else{
      tmpAneu <- "none"
    }
    tmpDisJoinAneuploidy[k, l] <- tmpAneu
    tmp_cna <- rbind(tmp_cna, tmpChr2)
    tmpCnaAneu <- rbind(tmpCnaAneu, tmpChr3)
  }
}

tmp_cna$seg.mean[which(abs(tmp_cna$seg.mean) < 0.2)] <- 0
tmp_cna$seg.mean[which(is.na(tmp_cna$seg.mean))] <- 0

tmpCnaAneu$seg.mean[which(abs(tmpCnaAneu$seg.mean) < 0.2)] <- 0
tmpCnaAneu$seg.mean[which(is.na(tmpCnaAneu$seg.mean))] <- 0

tmp_cna2 <- NULL
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

allMouseAneuploidy <- tmpDisJoinAneuploidy[,1:19]
allMouseCna <- tmp_cna2
allMouseAneu <- tmpCnaAneu

p53String <- "chr11 69588341 69588439"

expectedLogRCalc <- function(p, na, nb, pl){
  ### p is purity, pl is ploidy and nb and na are copies of parental allele
  res <- log2( (2* (1-p) + p* (na + nb)) / pl )
}

tmpRes <- expectedLogRCalc(0.3, 0, 0, 2)

ngsDisjoinP53 <- ngsDisjoin2[which(ngsDisjoin2$string == p53String),]

### cutoff if diploid, above shows loss of both copies at diploid and 30% tumor content is log2RR of -0.5
goodTcHgsc <- ngsDisjoinP53$sample[which(ngsDisjoinP53$seg.mean < -0.5)]




mm <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/cytoBand.txt.gz", 
            col.names = c("chrom","chromStart","chromEnd","name","gieStain"))

mm10ChromArms <- mm[ , .(length = sum(chromEnd - chromStart)), 
                     by = .(chrom, arm = substring(name, 1, 1)) ]

tmpChrBoundaries <- data.frame("chr" = mm10ChromArms$chrom, "start" = 1,
                               "end" = mm10ChromArms$length)

tmpChrBoundaries <- tmpChrBoundaries[which(tmpChrBoundaries$chr %in% paste0("chr", 1:19)), ]
tmpChrBoundaries <- tmpChrBoundaries[match(paste0("chr", 1:19), tmpChrBoundaries$chr), ]


allMouseAneuploidyMelt <- melt(allMouseAneuploidy)
allMouseAneuploidyNum <- NULL
for (i in unique(allMouseAneuploidyMelt$Var1)) {
  tmp <- allMouseAneuploidyMelt[which(allMouseAneuploidyMelt$Var1 == i), ]
  tmpRes <- data.frame("sample" = i, tmpChrBoundaries)
  tmpRes$seg.mean <- 0
  tmpRes$seg.mean[which(tmp$value == "loss")] <- -1
  tmpRes$seg.mean[which(tmp$value == "gain")] <- 1
  allMouseAneuploidyNum <- rbind(allMouseAneuploidyNum, tmpRes)
}


allMouseAneu_hgsc <- allMouseAneuploidyNum[-grep("efd|pp|umm|ef", allMouseAneuploidyNum$sample),]
allMouseAneu_hgsc2 <- data.frame("sampleID" = allMouseAneu_hgsc$sample, "chrom" = as.numeric(str_remove(allMouseAneu_hgsc$chr, "chr")),
                                 "start.pos" = allMouseAneu_hgsc$start, "end.pos" = allMouseAneu_hgsc$end,
                                 "n.probes" = rep(NA , length(allMouseAneu_hgsc$sample)), "mean" = allMouseAneu_hgsc$seg.mean,
                                 "str" = paste0(allMouseAneu_hgsc$sample, str_remove("chr", allMouseAneu_hgsc$chr), allMouseAneu_hgsc$start, allMouseAneu_hgsc$end),
                                 "length" = allMouseAneu_hgsc$end - allMouseAneu_hgsc$start)


### read in different colorectal annotation and fix pathology

goodTcHgsc2 <- goodTcHgsc
goodTcHgsc2[goodTcHgsc2 == "mg1"] <- "1628lt"
goodTcHgsc2[goodTcHgsc2 == "mg4"] <- "2027lte"
goodTcHgsc2[goodTcHgsc2 == "mg5"] <- "2163lt"
goodTcHgsc2[goodTcHgsc2 == "mg7"] <- "2405rt"
goodTcHgsc2[goodTcHgsc2 == "mg12"] <- "3807lt"
goodTcHgsc2[goodTcHgsc2 == "mg15"] <- "13085rt"
goodTcHgsc2[goodTcHgsc2 == "mg16"] <- "14085lt"
goodTcHgsc2[goodTcHgsc2 == "mg19"] <- "14396lt"
goodTcHgsc2[goodTcHgsc2 == "mg20"] <- "14399rt"
goodTcHgsc2[goodTcHgsc2 == "mg22"] <- "14433lte"


hgscHistology <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/2021078kathyMouseSelfAnno.xlsx", sheet = 1)

## converting mgids from 20201207annotation

### from quick check below a few samples to add for annotation from first histology 
hgscHistology2 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220110allHgscNewPanelAnno.xlsx", sheet = 1)
hgscHistology2$Sample2 <- str_remove_all(nameStripper(hgscHistology2$Sample), "-")

annoTableHgsc <- hgscHistology2
annoTableHgsc <- rbind(annoTableHgsc, data.frame("Sample" = c("12167rt","12175rt", "12402rt", "12402lt", "2007rt", "2519lt", "6805lt"),
                                                 "Type" = c("MMMT", "HGSC", "HGSC", "HGSC", "HGSC", "HGSC", "HGSC"), 
                                                 "Geno" = c("PRN", "PRN", "PRN", "PRN", "UNK", "BPRN", "PRN"),
                                                 "Sample2" =  c("12167rt","12175rt", "12402rt", "12402lt", "2007rt", "2519lt", "6805lt")))

annoTableHgsc$Geno2 <- annoTableHgsc$Geno
annoTableHgsc$Geno2[which(annoTableHgsc$Geno == "Aging (aged)")] <- "UNK"
annoTableHgsc$Geno2[which(annoTableHgsc$Geno == "Aging (control)")] <- "UNK"
annoTableHgsc$Geno2[which(annoTableHgsc$Geno == "Pregnancy (high parity)")] <- "UNK"
annoTableHgsc$Geno2[which(annoTableHgsc$Geno == "Pregnancy (control)")] <- "UNK"
annoTableHgsc$Geno2[which(annoTableHgsc$Geno == "Pregnancy (control)")] <- "UNK"
annoTableHgsc <- annoTableHgsc[which(annoTableHgsc$Geno2 %in% c("BPRN", "BPN", "BPR", "BPP", "UNK")),]
annoTableHgsc <- annoTableHgsc[-which(annoTableHgsc$Sample2 %in% c("14656peritnealmt", "14433mt")),]

### above is large cna, below aneuploidy
allMouseAneu_hgsc2 <- allMouseAneu_hgsc2[which(allMouseAneu_hgsc2$sampleID %in% annoTableHgsc$Sample2[which(annoTableHgsc$Type == "HGSC" | annoTableHgsc$Type == "MMMT")]), ]
allMouseAneu_hgsc2 <- allMouseAneu_hgsc2[which(allMouseAneu_hgsc2$sampleID %in% goodTcHgsc2),]

goodTcHgsc2[-which(goodTcHgsc2 %in% unique(allMouseAneu_hgsc2$sampleID))]

### before I got filtered out events < 10%; important to have it for permutation method

allMouseAneu_hgsc_freq <- getFreqData(allMouseAneu_hgsc2)
allMouseAneu_hgsc_freq_res <- ampsDels(allMouseAneu_hgsc_freq)
allMouseAneu_hgsc_amp_bed <- reducingFreqBed2(allMouseAneu_hgsc_freq_res[[1]])
allMouseAneu_hgsc_del_bed <- reducingFreqBed2(allMouseAneu_hgsc_freq_res[[3]])



### annotation of apc and other mutations for GISTIC
hgsc_anno <- read.table("/mnt/DATA5/tmp/kev/misc/tp53_gsc_PATIENT_DATA_oncoprint.tsv", sep = "\t",
                        header = TRUE)
hgsc_anno2 <- hgsc_anno


hg19cyto <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genome/hg19/cytoBand.txt", sep = "\t",
                       col.names = c("chrom","chromStart","chromEnd","name","gieStain"))

broadBySampleGisticHgsc <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/broad_values_by_arm.txt",
                                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

hg19ArmLocations <- NULL
for (i in unique(hg19cyto$chrom)) {
  tmp <- hg19cyto[which(hg19cyto$chrom == i),]
  tmpP <- tmp[grep("p", tmp$name), ]
  tmpQ <- tmp[grep("q", tmp$name), ]
  hg19ArmLocations <- rbind(hg19ArmLocations, data.frame("chrom" = rep(i, 2), "arm" = c("p", "q"),
                                                         "start" = c(min(tmpP$chromStart), min(tmpQ$chromStart)),                                                        "end" = c(max(tmpP$chromEnd), max(tmpQ$chromEnd))))
}

hg19ArmLocations$chrStripped <- as.numeric(str_remove(hg19ArmLocations$chrom, "chr"))
hg19ArmLocations <- hg19ArmLocations[which(hg19ArmLocations$chrStripped %in% 1:22),]
hg19ArmLocations <- hg19ArmLocations[order(hg19ArmLocations$chrStripped,decreasing = FALSE),]
hg19ArmLocations$str <- paste0(hg19ArmLocations$chrStripped, hg19ArmLocations$arm)
hg19ArmLocations <- hg19ArmLocations[which(hg19ArmLocations$str %in% broadBySampleGisticHgsc$`Chromosome Arm`), ]


broadBySampleGisticHgsc <- broadBySampleGisticHgsc[-grep("X", broadBySampleGisticHgsc$`Chromosome Arm`), ]
broadBySampleGisticMatHgsc <- broadBySampleGisticHgsc[, 2:ncol(broadBySampleGisticHgsc)]
broadBySampleGisticHgsc2 <- cbind(hg19ArmLocations[, c("chrStripped", "start", "end")], broadBySampleGisticMatHgsc)
broadBySampleGisticHgsc2_sample <- gsub("(^.*?-.{3}?)-.*", "\\1", colnames(broadBySampleGisticHgsc2))
broadBySampleGisticHgsc2 <- broadBySampleGisticHgsc2[c(1:3, which(broadBySampleGisticHgsc2_sample %in% paste0(hgsc_anno2$Sample.ID, "A")))]


broadBySampleGisticMeltHgsc <- melt(broadBySampleGisticHgsc2, id.vars = c("chrStripped", "start", "end"))
broadBySampleGisticMeltHgsc2 <- data.frame("sampleID" = broadBySampleGisticMeltHgsc$variable, "chrom" = broadBySampleGisticMeltHgsc$chrStripped,
                                           "start.pos" = broadBySampleGisticMeltHgsc$start, "end.pos" = broadBySampleGisticMeltHgsc$end, 
                                           "n.probes" = NA, "mean" = broadBySampleGisticMeltHgsc$value)
broadBySampleGisticMeltHgsc2$str <- paste0(broadBySampleGisticMeltHgsc2$sampleID, broadBySampleGisticMeltHgsc2$chrom, 
                                           broadBySampleGisticMeltHgsc2$start.pos, broadBySampleGisticMeltHgsc2$end.pos)
broadBySampleGisticMeltHgsc2$length <- broadBySampleGisticMeltHgsc2$end.pos - broadBySampleGisticMeltHgsc2$start.pos


armGisticHgsc_freq <- getFreqData(broadBySampleGisticMeltHgsc2)
armGisticHgsc_freq_res <- ampsDels(armGisticHgsc_freq)
armGisticHgsc_amp_bed <- reducingFreqBed2(armGisticHgsc_freq_res[[1]])
armGisticHgsc_del_bed <- reducingFreqBed2(armGisticHgsc_freq_res[[3]])


hgsc_arm_allSynTable <- circosFreq2(allMouseAneu_hgsc_amp_bed, allMouseAneu_hgsc_del_bed, armGisticHgsc_amp_bed,
                                             armGisticHgsc_del_bed, filename = "20231115_aneuAll_hgsc", ref = "human")

### start of permutation

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


hgsc_arm_allSynTable$str <- paste(hgsc_arm_allSynTable$h_chr,
                                           hgsc_arm_allSynTable$h_start,
                                           hgsc_arm_allSynTable$h_end)
hgscTp53Arms <- broadBySampleGisticHgsc2
tmpSynTable_hgsc <- hgsc_arm_allSynTable[-which(duplicated(hgsc_arm_allSynTable$str)), ]
tmpSynTable_hgsc$h_chr <- as.numeric(str_remove(tmpSynTable_hgsc$h_chr, "h_chr"))
tmpSynTable_hgsc$m_chr <- as.numeric(str_remove(tmpSynTable_hgsc$m_chr, "m_chr"))
tmpSynTable_hgsc <- tmpSynTable_hgsc[order(tmpSynTable_hgsc$h_chr, tmpSynTable_hgsc$h_start, decreasing = FALSE),]
hgscTp53ArmsGr <- GRanges(seqnames = hgscTp53Arms$chrStripped, IRanges(start = hgscTp53Arms$start, end = hgscTp53Arms$end))
hgscSyntenyGrange <- GRanges(seqnames = tmpSynTable_hgsc$h_chr, IRanges(start = tmpSynTable_hgsc$h_start, end = tmpSynTable_hgsc$h_end))
hgscMm10SyntenyGrange <- GRanges(seqnames = tmpSynTable_hgsc$m_chr, IRanges(start = tmpSynTable_hgsc$m_start, end = tmpSynTable_hgsc$m_end))

humanHgscRegionMat <- data.frame("chr" = tmpSynTable_hgsc$h_chr, "start" = tmpSynTable_hgsc$h_start, "end" = tmpSynTable_hgsc$h_end)
for (i in 4:ncol(hgscTp53Arms)) {
  
  tmpIdx <- subjectHits(findOverlaps(query = hgscSyntenyGrange, subject = hgscTp53ArmsGr))
  testQuery <- queryHits(findOverlaps(query = hgscSyntenyGrange, subject = hgscTp53ArmsGr))
  
  
  if (length(which(duplicated(testQuery))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated(testQuery))]
    testQuery <- testQuery[-which(duplicated(testQuery))]
  }
  
  ### error probably happens b/c there is no match for one of the queries
  humanHgscRegionMat[ ,i] <- NA
  humanHgscRegionMat[ ,i] <- hgscTp53Arms[tmpIdx, i]
}

humanHgscRegionMat2 <- humanHgscRegionMat[ , 4:ncol(humanHgscRegionMat)]

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
                           tmp <- apply(humanHgscRegionMat2, 2, function(x) sample(x, nrow(humanHgscRegionMat2)))
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

humanPermAmpHgscStat <- do.call(cbind, synHgscPerm[[1]])
humanPermDelHgscStat <- do.call(cbind, synHgscPerm[[2]])
humanPermAmpHgscFreq <- do.call(cbind, synHgscPerm[[3]])
humanPermDelHgscFreq <- do.call(cbind, synHgscPerm[[4]])



mouseHgscRegionMat <- data.frame("chr" = tmpSynTable_hgsc$m_chr, "start" = tmpSynTable_hgsc$m_start, "end" = tmpSynTable_hgsc$m_end)
for (i in seq_along(unique(allMouseAneu_hgsc2$sampleID))) {
  
  tmp <- allMouseAneu_hgsc2[which(allMouseAneu_hgsc2$sampleID == unique(allMouseAneu_hgsc2$sampleID)[i]),]
  tmpGrange <- GRanges(seqnames = tmp$chrom, IRanges(start = tmp$start.pos, end = tmp$end.pos))
  
  print(unique(allMouseAneu_hgsc2$sampleID)[i])
  
  tmpIdx <- subjectHits(findOverlaps(query = hgscMm10SyntenyGrange, subject = tmpGrange))
  testQuery <- queryHits(findOverlaps(query = hgscMm10SyntenyGrange, subject = tmpGrange))
  
  if (length(which(duplicated(testQuery))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated( testQuery))]
    testQuery <- testQuery[-which(duplicated( testQuery))]
  }
  ### error probably happens b/c there is no match for one of the queries
  
  mouseHgscRegionMat[ , i + 3] <- NA
  mouseHgscRegionMat[testQuery, i + 3] <- tmp$mean[tmpIdx]
}

mouseHgscRegionMat2 <- mouseHgscRegionMat[ , 4:ncol(mouseHgscRegionMat)]

cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
synMmHgscPerm <- foreach(i=1:25,
                         .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                           tmpAmpAll <- NULL
                           tmpDelAll <- NULL
                           tmpAmpFreq <- NULL
                           tmpDelFreq <- NULL
                           j <- 0
                           while (j < 400) {
                             tmp <- apply(mouseHgscRegionMat2, 2, function(x) sample(x, nrow(mouseHgscRegionMat2)))
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

mousePermAmpFreq_hgsc <- do.call(cbind, synMmHgscPerm[[3]])
mousePermDelFreq_hgsc <- do.call(cbind, synMmHgscPerm[[4]])



hgsc_arm_allSynTable$str <- paste(hgsc_arm_allSynTable$h_chr, hgsc_arm_allSynTable$h_start, hgsc_arm_allSynTable$h_end)
tmpSynTable_hgsc <- hgsc_arm_allSynTable[-which(duplicated(hgsc_arm_allSynTable$str)), ]
tmpSynTable_hgsc$h_chr <- as.numeric(str_remove(tmpSynTable_hgsc$h_chr, "h_chr"))
tmpSynTable_hgsc$m_chr <- as.numeric(str_remove(tmpSynTable_hgsc$m_chr, "m_chr"))
tmpSynTable_hgsc <- tmpSynTable_hgsc[order(tmpSynTable_hgsc$h_chr, tmpSynTable_hgsc$h_start, decreasing = FALSE),]

ampSynTable_hgsc <- data.frame("chr" = tmpSynTable_hgsc$m_chr, "start" = tmpSynTable_hgsc$m_start, "end" = tmpSynTable_hgsc$m_end)
delSynTable_hgsc <- data.frame("chr" = tmpSynTable_hgsc$m_chr, "start" = tmpSynTable_hgsc$m_start, "end" = tmpSynTable_hgsc$m_end)
ampSynTable_hgsc$h_freq <- apply(humanHgscRegionMat2, 1, function(x) length(which(x > 0.2))/length(x))
ampSynTable_hgsc$m_freq <- apply(mouseHgscRegionMat2, 1, function(x) length(which(x > 0.2))/length(x))
delSynTable_hgsc$h_freq <- apply(humanHgscRegionMat2, 1, function(x) length(which(x < -0.2))/length(x))
delSynTable_hgsc$m_freq <- apply(mouseHgscRegionMat2, 1, function(x) length(which(x < -0.2))/length(x))


humanPermAmp2_hgsc <- humanPermAmpHgscFreq
mousePermAmp2_hgsc <- mousePermAmpFreq_hgsc
humanPermDel2_hgsc <- humanPermDelHgscFreq
mousePermDel2_hgsc <- mousePermDelFreq_hgsc



permStatResTbl_hgsc <- NULL
for (i in 1:396) {
  
  
  ### alternatively I can look at them mouse and human separately and only nominate ones where they're both significant
  tmpAmpZ_hfreq <- (ampSynTable_hgsc$h_freq[i] - mean(humanPermAmp2_hgsc[i, ], na.rm = TRUE))/sd(humanPermAmp2_hgsc[i, ], na.rm = TRUE)
  tmpDelZ_hfreq <- (delSynTable_hgsc$h_freq[i] - mean(humanPermDel2_hgsc[i, ], na.rm = TRUE))/sd(humanPermDel2_hgsc[i, ], na.rm = TRUE)
  
  tmpAmpZ_mfreq <- (ampSynTable_hgsc$m_freq[i] - mean(mousePermAmp2_hgsc[i, ], na.rm = TRUE))/sd(mousePermAmp2_hgsc[i, ], na.rm = TRUE)
  tmpDelZ_mfreq <- (delSynTable_hgsc$m_freq[i] - mean(mousePermDel2_hgsc[i, ], na.rm = TRUE))/sd(mousePermDel2_hgsc[i, ], na.rm = TRUE)
  
  
  
  permStatResTbl_hgsc <- rbind(permStatResTbl_hgsc,
                                   data.frame("ampZ_hfreq" = tmpAmpZ_hfreq, "delZ_hfreq" = tmpDelZ_hfreq,
                                              "ampZ_mfreq" = tmpAmpZ_mfreq, "delZ_mfreq" = tmpDelZ_mfreq))
}

permStatResTbl_hgsc$ampZ_hfreq2 <- pnorm(q=permStatResTbl_hgsc$ampZ_hfreq, lower.tail=FALSE)
permStatResTbl_hgsc$delZ_hfreq2 <- pnorm(q=permStatResTbl_hgsc$delZ_hfreq, lower.tail=FALSE)
permStatResTbl_hgsc$ampZ_mfreq2 <- pnorm(q=permStatResTbl_hgsc$ampZ_mfreq, lower.tail=FALSE)
permStatResTbl_hgsc$delZ_mfreq2 <- pnorm(q=permStatResTbl_hgsc$delZ_mfreq, lower.tail=FALSE)


hgsc_arm_allSynTable_amp <- tmpSynTable_hgsc
hgsc_arm_allSynTable_amp$h_freq <- ampSynTable_hgsc$h_freq
hgsc_arm_allSynTable_amp$m_freq <- ampSynTable_hgsc$m_freq
hgsc_arm_allSynTable_amp$hfreqP <- permStatResTbl_hgsc$ampZ_hfreq2
hgsc_arm_allSynTable_amp$hfreqAdjP <- p.adjust(hgsc_arm_allSynTable_amp$hfreqP, method = "BH",
                                                   n = 4 * length(hgsc_arm_allSynTable_amp$hfreqP))
hgsc_arm_allSynTable_amp$mfreqP <- permStatResTbl_hgsc$ampZ_mfreq2
hgsc_arm_allSynTable_amp$mfreqAdjP <- p.adjust(hgsc_arm_allSynTable_amp$mfreqP, method = "BH",
                                                   n = 4 * length(hgsc_arm_allSynTable_amp$mfreqP))
hgsc_arm_allSynTable_amp$sigCheck <- ifelse(hgsc_arm_allSynTable_amp$hfreqAdjP < 0.05 & hgsc_arm_allSynTable_amp$mfreqAdjP < 0.05,
                                                "good", "bad")


hgsc_arm_allSynTable_del <- tmpSynTable_hgsc
hgsc_arm_allSynTable_del$h_freq <- delSynTable_hgsc$h_freq
hgsc_arm_allSynTable_del$m_freq <- delSynTable_hgsc$m_freq
hgsc_arm_allSynTable_del$hfreqP <- permStatResTbl_hgsc$delZ_hfreq2
hgsc_arm_allSynTable_del$hfreqAdjP <- p.adjust(hgsc_arm_allSynTable_del$hfreqP, method = "BH",
                                                   n = 4 * length(hgsc_arm_allSynTable_del$hfreqP))
hgsc_arm_allSynTable_del$mfreqP <- permStatResTbl_hgsc$delZ_mfreq2
hgsc_arm_allSynTable_del$mfreqAdjP <- p.adjust(hgsc_arm_allSynTable_del$mfreqP, method = "BH",
                                                   n = 4 * length(hgsc_arm_allSynTable_del$mfreqP))
hgsc_arm_allSynTable_del$sigCheck <- ifelse(hgsc_arm_allSynTable_del$hfreqAdjP < 0.05 & hgsc_arm_allSynTable_del$mfreqAdjP < 0.05,
                                                "good", "bad")


### works, now do the RNA portion - load rna then look at only significant regions

### possibly rerun the analysis using htseq pipeline instead of rsem values
# xenaHtseqExpression <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/xena/TCGA-OV.htseq_counts.tsv",
#                                   sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
# tmpDf <- xenaHtseqExpression[,1:2]
# tmpDf$Ensembl_ID <- gsub('\\..+$', '', tmpDf$Ensembl_ID)

queryTCGA <- GDCquery(project = "TCGA-OV",
                      data.category = 'Transcriptome Profiling',
                      experimental.strategy = 'RNA-Seq',
                      workflow.type = 'STAR - Counts',
                      access = 'Open',
                      sample.type = 'Primary Tumor')
tcga_ov_txn <- GDCprepare(queryTCGA, summarizedExperiment = TRUE,
                          directory = "/mnt/DATA5/tmp/kev/tmpDbs/tcgaBiolinks/")
ov_txn_matrix <- assay(tcga_ov_txn, 'unstranded')

tmpDf <- rownames(ov_txn_matrix)
tmpDf <- gsub('\\..+$', '', tmpDf)
library(biomaRt)
library(gprofiler2)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- tmpDf
symbol <- getBM(filters = "ensembl_gene_id",
                attributes = c("ensembl_gene_id","hgnc_symbol"),
                values = genes, 
                mart = mart)
symbol2 <- symbol[-which(symbol$hgnc_symbol == ""),]
ov_txn_matrix2 <- ov_txn_matrix
rownames(ov_txn_matrix2) <- gsub('\\..+$', '', rownames(ov_txn_matrix2))
ov_txn_matrix2 <- ov_txn_matrix2[which(rownames(ov_txn_matrix2) %in% symbol2$ensembl_gene_id), ]
rownames(ov_txn_matrix2) <- symbol2$hgnc_symbol[match(rownames(ov_txn_matrix2), symbol2$ensembl_gene_id)]
ov_txn_matrix2 <- ov_txn_matrix2[-which(duplicated(rownames(ov_txn_matrix2))), ] 

# hgscRsemExpression <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/gdac.broadinstitute.org_OV.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/OV.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt",
#                                  sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
# 
# which(colnames(xenaHtseqExpression) %in% colnames(hgscRsemExpression))
# 
# hgscRsemExpression <- hgscRsemExpression[, c(1, which(hgscRsemExpression[1, ] == "raw_count"))]
# colnames(hgscRsemExpression)[2:ncol(hgscRsemExpression)] <- gsub("(^.*?-.{3}?)-.*", "\\1",
#                                                                  colnames(hgscRsemExpression)[2:ncol(hgscRsemExpression)])
# hgscRsemExpression <- hgscRsemExpression[2:nrow(hgscRsemExpression), ]
# 
# combinedRsemRaw <- hgscRsemExpression
# combinedRsemRawMat <- combinedRsemRaw[, 2:ncol(combinedRsemRaw)]
# combinedRsemRawMat <- apply(combinedRsemRawMat, 2, as.numeric)
# rownames(combinedRsemRawMat) <- combinedRsemRaw$`Hybridization REF`

combinedRsemRawMat <- ov_txn_matrix2
colnames(combinedRsemRawMat) <- gsub("(^.*?-.{3}?)-.*", "\\1", colnames(combinedRsemRawMat))
colnames(combinedRsemRawMat) <- str_replace_all(colnames(combinedRsemRawMat), "01A", "01")
colnames(combinedRsemRawMat) <- str_replace_all(colnames(combinedRsemRawMat), "01B", "01")

### doing it by segment instead of arm from 20231102testingMultipleDiffExp.R
### using deseq - mainly for the shrunk or corrected logFC


hgscSegFile <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/gdac.broadinstitute.org_OV.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt",
                          sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
hgscSegFile$Sample <- gsub("(^.*?-.{3}?)-.*", "\\1",hgscSegFile$Sample)
hgscSegFile$Sample <- substr(hgscSegFile$Sample, 1, nchar(hgscSegFile$Sample) -1)
hgscSegFile <- hgscSegFile[which(hgscSegFile$Sample %in% hgsc_anno2$Sample.ID),]
hgscSegFile$lengthKb <- (hgscSegFile$End - hgscSegFile$Start)/1e4
hgscSegFile2 <- hgscSegFile[which(hgscSegFile$lengthKb > 24), ]

hgsc_arm_allSynTable_amp_good <- hgsc_arm_allSynTable_amp[which(hgsc_arm_allSynTable_amp$sigCheck == "good"), ]
hgsc_arm_allSynTable_del_good <- hgsc_arm_allSynTable_del[which(hgsc_arm_allSynTable_del$sigCheck == "good"), ]
hgsc_arm_allSynTable_amp_good$length <- (hgsc_arm_allSynTable_amp_good$h_end - hgsc_arm_allSynTable_amp_good$h_start)/1e6
hgsc_arm_allSynTable_del_good$length <- (hgsc_arm_allSynTable_del_good$h_end - hgsc_arm_allSynTable_del_good$h_start)/1e6
hgsc_arm_allSynTable_amp_good$block <- paste0("block", 1:nrow(hgsc_arm_allSynTable_amp_good))
hgsc_arm_allSynTable_del_good$block <- paste0("block", 1:nrow(hgsc_arm_allSynTable_del_good))


### find a way to iterate through each significant syntenic region - either 50 indivudal blocks or I combine
### should combine or too many variables for downstream analysis
### nevermind; too much work to combine, need to do individual



i <- unique(hgsc_arm_allSynTable_amp_good$h_chr)[1]
hGoodDiff <- NULL
for (i in unique(hgsc_arm_allSynTable_amp_good$h_chr)) {
  tmpDf <- hgsc_arm_allSynTable_amp_good[which(hgsc_arm_allSynTable_amp_good$h_chr == i), ]
  tmpDiff <- NULL
  if (nrow(tmpDf) > 1) {
    for (j in 2:nrow(tmpDf)) {
      tmpDiff <- c(tmpDiff, tmpDf$h_start[j] - tmpDf$h_end[j-1])
    }
    tmpDiff <- c(0, tmpDiff)
    hGoodDiff <- c(hGoodDiff, tmpDiff)
  } else{
    hGoodDiff <- c(hGoodDiff, 0)
  }
}
hGoodDiff <- hGoodDiff/1e6
hgsc_arm_allSynTable_amp_good$diff <- hGoodDiff

### for each region, which samples have these regions amplified or deleted in the same direction
### do differential expression analysis for those samples against tumor samples without;
### make list of all samples then iterate 

hgscSegFile_amp <- hgscSegFile[which(hgscSegFile$Segment_Mean > 0.2), ]
hgscSegFile_ampGr <- GRanges(seqnames = hgscSegFile_amp$Chromosome,
                             IRanges(start = hgscSegFile_amp$Start,
                                     end = hgscSegFile_amp$End))

allPercentDf <- NULL
listAllSampsWithGain <- list()
for (i in 1:nrow(hgsc_arm_allSynTable_amp_good)) {
  tmpHgGrange <- GRanges(seqnames = hgsc_arm_allSynTable_amp_good$h_chr[i],
                         IRanges(start = hgsc_arm_allSynTable_amp_good$h_start[i],
                                 end = hgsc_arm_allSynTable_amp_good$h_end[i]))
  tmpOverlap <- findOverlaps(tmpHgGrange, hgscSegFile_ampGr)
  tmpIntersect <- pintersect(tmpHgGrange[queryHits(tmpOverlap)],
                             hgscSegFile_ampGr[subjectHits(tmpOverlap)])
  ### getting percent of overlap with the syntenic contig
  tmpPercent <- width(tmpIntersect) / width(tmpHgGrange[queryHits(tmpOverlap)])
  
  tmpAmp <- hgscSegFile_amp[subjectHits(tmpOverlap),]
  tmpAmp$percentOverlap <- tmpPercent
  
  tmpPercentDf <- NULL
  for (j in unique(tmpAmp$Sample)) {
    tmpPercentDf <- rbind(tmpPercentDf, data.frame("sample" = j,
                                                   "percent" = sum(tmpAmp$percentOverlap[which(tmpAmp$Sample == j)])))
  }
  
  tmpKeepSamps <- tmpPercentDf$sample[which(tmpPercentDf$percent > 0.99)]
  listAllSampsWithGain[[i]] <- tmpKeepSamps
  
  
  # quantile(tmpPercentDf$percent, seq(0, 1, 0.01))
  ### from all of the overlaps, give a tiny bit of leeway
  ### only keep things >= 99%
  allPercentDf <- rbind(allPercentDf, tmpPercentDf)
}

quantile(allPercentDf$percent, seq(0, 1, 0.01))
hist(allPercentDf$percent, breaks = 99)
names(listAllSampsWithGain) <- paste0("block_amp_", 1:nrow(hgsc_arm_allSynTable_amp_good))
listAllSampsWithGain2 <- listAllSampsWithGain[-which(unlist(lapply(listAllSampsWithGain, length)) == 0)]
### iterate over rna comparisons. samps with vs without

# combinedRsemRawMat2 <- combinedRsemRawMat[, which(colnames(combinedRsemRawMat) %in% paste0(hgsc_anno2$Sample.ID, "A"))]
combinedRsemRawMat2 <- combinedRsemRawMat[, which(colnames(combinedRsemRawMat) %in% hgsc_anno2$Sample.ID)]



### for deletions
###
###

hgscSegFile_del <- hgscSegFile[which(hgscSegFile$Segment_Mean < -0.2), ]
hgscSegFile_delGr <- GRanges(seqnames = hgscSegFile_del$Chromosome,
                             IRanges(start = hgscSegFile_del$Start,
                                     end = hgscSegFile_del$End))


listAllSampsWithLoss <- list()
for (i in 1:nrow(hgsc_arm_allSynTable_del_good)) {
  tmpHgGrange <- GRanges(seqnames = hgsc_arm_allSynTable_del_good$h_chr[i],
                         IRanges(start = hgsc_arm_allSynTable_del_good$h_start[i],
                                 end = hgsc_arm_allSynTable_del_good$h_end[i]))
  tmpOverlap <- findOverlaps(tmpHgGrange, hgscSegFile_delGr)
  tmpIntersect <- pintersect(tmpHgGrange[queryHits(tmpOverlap)],
                             hgscSegFile_delGr[subjectHits(tmpOverlap)])
  ### getting percent of overlap with the syntenic contig
  tmpPercent <- width(tmpIntersect) / width(tmpHgGrange[queryHits(tmpOverlap)])
  
  tmpDel <- hgscSegFile_del[subjectHits(tmpOverlap),]
  tmpDel$percentOverlap <- tmpPercent
  
  tmpPercentDf <- NULL
  for (j in unique(tmpDel$Sample)) {
    tmpPercentDf <- rbind(tmpPercentDf, data.frame("sample" = j,
                                                   "percent" = sum(tmpDel$percentOverlap[which(tmpDel$Sample == j)])))
  }
  
  tmpKeepSamps <- tmpPercentDf$sample[which(tmpPercentDf$percent > 0.99)]
  listAllSampsWithLoss[[i]] <- tmpKeepSamps
  
}


names(listAllSampsWithLoss) <- paste0("block_loss_", 1:nrow(hgsc_arm_allSynTable_del_good))
listAllSampsWithLoss2 <- listAllSampsWithLoss[-which(unlist(lapply(listAllSampsWithLoss, length)) == 0)]


hg19GeneLocations <- read.table("/mnt/DATA5/tmp/kev/misc/20231018hg19GeneLocations.txt", sep = "\t",
                                header = TRUE, stringsAsFactors = FALSE)

hg19ExonGr <- GRanges(seqnames = hg19GeneLocations$chr, 
                      IRanges(start = hg19GeneLocations$start, hg19GeneLocations$end))
affectedGenesGains <- NULL
for (i in unique(hgsc_arm_allSynTable_amp_good$block)) {
  tmp <- hgsc_arm_allSynTable_amp_good[which(hgsc_arm_allSynTable_amp_good$block == i),]
  tmpGr <- GRanges(seqnames = tmp$h_chr, IRanges(start = tmp$h_start, end = tmp$h_end))
  tmpGene <- hg19GeneLocations[queryHits(findOverlaps(hg19ExonGr, tmpGr)),]
  affectedGenesGains <- rbind(affectedGenesGains, data.frame(tmpGene, "type" = rep(paste0("amp_", i), nrow(tmpGene))))
}


affectedGenesLosses <- NULL
for (i in unique(hgsc_arm_allSynTable_del_good$block)) {
  tmp <- hgsc_arm_allSynTable_del_good[which(hgsc_arm_allSynTable_del_good$block == i),]
  tmpGr <- GRanges(seqnames = tmp$h_chr, IRanges(start = tmp$h_start, end = tmp$h_end))
  tmpGene <- hg19GeneLocations[queryHits(findOverlaps(hg19ExonGr, tmpGr)),]
  affectedGenesLosses <- rbind(affectedGenesLosses, data.frame(tmpGene, "type" = rep(paste0("del_", i), nrow(tmpGene))))
}

affectedGenesDf <- rbind(affectedGenesGains, affectedGenesLosses)
affectedGenesDf$type2 <- str_replace_all(affectedGenesDf$type, "amp_block", "block_amp_")
affectedGenesDf$type2 <- str_replace_all(affectedGenesDf$type2, "del_block", "block_loss_")
### using either cox-proportional for all synteny blocks


library(TCGAbiolinks)
library(survminer)
library(survival)

clinical_hgsc <- GDCquery_clinic("TCGA-OV")

clinical_hgsc2 <- clinical_hgsc[, which(colnames(clinical_hgsc) %in% c("submitter_id", "vital_status", "days_to_last_follow_up", "days_to_death"))]
clinical_hgsc2$deceased <- ifelse(clinical_hgsc2$vital_status == "Alive", FALSE, TRUE)
clinical_hgsc2$overall_survival <- ifelse(clinical_hgsc2$vital_status == "Alive",
                                             clinical_hgsc2$days_to_last_follow_up,
                                             clinical_hgsc2$days_to_death)
clinical_hgsc2$submitter_id2 <- paste0(clinical_hgsc2$submitter_id, "-01")
clinical_hgsc2 <- clinical_hgsc2[match(colnames(combinedRsemRawMat2), paste0(clinical_hgsc2$submitter_id2, "A")),]


### trying xena's TCGA clinical data provides precalculated OS.time. DSS and PFI
tmpRsemColnames <- gsub("\\-01", "", colnames(combinedRsemRawMat2))
clinical_hgsc <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/xena/survivalOV_survival.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)
clinical_hgsc2 <- clinical_hgsc[match(tmpRsemColnames, clinical_hgsc$X_PATIENT),]
colnames(clinical_hgsc2)[1] <- "submitter_id2"
colnames(clinical_hgsc2) <- str_replace_all(colnames(clinical_hgsc2), "\\.", "_")
clinical_hgsc2$PFI_status <- ifelse(clinical_hgsc2$PFI == 1, TRUE, FALSE)

### doing it by arm status
hgscArmStatus <- cbind(hg19ArmLocations$str, broadBySampleGisticHgsc2[, 4:ncol(broadBySampleGisticHgsc2)])
colnames(hgscArmStatus)[1] <- "Chromosome Arm"
colnames(hgscArmStatus)[2:ncol(hgscArmStatus)] <- gsub("(^.*?-.{3}?)-.*", "\\1", colnames(hgscArmStatus)[2:ncol(hgscArmStatus)])
colnames(hgscArmStatus)[2:ncol(hgscArmStatus)] <- str_replace_all(colnames(hgscArmStatus)[2:ncol(hgscArmStatus)], "01A", "01")
hgscArmStatusMat <- hgscArmStatus[, 2:ncol(hgscArmStatus)]
rownames(hgscArmStatusMat) <- hgscArmStatus$`Chromosome Arm`



chromGainStatus <- list()
chromLossStatus <- list()
for (i in 1:nrow(hgscArmStatusMat)) {
  tmpGain <- colnames(hgscArmStatusMat)[which(hgscArmStatusMat[i, ] > 0.2)]
  tmpLoss <- colnames(hgscArmStatusMat)[which(hgscArmStatusMat[i, ] < -0.2)]
  
  chromGainStatus[[i]] <- tmpGain
  chromLossStatus[[i]] <- tmpLoss
  
}

chromGainStatusDf <- NULL
for (i in seq_along(chromGainStatus)) {
  tmpNames <- chromGainStatus[[i]]
  tmpVector <- ifelse(clinical_hgsc2$submitter_id2 %in% tmpNames, "Yes", "No")
  chromGainStatusDf <- cbind(chromGainStatusDf, tmpVector)
}

colnames(chromGainStatusDf) <- paste0("gain_", rownames(hgscArmStatusMat))


chromLossStatusDf <- NULL
for (i in seq_along(chromLossStatus)) {
  tmpNames <- chromLossStatus[[i]]
  tmpVector <- ifelse(clinical_hgsc2$submitter_id2 %in% tmpNames, "Yes", "No")
  chromLossStatusDf <- cbind(chromLossStatusDf, tmpVector)
}

colnames(chromLossStatusDf) <- paste0("loss_", rownames(hgscArmStatusMat))

clinical_hgsc_arm <- cbind(clinical_hgsc2[, c("submitter_id2", "X_PATIENT", "OS", "OS_time")], chromGainStatusDf, chromLossStatusDf)
# clinical_hgsc_arm <- cbind(clinical_hgsc2[, c("submitter_id2", "X_PATIENT", "DSS", "DSS_time")], chromGainStatusDf, chromLossStatusDf)
# clinical_hgsc_arm <- cbind(clinical_hgsc2[, c("submitter_id2", "X_PATIENT", "PFI", "PFI_time")], chromGainStatusDf, chromLossStatusDf)

### for os

### old with GDCquery
### new with xena browser data b/c they have DSS and PFI

# survivalarmsDf_arm <- NULL
# for (i in 8:ncol(clinical_hgsc_arm)) {
#   tmpFit <- survdiff(Surv(overall_survival, deceased) ~ clinical_hgsc_arm[,i], data = clinical_hgsc_arm)
#   tmpPVal <- pchisq(tmpFit$chisq, length(tmpFit$n)-1, lower.tail = FALSE)
#   survivalarmsDf_arm <- rbind(survivalarmsDf_arm, data.frame("var" = colnames(clinical_hgsc_arm)[i], "pval" = tmpPVal))
# }

# survivalarmsDf_arm <- NULL
# for (i in 8:ncol(clinical_hgsc_arm)) {
#   tmpFit <- survdiff(Surv(OS_time, OS) ~ clinical_hgsc_arm[,i], data = clinical_hgsc_arm)
#   tmpPVal <- pchisq(tmpFit$chisq, length(tmpFit$n)-1, lower.tail = FALSE)
#   survivalarmsDf_arm <- rbind(survivalarmsDf_arm, data.frame("var" = colnames(clinical_hgsc_arm)[i], "pval" = tmpPVal))
# }


survivalarmsDf_arm <- NULL
for (i in 5:ncol(clinical_hgsc_arm)) {
  tmpFit <- survdiff(Surv(OS_time, OS) ~ clinical_hgsc_arm[,i], data = clinical_hgsc_arm)
  tmpPVal <- pchisq(tmpFit$chisq, length(tmpFit$n)-1, lower.tail = FALSE)
  survivalarmsDf_arm <- rbind(survivalarmsDf_arm, data.frame("var" = colnames(clinical_hgsc_arm)[i], "pval" = tmpPVal))
}


# survivalarmsDf_arm <- NULL
# for (i in 8:ncol(clinical_hgsc_arm)) {
#   tmpFit <- survdiff(Surv(DSS_time, DSS) ~ clinical_hgsc_arm[,i], data = clinical_hgsc_arm)
#   tmpPVal <- pchisq(tmpFit$chisq, length(tmpFit$n)-1, lower.tail = FALSE)
#   survivalarmsDf_arm <- rbind(survivalarmsDf_arm, data.frame("var" = colnames(clinical_hgsc_arm)[i], "pval" = tmpPVal))
# }
# 
# 
# survivalarmsDf_arm <- NULL
# for (i in 8:ncol(clinical_hgsc_arm)) {
#   tmpFit <- survdiff(Surv(PFI_time, PFI) ~ clinical_hgsc_arm[,i], data = clinical_hgsc_arm)
#   tmpPVal <- pchisq(tmpFit$chisq, length(tmpFit$n)-1, lower.tail = FALSE)
#   survivalarmsDf_arm <- rbind(survivalarmsDf_arm, data.frame("var" = colnames(clinical_hgsc_arm)[i], "pval" = tmpPVal))
# }



survivalarmsDf_arm$adj.pval <- p.adjust(survivalarmsDf_arm$pval, method = "BH")

# fit <- survfit(Surv(overall_survival, deceased) ~ loss_16p, data = clinical_hgsc_arm)
fit <- survfit(Surv(OS_time, OS) ~ loss_16p, data = clinical_hgsc_arm)
ggsurvplot(fit,
           data = clinical_hgsc_arm,
           pval = T,
           risk.table = T)


### for blocks


blockStatusLossDf <- NULL
for (i in seq_along(names(listAllSampsWithLoss2))) {
  tmpNames <- listAllSampsWithLoss2[[i]]
  tmpVector <- ifelse(clinical_hgsc2$submitter_id2 %in% tmpNames, "Yes", "No")
  blockStatusLossDf <- cbind(blockStatusLossDf, tmpVector)
}
colnames(blockStatusLossDf) <- names(listAllSampsWithLoss2)

blockStatusGainDf <- NULL
for (i in seq_along(names(listAllSampsWithGain2))) {
  tmpNames <- listAllSampsWithGain2[[i]]
  tmpVector <- ifelse(clinical_hgsc2$submitter_id2 %in% tmpNames, "Yes", "No")
  blockStatusGainDf <- cbind(blockStatusGainDf, tmpVector)
}
colnames(blockStatusGainDf) <- names(listAllSampsWithGain2)

clinical_hgsc_block <- cbind(clinical_hgsc2[, c("submitter_id2", "X_PATIENT", "OS", "OS_time")], blockStatusGainDf, blockStatusLossDf)
# clinical_hgsc_block <- cbind(clinical_hgsc2[, c("submitter_id2", "X_PATIENT", "DSS", "DSS_time")], blockStatusGainDf, blockStatusLossDf)
# clinical_hgsc_block <- cbind(clinical_hgsc2[, c("submitter_id2", "X_PATIENT", "PFI", "PFI_time")], blockStatusGainDf, blockStatusLossDf)



# fit <- survfit(Surv(overall_survival, deceased) ~ block_amp_20, data = clinical_hgsc2)
# pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20231109logRankBlockAmp20", ".pdf"),
#     useDingbats = FALSE, width = 10, height = 5)
# ggsurvplot(fit,
#            data = clinical_coadread,
#            pval = T,
#            risk.table = T)
# dev.off()
### iterating for all arms

survivalBlocksDf_block <- NULL
for (i in 5:ncol(clinical_hgsc_block)) {
  tmpFit <- survdiff(Surv(OS_time, OS) ~ clinical_hgsc_block[,i], data = clinical_hgsc_block)
  tmpPVal <- pchisq(tmpFit$chisq, length(tmpFit$n)-1, lower.tail = FALSE)
  survivalBlocksDf_block <- rbind(survivalBlocksDf_block, data.frame("var" = colnames(clinical_hgsc_block)[i], "pval" = tmpPVal))
}

fit <- survfit(Surv(OS_time, OS) ~ block_amp_26, data = clinical_hgsc_block)
block_cox <- coxph(Surv(OS_time, OS) ~ block_amp_26, data = clinical_hgsc_block)
exp(coef(block_cox)) 
exp(confint(block_cox))

pdf("/mnt/DATA5/tmp/kev/misc/20240107hgscBlock26Coxph.pdf")
ggsurvplot(fit,
           data = clinical_hgsc_block,
           pval = T,
           risk.table = T)
dev.off()

# survivalBlocksDf_block <- NULL
# for (i in 5:ncol(clinical_hgsc_block)) {
#   tmpFit <- survdiff(Surv(DSS_time, DSS) ~ clinical_hgsc_block[,i], data = clinical_hgsc_block)
#   tmpPVal <- pchisq(tmpFit$chisq, length(tmpFit$n)-1, lower.tail = FALSE)
#   survivalBlocksDf_block <- rbind(survivalBlocksDf_block, data.frame("var" = colnames(clinical_hgsc_block)[i], "pval" = tmpPVal))
# }
# 
# 
# survivalBlocksDf_block <- NULL
# for (i in 5:ncol(clinical_hgsc_block)) {
#   tmpFit <- survdiff(Surv(PFI_time, PFI) ~ clinical_hgsc_block[,i], data = clinical_hgsc_block)
#   tmpPVal <- pchisq(tmpFit$chisq, length(tmpFit$n)-1, lower.tail = FALSE)
#   survivalBlocksDf_block <- rbind(survivalBlocksDf_block, data.frame("var" = colnames(clinical_hgsc_block)[i], "pval" = tmpPVal))
# }


survivalBlocksDf_block$adj.pval <- p.adjust(survivalBlocksDf_block$pval, method = "BH")


### similar to above, but instead of edgeR use DEseq2

library(foreach)
library(doParallel)

cl <- makeCluster(20)
registerDoParallel(cl)
start <- Sys.time()
diffArmAmpDeSeq <- foreach(i=seq_along(listAllSampsWithGain2),.combine = 'comb', .multicombine = TRUE,
                           .init = list(list()), .packages = c("DESeq2", "stringr")) %dopar% {
                             
                             tmpHgscArmNames <- unlist(listAllSampsWithGain2[i])
                             tmpMat <- combinedRsemRawMat2
                             colnames(tmpMat) <- substr(colnames(tmpMat), 1, nchar(colnames(tmpMat)) -1)
                             tmpRawCounts <- tmpMat
                             tmpGroup <- rep(1, ncol(tmpRawCounts))
                             tmpGroup[which(colnames(tmpRawCounts) %in% tmpHgscArmNames)] <- 2
                             
                             tmpColData <- data.frame("block" = factor(tmpGroup))
                             rownames(tmpColData) <- colnames(tmpRawCounts)
                             
                             dds <-  DESeqDataSetFromMatrix(countData = round(tmpRawCounts),
                                                            colData = tmpColData,
                                                            design = ~block)
                             smallestGroupSize <- 3
                             keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
                             dds <- dds[keep,]
                             
                             dds$block <- relevel(dds$block, ref = "1")
                             
                             dds <- DESeq(dds)
                             # res <- results(dds)
                             res2 <- lfcShrink(dds, coef = "block_2_vs_1",type="apeglm")
                             
                             # resTable <- data.frame(cbind("genes" = rownames(res), 
                             #                   res))
                             resTable2 <- data.frame(cbind("genes" = rownames(res2), 
                                                           res2))
                             # or to shrink log fold changes association with condition:
                             
                             resTable2$block <- names(listAllSampsWithGain2)[i]
                             
                             return(list(resTable2))
                           }


stopCluster(cl)
print( Sys.time() - start )


diffArmAmpDeSeqDf <- do.call(rbind, diffArmAmpDeSeq[[1]])
diffArmAmpDeSeqDf$genes <- str_remove(diffArmAmpDeSeqDf$genes, "\\|.*")

# write.table(diffArmAmpDeSeqDf, "/mnt/DATA5/tmp/kev/misc/20240105HgscDiffArmAmpDeSeq.txt",sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)



### same DEseq but for losses


cl <- makeCluster(20)
registerDoParallel(cl)
start <- Sys.time()
diffArmDelDeSeq <- foreach(i=seq_along(listAllSampsWithLoss2),.combine = 'comb', .multicombine = TRUE,
                           .init = list(list()), .packages = c("DESeq2", "stringr")) %dopar% {
                             
                             tmpHgscArmNames <- unlist(listAllSampsWithLoss2[i])
                             tmpMat <- combinedRsemRawMat2
                             colnames(tmpMat) <- substr(colnames(tmpMat), 1, nchar(colnames(tmpMat)) -1)
                             tmpRawCounts <- tmpMat
                             tmpGroup <- rep(1, ncol(tmpRawCounts))
                             tmpGroup[which(colnames(tmpRawCounts) %in% tmpHgscArmNames)] <- 2
                             
                             tmpColData <- data.frame("block" = factor(tmpGroup))
                             rownames(tmpColData) <- colnames(tmpRawCounts)
                             
                             dds <-  DESeqDataSetFromMatrix(countData = round(tmpRawCounts),
                                                            colData = tmpColData,
                                                            design = ~block)
                             smallestGroupSize <- 3
                             keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
                             dds <- dds[keep,]
                             
                             dds$block <- relevel(dds$block, ref = "1")
                             
                             dds <- DESeq(dds)
                             # res <- results(dds)
                             res2 <- lfcShrink(dds, coef = "block_2_vs_1",type="apeglm")
                             
                             # resTable <- data.frame(cbind("genes" = rownames(res), 
                             #                   res))
                             resTable2 <- data.frame(cbind("genes" = rownames(res2), 
                                                           res2))
                             # or to shrink log fold changes association with condition:
                             
                             resTable2$block <- names(listAllSampsWithLoss2)[i]
                             
                             return(list(resTable2))
                           }


stopCluster(cl)
print( Sys.time() - start )


diffArmDelDeSeqDf <- do.call(rbind, diffArmDelDeSeq[[1]])
diffArmDelDeSeqDf$genes <- str_remove(diffArmDelDeSeqDf$genes, "\\|.*")

# write.table(diffArmDelDeSeqDf, "/mnt/DATA5/tmp/kev/misc/20240105HgscDiffArmDelDeSeq.txt",sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)




### would make things easier if I match blocks with their chromosomes


diffArmAmpDeSeqDf$blockStripped <- str_remove(diffArmAmpDeSeqDf$block, "_amp_")
diffArmAmpDeSeqDf$h_chr <- hgsc_arm_allSynTable_amp_good$h_chr[match(diffArmAmpDeSeqDf$blockStripped, hgsc_arm_allSynTable_amp_good$block)]

diffArmDelDeSeqDf$blockStripped <- str_remove(diffArmDelDeSeqDf$block, "_loss_")
diffArmDelDeSeqDf$h_chr <- hgsc_arm_allSynTable_del_good$h_chr[match(diffArmDelDeSeqDf$blockStripped, hgsc_arm_allSynTable_del_good$block)]


allDiffArmDeSeqDf <- rbind(diffArmDelDeSeqDf, diffArmAmpDeSeqDf)
i <- unique(allDiffArmDeSeqDf$block)[51]

nominatedGenesDf <- NULL
for (i in unique(allDiffArmDeSeqDf$block)) {
  tmp <- allDiffArmDeSeqDf[which(allDiffArmDeSeqDf$block == i), ]
  if (grepl("amp", i)) {
    tmp$color <- "#000000"
    # tmp$color[which(tmp$padj < 0.05 & tmp$log2FoldChange > (log2(1.5 * 0.9)))] <- "#8b0000"
    tmp$color[which(tmp$padj < 0.05 & tmp$log2FoldChange > (log2(1.5)))] <- "#8b0000"
    tmp$on_block <- "no"
    tmp$on_block[which(tmp$genes %in% affectedGenesDf$gene[which(affectedGenesDf$type2 == i)])] <- "yes"
    # tmp2 <- tmp[which(tmp$on_block  == "yes"), ]
    tmp2 <- tmp[which(tmp$color != "#000000" & tmp$on_block  == "yes"), ]
    nominatedGenesDf <- rbind(nominatedGenesDf, tmp2)
  } else if(grepl("loss", i)){
    tmp$color <- "#000000"
    # tmp$color[which(tmp$padj < 0.05 & tmp$log2FoldChange < (log2(0.5/0.9)))] <- "#00008B"
    tmp$color[which(tmp$padj < 0.05 & tmp$log2FoldChange < (log2(0.5)))] <- "#00008B"
    tmp$on_block <- "no"
    tmp$on_block[which(tmp$genes %in% affectedGenesDf$gene[which(affectedGenesDf$type2 == i)])] <- "yes"
    tmp2 <- tmp[which(tmp$color != "#000000" & tmp$on_block  == "yes"), ]
    # tmp2 <- tmp[which(tmp$on_block  == "yes"), ]
    nominatedGenesDf <- rbind(nominatedGenesDf, tmp2)
  }
}

### depmap overlapping genes

depmapTableRnai <- depmap::depmap_rnai()
depmapTableRnai_ovary <- depmapTableRnai[grep("OVARY", depmapTableRnai$cell_line),]
depmapTableCrispr <- depmap::depmap_crispr()
depmapTableCrispr_ovary <- depmapTableCrispr[grep("OVARY", depmapTableCrispr$cell_line),]


depmapRnai_nom <- depmapTableRnai_ovary[which(depmapTableRnai_ovary$gene_name %in% nominatedGenesDf$genes), ]
depmapCripspr_nom <- depmapTableCrispr_ovary[which(depmapTableCrispr_ovary$gene_name %in% nominatedGenesDf$genes), ]
depmapRnai_nom2 <- depmapRnai_nom[which(depmapRnai_nom$gene_name %in% c("MED30", "PHB2", "USP5", "MCM2")), ]
depmapCripspr_nom2 <- depmapCripspr_nom[which(depmapCripspr_nom$gene_name %in% c("MED30", "PHB2", "USP5", "MCM2")), ]
depmapRnai_nom2$type <- "RNAi"
depmapCripspr_nom2$type <- "CRISPR"
depmapCombined <- rbind(depmapRnai_nom2, depmapCripspr_nom2)

pdf("/mnt/DATA5/tmp/kev/misc/20240106hgscBoxplotDepmap.pdf", width = 7, height = 5)
ggplot(depmapCombined) + geom_violin(aes(x = gene_name, y = dependency, color = type)) +
  scale_color_manual(values = c("darkblue", "purple2")) + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

### for whatever reason no normals can be found for ovarian cancer

combinedRsemRawMat3 <- combinedRsemRawMat2
rownames(combinedRsemRawMat3) <- str_remove(rownames(combinedRsemRawMat3), "\\|.*")
nominatedGenesZScore <- combinedRsemRawMat3[which(rownames(combinedRsemRawMat3) %in% c("MED30", "PHB2", "USP5", "MCM2")),]
nominatedGenesZScore2 <- data.frame(apply(t(nominatedGenesZScore), 2, function(x) (x - mean(x))/sd(x)))

### don't know if i should make 3 groups >1 zscore, 0 and <1 zscore, or just high vs everything 
mcm2Status <- ifelse(nominatedGenesZScore2$MCM2 > 0.99, "High", "Low")
med30Status <- ifelse(nominatedGenesZScore2$MED30 > 0.99, "High", "Low")
phb2Status <- ifelse(nominatedGenesZScore2$PHB2 > 0.99, "High", "Low")
usp5Status <- ifelse(nominatedGenesZScore2$USP5 > 0.99, "High", "Low")


clinical_hgsc_gene <- clinical_hgsc2[, c("submitter_id2", "X_PATIENT", "OS", "OS_time")]
clinical_hgsc_gene <- clinical_hgsc_gene[which(clinical_hgsc_gene$X_PATIENT %in% str_remove(rownames(nominatedGenesZScore2), "\\-01")), ]
clinical_hgsc_gene <- clinical_hgsc_gene[match(str_remove(rownames(nominatedGenesZScore2), "\\-01"), clinical_hgsc_gene$X_PATIENT), ]
clinical_hgsc_gene$mcm2 <- mcm2Status
clinical_hgsc_gene$med30 <-med30Status
clinical_hgsc_gene$phb2 <- phb2Status
clinical_hgsc_gene$usp5 <- usp5Status

survivalarmsDf_gene <- NULL
for (i in 5:ncol(clinical_hgsc_gene)) {
  tmpFit <- survdiff(Surv(OS_time, OS) ~ clinical_hgsc_gene[,i], data = clinical_hgsc_gene)
  tmpPVal <- pchisq(tmpFit$chisq, length(tmpFit$n)-1, lower.tail = FALSE)
  survivalarmsDf_gene <- rbind(survivalarmsDf_gene, data.frame("var" = colnames(clinical_hgsc_gene)[i], "pval" = tmpPVal))
}


### volcano plots
###
###

lrtTab_bl10 <- allDiffArmDeSeqDf[which(allDiffArmDeSeqDf$block == "block_amp_10"),]
lrtTab_bl10$log10Q <- -1 * log10(lrtTab_bl10$padj)
lrtTab_bl10$color <- "#000000"
lrtTab_bl10$color[which(lrtTab_bl10$padj < 0.05 & lrtTab_bl10$log2FoldChange > (log2(1.5)))] <- "#8b0000"
lrtTab_bl10$color[which(lrtTab_bl10$padj < 0.05 & lrtTab_bl10$log2FoldChange < (log2(0.5)))] <- "#00008B"
lrtTab_bl10$on_block <- "no"
lrtTab_bl10$on_block[which(lrtTab_bl10$genes %in% affectedGenesDf$gene[which(affectedGenesDf$type == "amp_block10")])] <- "yes"
lrtTab_bl10$color2 <- lrtTab_bl10$color
lrtTab_bl10$color2[which(lrtTab_bl10$on_block == "no")] <- "#000000"
lrtTab_bl10$color2[which(lrtTab_bl10$color == "#00008B")] <- "#000000"
lrtTab_bl10$alpha <- 1
lrtTab_bl10$alpha[which(lrtTab_bl10$color2 == "#000000")] <- 0.1



lrtTab_bl26 <- allDiffArmDeSeqDf[which(allDiffArmDeSeqDf$block == "block_amp_26"),]
lrtTab_bl26$log10Q <- -1 * log10(lrtTab_bl26$padj)
lrtTab_bl26$color <- "#000000"
lrtTab_bl26$color[which(lrtTab_bl26$padj < 0.05 & lrtTab_bl26$log2FoldChange > (log2(1.5)))] <- "#8b0000"
lrtTab_bl26$color[which(lrtTab_bl26$padj < 0.05 & lrtTab_bl26$log2FoldChange < (log2(0.5)))] <- "#00008B"
lrtTab_bl26$on_block <- "no"
lrtTab_bl26$on_block[which(lrtTab_bl26$genes %in% affectedGenesDf$gene[which(affectedGenesDf$type == "amp_block26")])] <- "yes"
lrtTab_bl26$color2 <- lrtTab_bl26$color
lrtTab_bl26$color2[which(lrtTab_bl26$on_block == "no")] <- "#000000"
lrtTab_bl26$color2[which(lrtTab_bl26$color == "#00008B")] <- "#000000"
lrtTab_bl26$alpha <- 1
lrtTab_bl26$alpha[which(lrtTab_bl26$color2 == "#000000")] <- 0.1


lrtTab_bl38 <- allDiffArmDeSeqDf[which(allDiffArmDeSeqDf$block == "block_amp_38"),]
lrtTab_bl38$log10Q <- -1 * log10(lrtTab_bl38$padj)
lrtTab_bl38$color <- "#000000"
lrtTab_bl38$color[which(lrtTab_bl38$padj < 0.05 & lrtTab_bl38$log2FoldChange > (log2(1.5)))] <- "#8b0000"
lrtTab_bl38$color[which(lrtTab_bl38$padj < 0.05 & lrtTab_bl38$log2FoldChange < (log2(0.5)))] <- "#00008B"
lrtTab_bl38$on_block <- "no"
lrtTab_bl38$on_block[which(lrtTab_bl38$genes %in% affectedGenesDf$gene[which(affectedGenesDf$type == "amp_block38")])] <- "yes"
lrtTab_bl38$color2 <- lrtTab_bl38$color
lrtTab_bl38$color2[which(lrtTab_bl38$on_block == "no")] <- "#000000"
lrtTab_bl38$color2[which(lrtTab_bl38$color == "#00008B")] <- "#000000"
lrtTab_bl38$alpha <- 1
lrtTab_bl38$alpha[which(lrtTab_bl38$color2 == "#000000")] <- 0.1



### the pathway analysis was done in brahm since avatar's 20230102 R version was too old
###
###



library(clusterProfiler)
library(msigdbr)
packageVersion("clusterProfiler")
library(gridExtra)
library(ggplot2)
library(ggpubr)

hg19GeneLocations <- read.table("/avatar_data5/tmp/kev/misc/20231018hg19GeneLocations.txt", sep = "\t",
                                header = TRUE, stringsAsFactors = FALSE)

msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
msig_hallmark<- msig_hallmark[, c("gs_name", "gene_symbol")]

### using DeSeq results

diffArmAmp <- read.table("/avatar_data5/tmp/kev/misc/20231109diffArmAmpDeSeq.txt", sep = "\t",
                         header = TRUE, stringsAsFactors = FALSE)

i <- unique(diffArmAmp$block)[1]
allGainsGseaDf <- NULL
for (i in unique(diffArmAmp$block)) {
  tmpDiffDf <- diffArmAmp[which(diffArmAmp$block == i), ]
  
  if (length(which(tmpDiffDf$genes == "?")) > 0) {
    tmpDiffDf <- tmpDiffDf[-which(tmpDiffDf$genes == "?"),]
  }
  
  if (length(which(duplicated(tmpDiffDf$genes))) > 0) {
    tmpDiffDf <- tmpDiffDf[-which(duplicated(tmpDiffDf$genes)),]
  }
  
  
  # gene_list <- -log10(tmpDiffDf$PValue) * sign(tmpDiffDf$logFC)
  gene_list <- sign(tmpDiffDf$log2FoldChange) * -log10(tmpDiffDf$pvalue)
  names(gene_list) <- tmpDiffDf$genes
  gene_list = sort(gene_list, decreasing = TRUE)
  
  with1 <- GSEA(gene_list, TERM2GENE = msig_hallmark, nPermSimple = 10000)
  if(nrow(with1@result) == 0){
    next()
  }
  # with1 <- GSEA(gene_list, TERM2GENE = msig_hallmark, pvalueCutoff = 0.05, eps = 0)
  with1@result$IDsign <- paste(with1@result$ID, sign(with1@result$NES))
  with1@result$count <- unname(sapply(with1@result$core_enrichment,
                                      function(x) length(unlist(strsplit(x, "/")))))
  with1@result$geneRatio <- unname(sapply(with1@result$core_enrichment,
                                          function(x) length(unlist(strsplit(x, "/")))))/with1@result$setSize
  res <- with1@result
  res$block <- i
  allGainsGseaDf <- rbind(allGainsGseaDf, res)
}


diffArmDel <- read.table("/avatar_data5/tmp/kev/misc/20231109diffArmDelDeSeq.txt", sep = "\t",
                         header = TRUE, stringsAsFactors = FALSE)


i <- unique(diffArmDel$block)[24]
allLossGseaDf <- NULL
for (i in unique(diffArmDel$block)) {
  tmpDiffDf <- diffArmDel[which(diffArmDel$block == i), ]
  
  if (length(which(tmpDiffDf$genes == "?")) > 0) {
    tmpDiffDf <- tmpDiffDf[-which(tmpDiffDf$genes == "?"),]
  }
  
  if (length(which(duplicated(tmpDiffDf$genes))) > 0) {
    tmpDiffDf <- tmpDiffDf[-which(duplicated(tmpDiffDf$genes)),]
  }
  
  
  # gene_list <- -log10(tmpDiffDf$PValue) * sign(tmpDiffDf$logFC)
  gene_list <- sign(tmpDiffDf$log2FoldChange) * -log10(tmpDiffDf$pvalue)
  names(gene_list) <- tmpDiffDf$genes
  gene_list = sort(gene_list, decreasing = TRUE)
  
  with1 <- GSEA(gene_list, TERM2GENE = msig_hallmark, nPermSimple = 10000)
  # with1 <- GSEA(gene_list, TERM2GENE = msig_hallmark, pvalueCutoff = 0.05, eps = 0)
  if(nrow(with1@result) == 0){
    next()
  }
  with1@result$IDsign <- paste(with1@result$ID, sign(with1@result$NES))
  with1@result$count <- unname(sapply(with1@result$core_enrichment,
                                      function(x) length(unlist(strsplit(x, "/")))))
  with1@result$geneRatio <- unname(sapply(with1@result$core_enrichment,
                                          function(x) length(unlist(strsplit(x, "/")))))/with1@result$setSize
  res <- with1@result
  res$block <- i
  allLossGseaDf <- rbind(allLossGseaDf, res)
}


allBlocksGsea <- rbind(allGainsGseaDf, allLossGseaDf)
allBlocksGsea$str <- paste(allBlocksGsea$ID, allBlocksGsea$block)
allBlocksGsea$str2 <- paste(allBlocksGsea$IDsign, allBlocksGsea$block)
### 10 blocks with unique alterations
table(allBlocksGsea$IDsign)[order(table(allBlocksGsea$IDsign))]



### will do figure based on order
###
###

### freqplots
library(gridExtra)

e <- freqPlotv2(armGisticHgsc_amp_bed, armGisticHgsc_del_bed,
                main = "hgsc n = 356", speciesType = "human")
f <- freqPlotv2(allMouseAneu_hgsc_amp_bed, allMouseAneu_hgsc_del_bed,
                main = "mouse hgsoc n = 84", speciesType = "mouse")

grid.arrange(e, f, nrow = 2)

pdf("/mnt/DATA5/tmp/kev/misc/20231117hgscAneuploidyFreqPlot.pdf", useDingbats = FALSE,
    width = 12, height = 10)
grid.arrange(e, f, nrow = 2)
dev.off()


### creating a different heatmap for each type i.e met, adenoma and adenocarcinoma
library(pheatmap)
annoTableHgsc2 <- annoTableHgsc[which(annoTableHgsc$Type %in% c("MMMT", "HGSC")), ]
allMouseAneuploidy2 <- allMouseAneuploidy
rownames(allMouseAneuploidy2) <- str_remove(rownames(allMouseAneuploidy2), "x.*")
annoCol2 <- list("Genotype" = c("BPRN" = "darkred", "BPN" = "darkorange",
                                "BPP" = "purple", "BPR" = "darkgreen", "UNK" = "grey"),
                 "Type" = c("MMMT" = "darkblue", "HGSC" = "darkred"))

allMouseAneuploidyHgsc <- allMouseAneuploidy2[which(rownames(allMouseAneuploidy2) %in% annoTableHgsc2$Sample2),]
allMouseAneuploidyHgsc[which(allMouseAneuploidyHgsc == "none")] <- 0
allMouseAneuploidyHgsc[which(allMouseAneuploidyHgsc == "gain")] <- 1
allMouseAneuploidyHgsc[which(allMouseAneuploidyHgsc == "loss")] <- -1
allMouseAneuploidyHgscMat <- apply(allMouseAneuploidyHgsc, 2, as.numeric)
rownames(allMouseAneuploidyHgscMat ) <- rownames(allMouseAneuploidyHgsc)
allMouseAneuploidyHgscMat2 <- allMouseAneuploidyHgscMat[which(rownames(allMouseAneuploidyHgscMat) %in% goodTcHgsc2), ]

heatmapAnnoDfHgsc <- data.frame(annoTableHgsc2[, c("Type", "Geno2")])
rownames(heatmapAnnoDfHgsc) <- annoTableHgsc2$Sample2
colnames(heatmapAnnoDfHgsc) <- c("Type", "Genotype")
heatmapAnnoDfHgsc2 <- heatmapAnnoDfHgsc[which(rownames(heatmapAnnoDfHgsc) %in% rownames(allMouseAneuploidyHgscMat2)),]


heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-1,1,2/1000)


pdf(paste0("/mnt/DATA5/tmp/kev/misc/20231116mm10HgscHeatmap", ".pdf"), useDingbats = FALSE, width = 10)
pheatmap(allMouseAneuploidyHgscMat2, cluster_rows = TRUE, cluster_cols = FALSE, cellwidth = 20,
         cellheight = 5, fontsize_row = 5, color = heatMapCol, breaks = colors.breaks,
         annotation_colors = annoCol2, annotation_row = heatmapAnnoDfHgsc2, main = "Brca/Tp53/Kras/Nf1 mice")
dev.off()


### fraction of aneuploidy altered

### need to do parallel loop to calculate all fga for graph

dir <- "/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/"
setwd(dir)
allDirs <- list.dirs()
allDirs <- allDirs[grep("TP", allDirs)]

listOfFiles <- allDirs

allCans <- c("BLCA", "BRCA", "CESC", "COADREAD", "ESCA", "GBM",
             "GBMLGG", "HNSC", "KIPAN", "KIRC", "LGG", "LIHC",
             "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD",
             "SARC", "STES", "TGCT", "THYM")


allCancerFga <- NULL
for (j in seq_along(listOfFiles)) {
  setwd(dir)
  setwd(listOfFiles[j])
  
  tmpBroad <- read.table("broad_values_by_arm.txt", sep = "\t", stringsAsFactors = FALSE,
                         header = TRUE, check.names = FALSE)
  colnames(tmpBroad)[2:ncol(tmpBroad)] <- gsub("(^.*?-.{3}?)-.*", "\\1",
                                               colnames(tmpBroad)[2:ncol(tmpBroad)])
  if (allCans[j] == "OV") {
    tmpBroad <- tmpBroad[, c(1, grep(paste0(hgsc_anno2$Sample.ID, collapse = "|"), colnames(tmpBroad)))] 
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
  tmpRes <- fgaCalculator_tcga_aneuploidy(tmpBroadMelt2)
  tmpRes$Cancer <- allCans[j]
  allCancerFga <- rbind(allCancerFga, tmpRes)
}


hgscFga <- fgaCalculator_amp_aneuploidy(allMouseAneu_hgsc2)
hgscFga$Cancer <- "mm_hgsc"

allFga <- rbind(allCancerFga, hgscFga)
allFga$Cancer <- reorder(allFga$Cancer, allFga$fga, median)
allFga$Cancer <- factor(allFga$Cancer, rev(levels(allFga$Cancer)))
fgaColorVector <- rep("#000000", length(unique(allFga$Cancer)))
fgaColorVector[which(levels(allFga$Cancer) %in% c("OV", "ESCA", "LUSC"))] <- "#E0115F"
fgaColorVector[which(levels(allFga$Cancer) %in% c("mm_hgsc"))] <- "#8B0000"

### pairwise ks test for mouse tumor for most similar
hgscKsStat <- NULL
for (i in unique(allFga$Cancer)) {
  if (i == "mm_hgsc") {
    next()
  } else{
    tmpDf <- allFga[which(allFga$Cancer == i), ]
    tmpKs <- ks.test(tmpDf$fga, allFga$fga[which(allFga$Cancer == "mm_hgsc")], alternative = "two.sided")
    hgscKsStat <- rbind(hgscKsStat, data_frame("sample" = i, "Dstat" = tmpKs$statistic, "pval" = tmpKs$p.value))
  }
}
hgscKsStat$padj <- p.adjust(hgscKsStat$pval)


# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231117fgaCompHgsc.pdf", useDingbats = TRUE,
#     width = 7, height = 4)
ggplot(allFga, aes(x = Cancer, y = fga)) + geom_boxplot(color = fgaColorVector) + 
  ylab("fraction of genome altered (aneuploidy)") + xlab("cancer type/model") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()

### doing for just mouse adeno +  adenocarcinoma + hg coadread

allFgaRed <- allFga[which(allFga$Cancer %in% c("OV", "mm_hgsc")), ]

library(ggpubr)
my_comparisons <- list( c("OV", "mm_hgsc") )
# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231024fgaCompRed.pdf", useDingbats = TRUE,
#     width = 7, height = 4)
ggboxplot(allFgaRed, x = "Cancer", y = "fga") + 
  ylab("fraction of genome altered (aneuploidy)") + xlab("cancer type/model") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+ # Add pairwise comparisons p-value
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()

### for methodlogy figure
### 20231024: creating dummy frequency plot for methodology figure

dummyFreqAmp <- armGisticCoadApc_amp_bed
dummyFreqAmp$Freq[-which(dummyFreqAmp$Chr == "13")] <- 0
dummyFreqDel <- armGisticCoadApc_del_bed
dummyFreqDel$Freq[-which(dummyFreqDel$Chr == "18")] <- 0

# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231024exampleFreq.pdf", useDingbats = TRUE,
#     width = 10, height = 4)
freqPlotv2(dummyFreqAmp, dummyFreqDel,
           main = "example", speciesType = "human")
# dev.off()

### 20231024: single chromosome for method
circosFreqSingleChr(allMouseAneu_adenoCar_amp_bed, allMouseAneu_adenoCar_del_bed, armGisticCoadApc_amp_bed,
                    armGisticCoadApc_del_bed, filename = "20231024coadreadSingleChr13", ref = "human",
                    chromosome = "h_chr13")


### 20231024: human and mouse freq histogram for method

human_freq_10000 <- rnorm(n = 10000, mean = 0.5, sd = 0.1)
mouse_freq_10000 <- rnorm(n = 10000, mean = 0.3, sd = 0.1)

human_freq_10000 <- c(human_freq_10000, 0, 1)
mouse_freq_10000 <- c(mouse_freq_10000, 0, 1)


# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231107exampleHumanPerm.pdf", useDingbats = TRUE,
#     width = 10, height = 4)
ggplot() + geom_histogram(aes(human_freq_10000), bins = 100, color="#cd5c5c", fill="white") +
  geom_vline(xintercept = c(0.9), linetype = "dashed") + 
  geom_vline(xintercept = c(0.8), linetype = "dashed", color = "darkred") + scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(breaks = seq(0, 500, 100)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")) + 
  xlab("Human Frequency 10,000 Permutations") + ylab("Frequency Count")
# dev.off()

# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231107exampleMousePerm.pdf", useDingbats = TRUE,
#     width = 10, height = 4)
ggplot() + geom_histogram(aes(mouse_freq_10000), bins = 100, color="#cd5c5c", fill="white") +
  geom_vline(xintercept = c(0.7), linetype = "dashed") +
  geom_vline(xintercept = c(0.6), linetype = "dashed", color = "darkred") +scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(breaks = seq(0, 500, 100)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")) + 
  xlab("Mouse Freqeuncy 10,000 Permutations") + ylab("Frequency Count")
# dev.off()


### barplot charts for synteny - 20231026 need to redo for first plot
### also add some quick analysis for paper stats i.e number of genes, segments 
### before and after method 

### might be able to replace this input with permutation resutls since those data already have freqs


syntenyDelIdx <- which(hgsc_arm_allSynTable_del$m_freq >  0.1 & hgsc_arm_allSynTable_del$h_freq > 0.1)
syntenyAmpIdx <- which(hgsc_arm_allSynTable_amp$m_freq >  0.1 & hgsc_arm_allSynTable_amp$h_freq > 0.1)

hgsc_arm_allSynTableMouseRef <- hgsc_arm_allSynTable_del
hgsc_arm_allSynTableMouseRef$str <- paste( paste0("m_chr", hgsc_arm_allSynTableMouseRef$m_chr),
                                                    paste0("h_chr", hgsc_arm_allSynTableMouseRef$h_chr))

barplotFractionMouse <- NULL
i <- unique(hgsc_arm_allSynTableMouseRef$m_chr)[1]
for (i in unique(hgsc_arm_allSynTableMouseRef$m_chr)) {
  tmpDf <- hgsc_arm_allSynTableMouseRef[which(hgsc_arm_allSynTableMouseRef$m_chr == i), ]
  tmpDf$length <- tmpDf$m_end - tmpDf$m_start
  for (j in unique(tmpDf$h_chr)) {
    tmpRes <- sum(tmpDf$length[which(tmpDf$h_chr == j)])/sum(tmpDf$length)
    barplotFractionMouse <- rbind(barplotFractionMouse, data.frame("h_chr" = j, "m_chr" = i,"fraction" = tmpRes,
                                                                   "numberRegions" = length(which(tmpDf$h_chr == j))))
  }
}

barplotMChrAmp <- unique(hgsc_arm_allSynTableMouseRef$str[syntenyAmpIdx])
barplotMChrDel <- unique(hgsc_arm_allSynTableMouseRef$str[syntenyDelIdx])


### 15 in both which will be purple

barplotFractionFiltMouseHgsc <- barplotFractionMouse
barplotFractionFiltMouseHgsc$h_chr <- paste0("h_chr",barplotFractionFiltMouseHgsc$h_chr)
barplotFractionFiltMouseHgsc$m_chr <- paste0("m_chr",barplotFractionFiltMouseHgsc$m_chr)
barplotFractionFiltMouseHgsc$str <- paste(barplotFractionFiltMouseHgsc$m_chr,
                                          barplotFractionFiltMouseHgsc$h_chr)
barplotFractionFiltMouseHgsc$color <- "#FFFFFF"
barplotFractionFiltMouseHgsc$color[which(barplotFractionFiltMouseHgsc$str %in% barplotMChrAmp)] <- "#8B0000"
barplotFractionFiltMouseHgsc$color[which(barplotFractionFiltMouseHgsc$str %in% barplotMChrDel)] <- "#00008B"
barplotFractionFiltMouseHgsc$color[which(barplotFractionFiltMouseHgsc$str %in% intersect(barplotMChrAmp, barplotMChrDel))] <- "#800080"
barplotFractionFiltMouseHgsc$alpha <- ifelse(barplotFractionFiltMouseHgsc$color == "#FFFFFF", 0.1, 0.5)

barplotFractionFiltMouseHgsc$h_chr <- factor(barplotFractionFiltMouseHgsc$h_chr, levels = paste0("h_chr", 1:22))
barplotFractionFiltMouseHgsc$m_chr <- factor(barplotFractionFiltMouseHgsc$m_chr, levels = paste0("m_chr", 1:19))



# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231117fractionSyntenyReccurentSignHgsc.pdf", useDingbats = TRUE,
#     width = 10, height = 5)
ggplot(barplotFractionFiltMouseHgsc, aes(fill= h_chr, y= fraction, x = m_chr)) +
  geom_bar(position="stack", stat="identity", color = barplotFractionFiltMouseHgsc$color,
           alpha =  barplotFractionFiltMouseHgsc$alpha, linewidth = 1.5) +
  scale_fill_manual(values = rep("#023020", 22)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) +
  ylab("human chromosome synteny fraction") + xlab("mouse chromsome")
# dev.off()


# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231027fractionSyntenyReccurentSign.pdf", useDingbats = TRUE,
#     width = 10, height = 5)
# ggplot(barplotFractionFiltMouseHgsc, aes(fill= h_chr, y= fraction, x = m_chr)) + 
#   geom_bar(position="stack", stat="identity", color = barplotFractionFiltMouseHgsc$color,
#            alpha =  barplotFractionFiltMouseHgsc$alpha) + 
#   scale_fill_manual(values = colorVector) +
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
#   theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + 
#   ylab("human chromosome synteny fraction") + xlab("mouse chromsome")
# dev.off()



### below is the filtered version after method
hgsc_arm_allSynTable_amp_good$str2 <- apply(hgsc_arm_allSynTable_amp_good[, c("h_chr", "m_chr")],
                                           1, function(x) paste0(x, collapse = "_"))
unique(hgsc_arm_allSynTable_amp_good$str2)

hgsc_arm_allSynTable_del_good$str2 <- apply(hgsc_arm_allSynTable_del_good[, c("h_chr", "m_chr")],
                                            1, function(x) paste0(x, collapse = "_"))
unique(hgsc_arm_allSynTable_del_good$str2)


barplotFractionFiltMouseHgsc2 <- barplotFractionFiltMouseHgsc
barplotFractionFiltMouseHgsc2$str <- paste(barplotFractionFiltMouseHgsc2$h_chr, barplotFractionFiltMouseHgsc2$m_chr)
barplotFractionFiltMouseHgsc2$alpha <- 0.1
barplotFractionFiltMouseHgsc2$alpha[which(barplotFractionFiltMouseHgsc2$str == "h_chr1 m_chr6" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr1 m_chr3" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr1 m_chr8" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr2 m_chr6" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr3 m_chr6" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr3 m_chr3" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr5 m_chr15" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr7 m_chr6" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr8 m_chr3" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr8 m_chr15" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr10 m_chr18" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr10 m_chr8" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr12 m_chr6" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr12 m_chr15" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr9 m_chr19" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr9 m_chr13" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr15 m_chr17" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr16 m_chr7" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr17 m_chr11" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr22 m_chr11")] <- 1

barplotFractionFiltMouseHgsc2$color <- "#FFFFFF"
barplotFractionFiltMouseHgsc2$color[which(barplotFractionFiltMouseHgsc2$str == "h_chr1 m_chr6" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr1 m_chr3" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr1 m_chr8" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr2 m_chr6" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr3 m_chr6" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr3 m_chr3" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr5 m_chr15" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr7 m_chr6" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr8 m_chr3" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr8 m_chr15" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr10 m_chr18" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr10 m_chr8" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr12 m_chr6" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr12 m_chr15")] <- "#8B0000"

barplotFractionFiltMouseHgsc2$color[which(barplotFractionFiltMouseHgsc2$str == "h_chr9 m_chr19" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr9 m_chr13" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr15 m_chr17" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr16 m_chr7" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr17 m_chr11" |
                                            barplotFractionFiltMouseHgsc2$str == "h_chr22 m_chr11")] <- "#00008B"



# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231117fractionSyntenyConMouseRefPermHgsc.pdf", useDingbats = TRUE,
#     width = 10, height = 5)
ggplot(barplotFractionFiltMouseHgsc2, aes(fill= h_chr, y= fraction, x = m_chr)) + 
  geom_bar(position="stack", stat="identity", color = barplotFractionFiltMouseHgsc2$color,
           alpha =  barplotFractionFiltMouseHgsc2$alpha, linewidth = 1.5) + 
  scale_fill_manual(values = colorVector) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + 
  ylab("human chromosome synteny fraction") + xlab("mouse chromsome")
# dev.off()




### specific gene names for nominated regions
### volcano plots with before and after filtering the genes

library(ggrepel)
pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20240107hgscVolBlockAmp10.pdf"),
    useDingbats = FALSE, width = 10, height = 10)
ggplot(lrtTab_bl10, aes(label = ifelse(lrtTab_bl10$alpha == 1, lrtTab_bl10$genes, ""),
                        x = log2FoldChange, y = log10Q)) +
  geom_point(color = lrtTab_bl10$color2, alpha = lrtTab_bl10$alpha) +
  geom_text_repel(max.overlaps = 30) + 
  ylab("-log10QValue") + xlab("log2FC") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  xlim(c(-3, 3))
dev.off()

pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20240107hgscVolBlockAmp26.pdf"),
    useDingbats = FALSE, width = 10, height = 10)
ggplot(lrtTab_bl26, aes(label = ifelse(lrtTab_bl26$alpha == 1, lrtTab_bl26$genes, ""),
                        x = log2FoldChange, y = log10Q)) +
  geom_point(color = lrtTab_bl26$color2, alpha = lrtTab_bl26$alpha) +
  geom_text_repel(max.overlaps = 30) + 
  ylab("-log10QValue") + xlab("log2FC") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  xlim(c(-3, 3))
dev.off()


pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20240107hgscVolBlockAmp38.pdf"),
    useDingbats = FALSE, width = 10, height = 10)
ggplot(lrtTab_bl38, aes(label = ifelse(lrtTab_bl38$alpha == 1, lrtTab_bl38$genes, ""),
                        x = log2FoldChange, y = log10Q)) +
  geom_point(color = lrtTab_bl38$color2, alpha = lrtTab_bl38$alpha) +
  geom_text_repel(max.overlaps = 30) + 
  ylab("-log10QValue") + xlab("log2FC") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  xlim(c(-3, 3))
dev.off()


### for each of these genes, look for overlaps in pathways affected
### done in brahm

library(ggrepel)



### pathways analysis chart -made on brahm
### from table with the most unique pathway changes

blocksOfInterest <- c("block_amp_12","block_amp_13","block_amp_20", "block_amp_22",
                      "block_loss_26","block_loss_27","block_loss_28", "block_loss_29","block_loss_30")

pathwaysOfInterest <- c("HALLMARK_NOTCH_SIGNALING", "HALLMARK_ADIPOGENESIS", "HALLMARK_APOPTOSIS",
                        "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS", "HALLMARK_P53_PATHWAY")
pathwaysOfInterest2 <- c("ACT. NOTCH_SIGNALING", "SUPP. ADIPOGENESIS", "SUPP. APOPTOSIS",
                         "SUPP. HYPOXIA", "SUPP. GLYCOLYSIS", "SUPP. P53_PATHWAY")

dotplotDf_blocks <- NULL
for (i in blocksOfInterest) {
  dotplotDf_blocks <- rbind(dotplotDf_blocks,
                            data.frame("blocks" = rep(i, length(pathwaysOfInterest)),
                                       "hallmarks" = pathwaysOfInterest,
                                       "hallmarks2" = pathwaysOfInterest2))
}
dotplotDf_blocks$str <- paste(dotplotDf_blocks$hallmarks, dotplotDf_blocks$blocks)
dotplotDf_blocks$GeneRatio <- NA
dotplotDf_blocks$p.adjust <- NA
dotplotDf_blocks$count <- NA

dotplotDf_blocks$GeneRatio <- allBlocksGsea$geneRatio[match(dotplotDf_blocks$str, allBlocksGsea$str)]
dotplotDf_blocks$p.adjust <- allBlocksGsea$p.adjust[match(dotplotDf_blocks$str, allBlocksGsea$str)]
dotplotDf_blocks$count <- allBlocksGsea$count[match(dotplotDf_blocks$str, allBlocksGsea$str)]

dotplotDf_blocks$blocks <- factor(dotplotDf_blocks$blocks, levels = blocksOfInterest)
dotplotDf_blocks$hallmarks2 <- factor(dotplotDf_blocks$hallmarks2, levels = pathwaysOfInterest2)
### make block labels same color as synteny block using illustrator

# pdf(file = paste0("/avatar_data5/tmp/kev/misc/20231109aneuploidyGsea", ".pdf"),
#     useDingbats = FALSE, width = 12, height = 5)
ggplot(data = dotplotDf_blocks, aes(x = blocks, y = hallmarks2, 
                                    color = `p.adjust`, size = GeneRatio)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("Gene set enrichment")
# dev.off()

